"""
domain_pore_analysis.py  -  roughness, domains and pores in JPK AFM height images

    python3 domain_pore_analysis.py pure   --preset deep    scan1.jpk [scan2.jpk ...]
    python3 domain_pore_analysis.py blend  --preset deep    scan1.jpk ...
    python3 domain_pore_analysis.py porous --preset shallow scan1.jpk ... [--deq 27.0 26.0 ...] [--png]

Two parameter sets are built in: "deep" for pores about 70 nm deep with raised rims around the
features, and "shallow" for pores about 20 nm deep. Every parameter can also be given on its own
(--margin, --rim, --seed, --split). Results are appended to domain_pore_summary_<film>.csv in the
current folder (one row per image); --png also writes <name>_segmentation.png with the objects
outlined.
Needs jpk_to_png.py from this repository and: pip install tifffile numpy scipy scikit-image matplotlib

Pipeline (one height image, z-sensor channel "measuredHeight", retrace; heights in nm):

 1. Leveling (blend and porous films). Plane fit and per-line median alignment, then five rounds of
    a per-line first-order fit and a plane fit computed only on the continuous matrix phase. The
    features (domains above, or pores below, the Otsu threshold) are excluded together with a
    margin around them (--margin, nm), which also keeps raised rims around the
    features out of the reference. The median height of the reference surface is set to 0 nm.
    Pure films: second-order polynomial background and per-line first-order fit on the surface
    without particles (higher than the median + 6 robust standard deviations, grown by 1 pixel).
 2. Roughness Sq = standard deviation of the height over all pixels (pure films: particles excluded;
    the value with particles is also reported).
 3. Objects.
      blend  - domains: pixels above the Otsu threshold, objects of <= 4 pixels removed.
      porous - pores: every pixel more than --rim nm below the film top that is connected to a pixel
               more than --seed nm below it (hysteresis threshold), objects of <= 4 pixels removed,
               holes filled. Touching pores are separated by a marker-controlled watershed on the
               height: a marker is set at every minimum of the (Gaussian-smoothed, sigma = 1 pixel)
               height that is at least --split nm deep relative to its surroundings (h-minima), so
               two pores are only split when each has its own floor.
    For every object of area A not touching the image border: D_eq = 2 (A/pi)^1/2; reported are the
    number of objects, the area-weighted mean of D_eq (weight A) and the median of D_eq.
 4. Pores.
      depth of a pore - height of the film top minus the 5th percentile of the heights in the pore.
      t               - median depth of the pores at least MIN_WIDTH_NM wide (D_eq), which are wide
                        enough for the tip to reach the floor (error: standard deviation over them).
      V_pore          - sum of the depth below the film top over the pore mask grown by 2 pixels (to
                        include the rim taper) x pixel area: pore volume within the image area A.
      porosity        - V_pore / (A t) and, with the UV-vis film thickness d_eq (--deq, nm),
                        V_pore / (A (d_eq + V_pore/A)).
 5. Errors of D_eq, V_pore and porosity: standard deviation over the 3 x 3 sub-areas of the image
    (objects assigned to a sub-area by their centroid).
"""
import os, sys, csv, argparse
import numpy as np
from scipy import ndimage as ndi
from skimage.filters import threshold_otsu, apply_hysteresis_threshold
from skimage.morphology import h_minima, remove_small_objects
from skimage.segmentation import watershed
from skimage.measure import label, regionprops

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from jpk_to_png import read_jpk, remove_scars

PRESETS = {                      # margin (nm), rim (nm), seed (nm), split (nm)
    "deep":    dict(margin=400, rim=10, seed=40, split=10),   # pores ~70 nm deep, raised rims
    "shallow": dict(margin=80,  rim=5,  seed=10, split=5),    # pores ~20 nm deep
}
CHANNEL = "measuredHeight"
MIN_PX = 4                       # objects of up to this many pixels are removed
MIN_WIDTH_NM = 150               # pores used for the depth t
PARTICLE_K, PARTICLE_GROW = 6, 1 # pure films: particle threshold (robust SD) and growth (px)
TILES = 3


# ------------------------------------------------------------------ input and leveling
def load_height(path):
    meta, ch = read_jpk(path)
    for name in (CHANNEL, "height"):
        c = [c for c in ch if c["name"] == name and c["retrace"]] or [c for c in ch if c["name"] == name]
        if c:
            return meta, c[0]["data"] * 1e9
    raise ValueError(f"{path}: no height channel")


def plane(z, mask=None):
    ny, nx = z.shape; y, x = np.mgrid[0:ny, 0:nx]
    m = np.ones_like(z, bool) if mask is None else mask
    c, *_ = np.linalg.lstsq(np.c_[x[m], y[m], np.ones(m.sum())], z[m], rcond=None)
    return z - (c[0] * x + c[1] * y + c[2])


def median_rows(z):
    return z - np.median(z, axis=1, keepdims=True)


def line_fit(z, mask=None):
    """first-order fit of every scan line on the pixels in mask (all pixels where fewer than 20)"""
    m = np.ones_like(z, bool) if mask is None else mask.copy()
    m[m.sum(1) < 20] = True
    w = m.astype(float); x = np.arange(z.shape[1], dtype=float)[None, :]
    S0, Sx, Sxx = w.sum(1), (w * x).sum(1), (w * x * x).sum(1)
    Sz, Sxz = (w * z).sum(1), (w * x * z).sum(1)
    b = (S0 * Sxz - Sx * Sz) / (S0 * Sxx - Sx ** 2); a = (Sz - b * Sx) / S0
    return z - (a[:, None] + b[:, None] * x)


def features(z, film, grow_px):
    zz = -z if film == "porous" else z
    return ndi.binary_dilation(zz > threshold_otsu(zz), iterations=grow_px) if grow_px else zz > threshold_otsu(zz)


def level_film(z0, film, grow_px, n_iter=5):
    z = median_rows(plane(z0))
    for _ in range(n_iter):
        ref = ~features(z, film, grow_px)
        z = plane(line_fit(z, ref), ref)
    ref = ~features(z, film, grow_px)
    return z - np.median(z[ref])


def poly_background(z, keep, order=2):
    ny, nx = z.shape; Y, X = np.mgrid[0:ny, 0:nx]; X = X / nx; Y = Y / ny
    terms = [X ** i * Y ** j for i in range(order + 1) for j in range(order + 1 - i)]
    c, *_ = np.linalg.lstsq(np.stack([t[keep] for t in terms], 1), z[keep], rcond=None)
    return z - sum(ci * t for ci, t in zip(c, terms))


def level_pure(z0):
    keep = np.ones_like(z0, bool)
    z = median_rows(plane(z0))
    for _ in range(3):
        med = np.median(z[keep]); sd = 1.4826 * np.median(np.abs(z[keep] - med))
        keep = ~ndi.binary_dilation(z > med + PARTICLE_K * sd, iterations=PARTICLE_GROW)
        z = remove_scars(line_fit(poly_background(z0, keep), keep))[0]
    return z - np.median(z[keep]), keep


# ------------------------------------------------------------------ objects
def segment(z, film, p):
    if film == "blend":
        mask = remove_small_objects(z > threshold_otsu(z), max_size=MIN_PX)
        return mask, label(mask)
    mask = apply_hysteresis_threshold(-z, p["rim"], p["seed"])
    mask = ndi.binary_fill_holes(remove_small_objects(mask, max_size=MIN_PX))
    markers = label(h_minima(ndi.gaussian_filter(z, 1), p["split"]) & mask)
    lab = watershed(z, markers, mask=mask)
    rest = label(mask & (lab == 0)); lab[rest > 0] = rest[rest > 0] + lab.max()   # pores without a marker
    return mask, lab


def measure_objects(z, lab, film, px_nm):
    """area (um^2), D_eq (um), depth or height (nm, relative to the film top at 0), border flag, centroid"""
    n, m = z.shape; out = []
    for r in regionprops(lab):
        v = z[r.coords[:, 0], r.coords[:, 1]]
        h = -np.percentile(v, 5) if film == "porous" else np.percentile(v, 95)
        edge = r.bbox[0] == 0 or r.bbox[1] == 0 or r.bbox[2] == n or r.bbox[3] == m
        out.append((r.area * px_nm ** 2 * 1e-6, r.equivalent_diameter_area * px_nm * 1e-3, h, edge, *r.centroid))
    return np.array(out, dtype=float).reshape(-1, 6)


def tiles(n):
    e = np.linspace(0, n, TILES + 1).astype(int)
    return [(e[i], e[i + 1], e[j], e[j + 1]) for i in range(TILES) for j in range(TILES)]


# ------------------------------------------------------------------ analysis
def analyse(path, film, p, deq_nm=None, png=False):
    meta, z0 = load_height(path)
    n = z0.shape[0]; px_nm = meta["ulen"] * 1e9 / z0.shape[1]; A_um2 = z0.size * (px_nm * 1e-3) ** 2
    res = dict(file=os.path.basename(path), film=film)
    if film == "pure":
        z, keep = level_pure(z0)
        res.update(Sq_nm=z[keep].std(), Sq_all_pixels_nm=z.std(), particle_area_pct=100 * (~keep).mean())
        return res
    z = level_film(z0, film, int(round(p["margin"] / px_nm)))
    res["Sq_nm"] = z.std()
    mask, lab = segment(z, film, p)
    o = measure_objects(z, lab, film, px_nm)
    inner = o[o[:, 3] == 0]; A, D = inner[:, 0], inner[:, 1]
    aw = []
    for a, b, c, d in tiles(n):
        s = (inner[:, 4] >= a) & (inner[:, 4] < b) & (inner[:, 5] >= c) & (inner[:, 5] < d)
        if s.any(): aw.append(np.average(D[s], weights=A[s]))
    res.update(n_objects=len(inner), D_eq_aw_um=np.average(D, weights=A), D_eq_aw_sd_um=np.std(aw, ddof=1),
               D_eq_median_um=np.median(D))
    if film == "porous":
        wide = inner[inner[:, 1] * 1e3 >= MIN_WIDTH_NM, 2]
        t = np.median(wide)
        dep = np.clip(-z, 0, None) * ndi.binary_dilation(mask, iterations=2)
        VA = dep.mean()                                                # V_pore / A, nm
        VA_t = np.array([dep[a:b, c:d].mean() for a, b, c, d in tiles(n)])
        res.update(pore_depth_t_nm=t, pore_depth_sd_nm=np.std(wide, ddof=1), n_pores_wide=len(wide),
                   V_pore_um3=VA * A_um2 * 1e-3, V_pore_sd_um3=np.std(VA_t, ddof=1) * A_um2 * 1e-3,
                   V_pore_over_A_nm=VA, porosity_vs_t_pct=100 * VA / t,
                   porosity_vs_t_sd_pct=np.std(100 * VA_t / t, ddof=1))
        if deq_nm:
            res.update(d_eq_nm=deq_nm, film_thickness_nm=deq_nm + VA,
                       porosity_vs_film_pct=100 * VA / (deq_nm + VA),
                       porosity_vs_film_sd_pct=np.std(100 * VA_t / (deq_nm + VA_t), ddof=1))
    if png:
        save_overlay(z, lab, os.path.splitext(os.path.basename(path))[0] + "_segmentation.png")
    return res


def save_overlay(z, lab, path):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from skimage.segmentation import find_boundaries
    lo, hi = np.percentile(z, [0.5, 99.5])
    fig, ax = plt.subplots(figsize=(4, 4), dpi=300)
    ax.imshow(z, cmap="afmhot", vmin=lo, vmax=hi, interpolation="nearest")
    ax.contour((lab > 0).astype(float), levels=[0.5], colors="cyan", linewidths=0.4)
    cut = find_boundaries(lab, mode="inner") & ~find_boundaries(lab > 0, mode="inner")
    ax.imshow(np.ma.masked_where(~cut, cut), cmap="gray", vmin=0, vmax=1, interpolation="nearest")
    ax.set_axis_off(); fig.savefig(path, bbox_inches="tight", pad_inches=0); plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description="Roughness, domains and pores in JPK AFM height images")
    ap.add_argument("film", choices=["pure", "blend", "porous"])
    ap.add_argument("files", nargs="+")
    ap.add_argument("--preset", choices=sorted(PRESETS), help="built-in parameter set: deep (pores ~70 nm, raised rims) or shallow (pores ~20 nm)")
    ap.add_argument("--margin", type=float, help="margin around the features excluded from leveling, nm")
    ap.add_argument("--rim", type=float, help="pore opening: depth below the film top, nm")
    ap.add_argument("--seed", type=float, help="a pore must reach this depth, nm")
    ap.add_argument("--split", type=float, help="watershed: depth of a pore's own floor, nm")
    ap.add_argument("--deq", type=float, nargs="+", help="UV-vis film thickness d_eq (nm), one per porous file")
    ap.add_argument("--png", action="store_true", help="write <name>_segmentation.png")
    a = ap.parse_args()
    p = dict(PRESETS.get(a.preset, PRESETS["shallow"]))
    p.update({k: getattr(a, k) for k in ("margin", "rim", "seed", "split") if getattr(a, k) is not None})
    if a.deq and len(a.deq) != len(a.files):
        ap.error("give one --deq value per file")
    rows = [analyse(f, a.film, p, a.deq[i] if a.deq else None, a.png) for i, f in enumerate(a.files)]
    keys = list(dict.fromkeys(k for r in rows for k in r))
    out = f"domain_pore_summary_{a.film}.csv"; new = not os.path.exists(out)
    with open(out, "a", newline="") as fh:
        w = csv.DictWriter(fh, keys, extrasaction="ignore")
        if new: w.writeheader()
        w.writerows({k: (round(float(v), 4) if isinstance(v, (float, np.floating)) else v) for k, v in r.items()} for r in rows)
    for r in rows:
        print(", ".join(f"{k}={float(v):.4g}" if isinstance(v, (float, np.floating)) else f"{k}={v}" for k, v in r.items()))


if __name__ == "__main__":
    main()
