"""
line_profiles.py  —  random line profiles from JPK AFM images

    python3 line_profiles.py [folder] [-n 3] [--match REGEX] [--raw]

For every .jpk image in the folder it picks N random horizontal rows
(whole scan lines, the way an AFM actually acquires them), and writes into
a Profiles/ subfolder:

    <name>_linemap.png    the image in greyscale with the chosen lines drawn
                          on it in colour and numbered - kept separate from
                          the height maps so the figure-ready images stay clean
    <name>_profiles.png   the N profiles, one colour per line
    AFM_line_profiles.csv          long table: one row per sampled point
    AFM_line_profiles_summary.csv  one row per line: mean, Rq, Ra, peak-to-valley

Row positions come from a seed derived from the file name, so re-running
reproduces exactly the same lines.

Heights use the same pipeline as jpk_to_png.py (plane + per-line levelling,
scar removal, then the display filters) and the same zero reference as that
script's colour bar, so a profile can be read directly against the height
map. Profiles are NOT clipped at zero - a clipped trace would show a flat
floor that is not in the data - so a few points can sit slightly below 0.
Pass --raw to skip the cosmetic filters and profile the levelled data.

Requires jpk_to_png.py next to this file.
"""
import os, re, sys, glob, csv, hashlib, argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import Rectangle
from matplotlib.colors import to_rgb, to_hex
from scipy.ndimage import gaussian_filter, median_filter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from jpk_to_png import read_jpk, plane_level, line_level, remove_scars

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Helvetica", "Nimbus Sans", "Arial", "Liberation Sans"]

# One base colour per image; the N lines inside an image are lightness steps
# of it. Replace with your own scheme, or map file names to colours in
# base_colour() below.
DEFAULT_COLOURS = ["#0072b2", "#d55e00", "#009e73", "#cc79a7", "#5b4ea8", "#b8860b"]

VMED_WEIGHT, SMOOTH_PX = 0.6, 0.8        # keep in step with jpk_to_png.py


# ---------------------------------------------------------------- colour
def _to_lab(hex_colour):
    rgb = np.array(to_rgb(hex_colour))
    lin = np.where(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)
    xyz = np.array([[0.4124, 0.3576, 0.1805],
                    [0.2126, 0.7152, 0.0722],
                    [0.0193, 0.1192, 0.9505]]) @ lin / np.array([0.95047, 1.0, 1.08883])
    f = np.where(xyz > 0.008856, np.cbrt(xyz), 7.787 * xyz + 16 / 116)
    return 116 * f[1] - 16, 500 * (f[0] - f[1]), 200 * (f[1] - f[2])


def shift_lightness(hex_colour, dL):
    """Same hue, L* moved by dL - gives distinguishable lines without dashes."""
    L, a, b = _to_lab(hex_colour)
    L = float(np.clip(L + dL, 8, 95))
    fy = (L + 16) / 116; fx = fy + a / 500; fz = fy - b / 200
    finv = lambda t: t ** 3 if t ** 3 > 0.008856 else (t - 16 / 116) / 7.787
    xyz = np.array([finv(fx), finv(fy), finv(fz)]) * np.array([0.95047, 1.0, 1.08883])
    lin = np.array([[3.2406, -1.5372, -0.4986],
                    [-0.9689, 1.8758, 0.0415],
                    [0.0557, -0.2040, 1.0570]]) @ xyz
    lin = np.clip(lin, 0, 1)
    srgb = np.where(lin <= 0.0031308, 12.92 * lin, 1.055 * lin ** (1 / 2.4) - 0.055)
    return to_hex(np.clip(srgb, 0, 1))


def base_colour(stem, index):
    """Override this to key the colour off your own sample naming."""
    return DEFAULT_COLOURS[index % len(DEFAULT_COLOURS)]


# ---------------------------------------------------------------- data
def height_map(path, raw=False):
    """Returns (metadata, height in nm with the colour-bar zero, colour-bar span)."""
    meta, channels = read_jpk(path)
    ch = ([c for c in channels if c["name"] == "height" and c["retrace"]]
          or [c for c in channels if c["name"] == "height"])
    z = line_level(plane_level(ch[0]["data"] * 1e9), 1)
    z, _ = remove_scars(z)
    lo, hi = np.percentile(z, [0.5, 99.5])
    if not raw:
        z = z + VMED_WEIGHT * (median_filter(z, size=(3, 1)) - z)
        z = gaussian_filter(z, SMOOTH_PX)
    return meta, z - lo, hi - lo


def pick_rows(stem, ny, n, margin=0.06, seed=None):
    key = stem if seed is None else "%s:%d" % (stem, seed)
    rng = np.random.default_rng(int(hashlib.md5(key.encode()).hexdigest()[:8], 16))
    lo, hi = int(margin * ny), int((1 - margin) * ny)
    return sorted(rng.choice(np.arange(lo, hi), size=n, replace=False).tolist())


# ---------------------------------------------------------------- figures
def _scalebar(ax, um, height_frac=0.017):
    bar = 5 if um > 10 else 1
    ax.add_patch(Rectangle((um * 0.06, um * 0.062), bar, um * height_frac,
                           facecolor="w", edgecolor="k", linewidth=1.8,
                           joinstyle="miter", zorder=6))
    ax.text(um * 0.06 + bar / 2, um * 0.125, "%g μm" % bar, color="w",
            ha="center", va="bottom", fontsize=11, weight="bold",
            path_effects=[pe.withStroke(linewidth=1.8, foreground="k")], zorder=6)


def draw_linemap(z, vmax, um, rows, colours, fname):
    ny = z.shape[0]
    fig, ax = plt.subplots(figsize=(4.2, 3.4), dpi=400)
    im = ax.imshow(np.clip(z, 0, None), cmap="Greys_r", vmin=0, vmax=vmax,
                   extent=[0, um, 0, um], origin="lower", interpolation="nearest")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_linewidth(0.8)
    for k, (r, c) in enumerate(zip(rows, colours), 1):
        y = (r + 0.5) * um / ny
        ax.axhline(y, color=c, lw=1.6, zorder=4)
        ax.text(um * 0.985, y, str(k), color=c, ha="right", va="bottom",
                fontsize=10, weight="bold", zorder=5,
                path_effects=[pe.withStroke(linewidth=1.8, foreground="k")])
    _scalebar(ax, um)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label("Height (nm)", fontsize=10); cb.ax.tick_params(labelsize=9)
    fig.tight_layout(); fig.savefig(fname); plt.close(fig)


def draw_profiles(x, profiles, rows, colours, um, ny, fname, label=""):
    fig, ax = plt.subplots(figsize=(4.6, 3.4), dpi=400)
    for k, (p, c, r) in enumerate(zip(profiles, colours, rows), 1):
        ax.plot(x, p, color=c, lw=1.1,
                label="%d   y = %.2f μm" % (k, (r + 0.5) * um / ny))
    ax.set_xlabel("x (μm)", fontsize=11)
    ax.set_ylabel("Height (nm)", fontsize=11)
    ax.tick_params(labelsize=10); ax.set_xlim(0, um)
    for s in ax.spines.values(): s.set_linewidth(0.8)
    ax.legend(fontsize=8.5, frameon=False, loc="upper left",
              bbox_to_anchor=(0, 1.30), handlelength=1.6)
    if label:
        ax.text(1.0, 1.02, label, transform=ax.transAxes,
                ha="right", va="bottom", fontsize=9)
    fig.tight_layout(); fig.savefig(fname); plt.close(fig)


# ---------------------------------------------------------------- batch
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folder", nargs="?", default=None, help="folder of .jpk files")
    ap.add_argument("-n", type=int, default=3, help="profiles per image (default 3)")
    ap.add_argument("--match", default=None,
                    help="only files whose name matches this regex (case-insensitive)")
    ap.add_argument("--seed", type=int, default=None,
                    help="change to draw a different set of rows")
    ap.add_argument("--raw", action="store_true",
                    help="profile the levelled data without the display filters")
    args = ap.parse_args()

    src = args.folder or os.path.dirname(os.path.abspath(__file__))
    out = os.path.join(src, "Profiles"); os.makedirs(out, exist_ok=True)

    files = sorted(glob.glob(os.path.join(src, "*.jpk")))
    if args.match:
        files = [f for f in files if re.search(args.match, os.path.basename(f), re.I)]
    if not files:
        print("no .jpk files matched"); return

    points, stats = [], []
    for idx, path in enumerate(files):
        stem = os.path.basename(path)[:-4]
        meta, z, vmax = height_map(path, args.raw)
        ny, nx = z.shape
        um, vm = meta["ulen"] * 1e6, meta["vlen"] * 1e6
        rows = pick_rows(stem, ny, args.n, seed=args.seed)
        base = base_colour(stem, idx)
        span = 16 * (args.n - 1) / 2 if args.n > 1 else 0
        colours = [shift_lightness(base, span - 32 * k / max(args.n - 1, 1))
                   for k in range(args.n)]
        x = (np.arange(nx) + 0.5) * um / nx
        profiles = [z[r] for r in rows]

        draw_linemap(z, vmax, um, rows, colours, os.path.join(out, stem + "_linemap.png"))
        draw_profiles(x, profiles, rows, colours, um, ny,
                      os.path.join(out, stem + "_profiles.png"),
                      "%g×%g μm" % (um, vm))

        for k, (r, p) in enumerate(zip(rows, profiles), 1):
            y = round((r + 0.5) * um / ny, 4)
            points += [[stem, round(um, 3), k, r, y, round(xi, 5), round(zi, 4)]
                       for xi, zi in zip(x, p)]
            stats.append([stem, round(um, 3), k, r, y,
                          round(p.mean(), 3), round(p.std(), 3),
                          round(np.abs(p - p.mean()).mean(), 3),
                          round(p.max() - p.min(), 3),
                          round(p.min(), 3), round(p.max(), 3)])
        print("%-64s %2gum  rows %s" % (stem[:64], um, rows))

    with open(os.path.join(out, "AFM_line_profiles.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sample", "scan_um", "line", "row_index", "y_um", "x_um", "height_nm"])
        w.writerows(points)
    with open(os.path.join(out, "AFM_line_profiles_summary.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sample", "scan_um", "line", "row_index", "y_um", "mean_nm",
                    "Rq_nm", "Ra_nm", "peak_to_valley_nm", "min_nm", "max_nm"])
        w.writerows(stats)
    print("\n%d images, %d profiles, %d points -> %s"
          % (len(files), len(stats), len(points), out))


if __name__ == "__main__":
    main()
