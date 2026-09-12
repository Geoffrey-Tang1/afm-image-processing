"""
jpk_to_png.py  —  JPK .jpk (AFM image) -> publication PNG
    python3 jpk_to_png.py [folder]        # default: folder of this script
Needs: pip install tifffile numpy scipy matplotlib

Pipeline (height, retrace):
    plane fit -> per-line 1st-order fit -> tip-scar / noise-line removal
    -> light Gaussian smoothing (display only)
Sq in the filename is measured BEFORE smoothing (after scar removal),
so the reported roughness is not lowered by the display filter.
Colour bar is zeroed at the 0.5 % height; 99.5 % is the top.
"""
from scipy.ndimage import gaussian_filter, median_filter

import tifffile, numpy as np

def _tv(pg,c,d=None):
    t=pg.tags.get(c); return t.value if t else d

def read_jpk(path):
    tf=tifffile.TiffFile(path)
    p0=tf.pages[0]
    meta=dict(ulen=float(_tv(p0,32834)), vlen=float(_tv(p0,32835)),
              nx=int(_tv(p0,32838)), ny=int(_tv(p0,32839)),
              scanline=_tv(p0,32843), date=_tv(p0,32771))
    chans=[]
    for pg in tf.pages[1:]:
        name=_tv(pg,32848)
        if name is None: continue
        retrace = "retrace : true" in str(_tv(pg,32851,""))
        default=_tv(pg,32897)
        slots={}
        for k in range(6):
            b=32912+48*k
            n=_tv(pg,b)
            if n is None: continue
            slots[n]=dict(parent=_tv(pg,b+2), unit=_tv(pg,b+18),
                          scaling=_tv(pg,b+19), mult=_tv(pg,b+20), off=_tv(pg,b+21))
        data=pg.asarray().astype(np.float64)
        c=slots.get(default)
        unit=None
        if c and c["scaling"]=="LinearScaling":
            data=data*float(c["mult"])+float(c["off"]); unit=c["unit"]
        elif c: unit=c["unit"]
        if meta["scanline"]=="topDown": data=data[::-1]
        chans.append(dict(name=name, retrace=retrace, unit=unit, slot=default, data=data))
    return meta, chans

def _mad(a):
    a=a[np.isfinite(a)]
    return np.median(np.abs(a-np.median(a)))*1.4826+1e-12

def _long_runs(bad, min_run):
    """keep only horizontally elongated runs (real scars / noise lines run along
    the fast-scan direction); isolated small features (pores, grains) are kept."""
    out=np.zeros_like(bad)
    for i in range(bad.shape[0]):
        row=bad[i]
        if not row.any(): continue
        idx=np.flatnonzero(np.diff(np.r_[0,row.view(np.int8),0]))
        for a,b in zip(idx[::2],idx[1::2]):
            if b-a>=min_run: out[i,a:b]=True
    return out

def remove_scars(z, thresh=1.6, max_thick=3, min_run=16):
    """tip scars / noise lines: a run of 1-max_thick lines that jumps away from
    BOTH neighbouring lines (while those agree) over at least min_run pixels
    along the fast-scan direction -> replaced by vertical interpolation."""
    z=z.copy(); ny=z.shape[0]; nfix=0
    for t in range(1,max_thick+1):
        up=z[:ny-t-1]; dn=z[t+1:]
        ref=_mad(up-dn)
        agree=np.abs(up-dn)<1.5*ref
        for k in range(t):
            blk=z[1+k:ny-t+k]
            w=(k+1)/(t+1)
            interp=up*(1-w)+dn*w
            bad=agree&(np.abs(blk-interp)>thresh*ref)
            bad=_long_runs(bad,min_run)
            blk[bad]=interp[bad]; nfix+=bad.sum()
    return z, nfix/z.size

def plane_level(z):
    ny,nx=z.shape; Y,X=np.mgrid[0:ny,0:nx]
    A=np.c_[X.ravel(),Y.ravel(),np.ones(X.size)]
    c,*_=np.linalg.lstsq(A,z.ravel(),rcond=None)
    return z-(A@c).reshape(ny,nx)

def line_level(z, order=1):
    z=z.copy(); x=np.arange(z.shape[1])
    for i in range(z.shape[0]):
        z[i]=z[i]-np.polyval(np.polyfit(x,z[i],order),x)
    return z

def process(z_nm, scar=True, vmed_weight=0.6, smooth=0.8):
    """Full display pipeline on one height channel, in nm.

    Returns (image, sq, scar_fraction). `sq` is the RMS roughness measured
    after artefact removal but BEFORE the cosmetic filters, so it is not
    lowered by them. Pass vmed_weight=0, smooth=0 for unfiltered output.
    """
    z = line_level(plane_level(z_nm), 1)
    frac = 0.0
    if scar:
        z, frac = remove_scars(z)
    sq = z.std()
    if vmed_weight:
        z = z + vmed_weight * (median_filter(z, size=(3, 1)) - z)
    if smooth:
        z = gaussian_filter(z, smooth)
    return z, sq, frac


# ---------------------------------------------------------------- batch
if __name__ == "__main__":
    import sys, os, glob, matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Helvetica", "Nimbus Sans", "Arial", "Liberation Sans"]

    # --- display filtering (does NOT affect the Sq written in the filename) ---
    SCAR_THRESH, SCAR_THICK, SCAR_MINRUN = 1.6, 3, 16   # tip scars / noise lines
    VMED_WEIGHT = 0.6      # 0..1, strength of the 3-px vertical median (kills 1-line stripes)
    SMOOTH_PX   = 0.8      # Gaussian sigma in pixels

    src = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
    out = os.path.join(src, "PNG"); os.makedirs(out, exist_ok=True)
    for p in sorted(glob.glob(os.path.join(src, "*.jpk"))):
        stem = os.path.basename(p)[:-4]
        meta, ch = read_jpk(p)
        h = [c for c in ch if c["name"] == "height" and c["retrace"]] \
            or [c for c in ch if c["name"] == "height"]
        z = line_level(plane_level(h[0]["data"] * 1e9), 1)       # nm
        z, frac = remove_scars(z, SCAR_THRESH, SCAR_THICK, SCAR_MINRUN)
        sq = z.std()                                   # quantitative: before display filter
        lo, hi = np.percentile(z, [0.5, 99.5])
        zf = z + VMED_WEIGHT * (median_filter(z, size=(3, 1)) - z)
        zs = np.clip(gaussian_filter(zf, SMOOTH_PX) - lo, 0, None)
        um = meta["ulen"] * 1e6
        bar = 5 if um > 10 else 1
        fig, ax = plt.subplots(figsize=(4.2, 3.4), dpi=300)
        im = ax.imshow(zs, cmap="afmhot", vmin=0, vmax=hi - lo,
                       extent=[0, um, 0, um], origin="lower", interpolation="nearest")
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values(): s.set_linewidth(0.8)
        ax.plot([um*0.06, um*0.06+bar], [um*0.07]*2, lw=4.5, color="w", solid_capstyle="butt")
        ax.text(um*0.06+bar/2, um*0.125, "%g \u03bcm" % bar, color="w",
                ha="center", va="bottom", fontsize=11, weight="bold")
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
        cb.set_label("Height (nm)", fontsize=10); cb.ax.tick_params(labelsize=9)
        fig.tight_layout()
        fn = "%s_%gx%gum_Sq%.1fnm.png" % (stem, um, meta["vlen"]*1e6, sq)
        fig.savefig(os.path.join(out, fn)); plt.close(fig)
        print("%-78s Sq=%5.2f nm  scar-fixed %4.1f%%" % (fn, sq, frac*100))
