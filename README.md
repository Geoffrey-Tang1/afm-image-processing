# afm-image-processing

Batch conversion of JPK (Bruker/JPK NanoWizard) AFM image files into
publication-ready height maps, with no GUI and no proprietary software.

A `.jpk` image file is a TIFF whose IFDs carry the channel data and whose
private tags (32768+) carry the scan geometry and the raw→physical
calibration chain. This reads those tags directly, so a full 12-channel
image stack can be decoded and rendered with nothing but `tifffile`,
`numpy`, `scipy` and `matplotlib`.

## Install

```bash
pip install tifffile numpy scipy matplotlib
```

Helvetica is used when available and falls back to Nimbus Sans / Arial /
Liberation Sans, so figures look the same on machines without it.

## Use

```bash
python3 jpk_to_png.py                 # every .jpk next to the script
python3 jpk_to_png.py path/to/folder  # every .jpk in that folder
```

PNGs are written to a `PNG/` subfolder, one per input file, named

```
<original name>_<scan size>_Sq<roughness>.png
```

so the scan size and RMS roughness are readable from the file listing
without opening anything. The scan size comes from the file's own
metadata, not from the file name — useful, because an operator who
re-scans a smaller area under the previous name is a common occurrence.

As a library:

```python
from jpk_to_png import read_jpk, plane_level, line_level, remove_scars

meta, channels = read_jpk("scan.jpk")
print(meta["ulen"], meta["nx"])          # scan width in m, pixels
for c in channels:
    print(c["name"], c["retrace"], c["unit"], c["data"].shape)
```

`read_jpk` returns every channel (height, measuredHeight, error,
amplitude, phase, and the auxiliary channels) for both trace and
retrace, each already converted to physical units by applying the
calibration slot the file itself marks as default.

## Pipeline

Applied to the height retrace channel:

1. **Plane fit** — removes sample tilt.
2. **Per-line first-order fit** — removes the line-to-line offset and
   slope that scanner drift leaves behind.
3. **Scar / noise-line removal** — a run of 1–3 scan lines that departs
   from *both* neighbouring lines while those two agree is a tip scar or
   a skipped line, and is replaced by interpolation between them. The
   flagged pixels must form a run of at least 16 pixels along the fast
   scan direction: real scars and line noise are horizontally elongated,
   whereas pores, grains and particles are not, so small round features
   survive a filter that would otherwise erase them.
4. **Vertical median (weight 0.6)** — suppresses single-line striping
   that is too shallow to trigger step 3.
5. **Gaussian blur (σ = 0.8 px)** — cosmetic only.

Steps 4 and 5 are display filters. **The Sq written into the file name
is measured after step 3 and before step 4**, so the reported roughness
is never lowered by cosmetic smoothing. Colour scale runs from the
0.5th percentile (drawn as 0 nm) to the 99.5th percentile.

Tune at the top of `__main__`:

| Constant | Default | Effect |
| --- | --- | --- |
| `SCAR_THRESH` | 1.6 | Lower catches fainter lines, at the risk of eating real texture |
| `SCAR_THICK` | 3 | Maximum scar thickness in scan lines |
| `SCAR_MINRUN` | 16 | Minimum horizontal run, in pixels, to count as a line artefact |
| `VMED_WEIGHT` | 0.6 | 0 disables the vertical median; 1 applies it fully |
| `SMOOTH_PX` | 0.8 | Gaussian σ in pixels; 0 disables |

Set `VMED_WEIGHT = 0` and `SMOOTH_PX = 0` for unfiltered output.

## Line profiles

`line_profiles.py` takes random line profiles across the same images:

```bash
python3 line_profiles.py                      # 3 profiles per .jpk here
python3 line_profiles.py data -n 5            # 5 per image
python3 line_profiles.py data --match porous  # only matching file names
python3 line_profiles.py data --raw           # skip the display filters
```

Each profile is a whole scan row, which is how an AFM acquires the data in
the first place — no interpolation, and no mixing of adjacent lines. Row
positions are drawn from a seed derived from the file name, so a re-run
reproduces the same lines; `--seed N` draws a different set.

Into a `Profiles/` subfolder it writes, per image, a `_linemap.png` showing
where the lines were taken (greyscale, with the lines drawn in colour and
numbered, so the height maps stay clean) and a `_profiles.png` of the traces
themselves. Lines within one image are lightness steps of one hue rather than
dashed variants of one colour. Two CSVs cover the whole batch: every sampled
point, and one row per line with mean, Rq, Ra and peak-to-valley.

Profiles carry the same zero reference as the height map's colour bar, so the
two can be read against each other, and they are deliberately not clipped at
zero — clipping would draw a flat floor that is not in the data.

## Limitations

Tip-shape artefacts are not corrected. A blunt or contaminated tip that
smears every feature in the fast-scan direction, or a groove ploughed by
debris dragged across the surface, cannot be undone by filtering — those
images need a new tip, not better processing. Filtering removes
artefacts narrower than the real features; it cannot recover information
the tip never resolved.

Only image files are handled. Force-curve files (`.jpk-force`,
`.jpk-qi-data`) are ZIP containers with a different layout and are not
read.

## License

MIT.
