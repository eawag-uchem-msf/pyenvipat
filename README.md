# pyEnvipat

Python port of the R package [`enviPat`](https://github.com/blosloos/enviPat)
(Loos, Gerber et al., EAWAG), licensed GPL-2.0-only, same as upstream --
calculation of isotope fine structure patterns and theoretical mass
spectrometry envelopes, plus tools for subsequent envelope generation,
centroidization and instrument resolution interpolation.

The distributed package is `pyEnvipat`; the importable module keeps the
original name, `envipat` (`import envipat as ep`), matching the R
package's own name.

This port reproduces the R package's algorithms and results (isotope
fine-structure calculation, Gaussian/Cauchy-Lorentz profile convolution,
centroid/intensoid/valley detection, and a from-scratch reimplementation
of R's smoothing-spline resolution interpolation) via a compiled C
extension built from the same, unmodified core C sources the R package's
own `.Call` entry points delegate to
(`envipat._core.enviPat_py` -- see that module's docstring). Parity with
the original R package is verified against fixtures frozen from real R
output; see `tests/` and the port's GitLab issue history for the exact
tolerances and any disclosed R-vs-Python divergences.

## Installation

Not yet published to PyPI -- install from a source checkout:

```bash
git clone git@gitlab.eawag.ch:joboog/pyenvipat.git
cd pyenvipat
pip install .
```

For development (editable install, so local edits take effect without
reinstalling):

```bash
pip install -e .
```

Requirements: Python >= 3.10, a C compiler (to build the bundled
`envipat._core.enviPat_py` extension at install time), and `numpy` at
build time. Runtime dependencies (`pandas`, `numpy`, `scipy`) are
installed automatically.

No R installation is required to install, import or use this package --
the bundled reference data (`envipat.data`) ships as plain CSV/TXT files,
and the compiled extension vendors the original algorithms' C sources
directly. (The R implementation was used once to freeze the reference
outputs under `tests/golden/` and has since been removed from this
repository; it remains in git history.)

## Quickstart

The typical workflow mirrors the R package's own
(`check_chemform` -> `isopattern` -> `envelope` -> `vdetect`): validate
and normalize a chemical formula, calculate its isotope fine structure,
convolute that into a theoretical instrument envelope, then detect
centroids (or intensoids/valleys) in that envelope.

```python
import envipat as ep
import envipat.data as data

# load the reference data bundled with the package
isotopes = data.isotopes

# validate/normalize formulas and get their monoisotopic mass
checked = ep.check_chemform(isotopes, ["C6H12O6", "C8H10N4O2"])
#    warning new_formula  monoisotopic_mass
# 0    False     C6H12O6          180.063388
# 1    False   C8H10N4O2          194.080376

# calculate isotope fine structure patterns
pattern = ep.isopattern(
    isotopes,
    checked["new_formula"],
    threshold=0.1,
    charge=False,
    algo=2,
)
pattern["C6H12O6"].head()
#           m/z   abundance  12C  13C  1H  2H  16O  18O  17O
# 0  180.063388  100.000000    6    0  12   0    6    0    0
# 1  181.066743    6.489437    5    1  12   0    6    0    0
# 2  181.067605    0.228555    6    0  12   0    5    0    1
# 3  181.069665    0.138016    6    0  11   1    6    0    0
# 4  182.067634    1.232996    6    0  12   0    5    1    0

# convert the patterns to envelope profiles at 1 million resolving power
profiles = ep.envelope(
    pattern,
    ppm=False,
    dmz=0.0001,
    frac=0.25,
    env="Gaussian",
    resolution=1e6,
)

# detect centroids in each envelope
centro = ep.vdetect(profiles, detect="centroid")
centro["C6H12O6"]
#           m/z   abundance
# 0  180.063388  100.000000
# 1  181.066743    6.525364
# 2  181.067605    0.229828
# 3  181.069665    0.138785
# 4  182.067634    1.246678
# 5  182.070098    0.177424
```

`isowrap` (`envipat.isowrap.isowrap`) composes the same three steps into
one call for a batch of formulas, given a `resolution` value (or a
per-formula sequence -- e.g. from `envipat.resolution.get_resolution`,
this port's reimplementation of R's `getR`, when interpolating resolution
from an instrument calibration curve such as one of
`envipat.data.resolution_list`'s 31 curves, rather than passing a fixed
value):

```python
resmass = data.resolution_list["Elite_R240000@400"]
resolution = ep.get_resolution(checked, resmass)  # per-formula, interpolated

result = ep.isowrap(
    isotopes,
    checked,
    resolution=resolution,
    threshold=0.1,
    charge=1,
    algo=2,
    detect="centroid",
)
```

## Package layout

The public pipeline functions listed below are also re-exported at the
top level (`import envipat as ep; ep.check_chemform(...)`, etc. -- see
the Quickstart above); the bundled datasets stay a submodule import
(`import envipat.data as data`) since `envipat.data` loads them lazily
on first attribute access, and re-exporting them at the top level would
defeat that.

- `envipat.data` -- lazily-loaded bundled reference datasets
  (`isotopes`, `chemforms`, `adducts`, `resolution_list`).
- `envipat.formula` -- `check_chemform`, `check_ded`, `mergeform`,
  `subform`, `multiform`: chemical formula validation and algebra.
- `envipat.isopattern` -- isotope fine structure calculation.
- `envipat.envelope` -- profile convolution (Gaussian / Cauchy-Lorentz).
- `envipat.vdetect` -- centroid / intensoid / valley detection.
- `envipat.check_several` -- overlap screening between compounds' peaks.
- `envipat.resolution` -- `get_resolution`, a from-scratch port of R's
  `getR` smoothing-spline resolution interpolation.
- `envipat.isowrap` -- convenience wrapper composing `isopattern`,
  `envelope` and `vdetect` for a batch of formulas.
- `envipat._core` -- the compiled C extension wrapping the vendored
  algorithm core; not intended for direct use (see its own docstring).

Every public function above has a full docstring (parameters, return
shape, and where relevant, exactly how it reproduces -- or deliberately
diverges from -- its R original); `help(...)` or your editor's docstring
popup is the fastest way to explore.

## Citation

If you use `envipat`/`enviPat` in your work, please cite the original
algorithm and R package:

Loos, M.; Gerber, C.; Corona, F.; Hollender, J.; Singer, H. (2015).
Accelerated isotope fine structure calculation using pruned transition
trees, *Analytical Chemistry*, 87(11), 5738-5744.
<https://doi.org/10.1021/acs.analchem.5b00941>

## License

GPL-2.0-only, same as the original R package -- see [`LICENSE`](LICENSE).
This is a derivative work of the original `enviPat` R package by Martin
Loos and Christian Gerber (EAWAG); see [`AUTHORS`](AUTHORS) for full
attribution and [`NOTICE`](NOTICE) for how the port relates to the
original C sources it vendors.
