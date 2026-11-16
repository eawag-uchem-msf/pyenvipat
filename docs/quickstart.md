# Quickstart

The typical workflow mirrors the R package:

1. {func}`~envipat.formula.check_chemform` validates and normalizes formulas and
   computes their monoisotopic mass,
2. {func}`~envipat.isopattern.isopattern` calculates the isotope fine structure,
3. {func}`~envipat.envelope.envelope` convolutes it into a theoretical instrument
   profile,
4. {func}`~envipat.vdetect.vdetect` detects centroids (or intensoids/valleys) in that
   profile.

```python
import envipat as ep
import envipat.data as data

isotopes = data.isotopes

checked = ep.check_chemform(isotopes, ["C6H12O6", "C8H10N4O2"])
pattern = ep.isopattern(isotopes, checked["new_formula"], threshold=0.1, algo=2)
profiles = ep.envelope(pattern, dmz=0.0001, resolution=1e6)
centroids = ep.vdetect(profiles, detect="centroid")
centroids["C6H12O6"]
```

{func}`~envipat.isowrap.isowrap` runs steps 2-4 in one call. Its `resolution` can be a
single value or one value per formula, for example interpolated from an
instrument calibration curve with {func}`~envipat.resolution.get_resolution`:

```python
resmass = data.resolution_list["Elite_R240000@400"]
resolution = ep.get_resolution(checked, resmass)

centroids = ep.isowrap(isotopes, checked, resolution=resolution, charge=1)
```

Every pipeline function returns a `dict` keyed by formula. A formula whose
calculation failed maps to the string `"error"` instead of a DataFrame, as in
the R package, so check with `isinstance(result[formula], str)` before using a
value.

The {doc}`tutorial` runs the same pipeline step by step with plots.
