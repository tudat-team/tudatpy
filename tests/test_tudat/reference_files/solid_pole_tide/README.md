The eight frozen native cases reproduce the IERS chapter 7 document dated
1 February 2018. `source.json` pins the exact primary PDF bytes reused from the
existing reference cache; the original retrieval timestamp was not retained.
The PDF itself identifies that document date. Its primary URL, SHA256 and
equations are recorded, so a later update at the same URL cannot silently
replace the fixture reference.

Reproduce and check offline, using a local copy matching the pinned hash:

```
python3 tests/test_tudat/reference_files/solid_pole_tide/generate_reference.py /path/to/icc7-2018.pdf --check
```

Omit `--check` to print the six ENU/Cartesian components for each original row.
The standalone scalar calculation uses the published secular pole (Eq.21),
wobble signs (Eq.25), original colatitude/south/east/radial formula (Eq.26) and
transposed Cartesian projection (Eq.27–28). It reads only the input and frozen
output literals from the native test, never the Tudat implementation or kernel.
The eight expected vectors and their existing 1e-14 metre tolerance are unchanged.
Python's standard double-precision trigonometric functions suffice for that
absolute tolerance; no independently executed full software model is claimed.
