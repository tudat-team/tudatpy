Solid Earth gravity tide references
===================================

The native tests contain four frozen published COST-G force samples and their
inputs, not predictions from Tudat. The immutable comparison criterion is
1e-11 m/s^2 (vector norm), from the COST-G software-comparison publication:
https://adgeo.copernicus.org/articles/55/1/2020/.
``sources.json`` pins URLs, SHA256 and sizes for every input. Tests perform no
downloads and need neither SPICE nor a Python kernel at runtime.

Sample extraction and clocks
-----------------------------

``costg_samples.txt`` selects zero-based rows 0, 720, 1440 and 2879 from the
published 2880-row 2008-07-03 files. Orbit, rotation, Doodson arguments and forces
are instrument-file lines 6, 726, 1446 and 2885. Corrected EOP matrix inputs are
lines 3, 723, 1443 and 2882. Each line retains its source sample index and GPS MJD.
The rotation direction is independently checked against both published ICRF and
ITRF positions. Sun/Moon vectors are geometric Earth-centered DE421 positions
in J2000, metres, with SPICE aberration setting NONE. GPS -> TT adds 51.184 s;
NAIF's standard DELTET constants convert TDT to ET before DE421 evaluation.

The benchmark's supplied Doodson arguments use UTC in its GMST helper, as
explicitly documented in its README. Injecting those published arguments tests
the harmonic synthesis; it does not establish production TT/UT1 clock wiring.
Production factories separately require the selected TDB-input Earth rotation
and explicit TT/UT1 conversion. The solid-pole test uses the published corrected
pole, with m1=xp-xbar and m2=-(yp-ybar) in arcseconds and the existing shared
2018 secular-pole helper evaluated at TT. The small GPS/TT epoch difference in
the benchmark pole description is immaterial to its published force budget;
no clock approximation is substituted into production.

The response source constants are explicitly GM=3.986004415e14 m^3/s^2 and
R=6378136.3 m (IERS section 6.1). Pinned GROOPS constants define these defaults;
its Earth-tide evaluator produces source coefficients before converting to the
requested field constants. The benchmark EIGEN-6C4 header supplies consuming
GM=3.986004415e14 and R=6378136.46. No benchmark configuration override has been
verified; this distinction is recorded rather than inferring source radius from
the static field. The frozen force comparisons test the resulting contract.
The benchmark nominal field is tide-free, so the complete permanent solid tide
is retained. A separate analytic control removes
4.4228e-8*(-0.31460)*0.30190 from C20 exactly once for zero-tide output.

Frequency-table transcription
-----------------------------

``iers2010SolidTideFrequencyTerms.inc`` contains only the six Doodson multipliers
and final ip/op amplitudes (units 1e-12) from IERS tables 6.5a/b/c: 48 diurnal,
21 zonal and two semidiurnal terms. ``extract_tables.py`` compares every row with
the primary chapter-6 PDF text, then checks the frozen C++ transcription.
The same 71 rows independently match the older ILR numerical files. Table 6.3
Love numbers are directly checked against the primary table; the degree-four
response uses degree-two forcing, not degree-four direct forcing.

The data transcription source is Orekit 13.1.8, under Apache License 2.0.
Original notices and licence are retained in ``third_parties/orekit-data`` and
identified in the repository NOTICE. No Orekit evaluator code is copied.
The small evaluation structure adapts Tudat's original ILR solid-tide code and
reuses the current SolidBodyTide geometry/normalized coefficient helper. It
corrects the old ILR diurnal S21 sign to equation 6.8b's C-iS convention, excludes
its implicit secular drifts, fixes the old pole constant to -1.333e-9 and returns
the declared 1x1 degree-2/order-1 pole block. GROOPS source is used only to verify
benchmark constants/conventions; no GPL source is copied.

Offline reproduction
--------------------

With the pinned source files already downloaded into a reference directory,
run from the repository root (Python with NumPy and SpiceyPy for DE421 only)::

    pdftotext -layout "$REF/iers_chapter6_2018.pdf" "$REF/iers_chapter6_2018.txt"
    python tests/test_tudat/reference_files/earth_gravity_solid/extract_tables.py "$REF" --check
    python tests/test_tudat/reference_files/earth_gravity_solid/generate_benchmark.py "$REF" --check

Omit ``--check`` to regenerate; ordinary native tests only read the frozen text.
No Tudat imports occur in either generator. The numerical controls retain source
and consumer constants independently and verify force invariance under output
GM/radius changes, normalized degree-2-to-4 signs, degree-three response,
permanent-tide treatment and pole coefficient/block conventions.
