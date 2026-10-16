Independent two-station temporal and source review
=================================================

Read-only inspection of all 306 raw points, their combined plot and source CRD
records from the completed March--April 2017 run. CSV SHA256:
``e9b3e9da081d5e12f8d10f7cf9579e250a9fe8c3ef2400722b1c36792816e908``.
No model rerun, rejection, reweighting, fitted correction or adjusted residual
was used. All 69 APOLLO rows reproduce the previous package exactly.

One isolated large discrepancy dominates the Grasse RMS
------------------------------------------------------

The full 237-point Grasse mean is +10.007 mm and RMS about zero 359.583 mm.
One A15 infrared record at 2017-03-09 18:34:15.297622200200 UTC has residual
+5.42994684 m and contributes 96.2155% of the full sum of squared residuals.
All other raw Grasse values lie between -160.75 and +175.51 mm. This identifies
the arithmetic source of the large aggregate RMS; it does not establish that
the observation is invalid or authorize its removal.

The exact record is ``apollo15_201703.npt`` line 171, session 29:
``11 66855.297622200200 2.490835100745 me09 2 ...``. The tag and RTT match the
provider bytes, also present in the earlier cached current-release file. Its
26 returns, 361.9 ps bin RMS, quality flag 0, me09 configuration and 1064.2 nm
wavelength do not identify a distinct correction regime. The 13 other normal
points in this same session span -71.30 to -36.70 mm; adjacent source lines
168/174 give -47.25/-67.86 mm. No station calibration was reapplied by the driver.

The source H4 interval is 18:09:48--18:21:08, while records continue through
22:29:05. This is a concrete metadata inconsistency affecting many following
records, not solely the large point. The driver uses the printed NP seconds
of day and the same session date; it does not use this end time as a tag.
No published specific reason for the large discrepancy was recovered. A model
failure at this isolated epoch is not excluded by this source-only review.

The `provider CRD file <https://edc.dgfi.tum.de/pub/slr/data/npt_crd/apollo15/2017/apollo15_201703.npt>`_
SHA256 is ``cd546e1d1f0947962366ab8441091a639233fae77ccd6dcb98e8053bf0963eb6``.
`EDC validation documentation <https://edc.dgfi.tum.de/en/oc/crd/v1/>`_
defines WGL021 as missing station XYZ for an orbit-fit check; this would not
diagnose a physically bad range. Bin RMS is distribution metadata, not a
provided formal normal-point standard error; dividing by sqrt(return count)
would require an unsupported noise assumption here.

All-point wavelength and reflector statistics
--------------------------------------------

Means and RMS about zero below are in mm and retain every point. Wavelength
groups sample different nights and reflectors, so their means are not a fitted
colour calibration.

==========  =========  =====  ==========  ==========
Wavelength  Reflector  Count  Mean        RMS
==========  =========  =====  ==========  ==========
1064.2 nm   A11        17      -53.156      75.286
1064.2 nm   A14        20      +50.205      72.334
1064.2 nm   A15        110     +17.631      521.338
1064.2 nm   L1         40      +66.285      85.545
1064.2 nm   L2         37      -55.354      74.883
532.1 nm    A14        2       -76.168      113.825
532.1 nm    A15        7       -41.123      42.000
532.1 nm    L1         4       +42.147      42.839
==========  =========  =====  ==========  ==========

The full infrared group has 224 points, mean +11.800 mm/RMS 369.595 mm; green
has 13 points, mean -20.893 mm/RMS 59.226 mm. On April 4, green/IR A15 means
are -42.17/-44.71 mm and L1 means +42.15/+37.55 mm. Green A14 has one -160.75 mm
record with six returns and another +8.42 mm record; same-night IR A14 spans
+3.67 to +39.57 mm. This is not evidence for a uniform green/IR offset, nor
does the low return count establish a formal uncertainty or rejection rule.

Shared temporal and reflector structure
--------------------------------------

The combined curves show similar curved temporal changes at both stations,
with A11/L2 lower, A14/L1 higher and A15 intermediate. On March 22, Grasse
minus APOLLO raw nightly means are +6.76, +19.39, +10.22, +6.93 and +14.99 mm
for A11, A14, A15, L1 and L2, respectively. Their absolute nightly means span
about -27 to +119 mm. On April 6/7, the corresponding A11/A14/A15 differences
are -26.56/-14.95/-18.37 mm; on April 14 they are -20.52/-13.00/-9.25 mm and
on April 15 -9.77/-15.30/-11.51 mm. These comparisons are separated by about
4.7--8 hours, so they are not simultaneous station biases.

The isolated metre-scale record and the smooth centimetre-scale structure
coexist. A station-specific constant alone cannot explain the common curved
variation and reflector ordering. The sparse, non-simultaneous geometry does
not uniquely identify a lunar or terrestrial physical cause. The published
fitted station/colour biases remain excluded, and intrinsic measurement noise
has not been reached or inferred.
