Grasse 1064.2 nm source and dispersion check
==========================================

This is a source/numerical diagnostic, with no native observation calculation,
model adjustment or substitution of a nominal wavelength.

The `INPOP17a report <https://www.imcce.fr/content/medias/recherche/equipes/asd/inpop/inpop17a_v3_090817.pdf>`_
Table 4, printed page 14, distinguishes the old Grasse 1064.0 nm Nd:YAG series
from the 2015--2017 MeO 1064.2 nm series. Section 2.1.1 gives reported wavelength
resolution as 0.1 nm. The frozen
`March 2017 Apollo 11 CRD <https://edc.dgfi.tum.de/pub/slr/data/npt_crd/apollo11/2017/apollo11_201703.npt>`_
likewise records 1064.2 in its C0/C1/C2 records. Neither source establishes a
formal plus/minus 0.5 nm model-domain rule. The proposed 1064.5 nm guard is an
engineering interpretation of the nominal 1064 nm limit, not a published
extension of its accuracy validation. An explicit vacuum-versus-air wavelength
definition was not recovered from these records.

`IERS chapter 9 <https://iers-conventions.obspm.fr/chapter9.php>`_ gives the
nominal 0.355--1.064 micrometre range. Its original
`FCUL_ZD_HPA.F <https://iers-conventions.obspm.fr/content/chapter9/software/FCUL_ZD_HPA.F>`_
evaluates the dispersion formula directly without a wavelength guard;
``FCUL_A.F`` has no wavelength argument. This supports evaluating the declared
1064.2 nm line rather than clamping it to 1064.0 nm; it does not independently
validate atmospheric accuracy outside the published nominal range.

An independent evaluation of IERS equations 9.3--9.10 used the actual
interpolated weather from all 237 frozen Grasse normal points. For the declared
1064.2 minus 1064.0 nm sensitivity, the per-leg path difference is:

====================  ======================
Elevation scenario    Path difference
====================  ======================
90 degrees            -11.34 to -11.56 um
20 degrees            -32.86 to -33.49 um
3 degrees             -0.1668 to -0.1699 mm
====================  ======================

These are weather-envelope scenarios, not the measured point elevations; the
3-degree value is a boundary diagnostic. Weather includes both green and IR
sessions, while the wavelength pair is held fixed for this sensitivity. The
first IR point uses 882.83345 hPa, 287.50701 K and relative humidity 0.340226;
its corresponding values are -11.549 um, -33.464 um and -0.169374 mm. No range
residual or source wavelength was changed.

Retained provenance
-------------------

Primary files remain in the external research cache. SHA256 values:

* INPOP17a PDF:
  ``1c73a95d11e1899f22a8b61a904eb1d200be1985d8d46b2b187927bc31cb318a``.
* Frozen ``apollo11_201703.npt``:
  ``9fae3c39f1e6563f0f371dff328a9e159b68bec886f0430239d1a24c939350d5``.
* IERS chapter 9 PDF:
  ``4bd0ef02008f02218af51eb1a293bda9849ca14f40e684a10fa80a6149d01903``.
* Official ``FCUL_ZD_HPA.F`` / ``FCUL_A.F``:
  ``92731affca053aad15a44be7db58dbf6df689e75cf2e1f3b39cb4d99a4da198b`` /
  ``fdeb39aee3c8d4c2d6eb6a7e743c420372e28da5b3e84942d09580a88847693a``.
* Independent formula source ``verification/llr_optical_actual_input_audit.py``:
  ``195727e384ad608055f191335cb8d3348a132fdcaf664218f2aa6b95f19fd1d2``.

Probe/output hashes and the active model/manifest hashes are recorded in the
adjacent runtime note. Both diagnostics used kernel
``80b924287707124a8efea66455b26233b5511b7856983d35c05d21f82734f215``.
