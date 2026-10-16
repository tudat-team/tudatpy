Remaining LLR physical terms: source audit
=========================================

A remaining centimetre-scale candidate is a difference between the integrated
lunar orientation and the full observation-reduction solution. This audit does
not establish that DE440 omits empirical libration terms and imports no fitted
correction or coefficient from another ephemeris version.

The `DE430 coordinate memo <https://ssd.jpl.nasa.gov/ftp/eph/planets/ioms/de430_moon_coord_iom.pdf>`_
(Williams, Boggs and Folkner 2013), printed p.9, explicitly identifies periodic
longitude-libration effects outside the integration. Its Eq.1 gives
``delta_tau = 5.0 cos(l') + 1.5 cos(2l-2D) - 3.6 cos(2F-2l)`` milliarcseconds,
with quoted fitted uncertainties. These are DE430 values only. The memo gives
4.5 mas RMS and 10 mas peak, equivalent to about 38/84 mm at the lunar equator;
range projections are smaller and reflector dependent. Printed p.12 describes
the supplied Euler angles as numerically integrated angles. Thus the distinction
is explicit for DE430, rather than an inference from its binary file format.

The `JPL LLR Model 2020 <https://ssd.jpl.nasa.gov/ftp/eph/planets/ioms/LLR_Model_2020_DR.pdf>`_
section 4 lists optional cosine/sine fit coefficients at 206 days, one year,
three years and six years, sometimes one month, with further empirical periods
of 1306/1643 days. It does not provide a DE440 coefficient set, exact phase/axis
mapping or instructions establishing whether those corrections are incorporated
in the released orientation arrays. Section 3.11 obtains lunar orientation from
the ephemeris file; listing possible solution parameters does not prove which
were used or how they were exported.

The `DE440 paper <https://ssd.jpl.nasa.gov/ftp/eph/planets/ioms/de440.aj.pdf>`_
(Park et al.2021), section 2.4, describes the stored Euler angles as integrated
simultaneously with the orbit; section 4.1 uses a delayed tidal response. The
official `header.440 <https://ssd.jpl.nasa.gov/ftp/eph/planets/ascii/de440/header.440>`_
contains 645 constants. Identified lunar rotation parameters include
PHI/THT/PSI, mantle/core angles and rates, K2M, TAUM, LBET/LGAM, IFAC, COBLAT,
KVC and PSIDOT. No named amplitudes/phases for the listed empirical periods were
identified. Header-name absence alone cannot prove that exported angles exclude
an already applied correction. The frozen PCK comments identify DE440 PA
orientation; comparison with the DE440 ASCII angle arrays tests representation,
not completeness of a separate LPRED reduction solution.

The missing information is a DE440-specific empirical coefficient/phase set and
an explicit statement of its application/export to the NIO/ASCII/PCK angles.
If it lies outside those arrays, it is a physical reduction/fit availability
limit. It cannot be repaired by changing binary interpolation, transplanting
the DE430 coefficients or estimating a correction from these residuals.

Other effects have different sizes in the same 2020 memo:

* Section 3.12 / Appendix A1: phase-shifted lunar surface tides are about 2 mm,
  with maximum vertical terms about 2.2 mm and horizontal terms a few tenths of
  a millimetre. They are distinct from the integrated tidal torque response.
* Section 3.13: expected Apollo thermal displacement is 1--1.5 mm; Lunokhod
  displacement is about 5 mm vertically, with a horizontal component. These
  uncertain models cannot individually explain several centimetres.
* Section 3.18: off-axis corner-cube path and wavelength-dependent refractive
  index effects are below 1 mm for the stated conditions. They do not justify
  a centimetre reflector offset or a second normal-point centroid correction.
* Sections 3.19 / 4: station/time-span biases and separate green/IR biases are
  fitted or input; typical interstation offsets are a few centimetres. Optical
  path metrology to the telescope-axis intersection is one possible cause.
  These are solution-specific terms, not transferable physical calibrations.

Published postfit residuals include jointly fitted parameters and therefore do
not establish a raw fixed-DE440 residual target. No effect in this note has been
added to the model. The next observational check should use raw time series
split by station, wavelength and reflector, retaining the individual points.

Local primary-source SHA256 identifiers:

* DE430 memo: ``c3e15f2033efaf45889674375da24acdc70eb71ab659a37f8715b01902be56c2``.
* JPL2020 memo: ``68947c087e0f6b8192c56d53d4b9026ff4d1c88621d1f8bd00e9b6a400ee709a``.
* Park2021 PDF: ``96c3945deb12e0f0843c2bf0703d3542d950cf85ff903b1fe11a501634f68ee0``.
* Header440: ``0f4636b663e6f00af9efa2b605ef23635489f56b23dd1a5e3821d3acc6ac3359``.

Original cached references remain external under ``research`` and
``research/llr-physical-terms``; the latter includes the download manifest and
the parsed header audit. No primary source or full software snapshot is copied
into this result package.
