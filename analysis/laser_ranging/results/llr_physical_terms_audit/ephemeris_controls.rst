Coherent modern ephemeris controls: proposal
===========================================

INPOP17a G+IR is the preferred first comparison: its station solution was
independently preselected, its orbit and PA rotation are publicly distributed
together, and the original reference table permits an explicit reconstruction
of the misprinted Lunokhod rows. This is a proposed fixed published solution,
not a new fit or an exact replay of its historical reduction. No implementation
or residual calculation was performed for this audit.

The subsequently reviewed implementation and both predeclared optical
sensitivities are retained in the `306-point results <../inpop17a_long_arc_20170301_20170430/README.md>`_.
They preserve the source limitations below and do not change the DE440 default.

INPOP17a reflector table repair for review
----------------------------------------

The `2017 report <https://www.imcce.fr/content/medias/recherche/equipes/asd/inpop/inpop17a_v3_090817.pdf>`_
Table 9, printed p.22, explicitly identifies INPOP13c as the reference values.
Its Lunokhod 1 baseline instead repeats terrestrial MLRS2 coordinates, and its
Lunokhod 2 baseline repeats the known Lunokhod 1 vector. These are visible in
the PDF, not text-extraction errors. The cited
`INPOP13c report <https://www.imcce.fr/content/medias/recherche/equipes/asd/inpop/inpop13c.pdf>`_
Table 13, printed p.27, provides all five original baselines unambiguously.

The proposed repair uses those original baselines and each named reflector's
printed G+IR correction; it does not swap correction rows or estimate offsets.
There is no recovered author erratum, so this reconstruction remains explicit
and requires design review. Vectors below are XYZ in metres::

  name  INPOP13c baseline                         INPOP17a G+IR correction
  A11   [1591924.511, 690802.582, 21003.774]       [42.1831, -102.9494, -0.1035]
  A14   [1652725.840, -520890.307, -109730.480]    [-36.0924, -107.0555, -0.1287]
  A15   [1554674.570, 98196.294, 765005.696]       [3.8898, -100.5489, -0.6184]
  L1    [1114345.496, -781226.597, 1076059.335]    [-52.9913, -71.6464, -0.8250]
  L2    [1339314.148, 801958.776, 756359.229]      [49.3269, -86.6622, -0.6982]

  name  reconstructed INPOP17a G+IR PA XYZ
  A11   [1591966.6941, 690699.6326, 21003.6705]
  A14   [1652689.7476, -520997.3625, -109730.6087]
  A15   [1554678.4598, 98095.7451, 765005.0776]
  L1    [1114292.5047, -781298.2434, 1076058.5100]
  L2    [1339363.4749, 801872.1138, 756358.5308]

The quoted 2017 correction uncertainties (XYZ, metres) are respectively
[.0181,.0052,.0118], [.0196,.0140,.0028], [.0155,.0068,.0067],
[.0149,.0042,.0021] and [.0145,.0015,.0083]. They are solution uncertainties,
not normal-point noise. These vectors must accompany INPOP17a PA rotation;
substituting DE440 coordinates or its ME/PA transform would mix realizations.

Actual provider kernels, constants and time convention
-----------------------------------------------------

The `official release <https://www.imcce.fr/recherche/equipes/asd/inpop/download17a>`_
provides the TDB archive ``inpop17a_TDB_m100_p100_spice.tar.gz`` (24189292 bytes),
retrieved after review. Its actual frame name is ``MOON_PA_INPOP``; the prose
mentions a versioned name, but the loaded assignments use this unversioned name.
Its PCK segment is type 2, class ID 1900301, base-frame ID 1 (SPICE J2000/ICRF).
SPK and PCK coverage is JD 2414105.0--2488985.0. Earth is a direct SSB segment;
Moon has a geocentric segment. Both are supplied by this same integration.

The official meta-kernel order is ``.tf, .tpc, .bpc, _time.bsp, .bsp``. A
comparison should clear the pool, load unchanged shape/leap-second kernels,
then those INPOP products and the identity ICRF alias. Do not retain DE440
translation, lunar rotation or mass kernels. The INPOP TPC supplies EMB, Moon
and Sun GMs but no BODY399_GM: Earth must explicitly use EMB minus Moon.

The original TDB header has AU=149597870.700 km and EMRAT=81.30057186659006.
Converting its au^3/day^2 mass parameters with 86400 seconds/day gives, in SI::

  GM_Sun   = 132712440042692823665.7546055 m^3/s^2
  GM_Earth = 398600436213176.9063798473   m^3/s^2
  GM_Moon  = 4902799907327.332799680577   m^3/s^2

For comparison, the provider TPC prints these values in km^3/s^2::

  BODY10_GM  = 132712440042.6928
  BODY3_GM   = 403503.2361205042
  BODY301_GM = 4902.799907327331
  Earth     := BODY3_GM - BODY301_GM = 398600.436213176869

Last-digit differences from the header conversion are rounding, at relative
order 1e-16; they are not alternative fitted constants. BODY399_GM is absent.

The header's RMOON is 1738 km. Table 7 gives h2=.0439 for G+IR; Table 6 fixes
l2=.0107. Its LB=1.550519768e-8 and LG=6.969290134e-10 are explicit header
constants. Select the TDB product; changing to TCB would also require compatible
coordinate/mass scaling. The time SPK (1000000000 to 1000000001) represents the
published time difference; loading it alone does not replace the native clock
converter. Keep the current observable clock path explicit in a first control.

The reduction reference `Manche (2011) <https://www.imcce.fr/content/medias/publications/publications-recherche/theses-habilitations/docs/Manche_These.pdf>`_
section 9.3.3, p.115, Eq.9.32 uses LB rather than JPL's LM in the lunar spatial
map. Their difference gives at most 1.183 mm at 1738 km. This is a documented
reference-model distinction, not evidence that every 2017 implementation detail
is unchanged. Section 9.3.1, pp.112--113, retains the full nearly permanent Earth
tide plus the Sun tide and uses actual reflector radius to the fourth power.
It gives a degree-2 surface model; degree-3 gravity in the orbital dynamics does
not establish h3/l3 surface terms. Section 9.3.2 includes spin-variation surface
motion below 0.2 mm. Missing exact 2017 h3/l3/spin details must remain declared;
blindly importing JPL's degree-3 surface constants would not establish a matched
reduction. No permanent-tide subtraction or newer TL convention is justified.

For an equation-faithful first control, use degree 2 and LB, explicitly omitting
the below-0.2-mm spin term and any unestablished h3/l3 terms. Actual radii are
1735472.585, 1736335.929, 1735477.069, 1734928.427 and 1734638.761 m for
A11/A14/A15/L1/L2. The respective ``(radius/1738000)^4`` factors are .99419584,
.99617565, .99420612, .99294951 and .99228653. At a declared Earth distance of
384400 km with aligned radial tide, using 1738 km would increase displacement
by 3.33, 2.19, 3.32, 4.04 and 4.42 mm. These are illustrative formula differences,
not computed range changes or bounds at every observation. The header's gravity
reference radius alone does not override the explicit surface equation.

Predeclared bounded observation comparison
-----------------------------------------

After source/design review, use the first chronological frozen A11, A14 and A15
point for each station on March 22, and APOLLO April 7 / Grasse April 6: twelve
points spanning the already observed crest/trough and opposite lunar longitudes
plus A15 latitude sensitivity. This rule precedes any INPOP residual evaluation.
Keep original identities, clocks, measured weather, all corrections and solver
tolerances; report individual changes, not an improved subset RMS. Only then
consider the complete unchanged 306-point arc, retaining the large source
discrepancy. Station/EOP/loading mismatches still prevent a postfit-noise claim.

Other modern controls and unresolved physical terms
--------------------------------------------------

* DE430 has public matched orbit/libration ASCII plus its coordinate memo and
  reflector table. The 1950--2050 ASCII block is about 29.3 MB; no DE430 lunar
  PCK was identified in the current/old generic NAIF directories. A faithful
  reader/export path and its own PA coordinates/constants are required. Its
  memo explicitly separates periodic longitude corrections from integrated
  Euler angles, so it cannot silently be treated as a complete raw reduction.
* DE441 is not a newer modern-data improvement: the
  `JPL export description <https://ssd.jpl.nasa.gov/planets/eph_export.html>`_
  removes core/mantle damping to support long historical integration. Its
  public header includes librations, but no matched DE441 PCK or reflector fit
  was identified. A DE441 SPK with DE440 PCK is not a coherent control.
* INPOP21a has public orbit/libration products, five sensible PA reflector
  coordinates and explicit longitude-correction amplitudes in Tables 6/7.
  Its station epoch and correction-export semantics remain unresolved; the
  2024 reference-system review still points back to Table 7 rather than fixing
  the missing station metadata. It is a useful later control once clarified.
* No DE440-specific public empirical libration coefficient/phase/export set was
  recovered beyond the sources already audited. The integrated angles agree
  through native rotation; the remaining question concerns reduction-solution
  availability. Slow empirical libration terms can produce shared temporal
  changes and reflector-dependent projection, but this is a hypothesis, not a
  correction inferred from the current residuals.
* `Kopeikin (2009) <https://arxiv.org/abs/0902.3489>`_ emphasizes consistent
  relativistic frames and gauge handling. Its criticism of older Newtonian
  rotation models does not prove a DE440 omission: JPL2020 explicitly says
  geodetic precession has been included in integrated librations since 2014.
* A concrete remaining light-propagation audit is the finite-endpoint moving
  monopole expression in `Soffel and Han (2015), Eq.33 <https://arxiv.org/html/1409.3743v2>`_.
  It changes both endpoint-relative vectors and the aberrated direction
  ``g=n-v_body/c``; current native Shapiro freezes one body position per leg.
  Prior endpoint-epoch and v/c diagnostics are not a completeness proof. A
  fixed-event Eq.33 comparison, with acceleration remainder and full two-leg
  accounting, would quantify the missing term without changing the model.
  No centimetre contribution has yet been established. Turyshev et al.2025
  Eqs.36/40 bound the spatial acceleration remainders at micrometres, and do
  not justify modifying lunar clock/frame conventions in a fixed JPL reduction.

Retained source provenance
-------------------------

External ``research/llr-ephemeris-controls`` retains source manifests, headers,
archive/member hashes and the reconstruction arithmetic. Visual controls are
``inpop13c-table13.png`` there and ``llr-physical-terms/inpop17a-table9.png``.
Key SHA256 identifiers:

* INPOP17a archive: ``52afaf265f3d1dd4d03b8e2239574af3d3d8c65c0b57e5a8c56fed0b137b98b2``.
* INPOP17a header: ``2082297553208d173d3dee0b76ed14cb9955e2a4590cb0ba63b7e09e37aee818``.
* 2017 report: ``1c73a95d11e1899f22a8b61a904eb1d200be1985d8d46b2b187927bc31cb318a``.
* INPOP13c report: ``976b2a0dac29960863d74f23af231f9d6c706452de318418462d04da0646c054``.
* Manche thesis: ``f36cf9602a79d4de0d98e1136e9cf5d45db3f9e319bfa0933b8f7a856cec1cb6``.
* Published 2018 article HTML: ``cc720f368ca8a192caba6833011855c799390be7b21c5e4751e1a06906768332``.
* Moving-multipole paper: ``0102ed0b52db706484432fe1482e260a8e15dd47cdc9f3ca58ad190aed9c637b``.

The `2018 published article <https://doi.org/10.1093/mnras/sty096>`_ explicitly
excludes DE430's extra longitude terms; its cached arXiv v1 predates that final
wording. Its empirical radiation/thermal adjustment is a separate few-mm
solution-level term and is not imported here. No new model, fit or default was
selected by residual reduction.
