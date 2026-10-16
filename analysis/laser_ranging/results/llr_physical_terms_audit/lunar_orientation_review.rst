Full DE440 lunar orientation convention control
===============================================

This source review and two matrix evaluations test whether a fixed lunar-frame
rotation could have entered the ASCII-to-PCK-to-native construction. No yaw,
offset or other correction was fitted; no observation calculation was run.

The `DE440 paper <https://ssd.jpl.nasa.gov/ftp/eph/planets/ioms/de440.aj.pdf>`_
equation 8 explicitly maps PA coordinates to ICRF3 as
``Rz(-phi) Rx(-theta) Rz(-psi)``, using the passive right-handed matrices
defined in equations 9--11. The original DE440 ASCII quantity 13 was evaluated
independently with 80-digit Chebyshev arithmetic and these matrices, retaining
the exact double epoch used by the native model.

The pinned binary lunar PCK has two type-2 segments, each with frame class ID
31008 and base-frame ID 1 (SPICE J2000), not ECLIPJ2000. The matched
`frame kernel <https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/satellites/a_old_versions/moon_de440_200625.tf>`_
defines MOON_PA_DE440 as class 2 / 31008. Its ICRF discussion explicitly says
that these binary lunar PCKs use the label J2000 for ICRF data and that users
normally should not add a J2000-to-ICRF transformation. The separate
MOON_ME_DE440_ME421 fixed rotation is not used by this analysis.

Native ``SpiceRotationalEphemeris::getRotationToBaseFrame`` requests the
target-to-base SPICE rotation directly. ``computeRotationQuaternionBetweenFrames``
calls ``pxform_c``, copies its matrix without transposition and converts it to
an Eigen quaternion. There is no added ecliptic or dynamical-J2000 bias in this
path. The analysis-only ICRF alias is an explicit identity relative to J2000.

At the first archived launch epochs of the frozen APOLLO and Matera cases:

========  ======================  =====================  ==================
Case      TDB seconds from J2000  Angular difference(rad) Surface bound (um)
========  ======================  =====================  ==================
APOLLO    542294424.1855052        5.7760e-14              0.1420
Matera    714322074.1868316        3.0143e-13              0.7409
========  ======================  =====================  ==================

Angular difference is the small-angle norm ``norm(R_ascii-R_native)/sqrt(2)``.
The conservative surface bound is the Frobenius norm times 1738000 m, rather
than a fitted line-of-sight projection. Native versus direct SPICE maximum
matrix-element differences are 1.11e-16 at both epochs; ICRF-alias versus
J2000 matrices are exactly identical. Native orthogonality errors are at most
3.33e-16. These differences are over 300000 times smaller than 20 mas.

Earlier independent full PCK-angle-to-SPICE matrix controls already gave
submicrometre A15 lever-arm agreement at all eleven APOLLO and two Matera
launch epochs. The separate original-ASCII/PCK angle comparison gave norms
below 4.17e-13 rad. The present two controls join both checks directly to the
actual native body rotation, rather than relying on Euler coefficient equality
alone. They close this implementation-frame hypothesis at the tested epochs;
they do not establish whether every empirical term of a published fitted LLR
reduction is included in the exported integrated orientation.

Provenance
----------

* Native kernel SHA256:
  ``cab49e738cfcf0ea5f251cb140ab8919d5372dbc3e5ba09314d0d5be4eedfcef``.
* Original ``ascp01950.440`` SHA256:
  ``e24f6356225f7d912e35af2ebc5f9bcf6fb0425c53aa92b5d7243cb063477286``.
* Original ``header.440`` SHA256:
  ``0f4636b663e6f00af9efa2b605ef23635489f56b23dd1a5e3821d3acc6ac3359``.
* PCK, frame-kernel and identity-alias hashes remain pinned in
  ``llr_dataset.json``; no provider kernel was modified.
* External probe ``verification/audit_de440_full_orientation.py`` SHA256:
  ``1a20a3fc55f0986643446e7c5edbad07e673a8c8ebfe7a891caa82927f66fad7``.

Native source locations are ``src/tudat/interface/spice/spiceRotationalEphemeris.cpp``
and ``src/tudat/interface/spice/spiceInterface.cpp``. The external retained
numeric report is ``research/llr-physical-terms/full_orientation_audit.json``.
Its SHA256 is
``e58f389e8d98486be429bc17a8f7dd1bb10e79782b1be0da0fe7a6f07f56a3a3``.
