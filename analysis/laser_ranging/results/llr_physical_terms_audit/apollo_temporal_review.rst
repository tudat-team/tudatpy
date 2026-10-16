Independent APOLLO temporal review
==================================

Read-only inspection of the 69 raw points and full-span plot in revision
``bd4778ec2``, ``results/apollo_long_arc_20170301_20170430``. The CSV SHA256 is
``2bd0b44740387095565d383dfaaf798a0372b23d6b3c67daef83acaba1114e6e``.
No model run, fit, correction, cut or replacement of a residual was performed.

The curves contain shared temporal variation and persistent reflector
separation. They do not support a single monotonic drift. Unweighted raw nightly
mean changes are:

==========  =====================  =====================
Reflector   March 9 to 22 (mm)      March 22 to April 3 (mm)
==========  =====================  =====================
A11         +78.17                 -63.17
A14         +81.72                 -61.78
A15         +86.63                 -56.17
L1          +86.78                 -57.20
==========  =====================  =====================

These parallel changes are much larger than the within-night scatter visible
in the plot. They cannot be explained by a static reflector bias alone. L2,
which was not observed on March 9, also falls 54.26 mm from March 22 to April 3.

Differences between raw nightly means, with no fitted offsets, show a separate
reflector structure:

===================  =====================  =============
Pair                 Difference range (mm)  Shared nights
===================  =====================  =============
L1 minus A15         +72.37 to +76.84        5
A11 minus A15        -50.15 to -28.52        7
A14 minus A15        +51.49 to +75.69        8
L2 minus A15         -66.30 to -46.79        4
===================  =====================  =============

These are descriptive comparisons of all raw points on shared nights, not
simultaneous measurements or adjusted data. Sparse counts, differing observing
times and provider uncertainties prevent interpreting them as exact offsets.
L1 and A15 nevertheless track particularly closely despite their large
separation in absolute residual.

The archived geometry spans mean Earth-direction PA longitude about -7.38 to
+5.51 degrees and latitude -6.65 to +5.67 degrees, lunar phase angle 35--135
degrees and transmit elevations about 27--69 degrees. These quantities covary
with the nine sampled nights. The data therefore do not distinguish lunar
libration geometry, station geometry or atmospheric dependence independently.
Similar phase angles alone do not determine the pattern: A15's March 22 and
April 21 means differ by 56.08 mm although their phase angles are about 110.5
and 114.2 degrees. The latter is a single point and is retained unchanged.

The warranted conclusion is shared curved variation plus reflector-dependent
structure. A monthly period, a unique physical attribution and an intrinsic
noise-level residual have not been established. Same-reflector observations
from Grasse will provide an independent station comparison when available.
