Fixed-event Shapiro audit
========================

This bounded audit retains all eleven saved APOLLO events from the historical
INPOP17a station-profile case. It reconstructs actual deformed surface endpoints
without another range solve. It uses the recorded ``10f8307b`` kernel and source
``a8d22d03a``; later ICRF/selected-clock experiments are separate configurations.
``summary.json`` records complete hashes; ``components.csv`` contains all 66
body/leg components. Delays are ``c * delta_t`` per leg, in metres; a round-trip
one-way range is half the sum of both legs.

The native formula agrees with `JPL2020 section 3.17, equation 19
<https://ssd.jpl.nasa.gov/ftp/eph/planets/ioms/LLR_Model_2020_DR.pdf>`_
and `IERS2010 section 11.2, equation 11.17
<https://iers-conventions.obspm.fr/content/chapter11/tn36_c11.pdf>`_. Sun position
is evaluated at each leg's midpoint, Earth at its station endpoint, and Moon at
the reflection event. The summed correction reconstructs saved native components
within 3.26 nm. The forward logarithm differs from an 80-digit exact-input
calculation by at most 0.91 pm. Squared-distance subtraction occurs in a partial
gradient, not this forward delay.

Sun, Earth and Moon delays are respectively 7.361--7.377 m, 36.99--37.55 mm and
0.589--0.608 mm per leg. Freezing a body at the alternative transmit/mid/receive
epochs changes the Earth term by at most 49.55 micrometres, the Moon term by
2.12 micrometres and the solar term by 0.244 nm. These epoch alternatives and the
4.14 micrometre summed ``v/c`` delay scale are diagnostics, not an exact
moving-monopole solution. No moving-body velocity factor or retarded epoch is
implemented by this native model. IERS equation 11.3 applies to an incoming VLBI
ray; it must not be substituted blindly for a finite Earth--Moon segment.
JPL2020 does not specify gravitating-body evaluation epochs in equation 19.
The omitted planets sum to at most 1.816 mm one-way here. This audit identifies
no missing centimetre-scale delay, and applies no residual subtraction or tuning.

``audit_executed.py.txt`` preserves the exact executed Python source bytes. It
requires the recorded local cache paths and historical environment script; it is
an execution record, not a new public test or current-default analysis entry
point. The larger per-event JSON and provider inputs remain external at the
paths recorded in ``summary.json``. Source code and selected derived results are
committed; provider products are not redistributed.
