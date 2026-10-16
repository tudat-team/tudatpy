Independent forward UTC calculation: bounded formula review
===========================================================

This is a proposal only; no second solver or additional observation runs were
implemented. Primary sources: Williams and Boggs, JPL LLR Model 2020,
https://ssd.jpl.nasa.gov/ftp/eph/planets/ioms/LLR_Model_2020_DR.pdf,
sections 3.7, 3.10-3.18, 3.20-3.22; IERS2010 chapter10,
https://iers-conventions.obspm.fr/content/chapter10/tn36_c10.pdf;
SOFA time-scale cookbook, https://www.iausofa.org/s/sofa_ts_c.pdf.
The cached JPL PDF/text were inspected; the public JPL URL was inaccessible
through the web reader during this review.

An independent forward control can use three frozen observations, unchanged
station/reflector conventions, and direct source kernels. Parse transmit UTC
independently, map it through ERFA UTC→TAI→TT→TDB with the selected station,
then solve the two finite-endpoint TDB light-time equations explicitly.
For each leg, tau = norm(x_receiver(t_receiver)-x_transmitter(t_transmitter))/c
+ atmospheric delay + Shapiro delay + corner-cube delay. Iterate the lunar
bounce first, then the Earth reception. The down-leg initial guess can simply
reuse the up-leg duration: observed RTT must not enter the final prediction.
JPL2020 section3.21 uses observed RTT only as an initial estimate.

The predicted atomic RTT is tau_up+tau_down+q(t3)-q(t1), where
q(t)=TAI(t)-TDB(t) at the terrestrial station. Equivalently convert the solved
reception TDB to TAI and subtract transmit TAI using split arithmetic.
For these 2017 records, no leap occurs during an RTT, so this equals the UTC
interval. This is the same physical identity as JPL2020 section3.22; neither
route calls for a second epoch conversion or a second station clock term.
If integrating geocentric time instead of using SOFA, derive the TDB/TT rate
from IERS10.1/10.3/10.4 with LC=(LB-LG)/(1-LG), external potential excluding
Earth, and the endpoint station v_E dot r/c² term. The ambiguous LB/external-U
wording of JPL2020 Eq4 should not be copied as an isolated rate formula.

A useful independence boundary is explicit: separate NumPy/high-precision
geometry, iteration and ERFA clocks can audit native model assembly. Reusing
native station displacement still makes this a hybrid assembly audit; a fully
independent physical model also requires independent IERS Earth displacements,
JPL lunar elastic displacement, geodetic zenith/optical delay, spatial scaling
and simultaneity, and consistent PA orientation. The reflector convention
must remain front-face plus nh/c per leg; replacing it by (n-1)h/c changes the
reference point. The moving-body Shapiro review is a separate bounded task.

Small event-time offsets and local vector differences should be preserved;
large BCRS epoch subtraction or direct double-SPICE queries can introduce
submillimetre representation floors. Compare intermediate endpoints, each
physical delay, both solved event times, and the final atomic interval before
interpreting any residual difference. A matching model would exclude plumbing
errors within shared inputs; it would not validate the orbit or fitted station
realization independently. The completed launch audit already excludes the
proposed omitted/doubled ~1.6ms map, so a full second solver is a larger check
rather than an immediate correction justified by that hypothesis.
