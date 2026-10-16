Published Hannover station sensitivity
======================================

``../../hofmann2017_station_profiles.json`` records the independently preselected
**standard** Hofmann 2017 solution, with source hashes and internal component
uncertainties. It imports no fitted range bias. The source fits station and
reflector coordinates together with its own lunar orbit/rotation; this is a
station-only sensitivity with DE440 retained, not a matched DE440 reduction.
Tables 5.4--5.5 and A.3 are used consistently. The separately constrained-velocity
solution in Tables 5.6/A.5 is excluded.

Both sites use telescope-axis intersection (SRP), J2000.0 and conventional
tide-free coordinates. APOLLO is explicitly not an official ITRF2008 station
solution: the a priori position came from Murphy (private communication, 2008),
with a nearby White Sands GPS velocity. Matera's standard velocity is fixed to
ITRF2008, while its position is estimated. Comparing the fitted network to ITRF
does not establish an exact frame tie; no Helmert transform is inferred.

Coordinates and velocities are in the JSON, copied at published precision.
TT Julian-year propagation is an explicit diagnostic interpretation of the
source's J2000.0 and metres/year labels. The station time scale is not explicitly
restated there. Full covariance, later local ties and post-2015 station history
are unavailable; the quoted internal uncertainties do not establish a complete
prediction uncertainty for these later observations.

Coordinate-only comparison
-------------------------

``station_position_differences.csv`` retains the actual eleven APOLLO tags and
two Matera tags. Values are Hofmann minus the existing undeformed station
position, with no frame adjustment, loading, rotation, light-time calculation
or residual selection. At the first tag of each case, in millimetres:

=========================  =========  =========  =========  =========
Reference                  Delta X    Delta Y    Delta Z    Norm
=========================  =========  =========  =========  =========
APOLLO survey2005           -189.388   +132.729   +262.311    349.703
APOLLO INPOP17a G+IR         +52.432    -13.718    -36.977     65.609
Matera SLRF2020               +5.649    +75.642    -47.824     89.670
=========================  =========  =========  =========  =========

The vector variation within either frozen arc is below one micrometre. These
norms are not range corrections. APOLLO's vertical velocity differs from INPOP
by -5.3 mm/year; the reference positions partly compensate its accumulated
effect. Thus a centimetre-scale sensitivity is plausible, without identifying
either published solution as correct for the frozen DE440 model.

Reproduction uses ``x0 + v * ((UTC - J2000 noon) + 69.184 s) / 31557600``
for both Hofmann and INPOP at these dates. The survey remains constant. Matera's
existing native linear history is evaluated algebraically as
``x2015 + vSI * (nominal UTC - 2015-01-01 UTC)``, excluding leap seconds from
that interval as documented by its adapter. Original CSV tag strings are
preserved; Python's microsecond parsing truncation is negligible for this
coordinate calculation. The existing profiles/adapter are in commit
``bb9d38085cd4a192c318d367d25e2aeaafe19f29``; their original frozen summaries
provide the exact survey and SLRF2020 reference vectors.
The latter is solution ``7941 A 1`` in the pinned SLRF2020 SINEX identified by
``matera_dataset.json``; its velocity is
``(-0.0186283100020115, 0.0190569144613420, 0.0146220896897415) m/year``.

Convention bounds and lunar limits
---------------------------------

IERS 2010 Eq. 7.14 gives conventional tide-free to mean-tide displacements
of 24.16 mm at APOLLO and 29.52 mm at Matera (geocentric latitude; radial/north
components +7.77/-22.87 mm and -15.87/-24.90 mm). The source and frozen Matera
profile both specify tide-free coordinates. An extra permanent-tide correction
is unsupported and this convention alone cannot explain 8--15 cm.

The Earth spatial scale ``LC * radius`` is about 94.4 mm; it is already applied
in the frozen model. Hofmann Eq. 2.13 includes this same scale, solar potential
and Lorentz map. ``LG * radius`` is only 4.44 mm. The older source's lunar LB
versus DE440 LM difference is about 1.18 mm at lunar radius; it does not justify
changing the matched DE440 map. Matera's documented optical reference point
already has zero eccentricity. No additional instrument term was recovered.

Hofmann section 5.5.1 explicitly pairs its Table 5.13 reflectors with station
Table 5.4 and its own integrated principal-axis system. It warns of correlations
with orbit, gravity, core/libration and Earth-Moon mass parameters. Its reflector
coordinates differ from DE430 by up to 0.6 m per component; a fitted network
transform reduces this to centimetres. This demonstrates solution dependence,
not a transferable reflector correction. Neither its PA/ME transform nor its
reflector coordinates should be mixed with the matched DE440 PA kernel.

Singh 2023 updates LUNAR using DE440 initial information, but continues joint
ephemeris/parameter estimation. Table 3.2 supplies station uncertainties, not
absolute coordinates; Table 4.1 supplies differences and explicitly omits
APOLLO because no ITRF2020 coordinates exist. Appendix A lists fitted parameter
types, not numerical station vectors. No recoverable matched DE440 absolute
station solution was found in this bounded literature review.

Primary sources
---------------

* `Hofmann 2017, DGK C797 <https://dgk.badw.de/fileadmin/user_upload/Files/DGK/docs/c-797.pdf>`_:
  sections 2.2.3, 5.3, 5.5.1; Tables 5.4--5.5/A.1/A.3. Hash in profile JSON.
* `Biskupek 2015, DGK C742 <https://dgk.badw.de/fileadmin/user_upload/Files/DGK/docs/c-742.pdf>`_:
  sections 4.2.2--4.2.3, printed pp.34--36. Hash in profile JSON.
* `Singh 2023, DGK C903 <https://dgk.badw.de/fileadmin/user_upload/Files/DGK/docs/c-903.pdf>`_:
  printed pp.24--25, 41--42, 89--90; SHA256
  ``4570bfd1701c2407d4da1958c7d312bb5ab77dcf9f828ccd969423e901a47f7b``.
* `IERS Conventions 2010 <https://iers-conventions.obspm.fr/conventions/content/tn36.pdf>`_:
  section 7.1.1.2, Eq. 7.14; cached PDF SHA256
  ``fdb5b74c5a6135c9d20fd99191741e8b62ac8d712eb96f046e7eeb8fb4236f51``.
  The three thesis PDFs/texts and download hash manifest are retained externally
  under ``research/hannover-llr``.
