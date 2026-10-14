# Optional RP03 loading-origin translation

IERS2010 Eq18 specifies X_CE=X_CM-dX. Eq20/Table7.6 supplies Cartesian XYZ
harmonic coefficients in metres, using the same midnight UT1 phase as Eq19.
Therefore a CE station loading vector converted toCM is deltaCM=deltaCE+dXYZ.
Rotate local source displacement toCartesian first, add dXYZ once, then project
into any requested station topocentric frame. The table itself supplies no grid
lookup latitude convention or local displacement basis.

The narrow pure helper computes only the published table. It does not convert
origins or identify a loading product. Its seven frozen references come from a
new independent Fortran implementation of the published equations using complex
Fourier factors and cross-product Cartesian projection; this is not original IERS
reference software, a provider CM grid, or a claim of complete model conformance.
Quarter-day C++ cases separately freeze exact arithmetic on published coefficients.

Proposed full-model guard (default preserves declared source origin): explicit
source profile iers2010Rp03NoIbFarrell, matching model/Green metadata, and known
source origin are required for opt-in CE-toCM. CE adds the table once; CM retains
the supplied vector. Unspecified origin/profile, incompatible model or Vienna
product rejects conversion. Source profile is declared/identified by the caller,
not inferred from filename or latitude. The untouched neutral provider records
retain original origins and headers. A converted model result declares CM.

Full model tests should cover preserveCE, opt-inCE->CM, alreadyCM unchanged,
unknown/profile/model rejection, equality between Cartesian addition and local
projection, plus independent grid lookup latitude and displacement basis choices.
Continuous UT1 comes from the shared EarthDeformationEpoch; velocity requires
its dUT1/dTDB factor and must not silently substitute a zero rate.

Primary source: https://iers-conventions.obspm.fr/content/chapter7/icc7.pdf
Provenance, source PDF/reference-generator hashes and compiler are pinned in
origin_translation_provenance.json. Neither the unresolved Bernese station
comparison nor another Vienna CM model validates this conversion's applicability.
