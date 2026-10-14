Exact numeric records from three potential-format products are pinned in
`manifest.json`: full-source SHA256, inclusive one-based source ranges, extract
SHA256 and edits (none). These are data records, not external program source.
The legacy FES2004 source is the user-requested local ILR data file; its complete
SHA256 matches the coordinator's source inventory. No old ILR program is built.

Each extract retains every source constituent in provider order and all rows
through degree 3. EOT11a additionally retains the full degree/order-120 M2 field.
Named-profile extents come from the pinned complete products, not the observed
extract maximum. Thus a request inside a published wave extent must find every
selected term; outside that extent the named product's contribution is zero.
The intentional single C20 Om1/Om2 rows are preserved. The raw layer retains
missing-entry identity and published degree-zero/one values. Body-centered
synthesis currently selects degrees >=2.

FES2004 units are 1e-12 normalized Stokes coefficients; FES2014b-v1 and EOT11a
units are 1e-11. EOT11a's leftover 'in cm' text is contradictory and rejects
unless its identified COST-G potential profile is explicitly selected.
FES2014b-v1 retains radiational S1=164.555 and Sa=056.555. Six-character source
IDs are zero-padded with correct digit offsets. Source GM/R are absent and
remain mandatory consuming declarations; the reader never guesses them.

Equation 6.15 gives Ccos=Cp+Cm, Csin=Sp+Sm, Scos=Sp-Sm, Ssin=Cm-Cp.
Physical Sn0 is zero. Potential products already contain their published phase
and loading conversion, and prograde zonal amplitudes require no extra doubling.
The old ILR phase/matrix accumulation structure is reused with corrected S signs,
Doodson digit offsets, checked EOF/duplicates and degree/order dimensions.

Four exact COST-G M2 force/argument samples reuse the ocean-pole fixtures'
ITRF positions and ICRF-to-ITRF rotations. The supplied benchmark Doodson
arguments are injected unchanged, so its documented UTC/UT1 distinction is not
recomputed by an unrelated production clock. The explicitly labelled matched
benchmark convention uses GM=3.986004415e14 m^3/s^2 and R=6378136.460 m from its
EIGEN-6C4 baseline header. These constants are not asserted as provider EOT11a
metadata. The largest force residual is 1.57e-13 m/s^2 against the fixed
1e-11 m/s^2 target. No claim of identical unrounded coefficient conversion.

Admittance takes an exact ordered major inventory, target inventory and supplied
matrix. Targets replace the input inventory and contribute once; no extra nodal
modulation is applied. This general evaluator does not silently select a
published interpolation profile or infer a matrix's major-wave ordering.

Offline reproduction, using the shared extraction script:

```
python ../ocean_pole_gravity/regenerate_provider_extracts.py /path/to/references/gravity .
```

The cache includes `legacy-fes2004-potential.txt`, an unchanged local copy of
ILR's `DataFiles/GravityFieldVariations/fes2004_Cnm-Snm.txt`. All complete source
and extract hashes are verified before any writes. Source URLs are in the
manifest; provider explanations are pinned separately in the coordinator cache.
