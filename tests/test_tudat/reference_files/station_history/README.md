Minimal fixed-width numerical extracts of the hash-pinned original IGN/IERS ITRF2020 products.
The source URLs, hashes, retrieval dates and exact original line/index mapping are in
`sources.json`. Only selected station records are retained. Fixture headers have the
parameter count adjusted and estimate indices renumbered sequentially; every other
selected record field remains byte-for-byte unchanged. Omitted covariance blocks are
not represented as a measured zero covariance. These fixtures are offline.

AB07 contains six position solutions but four independent velocity-discontinuity
identifiers. Observation spans intentionally differ from the discontinuity validity
intervals. The PSD and seasonal fixtures are parameter-only products and do not need
Cartesian positions. IGN relaxation parameters report unit `m` despite representing
years in the original IGN model: that conflict is deliberately retained.

No original reference implementation is redistributed. These are small numerical
product extracts, with source identity preserved, rather than a software licence grant.
The primary format is SINEX 2.02:
https://ivscc.gsfc.nasa.gov/products-data/sinex_v202.pdf

`repeated_psd.snx` retains real AB48 EXP and ANTC LOG multiplicities at one event.
`duplicate_sites.snx` retains GALA/A with distinct GNSS/DORIS DOMES identifiers
and identical repeated 7894/A metadata; no DOMES association is guessed.

`abmf_equipment.snx` retains the original lower-case `abmf` identifier, 25
receiver intervals, three antenna intervals and three UNE eccentricities from
the pinned IGS master SINEX template. Receiver firmware and antenna radome names
remain verbatim within their standard fields. The IGS antenna azimuth extension
after column 68 remains in each record's raw provenance. Header zero epochs stay
unspecified and equipment zero endpoints retain their header-default meaning.
The IGS describes this template at https://igs.org/acc/reprocessing/ ; source
bytes and extraction line mappings are recorded in `sources.json`.

`cedu_equipment.snx` retains the template's approximate latitude printed as
`-31 51 60.0`: exactly 60.0 arcseconds carries to the next minute when converted
to radians, with the original rounded F4.1 field retained in provenance. Values
above 60 are rejected; no numerical tolerance or fitted position is changed.

`igs_malformed_site.snx` deliberately retains malformed `qs07` SITE/ID source
line 389: its station description has 21 instead of 22 columns, shifting the
coordinates left by one. The strict reader rejects this unchanged record with
source/line information. Consequently the complete pinned IGS template is not
claimed to parse successfully; the eight complete pinned IGN products do parse.

Reproduce offline with:

```
python3 extract_fixtures.py /path/to/hash-pinned/station-history-cache
```

All source checksums are verified before extraction.
