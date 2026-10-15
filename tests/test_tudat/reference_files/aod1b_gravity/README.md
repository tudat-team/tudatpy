This offline fixture is a **labelled reduction**, not an unchanged GFZ daily product.
It retains degree 0--2 and all four components at 00:00, 03:00 and 06:00 GPS on
2026-01-01. Every coefficient line is copied verbatim; `sources.json` maps its
fixture line to the full source line. `original_provider_header.txt` is the full,
unmodified original header.

The reduction rewrites only `MAXIMUM DEGREE`, `NUMBER OF DATA RECORDS`,
`NUMBER OF DATA SETS`, `TIME LAST OBS(SEC PAST EPOCH)`, `FILENAME` and `FILESIZE`.
They describe the extract honestly. The provider constants, software revision,
normalization, creation timestamps and format declarations remain unchanged;
creation timestamps describe the provider source, not this reduction.

Source: `AOD1B_2026-01-01_X_07.asc.gz`, 8,860,900 bytes, SHA256
`bdb32469e4f9e9408886ebe2d4ccef9c9f56e90d27a8cd56a9ad5c6ba35d1b37`.
The source URL and output hash are recorded in `sources.json`. Reproduce from
this directory with:

```
python extract_fixture.py /path/to/AOD1B_2026-01-01_X_07.asc.gz
```

Attribution: Shihora, Linus; Balidakis, Kyriakos; Dill, Robert; Dobslaw, Henryk
(2022): *Atmosphere and Ocean Non-Tidal Dealiasing Level-1B (AOD1B) Product RL07*.
GFZ Data Services. https://doi.org/10.5880/GFZ.1.3.2022.003.
The dataset's DataCite rights record declares Creative Commons Attribution 4.0
International (CC BY 4.0): https://creativecommons.org/licenses/by/4.0/legalcode.
The reduction and its changed header fields are the modifications made here.
Rights metadata was checked on 2026-10-15 through
https://api.datacite.org/dois/10.5880/GFZ.1.3.2022.003.

Scientific conventions follow the GFZ Product Description Document, revision 7.0
(2023-10-13), especially equations 2.20--2.21 and sections 6, 8.1 and Appendix A.
These are already-converted potential anomalies relative to the 2007--2014 mean;
no extra loading multiplier is applied. GPS equals TAI minus 19 seconds. GLO is
ATM plus OCN, within printed coefficient rounding; OBA has a distinct domain.
Runtime consumption begins at degree two, leaving the published degree-zero/one
records available for inspection without assuming their origin semantics.
