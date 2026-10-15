# ICGEM fixtures

`gfz_rl061_ddk4_202104_degree2.gfc` retains the complete original header and the
degree 0–2 coefficient records of GFZ's DDK4-filtered April 2021 GRACE-FO RL06
v6.1 snapshot. The original degree-60 declaration, error columns, preamble and
CC BY 4.0 licence/citation remain unchanged. The preamble includes a
`max_degree` line before `begin_of_head`; it is a comment, not a duplicate header.
This extract exercises static `gfc` syntax, not time-variable `gfct` synthesis.
Use an explicit degree-two truncation when constructing a field from the extract.

`costg_eigen6_c4_duplicate_c00.gfc` preserves the first 16 lines of the COST-G
benchmark input, including its two numerically identical C00 records. It is a
negative fixture: duplicate records reject without silently repairing the source.
The extract contains metadata facts and numerical coefficient records.

`eigen6s4_v2_c20_1950_2006.gfc` preserves the full original EIGEN-6S4 Version 2
header, static C00 and verbatim C20 records through the interval ending in 2006.
The original `max_degree 300` declaration remains unchanged. This sparse extract
does not manufacture omitted coefficients. Its annual/semiannual records for
1950--2003 have independent reference/validity dates spanning many adjacent
`gfct` pieces. The later printed `.0060` suffix is retained without interpreting
it as a time of day or silently rolling it over. The reader validates Gregorian
dates and fixed-width token ordering only; a consuming adapter must explicitly
choose the time scale, subday interpretation and definition of a year.
The raw reader retains these temporal records; the static settings adapter rejects
them explicitly, including when only a static degree-zero truncation is requested.

The EIGEN extract is adapted by selection only from Foerste et al. (2016),
*EIGEN-6S4 A time-variable satellite-only gravity field model to d/o 300 based on
LAGEOS, GRACE and GOCE data from the collaboration of GFZ Potsdam and GRGS
Toulouse*, Version 2. GFZ Data Services, DOI `10.5880/icgem.2016.008`.
Its DataCite rights metadata declares CC BY 4.0:
https://creativecommons.org/licenses/by/4.0/legalcode. Exact source ranges,
source/fixture hashes and the rights-metadata URL are recorded in `sources.json`.

`sources.json` pins URLs, whole-source hashes, exact selections and extract hashes.
Reproduce offline from the coordinator's cache:

```sh
python extract_fixtures.py /path/to/references/gravity
```

The provider snapshot is attributed to Dahle et al. (2019), *GRACE-FO Geopotential
GSM Coefficients GFZ RL06*, version 6.1, DOI `10.5880/GFZ.GRACEFO_06_GSM`.
