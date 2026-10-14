# Static ICGEM fixtures

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

`sources.json` pins URLs, whole-source hashes, exact selections and extract hashes.
Reproduce offline from the coordinator's cache:

```sh
python extract_fixtures.py /path/to/references/gravity
```

The provider snapshot is attributed to Dahle et al. (2019), *GRACE-FO Geopotential
GSM Coefficients GFZ RL06*, version 6.1, DOI `10.5880/GFZ.GRACEFO_06_GSM`.
