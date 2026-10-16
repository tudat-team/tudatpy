Independent saved-output validation passes for all three annual cases:
DE440, INPOP17a with explicit JPL glass, and INPOP17a without separate glass.
The [compact summary](annual_output_verification.json) pins each annual CSV
and its corresponding current 305-point baseline by SHA256.

Each CSV contains exactly 1,579 eligible source identities: 311 APOLLO and
1,268 Grasse. The requested March 9 exclusion and both impossible February
flight times are absent. Both August source occurrences at lines 102/105
remain. Only the three declared June 22 points carry the long-weather flag.
All computed ranges and residuals are finite; epoch events remain transmit.
The coverage exactly matches eligible source rows: APOLLO January 6–December 31
on 42 UTC dates and Grasse January 6–December 21 on 113 dates. The source
manifest's 114 Grasse dates include February 20, whose only two source records
are the impossible flight times; raw and final coverage are distinct.

All 305 March–April source identities in each annual case match its prior
current-305 control. All 915 computed full two-way ranges are bit-identical
as binary64 values. APOLLO's changed subset filename is normalized only for
identity matching; original full-release line numbers and printed UTC remain
part of every key. No residual, correction, alignment or numerical tolerance
is introduced by this verification.

The portable [validator](../verify_llr_year_outputs.py) has no native imports.
For each case, supply the saved annual CSV, corresponding current-305 CSV,
frozen manifest, research data root and a JSON output path:

```sh
python analysis/laser_ranging/verify_llr_year_outputs.py \
  --annual /path/to/annual/case/residuals.csv \
  --baseline /path/to/current_305/case/residuals.csv \
  --manifest analysis/laser_ranging/llr_year_2017_dataset.json \
  --data-root /path/to/llr-year-2017 \
  --output /path/to/verification.json
```

The annual `inpop17a_no_separate_glass` case corresponds to baseline
`inpop17a_zero_explicit_glass`; the other case names match directly. No
observation solves or native builds were performed for this audit.
