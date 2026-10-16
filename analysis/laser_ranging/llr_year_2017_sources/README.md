Calendar 2017 sources were selected before annual residual evaluation: UTC
[2017-01-01, 2018-01-01), APOLLO and Grasse 7845, all five reflectors and both
reported Grasse wavelengths. The [manifest](../llr_year_2017_dataset.json)
pins 311 APOLLO Group-e rows and 1,271 Grasse CRDv1 rows. Exactly one Grasse
observation is excluded at Dominic's request: `apollo15_201703.npt:171`,
2017-03-09T18:34:15.297622200200 UTC. Thus 311 + 1,270 observations remain.
No quality, uncertainty, weather, or residual selection is applied.

The provider repeats the complete calibration/weather/NP block at
`luna21_201708.npt:100–102` and `103–105` within one session. The two NP
lines are byte-identical: L2, 2017-08-28T19:50:52.466626381600 UTC,
RTT 2.663180684343 s. Both original source-line identities remain in the
1,271 source / 1,270 retained counts. This is a provider-record repetition,
not an independently repeated measurement; no deduplication is applied.

All 60 official EDC monthly files are available. The ten March/April inputs
are byte-identical to the previous frozen two-month inputs. These products
also contain 33 Matera 7941 observations, outside the declared station scope.
The 429 Grasse sessions cover 114 UTC dates; APOLLO covers 42 dates.
[Monthly counts](monthly_coverage.csv) retain the source and retained totals.
Grasse has 93 green (532.1 nm) and 1,178 IR (1064.2 nm) source observations;
the sole requested exclusion is IR. APOLLO reports 532.0 nm throughout.

The [native reader audit](native_reader_validation.json) checks every source
product and APOLLO row without solving any observations. All Grasse headers
retain range type 2, transmit epoch event 2, UTC code 4, system 78/occupancy 1,
quality 0, troposphere/CoM/receive-amplitude/spacecraft-delay corrections 0,
and station-system-delay correction 1. Calibration is already applied.
All 1,271 Grasse weather records have measured origin 0. Bin RMS remains
metadata, not a normal-point uncertainty. The same published station profiles,
ephemeris choices, physical corrections and solver tolerances remain specified.

Three June 22 points outlast their measured weather: L1 at 09:47:06.763983741400
and L2 at 10:00:45.356444133700 and 11:19:16.531986743800 UTC. The last actual
measurement is at 09:29:16 UTC: 880.13 hPa, 295.95 K, 65% humidity. Earlier
09:14:48 weather is 880.11 hPa, 296.38 K, 64%; this does not establish later
weather stability. Official daily A15/L1/L2 products repeat the monthly weather
and supply no later bracket. Their URLs and hashes are in the manifest.
These three exact source identities are retained with predeclared bounds of
1,076, 1,895 and 6,606 s, including the existing five-second endpoint margin.
Holding the final measured state for those spans is an explicit, unvalidated
weather assumption. It does not establish intrinsic observational noise or
equal weather accuracy for every point. The ordinary 1,000 s guard remains
unchanged for all other points.

Provider bytes and the exact APOLLO subset remain in the research cache
`laser-ranging-work/research/llr-year-2017/`; large raw products are not
copied into the repository. The selected Grasse inventory preserves original
monthly filenames, line numbers, printed UTC/RTT fields, wavelengths, session
identities, weather ages and the requested exclusion marker. Recreate it with:

```sh
python analysis/laser_ranging/prepare_llr_year_2017.py \
  --data-root /path/to/llr-year-2017 \
  --apollo-release /path/to/160912e_apollo_np2.txt
```

The preparation script verifies the pinned bytes and source-row selection;
it does not calculate residuals, fit parameters, or change observations.

[CRDv1.01 §3.4](https://ilrs.gsfc.nasa.gov/docs/crd_v1.01.pdf) defines origin 0
as measured meteorology and recommends interpolation at NP epochs. It does not
specify the analysis's 1,000 s expiry or validate terminal holding. The long-hold
uncertainty therefore remains unquantified rather than absorbed into a fitted bias.
