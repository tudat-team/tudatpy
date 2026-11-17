Interval in test:
- Start epoch (TDB): 387314833.0 -> 2012-04-10T07:26:06.814 UTC
- End epoch (TDB): 387318433.0 -> 2012-04-10T08:26:06.814 UTC

Interval with margin (+-1h):
- Start epoch (UTC): 2012-04-10T06:26:06.814
- End epoch (UTC): 2012-04-10T09:26:06.814

Trim commands:

SPK (see grail_120301_120529_sci_v02-trimmed.log for contents of `spkmerge-grail.txt`):
```
spkmerge spkmerge-grail.txt
```

CK:
```
ckslicer -LSK /Users/lhinuber/.tudat/resource/spice_kernels/naif0012.tls -SCLK gra_sclkscet_00014.tsc -INPUTCK gra_rec_120409_120415.bc -OUTPUTCK gra_rec_120409_120415-trimmed.bc -START 2012 APR 10 06:00:00.000 -STOP 2012 APR 10 10:00:00.000
```
