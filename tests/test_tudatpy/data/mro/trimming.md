Interval in test:
- Start epoch (TDB): 385256700.0 -> 2012-03-17T11:43:53.814 UTC
- End epoch (TDB): 385260300.0 -> 2012-03-17T12:43:53.814 UTC

Interval with margin (+-1h):
- Start epoch (UTC): 2012-03-17T10:43:53.814
- End epoch (UTC): 2012-03-17T13:43:53.814

Trim commands:

SPK (see mro_psp22-trimmed.log for contents of `spkmerge-mro.txt`):
```
spkmerge spkmerge-mro.txt
```

CK:
```
ckslicer -LSK /Users/lhinuber/.tudat/resource/spice_kernels/naif0012.tls -SCLK mro_sclkscet_00112_65536.tsc -INPUTCK mro_hga_psp_120313_120319.bc -OUTPUTCK mro_hga_psp_120313_120319-trimmed.bc -START 2012 MAR 17 10:00:00.000 -STOP 2012 MAR 17 14:00:00.000
ckslicer -LSK /Users/lhinuber/.tudat/resource/spice_kernels/naif0012.tls -SCLK mro_sclkscet_00112_65536.tsc -INPUTCK mro_sc_psp_120313_120319.bc -OUTPUTCK mro_sc_psp_120313_120319-trimmed.bc -START 2012 MAR 17 10:00:00.000 -STOP 2012 MAR 17 14:00:00.000
```
