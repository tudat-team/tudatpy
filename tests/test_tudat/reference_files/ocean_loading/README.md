# Independent HARDISP fixtures

These fixtures come from the corrected IERS HARDISP Fortran package (2015 numerical
corrections and 2016 leap-second update), not from the production C++ expressions.
The twelve unmodified upstream files retain their original notices. Their SHA256
checksums, URLs, station sources and compiler provenance are pinned in
`reference_provenance.json`. The generator verifies every upstream checksum.

From this directory, with GNU Fortran 15.2 available:

```sh
python3 generate_hardisp_fixtures.py --compiler gfortran
python3 generate_hardisp_fixtures.py --compiler gfortran --check
```

The generator renames only the executable main program and increases its printed
precision; numerical reference routines remain unchanged. `--check` runs the same
independent reference and compares the frozen header without changing any files.
The generated header has clang-format guards so formatting preserves its bytes.
Fortran is optional tooling for regeneration; production and C++ tests need none.
The 36 cases cover six stations, nonmidnight epochs and 25 years. The one micrometre
C++ tolerance accounts for the reference's single-precision synthesis.

Pinned upstream SHA256 checksums:

- `ADMINT.F`: `478e3f0c001c09dd4d9e2e9920033d29e76eb60b808238e416e6db9ffaeae6c8`
- `ETUTC.F`: `a067c10de2c63269a54b6a14dce2daa49c30bfd52e66086b76492cd438414e1f`
- `EVAL.F`: `b4056a09ec5d77674cab8eac0fa32b5d1be8f471fbe33e753e81ab97a7767add`
- `HARDISP.F`: `df8556e0db396fe5358a1bb2edd8ae2ea19253993a5f8a8641495b9ca0224f75`
- `JULDAT.F`: `d1f39f83503178711845532e1a4ba46d0824b896fa961d81c983d694d48dde73`
- `LEAP.F`: `31d18153242823606beb9690ab0c685dd5403fe1e9bb5f214e22d21dd5e6a771`
- `MDAY.F`: `d4eae8a9a0b866a63f22134fb57bad1c57fc1f3baf98450892ddeee4f7ee8aec`
- `RECURS.F`: `8a3c88d69cebd130981887dc9c8c2f9e3d2af26fded3e052385c022126f9d44b`
- `SHELLS.F`: `6c69dede79f9adfb16d38dcbd890bbbfc3e4e47b2cbacd5fabfb520deba86701`
- `SPLINE.F`: `a13cc2b405079a2b86872341ca50d9b3f0cc219e43d1ce301b461af15df857f0`
- `TDFRPH.F`: `24068e0cd1e2e210fab7dd8d4473bb60ff8335bfcb53a2fc12dad7e6d1cccd19`
- `TOYMD.F`: `9b69b6a27d544215c516da60351e4ec51d03b04ce9bf71b6608ebf55080bf70a`
