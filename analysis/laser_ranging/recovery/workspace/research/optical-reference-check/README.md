# Independent optical zenith-delay control

The unmodified `FCUL_ZD_HPA.F` was retrieved on 2026-10-15 from
https://iers-conventions.obspm.fr/content/chapter9/software/FCUL_ZD_HPA.F .
SHA-256: `92731affca053aad15a44be7db58dbf6df689e75cf2e1f3b39cb4d99a4da198b`.
Its original IERS Conventions Software License notice is retained intact in the
source. This reference copy and the original Tudat driver are neither distributed
by nor endorsed by the IERS Conventions Center.

`check_fcul_zenith.f` calls that routine with the exact inputs printed in its
header: latitude 30.67166667 degrees, ellipsoidal height 2010.344 m, pressure
798.4188 hPa, water-vapour pressure 14.322 hPa, wavelength 0.532 micrometres.
It does not call Tudat or any comparison package.

Reproduce in this directory with a Fortran compiler supporting the original
fixed-form source, preserving default-real literal semantics:

```sh
gfortran -std=legacy -ffixed-form -O0 FCUL_ZD_HPA.F check_fcul_zenith.f -o check_fcul_zenith
./check_fcul_zenith
```

The independent execution used GNU Fortran (conda-forge gcc 15.2.0-19) 15.2.0
with the flags above and no default-real-width flags. Compilation and execution
both exited successfully; source and driver hashes were checked unchanged. Exact
command, compiler version, hashes and output are retained in the accompanying
provenance files. The outputs are total `1.93522972487606770`, hydrostatic
`1.93299597223469832`, and wet `0.00223375264136944527` metres.

## Inconsistent printed header control

The original header prints total `1.935225924846803114`, hydrostatic
`1.932992176591644462`, and wet `0.002233748255158703871` metres. At the stated
inputs, an independent direct evaluation of the published equations and the
Tudat implementation both give a total about 3.800121 micrometres larger. This
is not a Pa/hPa conversion, wavelength conversion, or floating-point tolerance
issue. The hydrostatic header value implies a height near 2003.344 m, seven
metres below the stated input; that is an inference, not an independently
published correction to the input. The driver therefore retains 2010.344 m.
The source's few single-precision literals account only for picometre/sub-nm
differences and cannot explain the micrometre discrepancy.

The test must retain its 1 nm arithmetic tolerance. Its independent reference
is the unmodified official routine evaluated at the printed inputs, not
the Tudat implementation's output. The inconsistent header values remain
recorded here instead of being silently relabelled as consistent controls.

## Height and coefficient convention

The official chapter https://iers-conventions.obspm.fr/content/chapter9/icc9.pdf
defines geodetic height in both Eq. 9.4 (zenith) and Eq. 9.10 (FCULa mapping).
Table 9.1 gives the temperature coefficient `a21 = 0.2346e-5`, as used here.
`FCUL_A.F` calls its height mean-sea-level height, which conflicts with that
chapter's Eq. 9.10 wording. Tudat follows the chapter and supplies the station's
ellipsoidal geodetic height to both functions. No geoid conversion is inferred.

For scale only, at latitude 32.78035394 degrees, height 2786.62 m, temperature
280.65 K, a +100 m change in the mapping height alone changes a 2 m zenith delay's
slant delay by 0.07736 mm at 20 degrees elevation, 0.02157 mm at 30 degrees,
0.00514 mm at 45 degrees, and zero at zenith. At 10 degrees the change is
0.57570 mm, and at 3.1 degrees it is 8.33529 mm. These are sensitivities of the
published continued fraction, not the actual APOLLO elevations or a geoid model;
the latter must be logged during the first observation run. Low-elevation
atmospheric gradients remain a separate physical limitation.
