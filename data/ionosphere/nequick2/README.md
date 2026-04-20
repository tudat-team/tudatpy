# NeQuick-2 Data Files

CCIR coefficient files and MODIP grid for the NeQuick-2 ionospheric model.

## Source

ITU-R Recommendation P.531-15, NeQuick2 P.531-15 reference implementation (release 2.1.0, October 2011).

Developed by the Aeronomy and Radiopropagation Laboratory, Abdus Salam International Centre for
Theoretical Physics (ICTP), Trieste, Italy.

## Files

- `ccir11.asc` ... `ccir22.asc` — Monthly CCIR coefficients for foF2 and M(3000)F2 (January through December)
- `modip.asc` — Modified dip latitude grid (1° lat × 2° lon)

## Installation

These files must be placed at `~/.tudat/resource/ionosphere/nequick2/` for runtime access by the
NeQuick-2 ionospheric model. This is done automatically when installing the tudat-resources package,
or can be done manually:

```bash
mkdir -p ~/.tudat/resource/ionosphere/nequick2
cp data/ionosphere/nequick2/*.asc ~/.tudat/resource/ionosphere/nequick2/
```

Alternatively, pass the `ccir_data_path` parameter to `nequick2_ionospheric_light_time_correction()`
to specify a custom path.
