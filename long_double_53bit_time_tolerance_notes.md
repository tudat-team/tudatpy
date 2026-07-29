# Tolerance handling for 53-bit `long double` time storage

These changes are implemented in the high-precision state-scalar test. The test
also prints the detected `long double` mantissa width and the relevant range,
time-conversion, and Doppler residuals so platform-specific failures provide
enough information for further adjustment.

## Rationale

`tudat::Time` stores its fractional component in `long double`. On x86 Linux,
`long double` normally has a 64-bit significand, but on systems where it is
equivalent to binary64 it has a 53-bit significand. Since `Time` normalizes its
fractional component within a 3600-second interval, the latter can have a
spacing of approximately 0.5 ps near the top of that interval.

The implementation detects this using
`std::numeric_limits<long double>::digits >
std::numeric_limits<double>::digits`, rather than `sizeof(long double)`.

## Tolerance changes

- TDB-to-UTC and UTC-to-TDB picosecond increment agreement:
  - retain `2e-16 s` when `long double` has more precision than `double`;
  - use `7.5e-13 s` for a 53-bit `long double`;
  - additionally require both converted increments to remain strictly positive,
    so a lost picosecond separation cannot pass under the wider tolerance.

- Two-way Jupiter range reference:
  - retain `1e-11 m` with extended `long double`;
  - use `1e-9 m` with a 53-bit `long double`, because the iterated reflection
    epoch is `Time`-valued and can move by approximately 0.5 ps.
  - The one-way range tolerance remains unchanged at `1e-18 m`.

- DSN Doppler sequences:
  - retain the `1e-9 Hz` increment and trend-residual bounds with extended
    `long double`;
  - use `2e-6 Hz` for individual increments and increment residuals on a 53-bit
    `long double`;
  - use `1e-6 Hz` for trend residuals on a 53-bit `long double`.

- Aggregate 100 ns Doppler sweep:
  - retain the same-sign and 2% relative-error checks with extended
    `long double`;
  - use an absolute `1e-6 Hz` residual check on a 53-bit `long double`, because
    the expected physical change is below that platform's `Time` noise floor.

- Double-to-quad Doppler comparison:
  - retain the historical `1e5` ratio with extended `long double`;
  - omit this ratio for a 53-bit `long double`, while still requiring a non-zero
    quad response, a zero double response, and the microhertz bounds above.
  - This ratio compares a baseline double-versus-quad offset with an epoch
    response and is not a direct resolution metric.
