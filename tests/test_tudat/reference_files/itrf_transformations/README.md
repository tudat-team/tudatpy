# Official ITRF transformation input and independent equation fixtures

The complete, small original IGN 2014/2020 transformation tables are retained
unchanged with their source notices and URL/SHA256 sidecars. Identical raw text is
pinned in the independent input target and parsed through its strict reader.
The 2020 table explicitly supersedes earlier published values; its catalogue is
kept distinct from the 2014 table used by compatibility functions.

The table prescribes a **linearized** position-vector matrix: identity plus
scale on its diagonal, with the published signed rotation cross-product terms.
It does not prescribe a finite-angle rotation. Each parameter propagates by
its published rate times the decimal-year difference. SI-normalized rates must
be multiplied by 31557600 to reconstruct that propagation; physical derivatives
add the explicit dDecimalYear/dInputSecond chain. No calendar/TT/TDB mapping is
inferred from the year printed in the table.

`generate_itrf_reference.py` independently evaluates these equations at 80-digit
Decimal precision. It uses scalar Gaussian elimination for the inverse, and
propagates original per-year parameters directly, independently of production
SI normalization/Eigen expressions. It covers every published row at the table
epoch and +/-10.25 years: 75 forward/inverse moving-state fixtures. The explicit
physical year rate is 0.75/31557600, exposing an erroneous assumption of unit
Julian-year time chaining. This is an independent implementation of published
equations, not original IGN/IERS reference software.

```sh
python3 generate_itrf_reference.py
python3 generate_itrf_reference.py --check
```

The check verifies both original table hashes and the embedded production raw
text, then compares the frozen generated header. Clang-format guards preserve
fixture reproducibility. Native tolerance is 2 nm for published-table positions
and 2e-17 m/input-second for velocities. The separate artificial large-rate
round-trip uses a bound derived from machine epsilon, matrix conditioning and
the magnitudes of the cancelled terms; its residual does not set the tolerance.

Legacy functions retain their explicitly documented uniform-Julian-year seconds
since J2000 (year=2000+seconds/31557600), with the 2014 table. Native structured
epochs/model clocks make a chosen physical-to-decimal-year mapping explicit.
Source release metadata is distinct from realization identity; aligned catalogue
updates do not acquire invented datum transformations. Composition uses the
catalogue's declared reference hub, with full affine values and derivatives.
Equal nonempty realization labels give identity before catalogue lookup, even
for unknown identical labels; identity is not a catalogue-membership test.

New C++/reference code remains under the repository licence. No complete upstream
software was translated. Production/tests require no network or Python runtime;
Python is optional for regenerating these independent fixtures.
