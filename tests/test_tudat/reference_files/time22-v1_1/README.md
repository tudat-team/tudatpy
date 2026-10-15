TiME22 revision 1.1, degree-two extracts
======================================

Source: Sulzbach, Roman; Balidakis, Kyriakos; Dobslaw, Henryk; Thomas, Maik
(2022), *TiME22: Periodic Disturbances of the Terrestrial Gravity Potential
Induced by Oceanic and Atmospheric Tides*, GFZ Data Services,
https://doi.org/10.5880/GFZ.1.3.2022.006. Revision 1.1 is dated March 22, 2025.
The provider licenses these data under Creative Commons Attribution 4.0
International, https://creativecommons.org/licenses/by/4.0/.

`ATM/` and `OCN/` retain all 16 and 57 paired constituents, respectively, but
only original GFC records of degree zero through two. These are sparse test
extracts, not complete gravity fields. Original physical/product header values
(including maximum source degree 360) and selected explanatory declarations
are unchanged. `manifest.json` identifies every retained source line, full
source SHA256, extract SHA256 and the exact official archive URL/hash. All six
supplementary tables are copied unchanged, including the full 57-by-435 ocean
admittance and 16-by-16 atmospheric identity matrix. No coefficients were
rounded or edited. Degree zero/one are preserved only to check the raw record;
Earth consumption starts at degree two.

Reproduce offline from the official archive extracted to `/path/to/TiME22`:

```
python tests/test_tudat/reference_files/time22-v1_1/regenerate.py /path/to/TiME22
```

The script verifies the pinned full-file hashes, rewrites the extracts and the
small compiled ordered profile inventories, and freezes an independent
80-digit Decimal synthesis. It imports no Tudat or numerical Python package.
Reference phases use the six binary64 arguments in the existing unmodified
COST-G `harmonic_potential/costg_90doodsonArguments.txt` fixture. References
follow the provider's contraction order: multiply the input-by-target matrix
by target cosines/sines, then contract these weights with input Stokes fields.
The C++ evaluator instead sums the already mapped canonical target fields.
Both optional-admittance settings are covered at all four supplied epochs.
The atmospheric S1 is `164.556`; the ocean radiational S1 is `164.555`.
Loading and Doodson-Warburg phase corrections are included by the provider.
No additional loading, nodal modulation, phase shift or zonal doubling applies.

The per-coefficient absolute comparison budget is derived from the actual
summation magnitudes, not a fitted test tolerance. With conservative binary64
`epsilon = 2^-52` and `gamma(k) = k*epsilon/(1-k*epsilon)`, each target contributes

```
B_t * [gamma(12)*sum_i(abs(q_i*argument_i)) + 4*epsilon
       + gamma(2*N+2) + gamma(2*M+4)]
```

where `B_t` is the absolute mapped cosine-plus-sine source magnitude, `N` is
the input count and `M` the target count. This includes phase products/additions,
input parsing, admittance multiply/adds, two trigonometric evaluations and
canonical synthesis accumulation. The Decimal reference error is negligible
at these budgets. Exactly zero zonal sine rows retain an exact zero budget.
