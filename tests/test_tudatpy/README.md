# Python tests and optional features

Build the kernel with the chosen development environment, then run from the
repository root:

```bash
PYTHONPATH=build-python-arcs/src python -m pytest tests/test_tudatpy --no-remote-data -q
```

This collects the complete Python suite. Do not manually ignore serialization
modules or deselect the GRAIL/MRO download tests. Selection follows the loaded
binary and the requested network policy:

- `requires_serialization` tests run when the compiled
  `tudatpy.kernel._build_with_serialization` boolean is true. They are skipped
  with an explicit reason when false. A missing or non-boolean flag is a
  configuration error: rebuild and check `tudatpy.kernel.__file__`. Missing
  serialization methods never imply that the feature was intentionally disabled.
- The Cereal file-I/O module and Cereal-dependent pickle classes carry this
  marker. `Time`'s independent tuple-based pickle tests run in both builds.
- `remote_data` tests are skipped with `--no-remote-data`. This includes the
  GRAIL and MRO residual tests, which download external SPICE/tracking files.
  Tests requiring network access remain enabled by default or with
  `--remote-data`; existing credential requirements and outage handling apply.

For a focused policy check:

```bash
PYTHONPATH=build-python-arcs/src python -m pytest \
  tests/test_tudatpy/test_serialization_build_policy.py --no-remote-data -q
```

On 2026-11-12, the rebuilt serialization-OFF kernel passed the full offline
command with **694 passed, 329 skipped**, without manual exclusions. The public
short independent-arc example also passed after rebuilding. These counts include
the uncommitted Python scaling-benchmark tests in this working tree.

A separate kernel built with `TUDAT_BUILD_WITH_SERIALIZATION=ON` must run the
Cereal roundtrip tests to validate serialization itself. The policy tests check
both values of the selection flag; they do not validate an ON binary's actual
roundtrips. The subsequent separate ON build was validated on 2026-11-12:
**915 passed, 108 skipped**, plus four native serialization CTests. The targeted
serialization/policy run passed 235 tests with two explicit refactoring skips.
See [the validation record](../../docs/serialization_validation_20261112.md)
for the fixture and value-equality repairs found by enabling real propagation
result tests. Eight accidental skips for an unrelated `spice` package and a
stale expected-failure mark were removed. The two deliberate observation-file
refactoring skips remain.

Both full-suite totals above include the 28 scaling-driver tests already present
in the validation working tree, whose benchmark changes are committed separately.

These test-policy changes and this documentation were AI-assisted. Review the
feature-marker boundaries when introducing more optional serialization tests.
