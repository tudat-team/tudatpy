# PR 905 review notes

PR: <https://github.com/tudat-team/tudatpy/pull/905>

Reviewed integration branch: `feature/observation_redesign_refactor`

Latest reviewed/pushed commit: `29ea922b3d701a39b1b282e48b6eb8ae707163e9`

Latest full CI run: <https://github.com/tudat-team/tudatpy/actions/runs/33860641446>

## Summary

The observation-data refactor had one immediate cross-platform CI problem and a number of correctness and compatibility problems in the new `ObservationDataset` implementation and its legacy `SingleObservationSet`/`ObservationCollection` facades. The corrections have been split into issue-oriented commits on the PR branch. Small C++ or Python regression tests were added where practical.

No local build was performed. Validation was delegated to the repository's full CI matrix. Linux, macOS, and Windows all completed successfully for the final pushed commit.

## Immediate CI failure

### Mixed Boost.Test linkage modes

Four test executables included the header-only Boost.Test entry point while also being linked through the project's compiled Boost.Test setup/PCH. The observable behavior differed by platform: the Linux observation-dataset test aborted after running, while four Windows test targets failed.

The affected tests now consistently use dynamic Boost.Test linkage:

```cpp
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
```

Commit: `1005ecb8e` (`Fix Boost.Test runtime linkage`)

## Correctness findings and corrections

### 1. Dataset invariants were incompletely enforced

The dataset could accept observation matrices whose row count did not match the declared observable type. Dependent-variable values could also be installed with incompatible link/type metadata, incorrect row widths or inconsistent per-observation counts. Appending rows could mix observations with and without dependent-variable data, and assigning an empty dependent-variable vector did not reliably clear existing values.

Corrections:

- validate observation dimensions against the declared observable type;
- validate dependent-variable observable type and link compatibility;
- validate every dependent-variable row and the total row count;
- reject mixed dependent-variable presence during append;
- make an empty dependent-variable assignment clear existing values;
- avoid returning fabricated zero matrices when only a dependent-variable layout exists;
- remove duplicate observations globally rather than only adjacent duplicates;
- use stable time sorting;
- guard null selection settings in dataset accessors.

Focused C++ regression tests cover the added invariants.

Commit: `5ba1836eb` (`Enforce observation dataset invariants`)

### 2. Residual and dependent-variable rows could be paired with the wrong observations

Simulation can sort observations internally. Residual code previously subtracted the simulated rows from the observed matrix by row index, even when the source dataset was intentionally unsorted. This silently associated residuals and dependent variables with the wrong observations. Duplicate epochs also made a simple time-to-index map insufficient. The legacy `ObservationCollection` path had a related ordering problem, and the Python wrapper's omitted `sort_observations=False` argument still went through a legacy path that sorted.

Corrections:

- map simulated rows back to source rows by epoch using a queue per epoch, preserving duplicate-epoch multiplicity and order;
- verify simulated and observed metadata and counts before residual calculation;
- use ordered flattening for legacy collection residuals;
- preserve the historical explicit sorting behavior only in tracking-data ingestion;
- make the Python default `sort_observations=False` use the native unsorted path.

C++ and Python regression tests cover unsorted observations and the legacy compatibility behavior.

Commit: `9f1856dac` (`Preserve observation ordering in residual calculations`)

### 3. Copies shared mutable metadata and conversion exposed unrelated observations

Dataset copies/reductions and the legacy facades could share mutable dependent-variable bookkeeping or ancillary settings. Mutating one result could therefore change another dataset unexpectedly. Self-copying into a dataset could invalidate the source while it was being traversed. Converting one `SingleObservationSet` back to a dataset could return its entire backing dataset instead of only the represented set.

Corrections:

- snapshot self-copy sources before insertion;
- clone mutable dependent-variable bookkeeping and ancillary settings for copies/reductions;
- apply copy-on-write before dependent-variable metadata mutation;
- make `createObservationDataset(SingleObservationSet)` return only the wrapped set.

Regression tests cover metadata isolation, self-copy, conversion, and datasets containing additional blocks.

Commit: `f8ccb4635` (`Isolate copied observation dataset metadata`)

### 4. Legacy transforms lost non-diagonal and set-level weights

Moving, splitting, filtering, or copying data through the legacy facade reconstructed only simple per-row weights. This discarded correlated weight entries and could also lose set-level or advanced weighting state. A move to the same source and target metadata was not rejected explicitly.

Corrections:

- preserve the full effective weight submatrix and row active/rejection state when moving observations;
- implement legacy split/filter/copy operations through dataset reductions;
- reject identical move endpoints and incompatible source/target metadata.

Regression tests cover set-block splitting, collection copying, and correlated-weight moves.

Commit: `63c95c62a` (`Preserve full weights in legacy observation transforms`)

### 5. Facade caches could become stale and viewers could dangle

Changing link or dependent-variable metadata did not always advance the structural version used by legacy facade caches. Cached collection getters could consequently return stale groupings. `ObservationDatasetViewer` held a raw dataset pointer with no way to detect destruction or replacement of the referenced dataset, allowing undefined behavior after stack destruction or assignment.

Corrections:

- bump the dataset structural version for link and dependent-variable metadata resets;
- refresh cached collection getters when the structural version changes;
- add a lifetime token checked by raw-pointer viewers so expired/replaced datasets produce a controlled exception;
- avoid applying condition-based metadata changes to empty sets.

Regression tests cover viewer use after destruction and assignment, cache refresh, and empty-set handling.

Commit: `85d178554` (`Keep observation facades synchronized and lifetime-safe`)

### 6. Public API validation was asymmetric or unsafe

Per-observation and set-level weight setters accepted different invalid shapes. Null dependent-variable condition settings could be dereferenced. The legacy facade could return fake dependent-variable zero matrices or dereference missing bookkeeping. A misspelling was exposed in a Python keyword.

Corrections:

- apply symmetric validation to per-observation and set-level weight matrices/settings;
- reject null dependent-variable condition settings cleanly;
- stop fabricating dependent-variable zeros and guard missing bookkeeping in `SingleObservationSet`;
- delegate legacy setters through the validated dataset implementation;
- correct `ancillary_settings_per_observatble` to `ancillary_settings_per_observable`.

Focused C++ tests cover invalid inputs.

Commit: `e835166a9` (`Tighten observation API validation`)

### 7. Legacy metadata could remain stale after a rebuild

Rebuilding a `SingleObservationSet` from its backing collection did not refresh cached dependent-variable bookkeeping. Asking an empty legacy collection for time bounds also took an unsafe legacy path.

Corrections:

- refresh dependent-variable bookkeeping when the single-set facade rebuilds;
- delegate empty-collection time-bound queries safely to the dataset implementation.

Regression tests cover both cases.

Commit: `c9dbc1986` (`Refresh legacy observation metadata`)

## CI-only follow-up corrections

The first full validation run reached the newly extended observation-dataset test target and exposed two mistakes in the tests themselves: a templated expression containing a comma needed parentheses inside a Boost.Test macro, and a test used a nonexistent legacy method name. These were corrected in `2ecb33241` (`Fix observation dataset regression tests`).

The next run reached the same target and showed that the legacy `setLinkEnds` facade unnecessarily required a mutable `LinkDefinition`. Making the argument a const reference is source-compatible and matches the underlying setter. This was corrected in `217dc64cd` (`Make legacy link metadata setter const-correct`).

The following macOS run built successfully and executed all 346 test executables. Its only failing executable reported two failures where older test fixtures supplied deliberately asymmetric values to APIs now correctly validating weight matrices as symmetric. Those tests exercised projection and precedence, not acceptance of asymmetric weight matrices. Their fixtures were made symmetric without reducing that coverage in `29ea922b3` (`Make weight matrix test fixtures valid`).

## Commit map

1. `1005ecb8e` — Fix Boost.Test runtime linkage
2. `5ba1836eb` — Enforce observation dataset invariants
3. `9f1856dac` — Preserve observation ordering in residual calculations
4. `f8ccb4635` — Isolate copied observation dataset metadata
5. `63c95c62a` — Preserve full weights in legacy observation transforms
6. `85d178554` — Keep observation facades synchronized and lifetime-safe
7. `e835166a9` — Tighten observation API validation
8. `c9dbc1986` — Refresh legacy observation metadata
9. `2ecb33241` — Fix observation dataset regression tests
10. `217dc64cd` — Make legacy link metadata setter const-correct
11. `29ea922b3` — Make weight matrix test fixtures valid

## Validation result

The Linux, macOS, and Windows jobs in run `33860641446` all completed successfully for commit `29ea922b3d701a39b1b282e48b6eb8ae707163e9`. The push and pull-request pre-commit workflows also passed.
