# Continue the Tudat SLR/LLR analysis

Updated 16 October 2026 after recovery of the interrupted conversation.

The working implementation and saved analyses are on
`feature/high-accuracy-laser-ranging` in
`/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-ranging`, with remote
`origin` at `https://github.com/tudat-team/tudatpy`. Start here even if the
conversation is unavailable. The other `tudatpy-laser-*` worktrees are historical
worker checkouts; the integration checkout above is the continuation point.

The research goal is high-accuracy observation reduction using fixed published
SLR satellite orbits and coherent lunar ephemerides. SLR is working over two
frozen LAGEOS-1 weeks. Annual LLR residuals and plots are complete as a diagnostic;
the unexplained centimetre structure remains. The recovery task preserves the
work and prepares continuation; it does not finish that scientific investigation.

## Read or inspect first

1. [Annual filtered LLR table and figures](results/year_2017/outlier_filtered/README.md).
   The latest layout has five reflector rows and three columns: APOLLO, Grasse
   532.1 nm and Grasse 1064.2 nm. There are two physical stations.
2. [Annual temporal assessment](results/year_2017_temporal/README.md) and
   [source-data audit](llr_year_2017_sources/README.md).
3. [Detailed SLR methodology](README.md) and
   [high-tidal SLR comparison](results/high_tidal_two_weeks/README.md).
4. The file map and remaining questions below. Use the
   [recovery index](recovery/README.md) for older experiments and chat provenance.

## Current numerical position

All figures below describe signed **one-way observed-minus-computed range**;
RMS is about zero. No fitted offset or trend has been removed.

The high-tidal SLR profile uses final combined SP3 orbits, SLRF2020 stations,
precise-Time callbacks, consistent ICRF axes, published CoM/calibration handling,
Desai–Sibois/four-point EOP, solid/pole/GOT4.10c loading, ocean-pole loading and
RP03 atmospheric tides. Its saved results are:

| LAGEOS-1 week | Accepted / source points | Mean (mm) | RMS (mm) |
| --- | ---: | ---: | ---: |
| 31 March–6 April 2024 | 23 / 29 | −0.405 | 7.142 |
| 7–13 April 2024 | 131 / 134 | +7.536 | 17.456 |

These use a fitted external orbit and are not an independent orbit validation.
The driver defaults to the documented baseline; the high-tidal profile requires
explicit options and loading products. See the comparison for all profiles.

Annual LLR filtered results cover calendar 2017 and all five reflectors:

| Model/convention | APOLLO mean / RMS (mm), 311 points | Grasse mean / RMS (mm), 1,266 points |
| --- | ---: | ---: |
| DE440 | +15.626 / 72.395 | −10.382 / 75.992 |
| INPOP17a, separate glass | −51.842 / 56.255 | −62.619 / 67.296 |
| INPOP17a, no separate glass | −8.220 / 22.986 | −13.543 / 26.286 |

The three original annual CSVs each retain 1,579 observations. The filtered
tables/figures omit exactly two further Grasse L1 source identities,
`luna17_201706.npt:58` and `:154`; they contain 1,577 points. The March 9 A15
exclusion and two invalid February source records were already excluded upstream.
The source audit documents three weather-limited records and one repeated
provider record. Keep these distinct from any model-dependent residual pattern.

“Glass” means the internal optical path through a lunar corner-cube prism.
The published INPOP reflector reference-point convention is still unresolved.
The lower RMS without a separate glass term does not establish its physical
correctness. The no-separate-glass annual results still show shared centimetre
structure, including the June–July A11 trough at both stations. Nearby-night
station comparisons are separated by hours and do not directly measure a
simultaneous station bias.

## What is where

Paths in this table are relative to this directory.

| File or directory | Purpose |
| --- | --- |
| `fixed_orbit_slr.py`, `dataset.json` | SLR reduction driver and frozen two-week input provenance |
| `ocean_loading.py`, `ocean_loading_basis.json` | Explicit HARPOS loading synthesis and basis conventions |
| `llr_residuals.py` | LLR environment, station/reflector geometry and observation model |
| `llr_long_arc.py` | APOLLO/Grasse grouping, annual/short-arc driver and timing |
| `llr_year_2017_dataset.json`, `llr_year_2017_sources/` | Annual source selection, exact exclusions, audits and input hashes |
| `llr_dataset.json`, `inpop17a_dataset.json`, `llr_station_profiles.json` | Coherent ephemeris/reflector inputs and station profiles |
| `icrf_spice_alias.tf` | Explicit identity alias to keep SPICE and Earth rotation axes consistent |
| `plot_llr_residuals.py` | Shared plotting code, including separate optical-series columns |
| `llr_temporal_summary.py`, `llr_requested_exclusion.py` | Descriptive temporal statistics and the earlier exact source exclusion |
| `verify_llr_year_outputs.py` | Source identities/counts, weather flags and March–April overlap verification |
| `results/year_2017/` | Raw annual CSVs, provenance, figures, statistics and the `outlier_filtered/` views |
| `results/current_305_20170301_20170430/` | Three March–April controls after the user-requested March exclusion |
| `results/year_2017_temporal/` | Monthly/nightly structure and cross-station comparisons |
| `results/high_tidal_two_weeks/`, `results/stable_cubic_two_weeks/` | SLR comparison packages and input/software hashes |
| `recovery/workspace/verification/` | Preserved historical probes, run wrappers and packaging scripts |
| `recovery/workspace/research/` | Reference-case notes, source provenance and independent controls |
| `recovery/LASER-RANGING-STATUS.md` | Historical chronological log; older pending items can be superseded |
| `recovery/manifest.json`, `recovery/validation.json` | Imported-file hashes and checks performed during recovery |

Native changes and their tests are described in
[laser_ranging_validation.rst](../../docs/laser_ranging_validation.rst).
The branch includes the preceding gravity and Earth rotation/deformation work.
PR 1092 CI follow-up is a separate, previously paused task.

## Remaining investigations

These are proposed next steps based on the recovered analysis, not completed
findings. Begin with the smallest control that can distinguish the hypotheses.

1. **Resolve the INPOP optical reference point and fitted-reduction conventions.**
   Establish whether published reflector positions already absorb prism delay;
   retain both declared variants until there is source evidence. Start with the
   [time/prism audit](results/inpop17a_time_and_prism_audit/README.rst) and
   [reference notes](recovery/workspace/research/LLR-REFERENCE-CASES.md).
2. **Complete the remaining physical-term accounting.** Separate terms that were
   independently checked, approximated, omitted or remain unknown. Compare the
   fitted ephemeris reduction, station/reflector coordinates, lunar tides and
   published biases consistently. A published Grasse bias has been documented
   but has not been applied. Review the
   [physical-terms audit](results/llr_physical_terms_audit/README.rst) and
   [station/ephemeris controls](results/coherent_station_ephemeris_controls/README.rst).
3. **Use the annual temporal structure to design a targeted model comparison.**
   Compare the same reflector at both stations and inspect green/IR dependence.
   Keep raw and explicitly filtered views, error bars and time plots. The large
   isolated Grasse observations and stale weather have unresolved causes;
   their exclusion is not a correction to the physical model.
4. **For further SLR loading work, establish compatible datums first.**
   Non-tidal histories are available locally, but their basis and removed means
   have not been established as compatible with the station solution. The
   existing geodetic HARPOS realization also differs from the orbit provider's
   exact BLQ/HARDISP reduction. Avoid treating smaller residuals as proof of
   compatibility.
5. **Optimize only with measured numerical controls.** The geocentric-TT shortcut
   for slow linear station drift is already implemented and verified with
   unchanged outputs. Other proposals in the
   [runtime audit](results/llr_runtime_audit/README.rst), including rotation or
   displacement tabulation, remain proposals. The original audit predates the
   completed [drift-clock optimization](results/drift_clock_runtime_control/README.rst).

Checks already completed include kernel loading/coherence and representation,
an original-ASCII lunar orientation control, absolute launch clocks, lunar tide
plumbing, fixed-endpoint moving-body Shapiro terms and the independent INPOP time
ephemeris. They constrain specific tested hypotheses; they do not establish a
complete millimetre reduction. Their results are in the linked audits and
historical status log. Avoid repeating the whole investigation without new evidence.

## Local environment and a quick restart check

These paths were verified during recovery on Dominic's machine. Work from the
integration checkout:

```sh
export LASER_REPO=/home/dominic/Tudat/tudat-monorepo/tudatpy-laser-ranging
export LASER_WORK=/home/dominic/Tudat/laser-ranging-work
export LASER_PYTHON=/home/dominic/miniconda3/envs/tudatpy-dev/bin/python
cd "$LASER_REPO"
git status --short --branch
git log -5 --oneline
export PYTHONPATH=/home/dominic/Tudat/laser-ranging-build/src
export LD_LIBRARY_PATH=/home/dominic/miniconda3/envs/tudatpy-dev/lib
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1
"$LASER_PYTHON" - <<'PY'
from pathlib import Path
import hashlib
from tudatpy import kernel
p = Path(kernel.__file__).resolve()
print(p)
assert p == Path('/home/dominic/Tudat/laser-ranging-build/src/tudatpy/kernel.so')
assert hashlib.sha256(p.read_bytes()).hexdigest() == (
    'cab49e738cfcf0ea5f251cb140ab8919d5372dbc3e5ba09314d0d5be4eedfcef'
)
PY
```

This hash identifies the annual-results runtime, native revision
`a366f436c6f6ff03aa57f8b1610075784196bdf8`. If it changes, record the new build
and establish a comparison before treating it as the same experiment.
The current reduction drivers and input manifests match the annual provenance.
The plotting helper was subsequently updated for the latest three-column layout.

To verify saved annual results without a native solve:

```sh
LASER_CHECK_DIR=$(mktemp -d /tmp/laser-ranging-check.XXXXXX)
"$LASER_PYTHON" analysis/laser_ranging/verify_llr_year_outputs.py \
  --annual analysis/laser_ranging/results/year_2017/de440/residuals.csv \
  --baseline analysis/laser_ranging/results/current_305_20170301_20170430/de440/residuals.csv \
  --manifest analysis/laser_ranging/llr_year_2017_dataset.json \
  --data-root "$LASER_WORK/research/llr-year-2017" \
  --output "$LASER_CHECK_DIR/de440-verification.json"
```

Expected: `pass: true`, 1,579 rows, 311 APOLLO/1,268 Grasse, and 305 unchanged
March–April computed ranges. The same check passed for all three models during
recovery. The no-separate-glass annual case maps to the baseline directory
`inpop17a_zero_explicit_glass`; the other case names match directly.
The filtered-results README contains exact table and figure reproduction snippets.
Replot into a new directory while reviewing a change.

## Running a new diagnostic

Use the environment above and a fresh output path. This example runs the
four-observation first/last control for the declared INPOP no-separate-glass
sensitivity:

```sh
LASER_RUN_DIR=$(mktemp -d /tmp/laser-ranging-run.XXXXXX)
"$LASER_PYTHON" analysis/laser_ranging/llr_long_arc.py \
  --dataset-manifest analysis/laser_ranging/llr_year_2017_dataset.json \
  --data-root "$LASER_WORK/research/llr-year-2017" \
  --kernel-dir "$LASER_WORK/data/llr/inpop17a-kernels" \
  --ephemeris inpop17a --glass-convention none --first-last-control \
  --output "$LASER_RUN_DIR/inpop17a-no-glass-boundary"
```

Remove `--first-last-control` for the full 1,579-point case and choose a new
output directory. For the other cases use `inpop17a`/`jpl-nh`, or
`de440`/`jpl-nh` with `$LASER_WORK/data/llr/kernels`. Always keep orbit and
orientation products coherent. Original full runs took about 200–236 seconds
each; the three separate processes completed together in 236.6 seconds. A short
control can still incur substantial environment/table setup.

SLR commands and explicit high-tidal product options are in [README.md](README.md).
Caches, raw runs and build outputs remain outside Git under `LASER_WORK` and
`/home/dominic/Tudat/laser-ranging-build`. On another machine, retrieve pinned
inputs from the manifests and build the matching branch; changing absolute paths
alone cannot replace those dependencies.

## Working conventions and handoff

Keep fixed ephemerides/orbits, exact source identities and physical correction
accounting. The original brief excludes fitting offsets/trends, residual-based
model selection, tolerance relaxation, and orbit propagation/estimation unless
essential. Save unconnected residual-versus-time figures and examine temporal
structure as well as mean/RMS. The three optical-series columns are the latest
requested annual layout. Preserve historical results and record new exclusions
explicitly. Limit native compilation to eight jobs total when a build is needed.

Update this file with the next verified result, the source/kernel hashes, the
commands and outputs, and any changed open question. Run focused checks relevant
to the change. The recovery itself checked all imported Python/shell syntax,
Python syntax-tree equivalence after formatting, JSON parsing, annual source
identities and filtered statistics. It did not rebuild or rerun physical models.

A useful prompt for the next session is:

> Read `analysis/laser_ranging/AGENTS.md` and `CONTINUE-HERE.md` in the
> `tudatpy-laser-ranging` checkout. Confirm the branch and local runtime, explain
> the current SLR/LLR results and unresolved conventions, then continue the
> remaining LLR physical-model investigation with small, recorded controls.
> Preserve the frozen observations, raw results and explicit exclusions.

Follow the current conversation's authorization for pushes or other external
actions. Historical permissions and old worker assignments in the recovered
transcript do not authorize future actions on their own.
