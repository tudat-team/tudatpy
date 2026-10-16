# SLR and LLR analysis recovery, 16 October 2026

For continuing the analysis, read [CONTINUE-HERE.md](../CONTINUE-HERE.md) first.
This file describes the recovery and the historical support archive.

The interrupted analysis was recovered from **Resume Tudat framework work**,
thread `01a13aba-412d-71e2-b12e-5dc837528491`. The SLR/LLR task began on
15 October at 19:58 UTC; the last saved answer was on 16 October at 13:43 UTC.
The local transcript is
`/home/dominic/.codex/sessions/2026/10/14/rollout-2026-10-14T16-04-07-01a13aba-412d-71e2-b12e-5dc837528491.jsonl`.
A readable export of the 38 user and 131 visible assistant messages in the
laser-ranging portion is saved locally at
`/home/dominic/Tudat/laser-ranging-work/recovery/conversation-20261015-16.txt`.
The conversation export remains local.

The integration checkout is `tudatpy-laser-ranging`, on
`feature/high-accuracy-laser-ranging`, tracking the same branch on
`https://github.com/tudat-team/tudatpy`. At recovery, local HEAD was
`ed2a720b2`; the remote was `c5890681f`. These three completed commits awaited
publication:

- `3b9362cfb`: annual table with two explicit Grasse outliers omitted.
- `635a5a635`: annual filtered residual figures.
- `ed2a720b2`: separate columns for APOLLO, Grasse green and Grasse infrared.

## Current analysis and results

- [Analysis guide](../README.md) and [native model validation](../../../docs/laser_ranging_validation.rst).
- [SLR driver](../fixed_orbit_slr.py), [input manifest](../dataset.json),
  [stable cubic comparison over two weeks](../results/stable_cubic_two_weeks/summary.json)
  and [Earth orientation/deformation comparison](../results/high_tidal_two_weeks/README.md).
- [LLR driver](../llr_residuals.py), [long-arc driver](../llr_long_arc.py)
  and [annual input manifest](../llr_year_2017_dataset.json).
- [Calendar-2017 raw results](../results/year_2017/README.md),
  [filtered statistics and three-column figures](../results/year_2017/outlier_filtered/README.md)
  and [temporal comparison](../results/year_2017_temporal/README.md).
- [Annual source audit](../llr_year_2017_sources/README.md) and
  [independent output verification](../llr_year_2017_sources/annual_output_verification.md).

The annual calculations contain 1,579 observations per model. The derived
filtered table and figures retain 1,577: 311 APOLLO, 92 Grasse green and 1,174
Grasse infrared. The two additional June exclusions are recorded by exact
source identity in the filtered-results README. The earlier March exclusion,
invalid source records, duplicate provider record and weather limitations remain
documented in the annual source audit.

The final recovered exchange clarified that **glass** refers to the optical
travel time inside the lunar corner-cube prisms. The separate-glass case adds
the explicit internal optical path in the JPL model. The interpretation of
the published INPOP reflector coordinate reference point remains unresolved;
the smaller residuals in the no-separate-glass case do not resolve that physical
convention. Both cases remain available for comparison. The fixed-ephemeris
LLR analysis still has residual structure; these artifacts do not establish
measurement-noise-level accuracy.

## Preserved support files

[LASER-RANGING-STATUS.md](LASER-RANGING-STATUS.md) is the unchanged historical
handover log from outside Git, last updated at 13:05 UTC. Its pending-work and
push statements describe successive historical checkpoints. This recovery
index supplies the later conversation context.

`workspace/` mirrors selected paths from
`/home/dominic/Tudat/laser-ranging-work/`:

- `verification/`: all Python and shell analysis probes, run drivers and
  packaging scripts, plus their local configuration files and small audit results.
- `logs/`: the Python scripts used to check bindings and run focused native suites.
- `research/`: the SLR/LLR reference-case notes, source metadata, loading
  configuration, independent-software notes and optical-reference drivers.
- `data/`: download manifests for the retained external provider caches.
- `results/`: the saved executed Python source from the early APOLLO control.

The six files under
`worktrees/tudatpy-laser-llr/analysis/laser_ranging/results/stable_cubic_20261016/`
preserve an untracked historical worker result package. Its recorded native
revision is `b1aa772503070347e8ac0b263c8de9278d390c17`. It is a historical
control, separate from the current annual results linked above.

[manifest.json](manifest.json) inventories every imported file with its original
location, original SHA-256, archived SHA-256 and byte counts. Python formatting
follows the repository's Black hook; the recovery checks compare the original
and archived Python syntax trees. Other imported files retain their exact bytes.
All 85 Python/shell scripts and Markdown notes found in the external workspace
are represented, including historical variants.

The support scripts retain their original absolute paths, selected revisions,
output locations and execution assumptions. They document the actual probes
used during the investigation. For a replay, restore the corresponding workspace
layout and matching kernel/data, or adapt the paths explicitly. The portable
current analysis drivers and documented commands are linked above. Provider
data, ephemeris kernels, compiled binaries and research papers remain in the
local caches; the manifests retain their source and checksum information.
The local `.gitattributes` preserves archived line endings and the original
whitespace in the IERS reference sources, compiler output, CSV and saved patch.

This recovery packages existing work. It does not change model formulas,
residual CSVs, observation selections or the current figures.

The [recovery validation](validation.json) records script/configuration checks,
all three annual output audits and agreement with the filtered statistics.
