# PR 1092: handover for compilation and testing (15 October 2026)

Start here. This file is self-contained. Background documents (read only if needed):
`REVIEW-IMPLEMENTATION-STATUS.md` (earlier stages, verification log), `REVIEW-DONE.md` (the review),
`REVIEW-RISKIER-REDESIGNS.md` (larger redesigns not done).

## 0. Rules you must follow
- **Never push** to any remote, publish or upload unless Dominic explicitly instructs that exact action in your own
  conversation (see `/home/dominic/AGENTS.md`). The commits below are local only.
- Build with **at most 8 parallel jobs** (`make -j8`, never more).
- No `Co-Authored-By` trailers or AI attribution lines in commit messages.
- Do not use bare `git stash` (the stash is shared between worktrees).
- Nothing may change the behaviour of `feature/data-refactor` (`3070b20e83e18a749767fe0332a5280464f59b16`) code, and
  numerical quality must not be compromised. Changing APIs that are new in the PR is allowed; rewriting tests is allowed.
- The pre-commit hooks (clang-format, black) reformat files. When a commit fails with "files were modified by this
  hook", re-add the files and commit again.

## 1. Where everything is
- Worktree: `/home/dominic/Tudat/tudat-monorepo/tudatpy-pr1092-external-review`, local branch
  `review/pr1092-external-simplify`.
- Remote branch: `origin/fix/earth_deformation_rotation` is at `92b9c7b16` (last push). The worktree is 4 commits ahead:

  | Commit | Content |
  | --- | --- |
  | `023a5323f` | Items 1–5: atmospheric, ocean-pole and reference-point policies replaced by single behaviours |
  | `9024aa95b` | Item 6 (C4): stations/deformation use the body rotation model as it is; warning for non-TDB input scale |
  | `ceade07f5` | Item 7 (C5): base ordering of frame processing restored in `createSystemOfBodies` |
  | `c1afc61f8` | Documentation pass of the Python API added by the PR |

- Build directory: `/home/dominic/Tudat/pr1092-review-build` (already configured; same options as
  `/home/dominic/Tudat/earth-station-python-build`; serialization is OFF in this configuration).
- Resource files that must be in `~/.tudat/resource` (already installed locally; flagged for a tudat-resources PR, staged
  in `/home/dominic/Tudat/tudat-monorepo/tudat-resources-additions/`):
  `earth_orientation/desaiSibois2016*.txt` (4 files + README) and `station_locations/Transfo-ITRF2014_ITRFs.txt`,
  `station_locations/Transfo-ITRF2020_TRFs.txt`.

## 2. What was changed, and why

Dominic's decisions (15 October) for the policies the PR had introduced, and for two compatibility findings:

### Item 1 — atmospheric grid seam/pole consistency (`023a5323f`)
- Removed `AtmosphericTideEndpointPolicy` (C++ enum, settings field, Python enum and `endpoint_policy` attribute).
- `AtmosphericTideGridInterpolator( grid )` now always averages the duplicated 0/360° longitude columns and replaces each
  polar row by its Cartesian mean (on a model copy; the provider grid is retained).
- New sanity check: `AtmosphericTideGridInterpolator::maximumEndpointDiscrepancyMetres = 1e-5`. A grid whose duplicates
  differ by more than 10 µm is rejected. The official RP03 grid differs by at most 8e-7 m (seam) and ~4e-7 m (poles).
- Why: the old default rejected the official IERS RP03 grid. Dominic: averaging is fine, the effect is very minor.
- Files: `include/tudat/astro/ground_stations/atmosphericTideLoading.h`, `src/tudat/astro/ground_stations/atmosphericTideLoading.cpp`,
  `atmosphericTideDeformation.{h,cpp}`.

### Item 2 — RP03 default conventions (`023a5323f`)
- In the grid constructor of `AtmosphericTideDeformation`: if the grid metadata has `modelName == "Ray-Ponte 2003"` and
  `greensFunction == "Farrell elastic Green functions"` (both set by the reader from the file header), undeclared lookup
  latitude convention defaults to `geocentric`, and undeclared displacement basis to `geocentricRadialNorthEast`.
  Explicit declarations still take precedence; other grids still require explicit declarations.
- Why: users had to know an undocumented convention. The provider does not document it (checked online); geodetic and
  geocentric latitude differ by ≤ 0.19° on the 1° grid, i.e. micrometres in displacement. Geocentric needs no body shape.

### Item 3 — atmospheric output origin (`023a5323f`)
- `AtmosphericLoadingOriginPolicy` replaced by `bool AtmosphericTideEvaluationSettings::convertToCentreOfMass = false`
  (Python `convert_to_centre_of_mass`). Behaviour unchanged (the old default was `preserveSourceOrigin`).

### Item 4 — ocean pole tide latitudes outside the grid (`023a5323f`)
- Removed `OceanPoleTideLatitudePolicy` (C++ and Python, including the `latitude_policy` arguments of `ocean_pole_tidal`,
  `ocean_pole_tidal_from_data`, `interpolateOceanPoleTideCoefficients`, the `OceanPoleTideDeformation` constructor and
  the settings class).
- Latitudes up to **one grid step** beyond the first/last grid row now use that row (covers ±89.75°…±90° of the
  official grid). Latitudes further outside still raise `std::out_of_range`. The one-step limit was added so that a
  station far outside a regional extract is not silently clamped.
- `OceanPoleTideDeformation` constructor is now `( grid, polarMotionCalculator, converter = nullptr )`.

### Item 5 — station reference-point state at an offset switch (`023a5323f`)
- Removed `StationReferencePointRatePolicy` (C++ and Python, including `rate_policy` arguments/properties of
  `StationReferencePointOffsetHistory` and `sinex_antenna_reference_point`).
- At an equipment/offset switch, the state now uses the new (right-hand) record, as elsewhere in the PR. Previously the
  default threw.

### Item 6 — C4: no hidden TDB copies of the Earth rotation (`9024aa95b`)
- Removed the three places that built a second `GcrsToItrsRotationModel` with `tdb_scale`:
  `createReferencePointEphemeris` (`createGroundStations.h`), `deformationRotationInTdb` (`createBodyDeformationModel.cpp`)
  and the relativistic station motion factory (`createGroundStations.cpp`). All now use the body's rotation model as is,
  which restores the `feature/data-refactor` behaviour.
- Added a warning (stderr) in `createRotationModelEarthOrientation.cpp` when a GCRS-ITRS rotation model is created from
  settings with a non-TDB input time scale; documented in `createRotationModel.h`.
- Why: Dominic does not want separate rotation models just to support different input time scales. With TDB (the default,
  and the only option from Python) results are identical.
- Test changed: `unitTestRelativisticStationPosition.cpp` (`native_public_factory_uses_selected_external_bodies_and_tdb_rotation`)
  looped over all input time scales expecting the TDB copy; it now uses only `tdb_scale`.

### Item 7 — C5: frame-processing order (`ceade07f5`)
- `createSystemOfBodies` (`createBodiesFactory.h`) again calls `processBodyFrameDefinitions` **after** ground stations
  and cameras are created, as on `feature/data-refactor`.
- The relativistic station motion factory (`createGroundStations.cpp`) needed frame processing only for its check that
  external-body ephemeris origin conversions are configured. That check now runs on the first evaluation of each
  barycentric state function (a one-time check captured in the lambda). Axis checks remain at creation.
- Test to watch: `unitTestRelativisticStationPosition.cpp`,
  `body_settings_factory_retains_environment_and_converts_external_body_origin_before_station_creation`.

### Documentation pass (`c1afc61f8`)
- Dominic's request: tudatpy style; function/class discussion at docstring level; module rst pages only a top-level
  description (plus at most one or two small critical examples); all arguments and types documented as in the rest of
  tudat.
- Module rst pages (`docs/tudatpy/source/...`) now contain an overview (or `automodule` of the package `__init__.py`
  docstring, as `space_weather` does) plus `Functions` / `Enumerations` / `Classes` autosummary and autodoc listings.
- Docstrings rewritten with numpy-style `Parameters`/`Returns`, enum value descriptions, and attribute docstrings with
  `:type:` (and `**read-only**` where applicable) for: `data_input.environment_data.{earth_orientation,
  atmospheric_loading, ocean_loading, station_loading, station_history, terrestrial_frames}`,
  `dynamics.environment` station-model classes (`expose_station_models.cpp`), `dynamics.environment_setup.{ground_station,
  shape_deformation, rotation_model}` and `add_ground_station`/`convert_ground_station_state_between_itrf_frames`
  (`expose_environment_setup.cpp`), and the terrestrial frame classes in `astro.frame_conversion`.
- New functions/classes that were missing from the rst listings were added (for instance all new `ground_station`
  functions and `StationHistoryGroundStationMotionSettings`).
- `GcrsToItrsRotationModelSettings.set_eop_data` now has a named argument `eop_data` (positional calls still work).

## 3. What you must do (build and test)

Nothing of the above has been compiled. Expect some compile errors in the binding files (large docstring rewrites).

1. Build the native libraries and the selected native tests (at most 8 jobs):
   ```sh
   cd /home/dominic/Tudat/pr1092-review-build
   make -j8 $(cat /tmp/claude-1000/-home-dominic/da490dc0-d594-4ff0-aa23-705088333291/scratchpad/selected.txt) \
       test_reference_frames_ReferenceFrameTransformations test_io_TerrestrialFrameTransformationData \
       test_simulation_EnvironmentModelSetup
   ```
   If the scratchpad file is gone, the 26 targets are: `test_earth_orientation_{EarthOrientationCalculator,
   EarthOrientationDataEvaluator, EarthRotationDerivatives, EopReader, HighFrequencyEarthOrientation,
   PolarMotionCalculator, ShortPeriodEopCorrections, TimeScaleConverter}`, `test_ephemerides_ItrsToGcrsRotationModel`,
   `test_ground_stations_{AtmosphericTideDeformation, AtmosphericTideLoading, BodyDeformationVelocity,
   CompleteEarthStationVelocity, DefaultGroundStationSettings, EarthDeformationTimeRouting, GroundStationState,
   Iers2010DeformationModel, OceanPoleTideDisplacement, OceanTideDisplacement, PoleTideDisplacement,
   RelativisticStationPosition, StationHistoryFrameTransformation, StationHistoryTrajectory, StationLoadingDeformation,
   StationReferencePoint}`, `test_reference_frames_TerrestrialFrameTransformation`. Verify the exact target names with
   `make help | grep test_` if a name is not found.
2. Run each test executable serially; they are built as `/home/dominic/Tudat/pr1092-review-build/tests/<target name>`
   (for instance `tests/test_ground_stations_StationReferencePoint`). All passed before these 4 commits.
   Tests most likely to need attention (they were edited without compiling):
   - `unitTestAtmosphericTideLoading.cpp`, `unitTestAtmosphericTideDeformation.cpp` (new case
     `rp03GridDefaultsToGeocentricConventions`; seam/pole cases now check the 10 µm rejection),
   - `unitTestOceanPoleTideDisplacement.cpp` (clamp and regional-extract rejection), `unitTestEarthDeformationTimeRouting.cpp`,
   - `unitTestStationReferencePoint.cpp` (state at the switch: `getState( atSwitch ).head< 3 >( )` must equal
     `getPosition( atSwitch )`),
   - `unitTestRelativisticStationPosition.cpp` (items 6 and 7),
   - `test_simulation_EnvironmentModelSetup` and anything else using `createSystemOfBodies` with ground stations
     (item 7 changes the ordering for every environment).
   Rewriting tests to the new behaviour is allowed; changing behaviour that exists on `feature/data-refactor` is not.
   If a test exposes a real numerical problem, report it rather than loosening tolerances.
3. Build the Python kernel and run the Python tests:
   ```sh
   cd /home/dominic/Tudat/pr1092-review-build && make -j8 kernel
   cd /home/dominic/Tudat/tudat-monorepo/tudatpy-pr1092-external-review
   PYTHONPATH=/home/dominic/Tudat/pr1092-review-build/src python -c "import tudatpy.dynamics.environment_setup, tudatpy.data_input.environment_data"
   PYTHONPATH=/home/dominic/Tudat/pr1092-review-build/src python -m pytest -q -p no:cacheprovider \
       tests/test_tudatpy/test_data_input_loading_settings.py tests/test_tudatpy/test_station_reference_point_bindings.py \
       tests/test_tudatpy/test_data_input_station_history.py tests/test_tudatpy/test_data_input_earth_orientation.py \
       tests/test_tudatpy/test_data_input_resource_paths.py tests/architecture --no-remote-data
   ```
   Then the full suite: `... python -m pytest -q -p no:cacheprovider tests/test_tudatpy --no-remote-data`. Before these
   commits it had 214 failures/19 errors, all environmental (serialization OFF; missing tracking-data files). Compare
   against that; any new failure involving the changed modules must be fixed.
4. Check that every documented name exists: in Python, for each module touched in `c1afc61f8`, check that all names in
   the rst autosummary lists resolve (for instance with a short script using `getattr`). Optionally build the Sphinx
   docs (`docs/tudatpy`) and check for autodoc warnings on the touched pages.
5. Record the results in `REVIEW-IMPLEMENTATION-STATUS.md` (verification log) and in the repository copy
   `docs/pr1092_review_notes.md`, commit locally, and **ask Dominic before pushing**.

## 4. Not implemented (still open)
- S5 (station-history time model to plain `Time`), S6 (remove `TransformedStationHistoryTrajectory` /
  `StationReferencePointTrajectory` wrappers), and the Python exposure redesign (settings factories instead of mutable
  classes; one factory per deformation model; `ign_station_history` returning settings). Dominic has since said API
  redesigns and test rewrites are acceptable, so these are no longer blocked by compatibility, only by effort; they were
  not part of this round.
- S2 (single EOP data route): base API, so only via deprecation; first make the observation tests that need tracking data
  (ODF, PSF, troposphere) runnable, because the shared time-converter code is used there.
- S8 (shared reading helpers, per-record provenance), S9 (single latitude/longitude grid interpolator), S10 (time
  converter in deformation models), remaining S11 items, comment/format pass of the C++ comments.
- The tudat-resources PR for the new resource files (section 1).
