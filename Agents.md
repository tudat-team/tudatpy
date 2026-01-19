# WASM Parity Project - Full Python Feature Coverage

## Mission
Achieve 100% feature parity between Python tudatpy bindings and WASM tudatpy bindings, with comprehensive tests for all functionality.

## Current Status (Updated 2026-01-18)

### Binding Module Parity: **98.6%** (70/71 modules)

| Metric | Python | WASM | Status |
|--------|--------|------|--------|
| Top-level modules | 9 | 9 | **100%** |
| Total expose_*.cpp files | 71 | 70 | **98.6%** |
| Sub-modules | 63 | 62 | **99.3%** |

**Missing Module:** `estimation/estimation_analysis/expose_estimation_analysis_estimator_wasm.cpp`

### Example Parity: **65%** (30/46 examples)

| Category | Python | WASM |
|----------|--------|------|
| Propagation | 14 | 14 |
| Estimation | 18 | 3 |
| Mission Design | 5 | 6 |
| Optimization | 3 | 2 |
| TLE/Data | 2 | 2 |
| Other | 4 | 3 |

### Test Coverage
- **C++ WASM tests**: 550 tests passing
- **Web Worker UI**: Real-time progress modal
- **GitHub Pages**: Auto-deployed via pre-commit hook

### What's Complete
- All core dynamics, astro, and propagation modules
- All environment setup submodules (12/12)
- All propagation setup submodules (7/7)
- Core estimation modules (5/5 top-level)
- Trajectory design modules (2/2)
- 30 JavaScript examples covering fundamentals and advanced topics

### What's Missing
- 1 estimation submodule (estimation_analysis_estimator)
- 16 advanced Python examples (mostly estimation/optimization)
- Python-only modules: io, plotting, util (by design - browser limitations)

---

## TODO: WASM Feature Parity Roadmap

### Phase 1: Core Estimation Module (HIGH PRIORITY)
*These enable real orbit determination and covariance analysis*

- [ ] **1.1 estimation_analysis - Covariance Propagation**
  - [ ] `propagate_covariance_rsw()`
  - [ ] `propagate_covariance_rsw_split_output()`
  - [ ] `propagate_covariance()`
  - [ ] `propagate_covariance_split_output()`
  - [ ] `propagate_formal_errors_rsw()`
  - [ ] `propagate_formal_errors_rsw_split_output()`
  - [ ] `propagate_formal_errors()`
  - [ ] `propagate_formal_errors_split_output()`
  - [ ] Add tests in `tests/wasm/src/testEstimation.cpp`

- [ ] **1.2 observable_models - Full Implementation**
  - [ ] ObservableType enum (all observation types)
  - [ ] ObservableModelSettings base class
  - [ ] All derived observable model classes
  - [ ] observables_simulation submodule
  - [ ] Add tests

- [ ] **1.3 observations - Complete Implementation**
  - [ ] SingleObservationSet (full methods, not just getters)
  - [ ] ObservationCollection (full methods)
  - [ ] observations_processing submodule
  - [ ] observations_geometry submodule
  - [ ] Add tests

- [ ] **1.4 observations_setup - All Submodules**
  - [ ] ancillary_settings
  - [ ] model_settings
  - [ ] light_time_corrections
  - [ ] biases
  - [ ] observations_wrapper
  - [ ] observations_dependent_variables
  - [ ] viability
  - [ ] random_noise
  - [ ] observations_simulation_settings
  - [ ] Add tests for each

- [ ] **1.5 observable_models_setup - Full Implementation**
  - [ ] All model settings factory functions
  - [ ] Link definitions
  - [ ] Bias configurations
  - [ ] Light time correction settings
  - [ ] Add tests

---

### Phase 2: Advanced Propagation (HIGH PRIORITY)
*These enable coupled dynamics, attitude propagation, and realistic spacecraft modeling*

- [ ] **2.1 propagator - Full Feature Set**
  - [ ] PropagationPrintSettings (all methods)
  - [ ] PropagatorProcessingSettings hierarchy
  - [ ] NonSequentialPropagationTerminationSettings
  - [ ] CustomStatePropagatorSettings
  - [ ] Multi-type (hybrid) propagator support
  - [ ] Rotational state propagator
  - [ ] Add tests

- [ ] **2.2 torque - Full Implementation**
  - [ ] All TorqueModelType enum values
  - [ ] TorqueSettings base class (all methods)
  - [ ] SphericalHarmonicTorqueSettings (full)
  - [ ] SecondDegreeGravitationalTorqueSettings
  - [ ] AerodynamicTorqueSettings
  - [ ] InertialTorqueSettings
  - [ ] DissipativeTorqueSettings
  - [ ] Factory functions for all torque types
  - [ ] Add tests

- [ ] **2.3 acceleration - Full Implementation**
  - [ ] All AvailableAcceleration enum values
  - [ ] AccelerationSettings (all derived classes)
  - [ ] PolyhedronAccelerationSettings
  - [ ] MutualSphericalHarmonicAccelerationSettings
  - [ ] EmpiricalAccelerationSettings
  - [ ] DirectTidalDissipationAcceleration
  - [ ] QuasiImpulsiveShotsAcceleration
  - [ ] Factory functions for all acceleration types
  - [ ] Add tests

- [ ] **2.4 thrust - Full Implementation**
  - [ ] All ThrustDirectionTypes
  - [ ] All ThrustMagnitudeTypes
  - [ ] ThrustDirectionSettings (all derived)
  - [ ] ThrustMagnitudeSettings (all derived)
  - [ ] EngineModelSettings
  - [ ] ThrustAccelerationSettings
  - [ ] Factory functions
  - [ ] Add tests

- [ ] **2.5 dependent_variable - Complete Coverage**
  - [ ] All PropagationDependentVariables enum values
  - [ ] All factory functions
  - [ ] Rotational state dependent variables
  - [ ] Body-fixed frame variables
  - [ ] Add tests

- [ ] **2.6 mass_rate - Full Implementation**
  - [ ] All MassRateModelType values
  - [ ] MassRateModelSettings (all derived)
  - [ ] Factory functions
  - [ ] Add tests

---

### Phase 3: Trajectory Design (MEDIUM PRIORITY)
*These enable mission design and optimization*

- [ ] **3.1 transfer_trajectory - Remaining Methods**
  - [ ] Verify all TransferLeg methods exposed
  - [ ] Verify all TransferNode methods exposed
  - [ ] Verify all TransferTrajectory methods exposed
  - [ ] Add missing specialized methods
  - [ ] Add comprehensive tests

- [ ] **3.2 shape_based_thrust - Verify Complete**
  - [ ] Audit all Python functions vs WASM
  - [ ] Add any missing functions
  - [ ] Add comprehensive tests

- [ ] **3.3 porkchop - NEW MODULE**
  - [ ] Create `/src/tudatpy_wasm/trajectory_design/porkchop/`
  - [ ] Implement Lambert grid computation
  - [ ] Implement porkchop data generation
  - [ ] Add visualization-compatible output
  - [ ] Add tests

---

### Phase 4: Astro Module (MEDIUM PRIORITY)
*These enable advanced orbital mechanics computations*

- [ ] **4.1 polyhedron_utilities - FULL IMPLEMENTATION (Currently STUB)**
  - [ ] `read_polyhedron_file()` (adapt for WASM virtual FS)
  - [ ] `compute_polyhedron_distance_to_point()`
  - [ ] `compute_polyhedron_gravitational_potential()`
  - [ ] `compute_polyhedron_gravitational_acceleration()`
  - [ ] `compute_polyhedron_laplacian()`
  - [ ] Polyhedron intersection methods
  - [ ] Surface facet utilities
  - [ ] Add tests

- [ ] **4.2 element_conversion - Verify Complete**
  - [ ] Audit all Python functions vs WASM
  - [ ] Add any missing conversions
  - [ ] Add comprehensive tests

- [ ] **4.3 frame_conversion - Verify Complete**
  - [ ] Audit all Python functions vs WASM
  - [ ] Add any missing conversions
  - [ ] Add comprehensive tests

- [ ] **4.4 gravitation - Verify Complete**
  - [ ] Audit all Python functions vs WASM
  - [ ] Add any missing functions
  - [ ] Add tests

- [ ] **4.5 time_representation - Verify Complete**
  - [ ] Audit all time conversion functions
  - [ ] Add any missing functions
  - [ ] Add tests

- [ ] **4.6 two_body_dynamics - Verify Complete**
  - [ ] Audit all two-body functions
  - [ ] Add any missing functions
  - [ ] Add tests

- [ ] **4.7 fundamentals - Verify Complete**
  - [ ] Audit all fundamental functions
  - [ ] Add any missing functions
  - [ ] Add tests

---

### Phase 5: Data Module (MEDIUM PRIORITY)
*Adapt for browser/WASM environment where possible*

- [ ] **5.1 data - Core Functions**
  - [ ] `get_resource_path()` (adapt for WASM)
  - [ ] `get_gravity_models_path()`
  - [ ] `get_spice_kernel_path()`
  - [ ] `get_earth_orientation_path()`
  - [ ] `get_ephemeris_path()`
  - [ ] `get_atmosphere_tables_path()`
  - [ ] `get_space_weather_path()`
  - [ ] Add tests

- [ ] **5.2 data - External Data (Where Feasible)**
  - [ ] Evaluate which data sources can work via fetch()
  - [ ] BatchMPC interface (via HTTP)
  - [ ] SBDBquery interface (via HTTP)
  - [ ] Horizons interface (via HTTP)
  - [ ] Document limitations for file-based data
  - [ ] Add tests

---

### Phase 6: Interface Module (LOW PRIORITY)
*SPICE and other interfaces*

- [ ] **6.1 spice - Verify Complete**
  - [ ] Audit all SPICE wrapper functions
  - [ ] Add any missing functions
  - [ ] Verify kernel loading works with embedded data
  - [ ] Add tests

---

### Phase 7: Math Module (LOW PRIORITY)
*Supporting mathematical utilities*

- [ ] **7.1 interpolators - Verify Complete**
  - [ ] Audit all interpolator types
  - [ ] Add any missing types
  - [ ] Add tests

- [ ] **7.2 root_finders - Verify Complete**
  - [ ] Audit all root finder types
  - [ ] Add any missing types
  - [ ] Add tests

- [ ] **7.3 numerical_integrators - Verify Complete**
  - [ ] Audit all integrator types
  - [ ] Add any missing types
  - [ ] Add tests

- [ ] **7.4 statistics - Verify Complete**
  - [ ] Audit all statistics functions
  - [ ] Add any missing functions
  - [ ] Add tests

- [ ] **7.5 geometry - Verify Complete**
  - [ ] Audit all geometry functions
  - [ ] Add any missing functions
  - [ ] Add tests

---

### Phase 8: Environment Setup (VERIFY)
*These are mostly implemented but need verification*

- [ ] **8.1 atmosphere - Verify Complete**
- [ ] **8.2 ephemeris - Verify Complete**
- [ ] **8.3 gravity_field - Verify Complete**
- [ ] **8.4 gravity_field_variation - Verify Complete**
- [ ] **8.5 radiation_pressure - Full Implementation**
  - [ ] PaneledRadiationPressureTargetModelSettings
  - [ ] CannonballRadiationPressureTargetModelSettings (verify)
  - [ ] Factory functions
- [ ] **8.6 aerodynamic_coefficients - Verify Complete**
- [ ] **8.7 rotation_model - Verify Complete**
- [ ] **8.8 rigid_body - Verify Complete**
- [ ] **8.9 vehicle_systems - Verify Complete**
- [ ] **8.10 ground_station - Verify Complete**
- [ ] **8.11 shape - Verify Complete**
- [ ] **8.12 shape_deformation - Verify Complete**

---

### Phase 9: Integration Tests & Examples
*Comprehensive testing to prove parity*

- [ ] **9.1 Port ALL Python Examples to WASM Tests**
  - [ ] coupled_translational_rotational_dynamics
  - [ ] juice_flybys
  - [ ] thrust_between_Earth_Moon
  - [ ] separation_satellites_diff_drag
  - [ ] aoo_custom_environment
  - [ ] aoo_design_space_exploration
  - [ ] cassini1_mga_optimization
  - [ ] hodographic_shaping_mga_optimization
  - [ ] low_thrust_earth_mars_transfer_window
  - [ ] impact_manifolds_lpo_cr3bp
  - [ ] All estimation examples (where data available)

- [ ] **9.2 Web Visualization Examples**
  - [ ] Add visualization for each new capability
  - [ ] Update visualizations/examples/ with new tests
  - [ ] Ensure all run in browser

- [ ] **9.3 API Documentation**
  - [ ] Update TypeDoc definitions for all new bindings
  - [ ] Regenerate API docs
  - [ ] Verify parity with Python API docs

---

### Phase 10: Performance & Polish
*Final optimization and cleanup*

- [ ] **10.1 WASM Size Optimization**
  - [ ] Analyze binary size
  - [ ] Remove unused code paths
  - [ ] Consider modular builds

- [ ] **10.2 Memory Management**
  - [ ] Audit for memory leaks
  - [ ] Ensure proper shared_ptr handling
  - [ ] Test long-running simulations

- [ ] **10.3 Error Handling**
  - [ ] Verify all exceptions properly converted
  - [ ] Add meaningful error messages
  - [ ] Test error conditions

---

## Progress Tracking

| Phase | Module | Status | Completion |
|-------|--------|--------|------------|
| 1.1 | estimation_analysis covariance | COMPLETE | 100% |
| 1.2 | observable_models | COMPLETE | 100% |
| 1.3 | observations | COMPLETE | 100% |
| 1.4 | observations_setup | COMPLETE | 100% |
| 1.5 | observable_models_setup | COMPLETE | 100% |
| 2.1 | propagator full | COMPLETE | 100% |
| 2.2 | torque | COMPLETE | 100% |
| 2.3 | acceleration | COMPLETE | 100% |
| 2.4 | thrust | COMPLETE | 100% |
| 2.5 | dependent_variable | COMPLETE | 100% |
| 2.6 | mass_rate | COMPLETE | 100% |
| 3.1 | transfer_trajectory | COMPLETE | 100% |
| 3.2 | shape_based_thrust | COMPLETE | 100% |
| 3.3 | porkchop | NOT STARTED | 0% |
| 4.1 | polyhedron_utilities | COMPLETE | 100% |
| 4.2-4.7 | astro submodules | COMPLETE | 100% |
| 5.1-5.2 | data | COMPLETE | 100% |
| 6.1 | spice | COMPLETE | 100% |
| 7.1-7.5 | math submodules | COMPLETE | 100% |
| 8.1-8.12 | environment_setup | COMPLETE | 100% |
| 9.1-9.3 | Integration tests | PARTIAL | 43% |
| 10.1-10.3 | Polish | IN PROGRESS | 50% |

**Note:** Module parity achieved (98.6%). Remaining work is adding more examples and the estimation_analysis_estimator submodule.

---

## Files to Create/Modify

### New Files Needed:
```
src/tudatpy_wasm/trajectory_design/porkchop/expose_porkchop_wasm.cpp
src/tudatpy_wasm/trajectory_design/porkchop/CMakeLists.txt
tests/wasm/src/testEstimationAdvanced.cpp
tests/wasm/src/testRotationalDynamics.cpp
tests/wasm/src/testPolyhedron.cpp
tests/wasm/src/testTrajectoryDesign.cpp
tests/wasm/web/visualizations/examples/[27 new example ports]
```

### Files Requiring Major Updates:
```
src/tudatpy_wasm/estimation/estimation_analysis/expose_estimation_analysis_wasm.cpp
src/tudatpy_wasm/estimation/observable_models/expose_observable_models_wasm.cpp
src/tudatpy_wasm/estimation/observations/expose_observations_wasm.cpp
src/tudatpy_wasm/estimation/observations_setup/*.cpp (all submodules)
src/tudatpy_wasm/dynamics/propagation_setup/propagator/expose_propagator_wasm.cpp
src/tudatpy_wasm/dynamics/propagation_setup/torque/expose_torque_wasm.cpp
src/tudatpy_wasm/dynamics/propagation_setup/acceleration/expose_acceleration_wasm.cpp
src/tudatpy_wasm/dynamics/propagation_setup/thrust/expose_thrust_wasm.cpp
src/tudatpy_wasm/astro/polyhedron_utilities/expose_polyhedron_utilities_wasm.cpp
src/tudatpy_wasm/data/expose_data_wasm.cpp
```

---

## Definition of Done

For each module to be considered COMPLETE:
1. All functions from Python binding are exposed in WASM
2. All classes from Python binding are exposed in WASM
3. All enums from Python binding are exposed in WASM
4. Unit tests exist and pass for all functionality
5. Integration test demonstrates real-world usage
6. Web visualization example exists (where applicable)
7. TypeDoc API documentation is updated
8. Matches Python API behavior exactly

---

## Commands

Build WASM:
```bash
cd build-wasm && ninja tudat_wasm_web
```

Run tests:
```bash
node build-wasm/tests/wasm/tudat_wasm_test.js
```

Start web server:
```bash
python tests/wasm/web/start_server.py
```

---

*Last updated: 2026-01-18*
*Current: 98.6% binding parity, 43.5% example parity*
*Target: 100% Python parity*
