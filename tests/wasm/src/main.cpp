/*    Copyright (c) 2010-2024, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Comprehensive WASM test suite for Tudat library.
 *    Includes full propagation tests without requiring external SPICE kernels.
 *    Run with: node build-wasm/tests/wasm/tudat_wasm_test.js
 */

#include <iostream>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#endif

// Test framework counters (shared across all test files)
int testsRun = 0;
int testsPassed = 0;
int testsFailed = 0;

// Forward declarations of all test functions

// Basic astrodynamics and math tests (testBasicAstro.cpp)
void testUnitConversions();
void testPhysicalConstants();
void testOrbitalElementConversions();
void testAnomalyConversions();
void testCoordinateConversions();
void testEigenOperations();
void testKeplerFunctions();
void testTimeConversions();
void testLegendrePolynomials();
void testLinearInterpolation();
void testNumericalIntegration();
void testCubicSplineInterpolation();
void testReferenceFrameTransformations();
void testModifiedEquinoctialElements();
void testStatistics();
void testSphericalHarmonics();
void testLinearAlgebra();
void testResourcePaths();
void testClohessyWiltshirePropagation();
void testKeplerPropagation();

#ifdef __EMSCRIPTEN__
void testEmscriptenEnvironment();
#endif

// Propagation tests (testPropagation.cpp)
void testCR3BPPropagation();
void testCustomStatePropagation();
void testMassPropagation();
void testTwoBodyPropagation();
void testMultiBodyMassPropagation();
void testPropagationTermination();

// SPICE tests (testSpice.cpp)
void testSpiceTimeConversions();
void testSpiceFrameRotations();
void testSpiceErrorHandling();
void testSpiceTLEPropagation();
void testSpiceTemeFrameRotation();

// Gravitation tests (testGravitation.cpp)
void testThirdBodyPerturbation();
void testLibrationPoints();
void testJacobiEnergy();
void testSphericalHarmonicsGravity();
void testCentralGravityModel();
void testDegreeTwoGravitationalTorque();
void testSphericalHarmonicGravitationalTorque();
void testInertiaFromSphericalHarmonics();

// Aerodynamics tests (testAerodynamics.cpp)
void testExponentialAtmosphere();
void testNRLMSISE00Atmosphere();
void testAerodynamicForce();
void testAerodynamicMoment();

// Mission segments tests (testMissionSegments.cpp)
void testLambertTargetingIzzo();
void testGravityAssistRoutines();
void testUnpoweredGravityAssistPropagation();
void testPoweredGravityAssistPropagation();
void testEscapeAndCapture();

// Electromagnetism tests (testElectromagnetism.cpp)
void testRadiationPressureForce();
void testRadiationPressureAccelerationEarth();
void testRadiationPressureAccelerationVenus();
void testRadiationPressureForceUranus();
void testRadiationPressureAccelerationUlysses();
void testRadiationPressureInverseSquareLaw();
void testRadiationPressureRandomPosition();
void testRadiationPressureGiancoliData();
void testLuminosityModel();

// Additional integrator tests (testIntegrators.cpp)
void testRungeKutta78Integrator();
void testRungeKutta87DormandPrinceIntegrator();
void testRungeKuttaFehlberg45Integrator();
void testBulirschStoerIntegrator();
void testAdamsBashforthMoultonIntegrator();

// Ephemerides tests (testEphemerides.cpp)
void testSimpleRotationalEphemeris();
void testKeplerEphemerisElliptical();
void testKeplerEphemerisHyperbolic();
void testTabulatedEphemeris();
void testConstantEphemeris();

// Earth orientation tests (testEarthOrientation.cpp)
void testTimeScaleConversions();
void testTimeScaleConversionPrecision();
void testLeapSecondConversions();
void testEopReaderData();
void testShortPeriodLibrationalPolarMotion();
void testShortPeriodOceanTidePolarMotion();
void testShortPeriodLibrationalUt1();
void testShortPeriodOceanTideUt1();
void testPolarMotionCalculator();
void testEarthOrientationRotationSetup();
void testHistoricalEarthRotation();
void testLeapSecondIdentification();

// Example tests ported from Python tudatpy examples (testExamples.cpp)
// Propagation examples
void testKeplerianSatelliteOrbit();
void testPerturbedSatelliteOrbit();
void testThrustWithMassPropagation();
void testCoupledTranslationalRotational();
void testDifferentialDrag();
void testSolarSystemPropagation();
void testThrustBetweenEarthMoon();
void testTwoStageRocketAscent();
void testLinearSensitivityAnalysis();
void testHybridTerminationConditions();
void testLambertTargeting();
void testVariationalEquations();
void testReentryTrajectory();
void testMultiArcPropagation();
void testCR3BPIrregularBody();
void testCustomThrustGuidance();
// Mission design examples
void testMGATrajectory();
void testPorkchopPattern();
void testLowThrustTransfer();
// Estimation examples
void testCovarianceAnalysisPattern();
void testObservationModelSetup();
void testTLEEphemeris();
void testOptimizationProblemSetup();
void testGalileanMoonsPattern();

// Estimation module tests (testEstimation.cpp)
void testStateTransitionMatrix();
void testSimpleBatchOrbitDetermination();
void testCovariancePropagation();
void testEstimationConvergenceChecker();
void testObservationTypesAndLinks();
void testFormalErrorPropagation();
void testMultiBodyEstimationSetup();

// Edge case tests (testEdgeCases.cpp)
void testNaNInfinityHandling();
void testSubnormalNumbers();
void testEpsilonComparisons();
void testCircularOrbitEdgeCase();
void testNearParabolicOrbitEdgeCase();
void testHyperbolicOrbitEdgeCase();
void testEquatorialOrbitEdgeCase();
void testPolarOrbitEdgeCase();
void testZeroTimePropagation();
void testFullOrbitPropagation();
void testVeryLongPropagation();
void testSphericalCoordinateSingularities();
void testZeroRadiusHandling();
void testIntegratorSmallStepSize();
void testIntegratorStiffODE();
void testInterpolationAtBoundaries();
void testSinglePointInterpolation();
void testSingularMatrixOperations();
void testEmptyAndZeroVectors();
void testLargeVectorOperations();

int main()
{
    std::cout << "=== Tudat WASM Test Suite ===" << std::endl;

    try {
        // Basic astrodynamics and math tests
        testUnitConversions();
        testPhysicalConstants();
        testOrbitalElementConversions();
        testAnomalyConversions();
        testCoordinateConversions();
        testEigenOperations();
        testKeplerFunctions();
        testTimeConversions();
        testLegendrePolynomials();
        testLinearInterpolation();
        testNumericalIntegration();
        testCubicSplineInterpolation();
        testReferenceFrameTransformations();
        testModifiedEquinoctialElements();
        testStatistics();
        testSphericalHarmonics();
        testLinearAlgebra();
        testResourcePaths();
        testClohessyWiltshirePropagation();
        testKeplerPropagation();

        // Propagation tests (full dynamics simulation)
        std::cout << "\n=== PROPAGATION TESTS ===" << std::endl;

        testCR3BPPropagation();           // Circular Restricted 3-Body Problem
        testCustomStatePropagation();     // Custom ODE propagation
        testMassPropagation();            // Single body mass propagation
        testTwoBodyPropagation();         // Two-body orbit propagation
        testMultiBodyMassPropagation();   // Coupled multi-body mass propagation
        testPropagationTermination();     // Termination conditions

        // SPICE tests (functions that work without external kernel files)
        std::cout << "\n=== SPICE TESTS ===" << std::endl;

        testSpiceTimeConversions();       // Julian Date <-> Ephemeris Time
        testSpiceFrameRotations();        // J2000 <-> ECLIPJ2000 rotations
        testSpiceErrorHandling();         // SPICE error control functions
        testSpiceTLEPropagation();        // SGP4 propagation with EOP files
        testSpiceTemeFrameRotation();     // TEME <-> J2000 frame rotation

        // Additional gravitation tests (ported from native)
        std::cout << "\n=== GRAVITATION TESTS ===" << std::endl;

        testThirdBodyPerturbation();      // Third-body gravitational perturbation
        testLibrationPoints();            // Lagrange point computation
        testJacobiEnergy();               // Jacobi integral of motion
        testSphericalHarmonicsGravity();  // Spherical harmonics gravity field
        testCentralGravityModel();        // Central body gravity model
        testDegreeTwoGravitationalTorque();       // Degree-2 gravitational torque
        testSphericalHarmonicGravitationalTorque(); // SH gravitational torque
        testInertiaFromSphericalHarmonics();      // Inertia tensor <-> SH coefficients

        // Aerodynamics tests
        std::cout << "\n=== AERODYNAMICS TESTS ===" << std::endl;

        testExponentialAtmosphere();      // Exponential atmosphere model
        testNRLMSISE00Atmosphere();       // NRLMSISE-00 atmosphere model
        testAerodynamicForce();           // Aerodynamic force calculation
        testAerodynamicMoment();          // Aerodynamic moment calculation

        // Mission segment tests
        std::cout << "\n=== MISSION SEGMENTS TESTS ===" << std::endl;

        testLambertTargetingIzzo();               // Izzo Lambert algorithm
        testGravityAssistRoutines();              // Gravity assist delta-V calculation
        testUnpoweredGravityAssistPropagation();  // Unpowered swing-by propagation
        testPoweredGravityAssistPropagation();    // Powered swing-by propagation
        testEscapeAndCapture();                   // Escape/capture maneuver delta-V

        // Electromagnetism tests (radiation pressure)
        std::cout << "\n=== ELECTROMAGNETISM TESTS ===" << std::endl;

        testRadiationPressureForce();                 // Cannon-ball SRP force
        testRadiationPressureAccelerationEarth();     // SRP acceleration at 1 AU
        testRadiationPressureAccelerationVenus();     // SRP acceleration at Venus
        testRadiationPressureForceUranus();           // SRP force at Uranus
        testRadiationPressureAccelerationUlysses();   // Ulysses spacecraft benchmark
        testRadiationPressureInverseSquareLaw();      // Inverse square law verification
        testRadiationPressureRandomPosition();        // Random 3D position test
        testRadiationPressureGiancoliData();          // Giancoli textbook benchmark
        testLuminosityModel();                        // Luminosity model

        // Additional integrator tests
        std::cout << "\n=== ADDITIONAL INTEGRATOR TESTS ===" << std::endl;

        testRungeKutta78Integrator();             // RKF78 adaptive integrator
        testRungeKutta87DormandPrinceIntegrator();// RKDP87 adaptive integrator
        testRungeKuttaFehlberg45Integrator();     // RKF45 adaptive integrator
        testAdamsBashforthMoultonIntegrator();    // ABM multi-step integrator
        testBulirschStoerIntegrator();            // BS integrator (fixed-step to avoid stack overflow)

        // Ephemerides tests
        std::cout << "\n=== EPHEMERIDES TESTS ===" << std::endl;

        testSimpleRotationalEphemeris();          // Venus rotational ephemeris
        testKeplerEphemerisElliptical();          // Elliptical Kepler orbit (ODTBX)
        testKeplerEphemerisHyperbolic();          // Hyperbolic Kepler orbit (GTOP)
        testTabulatedEphemeris();                 // Interpolated state ephemeris
        testConstantEphemeris();                  // Constant state ephemeris

        // Earth orientation tests
        std::cout << "\n=== EARTH ORIENTATION TESTS ===" << std::endl;

        testTimeScaleConversions();               // SOFA cookbook time scale conversions
        testTimeScaleConversionPrecision();       // High-precision time conversion roundtrip
        testLeapSecondConversions();              // UTC/TAI across leap seconds
        testEopReaderData();                      // EOP data reader and interpolation
        testShortPeriodLibrationalPolarMotion();  // Libration polar motion corrections
        testShortPeriodOceanTidePolarMotion();    // Ocean tide polar motion corrections
        testShortPeriodLibrationalUt1();          // Libration UT1 corrections
        testShortPeriodOceanTideUt1();            // Ocean tide UT1 corrections
        testPolarMotionCalculator();              // Combined polar motion calculator
        testEarthOrientationRotationSetup();      // GCRS/ITRS rotation matrices
        testHistoricalEarthRotation();            // Pre-1962 Earth orientation
        testLeapSecondIdentification();           // Leap second detection in EOP

        // Example tests ported from Python tudatpy examples
        std::cout << "\n=== PROPAGATION EXAMPLE TESTS (Ported from Python) ===" << std::endl;

        testKeplerianSatelliteOrbit();        // Basic two-body orbit propagation
        testPerturbedSatelliteOrbit();        // J2 + third body perturbations
        testThrustWithMassPropagation();      // Coupled thrust and mass propagation
        testCoupledTranslationalRotational(); // Coupled translational-rotational dynamics
        testDifferentialDrag();               // Multi-satellite propagation
        testSolarSystemPropagation();         // Multi-body planetary propagation
        testThrustBetweenEarthMoon();         // Engine thrust with mass rate model
        testTwoStageRocketAscent();           // Multi-stage rocket dynamics
        testLinearSensitivityAnalysis();      // Variational equations / STM
        testHybridTerminationConditions();    // Multiple termination conditions
        testLambertTargeting();               // Interplanetary transfer design
        testVariationalEquations();           // State transition matrix foundation
        testReentryTrajectory();              // Reentry with aerodynamic forces
        testMultiArcPropagation();            // Multi-arc propagation (JUICE flybys)
        testCR3BPIrregularBody();             // CR3BP with irregular body (impact manifolds)
        testCustomThrustGuidance();           // Custom thrust guidance (JUICE engine)

        // Mission design example tests
        std::cout << "\n=== MISSION DESIGN EXAMPLE TESTS ===" << std::endl;

        testMGATrajectory();                  // Multiple gravity assist trajectory
        testPorkchopPattern();                // Porkchop plot / launch window
        testLowThrustTransfer();              // Low-thrust transfer (hodographic shaping)

        // Estimation example tests
        std::cout << "\n=== ESTIMATION EXAMPLE TESTS ===" << std::endl;

        testCovarianceAnalysisPattern();      // Covariance analysis setup
        testObservationModelSetup();          // Ground station / observation geometry
        testTLEEphemeris();                   // TLE-based ephemeris
        testOptimizationProblemSetup();       // Optimization problem (PyGMO pattern)
        testGalileanMoonsPattern();           // Galilean moons multi-body estimation

        // Estimation module tests (comprehensive orbit determination)
        std::cout << "\n=== ESTIMATION MODULE TESTS ===" << std::endl;

        testStateTransitionMatrix();          // State transition matrix computation
        testSimpleBatchOrbitDetermination();  // Batch OD setup and parameter estimation
        testCovariancePropagation();          // Covariance propagation through dynamics
        testEstimationConvergenceChecker();   // Convergence checking functionality
        testObservationTypesAndLinks();       // Observable type definitions
        testFormalErrorPropagation();         // Formal error computation
        testMultiBodyEstimationSetup();       // Multi-body estimation (Galilean moons)

        // Edge case and boundary condition tests
        std::cout << "\n=== EDGE CASE TESTS ===" << std::endl;

        testNaNInfinityHandling();            // NaN and infinity handling
        testSubnormalNumbers();               // Subnormal/denormalized numbers
        testEpsilonComparisons();             // Machine epsilon comparisons
        testCircularOrbitEdgeCase();          // Circular orbit (e=0)
        testNearParabolicOrbitEdgeCase();     // Near-parabolic orbit (e≈1)
        testHyperbolicOrbitEdgeCase();        // Hyperbolic orbit (e>1)
        testEquatorialOrbitEdgeCase();        // Equatorial orbit (i=0)
        testPolarOrbitEdgeCase();             // Polar orbit (i=90°)
        testZeroTimePropagation();            // Zero time interval propagation
        testFullOrbitPropagation();           // Full orbital period propagation
        testVeryLongPropagation();            // Many orbital periods
        testSphericalCoordinateSingularities(); // Spherical coordinate poles
        testZeroRadiusHandling();             // Zero radius in coordinates
        testIntegratorSmallStepSize();        // Very small integrator steps
        testIntegratorStiffODE();             // Stiff differential equations
        testInterpolationAtBoundaries();      // Interpolation at data boundaries
        testSinglePointInterpolation();       // Minimal data interpolation
        testSingularMatrixOperations();       // Singular/ill-conditioned matrices
        testEmptyAndZeroVectors();            // Zero vector operations
        testLargeVectorOperations();          // Large value vector operations

#ifdef __EMSCRIPTEN__
        testEmscriptenEnvironment();
#endif
    } catch (const std::exception& e) {
        std::cerr << "\n[ERROR] Exception caught: " << e.what() << std::endl;
        return 1;
    }

    std::cout << "\n=== Test Results ===" << std::endl;
    std::cout << "[INFO] Tests run:    " << testsRun << std::endl;
    std::cout << "[INFO] Tests passed: " << testsPassed << std::endl;
    std::cout << "[INFO] Tests failed: " << testsFailed << std::endl;

    if (testsFailed > 0) {
        std::cout << "[FAIL] *** SOME TESTS FAILED ***" << std::endl;
        return 1;
    } else {
        std::cout << "[PASS] *** ALL TESTS PASSED ***" << std::endl;
        return 0;
    }
}
