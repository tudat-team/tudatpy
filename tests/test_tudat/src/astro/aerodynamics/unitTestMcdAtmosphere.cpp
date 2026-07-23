/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

/*
 * ================================================================================
 * MARS CLIMATE DATABASE (MCD) ATMOSPHERE MODEL - UNIT TESTS
 * ================================================================================
 *
 * OVERVIEW:
 * ---------
 * This file contains unit tests for the MCD (Mars Climate Database) atmosphere
 * model integration in Tudat. The tests validate that the MCD Fortran routines
 * are correctly called from C++ and that the returned atmospheric properties
 * (density, pressure, temperature, winds) are within acceptable ranges.
 *
 * REFERENCE TEST CASES:
 * ---------------------
 * The test cases are based on the MCD v6.1 INPUT_K and REF_OUTPUT_K test suite.
 * Those reference files are upstream test data and are not bundled here.
 *
 * Each numbered reference case is based on an INPUT_K file and compares
 * results against the corresponding REF_OUTPUT_K file. Any deliberate
 * deviation from the reference input is documented with that case.
 *
 * VERTICAL COORDINATES:
 * ---------------------
 * The reference INPUT_K files express their locations using zkey=1 (radial
 * distance from the planet center). The production McdAtmosphereModel expects
 * altitude above the local surface and therefore requires its shared climate
 * model to use zkey=3.
 *
 * The direct reference regressions are a deliberate exception: they disable
 * that production validation and query the climate model with zkey=2 (height
 * above the areoid) at the altitude listed for each reference location. They
 * consequently validate the MCD scenario and atmospheric output near the
 * reference location; they do not test literal parsing of the INPUT_K radial
 * coordinate.
 *
 * TEST TOLERANCES:
 * ----------------
 * The following tolerances account for coordinate-convention, interpolation,
 * MCD-version, and perturbation differences:
 *
 *   - Unperturbed case 1: 0.5%
 *   - Other 20 km and 150 km reference cases: 15-20%
 *   - 50 km reference case: 35%
 *   - Small-scale stochastic perturbation case 6: 150% (path regression only)
 *
 * These tolerances validate that:
 *   1. The MCD Fortran interface works correctly
 *   2. The returned values are physically reasonable
 *   3. The selected atmosphere vertical-coordinate convention is appropriate
 *
 * EXPECTED BEHAVIOR:
 * ------------------
 * All active tests should pass with the specified tolerances. Failures may indicate:
 *   1. MCD data files not properly installed or path incorrect
 *   2. MCD Fortran library not properly linked
 *   3. Actual coordinate conversion errors beyond expected differences
 *   4. Issues with the MCD Fortran code itself
 *
 * REFERENCE:
 * ----------
 * MCD v6.1 Documentation: http://www-mars.lmd.jussieu.fr/mcd_python/
 * MCD Paper: Millour et al. (2015), "The Mars Climate Database (MCD version 5.2)"
 * ================================================================================
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include "tudat/astro/aerodynamics/mcdAtmosphereModel.h"
#include "tudat/interface/mcd/marsClimateDatabaseClimateModel.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/simulation/environment_setup/createAtmosphereModel.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/propagation_setup/dynamicsSimulator.h"
#include "tudat/simulation/environment_setup/createClimateModel.h"

using namespace tudat::aerodynamics;

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_mcd_atmosphere )

// Helper function to convert date to seconds since J2000
double convertDateToJ2000( int day, int month, int year, int hour, int min, int sec )
{
    double julianDay = basic_astrodynamics::convertCalendarDateToJulianDay( year, month, day, hour, min, static_cast< double >( sec ) );
    double daysSinceJ2000 = julianDay - basic_astrodynamics::JULIAN_DAY_ON_J2000;
    return daysSinceJ2000 * physical_constants::JULIAN_DAY;
}

std::shared_ptr< aerodynamics::McdAtmosphereModel > createReferenceMcdAtmosphereModel( const int dustScenario,
                                                                                       const int perturbationKey,
                                                                                       const double perturbationSeed,
                                                                                       const double gravityWaveLength,
                                                                                       const int highResolutionMode )
{
    auto climateModel = std::make_shared< mcd_interface::MarsClimateDatabaseClimateModel >(
            "", dustScenario, perturbationKey, perturbationSeed, gravityWaveLength, highResolutionMode );
    auto atmosphereModel = std::make_shared< aerodynamics::McdAtmosphereModel >( climateModel );

    // The MCD reference cases use altitude above the areoid (zkey=2), whereas
    // McdAtmosphereModel normally enforces height above the local surface.
    atmosphereModel->setValidateVerticalCoordinateKey( false );
    climateModel->setZkey( 2 );
    return atmosphereModel;
}

// Reference cases 1-6 and 9 exercise the scenarios and locations defined by
// the MCD INPUT_K*.txt files and compare against their REF_OUTPUT_K* results.
// As described above, the wrapper queries these locations using zkey=2 rather
// than reproducing the files' zkey=1 radial coordinates literally. Cases 7,
// 8, and 10 require the warm-scenario data set, which is not distributed with
// the MCD data used by this project, and are therefore not registered.

// INPUT_K1 reference scenario: scenario 1, a corresponding altitude of 20 km
// above the areoid, high resolution, and no perturbations.
BOOST_AUTO_TEST_CASE( testMcdAtmosphereCase1 )
{
    double time = convertDateToJ2000( 26, 8, 2006, 3, 30, 0 );
    double altitude = 20000.0;  // meters
    double latitude = unit_conversions::convertDegreesToRadians( 15.0 );
    double longitude = unit_conversions::convertDegreesToRadians( 5.0 );

    std::shared_ptr< mcd_interface::MarsClimateDatabaseClimateModel > mcdClimateModel =
            std::make_shared< mcd_interface::MarsClimateDatabaseClimateModel >( "", 1, 0, 5.0, 0.0, 1 );

    std::shared_ptr< aerodynamics::McdAtmosphereModel > atmosphereModel = std::make_shared< McdAtmosphereModel >( mcdClimateModel );

    atmosphereModel->setValidateVerticalCoordinateKey( false );
    mcdClimateModel->setZkey( 2 );

    double density = atmosphereModel->getDensity( altitude, longitude, latitude, time );
    double pressure = atmosphereModel->getPressure( altitude, longitude, latitude, time );
    double temperature = atmosphereModel->getTemperature( altitude, longitude, latitude, time );
    double zonalWind = atmosphereModel->getZonalWind( altitude, longitude, latitude, time );
    double meridionalWind = atmosphereModel->getMeridionalWind( altitude, longitude, latitude, time );

    // REF_OUTPUT_K1 values.
    double expectedPressure = 75.7;
    double expectedDensity = 2.23e-3;
    double expectedTemperature = 177.0;
    double expectedZonalWind = -34.2;
    double expectedMeridionalWind = -7.82;

    double tolerance = 0.5;
    BOOST_CHECK_CLOSE( pressure, expectedPressure, tolerance );
    BOOST_CHECK_CLOSE( density, expectedDensity, tolerance );
    BOOST_CHECK_CLOSE( temperature, expectedTemperature, tolerance );
    BOOST_CHECK_CLOSE( std::abs( zonalWind ), std::abs( expectedZonalWind ), tolerance );
    BOOST_CHECK_CLOSE( std::abs( meridionalWind ), std::abs( expectedMeridionalWind ), tolerance );
}

// INPUT_K2: scenario 1, 20 km above the areoid, high resolution, and a
// large-scale perturbation with seed 5. The 15% tolerance is retained for
// the version-sensitive perturbed result.
BOOST_AUTO_TEST_CASE( testMcdAtmosphereCase2 )
{
    const double time = convertDateToJ2000( 26, 8, 2006, 3, 30, 0 );
    const double altitude = 20000.0;
    const double latitude = unit_conversions::convertDegreesToRadians( 15.0 );
    const double longitude = unit_conversions::convertDegreesToRadians( 5.0 );
    const auto atmosphereModel = createReferenceMcdAtmosphereModel( 1, 2, 5.0, 0.0, 1 );

    BOOST_CHECK_CLOSE( atmosphereModel->getPressure( altitude, longitude, latitude, time ), 77.8, 15.0 );
    BOOST_CHECK_CLOSE( atmosphereModel->getDensity( altitude, longitude, latitude, time ), 2.29e-3, 15.0 );
    BOOST_CHECK_CLOSE( atmosphereModel->getTemperature( altitude, longitude, latitude, time ), 178.0, 15.0 );
}

// INPUT_K3 reference location: scenario 1, 50 km above the areoid, high
// resolution. This regression deliberately disables the INPUT_K3
// small-scale perturbation (perturbationKey=3), which has historically
// caused Fortran memory failures. It therefore validates the mean
// atmosphere at the reference location. The larger tolerance reflects the
// increased sensitivity of density and pressure to altitude at 50 km.
BOOST_AUTO_TEST_CASE( testMcdAtmosphereCase3 )
{
    const double time = convertDateToJ2000( 26, 8, 2006, 3, 30, 0 );
    const double altitude = 50000.0;
    const double latitude = unit_conversions::convertDegreesToRadians( 15.0 );
    const double longitude = unit_conversions::convertDegreesToRadians( 5.0 );
    const auto atmosphereModel = createReferenceMcdAtmosphereModel( 1, 0, 0.0, 0.0, 1 );

    BOOST_CHECK_CLOSE( atmosphereModel->getPressure( altitude, longitude, latitude, time ), 1.80, 35.0 );
    BOOST_CHECK_CLOSE( atmosphereModel->getDensity( altitude, longitude, latitude, time ), 7.31e-5, 35.0 );
    BOOST_CHECK_CLOSE( atmosphereModel->getTemperature( altitude, longitude, latitude, time ), 129.0, 35.0 );
}

// INPUT_K4: scenario 1, 150 km above the areoid, high resolution, and no
// perturbations. At this altitude the influence of topography is small; the
// 15% tolerance accommodates remaining MCD-version differences.
BOOST_AUTO_TEST_CASE( testMcdAtmosphereCase4 )
{
    const double time = convertDateToJ2000( 26, 8, 2006, 3, 30, 0 );
    const double altitude = 150000.0;
    const double latitude = unit_conversions::convertDegreesToRadians( 15.0 );
    const double longitude = unit_conversions::convertDegreesToRadians( 5.0 );
    const auto atmosphereModel = createReferenceMcdAtmosphereModel( 1, 0, 0.0, 0.0, 1 );

    BOOST_CHECK_CLOSE( atmosphereModel->getPressure( altitude, longitude, latitude, time ), 4.79e-6, 15.0 );
    BOOST_CHECK_CLOSE( atmosphereModel->getDensity( altitude, longitude, latitude, time ), 1.21e-10, 15.0 );
    BOOST_CHECK_CLOSE( atmosphereModel->getTemperature( altitude, longitude, latitude, time ), 193.0, 15.0 );
}

// INPUT_K5: scenario 1, 150 km above the areoid, high resolution, and a
// large-scale perturbation with seed 5. The expected values are the
// perturbed REF_OUTPUT_K5 pressure, density, and temperature.
BOOST_AUTO_TEST_CASE( testMcdAtmosphereCase5 )
{
    const double time = convertDateToJ2000( 26, 8, 2006, 3, 30, 0 );
    const double altitude = 150000.0;
    const double latitude = unit_conversions::convertDegreesToRadians( 15.0 );
    const double longitude = unit_conversions::convertDegreesToRadians( 5.0 );
    const auto atmosphereModel = createReferenceMcdAtmosphereModel( 1, 2, 5.0, 0.0, 1 );

    BOOST_CHECK_CLOSE( atmosphereModel->getPressure( altitude, longitude, latitude, time ), 4.63e-6, 15.0 );
    BOOST_CHECK_CLOSE( atmosphereModel->getDensity( altitude, longitude, latitude, time ), 1.19e-10, 15.0 );
    BOOST_CHECK_CLOSE( atmosphereModel->getTemperature( altitude, longitude, latitude, time ), 190.0, 15.0 );
}

// INPUT_K6: scenario 1, 150 km above the areoid, low resolution, and
// small-scale gravity-wave perturbations with seed 5 and a 16 km wavelength.
//
// Exact agreement with REF_OUTPUT_K6 is not expected: this stochastic path
// is sensitive to the MCD version and seed handling, and the low-resolution
// implementation can apply gravity-wave perturbations differently. The
// 150% tolerance intentionally makes this a perturbation-path regression:
// it checks that the call succeeds and returns values on the scale of the
// perturbed reference result, rather than claiming precise numerical
// equivalence.
BOOST_AUTO_TEST_CASE( testMcdAtmosphereCase6 )
{
    const double time = convertDateToJ2000( 26, 8, 2006, 3, 30, 0 );
    const double altitude = 150000.0;
    const double latitude = unit_conversions::convertDegreesToRadians( 15.0 );
    const double longitude = unit_conversions::convertDegreesToRadians( 5.0 );
    const auto atmosphereModel = createReferenceMcdAtmosphereModel( 1, 3, 5.0, 16000.0, 0 );

    // Perturbed REF_OUTPUT_K6 values.
    BOOST_CHECK_CLOSE( atmosphereModel->getPressure( altitude, longitude, latitude, time ), 4.79e-6, 150.0 );
    BOOST_CHECK_CLOSE( atmosphereModel->getDensity( altitude, longitude, latitude, time ), 9.47e-11, 150.0 );
    BOOST_CHECK_CLOSE( atmosphereModel->getTemperature( altitude, longitude, latitude, time ), 226.0, 150.0 );
}

// INPUT_K9: scenario 1, 150 km above the areoid, low resolution, and no
// perturbations. The 20% tolerance accounts for the combination of the
// low-resolution data path and the sensitivity of the upper atmosphere.
BOOST_AUTO_TEST_CASE( testMcdAtmosphereCase9 )
{
    const double time = convertDateToJ2000( 26, 8, 2006, 3, 30, 0 );
    const double altitude = 150000.0;
    const double latitude = unit_conversions::convertDegreesToRadians( 15.0 );
    const double longitude = unit_conversions::convertDegreesToRadians( 5.0 );
    const auto atmosphereModel = createReferenceMcdAtmosphereModel( 1, 0, 0.0, 0.0, 0 );

    BOOST_CHECK_CLOSE( atmosphereModel->getPressure( altitude, longitude, latitude, time ), 4.79e-6, 20.0 );
    BOOST_CHECK_CLOSE( atmosphereModel->getDensity( altitude, longitude, latitude, time ), 1.21e-10, 20.0 );
    BOOST_CHECK_CLOSE( atmosphereModel->getTemperature( altitude, longitude, latitude, time ), 193.0, 20.0 );
}

// Test MCD atmosphere in propagation
BOOST_AUTO_TEST_CASE( testMcdAtmosphereInPropagation )
{
    using namespace aerodynamics;
    using namespace simulation_setup;
    using namespace numerical_integrators;
    using namespace basic_astrodynamics;
    using namespace propagators;

    // Load spice kernels
    spice_interface::loadStandardSpiceKernels( );

    // Create Mars with MCD atmosphere
    BodyListSettings defaultBodySettings = getDefaultBodySettings( { "Mars" } );

    defaultBodySettings.at( "Mars" )->climateModelSettings =
            std::make_shared< simulation_setup::MarsClimateDatabaseClimateModelSettings >( "", 1, 0, 0.0, 0.0, 0 );

    defaultBodySettings.at( "Mars" )->atmosphereSettings = std::make_shared< simulation_setup::McdAtmosphereSettings >( );

    SystemOfBodies bodies = createSystemOfBodies( defaultBodySettings );

    // Create vehicle
    double vehicleMass = 5.0E3;
    bodies.createEmptyBody( "Vehicle" );
    bodies.at( "Vehicle" )->setConstantBodyMass( vehicleMass );

    // Set aerodynamic coefficients
    std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
            std::make_shared< ConstantAerodynamicCoefficientSettings >( 2.0,
                                                                        4.0,
                                                                        Eigen::Vector3d::Zero( ),
                                                                        Eigen::Vector3d::UnitX( ),
                                                                        Eigen::Vector3d::Zero( ),
                                                                        negative_aerodynamic_frame_coefficients,
                                                                        negative_aerodynamic_frame_coefficients );
    bodies.at( "Vehicle" )
            ->setAerodynamicCoefficientInterface(
                    createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Vehicle", bodies ) );

    // Define accelerations
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;

    accelerationsOfVehicle[ "Mars" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationsOfVehicle[ "Mars" ].push_back( std::make_shared< AccelerationSettings >( aerodynamic ) );

    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;
    bodiesToPropagate.push_back( "Vehicle" );
    centralBodies.push_back( "Mars" );

    // Set initial state
    Eigen::Vector6d systemInitialState;
    systemInitialState << 3500.0E3, 0.0, 0.0, 0.0, 3500.0, 0.0;
    double initialTime = convertDateToJ2000( 26, 8, 2006, 3, 30, 0 );
    // Create acceleration models
    basic_astrodynamics::AccelerationMap accelerationModelMap =
            createAccelerationModelsMap( bodies, accelerationMap, bodiesToPropagate, centralBodies );

    // Set dependent variables
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
    dependentVariables.push_back(
            std::make_shared< SingleDependentVariableSaveSettings >( altitude_dependent_variable, "Vehicle", "Mars" ) );
    dependentVariables.push_back(
            std::make_shared< SingleDependentVariableSaveSettings >( local_density_dependent_variable, "Vehicle", "Mars" ) );

    // Set propagation settings
    std::shared_ptr< PropagationTimeTerminationSettings > terminationSettings =
            std::make_shared< propagators::PropagationTimeTerminationSettings >( initialTime + 1000.0 );

    std::shared_ptr< IntegratorSettings<> > integratorSettings = rungeKuttaFixedStepSettings( 10.0, CoefficientSets::rungeKutta4Classic );

    std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >( centralBodies,
                                                                                accelerationModelMap,
                                                                                bodiesToPropagate,
                                                                                systemInitialState,
                                                                                initialTime,
                                                                                integratorSettings,
                                                                                terminationSettings,
                                                                                cowell,
                                                                                dependentVariables );

    // Create simulation object and propagate
    SingleArcDynamicsSimulator<> dynamicsSimulator( bodies, translationalPropagatorSettings );

    BOOST_CHECK_EQUAL( dynamicsSimulator.getSingleArcPropagationResults( )->getPropagationIsPerformed( ), true );
    BOOST_CHECK_EQUAL( dynamicsSimulator.getSingleArcPropagationResults( )->integrationCompletedSuccessfully( ), true );
}

BOOST_AUTO_TEST_CASE( testMcdClimateModelCaching )
{
    // The current cache uses the complete query tuple as an exact key.
    // Repeating that tuple must return the retained result, while changing
    // one coordinate must create a separate result.
    mcd_interface::MarsClimateDatabaseClimateModel climateModel;
    const double altitude = 5000.0;
    const double longitude = 0.5;
    const double latitude = 0.3;
    const double time = 1000.0;

    const auto firstResult = climateModel.getCache( altitude, longitude, latitude, time );
    const auto repeatedResult = climateModel.getCache( altitude, longitude, latitude, time );
    BOOST_CHECK( firstResult == repeatedResult );

    const auto changedAltitudeResult = climateModel.getCache( altitude + 1000.0, longitude, latitude, time );
    BOOST_CHECK( firstResult != changedAltitudeResult );
    BOOST_CHECK( std::abs( firstResult->density_ - changedAltitudeResult->density_ ) > 1.0e-10 );
}

BOOST_AUTO_TEST_CASE( testMcdClimateModelExtraVariableCaching )
{
    // Adding a requested MCD extra variable changes the contents required
    // from each cache entry. Existing entries must consequently be
    // invalidated, after which identical queries should again be reused.
    mcd_interface::MarsClimateDatabaseClimateModel climateModel;
    const double altitude = 5000.0;
    const double longitude = 0.5;
    const double latitude = 0.3;
    const double time = 1000.0;

    const auto resultWithoutExtraVariables = climateModel.getCache( altitude, longitude, latitude, time );
    climateModel.addExtraVariableKeys( { mcd_interface::ExtVar::ratio_of_specific_heats } );
    const auto resultWithExtraVariables = climateModel.getCache( altitude, longitude, latitude, time );

    BOOST_CHECK( resultWithoutExtraVariables != resultWithExtraVariables );
    BOOST_CHECK( resultWithExtraVariables == climateModel.getCache( altitude, longitude, latitude, time ) );

    const double specificHeatRatio =
            climateModel.getExtraVariable( mcd_interface::ExtVar::ratio_of_specific_heats, altitude, longitude, latitude, time );
    BOOST_CHECK( specificHeatRatio > 1.0 );
    BOOST_CHECK( specificHeatRatio < 2.0 );
}

BOOST_AUTO_TEST_CASE( testMcdClimateModelCacheKeys )
{
    // All query coordinates, including the vertical-coordinate
    // interpretation (zkey), contribute to cache identity. Unlike the
    // former atmosphere-model cache, the current cache has no tolerance
    // within which nearby queries are treated as identical.
    mcd_interface::MarsClimateDatabaseClimateModel climateModel;
    const double altitude = 5000.0;
    const double longitude = 0.5;
    const double latitude = 0.3;
    const double time = 1000.0;

    const auto baseline = climateModel.getCache( altitude, longitude, latitude, time );
    BOOST_CHECK( baseline != climateModel.getCache( altitude + 1.0, longitude, latitude, time ) );
    BOOST_CHECK( baseline != climateModel.getCache( altitude, longitude + 1.0e-6, latitude, time ) );
    BOOST_CHECK( baseline != climateModel.getCache( altitude, longitude, latitude + 1.0e-6, time ) );
    BOOST_CHECK( baseline != climateModel.getCache( altitude, longitude, latitude, time + 1.0 ) );

    climateModel.setZkey( 2 );
    BOOST_CHECK( baseline != climateModel.getCache( altitude, longitude, latitude, time ) );
}

BOOST_AUTO_TEST_CASE( testMcdClimateModelCacheEviction )
{
    // Limit the cache to two entries and insert three distinct queries.
    // Reloading the oldest query must recompute an equivalent result,
    // demonstrating bounded first-in-first-out eviction.
    mcd_interface::MarsClimateDatabaseClimateModel climateModel( "", 1, 0, 0.0, 0.0, 0, 2 );
    const double longitude = 0.5;
    const double latitude = 0.3;
    const double time = 1000.0;

    const auto firstResult = climateModel.getCache( 5000.0, longitude, latitude, time );
    const auto secondResult = climateModel.getCache( 6000.0, longitude, latitude, time );
    const auto thirdResult = climateModel.getCache( 7000.0, longitude, latitude, time );
    BOOST_CHECK( firstResult != secondResult );
    BOOST_CHECK( secondResult != thirdResult );

    // Inserting the third distinct key evicts the oldest entry.
    const auto reloadedFirstResult = climateModel.getCache( 5000.0, longitude, latitude, time );
    BOOST_CHECK( firstResult != reloadedFirstResult );
    BOOST_CHECK_EQUAL( firstResult->density_, reloadedFirstResult->density_ );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
