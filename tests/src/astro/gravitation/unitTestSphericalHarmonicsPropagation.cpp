/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Mathworks. gravitysphericalharmonic, Implement spherical harmonic representation of
 *        planetary gravity, help documentation of Aerospace Toolbox of MATLAB R2012a, 2012.
 *
 *    Notes
 *      In future, more tests should be added here to test the completeness of the functions
 *      implemented. In particular, tests should be added to ascertain the maximum degree and order
 *      to which the functions are still able to produce accelerations. Further, the tests are
 *      currently all based off of data generated with MATLAB (Mathworks, 2012). Ideally, at least
 *      one other source of benchmark data should be included to thoroughly test the code and
 *      minimize the risks of bugs being present. The runtime errors thrown are also not tested;
 *      this behaviour needs to be tested rigorously too.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "tudat/io/readHistoryFromFile.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"
#include "tudat/simulation/simulation.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"


namespace tudat
{
namespace unit_tests
{

using namespace tudat::orbit_determination;
using namespace tudat::estimatable_parameters;
using namespace tudat::interpolators;
using namespace tudat::numerical_integrators;
using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;
using namespace tudat::orbital_element_conversions;
using namespace tudat::ephemerides;
using namespace tudat::propagators;
using namespace tudat::basic_astrodynamics;
using namespace tudat::coordinate_conversions;
using namespace tudat::physical_constants;

BOOST_AUTO_TEST_SUITE( test_SphericalHarmonicsGravityPropagation )

// Check single harmonics term of degree = 2 and order = 0.
BOOST_AUTO_TEST_CASE( testSphericalHarmonicsGravityPropagation )
{
// Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Define bodies in simulation
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Moon" );
    bodyNames.push_back( "Sun" );

    // Specify initial time
    double initialEphemerisTime = 1.0E8;
    double finalEphemerisTime = initialEphemerisTime + 10000;

    // Create bodies needed in simulation
    BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, "Earth", "ECLIPJ2000" );

    bodySettings.get( "Earth" )->gravityFieldVariationSettings.push_back(
        fixedSingleDegreeLoveNumberGravityFieldVariationSettings( "Moon", 50.0, 2 ) );
    bodySettings.get( "Earth" )->gravityFieldVariationSettings.push_back(
        fixedSingleDegreeLoveNumberGravityFieldVariationSettings( "Sun", 100.0, 2 ) );

    std::shared_ptr< SphericalHarmonicsGravityFieldSettings > moonGravityField =
        std::dynamic_pointer_cast< SphericalHarmonicsGravityFieldSettings >( bodySettings.at( "Moon" )->gravityFieldSettings );

    Eigen::MatrixXd cosineCoefficients = moonGravityField->getCosineCoefficients( );
    cosineCoefficients *= 1000.0;
    cosineCoefficients( 0, 0 ) = 1.0;
    moonGravityField->resetCosineCoefficients( cosineCoefficients );

    Eigen::MatrixXd sineCoefficients = moonGravityField->getSineCoefficients( );
    sineCoefficients *= 1000.0;
    moonGravityField->resetSineCoefficients( sineCoefficients );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    bodies.createEmptyBody( "Vehicle" );

    // Set accelerations on Vehicle that are to be taken into account.
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 12, 12 ) );
    accelerationsOfVehicle[ "Moon" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 2, 2 ) );
    accelerationsOfVehicle[ "Sun" ].push_back( pointMassGravityAcceleration( ) );

    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;

    // Set bodies for which initial state is to be estimated and integrated.
    std::vector< std::string > bodiesToIntegrate;
    std::vector< std::string > centralBodies;
    bodiesToIntegrate.push_back( "Vehicle" );
    centralBodies.push_back( "Earth" );

    // Create acceleration models
    AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodies, accelerationMap, bodiesToIntegrate, centralBodies );

    // Set Keplerian elements for Asterix.
    Eigen::Vector6d asterixInitialStateInKeplerianElements;
    asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7200.0E3;
    asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.05;
    asterixInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
    asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = unit_conversions::convertDegreesToRadians( 235.7 );
    asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = unit_conversions::convertDegreesToRadians( 23.4 );
    asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );

    double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );

    // Set (perturbed) initial state.
    Eigen::Matrix< double, 6, 1 > systemInitialState =
        convertKeplerianToCartesianElements( asterixInitialStateInKeplerianElements, earthGravitationalParameter );

    // Get dependent variables

    std::vector< std::pair< int, int > > earthComponentIndices;
    for( int i = 0; i <= 12; i++ )
    {
        for( int j = 0; ( j <= i && j <= 12 ); j++ )
        {
            earthComponentIndices.push_back( std::make_pair( i, j ) );
        }
    }

    std::vector< std::pair< int, int > > moonComponentIndices;
    std::vector< std::pair< int, int > > tideComponentIndices;
    for( int i = 0; i <= 2; i++ )
    {
        for( int j = 0; ( j <= i && j <= 2 ); j++ )
        {
            moonComponentIndices.push_back( std::make_pair( i, j ) );
            if( i == 2 )
            {
                tideComponentIndices.push_back( std::make_pair( i, j ) );
            }
        }
    }



    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
    dependentVariables.push_back( singleAccelerationDependentVariable( basic_astrodynamics::spherical_harmonic_gravity, "Vehicle", "Earth" ) );
    dependentVariables.push_back( singleAccelerationDependentVariable( basic_astrodynamics::spherical_harmonic_gravity, "Vehicle", "Moon" ) );
    dependentVariables.push_back( sphericalHarmonicAccelerationTermsDependentVariable( "Vehicle", "Earth", earthComponentIndices ) );
    dependentVariables.push_back( sphericalHarmonicAccelerationTermsDependentVariable( "Vehicle", "Moon", moonComponentIndices ) );
    dependentVariables.push_back( singleGravityFieldVariationSeparateTermsAccelerationContributionVariable(
        "Vehicle", "Earth", tideComponentIndices, gravitation::basic_solid_body, "Sun" ) );
    dependentVariables.push_back( singleGravityFieldVariationSeparateTermsAccelerationContributionVariable(
        "Vehicle", "Earth", tideComponentIndices, gravitation::basic_solid_body, "Moon" ) );

    // Create propagator settings
    std::shared_ptr< IntegratorSettings< double > > integratorSettings =
        rungeKuttaFixedStepSettings( 40.0, CoefficientSets::rungeKuttaFehlberg78 );
    std::shared_ptr< TranslationalStatePropagatorSettings< double, double > > propagatorSettings =
        std::make_shared< TranslationalStatePropagatorSettings< double, double > >(
            centralBodies, accelerationModelMap, bodiesToIntegrate, systemInitialState, initialEphemerisTime, integratorSettings,
            propagationTimeTerminationSettings( finalEphemerisTime ), cowell, dependentVariables );

    propagatorSettings->getOutputSettings( )->setResultsSaveFrequencyInSeconds( 10000.0 );
    propagatorSettings->getOutputSettings( )->setResultsSaveFrequencyInSteps( 0 );
    propagatorSettings->getOutputSettings()->getPrintSettings()->setPrintDependentVariableData(true);
    // Define parameters
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames = getInitialStateParameterSettings< double, double >(
        propagatorSettings, bodies );
    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
        2, 0, 12, 12, "Earth", spherical_harmonics_cosine_coefficient_block ) );
    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
        2, 1, 12, 12, "Earth", spherical_harmonics_sine_coefficient_block ) );
    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
        2, 0, 2, 2, "Earth", spherical_harmonics_cosine_coefficient_block ) );
    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
        2, 1, 2, 2, "Earth", spherical_harmonics_sine_coefficient_block ) );

    // Create parameters
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
        createParametersToEstimate< double, double >( parameterNames, bodies );

    printEstimatableParameterEntries( parametersToEstimate );

    std::shared_ptr< SingleArcVariationalEquationsSolver< double, double > > variationalEquationsSolver =
        std::make_shared< SingleArcVariationalEquationsSolver< double, double > >( bodies, propagatorSettings, parametersToEstimate, true );

    std::map< double, Eigen::MatrixXd > stateTransitionMatrixHistory = variationalEquationsSolver->getStateTransitionMatrixSolution( );
//    input_output::writeDataMapToTextFile( stateTransitionMatrixHistory, "sphericalHarmonicsTestStateTransition.dat", tudat::paths::getTudatTestDataPath( ) + "", "", 16 );

    std::map< double, Eigen::MatrixXd > referenceStateTransitionMatrixHistory =
    input_output::readMatrixHistoryFromFile< double, double >( stateTransitionMatrixHistory.begin( )->second.rows( ),
                                                 stateTransitionMatrixHistory.begin( )->second.cols( ),
                                                 tudat::paths::getTudatTestDataPath( ) + "sphericalHarmonicsTestStateTransition.dat" );
    for( auto it : stateTransitionMatrixHistory )
    {
        BOOST_CHECK_EQUAL( referenceStateTransitionMatrixHistory.count( it.first ), 1 );
        Eigen::MatrixXd currentMatrix = stateTransitionMatrixHistory.at( it.first );
        Eigen::MatrixXd matrixDifference = stateTransitionMatrixHistory.at( it.first ) - referenceStateTransitionMatrixHistory.at( it.first );

        for( int i = 0; i < 3; i++ )
        {
            for( int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL( matrixDifference( i, j ), currentMatrix.block( 0, 0, 3, 3 ).norm( ) * 1.0E-14 );
                BOOST_CHECK_SMALL( matrixDifference( i + 3, j ), currentMatrix.block( 3, 0, 3, 3 ).norm( ) * 1.0E-14 );
                BOOST_CHECK_SMALL( matrixDifference( i, j + 3 ), currentMatrix.block( 0, 3, 3, 3 ).norm( ) * 1.0E-14 );
                BOOST_CHECK_SMALL( matrixDifference( i + 3, j + 3 ), currentMatrix.block( 3, 3, 3, 3 ).norm( ) * 1.0E-14 );

            }
        }
    }

    std::map< double, Eigen::MatrixXd > sensitivityMatrixHistory = variationalEquationsSolver->getSensitivityMatrixSolution( );
//    input_output::writeDataMapToTextFile( sensitivityMatrixHistory, "sphericalHarmonicsTestSensitivity.dat", tudat::paths::getTudatTestDataPath( ) + "", "", 16 );

    std::map< double, Eigen::MatrixXd > referenceSensitivityMatrixHistory =
        input_output::readMatrixHistoryFromFile< double, double >( sensitivityMatrixHistory.begin( )->second.rows( ),
                                                                   sensitivityMatrixHistory.begin( )->second.cols( ),
                                                                   tudat::paths::getTudatTestDataPath( ) + "sphericalHarmonicsTestSensitivity.dat" );
    for( auto it : sensitivityMatrixHistory )
    {
        BOOST_CHECK_EQUAL( referenceSensitivityMatrixHistory.count( it.first ), 1 );
        Eigen::MatrixXd currentMatrix = sensitivityMatrixHistory.at( it.first );
        Eigen::MatrixXd matrixDifference = sensitivityMatrixHistory.at( it.first ) - referenceSensitivityMatrixHistory.at( it.first );
        for( int i = 0; i < 3; i++ )
        {
            for( int j = 0; j < currentMatrix.cols( ); j++ )
            {
                BOOST_CHECK_SMALL( matrixDifference( i, j ), currentMatrix.block( 0, j, 3, 3 ).norm( ) * 1.0E-12 );
                BOOST_CHECK_SMALL( matrixDifference( i + 3, j ), currentMatrix.block( 3, j, 3, 3 ).norm( ) * 1.0E-12 );
            }
        }
    }

    std::map< double, Eigen::VectorXd > stateHistory = variationalEquationsSolver->getEquationsOfMotionSolution( );
//    input_output::writeDataMapToTextFile( stateHistory, "sphericalHarmonicsTestStates.dat", tudat::paths::getTudatTestDataPath( ) + "", "", 16 );

    std::map< double, Eigen::MatrixXd > referenceStateHistory =
        input_output::readMatrixHistoryFromFile< double, double >( stateHistory.begin( )->second.rows( ),
                                                                   stateHistory.begin( )->second.cols( ),
                                                                   tudat::paths::getTudatTestDataPath( ) + "sphericalHarmonicsTestStates.dat" );
    for( auto it : stateHistory )
    {
        BOOST_CHECK_EQUAL( referenceStateHistory.count( it.first ), 1 );
        Eigen::VectorXd currentMatrix = stateHistory.at( it.first );
        Eigen::VectorXd matrixDifference = stateHistory.at( it.first ) - referenceStateHistory.at( it.first );
        for( int i = 0; i < 3; i++ )
        {
            BOOST_CHECK_SMALL( matrixDifference( i ), currentMatrix.segment( 0, 3 ).norm( ) * 1.0E-14 );
            BOOST_CHECK_SMALL( matrixDifference( i + 3 ), currentMatrix.segment( 3, 3 ).norm( ) * 1.0E-14 );
        }
    }

    std::map< double, Eigen::VectorXd > dependentVariableHistory =
        variationalEquationsSolver->getSingleArcVariationalPropagationResults( )->getDynamicsResults( )->getDependentVariableHistoryDouble( );
//    input_output::writeDataMapToTextFile( dependentVariableHistory, "sphericalHarmonicsTestDependentVariables.dat", tudat::paths::getTudatTestDataPath( ) + "", "", 16 );

    std::map< double, Eigen::MatrixXd > referenceDependentVariableHistory =
        input_output::readMatrixHistoryFromFile< double, double >( dependentVariableHistory.begin( )->second.rows( ),
                                                                   dependentVariableHistory.begin( )->second.cols( ),
                                                                   tudat::paths::getTudatTestDataPath( ) + "sphericalHarmonicsTestDependentVariables.dat" );
    for( auto it : dependentVariableHistory )
    {
        BOOST_CHECK_EQUAL( referenceDependentVariableHistory.count( it.first ), 1 );
        Eigen::VectorXd currentMatrix = dependentVariableHistory.at( it.first );
        Eigen::VectorXd matrixDifference = dependentVariableHistory.at( it.first ) - referenceDependentVariableHistory.at( it.first );

        for( int i = 0; i < currentMatrix.rows( ) / 3; i++ )
        {
            BOOST_CHECK_SMALL( matrixDifference( 3 * i ), currentMatrix.segment( 3 * i, 3 ).norm( ) * 1.0E-12 );
            BOOST_CHECK_SMALL( matrixDifference( 3 * i + 1), currentMatrix.segment( 3 * i, 3 ).norm( ) * 1.0E-12 );
            BOOST_CHECK_SMALL( matrixDifference( 3 * i + 2), currentMatrix.segment( 3 * i, 3 ).norm( ) * 1.0E-12 );

        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
