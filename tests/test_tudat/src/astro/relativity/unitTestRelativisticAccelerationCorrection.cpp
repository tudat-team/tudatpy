/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/createMetric.h"
#include "tudat/simulation/propagation_setup/singleArcDynamicsSimulator.h"
#include "tudat/astro/relativity/relativisticEquationsOfMotion.h"
#include "tudat/astro/relativity/relativisticAccelerationCorrection.h"

#include "tudat/math/basic/leastSquaresEstimation.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/math/statistics/basicStatistics.h"

namespace tudat
{

namespace unit_tests
{

using namespace tudat::simulation_setup;
using namespace tudat::propagators;
using namespace tudat::numerical_integrators;
using namespace tudat::orbital_element_conversions;
using namespace tudat::basic_mathematics;
using namespace tudat::unit_conversions;

double computeLenseThirringPericenterPrecession( const double gravitationalParameter,
                                                 const double angularMomentum,
                                                 const double semiMajorAxis,
                                                 const double eccentricity,
                                                 const double inclination )
{
    return -6.0 * gravitationalParameter * angularMomentum * std::cos( inclination ) /
            ( physical_constants::SPEED_OF_LIGHT * physical_constants::SPEED_OF_LIGHT * semiMajorAxis * semiMajorAxis * semiMajorAxis *
              std::pow( 1.0 - eccentricity * eccentricity, 1.5 ) );
}

double computeLenseThirringNodePrecession( const double gravitationalParameter,
                                           const double angularMomentum,
                                           const double semiMajorAxis,
                                           const double eccentricity )
{
    return 2.0 * gravitationalParameter * angularMomentum /
            ( physical_constants::SPEED_OF_LIGHT * physical_constants::SPEED_OF_LIGHT * semiMajorAxis * semiMajorAxis * semiMajorAxis *
              std::pow( 1.0 - eccentricity * eccentricity, 1.5 ) );
}

double computeSchwarzschildPericenterPrecession( const double gravitationalParameter,
                                                 const double semiMajorAxis,
                                                 const double eccentricity )
{
    return 3.0 * std::pow( gravitationalParameter, 1.5 ) /
            ( physical_constants::SPEED_OF_LIGHT * physical_constants::SPEED_OF_LIGHT * std::pow( semiMajorAxis, 2.5 ) *
              ( 1.0 - eccentricity * eccentricity ) );
}

double computeDeSitterPericenterPrecession( const double meanDistanceEarthToSun, const double meanEccentricity )
{
    return 1.5 * 1.327124E20 / ( physical_constants::SPEED_OF_LIGHT * physical_constants::SPEED_OF_LIGHT * meanDistanceEarthToSun ) * 2.0 *
            mathematical_constants::PI / (physical_constants::JULIAN_YEAR)*std::sqrt( 1.0 - meanEccentricity * meanEccentricity );
}

std::map< double, Eigen::Vector3d > runPropagationAndRetrieveTotalAcceleration(
        const SelectedAccelerationMap& selectedAccelerationSettings,
        const double simulationStartEpoch,
        const double simulationEndEpoch,
        const std::vector< std::string >& bodiesToPropagate,
        const std::vector< std::string >& centralBodies,
        const Eigen::Vector6d& vehicleInitialStateInKeplerianElements )
{
    BodyListSettings bodySettings = getDefaultBodySettings( { "Earth" } );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    bodies.createEmptyBody( "Vehicle" );

    createBaseMetric( schwardschildSpaceTimeMetricSettings( "Earth" ), bodies );

    basic_astrodynamics::AccelerationMap accelerationModelMap =
            createAccelerationModelsMap( bodies, selectedAccelerationSettings, bodiesToPropagate, centralBodies );

    const double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
    const Eigen::VectorXd systemInitialState =
            convertKeplerianToCartesianElements( vehicleInitialStateInKeplerianElements, earthGravitationalParameter );

    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
    dependentVariables.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
            total_acceleration_dependent_variable, "Vehicle" ) );

    std::shared_ptr< IntegratorSettings<> > integratorSettings =
            std::make_shared< IntegratorSettings<> >( rungeKutta4, simulationStartEpoch, 20.0 );

    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >(
                    centralBodies,
                    accelerationModelMap,
                    bodiesToPropagate,
                    systemInitialState,
                    simulationStartEpoch,
                    integratorSettings,
                    std::make_shared< PropagationTimeTerminationSettings >( simulationEndEpoch ),
                    cowell,
                    dependentVariables );

    SingleArcDynamicsSimulator<> dynamicsSimulator( bodies, propagatorSettings );
    const std::map< double, Eigen::VectorXd > dependentVariableHistory = dynamicsSimulator.getDependentVariableHistory( );

    std::map< double, Eigen::Vector3d > accelerationHistory;
    for( const auto& dependentVariableEntry : dependentVariableHistory )
    {
        if( dependentVariableEntry.second.rows( ) != 3 )
        {
            throw std::runtime_error(
                    "Error in runPropagationAndRetrieveTotalAcceleration: expected 3 components for total acceleration dependent variable." );
        }
        accelerationHistory[ dependentVariableEntry.first ] = dependentVariableEntry.second.segment( 0, 3 );
    }

    return accelerationHistory;
}

BOOST_AUTO_TEST_SUITE( test_relativistic_acceleration_corrections )

void testControlPropagation( Eigen::Vector6d asterixInitialStateInKeplerianElements,
                             std::vector< std::map< double, double > > elementMaps,
                             double earthGravitationalParameter )
{
    std::vector< double > polynomialPowers = { 0, 1 };
    for( unsigned elementIndex = 0; elementIndex < 5; elementIndex++ )
    {
        std::vector< double > fitOutput = linear_algebra::getLeastSquaresPolynomialFit( elementMaps[ elementIndex ], polynomialPowers );
        BOOST_CHECK_CLOSE_FRACTION( asterixInitialStateInKeplerianElements( elementIndex ), fitOutput.at( 0 ), 1.0E-10 );
        if( elementIndex == 1 )
        {
            BOOST_CHECK_SMALL( fitOutput.at( 1 ), 1.0E-18 );
        }
        else
        {
            BOOST_CHECK_SMALL( fitOutput.at( 1 ), 1.0E-12 );
        }
    }
}

void testLenseThirringPropagation( Eigen::Vector6d asterixInitialStateInKeplerianElements,
                                   std::vector< std::map< double, double > > elementMaps,
                                   double earthGravitationalParameter )
{
    double theoreticalLenseThirringPericenterPrecession =
            computeLenseThirringPericenterPrecession( earthGravitationalParameter,
                                                      1.0E9,
                                                      asterixInitialStateInKeplerianElements( semiMajorAxisIndex ),
                                                      asterixInitialStateInKeplerianElements( eccentricityIndex ),
                                                      asterixInitialStateInKeplerianElements( inclinationIndex ) );
    double theoreticalLenseThirringNodePrecession =
            computeLenseThirringNodePrecession( earthGravitationalParameter,
                                                1.0E9,
                                                asterixInitialStateInKeplerianElements( semiMajorAxisIndex ),
                                                asterixInitialStateInKeplerianElements( eccentricityIndex ) );

    std::vector< double > polynomialPowers = { 0, 1 };
    for( unsigned elementIndex = 0; elementIndex < 5; elementIndex++ )
    {
        std::vector< double > fitOutput = linear_algebra::getLeastSquaresPolynomialFit( elementMaps[ elementIndex ], polynomialPowers );
        BOOST_CHECK_CLOSE_FRACTION( asterixInitialStateInKeplerianElements( elementIndex ), fitOutput.at( 0 ), 1.0E-10 );
        if( elementIndex == 1 )
        {
            BOOST_CHECK_SMALL( fitOutput.at( 1 ), 1.0E-18 );
        }
        else if( elementIndex == 3 )
        {
            BOOST_CHECK_CLOSE_FRACTION( fitOutput.at( 1 ), theoreticalLenseThirringPericenterPrecession, 1.0E-5 );
        }
        else if( elementIndex == 4 )
        {
            BOOST_CHECK_CLOSE_FRACTION( fitOutput.at( 1 ), theoreticalLenseThirringNodePrecession, 1.0E-5 );
        }
        else
        {
            BOOST_CHECK_SMALL( fitOutput.at( 1 ), 1.0E-12 );
        }
    }
}

void testSchwarzschildPropagation( Eigen::Vector6d asterixInitialStateInKeplerianElements,
                                   std::vector< std::map< double, double > > elementMaps,
                                   double earthGravitationalParameter )
{
    double theoreticalSchwarzschildPericenterPrecession =
            computeSchwarzschildPericenterPrecession( earthGravitationalParameter,
                                                      asterixInitialStateInKeplerianElements( semiMajorAxisIndex ),
                                                      asterixInitialStateInKeplerianElements( eccentricityIndex ) );

    std::vector< double > polynomialPowers = { 0, 1 };
    for( unsigned elementIndex = 0; elementIndex < 5; elementIndex++ )
    {
        std::vector< double > fitOutput = linear_algebra::getLeastSquaresPolynomialFit( elementMaps[ elementIndex ], polynomialPowers );
        if( elementIndex != 1 )
        {
            BOOST_CHECK_CLOSE_FRACTION( asterixInitialStateInKeplerianElements( elementIndex ), fitOutput.at( 0 ), 1.0E-8 );
        }
        else
        {
            BOOST_CHECK_CLOSE_FRACTION( asterixInitialStateInKeplerianElements( elementIndex ), fitOutput.at( 0 ), 1.0E-7 );
        }
        if( elementIndex == 1 )
        {
            BOOST_CHECK_SMALL( fitOutput.at( 1 ), 1.0E-16 );
        }
        else if( elementIndex == 3 )
        {
            BOOST_CHECK_CLOSE_FRACTION( fitOutput.at( 1 ), theoreticalSchwarzschildPericenterPrecession, 1.0E-5 );
        }
        else
        {
            BOOST_CHECK_SMALL( fitOutput.at( 1 ), 1.0E-10 );
        }
    }
}

void testDeSitterPropagation( Eigen::Vector6d asterixInitialStateInKeplerianElements,
                              std::vector< std::map< double, double > > elementMaps,
                              double meanDistanceEarthToSun,
                              double meanEarthEccentricity )
{
    std::vector< double > polynomialPowers = { 0, 1 };
    for( unsigned elementIndex = 0; elementIndex < 5; elementIndex++ )
    {
        std::vector< double > fitOutput = linear_algebra::getLeastSquaresPolynomialFit( elementMaps[ elementIndex ], polynomialPowers );
        if( elementIndex != 4 )
        {
            BOOST_CHECK_CLOSE_FRACTION( asterixInitialStateInKeplerianElements( elementIndex ), fitOutput.at( 0 ), 1.0E-10 );
        }
        else
        {
            BOOST_CHECK_CLOSE_FRACTION( asterixInitialStateInKeplerianElements( elementIndex ), fitOutput.at( 0 ), 1.0E-8 );
        }
        if( elementIndex == 1 )
        {
            BOOST_CHECK_SMALL( fitOutput.at( 1 ), 1.0E-18 );
        }
        else if( elementIndex == 4 )
        {
            BOOST_CHECK_CLOSE_FRACTION(
                    fitOutput.at( 1 ), computeDeSitterPericenterPrecession( meanDistanceEarthToSun, meanEarthEccentricity ), 2.5E-2 );
        }
        else
        {
            BOOST_CHECK_SMALL( fitOutput.at( 1 ), 1.0E-12 );
        }
    }
}
BOOST_AUTO_TEST_CASE( testLenseThirring )
{
    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation end epoch.
    const double simulationEndEpoch = 0.25 * tudat::physical_constants::JULIAN_YEAR;

    // Create body objects.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Sun" );

    BodyListSettings bodySettings = getDefaultBodySettings( bodiesToCreate );

    // Create Earth object
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // Create spacecraft object.
    bodies.createEmptyBody( "Asterix" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    for( unsigned int testCase = 0; testCase < 4; testCase++ )
    {
        // Define propagator settings variables.
        SelectedAccelerationMap accelerationMap;
        std::vector< std::string > bodiesToPropagate;
        std::vector< std::string > centralBodies;

        bodiesToPropagate.push_back( "Asterix" );
        centralBodies.push_back( "Earth" );

        // Define propagation settings.
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfAsterix;
        accelerationsOfAsterix[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );
        if( testCase == 1 )
        {
            accelerationsOfAsterix[ "Earth" ].push_back( std::make_shared< RelativisticAccelerationCorrectionSettings >(
                    false, true, false, "", 1.0E9 * Eigen::Vector3d::UnitZ( ) ) );
        }
        if( testCase == 2 )
        {
            accelerationsOfAsterix[ "Earth" ].push_back(
                    std::make_shared< RelativisticAccelerationCorrectionSettings >( true, false, false ) );
        }
        if( testCase == 3 )
        {
            accelerationsOfAsterix[ "Earth" ].push_back(
                    std::make_shared< RelativisticAccelerationCorrectionSettings >( false, false, true, "Sun" ) );
        }
        accelerationMap[ "Asterix" ] = accelerationsOfAsterix;

        // Create acceleration models and propagation settings.
        basic_astrodynamics::AccelerationMap accelerationModelMap =
                createAccelerationModelsMap( bodies, accelerationMap, bodiesToPropagate, centralBodies );

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Set initial conditions for the Asterix satellite that will be propagated in this simulation.
        // The initial conditions are given in Keplerian elements and later on converted to Cartesian
        // elements.

        // Set Keplerian elements for Asterix.
        Eigen::Vector6d asterixInitialStateInKeplerianElements;
        asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 5000.0E3;
        asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.2;
        asterixInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 65.3 );
        asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 235.7 );
        asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 23.4 );
        asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 0.0 );

        // Convert Asterix state from Keplerian elements to Cartesian elements.
        double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
        Eigen::VectorXd systemInitialState =
                convertKeplerianToCartesianElements( asterixInitialStateInKeplerianElements, earthGravitationalParameter );

        std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >(
                        centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, simulationEndEpoch, encke );

        // Create numerical integrator.
        std::shared_ptr< IntegratorSettings<> > integratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettings<> >(
                0.0, 10.0, rungeKuttaFehlberg78, 1.0E-3, 1.0E3, 1.0E-12, 1.0E-12 );

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Create simulation object and propagate dynamics.
        SingleArcDynamicsSimulator<> dynamicsSimulator( bodies, integratorSettings, propagatorSettings );
        std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
        std::map< double, Eigen::VectorXd > keplerianIntegrationResult;

        // Compute map of Kepler elements
        Eigen::Vector6d currentCartesianState;
        std::vector< std::map< double, double > > elementMaps;
        elementMaps.resize( 6 );

        std::vector< double > solarDistances;
        std::vector< double > earthSemiMajorAxes;
        std::vector< double > earthEccentricities;

        Eigen::Vector6d earthKeplerianState;
        Eigen::Vector6d earthCartesianState;

        for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = integrationResult.begin( );
             stateIterator != integrationResult.end( );
             stateIterator++ )
        {
            // Retrieve current Cartesian state (convert to Moon-centered frame if needed)
            currentCartesianState = stateIterator->second;
            keplerianIntegrationResult[ stateIterator->first ] =
                    convertCartesianToKeplerianElements( currentCartesianState, earthGravitationalParameter );
            for( unsigned elementIndex = 0; elementIndex < 6; elementIndex++ )
            {
                elementMaps[ elementIndex ][ stateIterator->first ] = keplerianIntegrationResult[ stateIterator->first ]( elementIndex );
            }

            if( testCase == 3 )
            {
                earthCartesianState =
                        spice_interface::getBodyCartesianStateAtEpoch( "Earth", "Sun", "ECLIPJ2000", "None", stateIterator->first );
                earthKeplerianState =
                        convertCartesianToKeplerianElements( earthCartesianState, spice_interface::getBodyGravitationalParameter( "Sun" ) );

                earthSemiMajorAxes.push_back( earthKeplerianState( 0 ) );
                earthEccentricities.push_back( earthKeplerianState( 1 ) );

                solarDistances.push_back( earthCartesianState.segment( 0, 3 ).norm( ) );
            }
        }

        if( testCase == 0 )
        {
            testControlPropagation( asterixInitialStateInKeplerianElements, elementMaps, earthGravitationalParameter );
        }
        else if( testCase == 1 )
        {
            testLenseThirringPropagation( asterixInitialStateInKeplerianElements, elementMaps, earthGravitationalParameter );
        }
        else if( testCase == 2 )
        {
            testSchwarzschildPropagation( asterixInitialStateInKeplerianElements, elementMaps, earthGravitationalParameter );
        }
        else if( testCase == 3 )
        {
            testDeSitterPropagation( asterixInitialStateInKeplerianElements,
                                     elementMaps,
                                     statistics::computeSampleMean( solarDistances ),
                                     statistics::computeSampleMean( earthEccentricities ) );
        }
    }
}

BOOST_AUTO_TEST_CASE( testDirectRelativisticAccelerationFromMetricAgainstSchwarzschildCorrection )
{
    spice_interface::loadStandardSpiceKernels( );

    const double simulationStartEpoch = 0.0;
    const double simulationEndEpoch = simulationStartEpoch + 20.0;

    std::vector< std::string > bodiesToPropagate{ "Vehicle" };
    std::vector< std::string > centralBodies{ "Earth" };

    Eigen::Vector6d vehicleInitialStateInKeplerianElements;
    vehicleInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7000.0E3;
    vehicleInitialStateInKeplerianElements( eccentricityIndex ) = 0.05;
    vehicleInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 30.0 );
    vehicleInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 40.0 );
    vehicleInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 15.0 );
    vehicleInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 10.0 );

    SelectedAccelerationMap correctionAccelerationSettings;
    correctionAccelerationSettings[ "Vehicle" ][ "Earth" ] = {
            std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ),
            std::make_shared< RelativisticAccelerationCorrectionSettings >( true, false, false ) };

    std::map< double, Eigen::Vector3d > correctionTotalAccelerationHistory =
            runPropagationAndRetrieveTotalAcceleration(
                    correctionAccelerationSettings,
                    simulationStartEpoch,
                    simulationEndEpoch,
                    bodiesToPropagate,
                    centralBodies,
                    vehicleInitialStateInKeplerianElements );

    SelectedAccelerationMap metricAccelerationSettings;
    metricAccelerationSettings[ "Vehicle" ][ "" ] = {
            relativisticAccelerationFromMetric( ) };

    std::map< double, Eigen::Vector3d > metricTotalAccelerationHistory =
            runPropagationAndRetrieveTotalAcceleration(
                    metricAccelerationSettings,
                    simulationStartEpoch,
                    simulationEndEpoch,
                    bodiesToPropagate,
                    centralBodies,
                    vehicleInitialStateInKeplerianElements );

    BOOST_REQUIRE_EQUAL( correctionTotalAccelerationHistory.size( ), metricTotalAccelerationHistory.size( ) );

    double maximumRelativeDifference = 0.0;
    auto correctionIterator = correctionTotalAccelerationHistory.begin( );
    auto metricIterator = metricTotalAccelerationHistory.begin( );
    for( ; correctionIterator != correctionTotalAccelerationHistory.end( ); ++correctionIterator, ++metricIterator )
    {
        BOOST_CHECK_SMALL( correctionIterator->first - metricIterator->first, 1.0E-14 );

        const Eigen::Vector3d correctionAcceleration = correctionIterator->second;
        const Eigen::Vector3d metricAcceleration = metricIterator->second;
        const double correctionNorm = correctionAcceleration.norm( );
        const double absoluteDifferenceNorm = ( metricAcceleration - correctionAcceleration ).norm( );
        const double relativeDifference =
                absoluteDifferenceNorm / std::max( correctionNorm, 1.0E-30 );

        if( relativeDifference > maximumRelativeDifference )
        {
            maximumRelativeDifference = relativeDifference;
        }
    }

    BOOST_CHECK_SMALL( maximumRelativeDifference, 1.0E-14 );
}

BOOST_AUTO_TEST_CASE( testDirectMetricSchwarzschildRelativisticPartAgainstCorrectionAtGenericState )
{
    spice_interface::loadStandardSpiceKernels( );

    BodyListSettings bodySettings = getDefaultBodySettings( { "Earth" } );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    createBaseMetric( schwardschildSpaceTimeMetricSettings( "Earth" ), bodies );

    std::shared_ptr< relativity::Metric > metric = bodies.getSpaceTimeProperties( )->getBaseMetric( );
    BOOST_REQUIRE( metric != nullptr );

    const double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
    const std::shared_ptr< relativity::PPNParameterSet > ppnParameterSet =
            bodies.getSpaceTimeProperties( )->getPpnParameterSet( );
    BOOST_REQUIRE( ppnParameterSet != nullptr );

    Eigen::Vector6d testState = Eigen::Vector6d::Zero( );
    testState << 7.0e6, -2.0e6, 1.0e6, 1200.0, 7300.0, -900.0;

    const double testEpoch = 0.0;
    const Eigen::Vector3d directTotalAcceleration =
            relativity::evaluateRelativisticEquationsOfMotionInCoordinateTimeWithUpdate( metric, testState, testEpoch );

    const Eigen::Vector3d relativePosition = testState.segment( 0, 3 );
    const Eigen::Vector3d relativeVelocity = testState.segment( 3, 3 );
    const double relativeDistance = relativePosition.norm( );

    const Eigen::Vector3d newtonianPointMassAcceleration =
            -earthGravitationalParameter * relativePosition / std::pow( relativeDistance, 3.0 );
    const Eigen::Vector3d directRelativisticPart = directTotalAcceleration - newtonianPointMassAcceleration;

    const double commonCorrectionTerm =
            relativity::calculateRelativisticAccelerationCorrectionsCommonterm(
                    earthGravitationalParameter, relativeDistance );
    const Eigen::Vector3d schwarzschildCorrection =
            relativity::calculateScharzschildGravitationalAccelerationCorrection(
                    earthGravitationalParameter,
                    relativePosition,
                    relativeVelocity,
                    relativeDistance,
                    commonCorrectionTerm,
                    ppnParameterSet->getParameterGamma( ),
                    ppnParameterSet->getParameterBeta( ) );

    const double absoluteDifference = ( directRelativisticPart - schwarzschildCorrection ).norm( );
    const double relativeDifference =
            absoluteDifference /
            std::max( schwarzschildCorrection.norm( ), 1.0e-30 );

    BOOST_CHECK_SMALL( absoluteDifference, 1.0e-14 );
    BOOST_CHECK_SMALL( relativeDifference, 1.0e-6 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
