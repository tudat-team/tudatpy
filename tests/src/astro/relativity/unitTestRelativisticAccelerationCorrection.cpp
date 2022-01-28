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
#include <boost/make_shared.hpp>
#include <boost/bind/bind.hpp>
using namespace boost::placeholders;


#include "tudat/math/basic/leastSquaresEstimation.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/simulation/simulation.h"
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

double computeLenseThirringPericenterPrecession(
        const double gravitationalParameter,
        const double angularMomentum,
        const double semiMajorAxis,
        const double eccentricity,
        const double inclination )
{
    return - 6.0 * gravitationalParameter * angularMomentum * std::cos(inclination  ) /
            ( physical_constants::SPEED_OF_LIGHT * physical_constants::SPEED_OF_LIGHT * semiMajorAxis * semiMajorAxis * semiMajorAxis *
              std::pow( 1.0 - eccentricity * eccentricity, 1.5 ) );
}

double computeLenseThirringNodePrecession(
        const double gravitationalParameter,
        const double angularMomentum,
        const double semiMajorAxis,
        const double eccentricity )
{
    return 2.0 * gravitationalParameter * angularMomentum /
            ( physical_constants::SPEED_OF_LIGHT * physical_constants::SPEED_OF_LIGHT * semiMajorAxis * semiMajorAxis * semiMajorAxis *
              std::pow( 1.0 - eccentricity * eccentricity, 1.5 ) );
}

double computeSchwarzschildPericenterPrecession(
        const double gravitationalParameter,
        const double semiMajorAxis,
        const double eccentricity )
{
    return 3.0 * std::pow( gravitationalParameter, 1.5 ) /
            ( physical_constants::SPEED_OF_LIGHT * physical_constants::SPEED_OF_LIGHT * std::pow(
                  semiMajorAxis, 2.5 ) *( 1.0 - eccentricity * eccentricity ) );
}

double computeDeSitterPericenterPrecession( const double meanDistanceEarthToSun,
                                            const double meanEccentricity )
{
    return 1.5 * 1.327124E20 / ( physical_constants::SPEED_OF_LIGHT * physical_constants::SPEED_OF_LIGHT * meanDistanceEarthToSun ) * 2.0 * mathematical_constants::PI /
            ( physical_constants::JULIAN_YEAR ) * std::sqrt( 1.0 - meanEccentricity * meanEccentricity );
}

BOOST_AUTO_TEST_SUITE( test_relativistic_acceleration_corrections )

void testControlPropagation(
        Eigen::Vector6d asterixInitialStateInKeplerianElements,
        std::vector< std::map< double, double > > elementMaps,
        double earthGravitationalParameter )
{
    std::vector< double > polynomialPowers = { 0, 1 };
    for( unsigned elementIndex = 0; elementIndex < 5; elementIndex++ )
    {
        std::vector< double > fitOutput = linear_algebra::getLeastSquaresPolynomialFit(
                    elementMaps[ elementIndex ], polynomialPowers );
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

void testLenseThirringPropagation(
        Eigen::Vector6d asterixInitialStateInKeplerianElements,
        std::vector< std::map< double, double > > elementMaps,
        double earthGravitationalParameter )
{
    double theoreticalLenseThirringPericenterPrecession =
            computeLenseThirringPericenterPrecession(
                earthGravitationalParameter, 1.0E9,  asterixInitialStateInKeplerianElements( semiMajorAxisIndex ),
                asterixInitialStateInKeplerianElements( eccentricityIndex ),
                asterixInitialStateInKeplerianElements( inclinationIndex ) );
    double theoreticalLenseThirringNodePrecession =
            computeLenseThirringNodePrecession(
                earthGravitationalParameter, 1.0E9,  asterixInitialStateInKeplerianElements( semiMajorAxisIndex ),
                asterixInitialStateInKeplerianElements( eccentricityIndex ) );

    std::vector< double > polynomialPowers = { 0, 1 };
    for( unsigned elementIndex = 0; elementIndex < 5; elementIndex++ )
    {
        std::vector< double > fitOutput = linear_algebra::getLeastSquaresPolynomialFit(
                    elementMaps[ elementIndex ], polynomialPowers );
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

void testSchwarzschildPropagation(
        Eigen::Vector6d asterixInitialStateInKeplerianElements,
        std::vector< std::map< double, double > > elementMaps,
        double earthGravitationalParameter )
{
    double theoreticalSchwarzschildPericenterPrecession =
            computeSchwarzschildPericenterPrecession(
                earthGravitationalParameter, asterixInitialStateInKeplerianElements( semiMajorAxisIndex ),
                asterixInitialStateInKeplerianElements( eccentricityIndex ) );

    std::vector< double > polynomialPowers = { 0, 1 };
    for( unsigned elementIndex = 0; elementIndex < 5; elementIndex++ )
    {
        std::vector< double > fitOutput = linear_algebra::getLeastSquaresPolynomialFit(
                    elementMaps[ elementIndex ], polynomialPowers );
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

void testDeSitterPropagation(
        Eigen::Vector6d asterixInitialStateInKeplerianElements,
        std::vector< std::map< double, double > > elementMaps,
        double meanDistanceEarthToSun,
        double meanEarthEccentricity )
{
    std::vector< double > polynomialPowers = { 0, 1 };
    for( unsigned elementIndex = 0; elementIndex < 5; elementIndex++ )
    {
        std::vector< double > fitOutput = linear_algebra::getLeastSquaresPolynomialFit(
                    elementMaps[ elementIndex ], polynomialPowers );
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
           BOOST_CHECK_CLOSE_FRACTION( fitOutput.at( 1 ), computeDeSitterPericenterPrecession(
                                           meanDistanceEarthToSun, meanEarthEccentricity  ), 2.5E-2 );
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

    BodyListSettings bodySettings =
            getDefaultBodySettings( bodiesToCreate );

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
        accelerationsOfAsterix[ "Earth" ].push_back( std::make_shared< AccelerationSettings >(
                                                         basic_astrodynamics::point_mass_gravity ) );
        if( testCase == 1 )
        {
            accelerationsOfAsterix[ "Earth" ].push_back( std::make_shared< RelativisticAccelerationCorrectionSettings >(
                                                             false, true, false, "", 1.0E9 * Eigen::Vector3d::UnitZ( ) ) );
        }
        if( testCase == 2 )
        {
            accelerationsOfAsterix[ "Earth" ].push_back( std::make_shared< RelativisticAccelerationCorrectionSettings >(
                                                             true, false, false ) );
        }
        if( testCase == 3 )
        {
            accelerationsOfAsterix[ "Earth" ].push_back( std::make_shared< RelativisticAccelerationCorrectionSettings >(
                                                             false, false, true, "Sun" ) );
        }
        accelerationMap[ "Asterix" ] = accelerationsOfAsterix;

        // Create acceleration models and propagation settings.
        basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodies, accelerationMap, bodiesToPropagate, centralBodies );



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
        asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
                = convertDegreesToRadians( 235.7 );
        asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
                = convertDegreesToRadians( 23.4 );
        asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 0.0 );

        // Convert Asterix state from Keplerian elements to Cartesian elements.
        double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
        Eigen::VectorXd systemInitialState = convertKeplerianToCartesianElements(
                    asterixInitialStateInKeplerianElements,
                    earthGravitationalParameter );

        std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, simulationEndEpoch, encke );


        // Create numerical integrator.
        std::shared_ptr< IntegratorSettings< > > integratorSettings =
                std::make_shared< RungeKuttaVariableStepSizeSettings< > >
                ( 0.0, 10.0,
                  RungeKuttaCoefficients::rungeKuttaFehlberg78, 1.0E-3, 1.0E3, 1.0E-12, 1.0E-12 );

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Create simulation object and propagate dynamics.
        SingleArcDynamicsSimulator< > dynamicsSimulator(
                    bodies, integratorSettings, propagatorSettings );
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
             stateIterator != integrationResult.end( ); stateIterator++ )
        {
            // Retrieve current Cartesian state (convert to Moon-centered frame if needed)
            currentCartesianState = stateIterator->second;
            keplerianIntegrationResult[ stateIterator->first ] =
                    convertCartesianToKeplerianElements(
                        currentCartesianState, earthGravitationalParameter );
            for( unsigned elementIndex = 0; elementIndex < 6; elementIndex++ )
            {
                elementMaps[ elementIndex ][ stateIterator->first ] = keplerianIntegrationResult[ stateIterator->first ]( elementIndex );
            }

            if( testCase == 3 )
            {
                earthCartesianState = spice_interface:: getBodyCartesianStateAtEpoch(
                            "Earth", "Sun", "ECLIPJ2000", "None", stateIterator->first );
                earthKeplerianState  = convertCartesianToKeplerianElements(
                            earthCartesianState, spice_interface::getBodyGravitationalParameter( "Sun" ) );

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
            testDeSitterPropagation(
                        asterixInitialStateInKeplerianElements, elementMaps, statistics::computeSampleMean( solarDistances ),
                        statistics::computeSampleMean( earthEccentricities ) );
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}

}
