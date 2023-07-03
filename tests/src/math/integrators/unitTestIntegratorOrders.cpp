/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <cmath>
#include <iostream>
#include <limits>

#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "tudat/basics/utilityMacros.h"
#include "tudat/math/integrators/createNumericalIntegrator.h"
#include "tudat/simulation/environment_setup.h"
#include "tudat/simulation/propagation_setup.h"

using namespace tudat;
using namespace tudat::numerical_integrators;
using namespace tudat::simulation_setup;
using namespace tudat::basic_astrodynamics;
using namespace tudat::orbital_element_conversions;
using namespace tudat::propagators;
using namespace tudat::numerical_integrators;

namespace tudat
{
namespace unit_tests
{
BOOST_AUTO_TEST_SUITE( test_numerical_integrator_orders )

Eigen::Vector6d getFinalIntegrationError(
    const std::shared_ptr< IntegratorSettings< > > integratorSettings,
    const double numberOfOrbits,
    const double eccentricity)
{
    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Define bodies in simulation.
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );

    // Create bodies needed in simulation
    BodyListSettings bodySettings =
        getDefaultBodySettings(
            bodyNames, "Earth", "ECLIPJ2000" );
    SystemOfBodies bodies = createSystemOfBodies< double, double >( bodySettings );
    bodies.createEmptyBody( "Vehicle" );


    // Propagate the moon only
    std::vector< std::string > bodiesToIntegrate;
    std::vector< std::string > centralBodies;
    bodiesToIntegrate.push_back( "Vehicle" );
    centralBodies.push_back( "Earth" );

    // Set accelerations between bodies that are to be taken into account.
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;
    AccelerationMap accelerationModelMap = createAccelerationModelsMap(
        bodies, accelerationMap, bodiesToIntegrate, centralBodies );

    Eigen::Vector6d initialKeplerElements = Eigen::Vector6d::Zero( );
    initialKeplerElements[ semiMajorAxisIndex ] = 7.50E6;
    initialKeplerElements[ eccentricityIndex ] = eccentricity;
    initialKeplerElements[ inclinationIndex ] = 85.3 * mathematical_constants::PI / 180.0;
    initialKeplerElements[ argumentOfPeriapsisIndex ] = 235.7 * mathematical_constants::PI / 180.0;
    initialKeplerElements[ longitudeOfAscendingNodeIndex ] = 23.4 * mathematical_constants::PI / 180.0;
    initialKeplerElements[ trueAnomalyIndex ] = 139.87 * mathematical_constants::PI / 180.0;

    // Specify initial time
    double initialEphemerisTime = 0.0;
    double finalEphemerisTime = numberOfOrbits * mathematical_constants::PI * std::sqrt(
        std::pow( initialKeplerElements[ semiMajorAxisIndex ], 3.0 ) /
            bodies.getBody( "Earth" )->getGravitationalParameter( ) );

    Eigen::VectorXd systemInitialState = convertKeplerianToCartesianElements(
        initialKeplerElements, bodies.getBody( "Earth" )->getGravitationalParameter( ) );

    std::shared_ptr<TranslationalStatePropagatorSettings<double, double> > propagatorSettings =
        std::make_shared<TranslationalStatePropagatorSettings<double, double> >
        ( centralBodies, accelerationModelMap, bodiesToIntegrate, systemInitialState,
            initialEphemerisTime, integratorSettings,
            std::make_shared<PropagationTimeTerminationSettings>( finalEphemerisTime, true ) );

    // Create dynamics simulation object.
    SingleArcDynamicsSimulator<double, double> dynamicsSimulator(
        bodies, propagatorSettings );

    return ( dynamicsSimulator.getEquationsOfMotionNumericalSolution( ).rbegin( )->second - systemInitialState );


}

//! This tests the order of the fixed-step Runge-Kutta methods, for those with a scheme that only
//! permits fixed time steps.
//!
//! The method test computes the order p from a set of subsequent runs, and then computes the mean (to
//! verify whether the method behaves as would be expected) and standard deviation (to check the order
//! behaviour is robust)
//!
//! The exact values of the time steps have been tweaked to make sure the solutions fall in the range
//! where the time step is expected to behave well.
BOOST_AUTO_TEST_CASE( testPureFixedMultiStageNumericalIntegratorOrder )
{
    // List of schemes for which to test
    std::vector< CoefficientSets > coefficients = {
        rungeKutta4Classic,
        explicitMidPoint,
        explicitTrapezoidRule,
        ralston,
        rungeKutta3,
        ralston3,
        SSPRK3,
        ralston4,
        threeEighthRuleRK4 };

    // Multiplication of nominal time step for each method
    std::vector< double > timeStepMultiplications = {
        1.0, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 1.0, 1.0 };

    // Ideal order of each method
    std::vector< double > expectedOrders = {
        4.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 4.0, 4.0 };

    // Calculate final state error for each method
    for( unsigned int j = 0; j < coefficients.size( ); j++ )
    {
        std::vector< double > timeSteps = { 30.0, 60.0, 120.0, 240.0 };
        std::vector< double > errors;
        std::vector< double > calculatedOrder;

        for ( unsigned int i = 0; i < timeSteps.size( ); i++ )
        {
            Eigen::VectorXd finalState = getFinalIntegrationError(
                rungeKuttaFixedStepSettings(
                    timeStepMultiplications.at( j ) * timeSteps.at( i ), coefficients.at( j ) ), 10.0, 0.01 );
            errors.push_back( finalState.segment( 0, 3 ).norm( ));

            // Calculate the order p from two subsequent errors, using the fact that
            // error_{i} / error_{i-1} = ( step_{i} / step_{i-1} )^p
            if ( i > 0 )
            {
                calculatedOrder.push_back( std::log2( errors.at( i ) / errors.at( i - 1 )));
            }
        }

        // Calculate mean and standard deviation of values of order computed from eacg two subsequent steps
        double meanOrder = std::accumulate( calculatedOrder.begin( ), calculatedOrder.end( ), 0.0 ) / calculatedOrder.size( );
        double standardDeviationOrder = std::sqrt(
            std::inner_product( calculatedOrder.begin( ), calculatedOrder.end( ), calculatedOrder.begin( ), 0.0 ) /
            calculatedOrder.size( ) - meanOrder * meanOrder );

        // Check that the value of the order is correct, and reasonably constant
        // Note that, in some cases, the order computed from the subsequent time steps is
        // somewhat higher than theoretically expected from purely looking at the order of the method
        // This indicates that higher-order effects play a substantial role in the final errors.
        BOOST_CHECK( meanOrder > expectedOrders.at( j ) - 0.01 );
        BOOST_CHECK( meanOrder < expectedOrders.at( j ) + 1.0 );
        BOOST_CHECK( standardDeviationOrder < 0.2 );

    }
}

//! This tests the order of the fixed-step Runge-Kutta methods, for those with a scheme that also
//! permits variable time steps
//!
//! The method test computes the order p from a set of subsequent runs, and then computes the mean (to
//! verify whether the method behaves as would be expected) and standard deviation (to check the order
//! behaviour is robust)
//!
//! The exact values of the time steps have been tweaked to make sure the solutions fall in the range
//! where the time step is expected to behave well.
BOOST_AUTO_TEST_CASE( testFixedMultiStageNumericalIntegratorOrder )
{
    std::vector< CoefficientSets > coefficients = {
//        heunEuler,
//        rungeKuttaFehlberg12,
        rungeKuttaFehlberg45,
        rungeKuttaFehlberg56,
        rungeKuttaFehlberg78,
        rungeKutta87DormandPrince,
        rungeKuttaFehlberg89,
        rungeKuttaVerner89,
        rungeKuttaFeagin108,
        rungeKuttaFeagin1210,
        rungeKuttaFeagin1412 };

    // Multiplication of nominal time step for each method
    std::vector< double > timeStepMultiplications = {
//        0.25, 0.25,
        0.25, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.1, 1.5 };

    // Ideal order of each method (lower)
    std::vector< double > expectedLowerOrders = {
//        1.0, 1.0,
        4.0, 5.0, 7.0, 7.0, 8.0, 8.0, 8.0, 10.0, 12.0 };

    // Ideal order of each method (higher)
    std::vector< double > expectedHigherOrders = {
//        2.0, 2.0,
        5.0, 6.0, 8.0, 8.0, 9.0, 9.0, 10.0, 12.0, 14.0 };

    // Time steps that are taken form methods
    std::vector<double> timeSteps = { 316.22776602, 237.13737057, 177.827941, 133.35214322, 100. };
    std::vector<double> timeStepsHigh = { 316.22776602, 273.84196343, 237.13737057, 205.35250265, 177.827941};

    // Run test for higher and lower coefficient set
    for( unsigned int k = 0; k < 2; k++ )
    {
        // Calculate final state error for each method
        for ( unsigned int j = 0; j < coefficients.size( ); j++ )
        {
            std::vector<double> timeStepsToUse = ( j < coefficients.size( ) - 2 ) ? timeSteps : timeStepsHigh;
            std::vector<double> errors;
            std::vector<double> calculatedOrder;

            for ( unsigned int i = 0; i < timeStepsToUse.size( ); i++ )
            {
                Eigen::VectorXd finalState = getFinalIntegrationError(
                    rungeKuttaFixedStepSettings(
                        timeStepMultiplications.at( j ) * timeStepsToUse.at( i ), coefficients.at( j ),
                        static_cast< RungeKuttaCoefficients::OrderEstimateToIntegrate >( k ) ), 10.0, 0.01 );
                errors.push_back( finalState.segment( 0, 3 ).norm( ));

                // Calculate the order p from two subsequent errors, using the fact that
                // error_{i} / error_{i-1} = ( step_{i} / step_{i-1} )^p
                // NOTE, the first run on the 10(8), 12(10) and 14(12) method is omitted due to error stability issues
                if ( i > 0 && !( ( j > 5 ) && ( i == 1 ) ) )
                {
                    calculatedOrder.push_back( std::log( errors.at( i ) / errors.at( i - 1 )) /
                                            std::log( timeStepsToUse.at( i ) / timeStepsToUse.at( i - 1 )));
                }
            }

            // Calculate mean and standard deviation of values of order computed from eacg two subsequent steps
            double meanOrder = std::accumulate( calculatedOrder.begin( ), calculatedOrder.end( ), 0.0 ) / calculatedOrder.size( );
            double standardDeviationOrder = std::sqrt(
                std::inner_product( calculatedOrder.begin( ), calculatedOrder.end( ), calculatedOrder.begin( ), 0.0 ) /
                calculatedOrder.size( ) - meanOrder * meanOrder );

            // Check that the value of the order is correct, and reasonably constant
            // Note that, in some cases, the order computed from the subsequent time steps is
            // somewhat higher than theoretically expected from purely looking at the order of the method
            // This indicates that higher-order effects play a substantial role in the final errors.
            //
            if( k == 0 )
            {

                BOOST_CHECK( meanOrder > expectedLowerOrders.at( j ) - 0.5 );
                BOOST_CHECK( meanOrder < expectedLowerOrders.at( j ) + ( j == 8 ) ? 2.0 : 1.0 );
                BOOST_CHECK( standardDeviationOrder < 0.4 );
            }
            else
            {
                BOOST_CHECK( meanOrder > expectedHigherOrders.at( j ) - 0.5 );
                BOOST_CHECK( meanOrder < expectedHigherOrders.at( j ) + 1.0 );
                BOOST_CHECK( standardDeviationOrder < ( j == 6 ) ? 1.0 : 0.4 );

            }
        }
    }
}


//! This tests the order of the fixed-step Bulirsch-Stoer integrtor, using an unperturbed Earth orbiter.
//! For a method using k substeps, the order should ideally be 2*k + 1.
//! The method test computes the order p from a set of subsequent runs, and then computes the mean (to
//! verify whether the method behaves as would be expected) and standard deviation (to check the order
//! behaviour is robust)
//! The exact values of the time steps have been tweaked to make sure the solutions fall in the range
//! where the time step is expected to behave well. For high-order methods in particular, this has required
//! substantial tuning (for a 10th order method, a factor 2 in time step will result in a factor 1024 in error)
BOOST_AUTO_TEST_CASE( testFixedBulirschStoerNumericalIntegratorOrder )
{
    std::vector<double> timeSteps = { 400.0, 350.0, 300.0, 250.0, 200.0 };

    // Run test for different sequences
    for( unsigned int k = 0; k < 2; k++ )
    {
        // Run test for different number of substeps
        for ( unsigned int j = 3; j < 7; j++ )
        {
            // Use higher time step for highest order method
            double multiplicationFactor = ( j == 6 ) ? 1.8 : 1.0;

            // Calculate final state error for each method
            std::vector<double> errors;
            std::vector<double> calculatedOrder;
            for ( unsigned int i = 0; i < timeSteps.size( ); i++ )
            {
                Eigen::VectorXd finalState = getFinalIntegrationError(
                    bulirschStoerFixedStepIntegratorSettings(
                        multiplicationFactor * timeSteps.at( i ), static_cast< ExtrapolationMethodStepSequences >( k ), j ), 10.0, 0.01 );
                errors.push_back( finalState.segment( 0, 3 ).norm( ));

                // Calculate the order p from two subsequent errors, using the fact that
                // error_{i} / error_{i-1} = ( step_{i} / step_{i-1} )^p
                if ( i > 0 )
                {
                    calculatedOrder.push_back( std::log( errors.at( i ) / errors.at( i - 1 ) ) /
                                            std::log( timeSteps.at( i ) / timeSteps.at( i - 1 ) ) );
                }
            }

            // Calculate mean and standard deviation of values of order computed from eacg two subsequent steps
            double meanOrder = std::accumulate( calculatedOrder.begin( ), calculatedOrder.end( ), 0.0 ) / calculatedOrder.size( );
            double standardDeviationOrder = std::sqrt(
                std::inner_product( calculatedOrder.begin( ), calculatedOrder.end( ), calculatedOrder.begin( ), 0.0 ) /
                calculatedOrder.size( ) - meanOrder * meanOrder );

            std::cout<<j<<" "<<k<<" "<<meanOrder<<" "<<standardDeviationOrder<<std::endl;

            //TODO: reinstate testst
            // Check that the order is in the right range.
            //BOOST_CHECK( meanOrder > ( 2.0 * static_cast< double >( j ) ) );
            //BOOST_CHECK( meanOrder < ( 2.0 * static_cast< double >( j ) ) + 2 );

            // Check that the value of the order is reasonable constant; the high value for the test here is only needed for
            // the 13th order method.
            //BOOST_CHECK( standardDeviationOrder < 1.0 );
        }
    }
}
//
//BOOST_AUTO_TEST_CASE( testFixedAbmNumericalIntegratorOrder )
//{
//
//    std::vector<double> multiplicationFactors =
//        { 1.0, 1.0,
//        1.0, 1.0, 2.0, 2.0, 4.0, 6.0, 10.0, 12.0 };
//
//    // Run test for different number of substeps
//    for ( unsigned int j = 3; j < 10; j++ )
//    {    std::vector<double> timeSteps = { 2.0, 4.0, 8.0, 16.0 };
//
//        // Use higher time step for highest order method
//        double multiplicationFactor = multiplicationFactors.at( j );// ( j == 6 ) ? 1.8 : 1.0;
//
//        // Calculate final state error for each method
//        std::vector<double> errors;
//        std::vector<double> calculatedOrder;
//        for ( unsigned int i = 0; i < timeSteps.size( ); i++ )
//        {
//            Eigen::VectorXd finalState = getFinalIntegrationError(
//                adamsBashforthMoultonSettingsFixedStepFixedOrder(
//                    multiplicationFactor * timeSteps.at( i ), j ), 10.0, 0.01 );
//            errors.push_back( finalState.segment( 0, 3 ).norm( ));
//
//            // Calculate the order p from two subsequent errors, using the fact that
//            // error_{i} / error_{i-1} = ( step_{i} / step_{i-1} )^p
//            if ( i > 0 )
//            {
//                calculatedOrder.push_back( std::log( errors.at( i ) / errors.at( i - 1 )) /
//                                           std::log( timeSteps.at( i ) / timeSteps.at( i - 1 )));
//                std::cout << calculatedOrder.at( i - 1 ) << " " << errors.at( i ) << std::endl;
//            }
//        }
//        std::cout << std::endl;
//
//        // Calculate mean and standard deviation of values of order computed from eacg two subsequent steps
//        double meanOrder =
//            std::accumulate( calculatedOrder.begin( ), calculatedOrder.end( ), 0.0 ) / calculatedOrder.size( );
//        double standardDeviationOrder = std::sqrt(
//            std::inner_product( calculatedOrder.begin( ), calculatedOrder.end( ), calculatedOrder.begin( ), 0.0 ) /
//            calculatedOrder.size( ) - meanOrder * meanOrder );
//
////            // Check that the order is in the right range.
////            BOOST_CHECK( meanOrder > ( 2.0 * static_cast< double >( j ) ) );
////            BOOST_CHECK( meanOrder < ( 2.0 * static_cast< double >( j ) ) + 2 );
////
////            // Check that the value of the order is reasonable constant; the high value for the test here is only needed for
////            // the 13th order method.
////            BOOST_CHECK( standardDeviationOrder < 0.5 );
//    }
//}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
