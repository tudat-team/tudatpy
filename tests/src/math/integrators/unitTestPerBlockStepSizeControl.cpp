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
 *      Burden, R.L., Faires, J.D. Numerical Analysis, 7th Edition, Books/Cole, 2001.
 *      The MathWorks, Inc. Symbolic Math Toolbox, 2012.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <string>


#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "tudat/io/basicInputOutput.h"
#include "tudat/math/integrators/createNumericalIntegrator.h"
#include "tudat/simulation/environment_setup.h"
#include "tudat/simulation/propagation_setup.h"

namespace tudat
{
namespace unit_tests
{

Eigen::VectorXd circleState(
    const double currentTime,
    const double radius, const double angularRate )
{
    Eigen::VectorXd state = Eigen::VectorXd::Zero( 4, 1 );
    state( 0 ) = std::sin( angularRate * currentTime );
    state( 1 ) = std::cos( angularRate * currentTime );
    state( 2 ) = angularRate * std::cos( angularRate * currentTime );
    state( 3 ) = -angularRate * std::sin( angularRate * currentTime );
    state.segment( 0, 4 ) *= radius;
    return state;
}


Eigen::VectorXd circleStateDerivative(
    const double currentTime, const Eigen::VectorXd& currentState,
    const double radius, const double angularRate )
{
    Eigen::VectorXd stateDerivative = Eigen::VectorXd::Zero( 4, 1 );
    stateDerivative.segment( 0, 2 ) = currentState.segment( 2, 2 );
    stateDerivative( 2 ) = -std::sin( angularRate * currentTime );
    stateDerivative( 3 ) = -std::cos( angularRate * currentTime );
    stateDerivative.segment( 2, 2 ) *= radius * angularRate * angularRate;
    return stateDerivative;
}

using namespace numerical_integrators;

BOOST_AUTO_TEST_SUITE( test_block_step_size_control )

BOOST_AUTO_TEST_CASE( testPerBlockCircleStepSizeControl )
{

    double radius = 2.0;
    double period = 1.0;
    double angularRate = 2.0 * mathematical_constants::PI / period;
    double initiaStep = 1.0;
    double initialTime = 0.0;
    double tolerances = 1.0E-14;

    Eigen::VectorXd initialState =
        circleState( initialTime, radius, angularRate );

    std::function<Eigen::VectorXd( const double, const Eigen::VectorXd & )> stateDerivativeFunction =
        std::bind( &circleStateDerivative, std::placeholders::_1, std::placeholders::_2,
                   radius, angularRate );
    for ( unsigned int test = 0; test < 2; test++ )
    {

        std::shared_ptr<IntegratorSettings<> > integratorSettings;
        if( test == 0 )
        {
            integratorSettings = std::make_shared<MultiStageVariableStepSizeSettings<> >
                ( initiaStep, rungeKuttaFehlberg45,
                  std::make_shared<PerElementIntegratorStepSizeControlSettings<double> >( tolerances, 0.0 ),
                  std::make_shared<IntegratorStepSizeValidationSettings>( std::numeric_limits<double>::min( ),
                                                                          std::numeric_limits<double>::max( ),
                                                                          set_to_minimum_step_silently ));
        }
        else
        {
            std::vector< std::pair< int, int > > blocks;
            blocks.push_back( std::make_pair( 0, 2 ) );
            blocks.push_back( std::make_pair( 2, 2 ) );

            integratorSettings = std::make_shared<MultiStageVariableStepSizeSettings<> >
                ( initiaStep, rungeKuttaFehlberg45,
                  std::make_shared<PerBlockIntegratorStepSizeControlSettings<double> >(
                      [=](const int, const int){ return blocks; }, tolerances, 0.0 ),
                  std::make_shared<IntegratorStepSizeValidationSettings>( std::numeric_limits<double>::min( ),
                                                                          std::numeric_limits<double>::max( ),
                                                                          set_to_minimum_step_silently ));
        }

        std::shared_ptr<numerical_integrators::NumericalIntegrator<> >
            integrator = createIntegrator<double, Eigen::VectorXd>(
            stateDerivativeFunction, initialState, initialTime, integratorSettings );

        double finalTime = 10.0;
        std::map<double, Eigen::VectorXd> stateHistory;
        stateHistory[ integrator->getCurrentIndependentVariable( ) ] =
            integrator->getCurrentState( );
        double minimumStep = std::numeric_limits<double>::infinity( );
        double maximumStep = 0.0;

        Eigen::VectorXd stateAtMinimumStep;
        Eigen::VectorXd stateAtMaximumStep;

        double timeOfMinimumStep = TUDAT_NAN;
        double timeOfMaximumStep = TUDAT_NAN;

        while ( integrator->getCurrentIndependentVariable( ) < finalTime )
        {
            integrator->performIntegrationStep( integrator->getNextStepSize( ));
            double timeStep = integrator->getCurrentIndependentVariable( ) - stateHistory.rbegin( )->first;
            if ( timeStep > maximumStep )
            {
                maximumStep = timeStep;
                stateAtMaximumStep = integrator->getCurrentState( );
                timeOfMaximumStep = integrator->getCurrentIndependentVariable( );
            }

            if ( timeStep < minimumStep )
            {
                minimumStep = timeStep;
                stateAtMinimumStep = integrator->getCurrentState( );
                timeOfMinimumStep = integrator->getCurrentIndependentVariable( );
            }

            stateHistory[ integrator->getCurrentIndependentVariable( ) ] =
                integrator->getCurrentState( );
        }

        Eigen::VectorXd finalError = circleState(
            integrator->getCurrentIndependentVariable( ), radius, angularRate ) -
                integrator->getCurrentState( );
        std::cout<<finalError.transpose( )<<std::endl;

        double timeStepRatio = maximumStep / minimumStep;
        if( test == 1 )
        {
            BOOST_CHECK( ( timeStepRatio - 1.0 ) < 0.2 );
        }

        if( test == 0 )
        {
            BOOST_CHECK( timeStepRatio > 8 );
            BOOST_CHECK( ( std::fabs( stateAtMinimumStep( 0 ) ) < 1.0E-2 ) || ( std::fabs( stateAtMinimumStep( 1 ) ) < 1.0E-2 ) );
        }
    }
}

BOOST_AUTO_TEST_CASE( testCowellPropagatorKeplerCompare )
{
    using namespace tudat::numerical_integrators;
    using namespace tudat::spice_interface;
    using namespace tudat::simulation_setup;
    using namespace tudat::basic_astrodynamics;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::propagators;
    using namespace tudat;

    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Define bodies in simulation.
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );

    // Specify initial time
    double initialEphemerisTime = 1.0E7;
    double finalEphemerisTime = initialEphemerisTime + 86400.0;

    // Create bodies needed in simulation
    BodyListSettings bodySettings =
        getDefaultBodySettings(
            bodyNames, "Earth", "ECLIPJ2000" );
    SystemOfBodies bodies = createSystemOfBodies< double, double >( bodySettings );
    bodies.createEmptyBody( "Vehicle" );

    // Set accelerations between bodies that are to be taken into account.
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;

    // Propagate the moon only
    std::vector< std::string > bodiesToIntegrate;
    std::vector< std::string > centralBodies;
    bodiesToIntegrate.push_back( "Vehicle" );
    centralBodies.push_back( "Earth" );

    //Define initial position of satellite at the perigee
    Eigen::Vector6d initialKeplerElements = Eigen::Vector6d::Zero( );
    initialKeplerElements[ semiMajorAxisIndex ] = 8.0E6;
    initialKeplerElements[ inclinationIndex ] = 85.0 * mathematical_constants::PI / 180.0;
    initialKeplerElements[ trueAnomalyIndex ] = 10.0 * mathematical_constants::PI / 180.0;

    Eigen::VectorXd systemInitialState = convertKeplerianToCartesianElements(
        initialKeplerElements, bodies.getBody( "Earth" )->getGravitationalParameter( ) );

    // Define settings for numerical integrator.
    for ( unsigned int test = 0; test < 2; test++ )
    {
        double initialStep = 10.0;
        double tolerance = 1.0E-14;
        std::shared_ptr<IntegratorSettings<> > integratorSettings;
        if ( test == 0 )
        {
            integratorSettings = std::make_shared<MultiStageVariableStepSizeSettings<> >
                ( initialStep, rungeKuttaFehlberg45,
                  std::make_shared<PerElementIntegratorStepSizeControlSettings<double> >( tolerance, tolerance ),
                  std::make_shared<IntegratorStepSizeValidationSettings>( std::numeric_limits<double>::min( ),
                                                                          std::numeric_limits<double>::max( ),
                                                                          set_to_minimum_step_silently ));
        }
        else
        {
            std::vector<std::pair<int, int> > blocks;
            blocks.push_back( std::make_pair( 0, 3 ));
            blocks.push_back( std::make_pair( 3, 3 ));

            integratorSettings = std::make_shared<MultiStageVariableStepSizeSettings<> >
                ( initialStep, rungeKuttaFehlberg45,
                  std::make_shared<PerBlockIntegratorStepSizeControlSettings<double> >(
                      [ = ]( const int, const int )
                      { return blocks; }, tolerance, tolerance ),
                  std::make_shared<IntegratorStepSizeValidationSettings>( std::numeric_limits<double>::min( ),
                                                                          std::numeric_limits<double>::max( ),
                                                                          set_to_minimum_step_silently ));
        }


        // Create acceleration models and propagation settings.
        AccelerationMap accelerationModelMap = createAccelerationModelsMap(
            bodies, accelerationMap, bodiesToIntegrate, centralBodies );
        std::shared_ptr<TranslationalStatePropagatorSettings<double, double> > propagatorSettings =
            std::make_shared<TranslationalStatePropagatorSettings<double, double> >
                ( centralBodies, accelerationModelMap, bodiesToIntegrate, systemInitialState,
                  initialEphemerisTime, integratorSettings,
                  std::make_shared<PropagationTimeTerminationSettings>( finalEphemerisTime ));

        // Create dynamics simulation object.
        SingleArcDynamicsSimulator<double, double> dynamicsSimulator(
            bodies, propagatorSettings );

        double minimumStep = std::numeric_limits<double>::infinity( );
        double maximumStep = 0.0;

        Eigen::VectorXd stateAtMinimumStep;
        Eigen::VectorXd stateAtMaximumStep;

        double timeOfMinimumStep = TUDAT_NAN;
        double timeOfMaximumStep = TUDAT_NAN;

        std::map<double, Eigen::VectorXd> stateHistory = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
        auto firstIterator = stateHistory.begin( );
        auto secondIterator = stateHistory.begin( );
        secondIterator++;

        while ( secondIterator != stateHistory.end( ))
        {
            double timeStep = secondIterator->first - firstIterator->first;
            if ( timeStep > maximumStep )
            {
                maximumStep = timeStep;
                stateAtMaximumStep = firstIterator->second;
                timeOfMaximumStep = firstIterator->first;
            }

            if ( timeStep < minimumStep )
            {
                minimumStep = timeStep;
                stateAtMinimumStep = firstIterator->second;
                timeOfMinimumStep = firstIterator->first;;
            }
            firstIterator++;
            secondIterator++;
        }

        Eigen::Vector6d stateError = ( convertKeplerianToCartesianElements( propagateKeplerOrbit(
            initialKeplerElements, stateHistory.rbegin( )->first - initialEphemerisTime,
            bodies.getBody( "Earth" )->getGravitationalParameter( ) ) , bodies.getBody( "Earth" )->getGravitationalParameter( ) )
                                                           - stateHistory.rbegin( )->second );
//        std::cout<<stateError.segment( 0, 3 ).norm( )<<" "<<stateError.segment( 3, 3 ).norm( )<<std::endl;
//
//        std::cout<< dynamicsSimulator.getCumulativeNumberOfFunctionEvaluations( ).rbegin( )->second<<std::endl<<std::endl;

        double timeStepRatio = maximumStep / minimumStep;
        if ( test == 1 )
        {
            BOOST_CHECK(( timeStepRatio - 1.0 ) < 0.2 );
        }

        if ( test == 0 )
        {
            double minimumPositionRatio = stateAtMinimumStep.segment( 0, 3 ).cwiseAbs().minCoeff( ) / stateAtMinimumStep.segment( 0, 3 ).norm( );
            double minimumVelocityRatio = stateAtMinimumStep.segment( 3, 3 ).cwiseAbs().minCoeff( ) / stateAtMinimumStep.segment( 3, 3 ).norm( );

            BOOST_CHECK( timeStepRatio > 8 );
            BOOST_CHECK( std::min( minimumPositionRatio, minimumVelocityRatio ) < 2.0E-4 );
        }
    }
}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
