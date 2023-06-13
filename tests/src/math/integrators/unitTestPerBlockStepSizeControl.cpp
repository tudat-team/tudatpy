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
    {
        double radius = 2.0;
        double period = 1.0;
        double angularRate = 2.0 * mathematical_constants::PI / period;
        double initiaStep = 0.01;
        double initialTime = 0.0;
        double tolerances = 1.0E-14;

        Eigen::VectorXd initialState =
            circleState( initialTime, radius, angularRate );

        std::function< Eigen::VectorXd( const double, const Eigen::VectorXd& ) > stateDerivativeFunction =
            std::bind( &circleStateDerivative, std::placeholders::_1, std::placeholders::_2,
                       radius, angularRate );

        std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< MultiStageVariableStepSizeSettings< > >
                ( initiaStep, rungeKuttaFehlberg45,
                  std::make_shared< PerElementIntegratorStepSizeControlSettings< double > >( tolerances, 0.0 ),
                  std::make_shared< IntegratorStepSizeValidationSettings >( std::numeric_limits< double >::min( ),
                                                                            std::numeric_limits< double >::max( ),
                                                                            set_to_minimum_step_silently ) );

        std::shared_ptr< numerical_integrators::NumericalIntegrator< > > integrator = createIntegrator< double, Eigen::VectorXd >(
            stateDerivativeFunction, initialState, initialTime, integratorSettings );

        double finalTime = 10.0;
        std::map< double, Eigen::VectorXd > stateHistory;
        stateHistory[ integrator->getCurrentIndependentVariable( ) ] =
            integrator->getCurrentState( );
        double minimumStep = std::numeric_limits< double >::infinity( );
        double maximumStep = 0.0;

        Eigen::VectorXd stateAtMinimumStep;
        Eigen::VectorXd stateAtMaximumStep;

        while( integrator->getCurrentIndependentVariable( ) < finalTime )
        {
            integrator->performIntegrationStep( integrator->getNextStepSize( ) );
            double timeStep = integrator->getCurrentIndependentVariable( ) - stateHistory.rbegin( )->first;
            if( timeStep > maximumStep )
            {
                maximumStep = timeStep;
                stateAtMaximumStep = integrator->getCurrentState( );
            }

            if( timeStep < minimumStep )
            {
                minimumStep = timeStep;
                stateAtMinimumStep = integrator->getCurrentState( );
            }

            stateHistory[ integrator->getCurrentIndependentVariable( ) ] =
                integrator->getCurrentState( );
        }

        std::cout<<minimumStep<<" "<<stateAtMinimumStep

    }

}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
