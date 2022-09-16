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
 *      Montenbruck, O., Gill, E. Satellite Orbits: Models, Methods, Applications, Springer, 2005.
 *      The MathWorks, Inc. RKF54b, Symbolic Math Toolbox, 2012.
 *
 *    Notes
 *      For the tests using data from the Symbolic Math Toolbox (MathWorks, 2012), the single step
 *      and full integration error tolerances were picked to be as small as possible, without
 *      causing the tests to fail. These values are not deemed to indicate any bugs in the code;
 *      however, it is important to take these discrepancies into account when using this numerical
 *      integrator.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <limits>
#include <cmath>

#include <Eigen/Core>

#include <boost/test/unit_test.hpp>

#include "tudat/math/integrators/rungeKuttaVariableStepSizeIntegrator.h"
#include "tudat/math/integrators/rungeKuttaCoefficients.h"
#include "tudat/math/integrators/tests/numericalIntegratorTestFunctions.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_runge_kutta_fehlberg_56_integrator )

using numerical_integrator_test_functions::computeNonAutonomousModelStateDerivative;
using numerical_integrator_test_functions::computeVanDerPolStateDerivative;
using numerical_integrator_test_functions::computeFehlbergLogirithmicTestODEStateDerivative;
using numerical_integrator_test_functions::computeAnalyticalStateFehlbergODE;

using namespace numerical_integrators;

//! Compare with analytical solution of Fehlberg
BOOST_AUTO_TEST_CASE( test_RungeKuttaFehlberg56_Integrator_Fehlberg_Benchmark )
{
    RungeKuttaCoefficients coeff56 =
        RungeKuttaCoefficients::get( CoefficientSets::rungeKuttaFehlberg56);

    // Integrator settings
    double minimumStepSize   = std::numeric_limits< double >::epsilon( );
    double maximumStepSize   = std::numeric_limits< double >::infinity( );
    double initialStepSize   = 1E-6; // Don't make this too small
    double relativeTolerance = 1E-16;
    double absoluteTolerance = 1E-16;

    // Initial conditions
    double initialTime = 0.0;
    double finalTime   = 5.0;
    Eigen::Vector2d initialState( exp( 1.0 ), 1.0);

    // Setup integrator
    RungeKuttaVariableStepSizeIntegratorXd integrator56(
                coeff56, computeFehlbergLogirithmicTestODEStateDerivative,
                initialTime, initialState, minimumStepSize,
                maximumStepSize, relativeTolerance, absoluteTolerance );


    // Obtain numerical solution
    Eigen::Vector2d numericalSolution = integrator56.integrateTo( finalTime, initialStepSize );

    // Analytical solution
    // (page 30, Fehlberg, E. (1968). Classical Fifth-, Sixth-, Seventh- and Eigth-Order Runge-Kutta
    // Formulas with Stepsize Control)
    Eigen::Vector2d analyticalSolution =
        computeAnalyticalStateFehlbergODE( finalTime, initialState );

    Eigen::Vector2d computedError = numericalSolution - analyticalSolution;
    BOOST_CHECK_SMALL( std::fabs( computedError( 0 ) ), 1E-12 );
    BOOST_CHECK_SMALL( std::fabs( computedError( 1 ) ), 1E-12 );

    // Error calculated by -> Fehlberg, E. (1968) page 30
    // Initial stepsize unknown..
    Eigen::VectorXd fehlbergError( 2 );
    fehlbergError << 0.1072E-12, -0.2190E-12;

    // Sign check
    // Not always same sign -> initial step size = 1 or 1E-2, failure: computedError( 1 ),fehlbergError( 1 ) not same sign
//    BOOST_CHECK_GE( computedError( 0 ) / fehlbergError( 0 ), 0.0 );
//    BOOST_CHECK_GE( computedError( 1 ) / fehlbergError( 1 ), 0.0 );

    // Check error is similar in magnitude
    BOOST_CHECK_SMALL( std::fabs( computedError( 0 ) / fehlbergError( 0 ) ), 3.0);
    BOOST_CHECK_SMALL( std::fabs( computedError( 1 ) / fehlbergError( 1 ) ), 1.0);
}


//! Test Compare with Runge Kutta 78
BOOST_AUTO_TEST_CASE( test_RungeKuttaFehlberg56_Integrator_Compare78 )
{
    // Setup integrator
    RungeKuttaCoefficients coeff56 =
            RungeKuttaCoefficients::get( CoefficientSets::rungeKuttaFehlberg56);

    RungeKuttaCoefficients coeff78 =
            RungeKuttaCoefficients::get( CoefficientSets::rungeKuttaFehlberg78);

    // Integrator settings
    double minimumStepSize = std::numeric_limits< double >::epsilon( );
    double maximumStepSize = std::numeric_limits< double >::infinity( );
    double initialStepSize = 1E-4; // Don't make this too small
    double relativeTolerance = 1E-15;
    double absoluteTolerance = 1E-15;

    // Initial conditions
    double initialTime = 0.5;
    Eigen::VectorXd InitialState( 1 );
    InitialState << 0.5; // 1 large error

    // Setup integrator
    RungeKuttaVariableStepSizeIntegratorXd integrator56(
                coeff56, computeNonAutonomousModelStateDerivative, initialTime, InitialState,
                minimumStepSize, maximumStepSize, relativeTolerance, absoluteTolerance );

    RungeKuttaVariableStepSizeIntegratorXd integrator78(
                coeff78, computeNonAutonomousModelStateDerivative, initialTime, InitialState,
                minimumStepSize, maximumStepSize, relativeTolerance, absoluteTolerance );

    double endTime = 1.5;
    Eigen::VectorXd solution56 = integrator56.integrateTo( endTime, initialStepSize );
    Eigen::VectorXd solution78 = integrator78.integrateTo( endTime, initialStepSize );

    Eigen::VectorXd difference = solution78 - solution56;

    BOOST_CHECK_SMALL( std::fabs( difference( 0 ) ), 1E-13 );
}

//! Test Compare with Runge Kutta 78
BOOST_AUTO_TEST_CASE( test_RungeKuttaFehlberg56_Integrator_Compare78_v2 )
{
    // Setup integrator
    RungeKuttaCoefficients coeff56 =
            RungeKuttaCoefficients::get(
                rungeKuttaFehlberg56 );

    RungeKuttaCoefficients coeff78 =
            RungeKuttaCoefficients::get(
                rungeKuttaFehlberg78 );

    // Integrator settings
    double minimumStepSize = std::numeric_limits< double >::epsilon( );
    double maximumStepSize = std::numeric_limits< double >::infinity( );
    double initialStepSize = 1.0; // Don't make this too small
    double relativeTolerance = 1E-10;
    double absoluteTolerance = 1E-10;

    // Initial conditions
    double initialTime = 0.2;
    Eigen::VectorXd InitialState( 1 );
    InitialState << -1.0;

    // Setup integrator
    RungeKuttaVariableStepSizeIntegratorXd integrator56(
                coeff56, computeNonAutonomousModelStateDerivative, initialTime, InitialState,
                minimumStepSize, maximumStepSize, relativeTolerance, absoluteTolerance );

    RungeKuttaVariableStepSizeIntegratorXd integrator78(
                coeff78, computeNonAutonomousModelStateDerivative, initialTime, InitialState,
                minimumStepSize, maximumStepSize, relativeTolerance, absoluteTolerance );

    double endTime = 2.0;
    Eigen::VectorXd solution56 = integrator56.integrateTo( endTime, initialStepSize );
    Eigen::VectorXd solution78 = integrator78.integrateTo( endTime, initialStepSize );

    Eigen::VectorXd difference = solution78 - solution56;

    BOOST_CHECK_SMALL( std::fabs( difference( 0 ) ), 1E-8 );
}

//! Test Compare with Runge Kutta 78
BOOST_AUTO_TEST_CASE( test_RungeKuttaFehlberg56_Integrator_Compare78_VanDerPol )
{
    // Setup integrator
    RungeKuttaCoefficients coeff56 =
            RungeKuttaCoefficients::get( CoefficientSets::rungeKuttaFehlberg56 );

    RungeKuttaCoefficients coeff78 =
            RungeKuttaCoefficients::get( CoefficientSets::rungeKuttaFehlberg78 );

    // Integrator settings
    double minimumStepSize = std::numeric_limits< double >::epsilon( );
    double maximumStepSize = std::numeric_limits< double >::infinity( );
    double initialStepSize = 1; // Don't make this too small
    double relativeTolerance = 1E-15;
    double absoluteTolerance = 1E-15;

    // Initial conditions
    double initialTime = 0.2;
    Eigen::VectorXd InitialState( 2 );
    InitialState << -1.0, 1.0;

    // Setup integrator
    RungeKuttaVariableStepSizeIntegratorXd integrator56(
                coeff56, computeVanDerPolStateDerivative, initialTime, InitialState, minimumStepSize,
                maximumStepSize, relativeTolerance, absoluteTolerance );

    RungeKuttaVariableStepSizeIntegratorXd integrator78(
                coeff78, computeVanDerPolStateDerivative, initialTime, InitialState, minimumStepSize,
                maximumStepSize, relativeTolerance, absoluteTolerance );

    double endTime = 1.4;
    Eigen::VectorXd solution56 = integrator56.integrateTo(endTime,initialStepSize);
    Eigen::VectorXd solution78 = integrator78.integrateTo(endTime,initialStepSize);

    Eigen::VectorXd difference = solution78 - solution56;

    BOOST_CHECK_SMALL( std::fabs( difference( 0 ) ), 1E-13 );
    BOOST_CHECK_SMALL( std::fabs( difference( 1 ) ), 1E-13 );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
