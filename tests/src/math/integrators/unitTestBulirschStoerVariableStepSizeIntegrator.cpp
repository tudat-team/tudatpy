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

#include "tudat/math/integrators/bulirschStoerVariableStepsizeIntegrator.h"
#include "tudat/math/integrators/rungeKuttaVariableStepSizeIntegrator.h"
#include "tudat/math/integrators/rungeKuttaCoefficients.h"
#include "tudat/math/integrators/numericalIntegratorTestFunctions.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_bulirsch_stoer_integrator )

using numerical_integrator_test_functions::computeNonAutonomousModelStateDerivative;
using numerical_integrator_test_functions::computeVanDerPolStateDerivative;
using numerical_integrator_test_functions::computeFehlbergLogirithmicTestODEStateDerivative;
using numerical_integrator_test_functions::computeAnalyticalStateFehlbergODE;

using namespace numerical_integrators;

//! Test Compare with Runge Kutta 78
BOOST_AUTO_TEST_CASE( test_BulirschStoer_Integrator_Compare78 )
{
    // Setup integrator
    RungeKuttaCoefficients coeff_rk78 =
            RungeKuttaCoefficients::get( CoefficientSets::rungeKuttaFehlberg78);

    // Integrator settings
    double minimumStepSize = std::numeric_limits< double >::epsilon( );
    double maximumStepSize = std::numeric_limits< double >::infinity( );
    double initialStepSize = 1E-4; // Don't make this too small
    Eigen::VectorXd relativeTolerance(1);
    relativeTolerance << 1E-12;
    Eigen::VectorXd absoluteTolerance = relativeTolerance;

    // Initial conditions
    double initialTime = 0.5;
    Eigen::VectorXd initialState( 1 );
    initialState << 0.5; // 1 large error


    // Setup integrator
    BulirschStoerVariableStepSizeIntegratorXd integrator_bs(
                getBulirschStoerStepSequence( bulirsch_stoer_sequence, 4 ),
                computeNonAutonomousModelStateDerivative, initialTime, initialState,
                minimumStepSize, maximumStepSize, relativeTolerance, absoluteTolerance );

    RungeKuttaVariableStepSizeIntegratorXd integrator_rk78(
                coeff_rk78, computeNonAutonomousModelStateDerivative, initialTime, initialState,
                minimumStepSize, maximumStepSize, relativeTolerance, absoluteTolerance );

    double endTime = 1.5;
    Eigen::VectorXd solution_bs = integrator_bs.integrateTo( endTime, initialStepSize );
    Eigen::VectorXd solution_rk78 = integrator_rk78.integrateTo( endTime, initialStepSize );

    Eigen::VectorXd difference = solution_bs - solution_rk78;

    BOOST_CHECK_SMALL( std::fabs( difference( 0 ) ), 1E-11 );
}

//! Test Compare with Runge Kutta 78
BOOST_AUTO_TEST_CASE( test_BulirschStoer_Integrator_Compare78_v2 )
{
    // Setup integrator
    RungeKuttaCoefficients coeff_rk78 =
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
    Eigen::VectorXd initialState( 1 );
    initialState << -1.0;

    // Setup integrator
    BulirschStoerVariableStepSizeIntegratorXd integrator_bs(
                getBulirschStoerStepSequence( bulirsch_stoer_sequence, 4 ),
                computeNonAutonomousModelStateDerivative, initialTime, initialState,
                minimumStepSize, maximumStepSize, relativeTolerance, absoluteTolerance );

    RungeKuttaVariableStepSizeIntegratorXd integrator_rk78(
                coeff_rk78, computeNonAutonomousModelStateDerivative, initialTime, initialState,
                minimumStepSize, maximumStepSize, relativeTolerance, absoluteTolerance );

    double endTime = 2.0;
    Eigen::VectorXd solution_bs = integrator_bs.integrateTo( endTime, initialStepSize );
    Eigen::VectorXd solution_rk78 = integrator_rk78.integrateTo( endTime, initialStepSize );

    Eigen::VectorXd difference = solution_bs - solution_rk78;

    BOOST_CHECK_SMALL( std::fabs( difference( 0 ) ), 1E-8 );
}

//! Test Compare with Runge Kutta 78
BOOST_AUTO_TEST_CASE( test_BulirschStoer_Integrator_Compare78_VanDerPol )
{
    // Setup integrator
    RungeKuttaCoefficients coeff_rk78 =
            RungeKuttaCoefficients::get( CoefficientSets::rungeKuttaFehlberg78 );

    // Integrator settings
    double minimumStepSize = std::numeric_limits< double >::epsilon( );
    double maximumStepSize = std::numeric_limits< double >::infinity( );
    double initialStepSize = 1; // Don't make this too small
    double relativeTolerance = 1E-15;
    double absoluteTolerance = 1E-15;

    // Initial conditions
    double initialTime = 0.2;
    Eigen::VectorXd initialState( 2 );
    initialState << -1.0, 1.0;

    // Setup integrator
    BulirschStoerVariableStepSizeIntegratorXd integrator_bs(
                getBulirschStoerStepSequence( bulirsch_stoer_sequence, 4 ),
                computeVanDerPolStateDerivative, initialTime, initialState,
                minimumStepSize, maximumStepSize, relativeTolerance, absoluteTolerance );

    RungeKuttaVariableStepSizeIntegratorXd integrator_rk78(
                coeff_rk78, computeVanDerPolStateDerivative, initialTime, initialState,
                minimumStepSize, maximumStepSize, relativeTolerance, absoluteTolerance );

    double endTime = 1.4;
    Eigen::VectorXd solution_bs = integrator_bs.integrateTo( endTime,initialStepSize );
    Eigen::VectorXd solution_rk78 = integrator_rk78.integrateTo( endTime,initialStepSize );

    Eigen::VectorXd difference = solution_rk78 - solution_bs;

    BOOST_CHECK_SMALL( std::fabs( difference( 0 ) ), 5E-12 );
    BOOST_CHECK_SMALL( std::fabs( difference( 1 ) ), 5E-12 );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
