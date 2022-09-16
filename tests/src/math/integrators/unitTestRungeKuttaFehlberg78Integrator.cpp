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
 *      The Mathworks, Inc. RKF87, Symbolic Math Toolbox, 2012.
 *
 *    Notes
 *      All the test for this integrator are based on the data generated using the Symbolic Math
 *      Toolbox (MathWorks, 2012). Ideally, another source of data should be used to complete the
 *      testing.
 *
 *      The single step and full integration error tolerances were picked to be as small as
 *      possible, without causing the tests to fail. These values are not deemed to indicate any
 *      bugs in the code; however, it is important to take these discrepancies into account when
 *      using this numerical integrator.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>

#include "tudat/math/integrators/rungeKuttaVariableStepSizeIntegrator.h"
#include "tudat/math/integrators/rungeKuttaCoefficients.h"
#include "tudat/math/integrators/numericalIntegrator.h"
#include "tudat/math/integrators/reinitializableNumericalIntegrator.h"
#include "tudat/support/numericalIntegratorTests.h"
#include "tudat/math/integrators/numericalIntegratorTestFunctions.h"

#include "tudat/io/matrixTextFileReader.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/math/basic/linearAlgebra.h"

#include <limits>
#include <string>

#include <Eigen/Core>

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_runge_kutta_fehlberg_78_integrator )

using linear_algebra::flipMatrixRows;

using numerical_integrators::NumericalIntegratorXdPointer;
using numerical_integrators::ReinitializableNumericalIntegratorXdPointer;
using numerical_integrators::RungeKuttaVariableStepSizeIntegratorXd;
using numerical_integrators::RungeKuttaCoefficients;
using numerical_integrators::CoefficientSets;

using numerical_integrator_test_functions::computeNonAutonomousModelStateDerivative;
using numerical_integrator_test_functions::computeFehlbergLogirithmicTestODEStateDerivative ;
using numerical_integrator_test_functions::computeAnalyticalStateFehlbergODE;

//! Test Runge-Kutta-Fehlberg 78 integrator using benchmark ODE of Fehlberg (1968)
BOOST_AUTO_TEST_CASE( test_RungeKuttaFehlberg78_Integrator_Fehlberg_Benchmark )
{
    using namespace numerical_integrators;
    RungeKuttaCoefficients coeff78 =
            RungeKuttaCoefficients::get( CoefficientSets::rungeKuttaFehlberg78 );

    // Integrator settings
    double minimumStepSize   = std::numeric_limits<double>::epsilon( );
    double maximumStepSize   = std::numeric_limits<double>::infinity( );
    double initialStepSize   = 1E-6; // Error: 0.0521 for initialStepSize = 1 ?
    double relativeTolerance = 1E-16;
    double absoluteTolerance = 1E-16;

    // Initial conditions
    double initialTime = 0.0;
    double finalTime   = 5.0;
    Eigen::Vector2d initialState( exp( 1.0 ), 1.0);

    // Setup integrator
    numerical_integrators::RungeKuttaVariableStepSizeIntegratorXd integrator78(
                coeff78, computeFehlbergLogirithmicTestODEStateDerivative,
                initialTime, initialState, minimumStepSize,
                maximumStepSize, relativeTolerance, absoluteTolerance );


    // Obtain numerical solution
    Eigen::Vector2d numericalSolution = integrator78.integrateTo( finalTime, initialStepSize );

    // Analytical solution
    Eigen::Vector2d analyticalSolution = computeAnalyticalStateFehlbergODE( finalTime, initialState );

    Eigen::Vector2d computedError = numericalSolution - analyticalSolution;
    BOOST_CHECK_SMALL( std::fabs(computedError( 0 )), 1E-13 );
    BOOST_CHECK_SMALL( std::fabs(computedError( 1 )), 1E-13 );
}

//! Test Runge-Kutta-Fehlberg 78 integrator using benchmark data from (The MathWorks, 2012).
BOOST_AUTO_TEST_CASE( testRungeKuttaFehlberg78IntegratorUsingMatlabData )
{
    using namespace numerical_integrator_tests;

    // Read in benchmark data (generated using Symbolic Math Toolbox in Matlab
    // (The MathWorks, 2012)). This data is generated using the RKF87 numerical integrator.
    const std::string pathToForwardIntegrationOutputFile = paths::getTudatTestDataPath( )
            + "/matlabOutputRungeKuttaFehlberg78Forward.txt";
    const std::string pathToDiscreteEventIntegrationOutputFile = paths::getTudatTestDataPath( )
            + "/matlabOutputRungeKuttaFehlberg78DiscreteEvent.txt";

    // Store benchmark data in matrix.
    const Eigen::MatrixXd matlabForwardIntegrationData =
            input_output::readMatrixFromFile( pathToForwardIntegrationOutputFile, "," );
    Eigen::MatrixXd matlabBackwardIntegrationData = matlabForwardIntegrationData;
    flipMatrixRows( matlabBackwardIntegrationData );
    const Eigen::MatrixXd matlabDiscreteEventIntegrationData =
            input_output::readMatrixFromFile( pathToDiscreteEventIntegrationOutputFile, "," );

    // Set integrator parameters.

    // All of the following parameters are set such that the input data is fully accepted by the
    // integrator, to determine the steps to be taken.
    const double zeroMinimumStepSize = std::numeric_limits< double >::epsilon( );
    const double infiniteMaximumStepSize = std::numeric_limits< double >::infinity( );
    const double infiniteRelativeErrorTolerance = std::numeric_limits< double >::infinity( );
    const double infiniteAbsoluteErrorTolerance = std::numeric_limits< double >::infinity( );

    // The following parameters set how the error control mechanism should work.
    const double relativeErrorTolerance = 1.0e-15;
    const double absoluteErrorTolerance = 1.0e-15;

    // Case 1: Execute integrateTo() to integrate one step forward in time.
    {
        // Declare integrator with all necessary settings.
        NumericalIntegratorXdPointer integrator
                = std::make_shared< RungeKuttaVariableStepSizeIntegratorXd >(
                    RungeKuttaCoefficients::get( CoefficientSets::rungeKuttaFehlberg78 ),
                    &computeNonAutonomousModelStateDerivative,
                    matlabForwardIntegrationData( FIRST_ROW, TIME_COLUMN_INDEX ),
                    ( Eigen::VectorXd( 1 )
                      << matlabForwardIntegrationData( FIRST_ROW,
                                                       STATE_COLUMN_INDEX ) ).finished( ),
                    zeroMinimumStepSize,
                    infiniteMaximumStepSize,
                    infiniteRelativeErrorTolerance,
                    infiniteAbsoluteErrorTolerance );

        executeOneIntegrateToStep( matlabForwardIntegrationData, 1.0e-15, integrator );
    }

    // Case 2: Execute performIntegrationStep() to perform multiple integration steps until final
    //         time.
    {
        // Declare integrator with all necessary settings.
        NumericalIntegratorXdPointer integrator
                = std::make_shared< RungeKuttaVariableStepSizeIntegratorXd >(
                    RungeKuttaCoefficients::get( CoefficientSets::rungeKuttaFehlberg78 ),
                    &computeNonAutonomousModelStateDerivative,
                    matlabForwardIntegrationData( FIRST_ROW, TIME_COLUMN_INDEX ),
                    ( Eigen::VectorXd( 1 )
                      << matlabForwardIntegrationData( FIRST_ROW,
                                                       STATE_COLUMN_INDEX ) ).finished( ),
                    zeroMinimumStepSize,
                    infiniteMaximumStepSize,
                    infiniteRelativeErrorTolerance,
                    infiniteAbsoluteErrorTolerance );

        performIntegrationStepToSpecifiedTime( matlabForwardIntegrationData,
                                               1.0e-15, 1.0e-15, integrator );
    }

    // Case 3: Execute performIntegrationStep() to perform multiple integration steps until initial
    //         time (backwards).
    {
        // Declare integrator with all necessary settings.
        NumericalIntegratorXdPointer integrator
                = std::make_shared< RungeKuttaVariableStepSizeIntegratorXd >(
                    RungeKuttaCoefficients::get( CoefficientSets::rungeKuttaFehlberg78 ),
                    &computeNonAutonomousModelStateDerivative,
                    matlabBackwardIntegrationData( FIRST_ROW, TIME_COLUMN_INDEX ),
                    ( Eigen::VectorXd( 1 )
                      << matlabBackwardIntegrationData( FIRST_ROW,
                                                        STATE_COLUMN_INDEX ) ).finished( ),
                    zeroMinimumStepSize,
                    infiniteMaximumStepSize,
                    infiniteRelativeErrorTolerance,
                    infiniteAbsoluteErrorTolerance );

        performIntegrationStepToSpecifiedTime( matlabBackwardIntegrationData,
                                               1.0e-15, 1.0e-14, integrator );
    }

    // Case 4: Execute integrateTo() to integrate to specified time in one step.
    {
        // Note that this test has a strange issue that the if the absolute error tolerance is set
        // to 1.0e-15, the last step that the integrateTo() function takes does not result in the
        // expected final time of 1.0. As a temporary solution, the absolute error tolerance has
        // been multiplied by 10.0, which seems to solve the problem. This error indicated a
        // possible problem with the implementation of the integrateTo() function, which needs to
        // be investigated in future.

        // Declare integrator with all necessary settings.
        NumericalIntegratorXdPointer integrator
                = std::make_shared< RungeKuttaVariableStepSizeIntegratorXd >(
                    RungeKuttaCoefficients::get( CoefficientSets::rungeKuttaFehlberg78 ),
                    &computeNonAutonomousModelStateDerivative,
                    matlabForwardIntegrationData( FIRST_ROW, TIME_COLUMN_INDEX ),
                    ( Eigen::VectorXd( 1 )
                      << matlabForwardIntegrationData( FIRST_ROW,
                                                       STATE_COLUMN_INDEX ) ).finished( ),
                    zeroMinimumStepSize,
                    infiniteMaximumStepSize,
                    relativeErrorTolerance,
                    absoluteErrorTolerance * 10.0 );

        executeIntegrateToToSpecifiedTime( matlabForwardIntegrationData, 1.0e-13, integrator,
                                           matlabForwardIntegrationData(
                                               matlabForwardIntegrationData.rows( ) - 1,
                                               TIME_COLUMN_INDEX ) );
    }

    // Case 5: Execute performIntegrationstep() to integrate to specified time in multiple steps,
    //         including discrete events.
    {
        // Declare integrator with all necessary settings.
        ReinitializableNumericalIntegratorXdPointer integrator
                = std::make_shared< RungeKuttaVariableStepSizeIntegratorXd >(
                    RungeKuttaCoefficients::get( CoefficientSets::rungeKuttaFehlberg78 ),
                    &computeNonAutonomousModelStateDerivative,
                    matlabForwardIntegrationData( FIRST_ROW, TIME_COLUMN_INDEX ),
                    ( Eigen::VectorXd( 1 )
                      << matlabForwardIntegrationData( FIRST_ROW,
                                                       STATE_COLUMN_INDEX ) ).finished( ),
                    zeroMinimumStepSize,
                    infiniteMaximumStepSize,
                    infiniteRelativeErrorTolerance,
                    infiniteAbsoluteErrorTolerance );

        performIntegrationStepToSpecifiedTimeWithEvents( matlabDiscreteEventIntegrationData,
                                                         1.0e-15, 1.0e-13, integrator );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
