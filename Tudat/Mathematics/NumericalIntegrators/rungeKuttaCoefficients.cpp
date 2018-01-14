/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      The Mathworks, Inc. RKF78, Symbolic Math Toolbox, 2012.
 *      Fehlberg, E. Classical Fifth-, Sixth-, Seventh-, and Eighth-Order Runge-Kutta Formulas With
 *          Stepsize Control, Marshall Spaceflight Center, NASA TR R-278, 1968.
 *      Montenbruck, O., Gill, E. Satellite Orbits: Models, Methods, Applications, Springer, 2005.
 *
 *    Notes
 *      The naming of the coefficient sets follows (Montenbruck and Gill, 2005).
 *
 */

#include <Eigen/Core>

#include "Tudat/Mathematics/NumericalIntegrators/rungeKuttaCoefficients.h"

namespace tudat
{
namespace numerical_integrators
{

//! Initialize RKF45 coefficients.
void initializeRungeKuttaFehlberg45Coefficients( RungeKuttaCoefficients&
                                                 rungeKuttaFehlberg45Coefficients )
{
    // Define characteristics of coefficient set.
    rungeKuttaFehlberg45Coefficients.lowerOrder = 4;
    rungeKuttaFehlberg45Coefficients.higherOrder = 5;
    rungeKuttaFehlberg45Coefficients.orderEstimateToIntegrate = RungeKuttaCoefficients::lower;

    // This coefficient set is taken from (Fehlberg, 1968).

    // Define a-coefficients for the Runge-Kutta-Fehlberg method of order 5
    // with an embedded 4th-order method for stepsize control and a total of 6 stages.
    rungeKuttaFehlberg45Coefficients.aCoefficients = Eigen::MatrixXd::Zero( 6, 5 );
    rungeKuttaFehlberg45Coefficients.aCoefficients( 1, 0 ) = 1.0 / 4.0;

    rungeKuttaFehlberg45Coefficients.aCoefficients( 2, 0 ) = 3.0 / 32.0;
    rungeKuttaFehlberg45Coefficients.aCoefficients( 2, 1 ) = 9.0 / 32.0;

    rungeKuttaFehlberg45Coefficients.aCoefficients( 3, 0 ) =  1932.0 / 2197.0;
    rungeKuttaFehlberg45Coefficients.aCoefficients( 3, 1 ) = -7200.0 / 2197.0;
    rungeKuttaFehlberg45Coefficients.aCoefficients( 3, 2 ) =  7296.0 / 2197.0;

    rungeKuttaFehlberg45Coefficients.aCoefficients( 4, 0 ) = 439.0 / 216.0;
    rungeKuttaFehlberg45Coefficients.aCoefficients( 4, 1 ) = -8.0;
    rungeKuttaFehlberg45Coefficients.aCoefficients( 4, 2 ) = 3680.0 / 513.0;
    rungeKuttaFehlberg45Coefficients.aCoefficients( 4, 3 ) = -845.0 / 4104.0;

    rungeKuttaFehlberg45Coefficients.aCoefficients( 5, 0 ) = -8.0 / 27.0;
    rungeKuttaFehlberg45Coefficients.aCoefficients( 5, 1 ) = 2.0;
    rungeKuttaFehlberg45Coefficients.aCoefficients( 5, 2 ) = -3544.0 / 2565.0;
    rungeKuttaFehlberg45Coefficients.aCoefficients( 5, 3 ) = 1859.0 / 4104.0;
    rungeKuttaFehlberg45Coefficients.aCoefficients( 5, 4 ) = -11.0 / 40.0;


    // Define c-coefficients for the Runge-Kutta-Fehlberg method of order 5
    // with an embedded 4th-order method for stepsize control and a total of 6 stages.
    rungeKuttaFehlberg45Coefficients.cCoefficients = Eigen::VectorXd::Zero( 6 );
    rungeKuttaFehlberg45Coefficients.cCoefficients( 1 ) = 1.0 / 4.0;
    rungeKuttaFehlberg45Coefficients.cCoefficients( 2 ) = 3.0 / 8.0;
    rungeKuttaFehlberg45Coefficients.cCoefficients( 3 ) = 12.0 / 13.0;
    rungeKuttaFehlberg45Coefficients.cCoefficients( 4 ) = 1.0;
    rungeKuttaFehlberg45Coefficients.cCoefficients( 5 ) = 1.0 / 2.0;


    // Define b-coefficients for the Runge-Kutta method of order 5
    // with an embedded 4th-order method for stepsize control and a total of 6 stages.
    rungeKuttaFehlberg45Coefficients.bCoefficients = Eigen::MatrixXd::Zero( 2, 6 );
    rungeKuttaFehlberg45Coefficients.bCoefficients( 0, 0 ) = 25.0 / 216.0;
    rungeKuttaFehlberg45Coefficients.bCoefficients( 0, 2 ) = 1408.0 / 2565.0;
    rungeKuttaFehlberg45Coefficients.bCoefficients( 0, 3 ) = 2197.0 / 4104.0;
    rungeKuttaFehlberg45Coefficients.bCoefficients( 0, 4 ) = -1.0 / 5.0;

    rungeKuttaFehlberg45Coefficients.bCoefficients( 1, 0 ) = 16.0 / 135.0;
    rungeKuttaFehlberg45Coefficients.bCoefficients( 1, 2 ) = 6656.0 / 12825.0;
    rungeKuttaFehlberg45Coefficients.bCoefficients( 1, 3 ) = 28561.0 / 56430.0;
    rungeKuttaFehlberg45Coefficients.bCoefficients( 1, 4 ) = -9.0 / 50.0;
    rungeKuttaFehlberg45Coefficients.bCoefficients( 1, 5 ) = 2.0 / 55.0;
}

//! Initialize RKF56 coefficients.
void initializeRungeKuttaFehlberg56Coefficients( RungeKuttaCoefficients&
                                                 rungeKuttaFehlberg56Coefficients )
{
    // Define characteristics of coefficient set.
    rungeKuttaFehlberg56Coefficients.lowerOrder = 5;
    rungeKuttaFehlberg56Coefficients.higherOrder = 6;
    rungeKuttaFehlberg56Coefficients.orderEstimateToIntegrate = RungeKuttaCoefficients::lower;

    // This coefficient set is taken from (Fehlberg, 1968).

    // Define a-coefficients for the Runge-Kutta-Fehlberg method of order 6
    // with an embedded 5th-order method for stepsize control and a total of 8 stages.
    rungeKuttaFehlberg56Coefficients.aCoefficients = Eigen::MatrixXd::Zero( 8, 7 );
    rungeKuttaFehlberg56Coefficients.aCoefficients( 1, 0 ) = 1.0 / 6.0;

    rungeKuttaFehlberg56Coefficients.aCoefficients( 2, 0 ) =  4.0 / 75.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 2, 1 ) = 16.0 / 75.0;

    rungeKuttaFehlberg56Coefficients.aCoefficients( 3, 0 ) =  5.0 / 6.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 3, 1 ) = -8.0 / 3.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 3, 2 ) =  5.0 / 2.0;

    rungeKuttaFehlberg56Coefficients.aCoefficients( 4, 0 ) =  -8.0 / 5.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 4, 1 ) = 144.0 / 25.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 4, 2 ) =  -4.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 4, 3 ) =  16.0 / 25.0;

    rungeKuttaFehlberg56Coefficients.aCoefficients( 5, 0 ) = 361.0 / 320.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 5, 1 ) = -18.0 / 5.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 5, 2 ) = 407.0 / 128.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 5, 3 ) = -11.0 / 80.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 5, 4 ) =  55.0 / 128.0;

    rungeKuttaFehlberg56Coefficients.aCoefficients( 6, 0 ) = -11.0 / 640.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 6, 2 ) =  11.0 / 256.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 6, 3 ) = -11.0 / 160.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 6, 4 ) =  11.0 / 256.0;

    rungeKuttaFehlberg56Coefficients.aCoefficients( 7, 0 ) = 93.0 / 640.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 7, 1 ) = -18.0 / 5.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 7, 2 ) = 803.0 / 256.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 7, 3 ) = -11.0 / 160.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 7, 4 ) = 99.0 / 256.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 7, 6 ) = 1.0;

    // Define c-coefficients for the Runge-Kutta-Fehlberg method of order 6
    // with an embedded 5th-order method for stepsize control and a total of 8 stages.
    rungeKuttaFehlberg56Coefficients.cCoefficients = Eigen::VectorXd::Zero( 8 );
    rungeKuttaFehlberg56Coefficients.cCoefficients( 1 ) = 1.0 / 6.0;
    rungeKuttaFehlberg56Coefficients.cCoefficients( 2 ) = 4.0 / 15.0;
    rungeKuttaFehlberg56Coefficients.cCoefficients( 3 ) = 2.0 / 3.0;
    rungeKuttaFehlberg56Coefficients.cCoefficients( 4 ) = 4.0 / 5.0;
    rungeKuttaFehlberg56Coefficients.cCoefficients( 5 ) = 1.0;
    rungeKuttaFehlberg56Coefficients.cCoefficients( 7 ) = 1.0;

    // Define b-coefficients for the Runge-Kutta method of order 6
    // with an embedded 5th-order method for stepsize control and a total of 8 stages.
    rungeKuttaFehlberg56Coefficients.bCoefficients = Eigen::MatrixXd::Zero( 2, 8 );
    rungeKuttaFehlberg56Coefficients.bCoefficients( 0, 0 ) = 31.0 / 384.0;
    rungeKuttaFehlberg56Coefficients.bCoefficients( 0, 2 ) = 1125.0 / 2816.0;
    rungeKuttaFehlberg56Coefficients.bCoefficients( 0, 3 ) = 9.0 / 32.0;
    rungeKuttaFehlberg56Coefficients.bCoefficients( 0, 4 ) = 125.0 / 768.0;
    rungeKuttaFehlberg56Coefficients.bCoefficients( 0, 5 ) = 5.0 / 66.0;

    rungeKuttaFehlberg56Coefficients.bCoefficients( 1, 0 ) = 7.0 / 1408.0;
    rungeKuttaFehlberg56Coefficients.bCoefficients( 1, 2 ) = 1125.0 / 2816.0;
    rungeKuttaFehlberg56Coefficients.bCoefficients( 1, 3 ) = 9.0 / 32.0;
    rungeKuttaFehlberg56Coefficients.bCoefficients( 1, 4 ) = 125.0 / 768.0;
    rungeKuttaFehlberg56Coefficients.bCoefficients( 1, 6 ) = 5.0 / 66.0;
    rungeKuttaFehlberg56Coefficients.bCoefficients( 1, 7 ) = 5.0 / 66.0;
}

//! Initialize RKF78 coefficients.
void initializeRungeKuttaFehlberg78Coefficients( RungeKuttaCoefficients&
                                                 rungeKuttaFehlberg78Coefficients )
{
    // Define characteristics of coefficient set.
    rungeKuttaFehlberg78Coefficients.lowerOrder = 7;
    rungeKuttaFehlberg78Coefficients.higherOrder = 8;
    rungeKuttaFehlberg78Coefficients.orderEstimateToIntegrate = RungeKuttaCoefficients::lower;

    // This coefficient set is taken from (Fehlberg, 1968).

    // Define a-coefficients for the Runge-Kutta-Fehlberg method of order 7
    // with an embedded 8th-order method for stepsize control and a total of 13 stages.
    rungeKuttaFehlberg78Coefficients.aCoefficients = Eigen::MatrixXd::Zero( 13, 12 );
    rungeKuttaFehlberg78Coefficients.aCoefficients( 1, 0 ) = 2.0 / 27.0;

    rungeKuttaFehlberg78Coefficients.aCoefficients( 2, 0 ) = 1.0 / 36.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 2, 1 ) = 1.0 / 12.0;

    rungeKuttaFehlberg78Coefficients.aCoefficients( 3, 0 ) = 1.0 / 24.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 3, 2 ) = 1.0 / 8.0;

    rungeKuttaFehlberg78Coefficients.aCoefficients( 4, 0 ) = 5.0 / 12.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 4, 2 ) = -25.0 / 16.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 4, 3 ) = 25.0 / 16.0;

    rungeKuttaFehlberg78Coefficients.aCoefficients( 5, 0 ) = 1.0 / 20.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 5, 3 ) = 1.0 / 4.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 5, 4 ) = 1.0 / 5.0;

    rungeKuttaFehlberg78Coefficients.aCoefficients( 6, 0 ) = -25.0 / 108.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 6, 3 ) = 125.0 / 108.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 6, 4 ) = -65.0 / 27.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 6, 5 ) = 125.0 / 54.0;

    rungeKuttaFehlberg78Coefficients.aCoefficients( 7, 0 ) = 31.0 / 300.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 7, 4 ) = 61.0 / 225.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 7, 5 ) = -2.0 / 9.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 7, 6 ) = 13.0 / 900.0;

    rungeKuttaFehlberg78Coefficients.aCoefficients( 8, 0 ) = 2.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 8, 3 ) = -53.0 / 6.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 8, 4 ) = 704.0 / 45.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 8, 5 ) = -107.0 / 9.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 8, 6 ) = 67.0 / 90.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 8, 7 ) = 3.0;

    rungeKuttaFehlberg78Coefficients.aCoefficients( 9, 0 ) = -91.0 / 108.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 9, 3 ) = 23.0 / 108.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 9, 4 ) = -976.0 / 135.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 9, 5 ) = 311.0 / 54.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 9, 6 ) = -19.0 / 60.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 9, 7 ) = 17.0 / 6.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 9, 8 ) = -1.0 / 12.0;

    rungeKuttaFehlberg78Coefficients.aCoefficients( 10, 0 ) = 2383.0 / 4100.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 10, 3 ) = -341.0 / 164.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 10, 4 ) = 4496.0 / 1025.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 10, 5 ) = -301.0 / 82.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 10, 6 ) = 2133.0 / 4100.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 10, 7 ) = 45.0 / 82.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 10, 8 ) = 45.0 / 164.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 10, 9 ) = 18.0 / 41.0;

    rungeKuttaFehlberg78Coefficients.aCoefficients( 11, 0 ) = 3.0 / 205.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 11, 5 ) = -6.0 / 41.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 11, 6 ) = -3.0 / 205.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 11, 7 ) = -3.0 / 41.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 11, 8 ) = 3.0 / 41.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 11, 9 ) = 6.0 / 41.0;

    rungeKuttaFehlberg78Coefficients.aCoefficients( 12, 0 ) = -1777.0 / 4100.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 12, 3 ) = -341.0 / 164.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 12, 4 ) = 4496.0 / 1025.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 12, 5 ) = -289.0 / 82.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 12, 6 ) = 2193.0 / 4100.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 12, 7 ) = 51.0 / 82.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 12, 8 ) = 33.0 / 164.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 12, 9 ) = 12.0 / 41.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 12, 11 ) = 1.0;

    // Define c-coefficients for the Runge-Kutta-Fehlberg method of order 8
    // with an embedded 7th-order method for stepsize control and a total of 13 stages.
    rungeKuttaFehlberg78Coefficients.cCoefficients = Eigen::VectorXd::Zero( 13 );
    rungeKuttaFehlberg78Coefficients.cCoefficients( 0 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.cCoefficients( 1 ) = 2.0 / 27.0;
    rungeKuttaFehlberg78Coefficients.cCoefficients( 2 ) = 1.0 / 9.0;
    rungeKuttaFehlberg78Coefficients.cCoefficients( 3 ) = 1.0 / 6.0;
    rungeKuttaFehlberg78Coefficients.cCoefficients( 4 ) = 5.0 / 12.0;
    rungeKuttaFehlberg78Coefficients.cCoefficients( 5 ) = 1.0 / 2.0;
    rungeKuttaFehlberg78Coefficients.cCoefficients( 6 ) = 5.0 / 6.0;
    rungeKuttaFehlberg78Coefficients.cCoefficients( 7 ) = 1.0 / 6.0;
    rungeKuttaFehlberg78Coefficients.cCoefficients( 8 ) = 2.0 / 3.0;
    rungeKuttaFehlberg78Coefficients.cCoefficients( 9 ) = 1.0 / 3.0;
    rungeKuttaFehlberg78Coefficients.cCoefficients( 10 ) = 1.0;
    rungeKuttaFehlberg78Coefficients.cCoefficients( 11 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.cCoefficients( 12 ) = 1.0;

    // Define b-coefficients for the Runge-Kutta method of order 8
    // with an embedded 7th-order method for stepsize control and a total of 13 stages.
    rungeKuttaFehlberg78Coefficients.bCoefficients = Eigen::MatrixXd::Zero( 2, 13 );
    rungeKuttaFehlberg78Coefficients.bCoefficients( 0, 0 ) = 41.0 / 840.0;
    rungeKuttaFehlberg78Coefficients.bCoefficients( 0, 5 ) = 34.0 / 105.0;
    rungeKuttaFehlberg78Coefficients.bCoefficients( 0, 6 ) = 9.0 / 35.0;
    rungeKuttaFehlberg78Coefficients.bCoefficients( 0, 7 ) =
            rungeKuttaFehlberg78Coefficients.bCoefficients( 0, 6 );
    rungeKuttaFehlberg78Coefficients.bCoefficients( 0, 8 ) = 9.0 / 280.0;
    rungeKuttaFehlberg78Coefficients.bCoefficients( 0, 9 ) =
            rungeKuttaFehlberg78Coefficients.bCoefficients( 0, 8 );
    rungeKuttaFehlberg78Coefficients.bCoefficients( 0, 10 ) = 41.0 / 840.0;

    rungeKuttaFehlberg78Coefficients.bCoefficients( 1, 5 ) = 34.0 / 105.0;
    rungeKuttaFehlberg78Coefficients.bCoefficients( 1, 6 ) = 9.0 / 35.0;
    rungeKuttaFehlberg78Coefficients.bCoefficients( 1, 7 ) =
            rungeKuttaFehlberg78Coefficients.bCoefficients( 1, 6 );
    rungeKuttaFehlberg78Coefficients.bCoefficients( 1, 8 ) = 9.0 / 280.0;
    rungeKuttaFehlberg78Coefficients.bCoefficients( 1, 9 ) =
            rungeKuttaFehlberg78Coefficients.bCoefficients( 1, 8 );
    rungeKuttaFehlberg78Coefficients.bCoefficients( 1, 11 ) = 41.0 / 840.0;
    rungeKuttaFehlberg78Coefficients.bCoefficients( 1, 12 ) =
            rungeKuttaFehlberg78Coefficients.bCoefficients( 1, 11 );
}

//! Initialize RK87 (Dormand and Prince) coefficients.
void initializerungeKutta87DormandPrinceCoefficients(
        RungeKuttaCoefficients& rungeKutta87DormandPrinceCoefficients )
{
    // Define characteristics of coefficient set.
    rungeKutta87DormandPrinceCoefficients.lowerOrder = 7;
    rungeKutta87DormandPrinceCoefficients.higherOrder = 8;
    rungeKutta87DormandPrinceCoefficients.orderEstimateToIntegrate
            = RungeKuttaCoefficients::higher;

    // This coefficient set is taken from (Montenbruck and Gill, 2005).

    // a-coefficients for the Runge-Kutta method of order 8
    // with an embedded 7th-order method for stepsize control and a total of 13 stages.
    rungeKutta87DormandPrinceCoefficients.aCoefficients = Eigen::MatrixXd::Zero( 13, 12 );

    rungeKutta87DormandPrinceCoefficients.aCoefficients( 1, 0 ) = 1.0 / 18.0;

    rungeKutta87DormandPrinceCoefficients.aCoefficients( 2, 0 ) = 1.0 / 48.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 2, 1 ) = 1.0 / 16.0;

    rungeKutta87DormandPrinceCoefficients.aCoefficients( 3, 0 ) = 1.0 / 32.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 3, 2 ) = 3.0 / 32.0;

    rungeKutta87DormandPrinceCoefficients.aCoefficients( 4, 0 ) = 5.0 / 16.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 4, 2 ) = -75.0 / 64.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 4, 3 ) = 75.0 / 64.0;

    rungeKutta87DormandPrinceCoefficients.aCoefficients( 5, 0 ) = 3.0 / 80.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 5, 3 ) = 3.0 / 16.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 5, 4 ) = 3.0 / 20.0;

    rungeKutta87DormandPrinceCoefficients.aCoefficients( 6, 0 ) = 29443841.0 / 614563906.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 6, 3 ) = 77736538.0 / 692538347.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 6, 4 ) = -28693883.0 / 1125000000.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 6, 5 ) = 23124283.0 / 1800000000.0;

    rungeKutta87DormandPrinceCoefficients.aCoefficients( 7, 0 ) = 16016141.0 / 946692911.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 7, 3 ) = 61564180.0 / 158732637.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 7, 4 ) = 22789713.0 / 633445777.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 7, 5 ) = 545815736.0 / 2771057229.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 7, 6 ) = -180193667.0 / 1043307555.0;

    rungeKutta87DormandPrinceCoefficients.aCoefficients( 8, 0 ) = 39632708.0 / 573591083.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 8, 3 ) = -433636366.0 / 683701615.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 8, 4 ) = -421739975.0 / 2616292301.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 8, 5 ) = 100302831.0 / 723423059.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 8, 6 ) = 790204164.0 / 839813087.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 8, 7 ) = 800635310.0 / 3783071287.0;

    rungeKutta87DormandPrinceCoefficients.aCoefficients( 9, 0 ) = 246121993.0 / 1340847787.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 9, 3 ) = -37695042795.0 / 15268766246.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 9, 4 ) = -309121744.0 / 1061227803.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 9, 5 ) = -12992083.0 / 490766935.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 9, 6 ) = 6005943493.0 / 2108947869.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 9, 7 ) = 393006217.0 / 1396673457.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 9, 8 ) = 123872331.0 / 1001029789.0;

    rungeKutta87DormandPrinceCoefficients.aCoefficients( 10, 0 ) = -1028468189.0 / 846180014.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 10, 3 ) = 8478235783.0 / 508512852.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 10, 4 ) = 1311729495.0 / 1432422823.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 10, 5 ) = -10304129995.0 / 1701304382.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 10, 6 ) = -48777925059.0 / 3047939560.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 10, 7 ) = 15336726248.0 / 1032824649.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 10, 8 ) = -45442868181.0 / 3398467696.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 10, 9 ) = 3065993473.0 / 597172653.0;

    rungeKutta87DormandPrinceCoefficients.aCoefficients( 11, 0 ) = 185892177.0 / 718116043.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 11, 3 ) = -3185094517.0 / 667107341.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 11, 4 ) = -477755414.0 / 1098053517.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 11, 5 ) = -703635378.0 / 230739211.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 11, 6 ) = 5731566787.0 / 1027545527.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 11, 7 ) = 5232866602.0 / 850066563.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 11, 8 ) = -4093664535.0 / 808688257.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 11, 9 ) = 3962137247.0 / 1805957418.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 11, 10 ) = 65686358.0 / 487910083.0;

    rungeKutta87DormandPrinceCoefficients.aCoefficients( 12, 0 ) = 403863854.0 / 491063109.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 12, 3 ) = -5068492393.0 / 434740067.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 12, 4 ) = -411421997.0 / 543043805.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 12, 5 ) = 652783627.0 / 914296604.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 12, 6 ) = 11173962825.0 / 925320556.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 12, 7 ) = -13158990841.0 / 6184727034.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 12, 8 ) = 3936647629.0 / 1978049680.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 12, 9 ) = -160528059.0 / 685178525.0;
    rungeKutta87DormandPrinceCoefficients.aCoefficients( 12, 10 ) = 248638103.0 / 1413531060.0;

    // c-coefficients for the Runge-Kutta method of order 8
    // with an embedded 7th-order method for stepsize control and a total of 13 stages.
    rungeKutta87DormandPrinceCoefficients.cCoefficients = Eigen::VectorXd::Zero( 13 );

    rungeKutta87DormandPrinceCoefficients.cCoefficients( 1 ) = 1.0 / 18.0;
    rungeKutta87DormandPrinceCoefficients.cCoefficients( 2 ) = 1.0 / 12.0;
    rungeKutta87DormandPrinceCoefficients.cCoefficients( 3 ) = 1.0 / 8.0;
    rungeKutta87DormandPrinceCoefficients.cCoefficients( 4 ) = 5.0 / 16.0;
    rungeKutta87DormandPrinceCoefficients.cCoefficients( 5 ) = 3.0 / 8.0;
    rungeKutta87DormandPrinceCoefficients.cCoefficients( 6 ) = 59.0 / 400.0;
    rungeKutta87DormandPrinceCoefficients.cCoefficients( 7 ) = 93.0 / 200.0;
    rungeKutta87DormandPrinceCoefficients.cCoefficients( 8 ) = 5490023248.0 / 9719169821.0;
    rungeKutta87DormandPrinceCoefficients.cCoefficients( 9 ) = 13.0 / 20.0;
    rungeKutta87DormandPrinceCoefficients.cCoefficients( 10 ) = 1201146811.0 / 1299019798.0;
    rungeKutta87DormandPrinceCoefficients.cCoefficients( 11 ) = 1.0;
    rungeKutta87DormandPrinceCoefficients.cCoefficients( 12 ) = 1.0;

    // b-coefficients for the Runge-Kutta method of order 8
    // with an embedded 7th-order method for stepsize control and a total of 13 stages.
    rungeKutta87DormandPrinceCoefficients.bCoefficients = Eigen::MatrixXd::Zero( 2, 13 );

    rungeKutta87DormandPrinceCoefficients.bCoefficients( 0, 0 ) = 13451932.0 / 455176623.0;
    rungeKutta87DormandPrinceCoefficients.bCoefficients( 0, 5 ) = -808719846.0 / 976000145.0;
    rungeKutta87DormandPrinceCoefficients.bCoefficients( 0, 6 ) = 1757004468.0 / 5645159321.0;
    rungeKutta87DormandPrinceCoefficients.bCoefficients( 0, 7 ) = 656045339.0 / 265891186.0;
    rungeKutta87DormandPrinceCoefficients.bCoefficients( 0, 8 ) = -3867574721.0 / 1518517206.0;
    rungeKutta87DormandPrinceCoefficients.bCoefficients( 0, 9 ) = 465885868.0 / 322736535.0;
    rungeKutta87DormandPrinceCoefficients.bCoefficients( 0, 10 ) = 53011238.0 / 667516719.0;
    rungeKutta87DormandPrinceCoefficients.bCoefficients( 0, 11 ) = 2.0 / 45.0;

    rungeKutta87DormandPrinceCoefficients.bCoefficients( 1, 0 ) = 14005451.0 / 335480064.0;
    rungeKutta87DormandPrinceCoefficients.bCoefficients( 1, 5 ) = -59238493.0 / 1068277825.0;
    rungeKutta87DormandPrinceCoefficients.bCoefficients( 1, 6 ) = 181606767.0 / 758867731.0;
    rungeKutta87DormandPrinceCoefficients.bCoefficients( 1, 7 ) = 561292985.0 / 797845732.0;
    rungeKutta87DormandPrinceCoefficients.bCoefficients( 1, 8 ) = -1041891430.0 / 1371343529.0;
    rungeKutta87DormandPrinceCoefficients.bCoefficients( 1, 9 ) = 760417239.0 / 1151165299.0;
    rungeKutta87DormandPrinceCoefficients.bCoefficients( 1, 10 ) = 118820643.0 / 751138087.0;
    rungeKutta87DormandPrinceCoefficients.bCoefficients( 1, 11 ) = -528747749.0 / 2220607170.0;
    rungeKutta87DormandPrinceCoefficients.bCoefficients( 1, 12 ) = 1.0 / 4.0;
}

//! Get coefficients for a specified coefficient set
const RungeKuttaCoefficients& RungeKuttaCoefficients::get(
        RungeKuttaCoefficients::CoefficientSets coefficientSet )
{
    static RungeKuttaCoefficients rungeKuttaFehlberg45Coefficients,
                                  rungeKuttaFehlberg56Coefficients,
                                  rungeKuttaFehlberg78Coefficients,
                                  rungeKutta87DormandPrinceCoefficients;

    switch ( coefficientSet )
    {
    case rungeKuttaFehlberg45:
        if ( rungeKuttaFehlberg45Coefficients.higherOrder != 5 )
        {
            initializeRungeKuttaFehlberg45Coefficients( rungeKuttaFehlberg45Coefficients );
        }
        return rungeKuttaFehlberg45Coefficients;

    case rungeKuttaFehlberg56:
        if ( rungeKuttaFehlberg56Coefficients.higherOrder != 6 )
        {
            initializeRungeKuttaFehlberg56Coefficients( rungeKuttaFehlberg56Coefficients );
        }
        return rungeKuttaFehlberg56Coefficients;

    case rungeKuttaFehlberg78:
        if ( rungeKuttaFehlberg78Coefficients.higherOrder != 8 )
        {
            initializeRungeKuttaFehlberg78Coefficients( rungeKuttaFehlberg78Coefficients );
        }
        return rungeKuttaFehlberg78Coefficients;

    case rungeKutta87DormandPrince:
        if ( rungeKutta87DormandPrinceCoefficients.higherOrder != 7 )
        {
            initializerungeKutta87DormandPrinceCoefficients(
                        rungeKutta87DormandPrinceCoefficients );
        }
        return rungeKutta87DormandPrinceCoefficients;

    default: // The default case will never occur because CoefficientsSet is an enum.
        throw RungeKuttaCoefficients( );
    }
}

} // namespace numerical_integrators
} // namespace tudat
