/*    Copyright (c) 2010-2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120203    B. Tong Minh      File created
 *
 *    References
 *      Burden, R.L., Faires, J.D. Numerical Analysis, 7th Edition, Books/Cole, 2001.
 *
 */

#include <Eigen/Core>
#include "Tudat/Mathematics/NumericalIntegrators/rungeKuttaCoefficients.h"

namespace tudat
{
namespace mathematics
{
namespace numerical_integrators
{

//! Initialize RKF45 coefficients.
void initializeRungeKuttaFehlberg45Coefficients( RungeKuttaCoefficients&
                                                 rungeKuttaFehlberg45Coefficients )
{
    rungeKuttaFehlberg45Coefficients.higherOrder = 5;
    rungeKuttaFehlberg45Coefficients.lowerOrder = 4;

    // Define a-coefficients for the Runge-Kutta-Fehlberg method of order 5
    // with an embedded 4th-order method for stepsize cont rol and a total of 6 stages.
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
    rungeKuttaFehlberg45Coefficients.cCoefficients( 0 ) = 0.0;
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
    rungeKuttaFehlberg56Coefficients.higherOrder = 6;
    rungeKuttaFehlberg56Coefficients.lowerOrder = 5;

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
    rungeKuttaFehlberg56Coefficients.aCoefficients( 6, 1 ) =   0.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 6, 2 ) =  11.0 / 256.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 6, 3 ) = -11.0 / 160.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 6, 4 ) =  11.0 / 256.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 6, 5 ) =   0.0;

    rungeKuttaFehlberg56Coefficients.aCoefficients( 7, 0 ) = 93.0 / 640.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 7, 1 ) = -18.0 / 5.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 7, 2 ) = 803.0 / 256.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 7, 3 ) = -11.0 / 160.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 7, 4 ) = 99.0 / 256.0;
    rungeKuttaFehlberg56Coefficients.aCoefficients( 7, 5 ) = 0.0;
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
    rungeKuttaFehlberg78Coefficients.higherOrder = 8;
    rungeKuttaFehlberg78Coefficients.lowerOrder = 7;

    // Define a-coefficients for the Runge-Kutta-Fehlberg method of order 8
    // with an embedded 7th-order method for stepsize control and a total of 13 stages.
    rungeKuttaFehlberg78Coefficients.aCoefficients = Eigen::MatrixXd::Zero( 13, 12 );
    rungeKuttaFehlberg78Coefficients.aCoefficients( 1, 0 ) = 2.0 / 27.0;

    rungeKuttaFehlberg78Coefficients.aCoefficients( 2, 0 ) = 1.0 / 36.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 2, 1 ) = 1.0 / 12.0;

    rungeKuttaFehlberg78Coefficients.aCoefficients( 3, 0 ) = 1.0 / 24.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 3, 1 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 3, 2 ) = 1.0 / 8.0;

    rungeKuttaFehlberg78Coefficients.aCoefficients( 4, 0 ) = 5.0 / 12.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 4, 1 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 4, 2 ) = -25.0 / 16.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 4, 3 ) = 25.0 / 16.0;

    rungeKuttaFehlberg78Coefficients.aCoefficients( 5, 0 ) = 1.0 / 20.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 5, 1 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 5, 2 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 5, 3 ) = 1.0 / 4.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 5, 4 ) = 1.0 / 5.0;

    rungeKuttaFehlberg78Coefficients.aCoefficients( 6, 0 ) = -25.0 / 108.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 6, 1 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 6, 2 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 6, 3 ) = 125.0 / 108.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 6, 4 ) = -65.0 / 27.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 6, 5 ) = 125.0 / 54.0;

    rungeKuttaFehlberg78Coefficients.aCoefficients( 7, 0 ) = 31.0 / 300.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 7, 1 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 7, 2 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 7, 3 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 7, 4 ) = 61.0 / 225.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 7, 5 ) = -2.0 / 9.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 7, 6 ) = 13.0 / 900.0;

    rungeKuttaFehlberg78Coefficients.aCoefficients( 8, 0 ) = 2.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 8, 1 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 8, 2 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 8, 3 ) = -53.0 / 6.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 8, 4 ) = 704.0 / 45.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 8, 5 ) = -107.0 / 9.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 8, 6 ) = 67.0 / 90.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 8, 7 ) = 3.0;

    rungeKuttaFehlberg78Coefficients.aCoefficients( 9, 0 ) = -91.0 / 108.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 9, 1 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 9, 2 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 9, 3 ) = 23.0 / 108.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 9, 4 ) = -976.0 / 135.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 9, 5 ) = 311.0 / 54.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 9, 6 ) = -19.0 / 60.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 9, 7 ) = 17.0 / 6.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 9, 8 ) = -1.0 / 12.0;

    rungeKuttaFehlberg78Coefficients.aCoefficients( 10, 0 ) = 2383.0 / 4100.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 10, 1 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 10, 2 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 10, 3 ) = -341.0 / 164.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 10, 4 ) = 4496.0 / 1025.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 10, 5 ) = -301.0 / 82.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 10, 6 ) = 2133.0 / 4100.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 10, 7 ) = 45.0 / 82.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 10, 8 ) = 45.0 / 164.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 10, 9 ) = 18.0 / 41.0;

    rungeKuttaFehlberg78Coefficients.aCoefficients( 11, 0 ) = 3.0 / 205.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 11, 1 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 11, 2 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 11, 3 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 11, 4 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 11, 5 ) = -6.0 / 41.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 11, 6 ) = -3.0 / 205.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 11, 7 ) = -3.0 / 41.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 11, 8 ) = 3.0 / 41.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 11, 9 ) = 6.0 / 41.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 11, 10 ) = 0.0;

    rungeKuttaFehlberg78Coefficients.aCoefficients( 12, 0 ) = -1777.0 / 4100.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 12, 1 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 12, 2 ) = 0.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 12, 3 ) = -341.0 / 164.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 12, 4 ) = 4496.0 / 1025.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 12, 5 ) = -289.0 / 82.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 12, 6 ) = 2193.0 / 4100.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 12, 7 ) = 51.0 / 82.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 12, 8 ) = 33.0 / 164.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 12, 9 ) = 12.0 / 41.0;
    rungeKuttaFehlberg78Coefficients.aCoefficients( 12, 10 ) = 0.0;
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

//! Get coefficients for a specified coefficient set
const RungeKuttaCoefficients& RungeKuttaCoefficients::get(
        RungeKuttaCoefficients::CoefficientSets coefficientSet )
{
    static RungeKuttaCoefficients rungeKuttaFehlberg45Coefficients,
            rungeKuttaFehlberg56Coefficients, rungeKuttaFehlberg78Coefficients;

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

        default: // The default case will never occur because CoefficientsSet is an enum
            throw RungeKuttaCoefficients( );
    }
}

} // namespace integrators
} // namespace mathematics
} // namespace tudat
