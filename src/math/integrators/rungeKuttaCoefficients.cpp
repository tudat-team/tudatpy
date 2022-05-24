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

#include "tudat/math/integrators/rungeKuttaCoefficients.h"

namespace tudat
{
namespace numerical_integrators
{

//! Initialize forward Euler coefficients.
void initializeForwardEulerCoefficients( RungeKuttaCoefficients& forwardEulerCoefficients )
{
    // Define characteristics of coefficient set.
    forwardEulerCoefficients.isFixedStepSize = true;

    // Define a-coefficients for the forward Euler method of order 1 and with 1 stage.
    forwardEulerCoefficients.aCoefficients = Eigen::MatrixXd::Zero( 1, 1 );
    
    // Define c-coefficients for the forward Euler method of order 1 and with 1 stage.
    forwardEulerCoefficients.cCoefficients = Eigen::VectorXd::Zero( 1 );
    
    // Define b-coefficients for the forward Euler method of order 1 and with 1 stage.
    forwardEulerCoefficients.bCoefficients = Eigen::MatrixXd::Zero( 1, 1 );
    forwardEulerCoefficients.bCoefficients( 0, 0 ) = 1.0;

    // Set the name of these coefficients.
    forwardEulerCoefficients.name = "Forward Euler";
}

//! Initialize Heun-Euler coefficients.
void initializeHeunEulerCoefficients( RungeKuttaCoefficients& heunEulerCoefficients )
{    
    // Define characteristics of coefficient set.
    heunEulerCoefficients.lowerOrder = 1;
    heunEulerCoefficients.higherOrder = 2;
    heunEulerCoefficients.orderEstimateToIntegrate
            = RungeKuttaCoefficients::lower;
    
    // Define a-coefficients for the Heun-Euler method of order 2
    // with an embedded 1st-order method for stepsize control and a total of 2 stages.
    heunEulerCoefficients.aCoefficients = Eigen::MatrixXd::Zero( 2, 2 );
    heunEulerCoefficients.aCoefficients( 1, 0 ) = 1.0;
    
    // Define c-coefficients for the Heun-Euler method of order 2
    // with an embedded 1st-order method for stepsize control and a total of 2 stages.
    heunEulerCoefficients.cCoefficients = Eigen::VectorXd::Zero( 2 );
    heunEulerCoefficients.cCoefficients( 1 ) = 1.0;
    
    // Define b-coefficients for the Heun-Euler method of order 2
    // with an embedded 1st-order method for stepsize control and a total of 2 stages.
    heunEulerCoefficients.bCoefficients = Eigen::MatrixXd::Zero( 2, 2 );
    heunEulerCoefficients.bCoefficients( 0, 0 ) = 1.0;
    heunEulerCoefficients.bCoefficients( 1, 0 ) = 1.0 / 2.0;
    heunEulerCoefficients.bCoefficients( 1, 1 ) = 1.0 / 2.0;

    // Set the name of these coefficients.
    heunEulerCoefficients.name = "Heun-Euler";
}

//! Initialize RKF12 coefficients.
void initializeRungeKuttaFehlberg12Coefficients( RungeKuttaCoefficients& rungeKuttaFehlberg12Coefficients )
{    
    // Define characteristics of coefficient set.
    rungeKuttaFehlberg12Coefficients.lowerOrder = 1;
    rungeKuttaFehlberg12Coefficients.higherOrder = 2;
    rungeKuttaFehlberg12Coefficients.orderEstimateToIntegrate
            = RungeKuttaCoefficients::lower;
    
    // Define a-coefficients for the Runge-Kutta-Fehlberg method of order 1
    // with an embedded 2nd-order method for stepsize control and a total of 3 stages.
    rungeKuttaFehlberg12Coefficients.aCoefficients = Eigen::MatrixXd::Zero( 3, 3 );
    rungeKuttaFehlberg12Coefficients.aCoefficients( 1, 0 ) = 1.0 / 2.0;
    rungeKuttaFehlberg12Coefficients.aCoefficients( 2, 0 ) = 1.0 / 256.0;
    rungeKuttaFehlberg12Coefficients.aCoefficients( 2, 1 ) = 255.0 / 256.0;
    
    // Define c-coefficients for the Runge-Kutta-Fehlberg method of order 1
    // with an embedded 2nd-order method for stepsize control and a total of 3 stages.
    rungeKuttaFehlberg12Coefficients.cCoefficients = Eigen::VectorXd::Zero( 3 );
    rungeKuttaFehlberg12Coefficients.cCoefficients( 1 ) = 1.0 / 2.0;
    rungeKuttaFehlberg12Coefficients.cCoefficients( 2 ) = 1.0;
    
    // Define b-coefficients for the Runge-Kutta-Fehlberg method of order 1
    // with an embedded 2nd-order method for stepsize control and a total of 3 stages.
    rungeKuttaFehlberg12Coefficients.bCoefficients = Eigen::MatrixXd::Zero( 2, 3 );
    rungeKuttaFehlberg12Coefficients.bCoefficients( 0, 0 ) = 1.0 / 256.0;
    rungeKuttaFehlberg12Coefficients.bCoefficients( 0, 1 ) = 255.0 / 256.0;
    rungeKuttaFehlberg12Coefficients.bCoefficients( 1, 0 ) = 1.0 / 512.0;
    rungeKuttaFehlberg12Coefficients.bCoefficients( 1, 1 ) = 255.0 / 256.0;
    rungeKuttaFehlberg12Coefficients.bCoefficients( 1, 2 ) = 1.0 / 512.0;

    // Set the name of these coefficients.
    rungeKuttaFehlberg12Coefficients.name = "Runge-Kutta-Fehlberg 1/2";
}

//! Initialize RK4 coefficients.
void initializeRungeKutta4Coefficients( RungeKuttaCoefficients&
                                        rungeKutta4Coefficients )
{
    // Define characteristics of coefficient set.
    rungeKutta4Coefficients.isFixedStepSize = true;

    // Define a-coefficients for the Runge-Kutta method of order 4 and with 4 stages.
    rungeKutta4Coefficients.aCoefficients = Eigen::MatrixXd::Zero( 4, 4 );
    rungeKutta4Coefficients.aCoefficients( 1, 0 ) = 1.0 / 2.0;

    rungeKutta4Coefficients.aCoefficients( 2, 1 ) = 1.0 / 2.0;

    rungeKutta4Coefficients.aCoefficients( 3, 2 ) = 1.0;


    // Define c-coefficients for the Runge-Kutta method of order 4 and with 4 stages.
    rungeKutta4Coefficients.cCoefficients = Eigen::VectorXd::Zero( 4 );
    rungeKutta4Coefficients.cCoefficients( 1 ) = 1.0 / 2.0;
    rungeKutta4Coefficients.cCoefficients( 2 ) = 1.0 / 2.0;
    rungeKutta4Coefficients.cCoefficients( 3 ) = 1.0;


    // Define b-coefficients for the Runge-Kutta method of order 4 and with 4 stages.
    rungeKutta4Coefficients.bCoefficients = Eigen::MatrixXd::Zero( 1, 4 );
    rungeKutta4Coefficients.bCoefficients( 0, 0 ) = 1.0 / 6.0;
    rungeKutta4Coefficients.bCoefficients( 0, 1 ) = 2.0 / 6.0;
    rungeKutta4Coefficients.bCoefficients( 0, 2 ) = 2.0 / 6.0;
    rungeKutta4Coefficients.bCoefficients( 0, 3 ) = 1.0 / 6.0;

    // Set the name of these coefficients.
    rungeKutta4Coefficients.name = "Runge-Kutta 4";
}

//! Initialize explicit midpoint coefficients.
void initializeExplicitMidpointCoefficients( RungeKuttaCoefficients&
                                             explicitMidpointCoefficients )
{
    // Define characteristics of coefficient set.
    explicitMidpointCoefficients.isFixedStepSize = true;

    // Define a-coefficients for the explicit midpoint method of order 2 and with 2 stages.
    explicitMidpointCoefficients.aCoefficients = Eigen::MatrixXd::Zero( 2, 2 );
    explicitMidpointCoefficients.aCoefficients( 1, 0 ) = 1.0 / 2.0;

    // Define c-coefficients for the explicit midpoint method of order 2 and with 2 stages.
    explicitMidpointCoefficients.cCoefficients = Eigen::VectorXd::Zero( 2 );
    explicitMidpointCoefficients.cCoefficients( 1 ) = 1.0 / 2.0;

    // Define b-coefficients for the explicit midpoint method of order 2 and with 2 stages.
    explicitMidpointCoefficients.bCoefficients = Eigen::MatrixXd::Zero( 1, 2 );
    explicitMidpointCoefficients.bCoefficients( 0, 1 ) = 1.0;

    // Set the name of these coefficients.
    explicitMidpointCoefficients.name = "Explicit Midpoint";
}

//! Initialize explicit trapezoid coefficients.
void initializeExplicitTrapezoidRuleCoefficients( RungeKuttaCoefficients&
                                              explicitTrapezoidRuleCoefficients )
{
    // Define characteristics of coefficient set.
    explicitTrapezoidRuleCoefficients.isFixedStepSize = true;

    // Define a-coefficients for the explicit trapezoidal method of order 2 and with 2 stages.
    explicitTrapezoidRuleCoefficients.aCoefficients = Eigen::MatrixXd::Zero( 2, 2 );
    explicitTrapezoidRuleCoefficients.aCoefficients( 1, 0 ) = 1.0;

    // Define c-coefficients for the explicit trapezoidal method of order 2 and with 2 stages.
    explicitTrapezoidRuleCoefficients.cCoefficients = Eigen::VectorXd::Zero( 2 );
    explicitTrapezoidRuleCoefficients.cCoefficients( 1 ) = 1.0;

    // Define b-coefficients for the explicit trapezoidal method of order 2 and with 2 stages.
    explicitTrapezoidRuleCoefficients.bCoefficients = Eigen::MatrixXd::Zero( 1, 2 );
    explicitTrapezoidRuleCoefficients.bCoefficients( 0, 0 ) = 1.0 / 2.0;
    explicitTrapezoidRuleCoefficients.bCoefficients( 0, 1 ) = 1.0 / 2.0;

    // Set the name of these coefficients.
    explicitTrapezoidRuleCoefficients.name = "Explicit Trapezoid Rule";
}

//! Initialize Ralston coefficients.
void initializeRalstonCoefficients( RungeKuttaCoefficients&
                                    ralstonCoefficients )
{
    // Define characteristics of coefficient set.
    ralstonCoefficients.isFixedStepSize = true;

    // Define a-coefficients for the Ralston method of order 2 and with 2 stages.
    ralstonCoefficients.aCoefficients = Eigen::MatrixXd::Zero( 2, 2 );
    ralstonCoefficients.aCoefficients( 1, 0 ) = 2.0 / 3.0;

    // Define c-coefficients for the Ralston method of order 2 and with 2 stages.
    ralstonCoefficients.cCoefficients = Eigen::VectorXd::Zero( 2 );
    ralstonCoefficients.cCoefficients( 1 ) = 2.0 / 3.0;

    // Define b-coefficients for the Ralston method of order 2 and with 2 stages.
    ralstonCoefficients.bCoefficients = Eigen::MatrixXd::Zero( 1, 2 );
    ralstonCoefficients.bCoefficients( 0, 0 ) = 1.0 / 4.0;
    ralstonCoefficients.bCoefficients( 0, 1 ) = 3.0 / 4.0;

    // Set the name of these coefficients.
    ralstonCoefficients.name = "Ralston";
}

//! Initialize RK3 coefficients.
void initializeRungeKutta3Coefficients( RungeKuttaCoefficients&
                                rk3Coefficients )
{
    // Define characteristics of coefficient set.
    rk3Coefficients.isFixedStepSize = true;

    // Define a-coefficients for the RK3 method of order 3 and with 3 stages.
    rk3Coefficients.aCoefficients = Eigen::MatrixXd::Zero( 3, 3 );
    rk3Coefficients.aCoefficients( 1, 0 ) = 1.0 / 2.0;
    rk3Coefficients.aCoefficients( 2, 0 ) = -1.0;
    rk3Coefficients.aCoefficients( 2, 1 ) = 2.0;

    // Define c-coefficients for the RK3 method of order 3 and with 3 stages.
    rk3Coefficients.cCoefficients = Eigen::VectorXd::Zero( 3 );
    rk3Coefficients.cCoefficients( 1 ) = 1.0 / 2.0;
    rk3Coefficients.cCoefficients( 2 ) = 1.0;

    // Define b-coefficients for the RK3 method of order 3 and with 3 stages.
    rk3Coefficients.bCoefficients = Eigen::MatrixXd::Zero( 1, 3 );
    rk3Coefficients.bCoefficients( 0, 0 ) = 1.0 / 6.0;
    rk3Coefficients.bCoefficients( 0, 1 ) = 2.0 / 3.0;
    rk3Coefficients.bCoefficients( 0, 2 ) = 1.0 / 6.0;

    // Set the name of these coefficients.
    rk3Coefficients.name = "Runge-Kutta 3";
}

//! Initialize Ralston3 coefficients.
void initializeRalston3Coefficients( RungeKuttaCoefficients&
                                     ralston3Coefficients )
{
    // Define characteristics of coefficient set.
    ralston3Coefficients.isFixedStepSize = true;

    // Define a-coefficients for the Ralston3 method of order 3 and with 3 stages.
    ralston3Coefficients.aCoefficients = Eigen::MatrixXd::Zero( 3, 3 );
    ralston3Coefficients.aCoefficients( 1, 0 ) = 1.0 / 2.0;
    ralston3Coefficients.aCoefficients( 2, 1 ) = 3.0 / 4.0;

    // Define c-coefficients for the Ralston3 method of order 3 and with 3 stages.
    ralston3Coefficients.cCoefficients = Eigen::VectorXd::Zero( 3 );
    ralston3Coefficients.cCoefficients( 1 ) = 1.0 / 2.0;
    ralston3Coefficients.cCoefficients( 2 ) = 3.0 / 4.0;

    // Define b-coefficients for the Ralston3 method of order 3 and with 3 stages.
    ralston3Coefficients.bCoefficients = Eigen::MatrixXd::Zero( 1, 3 );
    ralston3Coefficients.bCoefficients( 0, 0 ) = 2.0 / 9.0;
    ralston3Coefficients.bCoefficients( 0, 1 ) = 1.0 / 3.0;
    ralston3Coefficients.bCoefficients( 0, 2 ) = 4.0 / 9.0;

    // Set the name of these coefficients.
    ralston3Coefficients.name = "Ralston 3";
}

//! Initialize SSPRK3 coefficients.
void initializeSSPRK3Coefficients( RungeKuttaCoefficients&
                                   ssprk3Coefficients )
{
    // Define characteristics of coefficient set.
    ssprk3Coefficients.isFixedStepSize = true;

    // Define a-coefficients for the SSPRK3 method of order 3 and with 3 stages.
    ssprk3Coefficients.aCoefficients = Eigen::MatrixXd::Zero( 3, 3 );
    ssprk3Coefficients.aCoefficients( 1, 0 ) = 1.0;
    ssprk3Coefficients.aCoefficients( 2, 0 ) = 1.0 / 4.0;
    ssprk3Coefficients.aCoefficients( 2, 1 ) = 1.0 / 4.0;

    // Define c-coefficients for the SSPRK3 method of order 3 and with 3 stages.
    ssprk3Coefficients.cCoefficients = Eigen::VectorXd::Zero( 3 );
    ssprk3Coefficients.cCoefficients( 1 ) = 1.0 ;
    ssprk3Coefficients.cCoefficients( 2 ) = 1.0 / 2.0;

    // Define b-coefficients for the SSPRK3 method of order 3 and with 3 stages.
    ssprk3Coefficients.bCoefficients = Eigen::MatrixXd::Zero( 1, 3 );
    ssprk3Coefficients.bCoefficients( 0, 0 ) = 1.0 / 6.0;
    ssprk3Coefficients.bCoefficients( 0, 1 ) = 1.0 / 6.0;
    ssprk3Coefficients.bCoefficients( 0, 2 ) = 2.0 / 3.0;

    // Set the name of these coefficients.
    ssprk3Coefficients.name = "Strong Stability Preserving Runge-Kutta 3";
}

//! Initialize Ralston 4 coefficients.
void initializeRalston4Coefficients( RungeKuttaCoefficients&
                                     ralston4Coefficients )
{
    // Define characteristics of coefficient set.
    ralston4Coefficients.isFixedStepSize = true;

    // Compute the square root of 5 to re-use later on.
    double SQRT5 = std::sqrt( 5.0 );

    // Define a-coefficients for the Ralston4 method of order 4 and with 4 stages.
    ralston4Coefficients.aCoefficients = Eigen::MatrixXd::Zero( 4, 4 );
    ralston4Coefficients.aCoefficients( 1, 0 ) = 2.0 / 5.0;
    ralston4Coefficients.aCoefficients( 2, 0 ) = (-2889.0 + 1428.0 * SQRT5) / 1024.0;
    ralston4Coefficients.aCoefficients( 2, 1 ) = (3785.0 - 1620.0 * SQRT5) / 1024.0;
    ralston4Coefficients.aCoefficients( 3, 0 ) = (-3365.0 + 2094.0 * SQRT5) / 6040.0;
    ralston4Coefficients.aCoefficients( 3, 1 ) = (-975.0 - 3046.0 * SQRT5) / 2552.0;
    ralston4Coefficients.aCoefficients( 3, 2 ) = (467040.0 + 203968.0 * SQRT5) / 240845.0;

    // Define c-coefficients for the Ralston4 method of order 4 and with 4 stages.
    ralston4Coefficients.cCoefficients = Eigen::VectorXd::Zero( 4 );
    ralston4Coefficients.cCoefficients( 1 ) = 2.0 / 5.0;
    ralston4Coefficients.cCoefficients( 2 ) = (14.0 - 3.0 * SQRT5) / 16.0;
    ralston4Coefficients.cCoefficients( 3 ) = 1.0;

    // Define b-coefficients for the Ralston4 method of order 4 and with 4 stages.
    ralston4Coefficients.bCoefficients = Eigen::MatrixXd::Zero( 1, 4 );
    ralston4Coefficients.bCoefficients( 0, 0 ) = (263.0 + 24.0 * SQRT5) / 1812.0;
    ralston4Coefficients.bCoefficients( 0, 1 ) = (125.0 - 1000.0 * SQRT5) / 3828.0;
    ralston4Coefficients.bCoefficients( 0, 2 ) = 1024.0 * (3346.0 + 1623.0 * SQRT5) / 5924787.0;
    ralston4Coefficients.bCoefficients( 0, 3 ) = (30.0 - 4.0 * SQRT5) / 123.0;

    // Set the name of these coefficients.
    ralston4Coefficients.name = "Ralston 4";
}

//! Initialize 3/8-rule RK4 coefficients.
void initializeThreeEighthRuleRK4Coefficients( RungeKuttaCoefficients&
                                               threeEighthRuleRK4Coefficients )
{
    // Define characteristics of coefficient set.
    threeEighthRuleRK4Coefficients.isFixedStepSize = true;

    // Define a-coefficients for the 3/8-rule RK4 method of order 4 and with 4 stages.
    threeEighthRuleRK4Coefficients.aCoefficients = Eigen::MatrixXd::Zero( 4, 4 );
    threeEighthRuleRK4Coefficients.aCoefficients( 1, 0 ) = 1.0 / 3.0;
    threeEighthRuleRK4Coefficients.aCoefficients( 2, 0 ) = -1.0 / 3.0;
    threeEighthRuleRK4Coefficients.aCoefficients( 2, 1 ) = 1.0;
    threeEighthRuleRK4Coefficients.aCoefficients( 3, 0 ) = 1.0;
    threeEighthRuleRK4Coefficients.aCoefficients( 3, 1 ) = -1.0;
    threeEighthRuleRK4Coefficients.aCoefficients( 3, 2 ) = 1.0;

    // Define c-coefficients for the 3/8-rule RK4 method of order 4 and with 4 stages.
    threeEighthRuleRK4Coefficients.cCoefficients = Eigen::VectorXd::Zero( 4 );
    threeEighthRuleRK4Coefficients.cCoefficients( 1 ) = 1.0 / 3.0;
    threeEighthRuleRK4Coefficients.cCoefficients( 2 ) = 2.0 / 3.0;
    threeEighthRuleRK4Coefficients.cCoefficients( 3 ) = 1.0;

    // Define b-coefficients for the 3/8-rule RK4 method of order 4 and with 4 stages.
    threeEighthRuleRK4Coefficients.bCoefficients = Eigen::MatrixXd::Zero( 1, 4 );
    threeEighthRuleRK4Coefficients.bCoefficients( 0, 0 ) = 1.0 / 8.0;
    threeEighthRuleRK4Coefficients.bCoefficients( 0, 1 ) = 3.0 / 8.0;
    threeEighthRuleRK4Coefficients.bCoefficients( 0, 2 ) = 3.0 / 8.0;
    threeEighthRuleRK4Coefficients.bCoefficients( 0, 3 ) = 1.0 / 8.0;

    // Set the name of these coefficients.
    threeEighthRuleRK4Coefficients.name = "3/8-rule RK4";
}

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

    //// COMMENTED OUT BELOW : alternative coefficients that give slightly different results
    //// but give orders closer to what is expected (4.6~4 and 5=5).

    // // Define a-coefficients for the Runge-Kutta-Fehlberg method of order 5
    // // with an embedded 4th-order method for stepsize control and a total of 6 stages.
    // rungeKuttaFehlberg45Coefficients.aCoefficients = Eigen::MatrixXd::Zero( 6, 5 );
    // rungeKuttaFehlberg45Coefficients.aCoefficients( 1, 0 ) = 2.0 / 9.0;

    // rungeKuttaFehlberg45Coefficients.aCoefficients( 2, 0 ) = 1.0 / 12.0;
    // rungeKuttaFehlberg45Coefficients.aCoefficients( 2, 1 ) = 1.0 / 4.0;

    // rungeKuttaFehlberg45Coefficients.aCoefficients( 3, 0 ) =  69.0 / 128.0;
    // rungeKuttaFehlberg45Coefficients.aCoefficients( 3, 1 ) = -243.0 / 128.0;
    // rungeKuttaFehlberg45Coefficients.aCoefficients( 3, 2 ) =  135.0 / 64.0;

    // rungeKuttaFehlberg45Coefficients.aCoefficients( 4, 0 ) = -17.0 / 12.0;
    // rungeKuttaFehlberg45Coefficients.aCoefficients( 4, 1 ) = 27.0 / 4.0;
    // rungeKuttaFehlberg45Coefficients.aCoefficients( 4, 2 ) = -27.0 / 5.0;
    // rungeKuttaFehlberg45Coefficients.aCoefficients( 4, 3 ) = 16.0 / 15.0;

    // rungeKuttaFehlberg45Coefficients.aCoefficients( 5, 0 ) = 65.0 / 432.0;
    // rungeKuttaFehlberg45Coefficients.aCoefficients( 5, 1 ) = -5.0 / 16.0;
    // rungeKuttaFehlberg45Coefficients.aCoefficients( 5, 2 ) = 13.0 / 16.0;
    // rungeKuttaFehlberg45Coefficients.aCoefficients( 5, 3 ) = 4.0 / 27.0;
    // rungeKuttaFehlberg45Coefficients.aCoefficients( 5, 4 ) = 5.0 / 144.0;


    // // Define c-coefficients for the Runge-Kutta-Fehlberg method of order 5
    // // with an embedded 4th-order method for stepsize control and a total of 6 stages.
    // rungeKuttaFehlberg45Coefficients.cCoefficients = Eigen::VectorXd::Zero( 6 );
    // rungeKuttaFehlberg45Coefficients.cCoefficients( 1 ) = 2.0 / 9.0;
    // rungeKuttaFehlberg45Coefficients.cCoefficients( 2 ) = 1.0 / 3.0;
    // rungeKuttaFehlberg45Coefficients.cCoefficients( 3 ) = 3.0 / 4.0;
    // rungeKuttaFehlberg45Coefficients.cCoefficients( 4 ) = 1.0;
    // rungeKuttaFehlberg45Coefficients.cCoefficients( 5 ) = 5.0 / 6.0;


    // // Define b-coefficients for the Runge-Kutta method of order 5
    // // with an embedded 4th-order method for stepsize control and a total of 6 stages.
    // rungeKuttaFehlberg45Coefficients.bCoefficients = Eigen::MatrixXd::Zero( 2, 6 );
    // rungeKuttaFehlberg45Coefficients.bCoefficients( 0, 0 ) = 1.0 / 9.0;
    // rungeKuttaFehlberg45Coefficients.bCoefficients( 0, 2 ) = 9.0 / 20.0;
    // rungeKuttaFehlberg45Coefficients.bCoefficients( 0, 3 ) = 16.0 / 45.0;
    // rungeKuttaFehlberg45Coefficients.bCoefficients( 0, 4 ) = 1.0 / 12.0;

    // rungeKuttaFehlberg45Coefficients.bCoefficients( 1, 0 ) = 47.0 / 450.0;
    // rungeKuttaFehlberg45Coefficients.bCoefficients( 1, 2 ) = 12.0 / 25.0;
    // rungeKuttaFehlberg45Coefficients.bCoefficients( 1, 3 ) = 32.0 / 225.0;
    // rungeKuttaFehlberg45Coefficients.bCoefficients( 1, 4 ) = 1.0 / 30.0;
    // rungeKuttaFehlberg45Coefficients.bCoefficients( 1, 5 ) = 6.0 / 25.0;

    // Set the name of these coefficients.
    rungeKuttaFehlberg45Coefficients.name = "Runge-Kutta-Fehlberg 4/5";
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

    // Set the name of these coefficients.
    rungeKuttaFehlberg56Coefficients.name = "Runge-Kutta-Fehlberg 5/6";
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

    // Set the name of these coefficients.
    rungeKuttaFehlberg78Coefficients.name = "Runge-Kutta-Fehlberg 7/8";
}

//! Initialize RK87 (Dormand and Prince) coefficients.
void initializeRungeKutta87DormandPrinceCoefficients(
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

    // Set the name of these coefficients.
    rungeKutta87DormandPrinceCoefficients.name = "Runge-Kutta 8/7 Dormand-Prince";
}

//! Initialize RKF89 coefficients.
void initializeRungeKuttaFehlberg89Coefficients( RungeKuttaCoefficients&
                                                 rungeKuttaFehlberg89Coefficients )
{
    // Define characteristics of coefficient set.
    rungeKuttaFehlberg89Coefficients.lowerOrder = 8;
    rungeKuttaFehlberg89Coefficients.higherOrder = 9;
    rungeKuttaFehlberg89Coefficients.orderEstimateToIntegrate = RungeKuttaCoefficients::lower;

    // This coefficient set is taken from
    // Fehlberg, E. (1968).
    // Classical fifth-, sixth-, seventh-, and eighth-order Runge-Kutta formulas with stepsize control.
    // National Aeronautics and Space Administration.

    // Define a-coefficients for the Runge-Kutta method of order 8 and a total of 17 stages.
    rungeKuttaFehlberg89Coefficients.aCoefficients = Eigen::MatrixXd::Zero( 17, 16 );

    // Assign all the a coefficients values
    rungeKuttaFehlberg89Coefficients.aCoefficients( 1, 0 ) = 0.44368940376498183109599404281370;

    rungeKuttaFehlberg89Coefficients.aCoefficients( 2, 0 ) = 0.16638352641186818666099776605514;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 2, 1 ) = 0.49915057923560455998299329816541;

    rungeKuttaFehlberg89Coefficients.aCoefficients( 3, 0 ) = 0.24957528961780227999149664908271;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 3, 2 ) = 0.74872586885340683997448994724812;

    rungeKuttaFehlberg89Coefficients.aCoefficients( 4, 0 ) = 0.20661891163400602426556710393185;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 4, 2 ) = 0.17707880377986347040380997288319;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 4, 3 ) = -0.68197715413869494669377076815048e-1;

    rungeKuttaFehlberg89Coefficients.aCoefficients( 5, 0 ) = 0.10927823152666408227903890926157;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 5, 3 ) = 0.40215962642367995421990563690087e-2;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 5, 4 ) = 0.39214118169078980444392330174325;

    rungeKuttaFehlberg89Coefficients.aCoefficients( 6, 0 ) = 0.98899281409164665304844765434355e-1;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 6, 3 ) = 0.35138370227963966951204487356703e-2;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 6, 4 ) = 0.12476099983160016621520625872489;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 6, 5 ) = -0.55745546834989799643742901466348e-1;

    rungeKuttaFehlberg89Coefficients.aCoefficients( 7, 0 ) = -0.36806865286242203724153101080691;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 7, 4 ) = -0.22273897469476007645024020944166e+1;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 7, 5 ) = 0.13742908256702910729565691245744e+1;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 7, 6 ) = 0.20497390027111603002159354092206e+1;

    rungeKuttaFehlberg89Coefficients.aCoefficients( 8, 0 ) = 0.45467962641347150077351950603349e-1;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 8, 5 ) = 0.32542131701589147114677469648853;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 8, 6 ) = 0.28476660138527908888182420573687;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 8, 7 ) = 0.97837801675979152435868397271099e-2;

    rungeKuttaFehlberg89Coefficients.aCoefficients( 9, 0 ) = 0.60842071062622057051094145205182e-1;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 9, 5 ) = -0.21184565744037007526325275251206e-1;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 9, 6 ) = 0.19596557266170831957464490662983;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 9, 7 ) = -0.42742640364817603675144835342899e-2;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 9, 8 ) = 0.17434365736814911965323452558189e-1;

    rungeKuttaFehlberg89Coefficients.aCoefficients( 10, 0 ) = 0.54059783296931917365785724111182e-1;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 10, 6 ) = 0.11029825597828926530283127648228;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 10, 7 ) = -0.12565008520072556414147763782250e-2;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 10, 8 ) = 0.36790043477581460136384043566339e-2;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 10, 9 ) = -0.57780542770972073040840628571866e-1;

    rungeKuttaFehlberg89Coefficients.aCoefficients( 11, 0 ) = 0.12732477068667114646645181799160;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 11, 7 ) = 0.11448805006396105323658875721817;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 11, 8 ) = 0.28773020709697992776202201849198;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 11, 9 ) = 0.50945379459611363153735885079465;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 11, 10 ) = -0.14799682244372575900242144449640;

    rungeKuttaFehlberg89Coefficients.aCoefficients( 12, 0 ) = -0.36526793876616740535848544394333e-2;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 12, 5 ) = 0.81629896012318919777819421247030e-1;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 12, 6 ) = -0.38607735635693506490517694343215;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 12, 7 ) = 0.30862242924605106450474166025206e-1;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 12, 8 ) = -0.58077254528320602815829374733518e-1;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 12, 9 ) = 0.33598659328884971493143451362322;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 12, 10 ) = 0.41066880401949958613549622786417;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 12, 11 ) = -0.11840245972355985520633156154536e-1;

    rungeKuttaFehlberg89Coefficients.aCoefficients( 13, 0 ) = -0.12375357921245143254979096135669e+1;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 13, 5 ) = -0.24430768551354785358734861366763e+2;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 13, 6 ) = 0.54779568932778656050436528991173;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 13, 7 ) = -0.44413863533413246374959896569346e+1;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 13, 8 ) = 0.10013104813713266094792617851022e+2;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 13, 9 ) = -0.14995773102051758447170985073142e+2;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 13, 10 ) = 0.58946948523217013620824539651427e+1;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 13, 11 ) = 0.17380377503428984877616857440542e+1;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 13, 12 ) = 0.27512330693166730263758622860276e+2;

    rungeKuttaFehlberg89Coefficients.aCoefficients( 14, 0 ) = -0.35260859388334522700502958875588;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 14, 5 ) = -0.18396103144848270375044198988231;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 14, 6 ) = -0.65570189449741645138006879985251;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 14, 7 ) = -0.39086144880439863435025520241310;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 14, 8 ) = 0.26794646712850022936584423271209;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 14, 9 ) = -0.10383022991382490865769858507427e+1;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 14, 10 ) = 0.16672327324258671664727346168501e+1;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 14, 11 ) = 0.49551925855315977067732967071441;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 14, 12 ) = 0.11394001132397063228586738141784e+1;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 14, 13 ) = 0.51336696424658613688199097191534e-1;

    rungeKuttaFehlberg89Coefficients.aCoefficients( 15, 0 ) = 0.10464847340614810391873002406755e-2;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 15, 8 ) = -0.67163886844990282237778446178020e-2;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 15, 9 ) = 0.81828762189425021265330065248999e-2;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 15 ,10 ) = -0.42640342864483347277142138087561e-2;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 15, 11 ) = 0.28009029474168936545976331153703e-3;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 15, 12 ) = -0.87835333876238676639057813145633e-2;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 15, 13 ) = 0.10254505110825558084217769664009e-1;

    rungeKuttaFehlberg89Coefficients.aCoefficients( 16, 0 ) = -0.13536550786174067080442168889966e+1;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 16, 5 ) = -0.18396103144848270375044198988231;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 16, 6 ) = -0.65570189449741645138006879985251;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 16, 7 ) = -0.39086144880439863435025520241310;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 16, 8 ) = 0.27466285581299925758962207732989;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 16, 9 ) = -0.10464851753571915887035188572676e+1;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 16, 10 ) = 0.16714967667123155012004488306588e+1;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 16, 11 ) = 0.49523916825841808131186990740287;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 16, 12 ) = 0.11481836466273301905225795954930e+1;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 16, 13 ) = 0.41082191313833055603981327527525e-1;
    rungeKuttaFehlberg89Coefficients.aCoefficients( 16, 15 ) = 1.0;
    

    // Define c-coefficients for the Runge-Kutta method of order 8 and a total of 17 stages.
    rungeKuttaFehlberg89Coefficients.cCoefficients = Eigen::VectorXd::Zero( 17 );
    rungeKuttaFehlberg89Coefficients.cCoefficients( 1 ) = 0.44368940376498183109599404281370;
    rungeKuttaFehlberg89Coefficients.cCoefficients( 2 ) = 0.66553410564747274664399106422055;
    rungeKuttaFehlberg89Coefficients.cCoefficients( 3 ) = 0.99830115847120911996598659633083;
    rungeKuttaFehlberg89Coefficients.cCoefficients( 4 ) = 0.31550000000000000000000000000000;
    rungeKuttaFehlberg89Coefficients.cCoefficients( 5 ) = 0.50544100948169068626516126737384;
    rungeKuttaFehlberg89Coefficients.cCoefficients( 6 ) = 0.17142857142857142857142857142857;
    rungeKuttaFehlberg89Coefficients.cCoefficients( 7 ) = 0.82857142857142857142857142857143;
    rungeKuttaFehlberg89Coefficients.cCoefficients( 8 ) = 0.66543966121011562534953769255586;
    rungeKuttaFehlberg89Coefficients.cCoefficients( 9 ) = 0.24878317968062652069722274560771;
    rungeKuttaFehlberg89Coefficients.cCoefficients( 10 ) = 0.10900000000000000000000000000000;
    rungeKuttaFehlberg89Coefficients.cCoefficients( 11 ) = 0.89100000000000000000000000000000;
    rungeKuttaFehlberg89Coefficients.cCoefficients( 12 ) = 0.39950000000000000000000000000000;
    rungeKuttaFehlberg89Coefficients.cCoefficients( 13 ) = 0.60050000000000000000000000000000;
    rungeKuttaFehlberg89Coefficients.cCoefficients( 14 ) = 1.0;
    rungeKuttaFehlberg89Coefficients.cCoefficients( 16 ) = 1.0;

    // Define b-coefficients for the Runge-Kutta method of order 8 and a total of 17 stages.
    rungeKuttaFehlberg89Coefficients.bCoefficients = Eigen::MatrixXd::Zero( 2, 17 );
    rungeKuttaFehlberg89Coefficients.bCoefficients( 0, 0 ) = 0.32256083500216249913612900960247e-1;
    rungeKuttaFehlberg89Coefficients.bCoefficients( 0, 8 ) = 0.25983725283715403018887023171963;
    rungeKuttaFehlberg89Coefficients.bCoefficients( 0, 9 ) = 0.92847805996577027788063714302190e-1;
    rungeKuttaFehlberg89Coefficients.bCoefficients( 0, 10 ) = 0.16452339514764342891647731842800;
    rungeKuttaFehlberg89Coefficients.bCoefficients( 0, 11 ) = 0.17665951637860074367084298397547;
    rungeKuttaFehlberg89Coefficients.bCoefficients( 0, 12 ) = 0.23920102320352759374108933320941;
    rungeKuttaFehlberg89Coefficients.bCoefficients( 0, 13 ) = 0.39484274604202853746752118829325e-2;
    rungeKuttaFehlberg89Coefficients.bCoefficients( 0, 14 ) = 0.30726495475860640406368305522124e-1;
    
    rungeKuttaFehlberg89Coefficients.bCoefficients( 1, 0 ) = 0.062982578976076890319981206482371;
    rungeKuttaFehlberg89Coefficients.bCoefficients( 1, 8 ) = 0.25983725283715403018887023171963;
    rungeKuttaFehlberg89Coefficients.bCoefficients( 1, 9 ) = 0.92847805996577027788063714302190e-1;
    rungeKuttaFehlberg89Coefficients.bCoefficients( 1, 10 ) = 0.16452339514764342891647731842800;
    rungeKuttaFehlberg89Coefficients.bCoefficients( 1, 11 ) = 0.17665951637860074367084298397547;
    rungeKuttaFehlberg89Coefficients.bCoefficients( 1, 12 ) = 0.23920102320352759374108933320941;
    rungeKuttaFehlberg89Coefficients.bCoefficients( 1, 13 ) = 0.39484274604202853746752118829325e-2;
    rungeKuttaFehlberg89Coefficients.bCoefficients( 1, 14 ) = 0.061452990951721280812736611044248;
    rungeKuttaFehlberg89Coefficients.bCoefficients( 1, 15 ) = -0.30726495475860640406368305522124e-1;
    rungeKuttaFehlberg89Coefficients.bCoefficients( 1, 16 ) = -0.30726495475860640406368305522124e-1;

    // Set the name of these coefficients.
    rungeKuttaFehlberg89Coefficients.name = "Runge-Kutta-Fehlberg 8/9";
}

//! Initialize RKV89 coefficients.
void initializeRungeKuttaVerner89Coefficients(
        RungeKuttaCoefficients& rungeKuttaVerner89Coefficients )
{
    // Define characteristics of coefficient set.
    rungeKuttaVerner89Coefficients.lowerOrder = 8;
    rungeKuttaVerner89Coefficients.higherOrder = 9;
    rungeKuttaVerner89Coefficients.orderEstimateToIntegrate
            = RungeKuttaCoefficients::lower;

    // This coefficient set is taken from
    // Cooper, G. J., & Verner, J. H. (1972).
    // Some Explicit Runge-Kutta Methods of High Order.
    // SIAM Journal on Numerical Analysis, 9(3), 389–405.
    // http://www.jstor.org/stable/2156139

    // Compute the square root of 6 to re-use later on.
    double SQRT6 = std::sqrt( 6.0 );

    // a-coefficients for the Runge-Kutta-Verner method of order 8
    // with an embedded 9th-order method for stepsize control and a total of 16 stages.
    rungeKuttaVerner89Coefficients.aCoefficients = Eigen::MatrixXd::Zero( 16, 16 );
    rungeKuttaVerner89Coefficients.aCoefficients( 1, 0 ) = 1.0 / 12.0;

    rungeKuttaVerner89Coefficients.aCoefficients( 2, 0 ) = 1.0 / 27.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 2, 1 ) = 2.0 / 27.0;

    rungeKuttaVerner89Coefficients.aCoefficients( 3, 0 ) = 1.0 / 24.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 3, 2 ) = 1.0 / 8.0;

    rungeKuttaVerner89Coefficients.aCoefficients( 4, 0 ) = (4.0 + 94.0 * SQRT6) / 375.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 4, 2 ) = (-94.0 - 84.0 * SQRT6) / 125.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 4, 3 ) = (328.0 + 208.0 * SQRT6) / 375.0;

    rungeKuttaVerner89Coefficients.aCoefficients( 5, 0 ) = (9.0 - SQRT6) / 150.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 5, 3 ) = (312.0 + 32.0 * SQRT6) / 1425.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 5, 4 ) = (69.0 + 29.0 * SQRT6) / 570.0;

    rungeKuttaVerner89Coefficients.aCoefficients( 6, 0 ) = (927.0 - 347.0 * SQRT6) / 1250.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 6, 3 ) = (-16248.0 + 7328.0 * SQRT6) / 9375.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 6, 4 ) = (-489.0 + 179.0 * SQRT6) / 3750.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 6, 5 ) = (14268.0 - 5798.0 * SQRT6) / 9375.0;

    rungeKuttaVerner89Coefficients.aCoefficients( 7, 0 ) = 2.0 / 27.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 7, 5 ) = (16.0 - SQRT6) / 54.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 7, 6 ) = (16.0 + SQRT6) / 54.0;

    rungeKuttaVerner89Coefficients.aCoefficients( 8, 0 ) = 19.0 / 256.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 8, 5 ) = (118.0 - 23.0 * SQRT6) / 512.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 8, 6 ) = (118.0 + 23.0 * SQRT6) / 512.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 8, 7 ) = -9.0 / 256.0;

    rungeKuttaVerner89Coefficients.aCoefficients( 9, 0 ) = 11.0 / 144.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 9, 5 ) = (266.0 - SQRT6) / 864.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 9, 6 ) = (266.0 + SQRT6) / 864.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 9, 7 ) = -1.0 / 16.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 9, 8 ) = -8.0 / 27.0;

    rungeKuttaVerner89Coefficients.aCoefficients( 10, 0 ) = (5034.0 - 271.0 * SQRT6) / 61440.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 10, 6 ) = (7859.0 - 1626.0 * SQRT6) / 10240.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 10, 7 ) = (-2232.0 + 813.0 * SQRT6) / 20480.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 10, 8 ) = (-594.0 + 271.0 * SQRT6) / 960.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 10, 9 ) = (657.0 - 813.0 * SQRT6) / 5120.0;

    rungeKuttaVerner89Coefficients.aCoefficients( 11, 0 ) = (5996.0 - 3794.0 * SQRT6) / 405.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 11, 5 ) = (-4342.0 - 338.0 * SQRT6) / 9.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 11, 6 ) = (154922.0 - 40458.0 * SQRT6) / 135.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 11, 7 ) = (-4176.0 + 3794.0 * SQRT6) / 45.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 11, 8 ) = (-340864.0 + 242816.0 * SQRT6) / 405.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 11, 9 ) = (26304.0 - 15176.0 * SQRT6) / 45.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 11, 10 ) = -26624.0 / 81.0;

    rungeKuttaVerner89Coefficients.aCoefficients( 12, 0 ) = (3793.0 + 2168.0 * SQRT6) / 103680.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 12, 5 ) = (4042.0 + 2263.0 * SQRT6) / 13824.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 12, 6 ) = (-231278.0 + 40717.0 * SQRT6) / 69120.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 12, 7 ) = (7947.0 - 2168.0 * SQRT6) / 11520.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 12, 8 ) = (1048.0 - 542.0 * SQRT6) / 405.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 12, 9 ) = (-1383.0 + 542.0 * SQRT6) / 720.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 12, 10 ) = 2624.0 / 1053.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 12, 11 ) = 3.0 / 1664.0;

    rungeKuttaVerner89Coefficients.aCoefficients( 13, 0 ) = - 137.0 / 1296.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 13, 5 ) = (5642.0 - 337.0 * SQRT6) / 864.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 13, 6 ) = (5642.0 + 337.0 * SQRT6) / 864.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 13, 7 ) = -299.0 / 48.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 13, 8 ) = 184.0 / 81.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 13, 9 ) = -44.0 / 9.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 13, 10 ) = -5120.0 / 1053.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 13, 11 ) = -11.0 / 468.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 13, 12 ) = 16.0 / 9.0;

    rungeKuttaVerner89Coefficients.aCoefficients( 14, 0 ) = (33617.0 - 2168.0 * SQRT6) / 518400.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 14, 5 ) = (-3846.0 + 31.0 * SQRT6) / 13824.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 14, 6 ) = (155338.0 - 52807.0 * SQRT6) / 345600.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 14, 7 ) = (-12537.0 + 2168.0 * SQRT6) / 57600.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 14, 8 ) = (92.0 + 542.0 * SQRT6) / 2025.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 14, 9 ) = (-1797.0 - 542.0 * SQRT6) / 3600.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 14, 10 ) = 320.0 / 567.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 14, 11 ) = -1.0 / 1920.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 14, 12 ) = 4.0 / 105.0;

    rungeKuttaVerner89Coefficients.aCoefficients( 15, 0 ) = (-36487.0 - 30352.0 * SQRT6) / 279600.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 15, 5 ) = (-29666.0 - 4499.0 * SQRT6) / 7456.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 15, 6 ) = (2779182.0 - 615973.0 * SQRT6) / 186400.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 15, 7 ) = (-94329.0 + 91056.0 * SQRT6) / 93200.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 15, 8 ) = (-232192.0 + 121408.0 * SQRT6) / 17475.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 15, 9 ) = (101226.0 - 22764.0 * SQRT6) / 5825.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 15, 10 ) = -169984.0 / 9087.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 15, 11 ) = -87.0 / 30290.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 15, 12 ) = 492.0 / 1165.0;
    rungeKuttaVerner89Coefficients.aCoefficients( 15, 14 ) = 1260.0 / 233.0;

    // c-coefficients for the Runge-Kutta-Verner method of order 8
    // with an embedded 9th-order method for stepsize control and a total of 16 stages.
    rungeKuttaVerner89Coefficients.cCoefficients = Eigen::VectorXd::Zero( 16 );
    rungeKuttaVerner89Coefficients.cCoefficients( 1 ) = 1.0 / 12.0;
    rungeKuttaVerner89Coefficients.cCoefficients( 2 ) = 1.0 / 9.0;
    rungeKuttaVerner89Coefficients.cCoefficients( 3 ) = 1.0 / 6.0;
    rungeKuttaVerner89Coefficients.cCoefficients( 4 ) = (2.0 + 2.0 * SQRT6) / 15.0;
    rungeKuttaVerner89Coefficients.cCoefficients( 5 ) = (6.0 + 1.0 * SQRT6) / 15.0;
    rungeKuttaVerner89Coefficients.cCoefficients( 6 ) = (6.0 - 1.0 * SQRT6) / 15.0;
    rungeKuttaVerner89Coefficients.cCoefficients( 7 ) = 2.0 / 3.0;
    rungeKuttaVerner89Coefficients.cCoefficients( 8 ) = 1.0 / 2.0;
    rungeKuttaVerner89Coefficients.cCoefficients( 9 ) = 1.0 / 3.0;
    rungeKuttaVerner89Coefficients.cCoefficients( 10 ) = 1.0 / 4.0;
    rungeKuttaVerner89Coefficients.cCoefficients( 11 ) = 4.0 / 3.0;
    rungeKuttaVerner89Coefficients.cCoefficients( 12 ) = 5.0 / 6.0;
    rungeKuttaVerner89Coefficients.cCoefficients( 13 ) = 1.0;
    rungeKuttaVerner89Coefficients.cCoefficients( 14 ) = 1.0 / 6.0;
    rungeKuttaVerner89Coefficients.cCoefficients( 15 ) = 1.0;

    // b-coefficients for the Runge-Kutta-Verner method of order 8
    // with an embedded 9th-order method for stepsize control and a total of 16 stages.
    rungeKuttaVerner89Coefficients.bCoefficients = Eigen::MatrixXd::Zero( 2, 16 );
    rungeKuttaVerner89Coefficients.bCoefficients( 0, 0 ) = 103.0 / 1680.0;
    rungeKuttaVerner89Coefficients.bCoefficients( 0, 7 ) = -27.0 / 140.0;
    rungeKuttaVerner89Coefficients.bCoefficients( 0, 8 ) = 76.0 / 105.0;
    rungeKuttaVerner89Coefficients.bCoefficients( 0, 9 ) = -201.0 / 280.0;
    rungeKuttaVerner89Coefficients.bCoefficients( 0, 10 ) = 1024.0 / 1365.0;
    rungeKuttaVerner89Coefficients.bCoefficients( 0, 11 ) = 3.0 / 7280.0;
    rungeKuttaVerner89Coefficients.bCoefficients( 0, 12 ) = 12.0 / 35.0;
    rungeKuttaVerner89Coefficients.bCoefficients( 0, 13 ) = 9.0 / 280.0;

    rungeKuttaVerner89Coefficients.bCoefficients( 1, 0 ) = 23.0 / 525.0;
    rungeKuttaVerner89Coefficients.bCoefficients( 1, 7 ) = 171.0 / 1400.0;
    rungeKuttaVerner89Coefficients.bCoefficients( 1, 8 ) = 86.0 / 525.0;
    rungeKuttaVerner89Coefficients.bCoefficients( 1, 9 ) = 93.0 / 280.0;
    rungeKuttaVerner89Coefficients.bCoefficients( 1, 10 ) = -2048.0 / 6825.0;
    rungeKuttaVerner89Coefficients.bCoefficients( 1, 11 ) = -3.0 / 18200.0;
    rungeKuttaVerner89Coefficients.bCoefficients( 1, 12 ) = 39.0 / 175.0;
    rungeKuttaVerner89Coefficients.bCoefficients( 1, 14 ) = 9.0 / 25.0;
    rungeKuttaVerner89Coefficients.bCoefficients( 1, 15 ) = 233.0 / 4200.0;

    // Set the name of these coefficients.
    rungeKuttaVerner89Coefficients.name = "Runge-Kutta-Verner 8/9";
}

//! Initialize Runge Kutta Feagin 10(8) coefficients.
void initializeRungeKuttaFeagin108Coefficients(
        RungeKuttaCoefficients& rungeKutta108Coefficients )
{
    // Define characteristics of coefficient set.
    rungeKutta108Coefficients.lowerOrder = 8;
    rungeKutta108Coefficients.higherOrder = 10;
    rungeKutta108Coefficients.orderEstimateToIntegrate
            = RungeKuttaCoefficients::higher;

    // This coefficient set is taken from
    // Feagin, T. (2012). High-order explicit Runge-Kutta methods using m-symmetry.
    // Neural, Parallel & Scientific Computations, 20(3-4), 437-458.

    // Define a-coefficients for the Runge-Kutta method.
    rungeKutta108Coefficients.aCoefficients = Eigen::MatrixXd::Zero( 17, 16 );
    rungeKutta108Coefficients.aCoefficients( 1, 0 ) = 0.100000000000000000000000000000000000000000000000000000000000;

    rungeKutta108Coefficients.aCoefficients( 2, 0 ) = -0.915176561375291440520015019275342154318951387664369720564660;
    rungeKutta108Coefficients.aCoefficients( 2, 1 ) = 1.45453440217827322805250021715664459117622483736537873607016;

    rungeKutta108Coefficients.aCoefficients( 3, 0 ) = 0.202259190301118170324681949205488413821477543637878380814562;
    rungeKutta108Coefficients.aCoefficients( 3, 1 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 3, 2 ) = 0.606777570903354510974045847616465241464432630913635142443687;

    rungeKutta108Coefficients.aCoefficients( 4, 0 ) = 0.184024714708643575149100693471120664216774047979591417844635;
    rungeKutta108Coefficients.aCoefficients( 4, 1 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 4, 2 ) = 0.197966831227192369068141770510388793370637287463360401555746;
    rungeKutta108Coefficients.aCoefficients( 4, 3 ) = -0.0729547847313632629185146671595558023015011608914382961421311;

    rungeKutta108Coefficients.aCoefficients( 5, 0 ) = 0.0879007340206681337319777094132125475918886824944548534041378;
    rungeKutta108Coefficients.aCoefficients( 5, 1 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 5, 2 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 5, 3 ) = 0.410459702520260645318174895920453426088035325902848695210406;
    rungeKutta108Coefficients.aCoefficients( 5, 4 ) = 0.482713753678866489204726942976896106809132737721421333413261;

    rungeKutta108Coefficients.aCoefficients( 6, 0 ) = 0.0859700504902460302188480225945808401411132615636600222593880;
    rungeKutta108Coefficients.aCoefficients( 6, 1 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 6, 2 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 6, 3 ) = 0.330885963040722183948884057658753173648240154838402033448632;
    rungeKutta108Coefficients.aCoefficients( 6, 4 ) = 0.489662957309450192844507011135898201178015478433790097210790;
    rungeKutta108Coefficients.aCoefficients( 6, 5 ) = -0.0731856375070850736789057580558988816340355615025188195854775;

    rungeKutta108Coefficients.aCoefficients( 7, 0 ) = 0.120930449125333720660378854927668953958938996999703678812621;
    rungeKutta108Coefficients.aCoefficients( 7, 1 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 7, 2 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 7, 3 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 7, 4 ) = 0.260124675758295622809007617838335174368108756484693361887839;
    rungeKutta108Coefficients.aCoefficients( 7, 5 ) = 0.0325402621549091330158899334391231259332716675992700000776101;
    rungeKutta108Coefficients.aCoefficients( 7, 6 ) = -0.0595780211817361001560122202563305121444953672762930724538856;

    rungeKutta108Coefficients.aCoefficients( 8, 0 ) = 0.110854379580391483508936171010218441909425780168656559807038;
    rungeKutta108Coefficients.aCoefficients( 8, 1 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 8, 2 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 8, 3 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 8, 4 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 8, 5 ) = -0.0605761488255005587620924953655516875526344415354339234619466;
    rungeKutta108Coefficients.aCoefficients( 8, 6 ) = 0.321763705601778390100898799049878904081404368603077129251110;
    rungeKutta108Coefficients.aCoefficients( 8, 7 ) = 0.510485725608063031577759012285123416744672137031752354067590;

    rungeKutta108Coefficients.aCoefficients( 9, 0 ) = 0.112054414752879004829715002761802363003717611158172229329393;
    rungeKutta108Coefficients.aCoefficients( 9, 1 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 9, 2 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 9, 3 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 9, 4 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 9, 5 ) = -0.144942775902865915672349828340980777181668499748506838876185;
    rungeKutta108Coefficients.aCoefficients( 9, 6 ) = -0.333269719096256706589705211415746871709467423992115497968724;
    rungeKutta108Coefficients.aCoefficients( 9, 7 ) = 0.499269229556880061353316843969978567860276816592673201240332;
    rungeKutta108Coefficients.aCoefficients( 9, 8 ) = 0.509504608929686104236098690045386253986643232352989602185060;

    rungeKutta108Coefficients.aCoefficients( 10, 0 ) = 0.113976783964185986138004186736901163890724752541486831640341;
    rungeKutta108Coefficients.aCoefficients( 10, 1 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 10, 2 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 10, 3 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 10, 4 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 10, 5 ) = -0.0768813364203356938586214289120895270821349023390922987406384;
    rungeKutta108Coefficients.aCoefficients( 10, 6 ) = 0.239527360324390649107711455271882373019741311201004119339563;
    rungeKutta108Coefficients.aCoefficients( 10, 7 ) = 0.397774662368094639047830462488952104564716416343454639902613;
    rungeKutta108Coefficients.aCoefficients( 10, 8 ) = 0.0107558956873607455550609147441477450257136782823280838547024;
    rungeKutta108Coefficients.aCoefficients( 10, 9 ) = -0.327769124164018874147061087350233395378262992392394071906457;

    rungeKutta108Coefficients.aCoefficients( 11, 0 ) = 0.0798314528280196046351426864486400322758737630423413945356284;
    rungeKutta108Coefficients.aCoefficients( 11, 1 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 11, 2 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 11, 3 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 11, 4 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 11, 5 ) = -0.0520329686800603076514949887612959068721311443881683526937298;
    rungeKutta108Coefficients.aCoefficients( 11, 6 ) = -0.0576954146168548881732784355283433509066159287152968723021864;
    rungeKutta108Coefficients.aCoefficients( 11, 7 ) = 0.194781915712104164976306262147382871156142921354409364738090;
    rungeKutta108Coefficients.aCoefficients( 11, 8 ) = 0.145384923188325069727524825977071194859203467568236523866582;
    rungeKutta108Coefficients.aCoefficients( 11, 9 ) = -0.0782942710351670777553986729725692447252077047239160551335016;
    rungeKutta108Coefficients.aCoefficients( 11, 10 ) = -0.114503299361098912184303164290554670970133218405658122674674;

    rungeKutta108Coefficients.aCoefficients( 12, 0 ) = 0.985115610164857280120041500306517278413646677314195559520529;
    rungeKutta108Coefficients.aCoefficients( 12, 1 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 12, 2 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 12, 3 ) = 0.330885963040722183948884057658753173648240154838402033448632;
    rungeKutta108Coefficients.aCoefficients( 12, 4 ) = 0.489662957309450192844507011135898201178015478433790097210790;
    rungeKutta108Coefficients.aCoefficients( 12, 5 ) = -1.37896486574843567582112720930751902353904327148559471526397;
    rungeKutta108Coefficients.aCoefficients( 12, 6 ) = -0.861164195027635666673916999665534573351026060987427093314412;
    rungeKutta108Coefficients.aCoefficients( 12, 7 ) = 5.78428813637537220022999785486578436006872789689499172601856;
    rungeKutta108Coefficients.aCoefficients( 12, 8 ) = 3.28807761985103566890460615937314805477268252903342356581925;
    rungeKutta108Coefficients.aCoefficients( 12, 9 ) = -2.38633905093136384013422325215527866148401465975954104585807;
    rungeKutta108Coefficients.aCoefficients( 12, 10 ) = -3.25479342483643918654589367587788726747711504674780680269911;
    rungeKutta108Coefficients.aCoefficients( 12, 11 ) = -2.16343541686422982353954211300054820889678036420109999154887;

    rungeKutta108Coefficients.aCoefficients( 13, 0 ) = 0.895080295771632891049613132336585138148156279241561345991710;
    rungeKutta108Coefficients.aCoefficients( 13, 1 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 13, 2 ) = 0.197966831227192369068141770510388793370637287463360401555746;
    rungeKutta108Coefficients.aCoefficients( 13, 3 ) = -0.0729547847313632629185146671595558023015011608914382961421311;
    rungeKutta108Coefficients.aCoefficients( 13, 4 ) = 0.0000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 13, 5 ) = -0.851236239662007619739049371445966793289359722875702227166105;
    rungeKutta108Coefficients.aCoefficients( 13, 6 ) = 0.398320112318533301719718614174373643336480918103773904231856;
    rungeKutta108Coefficients.aCoefficients( 13, 7 ) = 3.63937263181035606029412920047090044132027387893977804176229;
    rungeKutta108Coefficients.aCoefficients( 13, 8 ) = 1.54822877039830322365301663075174564919981736348973496313065;
    rungeKutta108Coefficients.aCoefficients( 13, 9 ) = -2.12221714704053716026062427460427261025318461146260124401561;
    rungeKutta108Coefficients.aCoefficients( 13, 10 ) = -1.58350398545326172713384349625753212757269188934434237975291;
    rungeKutta108Coefficients.aCoefficients( 13, 11 ) = -1.71561608285936264922031819751349098912615880827551992973034;
    rungeKutta108Coefficients.aCoefficients( 13, 12 ) = -0.0244036405750127452135415444412216875465593598370910566069132;

    rungeKutta108Coefficients.aCoefficients( 14, 0 ) = -0.915176561375291440520015019275342154318951387664369720564660;
    rungeKutta108Coefficients.aCoefficients( 14, 1 ) = 1.45453440217827322805250021715664459117622483736537873607016;
    rungeKutta108Coefficients.aCoefficients( 14, 2 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 14, 3 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 14, 4 ) = -0.777333643644968233538931228575302137803351053629547286334469;
    rungeKutta108Coefficients.aCoefficients( 14, 5 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 14, 6 ) = -0.0910895662155176069593203555807484200111889091770101799647985;
    rungeKutta108Coefficients.aCoefficients( 14, 7 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 14, 8 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 14, 9 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 14, 10 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 14, 11 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 14, 12 ) = 0.0910895662155176069593203555807484200111889091770101799647985;
    rungeKutta108Coefficients.aCoefficients( 14, 13 ) = 0.777333643644968233538931228575302137803351053629547286334469;

    rungeKutta108Coefficients.aCoefficients( 15, 0 ) = 0.100000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 15, 1 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 15, 2 ) = -0.157178665799771163367058998273128921867183754126709419409654;
    rungeKutta108Coefficients.aCoefficients( 15, 3 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 15, 4 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 15, 5 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 15, 6 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 15, 7 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 15, 8 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 15, 9 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 15, 10 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 15, 11 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 15, 12 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 15, 13 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 15, 14 ) = 0.157178665799771163367058998273128921867183754126709419409654;

    rungeKutta108Coefficients.aCoefficients( 16, 0 ) = 0.181781300700095283888472062582262379650443831463199521664945;
    rungeKutta108Coefficients.aCoefficients( 16, 1 ) = 0.675000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 16, 2 ) = 0.342758159847189839942220553413850871742338734703958919937260;
    rungeKutta108Coefficients.aCoefficients( 16, 3 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 16, 4 ) = 0.259111214548322744512977076191767379267783684543182428778156;
    rungeKutta108Coefficients.aCoefficients( 16, 5 ) = -0.358278966717952089048961276721979397739750634673268802484271;
    rungeKutta108Coefficients.aCoefficients( 16, 6 ) = -1.04594895940883306095050068756409905131588123172378489286080;
    rungeKutta108Coefficients.aCoefficients( 16, 7 ) = 0.930327845415626983292300564432428777137601651182965794680397;
    rungeKutta108Coefficients.aCoefficients( 16, 8 ) = 1.77950959431708102446142106794824453926275743243327790536000;
    rungeKutta108Coefficients.aCoefficients( 16, 9 ) = 0.100000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.aCoefficients( 16, 10 ) = -0.282547569539044081612477785222287276408489375976211189952877;
    rungeKutta108Coefficients.aCoefficients( 16, 11 ) = -0.159327350119972549169261984373485859278031542127551931461821;
    rungeKutta108Coefficients.aCoefficients( 16, 12 ) = -0.145515894647001510860991961081084111308650130578626404945571;
    rungeKutta108Coefficients.aCoefficients( 16, 13 ) = -0.259111214548322744512977076191767379267783684543182428778156;
    rungeKutta108Coefficients.aCoefficients( 16, 14 ) = -0.342758159847189839942220553413850871742338734703958919937260;
    rungeKutta108Coefficients.aCoefficients( 16, 15 ) = -0.675000000000000000000000000000000000000000000000000000000000;
    
    // Define c-coefficients for the Runge-Kutta method.
    rungeKutta108Coefficients.cCoefficients = Eigen::VectorXd::Zero( 17 );
    rungeKutta108Coefficients.cCoefficients( 0 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.cCoefficients( 1 ) = 0.100000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.cCoefficients( 2 ) = 0.539357840802981787532485197881302436857273449701009015505500;
    rungeKutta108Coefficients.cCoefficients( 3 ) = 0.809036761204472681298727796821953655285910174551513523258250;
    rungeKutta108Coefficients.cCoefficients( 4 ) = 0.309036761204472681298727796821953655285910174551513523258250;
    rungeKutta108Coefficients.cCoefficients( 5 ) = 0.981074190219795268254879548310562080489056746118724882027805;
    rungeKutta108Coefficients.cCoefficients( 6 ) = 0.833333333333333333333333333333333333333333333333333333333333;
    rungeKutta108Coefficients.cCoefficients( 7 ) = 0.354017365856802376329264185948796742115824053807373968324184;
    rungeKutta108Coefficients.cCoefficients( 8 ) = 0.882527661964732346425501486979669075182867844268052119663791;
    rungeKutta108Coefficients.cCoefficients( 9 ) = 0.642615758240322548157075497020439535959501736363212695909875;
    rungeKutta108Coefficients.cCoefficients( 10 ) = 0.357384241759677451842924502979560464040498263636787304090125;
    rungeKutta108Coefficients.cCoefficients( 11 ) = 0.117472338035267653574498513020330924817132155731947880336209;
    rungeKutta108Coefficients.cCoefficients( 12 ) = 0.833333333333333333333333333333333333333333333333333333333333;
    rungeKutta108Coefficients.cCoefficients( 13 ) = 0.309036761204472681298727796821953655285910174551513523258250;
    rungeKutta108Coefficients.cCoefficients( 14 ) = 0.539357840802981787532485197881302436857273449701009015505500;
    rungeKutta108Coefficients.cCoefficients( 15 ) = 0.100000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.cCoefficients( 16 ) = 1.00000000000000000000000000000000000000000000000000000000000;

    // Define b-coefficients for the Runge-Kutta method.
    rungeKutta108Coefficients.bCoefficients = Eigen::MatrixXd::Zero( 2, 17 );
    rungeKutta108Coefficients.bCoefficients( 0, 0 ) = 0.0333333333333333333333333333333333333333333333333333333333333;
    rungeKutta108Coefficients.bCoefficients( 0, 1 ) = 1.0 / 36.0;
    rungeKutta108Coefficients.bCoefficients( 0, 2 ) = 0.0333333333333333333333333333333333333333333333333333333333333;
    rungeKutta108Coefficients.bCoefficients( 0, 3 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.bCoefficients( 0, 4 ) = 0.0500000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.bCoefficients( 0, 5 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.bCoefficients( 0, 6 ) = 0.0400000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.bCoefficients( 0, 7 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.bCoefficients( 0, 8 ) = 0.189237478148923490158306404106012326238162346948625830327194;
    rungeKutta108Coefficients.bCoefficients( 0, 9 ) = 0.277429188517743176508360262560654340428504319718040836339472;
    rungeKutta108Coefficients.bCoefficients( 0, 10 ) = 0.277429188517743176508360262560654340428504319718040836339472;
    rungeKutta108Coefficients.bCoefficients( 0, 11 ) = 0.189237478148923490158306404106012326238162346948625830327194;
    rungeKutta108Coefficients.bCoefficients( 0, 12 ) = -0.0400000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.bCoefficients( 0, 13 ) = -0.0500000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.bCoefficients( 0, 14 ) = -0.0333333333333333333333333333333333333333333333333333333333333;
    rungeKutta108Coefficients.bCoefficients( 0, 15 ) = -1.0 / 36.0;
    rungeKutta108Coefficients.bCoefficients( 0, 16 ) = 0.0333333333333333333333333333333333333333333333333333333333333;

    rungeKutta108Coefficients.bCoefficients( 1, 0 ) = 0.0333333333333333333333333333333333333333333333333333333333333;
    rungeKutta108Coefficients.bCoefficients( 1, 1 ) = 0.0250000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.bCoefficients( 1, 2 ) = 0.0333333333333333333333333333333333333333333333333333333333333;
    rungeKutta108Coefficients.bCoefficients( 1, 3 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.bCoefficients( 1, 4 ) = 0.0500000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.bCoefficients( 1, 5 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.bCoefficients( 1, 6 ) = 0.0400000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.bCoefficients( 1, 7 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.bCoefficients( 1, 8 ) = 0.189237478148923490158306404106012326238162346948625830327194;
    rungeKutta108Coefficients.bCoefficients( 1, 9 ) = 0.277429188517743176508360262560654340428504319718040836339472;
    rungeKutta108Coefficients.bCoefficients( 1, 10 ) = 0.277429188517743176508360262560654340428504319718040836339472;
    rungeKutta108Coefficients.bCoefficients( 1, 11 ) = 0.189237478148923490158306404106012326238162346948625830327194;
    rungeKutta108Coefficients.bCoefficients( 1, 12 ) = -0.0400000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.bCoefficients( 1, 13 ) = -0.0500000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.bCoefficients( 1, 14 ) = -0.0333333333333333333333333333333333333333333333333333333333333;
    rungeKutta108Coefficients.bCoefficients( 1, 15 ) = -0.0250000000000000000000000000000000000000000000000000000000000;
    rungeKutta108Coefficients.bCoefficients( 1, 16 ) = 0.0333333333333333333333333333333333333333333333333333333333333;
    
    // Define the name of the coefficient set.
    rungeKutta108Coefficients.name = "Runge-Kutta-Feagin 10/8";
}

//! Initialize Runge Kutta Feagin 12(10) coefficients.
void initializeRungeKuttaFeagin1210Coefficients(
        RungeKuttaCoefficients& rungeKutta1210Coefficients )
{
    // Define characteristics of coefficient set.
    rungeKutta1210Coefficients.lowerOrder = 10;
    rungeKutta1210Coefficients.higherOrder = 12;
    rungeKutta1210Coefficients.orderEstimateToIntegrate
            = RungeKuttaCoefficients::higher;

    // This coefficient set is taken from
    // Feagin, T. (2012). High-order explicit Runge-Kutta methods using m-symmetry.
    // Neural, Parallel & Scientific Computations, 20(3-4), 437-458.

    // Define a-coefficients for the Runge-Kutta method.
    rungeKutta1210Coefficients.aCoefficients = Eigen::MatrixXd::Zero( 25, 24 );
    rungeKutta1210Coefficients.aCoefficients( 1, 0 ) = 0.200000000000000000000000000000000000000000000000000000000000;

    rungeKutta1210Coefficients.aCoefficients( 2, 0 ) = -0.216049382716049382716049382716049382716049382716049382716049;
    rungeKutta1210Coefficients.aCoefficients( 2, 1 ) = 0.771604938271604938271604938271604938271604938271604938271605;

    rungeKutta1210Coefficients.aCoefficients( 3, 0 ) = 0.208333333333333333333333333333333333333333333333333333333333;
    rungeKutta1210Coefficients.aCoefficients( 3, 2 ) = 0.625000000000000000000000000000000000000000000000000000000000;

    rungeKutta1210Coefficients.aCoefficients( 4, 0 ) = 0.193333333333333333333333333333333333333333333333333333333333;
    rungeKutta1210Coefficients.aCoefficients( 4, 2 ) = 0.220000000000000000000000000000000000000000000000000000000000;
    rungeKutta1210Coefficients.aCoefficients( 4, 3 ) = -0.0800000000000000000000000000000000000000000000000000000000000;

    rungeKutta1210Coefficients.aCoefficients( 5, 0 ) = 0.100000000000000000000000000000000000000000000000000000000000;
    rungeKutta1210Coefficients.aCoefficients( 5, 3 ) = 0.400000000000000000000000000000000000000000000000000000000000;
    rungeKutta1210Coefficients.aCoefficients( 5, 4 ) = 0.500000000000000000000000000000000000000000000000000000000000;

    rungeKutta1210Coefficients.aCoefficients( 6, 0 ) = 0.103364471650010477570395435690481791543342708330349879244197;
    rungeKutta1210Coefficients.aCoefficients( 6, 3 ) = 0.124053094528946761061581889237115328211074784955180298044074;
    rungeKutta1210Coefficients.aCoefficients( 6, 4 ) = 0.483171167561032899288836480451962508724109257517289177302380;
    rungeKutta1210Coefficients.aCoefficients( 6, 5 ) = -0.0387530245694763252085681443767620580395733302341368038804290;

    rungeKutta1210Coefficients.aCoefficients( 7, 0 ) = 0.124038261431833324081904585980175168140024670698633612292480;
    rungeKutta1210Coefficients.aCoefficients( 7, 4 ) = 0.217050632197958486317846256953159942875916353757734167684657;
    rungeKutta1210Coefficients.aCoefficients( 7, 5 ) = 0.0137455792075966759812907801835048190594443990939408530842918;
    rungeKutta1210Coefficients.aCoefficients( 7, 6 ) = -0.0661095317267682844455831341498149531672668252085016565917546;

    rungeKutta1210Coefficients.aCoefficients( 8, 0 ) = 0.0914774894856882983144991846980432197088832099976660100090486;
    rungeKutta1210Coefficients.aCoefficients( 8, 5 ) = -0.00544348523717469689965754944144838611346156873847009178068318;
    rungeKutta1210Coefficients.aCoefficients( 8, 6 ) = 0.0680716801688453518578515120895103863112751730758794372203952;
    rungeKutta1210Coefficients.aCoefficients( 8, 7 ) = 0.408394315582641046727306852653894780093303185664924644551239;

    rungeKutta1210Coefficients.aCoefficients( 9, 0 ) = 0.0890013652502551018954509355423841780143232697403434118692699;
    rungeKutta1210Coefficients.aCoefficients( 9, 5 ) = 0.00499528226645532360197793408420692800405891149406814091955810;
    rungeKutta1210Coefficients.aCoefficients( 9, 6 ) = 0.397918238819828997341739603001347156083435060931424970826304;
    rungeKutta1210Coefficients.aCoefficients( 9, 7 ) = 0.427930210752576611068192608300897981558240730580396406312359;
    rungeKutta1210Coefficients.aCoefficients( 9, 8 ) = -0.0865117637557827005740277475955029103267246394128995965941585;

    rungeKutta1210Coefficients.aCoefficients( 10, 0 ) = 0.0695087624134907543112693906409809822706021061685544615255758;
    rungeKutta1210Coefficients.aCoefficients( 10, 5 ) = 0.129146941900176461970759579482746551122871751501482634045487;
    rungeKutta1210Coefficients.aCoefficients( 10, 6 ) = 1.53073638102311295076342566143214939031177504112433874313011;
    rungeKutta1210Coefficients.aCoefficients( 10, 7 ) = 0.577874761129140052546751349454576715334892100418571882718036;
    rungeKutta1210Coefficients.aCoefficients( 10, 8 ) = -0.951294772321088980532340837388859453930924498799228648050949;
    rungeKutta1210Coefficients.aCoefficients( 10, 9 ) = -0.408276642965631951497484981519757463459627174520978426909934;

    rungeKutta1210Coefficients.aCoefficients( 11, 0 ) = 0.0444861403295135866269453507092463581620165501018684152933313;
    rungeKutta1210Coefficients.aCoefficients( 11, 5 ) = -0.00380476867056961731984232686574547203016331563626856065717964;
    rungeKutta1210Coefficients.aCoefficients( 11, 6 ) = 0.0106955064029624200721262602809059154469206077644957399593972;
    rungeKutta1210Coefficients.aCoefficients( 11, 7 ) = 0.0209616244499904333296674205928919920806734650660039898074652;
    rungeKutta1210Coefficients.aCoefficients( 11, 8 ) = -0.0233146023259321786648561431551978077665337818756053603898847;
    rungeKutta1210Coefficients.aCoefficients( 11, 9 ) = 0.00263265981064536974369934736325334761174975280887405725010964;
    rungeKutta1210Coefficients.aCoefficients( 11, 10 ) = 0.00315472768977025060103545855572111407955208306374459723959783;

    rungeKutta1210Coefficients.aCoefficients( 12, 0 ) = 0.0194588815119755475588801096525317761242073762016273186231215;
    rungeKutta1210Coefficients.aCoefficients( 12, 8 ) = 0.0000678512949171812509306121653452367476194364781259165332321534;
    rungeKutta1210Coefficients.aCoefficients( 12, 9 ) = -0.0000429795859049273623271005330230162343568863387724883603675550;
    rungeKutta1210Coefficients.aCoefficients( 12, 10 ) = 0.0000176358982260285155407485928953302139937553442829975734148981;
    rungeKutta1210Coefficients.aCoefficients( 12, 11 ) = 0.0653866627415027051009595231385181033549511358787382098351924;

    rungeKutta1210Coefficients.aCoefficients( 13, 0 ) = 0.206836835664277105916828174798272361078909196043446411598231;
    rungeKutta1210Coefficients.aCoefficients( 13, 8 ) = 0.0166796067104156472828045866664696450306326505094792505215514;
    rungeKutta1210Coefficients.aCoefficients( 13, 9 ) = -0.00879501563200710214457024178249986591130234990219959208704979;
    rungeKutta1210Coefficients.aCoefficients( 13, 10 ) = 0.00346675455362463910824462315246379209427513654098596403637231;
    rungeKutta1210Coefficients.aCoefficients( 13, 11 ) = -0.861264460105717678161432562258351242030270498966891201799225;
    rungeKutta1210Coefficients.aCoefficients( 13, 12 ) = 0.908651882074050281096239478469262145034957129939256789178785;

    rungeKutta1210Coefficients.aCoefficients( 14, 0 ) = 0.0203926084654484010091511314676925686038504449562413004562382;
    rungeKutta1210Coefficients.aCoefficients( 14, 8 ) = 0.0869469392016685948675400555583947505833954460930940959577347;
    rungeKutta1210Coefficients.aCoefficients( 14, 9 ) = -0.0191649630410149842286436611791405053287170076602337673587681;
    rungeKutta1210Coefficients.aCoefficients( 14, 10 ) = 0.00655629159493663287364871573244244516034828755253746024098838;
    rungeKutta1210Coefficients.aCoefficients( 14, 11 ) = 0.0987476128127434780903798528674033899738924968006632201445462;
    rungeKutta1210Coefficients.aCoefficients( 14, 12 ) = 0.00535364695524996055083260173615567408717110247274021056118319;
    rungeKutta1210Coefficients.aCoefficients( 14, 13 ) = 0.301167864010967916837091303817051676920059229784957479998077;

    rungeKutta1210Coefficients.aCoefficients( 15, 0 ) = 0.228410433917778099547115412893004398779136994596948545722283;
    rungeKutta1210Coefficients.aCoefficients( 15, 8 ) = -0.498707400793025250635016567442511512138603770959682292383042;
    rungeKutta1210Coefficients.aCoefficients( 15, 9 ) = 0.134841168335724478552596703792570104791700727205981058201689;
    rungeKutta1210Coefficients.aCoefficients( 15, 10 ) = -0.0387458244055834158439904226924029230935161059142806805674360;
    rungeKutta1210Coefficients.aCoefficients( 15, 11 ) = -1.27473257473474844240388430824908952380979292713250350199641;
    rungeKutta1210Coefficients.aCoefficients( 15, 12 ) = 1.43916364462877165201184452437038081875299303577911839630524;
    rungeKutta1210Coefficients.aCoefficients( 15, 13 ) = -0.214007467967990254219503540827349569639028092344812795499026;
    rungeKutta1210Coefficients.aCoefficients( 15, 14 ) = 0.958202417754430239892724139109781371059908874605153648768037;

    rungeKutta1210Coefficients.aCoefficients( 16, 0 ) = 2.00222477655974203614249646012506747121440306225711721209798;
    rungeKutta1210Coefficients.aCoefficients( 16, 8 ) = 2.06701809961524912091954656438138595825411859673341600679555;
    rungeKutta1210Coefficients.aCoefficients( 16, 9 ) = 0.623978136086139541957471279831494466155292316167021080663140;
    rungeKutta1210Coefficients.aCoefficients( 16, 10 ) = -0.0462283685500311430283203554129062069391947101880112723185773;
    rungeKutta1210Coefficients.aCoefficients( 16, 11 ) = -8.84973288362649614860075246727118949286604835457092701094630;
    rungeKutta1210Coefficients.aCoefficients( 16, 12 ) = 7.74257707850855976227437225791835589560188590785037197433615;
    rungeKutta1210Coefficients.aCoefficients( 16, 13 ) = -0.588358519250869210993353314127711745644125882130941202896436;
    rungeKutta1210Coefficients.aCoefficients( 16, 14 ) = -1.10683733362380649395704708016953056176195769617014899442903;
    rungeKutta1210Coefficients.aCoefficients( 16, 15 ) = -0.929529037579203999778397238291233214220788057511899747507074;

    rungeKutta1210Coefficients.aCoefficients( 17, 0 ) = 3.13789533412073442934451608989888796808161259330322100268310;
    rungeKutta1210Coefficients.aCoefficients( 17, 5 ) = 0.129146941900176461970759579482746551122871751501482634045487;
    rungeKutta1210Coefficients.aCoefficients( 17, 6 ) = 1.53073638102311295076342566143214939031177504112433874313011;
    rungeKutta1210Coefficients.aCoefficients( 17, 7 ) = 0.577874761129140052546751349454576715334892100418571882718036;
    rungeKutta1210Coefficients.aCoefficients( 17, 8 ) = 5.42088263055126683050056840891857421941300558851862156403363;
    rungeKutta1210Coefficients.aCoefficients( 17, 9 ) = 0.231546926034829304872663800877643660904880180835945693836936;
    rungeKutta1210Coefficients.aCoefficients( 17, 10 ) = 0.0759292995578913560162301311785251873561801342333194895292058;
    rungeKutta1210Coefficients.aCoefficients( 17, 11 ) = -12.3729973380186513287414553402595806591349822617535905976253;
    rungeKutta1210Coefficients.aCoefficients( 17, 12 ) = 9.85455883464769543935957209317369202080367765721777101906955;
    rungeKutta1210Coefficients.aCoefficients( 17, 13 ) = 0.0859111431370436529579357709052367772889980495122329601159540;
    rungeKutta1210Coefficients.aCoefficients( 17, 14 ) = -5.65242752862643921117182090081762761180392602644189218673969;
    rungeKutta1210Coefficients.aCoefficients( 17, 15 ) = -1.94300935242819610883833776782364287728724899124166920477873;
    rungeKutta1210Coefficients.aCoefficients( 17, 16 ) = -0.128352601849404542018428714319344620742146491335612353559923;

    rungeKutta1210Coefficients.aCoefficients( 18, 0 ) = 1.38360054432196014878538118298167716825163268489922519995564;
    rungeKutta1210Coefficients.aCoefficients( 18, 5 ) = 0.00499528226645532360197793408420692800405891149406814091955810;
    rungeKutta1210Coefficients.aCoefficients( 18, 6 ) = 0.397918238819828997341739603001347156083435060931424970826304;
    rungeKutta1210Coefficients.aCoefficients( 18, 7 ) = 0.427930210752576611068192608300897981558240730580396406312359;
    rungeKutta1210Coefficients.aCoefficients( 18, 8 ) = -1.30299107424475770916551439123047573342071475998399645982146;
    rungeKutta1210Coefficients.aCoefficients( 18, 9 ) = 0.661292278669377029097112528107513072734573412294008071500699;
    rungeKutta1210Coefficients.aCoefficients( 18, 10 ) = -0.144559774306954349765969393688703463900585822441545655530145;
    rungeKutta1210Coefficients.aCoefficients( 18, 11 ) = -6.96576034731798203467853867461083919356792248105919255460819;
    rungeKutta1210Coefficients.aCoefficients( 18, 12 ) = 6.65808543235991748353408295542210450632193197576935120716437;
    rungeKutta1210Coefficients.aCoefficients( 18, 13 ) = -1.66997375108841486404695805725510845049807969199236227575796;
    rungeKutta1210Coefficients.aCoefficients( 18, 14 ) = 2.06413702318035263832289040301832647130604651223986452170089;
    rungeKutta1210Coefficients.aCoefficients( 18, 15 ) = -0.674743962644306471862958129570837723192079875998405058648892;
    rungeKutta1210Coefficients.aCoefficients( 18, 16 ) = -0.00115618834794939500490703608435907610059605754935305582045729;
    rungeKutta1210Coefficients.aCoefficients( 18, 17 ) = -0.00544057908677007389319819914241631024660726585015012485938593;

    rungeKutta1210Coefficients.aCoefficients( 19, 0 ) = 0.951236297048287669474637975894973552166903378983475425758226;
    rungeKutta1210Coefficients.aCoefficients( 19, 4 ) = 0.217050632197958486317846256953159942875916353757734167684657;
    rungeKutta1210Coefficients.aCoefficients( 19, 5 ) = 0.0137455792075966759812907801835048190594443990939408530842918;
    rungeKutta1210Coefficients.aCoefficients( 19, 6 ) = -0.0661095317267682844455831341498149531672668252085016565917546;
    rungeKutta1210Coefficients.aCoefficients( 19, 8 ) = 0.152281696736414447136604697040747131921486432699422112099617;
    rungeKutta1210Coefficients.aCoefficients( 19, 9 ) = -0.337741018357599840802300793133998004354643424457539667670080;
    rungeKutta1210Coefficients.aCoefficients( 19, 10 ) = -0.0192825981633995781534949199286824400469353110630787982121133;
    rungeKutta1210Coefficients.aCoefficients( 19, 11 ) = -3.68259269696866809932409015535499603576312120746888880201882;
    rungeKutta1210Coefficients.aCoefficients( 19, 12 ) = 3.16197870406982063541533528419683854018352080342887002331312;
    rungeKutta1210Coefficients.aCoefficients( 19, 13 ) = -0.370462522106885290716991856022051125477943482284080569177386;
    rungeKutta1210Coefficients.aCoefficients( 19, 14 ) = -0.0514974200365440434996434456698127984941168616474316871020314;
    rungeKutta1210Coefficients.aCoefficients( 19, 15 ) = -0.000829625532120152946787043541792848416659382675202720677536554;
    rungeKutta1210Coefficients.aCoefficients( 19, 16 ) = 0.00000279801041419278598986586589070027583961355402640879503213503;
    rungeKutta1210Coefficients.aCoefficients( 19, 17 ) = 0.0418603916412360287969841020776788461794119440689356178942252;
    rungeKutta1210Coefficients.aCoefficients( 19, 18 ) = 0.279084255090877355915660874555379649966282167560126269290222;

    rungeKutta1210Coefficients.aCoefficients( 20, 0 ) = 0.103364471650010477570395435690481791543342708330349879244197;
    rungeKutta1210Coefficients.aCoefficients( 20, 3 ) = 0.124053094528946761061581889237115328211074784955180298044074;
    rungeKutta1210Coefficients.aCoefficients( 20, 4 ) = 0.483171167561032899288836480451962508724109257517289177302380;
    rungeKutta1210Coefficients.aCoefficients( 20, 5 ) = -0.0387530245694763252085681443767620580395733302341368038804290;
    rungeKutta1210Coefficients.aCoefficients( 20, 7 ) = -0.438313820361122420391059788940960176420682836652600698580091;
    rungeKutta1210Coefficients.aCoefficients( 20, 9 ) = -0.218636633721676647685111485017151199362509373698288330593486;
    rungeKutta1210Coefficients.aCoefficients( 20, 10 ) = -0.0312334764394719229981634995206440349766174759626578122323015;
    rungeKutta1210Coefficients.aCoefficients( 20, 17 ) = 0.0312334764394719229981634995206440349766174759626578122323015;
    rungeKutta1210Coefficients.aCoefficients( 20, 18 ) = 0.218636633721676647685111485017151199362509373698288330593486;
    rungeKutta1210Coefficients.aCoefficients( 20, 19 ) = 0.438313820361122420391059788940960176420682836652600698580091;

    rungeKutta1210Coefficients.aCoefficients( 21, 0 ) = 0.193333333333333333333333333333333333333333333333333333333333;
    rungeKutta1210Coefficients.aCoefficients( 21, 2 ) = 0.220000000000000000000000000000000000000000000000000000000000;
    rungeKutta1210Coefficients.aCoefficients( 21, 3 ) = -0.0800000000000000000000000000000000000000000000000000000000000;
    rungeKutta1210Coefficients.aCoefficients( 21, 6 ) = 0.0984256130499315928152900286856048243348202521491288575952143;
    rungeKutta1210Coefficients.aCoefficients( 21, 7 ) = -0.196410889223054653446526504390100417677539095340135532418849;
    rungeKutta1210Coefficients.aCoefficients( 21, 9 ) = 0.436457930493068729391826122587949137609670676712525034763317;
    rungeKutta1210Coefficients.aCoefficients( 21, 10 ) = 0.0652613721675721098560370939805555698350543810708414716730270;
    rungeKutta1210Coefficients.aCoefficients( 21, 17 ) = -0.0652613721675721098560370939805555698350543810708414716730270;
    rungeKutta1210Coefficients.aCoefficients( 21, 18 ) = -0.436457930493068729391826122587949137609670676712525034763317;
    rungeKutta1210Coefficients.aCoefficients( 21, 19 ) = 0.196410889223054653446526504390100417677539095340135532418849;
    rungeKutta1210Coefficients.aCoefficients( 21, 20 ) = -0.0984256130499315928152900286856048243348202521491288575952143;

    rungeKutta1210Coefficients.aCoefficients( 22, 0 ) = -0.216049382716049382716049382716049382716049382716049382716049;
    rungeKutta1210Coefficients.aCoefficients( 22, 1 ) = 0.771604938271604938271604938271604938271604938271604938271605;
    rungeKutta1210Coefficients.aCoefficients( 22, 4 ) = -0.666666666666666666666666666666666666666666666666666666666667;
    rungeKutta1210Coefficients.aCoefficients( 22, 6 ) = -0.390696469295978451446999802258495981249099665294395945559163;
    rungeKutta1210Coefficients.aCoefficients( 22, 20 ) = 0.390696469295978451446999802258495981249099665294395945559163;
    rungeKutta1210Coefficients.aCoefficients( 22, 21 ) = 0.666666666666666666666666666666666666666666666666666666666667;

    rungeKutta1210Coefficients.aCoefficients( 23, 0 ) = 0.200000000000000000000000000000000000000000000000000000000000;
    rungeKutta1210Coefficients.aCoefficients( 23, 2 ) = -0.164609053497942386831275720164609053497942386831275720164609;
    rungeKutta1210Coefficients.aCoefficients( 23, 22 ) = 0.164609053497942386831275720164609053497942386831275720164609;

    rungeKutta1210Coefficients.aCoefficients( 24, 0 ) = 1.47178724881110408452949550989023611293535315518571691939396;
    rungeKutta1210Coefficients.aCoefficients( 24, 1 ) = 0.787500000000000000000000000000000000000000000000000000000000;
    rungeKutta1210Coefficients.aCoefficients( 24, 2 ) = 0.421296296296296296296296296296296296296296296296296296296296;
    rungeKutta1210Coefficients.aCoefficients( 24, 4 ) = 0.291666666666666666666666666666666666666666666666666666666667;
    rungeKutta1210Coefficients.aCoefficients( 24, 6 ) = 0.348600717628329563206854421629657569274689947367847465753757;
    rungeKutta1210Coefficients.aCoefficients( 24, 7 ) = 0.229499544768994849582890233710555447073823569666506700662510;
    rungeKutta1210Coefficients.aCoefficients( 24, 8 ) = 5.79046485790481979159831978177003471098279506036722411333192;
    rungeKutta1210Coefficients.aCoefficients( 24, 9 ) = 0.418587511856506868874073759426596207226461447604248151080016;
    rungeKutta1210Coefficients.aCoefficients( 24, 10 ) = 0.307039880222474002649653817490106690389251482313213999386651;
    rungeKutta1210Coefficients.aCoefficients( 24, 11 ) = -4.68700905350603332214256344683853248065574415794742040470287;
    rungeKutta1210Coefficients.aCoefficients( 24, 12 ) = 3.13571665593802262152038152399873856554395436199962915429076;
    rungeKutta1210Coefficients.aCoefficients( 24, 13 ) = 1.40134829710965720817510506275620441055845017313930508348898;
    rungeKutta1210Coefficients.aCoefficients( 24, 14 ) = -5.52931101439499023629010306005764336421276055777658156400910;
    rungeKutta1210Coefficients.aCoefficients( 24, 15 ) = -0.853138235508063349309546894974784906188927508039552519557498;
    rungeKutta1210Coefficients.aCoefficients( 24, 16 ) = 0.103575780373610140411804607167772795518293914458500175573749;
    rungeKutta1210Coefficients.aCoefficients( 24, 17 ) = -0.140474416950600941142546901202132534870665923700034957196546;
    rungeKutta1210Coefficients.aCoefficients( 24, 18 ) = -0.418587511856506868874073759426596207226461447604248151080016;
    rungeKutta1210Coefficients.aCoefficients( 24, 19 ) = -0.229499544768994849582890233710555447073823569666506700662510;
    rungeKutta1210Coefficients.aCoefficients( 24, 20 ) = -0.348600717628329563206854421629657569274689947367847465753757;
    rungeKutta1210Coefficients.aCoefficients( 24, 21 ) = -0.291666666666666666666666666666666666666666666666666666666667;
    rungeKutta1210Coefficients.aCoefficients( 24, 22 ) = -0.421296296296296296296296296296296296296296296296296296296296;
    rungeKutta1210Coefficients.aCoefficients( 24, 23 ) = -0.787500000000000000000000000000000000000000000000000000000000;
    
    // Define c-coefficients for the Runge-Kutta method.
    rungeKutta1210Coefficients.cCoefficients = Eigen::VectorXd::Zero( 25 );
    rungeKutta1210Coefficients.cCoefficients( 1 ) = 0.200000000000000000000000000000000000000000000000000000000000;
    rungeKutta1210Coefficients.cCoefficients( 2 ) = 0.555555555555555555555555555555555555555555555555555555555556;
    rungeKutta1210Coefficients.cCoefficients( 3 ) = 0.833333333333333333333333333333333333333333333333333333333333;
    rungeKutta1210Coefficients.cCoefficients( 4 ) = 0.333333333333333333333333333333333333333333333333333333333333;
    rungeKutta1210Coefficients.cCoefficients( 5 ) = 1.00000000000000000000000000000000000000000000000000000000000;
    rungeKutta1210Coefficients.cCoefficients( 6 ) = 0.671835709170513812712245661002797570438953420568682550710222;
    rungeKutta1210Coefficients.cCoefficients( 7 ) = 0.288724941110620201935458488967024976908118598341806976469674;
    rungeKutta1210Coefficients.cCoefficients( 8 ) = 0.562500000000000000000000000000000000000000000000000000000000;
    rungeKutta1210Coefficients.cCoefficients( 9 ) = 0.833333333333333333333333333333333333333333333333333333333333;
    rungeKutta1210Coefficients.cCoefficients( 10 ) = 0.947695431179199287562380162101836721649589325892740646458322;
    rungeKutta1210Coefficients.cCoefficients( 11 ) = 0.0548112876863802643887753674810754475842153612931128785028369;
    rungeKutta1210Coefficients.cCoefficients( 12 ) = 0.0848880518607165350639838930162674302064148175640019542045934;
    rungeKutta1210Coefficients.cCoefficients( 13 ) = 0.265575603264642893098114059045616835297201264164077621448665;
    rungeKutta1210Coefficients.cCoefficients( 14 ) = 0.500000000000000000000000000000000000000000000000000000000000;
    rungeKutta1210Coefficients.cCoefficients( 15 ) = 0.734424396735357106901885940954383164702798735835922378551335;
    rungeKutta1210Coefficients.cCoefficients( 16 ) = 0.915111948139283464936016106983732569793585182435998045795407;
    rungeKutta1210Coefficients.cCoefficients( 17 ) = 0.947695431179199287562380162101836721649589325892740646458322;
    rungeKutta1210Coefficients.cCoefficients( 18 ) = 0.833333333333333333333333333333333333333333333333333333333333;
    rungeKutta1210Coefficients.cCoefficients( 19 ) = 0.288724941110620201935458488967024976908118598341806976469674;
    rungeKutta1210Coefficients.cCoefficients( 20 ) = 0.671835709170513812712245661002797570438953420568682550710222;
    rungeKutta1210Coefficients.cCoefficients( 21 ) = 0.333333333333333333333333333333333333333333333333333333333333;
    rungeKutta1210Coefficients.cCoefficients( 22 ) = 0.555555555555555555555555555555555555555555555555555555555556;
    rungeKutta1210Coefficients.cCoefficients( 23 ) = 0.200000000000000000000000000000000000000000000000000000000000;
    rungeKutta1210Coefficients.cCoefficients( 24 ) = 1.00000000000000000000000000000000000000000000000000000000000;

    // Define b-coefficients for the Runge-Kutta method.
    rungeKutta1210Coefficients.bCoefficients = Eigen::MatrixXd::Zero( 2, 25 );
    rungeKutta1210Coefficients.bCoefficients( 0, 0 ) = 0.0238095238095238095238095238095238095238095238095238095238095;
    rungeKutta1210Coefficients.bCoefficients( 0, 1 ) = 1.0 / 10.0;
    rungeKutta1210Coefficients.bCoefficients( 0, 2 ) = 0.0312500000000000000000000000000000000000000000000000000000000;
    rungeKutta1210Coefficients.bCoefficients( 0, 4 ) = 0.0416666666666666666666666666666666666666666666666666666666667;
    rungeKutta1210Coefficients.bCoefficients( 0, 6 ) = 0.0500000000000000000000000000000000000000000000000000000000000;
    rungeKutta1210Coefficients.bCoefficients( 0, 7 ) = 0.0500000000000000000000000000000000000000000000000000000000000;
    rungeKutta1210Coefficients.bCoefficients( 0, 9 ) = 0.100000000000000000000000000000000000000000000000000000000000;
    rungeKutta1210Coefficients.bCoefficients( 0, 10 ) = 0.0714285714285714285714285714285714285714285714285714285714286;
    rungeKutta1210Coefficients.bCoefficients( 0, 12 ) = 0.138413023680782974005350203145033146748813640089941234591267;
    rungeKutta1210Coefficients.bCoefficients( 0, 13 ) = 0.215872690604931311708935511140681138965472074195773051123019;
    rungeKutta1210Coefficients.bCoefficients( 0, 14 ) = 0.243809523809523809523809523809523809523809523809523809523810;
    rungeKutta1210Coefficients.bCoefficients( 0, 15 ) = 0.215872690604931311708935511140681138965472074195773051123019;
    rungeKutta1210Coefficients.bCoefficients( 0, 16 ) = 0.138413023680782974005350203145033146748813640089941234591267;
    rungeKutta1210Coefficients.bCoefficients( 0, 17 ) = -0.0714285714285714285714285714285714285714285714285714285714286;
    rungeKutta1210Coefficients.bCoefficients( 0, 18 ) = -0.100000000000000000000000000000000000000000000000000000000000;
    rungeKutta1210Coefficients.bCoefficients( 0, 19 ) = -0.0500000000000000000000000000000000000000000000000000000000000;
    rungeKutta1210Coefficients.bCoefficients( 0, 20 ) = -0.0500000000000000000000000000000000000000000000000000000000000;
    rungeKutta1210Coefficients.bCoefficients( 0, 21 ) = -0.0416666666666666666666666666666666666666666666666666666666667;
    rungeKutta1210Coefficients.bCoefficients( 0, 22 ) = -0.0312500000000000000000000000000000000000000000000000000000000;
    rungeKutta1210Coefficients.bCoefficients( 0, 23 ) = -1.0 / 10.0;
    rungeKutta1210Coefficients.bCoefficients( 0, 24 ) = 0.0238095238095238095238095238095238095238095238095238095238095;
    
    rungeKutta1210Coefficients.bCoefficients( 1, 0 ) = 0.0238095238095238095238095238095238095238095238095238095238095;
    rungeKutta1210Coefficients.bCoefficients( 1, 1 ) = 0.0234375000000000000000000000000000000000000000000000000000000;
    rungeKutta1210Coefficients.bCoefficients( 1, 2 ) = 0.0312500000000000000000000000000000000000000000000000000000000;
    rungeKutta1210Coefficients.bCoefficients( 1, 4 ) = 0.0416666666666666666666666666666666666666666666666666666666667;
    rungeKutta1210Coefficients.bCoefficients( 1, 6 ) = 0.0500000000000000000000000000000000000000000000000000000000000;
    rungeKutta1210Coefficients.bCoefficients( 1, 7 ) = 0.0500000000000000000000000000000000000000000000000000000000000;
    rungeKutta1210Coefficients.bCoefficients( 1, 9 ) = 0.100000000000000000000000000000000000000000000000000000000000;
    rungeKutta1210Coefficients.bCoefficients( 1, 10 ) = 0.0714285714285714285714285714285714285714285714285714285714286;
    rungeKutta1210Coefficients.bCoefficients( 1, 12 ) = 0.138413023680782974005350203145033146748813640089941234591267;
    rungeKutta1210Coefficients.bCoefficients( 1, 13 ) = 0.215872690604931311708935511140681138965472074195773051123019;
    rungeKutta1210Coefficients.bCoefficients( 1, 14 ) = 0.243809523809523809523809523809523809523809523809523809523810;
    rungeKutta1210Coefficients.bCoefficients( 1, 15 ) = 0.215872690604931311708935511140681138965472074195773051123019;
    rungeKutta1210Coefficients.bCoefficients( 1, 16 ) = 0.138413023680782974005350203145033146748813640089941234591267;
    rungeKutta1210Coefficients.bCoefficients( 1, 17 ) = -0.0714285714285714285714285714285714285714285714285714285714286;
    rungeKutta1210Coefficients.bCoefficients( 1, 18 ) = -0.100000000000000000000000000000000000000000000000000000000000;
    rungeKutta1210Coefficients.bCoefficients( 1, 19 ) = -0.0500000000000000000000000000000000000000000000000000000000000;
    rungeKutta1210Coefficients.bCoefficients( 1, 20 ) = -0.0500000000000000000000000000000000000000000000000000000000000;
    rungeKutta1210Coefficients.bCoefficients( 1, 21 ) = -0.0416666666666666666666666666666666666666666666666666666666667;
    rungeKutta1210Coefficients.bCoefficients( 1, 22 ) = -0.0312500000000000000000000000000000000000000000000000000000000;
    rungeKutta1210Coefficients.bCoefficients( 1, 23 ) = -0.0234375000000000000000000000000000000000000000000000000000000;
    rungeKutta1210Coefficients.bCoefficients( 1, 24 ) = 0.0238095238095238095238095238095238095238095238095238095238095;
    
    // Define the name of the coefficient set.
    rungeKutta1210Coefficients.name = "Runge-Kutta-Feagin 12/10";
}

//! Initialize Runge Kutta Feagin 14(12) coefficients.
void initializeRungeKuttaFeagin1412Coefficients(
        RungeKuttaCoefficients& rungeKutta1412Coefficients )
{
    // Define characteristics of coefficient set.
    rungeKutta1412Coefficients.lowerOrder = 12;
    rungeKutta1412Coefficients.higherOrder = 14;
    rungeKutta1412Coefficients.orderEstimateToIntegrate
            = RungeKuttaCoefficients::higher;

    // This coefficient set is taken from
    // Feagin, T. (2012). High-order explicit Runge-Kutta methods using m-symmetry.
    // Neural, Parallel & Scientific Computations, 20(3-4), 437-458.

    // Define a-coefficients for the Runge-Kutta Feagin method.
    rungeKutta1412Coefficients.aCoefficients = Eigen::MatrixXd::Zero( 35, 34 );
    rungeKutta1412Coefficients.aCoefficients( 1, 0 ) = 0.1111111111111111111111111111111111111111111111111111111111111111111111111111111111111;
    
    rungeKutta1412Coefficients.aCoefficients( 2, 0 ) = -0.8333333333333333333333333333333333333333333333333333333333333333333333333333333333333;
    rungeKutta1412Coefficients.aCoefficients( 2, 1 ) = 1.388888888888888888888888888888888888888888888888888888888888888888888888888888888889;
    
    rungeKutta1412Coefficients.aCoefficients( 3, 0 ) = 0.2083333333333333333333333333333333333333333333333333333333333333333333333333333333333;
    rungeKutta1412Coefficients.aCoefficients( 3, 2 ) = 0.625;
    
    rungeKutta1412Coefficients.aCoefficients( 4, 0 ) = 0.1933333333333333333333333333333333333333333333333333333333333333333333333333333333333;
    rungeKutta1412Coefficients.aCoefficients( 4, 2 ) = 0.22;
    rungeKutta1412Coefficients.aCoefficients( 4, 3 ) = -0.8e-1;
    
    rungeKutta1412Coefficients.aCoefficients( 5, 0 ) = 0.1;
    rungeKutta1412Coefficients.aCoefficients( 5, 3 ) = 0.4;
    rungeKutta1412Coefficients.aCoefficients( 5, 4 ) = 0.5;
    
    rungeKutta1412Coefficients.aCoefficients( 6, 0 ) = 0.1034845616366797766729935465119103444997447982019713166066629728281981965079290745983;
    rungeKutta1412Coefficients.aCoefficients( 6, 3 ) = 0.1220688873064072225896440828689620771395927148341621347412746563709055937325311521675;
    rungeKutta1412Coefficients.aCoefficients( 6, 4 ) = 0.4825744903312466224751347801256881128659190238501680496794015023696413273862321544150;
    rungeKutta1412Coefficients.aCoefficients( 6, 5 ) = -0.3814096000156069997308862400056202056641130724784114774219699240039767479629669855696e-1;
    
    rungeKutta1412Coefficients.aCoefficients( 7, 0 ) = 0.1243805266540944128815164208687993162684914663596714231632892354628068537117612942798;
    rungeKutta1412Coefficients.aCoefficients( 7, 4 ) = 0.2261202821975843014222386629792029011967523207426331439651447460281196206643404356021;
    rungeKutta1412Coefficients.aCoefficients( 7, 5 ) = 0.1378858876180808806076958370164778145309694174914933853635428709475288586061552782365e-1;
    rungeKutta1412Coefficients.aCoefficients( 7, 6 ) = -0.6722101339966844497493995074143058569500863415253821828561997825320849038679063596730e-1;
    
    rungeKutta1412Coefficients.aCoefficients( 8, 0 ) = 0.9369190656596738155308854560830059338663496952177500856556033862893464429241815101000e-1;
    rungeKutta1412Coefficients.aCoefficients( 8, 5 ) = -0.6134068434505109872294989956416647356209145071288588710070986068372475355320835997035e-2;
    rungeKutta1412Coefficients.aCoefficients( 8, 6 ) = 0.2160198256255030637088600976598665734909794332781173201886676706066128640340557614360;
    rungeKutta1412Coefficients.aCoefficients( 8, 7 ) = 0.4236950635157619373376190739609767532058674695441235326831157041055522397561196508237;
    
    rungeKutta1412Coefficients.aCoefficients( 9, 0 ) = 0.8384798124090526646169687913728140859805331392249111310693346670107922625197375034871e-1;
    rungeKutta1412Coefficients.aCoefficients( 9, 5 ) = -0.1179493671009738143197550560312957753679619605907361507776128268875265788248790903515e-1;
    rungeKutta1412Coefficients.aCoefficients( 9, 6 ) = -0.2472990205688126523394738387431945983259928403533401326974984247503501083158412965835;
    rungeKutta1412Coefficients.aCoefficients( 9, 7 ) = 0.9780808583677290122593130140812916655037406554767339407565991037499621093437371932341e-1;
    rungeKutta1412Coefficients.aCoefficients( 9, 8 ) = 0.2175906892434206313600086517678603183441681200247821768799893467069296630467914197921;
    
    rungeKutta1412Coefficients.aCoefficients( 10, 0 ) = 0.6152553597694282279545623896143147143334239690648211074539397569215087099333654844097e-1;
    rungeKutta1412Coefficients.aCoefficients( 10, 5 ) = 0.5922327803245033080429900057980465247383895604442571368349896773084347972825775455007e-2;
    rungeKutta1412Coefficients.aCoefficients( 10, 6 ) = 0.4703261599638411122172243032058941134553625307461088250108483236601604516650193568134;
    rungeKutta1412Coefficients.aCoefficients( 10, 7 ) = 0.2996888638486790008539818370961923991368311216717812791841936858888827504094204242461;
    rungeKutta1412Coefficients.aCoefficients( 10, 8 ) = -0.2476568775939949146899922763298108258539580692639470955481886317480090967647905771626;
    rungeKutta1412Coefficients.aCoefficients( 10, 9 ) = 0.1108950297714376828939998518390617145224451736006787182086245987785252503880550245038;
    
    rungeKutta1412Coefficients.aCoefficients( 11, 0 ) = 0.4197000733627825798617928647872777872134836565431046112459945389674655429048057710370e-1;
    rungeKutta1412Coefficients.aCoefficients( 11, 5 ) = -0.317987696266205093901912847692712407988609169703103952205634e-2;
    rungeKutta1412Coefficients.aCoefficients( 11, 6 ) = 0.8063977149061920772608217115203795063935431115674197501197468839656405367779525213500;
    rungeKutta1412Coefficients.aCoefficients( 11, 7 ) = 0.9759831264123889790935228506842888513146720480030545503571875185550549213299958241991e-1;
    rungeKutta1412Coefficients.aCoefficients( 11, 8 ) = 0.7785755781583989090275124464529272389997634605941819649588520345133050850477185489203;
    rungeKutta1412Coefficients.aCoefficients( 11, 9 ) = 0.2048904238315994281894992020981056033120292350814206535748293420400885242747823516625;
    rungeKutta1412Coefficients.aCoefficients( 11, 10 ) = -1.562615796274681883070709439505278252114628922364243608928053762634922556160297217820;
    
    rungeKutta1412Coefficients.aCoefficients( 12, 0 ) = 0.4377267822337301635744652424953398116882149670716141232569729223172939742940416733395e-1;
    rungeKutta1412Coefficients.aCoefficients( 12, 8 ) = 0.6243650275201952087943586285809336252816312169030959172012504609444028241438248581173e-2;
    rungeKutta1412Coefficients.aCoefficients( 12, 9 ) = 0.2000430971095773149944351654696478568290662322182649696087680691197048872391143823078;
    rungeKutta1412Coefficients.aCoefficients( 12, 10 ) = -0.8053283678049830368238571620489029119233928873370293148442058928084075077460302544840e-2;
    rungeKutta1412Coefficients.aCoefficients( 12, 11 ) = 0.2115175280673965219157119035233996013168778251575505730512208770404786743066139905871e-1;
    
    rungeKutta1412Coefficients.aCoefficients( 13, 0 ) = 0.2834992503635145630950235919207173122471376548964770977684956012393009143065795513785e-1;
    rungeKutta1412Coefficients.aCoefficients( 13, 8 ) = 0.2491632048558174075389491488059951494598846535854176800982219995075912885766744587193e-2;
    rungeKutta1412Coefficients.aCoefficients( 13, 9 ) = 0.2301387878545931496383998463737427687720871226381422342236583655735620108657836993957e-1;
    rungeKutta1412Coefficients.aCoefficients( 13, 10 ) = -0.3221559566929770987244760924671208781894636047606204610433085107190031098987004938258e-2;
    rungeKutta1412Coefficients.aCoefficients( 13, 11 ) = 0.9884425494476646689463354144878852560408199827860146481292993078049373245839618405001e-2;
    rungeKutta1412Coefficients.aCoefficients( 13, 12 ) = -0.2130107713288873513843076428759273848866345654295724666320922464722154754985568313136e-1;
    
    rungeKutta1412Coefficients.aCoefficients( 14, 0 ) = 0.3435118942902430010494322347351479430833531749807014262686507474123120416010457867571;
    rungeKutta1412Coefficients.aCoefficients( 14, 8 ) = 0.2104519120236273856090970119990106557888074052256267000419050051487632641518018732685;
    rungeKutta1412Coefficients.aCoefficients( 14, 9 ) = 1.034274520572304119364829268288257099386679996983247401666929134177931632176349026735;
    rungeKutta1412Coefficients.aCoefficients( 14, 10 ) = 0.6003036458644224870512404482066405749390780924061569454673075686417142117164254262878e-2;
    rungeKutta1412Coefficients.aCoefficients( 14, 11 ) = 0.8559381250996195375780121060024077289150626526164160058172684354881277648341960563008;
    rungeKutta1412Coefficients.aCoefficients( 14, 12 ) = -0.9772350050367668108722648523725256330131076568928396776974412446349105799705851506077;
    rungeKutta1412Coefficients.aCoefficients( 14, 13 ) = -0.6600269804792946946162250138563276937205739812199748747775581736879654453322759683463;
    
    rungeKutta1412Coefficients.aCoefficients( 15, 0 ) = -0.1435740016721680695382063999350763666577559543783998809757153672896315044426183882232e-1;
    rungeKutta1412Coefficients.aCoefficients( 15, 8 ) = -0.3662532700490399702936857968489747917331190817335522078657913621382824038988807796287e-1;
    rungeKutta1412Coefficients.aCoefficients( 15, 9 ) = 0.3502549756362136819768494069798465243467890824711035742020654749717518291597210559354e-1;
    rungeKutta1412Coefficients.aCoefficients( 15, 10 ) = 0.3609460163621135089317866587583352398236899298642376718895880083960486970547825683491e-1;
    rungeKutta1412Coefficients.aCoefficients( 15, 11 ) = -0.2652199675536811063515959468346019236496270124574642848667252606942739787160130682669e-1;
    rungeKutta1412Coefficients.aCoefficients( 15, 12 ) = 0.4456990113056981196389115375088399081043363230822267716707629092315111479614958673268e-1;
    rungeKutta1412Coefficients.aCoefficients( 15, 13 ) = 0.1243430933313582432862255957417864480389734088951067419167759990001419776217292554191;
    rungeKutta1412Coefficients.aCoefficients( 15, 14 ) = 0.4138296932394806944035124962043359604261929086744760344472227418812310333088685698274e-2;
    
    rungeKutta1412Coefficients.aCoefficients( 16, 0 ) = 0.3560324044251202909756091163980891762641062223797488026536968101501275805051289760823;
    rungeKutta1412Coefficients.aCoefficients( 16, 8 ) = -0.450192758947562595966821779075956175110645100214763601190349;
    rungeKutta1412Coefficients.aCoefficients( 16, 9 ) = 0.430527907083710898626656292808782917793030154094709462877146;
    rungeKutta1412Coefficients.aCoefficients( 16, 10 ) = 0.5119730290110222376685569603940716920771257870306513863906805244405042813755411512463;
    rungeKutta1412Coefficients.aCoefficients( 16, 11 ) = 0.9083036388864042603901591246381102139974962148199046305445452541870528153933236088278;
    rungeKutta1412Coefficients.aCoefficients( 16, 12 ) = -1.239210933719339317573724691515340288544138892486057261860887966510000755220957594942;
    rungeKutta1412Coefficients.aCoefficients( 16, 13 ) = -0.6490486616717614651416723488790625539054028319671910976544025456235491510559878435372;
    rungeKutta1412Coefficients.aCoefficients( 16, 14 ) = 0.2517089045868192922104805299489705414048878529314474912189256354259853776829630937658;
    rungeKutta1412Coefficients.aCoefficients( 16, 15 ) = 0.7799064703455863988107567952823344760235405934115501870206452879298798513199886085571;
    
    rungeKutta1412Coefficients.aCoefficients( 17, 0 ) = 0.1309356874065130664068812064188349801274704382131924878449566575565302965696195341197e-1;
    rungeKutta1412Coefficients.aCoefficients( 17, 12 ) = -0.9320530679851139459084619627671082378586315096846671421247697017556505173897578610165e-4;
    rungeKutta1412Coefficients.aCoefficients( 17, 13 ) = 0.5053743342622993596400904431385907267709423447161223817027456630856526555478831396014e-1;
    rungeKutta1412Coefficients.aCoefficients( 17, 14 ) = 0.8044703419444879791095791096101977976413118689308653610493721999399129417586629251430e-6;
    rungeKutta1412Coefficients.aCoefficients( 17, 15 ) = 0.5917260294941711905287557427777172598443409719243215281782302034071342229921661278343e-3;
    rungeKutta1412Coefficients.aCoefficients( 17, 16 ) = -0.4016147221545573370646916849063755877322642479500938046774565993013424294867398455789e-6;
    
    rungeKutta1412Coefficients.aCoefficients( 18, 0 ) = 0.2079264844660530125419445440007656521672552061443734079797586969853055549175505457737e-1;
    rungeKutta1412Coefficients.aCoefficients( 18, 12 ) = 0.5826959188000859151019026978372841089514061030298715701031065480360641416298102920851e-3;
    rungeKutta1412Coefficients.aCoefficients( 18, 13 ) = -0.8017007323588159390833421865258527466405584659196335246554992680506588169863285718822e-2;
    rungeKutta1412Coefficients.aCoefficients( 18, 14 ) = 0.4038476438471369403751708217435605704841172903308955066191655368223862388605213690921e-5;
    rungeKutta1412Coefficients.aCoefficients( 18, 15 ) = 0.8546099980555061442250561145675356025101146220336224918025961310211940592009621595606e-1;
    rungeKutta1412Coefficients.aCoefficients( 18, 16 ) = -0.2044864809358042427067075696910043079044428375526774562331430989116458814609927891477e-5;
    rungeKutta1412Coefficients.aCoefficients( 18, 17 ) = 0.1053285788244318933997994029790939973542409042351728431465827473723673651882417656762;
    
    rungeKutta1412Coefficients.aCoefficients( 19, 0 ) = 1.401534497957360214154462473557713067184864529175977331289881318884096354294079099114;
    rungeKutta1412Coefficients.aCoefficients( 19, 12 ) = -0.2302520009842212616162724103674156212611302982744556219175010157057031125814669239016;
    rungeKutta1412Coefficients.aCoefficients( 19, 13 ) = -7.211068404669129056595822371068742471658564935099615697324849532576890894506619405031;
    rungeKutta1412Coefficients.aCoefficients( 19, 14 ) = 0.3729015606948363352369953278521323402177595666786623882373057096229137360164435411243e-2;
    rungeKutta1412Coefficients.aCoefficients( 19, 15 ) = -4.714154957271250206787781793922247570113233732218200980194845522013711035054762664884;
    rungeKutta1412Coefficients.aCoefficients( 19, 16 ) = -0.1763676575453492420538419950327976735749038866956001340593194717236122233799126229446e-2;
    rungeKutta1412Coefficients.aCoefficients( 19, 17 ) = 7.641305480386987655630293108802376511851733678139370059818519661401442202665741111270;
    rungeKutta1412Coefficients.aCoefficients( 19, 18 ) = 3.506020436597518349898960829497447109682129498933753736341591881470708008233521976557;
    
    rungeKutta1412Coefficients.aCoefficients( 20, 0 ) = 11.95146506941206867993723858307164016744736108265535168242754934626543968357331742096;
    rungeKutta1412Coefficients.aCoefficients( 20, 12 ) = 7.794809321081759687835167002317643882202842795989809549197917776161588225206322580459;
    rungeKutta1412Coefficients.aCoefficients( 20, 13 ) = -56.45013938673257925235609911209042814404681000613405538635967763011214022629172907669;
    rungeKutta1412Coefficients.aCoefficients( 20, 14 ) = 0.9123763069306449013445304492902766457096074504036737047499704936582270274950128398912e-1;
    rungeKutta1412Coefficients.aCoefficients( 20, 15 ) = -12.73362799254348862019455243091992750381627175299189605168457824373779389828110581300;
    rungeKutta1412Coefficients.aCoefficients( 20, 16 ) = -0.3968959219047197123135428109397366747123830704331478729319411886202118671113516172493e-1;
    rungeKutta1412Coefficients.aCoefficients( 20, 17 ) = 54.43921418835708869962257651553077918614383784233053341001985423053366890118247056463;
    rungeKutta1412Coefficients.aCoefficients( 20, 18 ) = -3.644116379215692368464069903613506458067214784092667356589342345057374050114156075061;
    rungeKutta1412Coefficients.aCoefficients( 20, 19 ) = -0.8045032499105099108990307879585794993156949132107878807481027183961246894903442258757;
    
    rungeKutta1412Coefficients.aCoefficients( 21, 0 ) = -148.8094265071004884278388682686476255619306120821485965777899951377767737092911763254;
    rungeKutta1412Coefficients.aCoefficients( 21, 12 ) = -91.72952782912564843579356624023216234952287290363542836291360346578688265538801398361;
    rungeKutta1412Coefficients.aCoefficients( 21, 13 ) = 707.6561449715983598345757192863357161548211289666495623584804744987957677893379157809;
    rungeKutta1412Coefficients.aCoefficients( 21, 14 ) = -1.10563611857482440905296961311590930801338308942637769555540;
    rungeKutta1412Coefficients.aCoefficients( 21, 15 ) = 176.1345918838113725878598980760556604069995167623016865882869129962911416096097878945;
    rungeKutta1412Coefficients.aCoefficients( 21, 16 ) = 0.4913848242148806622688983451644545574168846314027647925019604519368994965045299923826;
    rungeKutta1412Coefficients.aCoefficients( 21, 17 ) = -684.2780004498149443582375356108950819560771678936002751371799726829821841834791232605;
    rungeKutta1412Coefficients.aCoefficients( 21, 18 ) = 27.99106049983982589842243321243804074460025184006686868209688958109916979926727384229;
    rungeKutta1412Coefficients.aCoefficients( 21, 19 ) = 13.19397100302823334436709643711532384350641596237449753683872220663989495376087330358;
    rungeKutta1412Coefficients.aCoefficients( 21, 20 ) = 1.251287812839804454501149741480560063172688300773964063605141347518040989702499199856;
    
    rungeKutta1412Coefficients.aCoefficients( 22, 0 ) = -9.673079469481967636441261184332193958399514085718772596349277868068021458303626779169;
    rungeKutta1412Coefficients.aCoefficients( 22, 12 ) = -4.469901508585055314438462277019603604978306814087514357488023393670679083633020106516;
    rungeKutta1412Coefficients.aCoefficients( 22, 13 ) = 45.51271286909526819682419504000527511789059078173984816890412459840121969200961260987;
    rungeKutta1412Coefficients.aCoefficients( 22, 14 ) = -0.713085086183826912791492024438246129930559805352394367050813e-1;
    rungeKutta1412Coefficients.aCoefficients( 22, 15 ) = 11.22736140684127415825906244799393842078268007767944830815221105133516977144595052189;
    rungeKutta1412Coefficients.aCoefficients( 22, 16 ) = 0.1262443767176227245162379129091388093617868898191054263714925416869147773104813482457;
    rungeKutta1412Coefficients.aCoefficients( 22, 17 ) = -43.54393395494833136058106249072421076238143044676214056937881652359375369765457150165;
    rungeKutta1412Coefficients.aCoefficients( 22, 18 ) = 0.7871743075430589783987929949965509020645460914432340378113766124779028133099797867162;
    rungeKutta1412Coefficients.aCoefficients( 22, 19 ) = 0.5322646967446842156693007086038866907853957768215038536520118921656033723449302296244;
    rungeKutta1412Coefficients.aCoefficients( 22, 20 ) = 0.4224227339963253260102251274713887725750865388096033468497941673910509540050957057177;
    rungeKutta1412Coefficients.aCoefficients( 22, 21 ) = 0.8591312495030671073084380314998594434411150562941549563989586466154235621165245563192e-1;
    
    rungeKutta1412Coefficients.aCoefficients( 23, 0 ) = -10.06640324470547024033966069004268914722028247579687652710623604380152449409080444899;
    rungeKutta1412Coefficients.aCoefficients( 23, 8 ) = -0.3662532700490399702936857968489747917331190817335522078657913621382824038988807796287e-1;
    rungeKutta1412Coefficients.aCoefficients( 23, 9 ) = 0.3502549756362136819768494069798465243467890824711035742020654749717518291597210559354e-1;
    rungeKutta1412Coefficients.aCoefficients( 23, 10 ) = 0.3609460163621135089317866587583352398236899298642376718895880083960486970547825683491e-1;
    rungeKutta1412Coefficients.aCoefficients( 23, 11 ) = -0.2652199675536811063515959468346019236496270124574642848667252606942739787160130682669e-1;
    rungeKutta1412Coefficients.aCoefficients( 23, 12 ) = -6.270889721814641435905531494788716038393561229573960230194057818533161624674313994502;
    rungeKutta1412Coefficients.aCoefficients( 23, 13 ) = 48.20792374425629890907021030081950639234925931416361161278899187780407980462426656808;
    rungeKutta1412Coefficients.aCoefficients( 23, 14 ) = -0.694471689136165640882395180583732834557754169149088630301342e-1;
    rungeKutta1412Coefficients.aCoefficients( 23, 15 ) = 12.68106902048502956983413709136098070661084838114121251454273060707937017246509534894;
    rungeKutta1412Coefficients.aCoefficients( 23, 16 ) = 0.119671168968323754838161435501011294100927813964199613229864e-1;
    rungeKutta1412Coefficients.aCoefficients( 23, 17 ) = -46.72497649924824080033582682426626955932013216597956070401309263301039263373634230581;
    rungeKutta1412Coefficients.aCoefficients( 23, 18 ) = 1.330296133266267113147100392982165913990335111912271192356479099067512051132965697343;
    rungeKutta1412Coefficients.aCoefficients( 23, 19 ) = 1.007667875033982983534389036199266577711627177936617199056121787956529680139072027935;
    rungeKutta1412Coefficients.aCoefficients( 23, 20 ) = 0.2095120519336650916641223884754807028927707538644872411247284065032940106679251005781e-1;
    rungeKutta1412Coefficients.aCoefficients( 23, 21 ) = 0.2101347063312641773177354243313964074244121884437574908902263894855162847478911411134e-1;
    rungeKutta1412Coefficients.aCoefficients( 23, 22 ) = 0.9521960144171217941751015424545759073763602336583562405468424451848266905185171865534e-2;
    
    rungeKutta1412Coefficients.aCoefficients( 24, 0 ) = -409.4780816777437087725890974093703576244243416067520683455326035855162023776088699896;
    rungeKutta1412Coefficients.aCoefficients( 24, 8 ) = 0.2104519120236273856090970119990106557888074052256267000419050051487632641518018732685;
    rungeKutta1412Coefficients.aCoefficients( 24, 9 ) = 1.034274520572304119364829268288257099386679996983247401666929134177931632176349026735;
    rungeKutta1412Coefficients.aCoefficients( 24, 10 ) = 0.6003036458644224870512404482066405749390780924061569454673075686417142117164254262878e-2;
    rungeKutta1412Coefficients.aCoefficients( 24, 11 ) = 0.8559381250996195375780121060024077289150626526164160058172684354881277648341960563008;
    rungeKutta1412Coefficients.aCoefficients( 24, 12 ) = -250.5169985474478604927776577293161303865840504207820779326393997812026874735614210230;
    rungeKutta1412Coefficients.aCoefficients( 24, 13 ) = 1946.424666523884277660537503282647585958298508957614274560610260899186136259514015246;
    rungeKutta1412Coefficients.aCoefficients( 24, 14 ) = -3.045038821023103655061058090868608827869505440976021016842196622317831446605499698935;
    rungeKutta1412Coefficients.aCoefficients( 24, 15 ) = 490.6263795282817135212082652991680838415985422740616633051003594128766152337185220086;
    rungeKutta1412Coefficients.aCoefficients( 24, 16 ) = 1.566475895312709071154840670135974457395956152459667753199388690841173424714434871921;
    rungeKutta1412Coefficients.aCoefficients( 24, 17 ) = -1881.974289940111733622172673770358706192159066384530557689275696031792911993357071098;
    rungeKutta1412Coefficients.aCoefficients( 24, 18 ) = 75.25922247248471752788377136433031498216206189142459440229301807516615379972994062700;
    rungeKutta1412Coefficients.aCoefficients( 24, 19 ) = 34.57343569803310676224343447365546896967286447935510158001529990937243976348724448442;
    rungeKutta1412Coefficients.aCoefficients( 24, 20 ) = 3.211476794409689614354173618470737551690229667488916278855754113243135684398993410117;
    rungeKutta1412Coefficients.aCoefficients( 24, 21 ) = -0.4604080417384143913072014042370588488672450952653828208427296561415079214017074427602;
    rungeKutta1412Coefficients.aCoefficients( 24, 22 ) = -0.8707183398418105224318841379579862457242520473889365722145748143125162133630944128398e-1;
    rungeKutta1412Coefficients.aCoefficients( 24, 23 ) = -7.393518141583030675670169521955210639991857732491329543926346613193825315394087286297;
    
    rungeKutta1412Coefficients.aCoefficients( 25, 0 ) = 3.433474758535508789210934962575967811206238910720084588712755786644583035514752699598;
    rungeKutta1412Coefficients.aCoefficients( 25, 8 ) = 0.2491632048558174075389491488059951494598846535854176800982219995075912885766744587193e-2;
    rungeKutta1412Coefficients.aCoefficients( 25, 9 ) = 0.2301387878545931496383998463737427687720871226381422342236583655735620108657836993957e-1;
    rungeKutta1412Coefficients.aCoefficients( 25, 10 ) = -0.3221559566929770987244760924671208781894636047606204610433085107190031098987004938258e-2;
    rungeKutta1412Coefficients.aCoefficients( 25, 11 ) = 0.9884425494476646689463354144878852560408199827860146481292993078049373245839618405001e-2;
    rungeKutta1412Coefficients.aCoefficients( 25, 12 ) = 2.162527993779225077883078419047573540457592253357327094851479956564246957314476133478;
    rungeKutta1412Coefficients.aCoefficients( 25, 13 ) = -16.26998645464574213280656406601394890069875520402288517985775075363232756881970486667;
    rungeKutta1412Coefficients.aCoefficients( 25, 14 ) = -0.1285345021205245528435834174709350105380290375426545062302651848844352856037884822181;
    rungeKutta1412Coefficients.aCoefficients( 25, 15 ) = -8.98915042666504253089307820833379330486511746063552853023189;
    rungeKutta1412Coefficients.aCoefficients( 25, 16 ) = -0.3485953632320253333870802018510136501924017672505137649688730136175086767654181319387e-2;
    rungeKutta1412Coefficients.aCoefficients( 25, 17 ) = 15.79361941133398075362351873886955741358533870251397376656158275266140525531011608606;
    rungeKutta1412Coefficients.aCoefficients( 25, 18 ) = -0.5744033309140950656281654820173358201483836631956754708231458398423255984252281047127;
    rungeKutta1412Coefficients.aCoefficients( 25, 19 ) = -0.3456020390213932966927224966081249825352372288276553067081833889419898565070467534157;
    rungeKutta1412Coefficients.aCoefficients( 25, 20 ) = -0.6622414902065850917316199913837577811330679927074186873906450413385445874036001388495e-2;
    rungeKutta1412Coefficients.aCoefficients( 25, 21 ) = -0.7777881292422041640325464586073643097593472096267591120155367761150273183248441708392e-2;
    rungeKutta1412Coefficients.aCoefficients( 25, 22 ) = -0.3560841924022749133388272326974373646752408187917065879526063406092336300493607300593e-2;
    rungeKutta1412Coefficients.aCoefficients( 25, 23 ) = 4.792825064499307996497977496298401894572969341393590555417712618624354747222657791607;
    rungeKutta1412Coefficients.aCoefficients( 25, 24 ) = 0.153725464873068577844576387402512082757034273069877432944621;
    
    rungeKutta1412Coefficients.aCoefficients( 26, 0 ) = 32.30385208719854423269947344400315350913649750477846297617061421719281146058139852238;
    rungeKutta1412Coefficients.aCoefficients( 26, 5 ) = -0.317987696266205093901912847692712407988609169703103952205634e-2;
    rungeKutta1412Coefficients.aCoefficients( 26, 6 ) = 0.8063977149061920772608217115203795063935431115674197501197468839656405367779525213500;
    rungeKutta1412Coefficients.aCoefficients( 26, 7 ) = 0.9759831264123889790935228506842888513146720480030545503571875185550549213299958241991e-1;
    rungeKutta1412Coefficients.aCoefficients( 26, 8 ) = 0.7785755781583989090275124464529272389997634605941819649588520345133050850477185489203;
    rungeKutta1412Coefficients.aCoefficients( 26, 9 ) = 0.2048904238315994281894992020981056033120292350814206535748293420400885242747823516625;
    rungeKutta1412Coefficients.aCoefficients( 26, 10 ) = -1.562615796274681883070709439505278252114628922364243608928053762634922556160297217820;
    rungeKutta1412Coefficients.aCoefficients( 26, 12 ) = 16.34298918823105706485042439739271747087533535041545512917666902744198799725970841669;
    rungeKutta1412Coefficients.aCoefficients( 26, 13 ) = -154.5445552935436212307301896314710363993166836696091165017078152549564923882084122674;
    rungeKutta1412Coefficients.aCoefficients( 26, 14 ) = 1.569710887033348726920342834176217614662635935824970859658624964687079589089479471888;
    rungeKutta1412Coefficients.aCoefficients( 26, 15 ) = 3.276855450872481313214298172699007311655224049747336000450385269517693130775985884604;
    rungeKutta1412Coefficients.aCoefficients( 26, 16 ) = -0.5034892451936531763480407271997836265340810956916323972462042700071863164675818955838e-1;
    rungeKutta1412Coefficients.aCoefficients( 26, 17 ) = 153.3211518580416650705937678859146940112243631025945564907021486707139114294996134941;
    rungeKutta1412Coefficients.aCoefficients( 26, 18 ) = 7.175681863277204958467664848147841435678263080348653386540185145833155908488128910568;
    rungeKutta1412Coefficients.aCoefficients( 26, 19 ) = -2.940367486753004819459176598969309892153205943807775979427615740476908865098135595635;
    rungeKutta1412Coefficients.aCoefficients( 26, 20 ) = -0.6658459460768031444707496760226288702819204931972568878708744783028558369468497032253e-1;
    rungeKutta1412Coefficients.aCoefficients( 26, 21 ) = -0.4623460549908436612292486685622172611769665140168592842374268449140643068786760618896e-1;
    rungeKutta1412Coefficients.aCoefficients( 26, 22 ) = -0.2041987335856794015393882286172697788485797748215817776751235910664984352284968100100e-1;
    rungeKutta1412Coefficients.aCoefficients( 26, 23 ) = -53.35231064387358505159534411659981079740450904957915977996876390672711239156977103431;
    rungeKutta1412Coefficients.aCoefficients( 26, 24 ) = -1.355487147150786549787321867059964040175545016141913251148206738329360142936656282958;
    rungeKutta1412Coefficients.aCoefficients( 26, 25 ) = -1.571962758012327518829017351714592491776872191144425834618663282570958684038698495739;
    
    rungeKutta1412Coefficients.aCoefficients( 27, 0 ) = -16.64514674863415128720312944039317587645603711308189782044257016154825923946758475845;
    rungeKutta1412Coefficients.aCoefficients( 27, 5 ) = 0.5922327803245033080429900057980465247383895604442571368349896773084347972825775455007e-2;
    rungeKutta1412Coefficients.aCoefficients( 27, 6 ) = 0.4703261599638411122172243032058941134553625307461088250108483236601604516650193568134;
    rungeKutta1412Coefficients.aCoefficients( 27, 7 ) = 0.2996888638486790008539818370961923991368311216717812791841936858888827504094204242461;
    rungeKutta1412Coefficients.aCoefficients( 27, 8 ) = -0.2476568775939949146899922763298108258539580692639470955481886317480090967647905771626;
    rungeKutta1412Coefficients.aCoefficients( 27, 9 ) = 0.1108950297714376828939998518390617145224451736006787182086245987785252503880550245038;
    rungeKutta1412Coefficients.aCoefficients( 27, 11 ) = -0.4917190438462291470706666287041940976780819072106730449888664749836403474888832394921;
    rungeKutta1412Coefficients.aCoefficients( 27, 12 ) = -11.47431544272894969683894925643525363508424541308531757856483965863898534849416840511;
    rungeKutta1412Coefficients.aCoefficients( 27, 13 ) = 80.25931665762302725417024858864844001527933666235899875893849400507278534931158408231;
    rungeKutta1412Coefficients.aCoefficients( 27, 14 ) = -0.3841323039800428476253125267590291037469268413420882192068133107492120348263618466046;
    rungeKutta1412Coefficients.aCoefficients( 27, 15 ) = 7.281476674681075834713269509261361157676125818628777243483988994104498714011047355205;
    rungeKutta1412Coefficients.aCoefficients( 27, 16 ) = -0.1326993846122483795105717081760352748368273416167518843018178653526280269065470590467;
    rungeKutta1412Coefficients.aCoefficients( 27, 17 ) = -81.07998325257307266746792897522552400060707166336329885641562357237166810196760593013;
    rungeKutta1412Coefficients.aCoefficients( 27, 18 ) = -1.250374928356206395217681856561791199622537474924031863192434629401819729868852090550;
    rungeKutta1412Coefficients.aCoefficients( 27, 19 ) = 2.592635949695436810237763795043773249942264473592968880837586883560068434349818491911;
    rungeKutta1412Coefficients.aCoefficients( 27, 20 ) = -0.3014402983464045398301639972605268752644315372756414953420797074457552586137488110716;
    rungeKutta1412Coefficients.aCoefficients( 27, 21 ) = 0.2213844607898323374517064515727737916952468390573184143179573617704323166985265217363;
    rungeKutta1412Coefficients.aCoefficients( 27, 22 ) = 0.8275772747718929319559898709746931529962764354298098905497078729734353980896315305691e-1;
    rungeKutta1412Coefficients.aCoefficients( 27, 23 ) = 18.99606620406115204646724500372432639981751614122371589366718674999943569769696943522;
    rungeKutta1412Coefficients.aCoefficients( 27, 24 ) = 0.2692319464096396856234680151283341674600519103489128451211866688910668614577677735665;
    rungeKutta1412Coefficients.aCoefficients( 27, 25 ) = 1.626748274470665374629893649296289339881250292841836802790201430504847697803528636395;
    rungeKutta1412Coefficients.aCoefficients( 27, 26 ) = 0.4917190438462291470706666287041940976780819072106730449888664749836403474888832394921;
    
    rungeKutta1412Coefficients.aCoefficients( 28, 0 ) = 0.8384798124090526646169687913728140859805331392249111310693346670107922625197375034871e-1;
    rungeKutta1412Coefficients.aCoefficients( 28, 5 ) = -0.1179493671009738143197550560312957753679619605907361507776128268875265788248790903515e-1;
    rungeKutta1412Coefficients.aCoefficients( 28, 6 ) = -0.2472990205688126523394738387431945983259928403533401326974984247503501083158412965835;
    rungeKutta1412Coefficients.aCoefficients( 28, 7 ) = 0.9780808583677290122593130140812916655037406554767339407565991037499621093437371932341e-1;
    rungeKutta1412Coefficients.aCoefficients( 28, 8 ) = 0.2175906892434206313600086517678603183441681200247821768799893467069296630467914197921;
    rungeKutta1412Coefficients.aCoefficients( 28, 10 ) = 0.1375856067633252248656596321967877466474472229750848659754400903987833771639575727867;
    rungeKutta1412Coefficients.aCoefficients( 28, 11 ) = 0.4398702297150466850587900923415450260461038902942613590425808839943205635447284745074e-1;
    rungeKutta1412Coefficients.aCoefficients( 28, 13 ) = -0.5137008137681933419570044566186303037387573636419640300869712169933398305905931343468;
    rungeKutta1412Coefficients.aCoefficients( 28, 14 ) = 0.8263556911513155086442113083991534587014231586161685769224194977471882335420141183213;
    rungeKutta1412Coefficients.aCoefficients( 28, 15 ) = 25.70181397198118326258738829725199395111365563419600781824702737091645129169813134401;
    rungeKutta1412Coefficients.aCoefficients( 28, 23 ) = -25.70181397198118326258738829725199395111365563419600781824702737091645129169813134401;
    rungeKutta1412Coefficients.aCoefficients( 28, 24 ) = -0.8263556911513155086442113083991534587014231586161685769224194977471882335420141183213;
    rungeKutta1412Coefficients.aCoefficients( 28, 25 ) = 0.5137008137681933419570044566186303037387573636419640300869712169933398305905931343468;
    rungeKutta1412Coefficients.aCoefficients( 28, 26 ) = -0.4398702297150466850587900923415450260461038902942613590425808839943205635447284745074e-1;
    rungeKutta1412Coefficients.aCoefficients( 28, 27 ) = -0.1375856067633252248656596321967877466474472229750848659754400903987833771639575727867;
    
    rungeKutta1412Coefficients.aCoefficients( 29, 0 ) = 0.1243805266540944128815164208687993162684914663596714231632892354628068537117612942798;
    rungeKutta1412Coefficients.aCoefficients( 29, 4 ) = 0.2261202821975843014222386629792029011967523207426331439651447460281196206643404356021;
    rungeKutta1412Coefficients.aCoefficients( 29, 5 ) = 0.1378858876180808806076958370164778145309694174914933853635428709475288586061552782365e-1;
    rungeKutta1412Coefficients.aCoefficients( 29, 6 ) = -0.6722101339966844497493995074143058569500863415253821828561997825320849038679063596730e-1;
    rungeKutta1412Coefficients.aCoefficients( 29, 9 ) = -0.8562389750854283547553497698795017721121215974115638028550665385850612741040225222977;
    rungeKutta1412Coefficients.aCoefficients( 29, 10 ) = -1.963375228668589089282628500280938139881804405182674045535756631526916950083353845169;
    rungeKutta1412Coefficients.aCoefficients( 29, 11 ) = -0.2323328227241194012372462573089218472501081992304199949782180319905262045718872259601;
    rungeKutta1412Coefficients.aCoefficients( 29, 13 ) = 4.306607190864533494616689368765629477724325620534780926267640393608500758570100495873;
    rungeKutta1412Coefficients.aCoefficients( 29, 14 ) = -2.927229632494654826597879112023904466876873949506336126307786635262992367484998786517;
    rungeKutta1412Coefficients.aCoefficients( 29, 15 ) = -82.31316663978589444544923341054587077357619664281386893950601309356417181948645997040;
    rungeKutta1412Coefficients.aCoefficients( 29, 23 ) = 82.31316663978589444544923341054587077357619664281386893950601309356417181948645997040;
    rungeKutta1412Coefficients.aCoefficients( 29, 24 ) = 2.927229632494654826597879112023904466876873949506336126307786635262992367484998786517;
    rungeKutta1412Coefficients.aCoefficients( 29, 25 ) = -4.306607190864533494616689368765629477724325620534780926267640393608500758570100495873;
    rungeKutta1412Coefficients.aCoefficients( 29, 26 ) = 0.2323328227241194012372462573089218472501081992304199949782180319905262045718872259601;
    rungeKutta1412Coefficients.aCoefficients( 29, 27 ) = 1.963375228668589089282628500280938139881804405182674045535756631526916950083353845169;
    rungeKutta1412Coefficients.aCoefficients( 29, 28 ) = 0.8562389750854283547553497698795017721121215974115638028550665385850612741040225222977;
    
    rungeKutta1412Coefficients.aCoefficients( 30, 0 ) = 0.1034845616366797766729935465119103444997447982019713166066629728281981965079290745983;
    rungeKutta1412Coefficients.aCoefficients( 30, 3 ) = 0.1220688873064072225896440828689620771395927148341621347412746563709055937325311521675;
    rungeKutta1412Coefficients.aCoefficients( 30, 4 ) = 0.4825744903312466224751347801256881128659190238501680496794015023696413273862321544150;
    rungeKutta1412Coefficients.aCoefficients( 30, 5 ) = -0.3814096000156069997308862400056202056641130724784114774219699240039767479629669855696e-1;
    rungeKutta1412Coefficients.aCoefficients( 30, 7 ) = -0.5504995253108023241383885070205081774114143110000375617128363206424473498745141065969;
    rungeKutta1412Coefficients.aCoefficients( 30, 9 ) = -0.7119158115851892278876482620437943875782918824067455704957652139710574799878630163853;
    rungeKutta1412Coefficients.aCoefficients( 30, 10 ) = -0.5841296056715513404329887301584808720953353296452275957070524410065417676683463009109;
    rungeKutta1412Coefficients.aCoefficients( 30, 13 ) = 2.110463081258649321287173000466227503003750542789369878507182287710881470618943318741;
    rungeKutta1412Coefficients.aCoefficients( 30, 14 ) = -0.8374947367395721355257420230010379926952601753351235177405529298334532793741463162845e-1;
    rungeKutta1412Coefficients.aCoefficients( 30, 15 ) = 5.100214990723209140752959690433441131075450608628042491597346388445135412965217165555;
    rungeKutta1412Coefficients.aCoefficients( 30, 23 ) = -5.100214990723209140752959690433441131075450608628042491597346388445135412965217165555;
    rungeKutta1412Coefficients.aCoefficients( 30, 24 ) = 0.8374947367395721355257420230010379926952601753351235177405529298334532793741463162845e-1;
    rungeKutta1412Coefficients.aCoefficients( 30, 25 ) = -2.110463081258649321287173000466227503003750542789369878507182287710881470618943318741;
    rungeKutta1412Coefficients.aCoefficients( 30, 27 ) = 0.5841296056715513404329887301584808720953353296452275957070524410065417676683463009109;
    rungeKutta1412Coefficients.aCoefficients( 30, 28 ) = 0.7119158115851892278876482620437943875782918824067455704957652139710574799878630163853;
    rungeKutta1412Coefficients.aCoefficients( 30, 29 ) = 0.5504995253108023241383885070205081774114143110000375617128363206424473498745141065969;
    
    rungeKutta1412Coefficients.aCoefficients( 31, 0 ) = 0.1933333333333333333333333333333333333333333333333333333333333333333333333333333333333;
    rungeKutta1412Coefficients.aCoefficients( 31, 2 ) = 0.22;
    rungeKutta1412Coefficients.aCoefficients( 31, 3 ) = -0.8e-1;
    rungeKutta1412Coefficients.aCoefficients( 31, 6 ) = 0.1099934255807247039194624048650683408451190582958464264636524271459687549994002654752;
    rungeKutta1412Coefficients.aCoefficients( 31, 7 ) = -0.2542970480762701613840685069971531221418356269767039208462421656164179875269042982442;
    rungeKutta1412Coefficients.aCoefficients( 31, 9 ) = 0.8655707771166942543437703438210982818328474012330118593467368132762510892051242759318;
    rungeKutta1412Coefficients.aCoefficients( 31, 10 ) = 3.324164491140930831067995527865720183368600929369864071601998386039920635781409865040;
    rungeKutta1412Coefficients.aCoefficients( 31, 13 ) = -12.01022233159779338823523851486618412603019426339968151272769528462035002110216728101;
    rungeKutta1412Coefficients.aCoefficients( 31, 14 ) = 0.4766014662424932394304427768620618996029637820035802094825720242694315551196576125507;
    rungeKutta1412Coefficients.aCoefficients( 31, 15 ) = -29.02430112210363905258026232136540995962512213324709106915239870601916450708546744075;
    rungeKutta1412Coefficients.aCoefficients( 31, 23 ) = 29.02430112210363905258026232136540995962512213324709106915239870601916450708546744075;
    rungeKutta1412Coefficients.aCoefficients( 31, 24 ) = -0.4766014662424932394304427768620618996029637820035802094825720242694315551196576125507;
    rungeKutta1412Coefficients.aCoefficients( 31, 25 ) = 12.01022233159779338823523851486618412603019426339968151272769528462035002110216728101;
    rungeKutta1412Coefficients.aCoefficients( 31, 27 ) = -3.324164491140930831067995527865720183368600929369864071601998386039920635781409865040;
    rungeKutta1412Coefficients.aCoefficients( 31, 28 ) = -0.8655707771166942543437703438210982818328474012330118593467368132762510892051242759318;
    rungeKutta1412Coefficients.aCoefficients( 31, 29 ) = 0.2542970480762701613840685069971531221418356269767039208462421656164179875269042982442;
    rungeKutta1412Coefficients.aCoefficients( 31, 30 ) = -0.1099934255807247039194624048650683408451190582958464264636524271459687549994002654752;
    
    rungeKutta1412Coefficients.aCoefficients( 32, 0 ) = -0.8333333333333333333333333333333333333333333333333333333333333333333333333333333333333;
    rungeKutta1412Coefficients.aCoefficients( 32, 1 ) = 1.388888888888888888888888888888888888888888888888888888888888888888888888888888888889;
    rungeKutta1412Coefficients.aCoefficients( 32, 4 ) = -0.75;
    rungeKutta1412Coefficients.aCoefficients( 32, 6 ) = -0.4925295437180263044226820491140213202002146815806577847190740839644346370048749342561;
    rungeKutta1412Coefficients.aCoefficients( 32, 30 ) = 0.4925295437180263044226820491140213202002146815806577847190740839644346370048749342561;
    rungeKutta1412Coefficients.aCoefficients( 32, 31 ) = 0.75;
    
    rungeKutta1412Coefficients.aCoefficients( 33, 0 ) = 0.1111111111111111111111111111111111111111111111111111111111111111111111111111111111111;
    rungeKutta1412Coefficients.aCoefficients( 33, 2 ) = -0.2222222222222222222222222222222222222222222222222222222222222222222222222222222222222;
    rungeKutta1412Coefficients.aCoefficients( 33, 32 ) = 0.2222222222222222222222222222222222222222222222222222222222222222222222222222222222222;
    
    rungeKutta1412Coefficients.aCoefficients( 34, 0 ) = 0.2858351403889715587960888421638364148529275378945964668924322897553490152559792262023;
    rungeKutta1412Coefficients.aCoefficients( 34, 1 ) = 0.2916666666666666666666666666666666666666666666666666666666666666666666666666666666667;
    rungeKutta1412Coefficients.aCoefficients( 34, 2 ) = 0.21875;
    rungeKutta1412Coefficients.aCoefficients( 34, 4 ) = 0.1640625;
    rungeKutta1412Coefficients.aCoefficients( 34, 6 ) = 0.2181943549455566583271882415813521070932888243221879411415164327116967439531911272777;
    rungeKutta1412Coefficients.aCoefficients( 34, 7 ) = 0.1803928984786977668636352219467754377196200536418492285624347210514163759703679527180;
    rungeKutta1412Coefficients.aCoefficients( 34, 9 ) = 0.2057138394048450188591207551229295422775700949828089053939914789386228504942804843989;
    rungeKutta1412Coefficients.aCoefficients( 34, 10 ) = 0.2427157915817702399702829279594465157627459713866705419485763522859549196625913978401;
    rungeKutta1412Coefficients.aCoefficients( 34, 11 ) = 0.2464657808136293058336092911818914077992281038693057051370210135284213379790417930740;
    rungeKutta1412Coefficients.aCoefficients( 34, 12 ) = -3.449919407908908249798341546016226620603704606149316442883265523381128452524989278943;
    rungeKutta1412Coefficients.aCoefficients( 34, 13 ) = 0.2288755621600360817607290607384585842942203725527402184592948392511281334278617959957;
    rungeKutta1412Coefficients.aCoefficients( 34, 14 ) = 0.2832905997021514153215274190567333359784365954938557898314048426595070708424182066065;
    rungeKutta1412Coefficients.aCoefficients( 34, 15 ) = 3.210851258377666409601314905442367870055573203322387098512984999880577120008173123283;
    rungeKutta1412Coefficients.aCoefficients( 34, 16 ) = -0.2235387773648456999202337562141625079641252300836740320899016275445898395177373582441;
    rungeKutta1412Coefficients.aCoefficients( 34, 17 ) = -0.7071211572044190735187272862074872121300912319552061607910521928571247612111795934106;
    rungeKutta1412Coefficients.aCoefficients( 34, 18 ) = 3.211233451502870804081747292028565008932600344430223743249588034157195885590228893622;
    rungeKutta1412Coefficients.aCoefficients( 34, 19 ) = 1.409543483096697660304144743011231757690459455735489863573218752821178310978199657967;
    rungeKutta1412Coefficients.aCoefficients( 34, 20 ) = -0.1513620534437426131216022767425181110909630262036760559494590353712667648924754181285;
    rungeKutta1412Coefficients.aCoefficients( 34, 21 ) = 0.3723505745270142764547240802146199843971210282021482987373568243836683323798121465643;
    rungeKutta1412Coefficients.aCoefficients( 34, 22 ) = 0.2529787464063613367221999077621412859157757281294143192610824780367182739421617243696;
    rungeKutta1412Coefficients.aCoefficients( 34, 23 ) = -3.210851258377666409601314905442367870055573203322387098512984999880577120008173123283;
    rungeKutta1412Coefficients.aCoefficients( 34, 24 ) = -0.2832905997021514153215274190567333359784365954938557898314048426595070708424182066065;
    rungeKutta1412Coefficients.aCoefficients( 34, 25 ) = -0.2288755621600360817607290607384585842942203725527402184592948392511281334278617959957;
    rungeKutta1412Coefficients.aCoefficients( 34, 26 ) = -0.2464657808136293058336092911818914077992281038693057051370210135284213379790417930740;
    rungeKutta1412Coefficients.aCoefficients( 34, 27 ) = -0.2427157915817702399702829279594465157627459713866705419485763522859549196625913978401;
    rungeKutta1412Coefficients.aCoefficients( 34, 28 ) = -0.2057138394048450188591207551229295422775700949828089053939914789386228504942804843989;
    rungeKutta1412Coefficients.aCoefficients( 34, 29 ) = -0.1803928984786977668636352219467754377196200536418492285624347210514163759703679527180;
    rungeKutta1412Coefficients.aCoefficients( 34, 30 ) = -0.2181943549455566583271882415813521070932888243221879411415164327116967439531911272777;
    rungeKutta1412Coefficients.aCoefficients( 34, 31 ) = -0.1640625;
    rungeKutta1412Coefficients.aCoefficients( 34, 32 ) = -0.21875;
    rungeKutta1412Coefficients.aCoefficients( 34, 33 ) = -0.2916666666666666666666666666666666666666666666666666666666666666666666666666666666667;
    
    // Define c-coefficients for the Runge-Kutta method.
    rungeKutta1412Coefficients.cCoefficients = Eigen::VectorXd::Zero( 35 );
    rungeKutta1412Coefficients.cCoefficients( 1 ) = 0.1111111111111111111111111111111111111111111111111111111111111111111111111111111111111; 
    rungeKutta1412Coefficients.cCoefficients( 2 ) = 0.5555555555555555555555555555555555555555555555555555555555555555555555555555555555556; 
    rungeKutta1412Coefficients.cCoefficients( 3 ) = 0.8333333333333333333333333333333333333333333333333333333333333333333333333333333333333; 
    rungeKutta1412Coefficients.cCoefficients( 4 ) = 0.3333333333333333333333333333333333333333333333333333333333333333333333333333333333333; 
    rungeKutta1412Coefficients.cCoefficients( 5 ) = 1.0; 
    rungeKutta1412Coefficients.cCoefficients( 6 ) = 0.6699869792727729217646837855059985139388452296384603532851421391683474428303956826239; 
    rungeKutta1412Coefficients.cCoefficients( 7 ) = 0.2970683842138183573895847168082194132233320946989156873791682903324708698499266217383; 
    rungeKutta1412Coefficients.cCoefficients( 8 ) = 0.7272727272727272727272727272727272727272727272727272727272727272727272727272727272727; 
    rungeKutta1412Coefficients.cCoefficients( 9 ) = 0.1401527990421887652761874879669467176298064630825329362873230163439023340348096838456; 
    rungeKutta1412Coefficients.cCoefficients( 10 ) = 0.7007010397701507371510998548307493379414070492655464089692218490447945746638665522966; 
    rungeKutta1412Coefficients.cCoefficients( 11 ) = 0.3636363636363636363636363636363636363636363636363636363636363636363636363636363636364; 
    rungeKutta1412Coefficients.cCoefficients( 12 ) = 0.2631578947368421052631578947368421052631578947368421052631578947368421052631578947368; 
    rungeKutta1412Coefficients.cCoefficients( 13 ) = 0.392172246650270859125196642501208648863714315266128052078483e-1; 
    rungeKutta1412Coefficients.cCoefficients( 14 ) = 0.8129175029283767629833931592780365061896123726172385507744269795906758195776958783707; 
    rungeKutta1412Coefficients.cCoefficients( 15 ) = 0.1666666666666666666666666666666666666666666666666666666666666666666666666666666666667; 
    rungeKutta1412Coefficients.cCoefficients( 16 ) = 0.9; 
    rungeKutta1412Coefficients.cCoefficients( 17 ) = 0.6412992574519669233127711938966828094810966516150832254029235721305050295351572963693e-1; 
    rungeKutta1412Coefficients.cCoefficients( 18 ) = 0.2041499092834288489277446343010234050271495052413337516288702042649259099754335560687; 
    rungeKutta1412Coefficients.cCoefficients( 19 ) = 0.3953503910487605656156713698273243723522272974566594505545766538389345381768585023057; 
    rungeKutta1412Coefficients.cCoefficients( 20 ) = 0.6046496089512394343843286301726756276477727025433405494454233461610654618231414976943; 
    rungeKutta1412Coefficients.cCoefficients( 21 ) = 0.7958500907165711510722553656989765949728504947586662483711297957350740900245664439313; 
    rungeKutta1412Coefficients.cCoefficients( 22 ) = 0.9358700742548033076687228806103317190518903348384916774597076427869494970464842703631; 
    rungeKutta1412Coefficients.cCoefficients( 23 ) = 0.1666666666666666666666666666666666666666666666666666666666666666666666666666666666667; 
    rungeKutta1412Coefficients.cCoefficients( 24 ) = 0.8129175029283767629833931592780365061896123726172385507744269795906758195776958783707; 
    rungeKutta1412Coefficients.cCoefficients( 25 ) = 0.392172246650270859125196642501208648863714315266128052078483e-1; 
    rungeKutta1412Coefficients.cCoefficients( 26 ) = 0.3636363636363636363636363636363636363636363636363636363636363636363636363636363636364; 
    rungeKutta1412Coefficients.cCoefficients( 27 ) = 0.7007010397701507371510998548307493379414070492655464089692218490447945746638665522966; 
    rungeKutta1412Coefficients.cCoefficients( 28 ) = 0.1401527990421887652761874879669467176298064630825329362873230163439023340348096838456; 
    rungeKutta1412Coefficients.cCoefficients( 29 ) = 0.2970683842138183573895847168082194132233320946989156873791682903324708698499266217383; 
    rungeKutta1412Coefficients.cCoefficients( 30 ) = 0.6699869792727729217646837855059985139388452296384603532851421391683474428303956826239;
    rungeKutta1412Coefficients.cCoefficients( 31 ) = 0.3333333333333333333333333333333333333333333333333333333333333333333333333333333333333; 
    rungeKutta1412Coefficients.cCoefficients( 32 ) = 0.5555555555555555555555555555555555555555555555555555555555555555555555555555555555556;
    rungeKutta1412Coefficients.cCoefficients( 33 ) = 0.1111111111111111111111111111111111111111111111111111111111111111111111111111111111111;
    rungeKutta1412Coefficients.cCoefficients( 34 ) = 1.0;
    
    // Define b-coefficients for the Runge-Kutta method.
    rungeKutta1412Coefficients.bCoefficients = Eigen::MatrixXd::Zero( 2, 35 );    
    rungeKutta1412Coefficients.bCoefficients( 0, 0 ) = 0.1785714285714285714285714285714285714285714285714285714285714285714285714285714285714e-1;
    rungeKutta1412Coefficients.bCoefficients( 0, 1 ) = 439.0 / 64000.0;
    rungeKutta1412Coefficients.bCoefficients( 0, 2 ) = 0.1171875e-1;
    rungeKutta1412Coefficients.bCoefficients( 0, 4 ) = 0.17578125e-1;
    rungeKutta1412Coefficients.bCoefficients( 0, 6 ) = 0.234375e-1;
    rungeKutta1412Coefficients.bCoefficients( 0, 7 ) = 0.29296875e-1;
    rungeKutta1412Coefficients.bCoefficients( 0, 9 ) = 0.3515625e-1;
    rungeKutta1412Coefficients.bCoefficients( 0, 10 ) = 0.41015625e-1;
    rungeKutta1412Coefficients.bCoefficients( 0, 11 ) = 0.46875e-1;
    rungeKutta1412Coefficients.bCoefficients( 0, 13 ) = 0.52734375e-1;
    rungeKutta1412Coefficients.bCoefficients( 0, 14 ) = 0.5859375e-1;
    rungeKutta1412Coefficients.bCoefficients( 0, 15 ) = 0.64453125e-1;
    rungeKutta1412Coefficients.bCoefficients( 0, 17 ) = 0.1053521135717530196914960328878781622276730830805238840416702908213176249782427570033;
    rungeKutta1412Coefficients.bCoefficients( 0, 18 ) = 0.1705613462417521823821203385538740858875554878027908047375010369442754416180982144816;
    rungeKutta1412Coefficients.bCoefficients( 0, 19 ) = 0.2062293973293519407835264857011048947419142862595424540779715293772640762608018856579;
    rungeKutta1412Coefficients.bCoefficients( 0, 20 ) = 0.2062293973293519407835264857011048947419142862595424540779715293772640762608018856579;
    rungeKutta1412Coefficients.bCoefficients( 0, 21 ) = 0.1705613462417521823821203385538740858875554878027908047375010369442754416180982144816;
    rungeKutta1412Coefficients.bCoefficients( 0, 22 ) = 0.1053521135717530196914960328878781622276730830805238840416702908213176249782427570033;
    rungeKutta1412Coefficients.bCoefficients( 0, 23 ) = -0.64453125e-1;
    rungeKutta1412Coefficients.bCoefficients( 0, 24 ) = -0.5859375e-1;
    rungeKutta1412Coefficients.bCoefficients( 0, 25 ) = -0.52734375e-1;
    rungeKutta1412Coefficients.bCoefficients( 0, 26 ) = -0.46875e-1;
    rungeKutta1412Coefficients.bCoefficients( 0, 27 ) = -0.41015625e-1;
    rungeKutta1412Coefficients.bCoefficients( 0, 28 ) = -0.3515625e-1;
    rungeKutta1412Coefficients.bCoefficients( 0, 29 ) = -0.29296875e-1;
    rungeKutta1412Coefficients.bCoefficients( 0, 30 ) = -0.234375e-1;
    rungeKutta1412Coefficients.bCoefficients( 0, 31 ) = -0.17578125e-1;
    rungeKutta1412Coefficients.bCoefficients( 0, 32 ) = -0.1171875e-1;
    rungeKutta1412Coefficients.bCoefficients( 0, 33 ) = -439.0 / 64000.0;
    rungeKutta1412Coefficients.bCoefficients( 0, 34 ) = 0.1785714285714285714285714285714285714285714285714285714285714285714285714285714285714e-1;
    
    rungeKutta1412Coefficients.bCoefficients( 1, 0 ) = 0.1785714285714285714285714285714285714285714285714285714285714285714285714285714285714e-1;
    rungeKutta1412Coefficients.bCoefficients( 1, 1 ) = 0.5859375e-2;
    rungeKutta1412Coefficients.bCoefficients( 1, 2 ) = 0.1171875e-1;
    rungeKutta1412Coefficients.bCoefficients( 1, 4 ) = 0.17578125e-1;
    rungeKutta1412Coefficients.bCoefficients( 1, 6 ) = 0.234375e-1;
    rungeKutta1412Coefficients.bCoefficients( 1, 7 ) = 0.29296875e-1;
    rungeKutta1412Coefficients.bCoefficients( 1, 9 ) = 0.3515625e-1;
    rungeKutta1412Coefficients.bCoefficients( 1, 10 ) = 0.41015625e-1;
    rungeKutta1412Coefficients.bCoefficients( 1, 11 ) = 0.46875e-1;
    rungeKutta1412Coefficients.bCoefficients( 1, 13 ) = 0.52734375e-1;
    rungeKutta1412Coefficients.bCoefficients( 1, 14 ) = 0.5859375e-1;
    rungeKutta1412Coefficients.bCoefficients( 1, 15 ) = 0.64453125e-1;
    rungeKutta1412Coefficients.bCoefficients( 1, 17 ) = 0.1053521135717530196914960328878781622276730830805238840416702908213176249782427570033;
    rungeKutta1412Coefficients.bCoefficients( 1, 18 ) = 0.1705613462417521823821203385538740858875554878027908047375010369442754416180982144816;
    rungeKutta1412Coefficients.bCoefficients( 1, 19 ) = 0.2062293973293519407835264857011048947419142862595424540779715293772640762608018856579;
    rungeKutta1412Coefficients.bCoefficients( 1, 20 ) = 0.2062293973293519407835264857011048947419142862595424540779715293772640762608018856579;
    rungeKutta1412Coefficients.bCoefficients( 1, 21 ) = 0.1705613462417521823821203385538740858875554878027908047375010369442754416180982144816;
    rungeKutta1412Coefficients.bCoefficients( 1, 22 ) = 0.1053521135717530196914960328878781622276730830805238840416702908213176249782427570033;
    rungeKutta1412Coefficients.bCoefficients( 1, 23 ) = -0.64453125e-1;
    rungeKutta1412Coefficients.bCoefficients( 1, 24 ) = -0.5859375e-1;
    rungeKutta1412Coefficients.bCoefficients( 1, 25 ) = -0.52734375e-1;
    rungeKutta1412Coefficients.bCoefficients( 1, 26 ) = -0.46875e-1;
    rungeKutta1412Coefficients.bCoefficients( 1, 27 ) = -0.41015625e-1;
    rungeKutta1412Coefficients.bCoefficients( 1, 28 ) = -0.3515625e-1;
    rungeKutta1412Coefficients.bCoefficients( 1, 29 ) = -0.29296875e-1;
    rungeKutta1412Coefficients.bCoefficients( 1, 30 ) = -0.234375e-1;
    rungeKutta1412Coefficients.bCoefficients( 1, 31 ) = -0.17578125e-1;
    rungeKutta1412Coefficients.bCoefficients( 1, 32 ) = -0.1171875e-1;
    rungeKutta1412Coefficients.bCoefficients( 1, 33 ) = -0.5859375e-2;
    rungeKutta1412Coefficients.bCoefficients( 1, 34 ) = 0.1785714285714285714285714285714285714285714285714285714285714285714285714285714285714e-1;

    // Define the name of the coefficient set.
    rungeKutta1412Coefficients.name = "Runge-Kutta-Feagin 14/12";
}

const RungeKuttaCoefficients& RungeKuttaCoefficients::get(
        CoefficientSets coefficientSet )
{
    static RungeKuttaCoefficients forwardEulerCoefficients,
                                  rungeKutta4Coefficients,
                                  explicitMidPointCoefficients,
                                  explicitTrapezoidRuleCoefficients,
                                  ralstonCoefficients,
                                  rungeKutta3Coefficients,
                                  ralston3Coefficients,
                                  SSPRK3Coefficients,
                                  ralston4Coefficients,
                                  threeEighthRuleRK4Coefficients,
                                  heunEulerCoefficients,
                                  rungeKuttaFehlberg12Coefficients,
                                  rungeKuttaFehlberg45Coefficients,
                                  rungeKuttaFehlberg56Coefficients,
                                  rungeKuttaFehlberg78Coefficients,
                                  rungeKutta87DormandPrinceCoefficients,
                                  rungeKuttaFehlberg89Coefficients,
                                  rungeKuttaVerner89Coefficients,
                                  rungeKuttaFeagin108Coefficients,
                                  rungeKuttaFeagin1210Coefficients,
                                  rungeKuttaFeagin1412Coefficients;

    switch ( coefficientSet )
    {
        
    case forwardEuler:
        if ( !(forwardEulerCoefficients.isFixedStepSize) )
        {
            initializeForwardEulerCoefficients( forwardEulerCoefficients );
        }
        return forwardEulerCoefficients;

    case rungeKutta4Classic:
        if ( !(rungeKutta4Coefficients.isFixedStepSize) )
        {
            initializeRungeKutta4Coefficients( rungeKutta4Coefficients );
        }
        return rungeKutta4Coefficients;

    case explicitMidPoint:
        if ( !(explicitMidPointCoefficients.isFixedStepSize) )
        {
            initializeExplicitMidpointCoefficients( explicitMidPointCoefficients );
        }
        return explicitMidPointCoefficients;

    case explicitTrapezoidRule:
        if ( !(explicitTrapezoidRuleCoefficients.isFixedStepSize) )
        {
            initializeExplicitTrapezoidRuleCoefficients( explicitTrapezoidRuleCoefficients );
        }
        return explicitTrapezoidRuleCoefficients;

    case ralston:
        if ( !(ralstonCoefficients.isFixedStepSize) )
        {
            initializeRalstonCoefficients( ralstonCoefficients );
        }
        return ralstonCoefficients;

    case rungeKutta3:
        if ( !(rungeKutta3Coefficients.isFixedStepSize) )
        {
            initializeRungeKutta3Coefficients( rungeKutta3Coefficients );
        }
        return rungeKutta3Coefficients;

    case ralston3:
        if ( !(ralston3Coefficients.isFixedStepSize) )
        {
            initializeRalston3Coefficients( ralston3Coefficients );
        }
        return ralston3Coefficients;

    case SSPRK3:
        if ( !(SSPRK3Coefficients.isFixedStepSize) )
        {
            initializeSSPRK3Coefficients( SSPRK3Coefficients );
        }
        return SSPRK3Coefficients;

    case ralston4:
        if ( !(ralston4Coefficients.isFixedStepSize) )
        {
            initializeRalston4Coefficients( ralston4Coefficients );
        }
        return ralston4Coefficients;

    case threeEighthRuleRK4:
        if ( !(threeEighthRuleRK4Coefficients.isFixedStepSize) )
        {
            initializeThreeEighthRuleRK4Coefficients( threeEighthRuleRK4Coefficients );
        }
        return threeEighthRuleRK4Coefficients;

    case heunEuler:
        if ( heunEulerCoefficients.higherOrder != 2 )
        {
            initializeHeunEulerCoefficients( heunEulerCoefficients );
        }
        return heunEulerCoefficients;

    case rungeKuttaFehlberg12:
        if ( rungeKuttaFehlberg12Coefficients.higherOrder != 2 )
        {
            initializeRungeKuttaFehlberg12Coefficients( rungeKuttaFehlberg12Coefficients );
        }
        return rungeKuttaFehlberg12Coefficients;

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
            initializeRungeKutta87DormandPrinceCoefficients(
                        rungeKutta87DormandPrinceCoefficients );
        }
        return rungeKutta87DormandPrinceCoefficients;

    case rungeKuttaFehlberg89:
        if ( rungeKuttaFehlberg89Coefficients.higherOrder != 9 )
        {
            initializeRungeKuttaFehlberg89Coefficients( rungeKuttaFehlberg89Coefficients );
        }
        return rungeKuttaFehlberg89Coefficients;

    case rungeKuttaVerner89:
        if ( rungeKuttaVerner89Coefficients.higherOrder != 9 )
        {
            initializeRungeKuttaVerner89Coefficients( rungeKuttaVerner89Coefficients );
        }
        return rungeKuttaVerner89Coefficients;
    
    case rungeKuttaFeagin108:
        if ( rungeKuttaFeagin108Coefficients.higherOrder != 10 )
        {
            initializeRungeKuttaFeagin108Coefficients( rungeKuttaFeagin108Coefficients );
        }
        return rungeKuttaFeagin108Coefficients;

    case rungeKuttaFeagin1210:
        if ( rungeKuttaFeagin1210Coefficients.higherOrder != 12 )
        {
            initializeRungeKuttaFeagin1210Coefficients( rungeKuttaFeagin1210Coefficients );
        }
        return rungeKuttaFeagin1210Coefficients;

    case rungeKuttaFeagin1412:
        if ( rungeKuttaFeagin1412Coefficients.higherOrder != 14 )
        {
            initializeRungeKuttaFeagin1412Coefficients( rungeKuttaFeagin1412Coefficients );
        }
        return rungeKuttaFeagin1412Coefficients;

    default: // The default case will never occur because CoefficientsSet is an enum.
        throw RungeKuttaCoefficients( );
    }
}

// Function to print the Butcher tableau of a given coefficient set.
void printButcherTableau( CoefficientSets coefficientSet )
{
    RungeKuttaCoefficients coefficients;
    coefficients = coefficients.get( coefficientSet );

    std::cout << "Butcher tableau of the " << coefficients.name << " coefficients: " << std::endl;

    // Create a zero matrix of the same size as aCoefficients plus 1.
    Eigen::MatrixXd ButcherTable = Eigen::MatrixXd::Zero(
        coefficients.aCoefficients.rows( ) + coefficients.bCoefficients.rows( ),
        coefficients.bCoefficients.cols( ) + 1 );
    
    // Set the table first column to the values of cCoefficients.
    ButcherTable.block( 0, 0, coefficients.cCoefficients.rows( ), 1 ) = coefficients.cCoefficients;

    // Set the table last row(s) to the values of bCoefficients.
    ButcherTable.block(
        coefficients.aCoefficients.rows( ), 
        1,
        coefficients.bCoefficients.rows( ),
        coefficients.bCoefficients.cols( )
    ) = coefficients.bCoefficients;

    // Set the rest of the table to the values of aCoefficients.
    ButcherTable.block(
        0,
        1,
        coefficients.aCoefficients.rows( ),
        coefficients.aCoefficients.cols( )
    ) = coefficients.aCoefficients;
    
    // Feed the full Butcher tableau into a stringstream (to make use of the precision/formatting from IOFormat implemented in Eigen).
    std::stringstream tableStream;
    tableStream << ButcherTable;
    int line_i = 0;
    int first_column_end = -1;
    std::string line;
    // Go trough each of the table lines.
    while(std::getline(tableStream,line,'\n'))
    {
        // If the index at which the first column ends is not known yet, find it.
        if (first_column_end == -1)
        {
            bool has_encountered_non_space = false;
            // Go trough each of the line characters.
            for (unsigned int i = 0; i < line.length( ); i++)
            {
                // If the character is not a space, remember it.
                if (line[i] != ' ')
                {
                    has_encountered_non_space = true;
                }
                // If the character is a space and we have encountered a non-space character before...
                if (line[i] == ' ' && has_encountered_non_space)
                {
                    // Remember the index at which the first column ends.
                    first_column_end = i;
                    break;
                }
            }
        }

        // Add a vertical bar after the end of the first column.
        line.insert(first_column_end + 1, "| ");

        // If we are in the last row(s)...
        if (line_i >= ButcherTable.rows( ) - coefficients.bCoefficients.rows( ))
        {
            // Replace every character in the line before the first column end by a space (do not print 0 in bottom left tableau section).
            for (int i = 0; i < first_column_end; i++)
            {
                line[i] = ' ';
            }                
        }

        // Print the line.
        std::cout << line << std::endl;

        // Print a dash line to show the separation between a and b coefficients.
        if(line_i == ButcherTable.rows( ) - coefficients.bCoefficients.rows( ) - 1)
        {
            std::string dash_line(line.length( ) - 1, '-');
            dash_line.insert(first_column_end + 1, "|");
            std::cout << dash_line << std::endl;
        }
        line_i++;
    }
}

} // namespace numerical_integrators
} // namespace tudat
