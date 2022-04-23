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
    heunEulerCoefficients.bCoefficients( 0, 0 ) = 1.0 / 2.0;
    heunEulerCoefficients.bCoefficients( 1, 0 ) = 1.0;
    heunEulerCoefficients.bCoefficients( 0, 1 ) = 1.0 / 2.0;

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
    rungeKuttaFehlberg12Coefficients.bCoefficients( 0, 0 ) = 1.0 / 512.0;
    rungeKuttaFehlberg12Coefficients.bCoefficients( 0, 1 ) = 255.0 / 256.0;
    rungeKuttaFehlberg12Coefficients.bCoefficients( 0, 2 ) = 1.0 / 512.0;
    rungeKuttaFehlberg12Coefficients.bCoefficients( 1, 0 ) = 1.0 / 256.0;
    rungeKuttaFehlberg12Coefficients.bCoefficients( 1, 1 ) = 255.0 / 256.0;

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
    rungeKutta4Coefficients.bCoefficients( 0, 1 ) = 1.0 / 3.0;
    rungeKutta4Coefficients.bCoefficients( 0, 2 ) = 1.0 / 3.0;
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
    rungeKuttaFehlberg89Coefficients.bCoefficients( 1, 14 ) = 0.062982578976076890319981206482371;
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
    rungeKutta1412Coefficients.aCoefficients( 1, 0 ) = 0.111111111111111111111111111111111111111111111111111111111111;
    rungeKutta1412Coefficients.aCoefficients( 2, 0 ) = -0.833333333333333333333333333333333333333333333333333333333333;
    rungeKutta1412Coefficients.aCoefficients( 2, 1 ) = 1.38888888888888888888888888888888888888888888888888888888889;

    rungeKutta1412Coefficients.aCoefficients( 3, 0 ) = 0.208333333333333333333333333333333333333333333333333333333333;
    rungeKutta1412Coefficients.aCoefficients( 3, 2 ) = 0.625000000000000000000000000000000000000000000000000000000000;

    rungeKutta1412Coefficients.aCoefficients( 4, 0 ) = 0.193333333333333333333333333333333333333333333333333333333333;
    rungeKutta1412Coefficients.aCoefficients( 4, 2 ) = 0.220000000000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.aCoefficients( 4, 3 ) = -0.0800000000000000000000000000000000000000000000000000000000000;

    rungeKutta1412Coefficients.aCoefficients( 5, 0 ) = 0.100000000000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.aCoefficients( 5, 3 ) = 0.400000000000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.aCoefficients( 5, 4 ) = 0.500000000000000000000000000000000000000000000000000000000000;

    rungeKutta1412Coefficients.aCoefficients( 6, 0 ) = 0.103484561636679776672993546511910344499744798201971316606663;
    rungeKutta1412Coefficients.aCoefficients( 6, 3 ) = 0.122068887306407222589644082868962077139592714834162134741275;
    rungeKutta1412Coefficients.aCoefficients( 6, 4 ) = 0.482574490331246622475134780125688112865919023850168049679402;
    rungeKutta1412Coefficients.aCoefficients( 6, 5 ) = -0.0381409600015606999730886240005620205664113072478411477421970;

    rungeKutta1412Coefficients.aCoefficients( 7, 0 ) = 0.124380526654094412881516420868799316268491466359671423163289;
    rungeKutta1412Coefficients.aCoefficients( 7, 4 ) = 0.226120282197584301422238662979202901196752320742633143965145;
    rungeKutta1412Coefficients.aCoefficients( 7, 5 ) = 0.0137885887618080880607695837016477814530969417491493385363543;
    rungeKutta1412Coefficients.aCoefficients( 7, 6 ) = -0.0672210133996684449749399507414305856950086341525382182856200;

    rungeKutta1412Coefficients.aCoefficients( 8, 0 ) = 0.0936919065659673815530885456083005933866349695217750085655603;
    rungeKutta1412Coefficients.aCoefficients( 8, 5 ) = -0.00613406843450510987229498995641664735620914507128858871007099;
    rungeKutta1412Coefficients.aCoefficients( 8, 6 ) = 0.216019825625503063708860097659866573490979433278117320188668;
    rungeKutta1412Coefficients.aCoefficients( 8, 7 ) = 0.423695063515761937337619073960976753205867469544123532683116;

    rungeKutta1412Coefficients.aCoefficients( 9, 0 ) = 0.0838479812409052664616968791372814085980533139224911131069335;
    rungeKutta1412Coefficients.aCoefficients( 9, 5 ) = -0.0117949367100973814319755056031295775367961960590736150777613;
    rungeKutta1412Coefficients.aCoefficients( 9, 6 ) = -0.247299020568812652339473838743194598325992840353340132697498;
    rungeKutta1412Coefficients.aCoefficients( 9, 7 ) = 0.0978080858367729012259313014081291665503740655476733940756599;
    rungeKutta1412Coefficients.aCoefficients( 9, 8 ) = 0.217590689243420631360008651767860318344168120024782176879989;

    rungeKutta1412Coefficients.aCoefficients( 10, 0 ) = 0.0615255359769428227954562389614314714333423969064821107453940;
    rungeKutta1412Coefficients.aCoefficients( 10, 5 ) = 0.00592232780324503308042990005798046524738389560444257136834990;
    rungeKutta1412Coefficients.aCoefficients( 10, 6 ) = 0.470326159963841112217224303205894113455362530746108825010848;
    rungeKutta1412Coefficients.aCoefficients( 10, 7 ) = 0.299688863848679000853981837096192399136831121671781279184194;
    rungeKutta1412Coefficients.aCoefficients( 10, 8 ) = -0.247656877593994914689992276329810825853958069263947095548189;
    rungeKutta1412Coefficients.aCoefficients( 10, 9 ) = 0.110895029771437682893999851839061714522445173600678718208625;

    rungeKutta1412Coefficients.aCoefficients( 11, 0 ) = 0.0419700073362782579861792864787277787213483656543104611245994;
    rungeKutta1412Coefficients.aCoefficients( 11, 5 ) = -0.00317987696266205093901912847692712407988609169703103952205634;
    rungeKutta1412Coefficients.aCoefficients( 11, 6 ) = 0.806397714906192077260821711520379506393543111567419750119748;
    rungeKutta1412Coefficients.aCoefficients( 11, 7 ) = 0.0975983126412388979093522850684288851314672048003054550357187;
    rungeKutta1412Coefficients.aCoefficients( 11, 8 ) = 0.778575578158398909027512446452927238999763460594181964958853;
    rungeKutta1412Coefficients.aCoefficients( 11, 9 ) = 0.204890423831599428189499202098105603312029235081420653574829;
    rungeKutta1412Coefficients.aCoefficients( 11, 10 ) = -1.56261579627468188307070943950527825211462892236424360892806;

    rungeKutta1412Coefficients.aCoefficients( 12, 0 ) = 0.0437726782233730163574465242495339811688214967071614123256973;
    rungeKutta1412Coefficients.aCoefficients( 12, 8 ) = 0.00624365027520195208794358628580933625281631216903095917201250;
    rungeKutta1412Coefficients.aCoefficients( 12, 9 ) = 0.200043097109577314994435165469647856829066232218264969608768;
    rungeKutta1412Coefficients.aCoefficients( 12, 10 ) = -0.00805328367804983036823857162048902911923392887337029314844206;
    rungeKutta1412Coefficients.aCoefficients( 12, 11 ) = 0.0211517528067396521915711903523399601316877825157550573051221;

    rungeKutta1412Coefficients.aCoefficients( 13, 0 ) = 0.0283499250363514563095023591920717312247137654896477097768495;
    rungeKutta1412Coefficients.aCoefficients( 13, 8 ) = 0.00249163204855817407538949148805995149459884653585417680098222;
    rungeKutta1412Coefficients.aCoefficients( 13, 9 ) = 0.0230138787854593149638399846373742768772087122638142234223658;
    rungeKutta1412Coefficients.aCoefficients( 13, 10 ) = -0.00322155956692977098724476092467120878189463604760620461043308;
    rungeKutta1412Coefficients.aCoefficients( 13, 11 ) = 0.00988442549447664668946335414487885256040819982786014648129297;
    rungeKutta1412Coefficients.aCoefficients( 13, 12 ) = -0.0213010771328887351384307642875927384886634565429572466632092;

    rungeKutta1412Coefficients.aCoefficients( 14, 0 ) = 0.343511894290243001049432234735147943083353174980701426268122;
    rungeKutta1412Coefficients.aCoefficients( 14, 8 ) = 0.210451912023627385609097011999010655788807405225626700040882;
    rungeKutta1412Coefficients.aCoefficients( 14, 9 ) = 1.03427452057230411936482926828825709938667999698324740166559;
    rungeKutta1412Coefficients.aCoefficients( 14, 10 ) = 0.00600303645864422487051240448206640574939078092406156945568306;
    rungeKutta1412Coefficients.aCoefficients( 14, 11 ) = 0.855938125099619537578012106002407728915062652616416005816477;
    rungeKutta1412Coefficients.aCoefficients( 14, 12 ) = -0.977235005036766810872264852372525633013107656892839677696022;
    rungeKutta1412Coefficients.aCoefficients( 14, 13 ) = -0.660026980479294694616225013856327693720573981219974874776419;

    rungeKutta1412Coefficients.aCoefficients( 15, 0 ) = -0.0143574001672168069538206399935076366657755954378399880691949;
    rungeKutta1412Coefficients.aCoefficients( 15, 8 ) = -0.0366253270049039970293685796848974791733119081733552207318285;
    rungeKutta1412Coefficients.aCoefficients( 15, 9 ) = 0.0350254975636213681976849406979846524346789082471103574920148;
    rungeKutta1412Coefficients.aCoefficients( 15, 10 ) = 0.0360946016362113508931786658758335239823689929864237671348749;
    rungeKutta1412Coefficients.aCoefficients( 15, 11 ) = -0.0265219967553681106351595946834601923649627012457464284442911;
    rungeKutta1412Coefficients.aCoefficients( 15, 12 ) = 0.0445699011305698119638911537508839908104336323082226770910408;
    rungeKutta1412Coefficients.aCoefficients( 15, 13 ) = 0.124343093331358243286225595741786448038973408895106741855721;
    rungeKutta1412Coefficients.aCoefficients( 15, 14 ) = 0.00413829693239480694403512496204335960426192908674476033832967;

    rungeKutta1412Coefficients.aCoefficients( 16, 0 ) = 0.356032404425120290975609116398089176264106222379748802654822;
    rungeKutta1412Coefficients.aCoefficients( 16, 8 ) = -0.450192758947562595966821779075956175110645100214763601190349;
    rungeKutta1412Coefficients.aCoefficients( 16, 9 ) = 0.430527907083710898626656292808782917793030154094709462877146;
    rungeKutta1412Coefficients.aCoefficients( 16, 10 ) = 0.511973029011022237668556960394071692077125787030651386389972;
    rungeKutta1412Coefficients.aCoefficients( 16, 11 ) = 0.908303638886404260390159124638110213997496214819904630546596;
    rungeKutta1412Coefficients.aCoefficients( 16, 12 ) = -1.23921093371933931757372469151534028854413889248605726186520;
    rungeKutta1412Coefficients.aCoefficients( 16, 13 ) = -0.649048661671761465141672348879062553905402831967191097656668;
    rungeKutta1412Coefficients.aCoefficients( 16, 14 ) = 0.251708904586819292210480529948970541404887852931447491219418;
    rungeKutta1412Coefficients.aCoefficients( 16, 15 ) = 0.779906470345586398810756795282334476023540593411550187024263;

    rungeKutta1412Coefficients.aCoefficients( 17, 0 ) = 0.0130935687406513066406881206418834980127470438213192487844956;
    rungeKutta1412Coefficients.aCoefficients( 17, 12 ) = -0.0000932053067985113945908461962767108237858631509684667142124826;
    rungeKutta1412Coefficients.aCoefficients( 17, 13 ) = 0.0505374334262299359640090443138590726770942344716122381702746;
    rungeKutta1412Coefficients.aCoefficients( 17, 15 ) = 0.000591726029494171190528755742777717259844340971924321528178248;

    rungeKutta1412Coefficients.aCoefficients( 18, 0 ) = 0.0207926484466053012541944544000765652167255206144373407979758;
    rungeKutta1412Coefficients.aCoefficients( 18, 12 ) = 0.000582695918800085915101902697837284108951406103029871570103075;
    rungeKutta1412Coefficients.aCoefficients( 18, 13 ) = -0.00801700732358815939083342186525852746640558465919633524655451;
    rungeKutta1412Coefficients.aCoefficients( 18, 15 ) = 0.0854609998055506144225056114567535602510114622033622491802597;
    rungeKutta1412Coefficients.aCoefficients( 18, 17 ) = 0.105328578824431893399799402979093997354240904235172843146582;

    rungeKutta1412Coefficients.aCoefficients( 19, 0 ) = 1.40153449795736021415446247355771306718486452917597731683689;
    rungeKutta1412Coefficients.aCoefficients( 19, 12 ) = -0.230252000984221261616272410367415621261130298274455611733277;
    rungeKutta1412Coefficients.aCoefficients( 19, 13 ) = -7.21106840466912905659582237106874247165856493509961561958267;
    rungeKutta1412Coefficients.aCoefficients( 19, 14 ) = 0.00372901560694836335236995327852132340217759566678662385552634;
    rungeKutta1412Coefficients.aCoefficients( 19, 15 ) = -4.71415495727125020678778179392224757011323373221820091641216;
    rungeKutta1412Coefficients.aCoefficients( 19, 16 ) = -0.00176367657545349242053841995032797673574903886695600132759652;
    rungeKutta1412Coefficients.aCoefficients( 19, 17 ) = 7.64130548038698765563029310880237651185173367813936997648198;
    rungeKutta1412Coefficients.aCoefficients( 19, 18 ) = 3.50602043659751834989896082949744710968212949893375368243588;

    rungeKutta1412Coefficients.aCoefficients( 20, 0 ) = 11.9514650694120686799372385830716401674473610826553517297976;
    rungeKutta1412Coefficients.aCoefficients( 20, 12 ) = 7.79480932108175968783516700231764388220284279598980948538579;
    rungeKutta1412Coefficients.aCoefficients( 20, 13 ) = -56.4501393867325792523560991120904281440468100061340556540132;
    rungeKutta1412Coefficients.aCoefficients( 20, 14 ) = 0.0912376306930644901344530449290276645709607450403673704844997;
    rungeKutta1412Coefficients.aCoefficients( 20, 15 ) = -12.7336279925434886201945524309199275038162717529918963305155;
    rungeKutta1412Coefficients.aCoefficients( 20, 16 ) = -0.0396895921904719712313542810939736674712383070433147873009352;
    rungeKutta1412Coefficients.aCoefficients( 20, 17 ) = 54.4392141883570886996225765155307791861438378423305337073797;
    rungeKutta1412Coefficients.aCoefficients( 20, 18 ) = -3.64411637921569236846406990361350645806721478409266709351203;
    rungeKutta1412Coefficients.aCoefficients( 20, 19 ) = -0.804503249910509910899030787958579499315694913210787878260459;

    rungeKutta1412Coefficients.aCoefficients( 21, 0 ) = -148.809426507100488427838868268647625561930612082148597076690;
    rungeKutta1412Coefficients.aCoefficients( 21, 12 ) = -91.7295278291256484357935662402321623495228729036354276506427;
    rungeKutta1412Coefficients.aCoefficients( 21, 13 ) = 707.656144971598359834575719286335716154821128966649565194286;
    rungeKutta1412Coefficients.aCoefficients( 21, 14 ) = -1.10563611857482440905296961311590930801338308942637769555540;
    rungeKutta1412Coefficients.aCoefficients( 21, 15 ) = 176.134591883811372587859898076055660406999516762301689616841;
    rungeKutta1412Coefficients.aCoefficients( 21, 16 ) = 0.491384824214880662268898345164454557416884631402764792538746;
    rungeKutta1412Coefficients.aCoefficients( 21, 17 ) = -684.278000449814944358237535610895081956077167893600278300805;
    rungeKutta1412Coefficients.aCoefficients( 21, 18 ) = 27.9910604998398258984224332124380407446002518400668657974589;
    rungeKutta1412Coefficients.aCoefficients( 21, 19 ) = 13.1939710030282333443670964371153238435064159623744975073252;
    rungeKutta1412Coefficients.aCoefficients( 21, 20 ) = 1.25128781283980445450114974148056006317268830077396406361417;

    rungeKutta1412Coefficients.aCoefficients( 22, 0 ) = -9.67307946948196763644126118433219395839951408571877262880482;
    rungeKutta1412Coefficients.aCoefficients( 22, 12 ) = -4.46990150858505531443846227701960360497830681408751431146712;
    rungeKutta1412Coefficients.aCoefficients( 22, 13 ) = 45.5127128690952681968241950400052751178905907817398483534845;
    rungeKutta1412Coefficients.aCoefficients( 22, 14 ) = -0.0713085086183826912791492024438246129930559805352394367050813;
    rungeKutta1412Coefficients.aCoefficients( 22, 15 ) = 11.2273614068412741582590624479939384207826800776794485051540;
    rungeKutta1412Coefficients.aCoefficients( 22, 16 ) = 0.126244376717622724516237912909138809361786889819105426371393;
    rungeKutta1412Coefficients.aCoefficients( 22, 17 ) = -43.5439339549483313605810624907242107623814304467621407753424;
    rungeKutta1412Coefficients.aCoefficients( 22, 18 ) = 0.787174307543058978398792994996550902064546091443233850464377;
    rungeKutta1412Coefficients.aCoefficients( 22, 19 ) = 0.532264696744684215669300708603886690785395776821503851830821;
    rungeKutta1412Coefficients.aCoefficients( 22, 20 ) = 0.422422733996325326010225127471388772575086538809603346825334;
    rungeKutta1412Coefficients.aCoefficients( 22, 21 ) = 0.0859131249503067107308438031499859443441115056294154956487671;

    rungeKutta1412Coefficients.aCoefficients( 23, 0 ) = -10.0664032447054702403396606900426891472202824757968765569183;
    rungeKutta1412Coefficients.aCoefficients( 23, 8 ) = -0.0366253270049039970293685796848974791733119081733552207318285;
    rungeKutta1412Coefficients.aCoefficients( 23, 9 ) = 0.0350254975636213681976849406979846524346789082471103574920148;
    rungeKutta1412Coefficients.aCoefficients( 23, 10 ) = 0.0360946016362113508931786658758335239823689929864237671348749;
    rungeKutta1412Coefficients.aCoefficients( 23, 11 ) = -0.0265219967553681106351595946834601923649627012457464284442911;
    rungeKutta1412Coefficients.aCoefficients( 23, 12 ) = -6.27088972181464143590553149478871603839356122957396018530209;
    rungeKutta1412Coefficients.aCoefficients( 23, 13 ) = 48.2079237442562989090702103008195063923492593141636117832993;
    rungeKutta1412Coefficients.aCoefficients( 23, 14 ) = -0.0694471689136165640882395180583732834557754169149088630301342;
    rungeKutta1412Coefficients.aCoefficients( 23, 15 ) = 12.6810690204850295698341370913609807066108483811412127009785;
    rungeKutta1412Coefficients.aCoefficients( 23, 16 ) = 0.0119671168968323754838161435501011294100927813964199613229864;
    rungeKutta1412Coefficients.aCoefficients( 23, 17 ) = -46.7249764992482408003358268242662695593201321659795608950429;
    rungeKutta1412Coefficients.aCoefficients( 23, 18 ) = 1.33029613326626711314710039298216591399033511191227101321435;
    rungeKutta1412Coefficients.aCoefficients( 23, 19 ) = 1.00766787503398298353438903619926657771162717793661719708370;
    rungeKutta1412Coefficients.aCoefficients( 23, 20 ) = 0.0209512051933665091664122388475480702892770753864487241177616;
    rungeKutta1412Coefficients.aCoefficients( 23, 21 ) = 0.0210134706331264177317735424331396407424412188443757490871603;
    rungeKutta1412Coefficients.aCoefficients( 23, 22 ) = 0.00952196014417121794175101542454575907376360233658356240547761;

    rungeKutta1412Coefficients.aCoefficients( 24, 0 ) = -409.478081677743708772589097409370357624424341606752069725341;
    rungeKutta1412Coefficients.aCoefficients( 24, 8 ) = 0.210451912023627385609097011999010655788807405225626700040882;
    rungeKutta1412Coefficients.aCoefficients( 24, 9 ) = 1.03427452057230411936482926828825709938667999698324740166559;
    rungeKutta1412Coefficients.aCoefficients( 24, 10 ) = 0.00600303645864422487051240448206640574939078092406156945568306;
    rungeKutta1412Coefficients.aCoefficients( 24, 11 ) = 0.855938125099619537578012106002407728915062652616416005816477;
    rungeKutta1412Coefficients.aCoefficients( 24, 12 ) = -250.516998547447860492777657729316130386584050420782075966990;
    rungeKutta1412Coefficients.aCoefficients( 24, 13 ) = 1946.42466652388427766053750328264758595829850895761428240231;
    rungeKutta1412Coefficients.aCoefficients( 24, 14 ) = -3.04503882102310365506105809086860882786950544097602101685174;
    rungeKutta1412Coefficients.aCoefficients( 24, 15 ) = 490.626379528281713521208265299168083841598542274061671576230;
    rungeKutta1412Coefficients.aCoefficients( 24, 16 ) = 1.56647589531270907115484067013597445739595615245966775329993;
    rungeKutta1412Coefficients.aCoefficients( 24, 17 ) = -1881.97428994011173362217267377035870619215906638453056643641;
    rungeKutta1412Coefficients.aCoefficients( 24, 18 ) = 75.2592224724847175278837713643303149821620618914245864351135;
    rungeKutta1412Coefficients.aCoefficients( 24, 19 ) = 34.5734356980331067622434344736554689696728644793551014989002;
    rungeKutta1412Coefficients.aCoefficients( 24, 20 ) = 3.21147679440968961435417361847073755169022966748891627882572;
    rungeKutta1412Coefficients.aCoefficients( 24, 21 ) = -0.460408041738414391307201404237058848867245095265382820823055;
    rungeKutta1412Coefficients.aCoefficients( 24, 22 ) = -0.0870718339841810522431884137957986245724252047388936572215438;
    rungeKutta1412Coefficients.aCoefficients( 24, 23 ) = -7.39351814158303067567016952195521063999185773249132944724553;

    rungeKutta1412Coefficients.aCoefficients( 25, 0 ) = 3.43347475853550878921093496257596781120623891072008459930197;
    rungeKutta1412Coefficients.aCoefficients( 25, 8 ) = 0.00249163204855817407538949148805995149459884653585417680098222;
    rungeKutta1412Coefficients.aCoefficients( 25, 9 ) = 0.0230138787854593149638399846373742768772087122638142234223658;
    rungeKutta1412Coefficients.aCoefficients( 25, 10 ) = -0.00322155956692977098724476092467120878189463604760620461043308;
    rungeKutta1412Coefficients.aCoefficients( 25, 11 ) = 0.00988442549447664668946335414487885256040819982786014648129297;
    rungeKutta1412Coefficients.aCoefficients( 25, 12 ) = 2.16252799377922507788307841904757354045759225335732707916530;
    rungeKutta1412Coefficients.aCoefficients( 25, 13 ) = -16.2699864546457421328065640660139489006987552040228852402716;
    rungeKutta1412Coefficients.aCoefficients( 25, 14 ) = -0.128534502120524552843583417470935010538029037542654506231743;
    rungeKutta1412Coefficients.aCoefficients( 25, 15 ) = -8.98915042666504253089307820833379330486511746063552853023189;
    rungeKutta1412Coefficients.aCoefficients( 25, 16 ) = -0.00348595363232025333387080201851013650192401767250513765000963;
    rungeKutta1412Coefficients.aCoefficients( 25, 17 ) = 15.7936194113339807536235187388695574135853387025139738341334;
    rungeKutta1412Coefficients.aCoefficients( 25, 18 ) = -0.574403330914095065628165482017335820148383663195675408024658;
    rungeKutta1412Coefficients.aCoefficients( 25, 19 ) = -0.345602039021393296692722496608124982535237228827655306030152;
    rungeKutta1412Coefficients.aCoefficients( 25, 20 ) = -0.00662241490206585091731619991383757781133067992707418687587487;
    rungeKutta1412Coefficients.aCoefficients( 25, 21 ) = -0.00777788129242204164032546458607364309759347209626759111946150;
    rungeKutta1412Coefficients.aCoefficients( 25, 22 ) = -0.00356084192402274913338827232697437364675240818791706587952939;
    rungeKutta1412Coefficients.aCoefficients( 25, 23 ) = 4.79282506449930799649797749629840189457296934139359048988332;
    rungeKutta1412Coefficients.aCoefficients( 25, 24 ) = 0.153725464873068577844576387402512082757034273069877432944621;

    rungeKutta1412Coefficients.aCoefficients( 26, 0 ) = 32.3038520871985442326994734440031535091364975047784630088983;
    rungeKutta1412Coefficients.aCoefficients( 26, 5 ) = -0.00317987696266205093901912847692712407988609169703103952205634;
    rungeKutta1412Coefficients.aCoefficients( 26, 6 ) = 0.806397714906192077260821711520379506393543111567419750119748;
    rungeKutta1412Coefficients.aCoefficients( 26, 7 ) = 0.0975983126412388979093522850684288851314672048003054550357187;
    rungeKutta1412Coefficients.aCoefficients( 26, 8 ) = 0.778575578158398909027512446452927238999763460594181964958853;
    rungeKutta1412Coefficients.aCoefficients( 26, 9 ) = 0.204890423831599428189499202098105603312029235081420653574829;
    rungeKutta1412Coefficients.aCoefficients( 26, 10 ) = -1.56261579627468188307070943950527825211462892236424360892806;
    rungeKutta1412Coefficients.aCoefficients( 26, 12 ) = 16.3429891882310570648504243973927174708753353504154550405647;
    rungeKutta1412Coefficients.aCoefficients( 26, 13 ) = -154.544555293543621230730189631471036399316683669609116705323;
    rungeKutta1412Coefficients.aCoefficients( 26, 14 ) = 1.56971088703334872692034283417621761466263593582497085955201;
    rungeKutta1412Coefficients.aCoefficients( 26, 15 ) = 3.27685545087248131321429817269900731165522404974733504794135;
    rungeKutta1412Coefficients.aCoefficients( 26, 16 ) = -0.0503489245193653176348040727199783626534081095691632396802451;
    rungeKutta1412Coefficients.aCoefficients( 26, 17 ) = 153.321151858041665070593767885914694011224363102594556731397;
    rungeKutta1412Coefficients.aCoefficients( 26, 18 ) = 7.17568186327720495846766484814784143567826308034865369443637;
    rungeKutta1412Coefficients.aCoefficients( 26, 19 ) = -2.94036748675300481945917659896930989215320594380777597403592;
    rungeKutta1412Coefficients.aCoefficients( 26, 20 ) = -0.0665845946076803144470749676022628870281920493197256887985612;
    rungeKutta1412Coefficients.aCoefficients( 26, 21 ) = -0.0462346054990843661229248668562217261176966514016859284197145;
    rungeKutta1412Coefficients.aCoefficients( 26, 22 ) = -0.0204198733585679401539388228617269778848579774821581777675337;
    rungeKutta1412Coefficients.aCoefficients( 26, 23 ) = -53.3523106438735850515953441165998107974045090495791591218714;
    rungeKutta1412Coefficients.aCoefficients( 26, 24 ) = -1.35548714715078654978732186705996404017554501614191325114947;
    rungeKutta1412Coefficients.aCoefficients( 26, 25 ) = -1.57196275801232751882901735171459249177687219114442583461866;

    rungeKutta1412Coefficients.aCoefficients( 27, 0 ) = -16.6451467486341512872031294403931758764560371130818978459405;
    rungeKutta1412Coefficients.aCoefficients( 27, 5 ) = 0.00592232780324503308042990005798046524738389560444257136834990;
    rungeKutta1412Coefficients.aCoefficients( 27, 6 ) = 0.470326159963841112217224303205894113455362530746108825010848;
    rungeKutta1412Coefficients.aCoefficients( 27, 7 ) = 0.299688863848679000853981837096192399136831121671781279184194;
    rungeKutta1412Coefficients.aCoefficients( 27, 8 ) = -0.247656877593994914689992276329810825853958069263947095548189;
    rungeKutta1412Coefficients.aCoefficients( 27, 9 ) = 0.110895029771437682893999851839061714522445173600678718208625;
    rungeKutta1412Coefficients.aCoefficients( 27, 11 ) = -0.491719043846229147070666628704194097678081907210673044988866;
    rungeKutta1412Coefficients.aCoefficients( 27, 12 ) = -11.4743154427289496968389492564352536350842454130853175250727;
    rungeKutta1412Coefficients.aCoefficients( 27, 13 ) = 80.2593166576230272541702485886484400152793366623589989106256;
    rungeKutta1412Coefficients.aCoefficients( 27, 14 ) = -0.384132303980042847625312526759029103746926841342088219165648;
    rungeKutta1412Coefficients.aCoefficients( 27, 15 ) = 7.28147667468107583471326950926136115767612581862877764249646;
    rungeKutta1412Coefficients.aCoefficients( 27, 16 ) = -0.132699384612248379510571708176035274836827341616751884314074;
    rungeKutta1412Coefficients.aCoefficients( 27, 17 ) = -81.0799832525730726674679289752255240006070716633632990308935;
    rungeKutta1412Coefficients.aCoefficients( 27, 18 ) = -1.25037492835620639521768185656179119962253747492403205797494;
    rungeKutta1412Coefficients.aCoefficients( 27, 19 ) = 2.59263594969543681023776379504377324994226447359296887778718;
    rungeKutta1412Coefficients.aCoefficients( 27, 20 ) = -0.301440298346404539830163997260526875264431537275641495291993;
    rungeKutta1412Coefficients.aCoefficients( 27, 21 ) = 0.221384460789832337451706451572773791695246839057318414301020;
    rungeKutta1412Coefficients.aCoefficients( 27, 22 ) = 0.0827577274771892931955989870974693152996276435429809890551210;
    rungeKutta1412Coefficients.aCoefficients( 27, 23 ) = 18.9960662040611520464672450037243263998175161412237156872211;
    rungeKutta1412Coefficients.aCoefficients( 27, 24 ) = 0.269231946409639685623468015128334167460051910348912845121977;
    rungeKutta1412Coefficients.aCoefficients( 27, 25 ) = 1.62674827447066537462989364929628933988125029284183680279020;
    rungeKutta1412Coefficients.aCoefficients( 27, 26 ) = 0.491719043846229147070666628704194097678081907210673044988866;

    rungeKutta1412Coefficients.aCoefficients( 28, 0 ) = 0.0838479812409052664616968791372814085980533139224911131069335;
    rungeKutta1412Coefficients.aCoefficients( 28, 5 ) = -0.0117949367100973814319755056031295775367961960590736150777613;
    rungeKutta1412Coefficients.aCoefficients( 28, 6 ) = -0.247299020568812652339473838743194598325992840353340132697498;
    rungeKutta1412Coefficients.aCoefficients( 28, 7 ) = 0.0978080858367729012259313014081291665503740655476733940756599;
    rungeKutta1412Coefficients.aCoefficients( 28, 8 ) = 0.217590689243420631360008651767860318344168120024782176879989;
    rungeKutta1412Coefficients.aCoefficients( 28, 10 ) = 0.137585606763325224865659632196787746647447222975084865975440;
    rungeKutta1412Coefficients.aCoefficients( 28, 11 ) = 0.0439870229715046685058790092341545026046103890294261359042581;
    rungeKutta1412Coefficients.aCoefficients( 28, 13 ) = -0.513700813768193341957004456618630303738757363641964030086972;
    rungeKutta1412Coefficients.aCoefficients( 28, 14 ) = 0.826355691151315508644211308399153458701423158616168576922372;
    rungeKutta1412Coefficients.aCoefficients( 28, 15 ) = 25.7018139719811832625873882972519939511136556341960074626615;
    rungeKutta1412Coefficients.aCoefficients( 28, 23 ) = -25.7018139719811832625873882972519939511136556341960074626615;
    rungeKutta1412Coefficients.aCoefficients( 28, 24 ) = -0.826355691151315508644211308399153458701423158616168576922372;
    rungeKutta1412Coefficients.aCoefficients( 28, 25 ) = 0.513700813768193341957004456618630303738757363641964030086972;
    rungeKutta1412Coefficients.aCoefficients( 28, 26 ) = -0.0439870229715046685058790092341545026046103890294261359042581;
    rungeKutta1412Coefficients.aCoefficients( 28, 27 ) = -0.137585606763325224865659632196787746647447222975084865975440;

    rungeKutta1412Coefficients.aCoefficients( 29, 0 ) = 0.124380526654094412881516420868799316268491466359671423163289;
    rungeKutta1412Coefficients.aCoefficients( 29, 4 ) = 0.226120282197584301422238662979202901196752320742633143965145;
    rungeKutta1412Coefficients.aCoefficients( 29, 5 ) = 0.0137885887618080880607695837016477814530969417491493385363543;
    rungeKutta1412Coefficients.aCoefficients( 29, 6 ) = -0.0672210133996684449749399507414305856950086341525382182856200;
    rungeKutta1412Coefficients.aCoefficients( 29, 9 ) = -0.856238975085428354755349769879501772112121597411563802855067;
    rungeKutta1412Coefficients.aCoefficients( 29, 10 ) = -1.96337522866858908928262850028093813988180440518267404553576;
    rungeKutta1412Coefficients.aCoefficients( 29, 11 ) = -0.232332822724119401237246257308921847250108199230419994978218;
    rungeKutta1412Coefficients.aCoefficients( 29, 13 ) = 4.30660719086453349461668936876562947772432562053478092626764;
    rungeKutta1412Coefficients.aCoefficients( 29, 14 ) = -2.92722963249465482659787911202390446687687394950633612630592;
    rungeKutta1412Coefficients.aCoefficients( 29, 15 ) = -82.3131666397858944454492334105458707735761966428138676971041;
    rungeKutta1412Coefficients.aCoefficients( 29, 23 ) = 82.3131666397858944454492334105458707735761966428138676971041;
    rungeKutta1412Coefficients.aCoefficients( 29, 24 ) = 2.92722963249465482659787911202390446687687394950633612630592;
    rungeKutta1412Coefficients.aCoefficients( 29, 25 ) = -4.30660719086453349461668936876562947772432562053478092626764;
    rungeKutta1412Coefficients.aCoefficients( 29, 26 ) = 0.232332822724119401237246257308921847250108199230419994978218;
    rungeKutta1412Coefficients.aCoefficients( 29, 27 ) = 1.96337522866858908928262850028093813988180440518267404553576;
    rungeKutta1412Coefficients.aCoefficients( 29, 28 ) = 0.856238975085428354755349769879501772112121597411563802855067;

    rungeKutta1412Coefficients.aCoefficients( 30, 0 ) = 0.103484561636679776672993546511910344499744798201971316606663;
    rungeKutta1412Coefficients.aCoefficients( 30, 3 ) = 0.122068887306407222589644082868962077139592714834162134741275;
    rungeKutta1412Coefficients.aCoefficients( 30, 4 ) = 0.482574490331246622475134780125688112865919023850168049679402;
    rungeKutta1412Coefficients.aCoefficients( 30, 5 ) = -0.0381409600015606999730886240005620205664113072478411477421970;
    rungeKutta1412Coefficients.aCoefficients( 30, 7 ) = -0.550499525310802324138388507020508177411414311000037561712836;
    rungeKutta1412Coefficients.aCoefficients( 30, 9 ) = -0.711915811585189227887648262043794387578291882406745570495765;
    rungeKutta1412Coefficients.aCoefficients( 30, 10 ) = -0.584129605671551340432988730158480872095335329645227595707052;
    rungeKutta1412Coefficients.aCoefficients( 30, 13 ) = 2.11046308125864932128717300046622750300375054278936987850718;
    rungeKutta1412Coefficients.aCoefficients( 30, 14 ) = -0.0837494736739572135525742023001037992695260175335123517729291;
    rungeKutta1412Coefficients.aCoefficients( 30, 15 ) = 5.10021499072320914075295969043344113107545060862804249161191;
    rungeKutta1412Coefficients.aCoefficients( 30, 23 ) = -5.10021499072320914075295969043344113107545060862804249161191;
    rungeKutta1412Coefficients.aCoefficients( 30, 24 ) = 0.0837494736739572135525742023001037992695260175335123517729291;
    rungeKutta1412Coefficients.aCoefficients( 30, 25 ) = -2.11046308125864932128717300046622750300375054278936987850718;
    rungeKutta1412Coefficients.aCoefficients( 30, 27 ) = 0.584129605671551340432988730158480872095335329645227595707052;
    rungeKutta1412Coefficients.aCoefficients( 30, 28 ) = 0.711915811585189227887648262043794387578291882406745570495765;
    rungeKutta1412Coefficients.aCoefficients( 30, 29 ) = 0.550499525310802324138388507020508177411414311000037561712836;

    rungeKutta1412Coefficients.aCoefficients( 31, 0 ) = 0.193333333333333333333333333333333333333333333333333333333333;
    rungeKutta1412Coefficients.aCoefficients( 31, 2 ) = 0.220000000000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.aCoefficients( 31, 3 ) = -0.0800000000000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.aCoefficients( 31, 6 ) = 0.109993425580724703919462404865068340845119058295846426463652;
    rungeKutta1412Coefficients.aCoefficients( 31, 7 ) = -0.254297048076270161384068506997153122141835626976703920846242;
    rungeKutta1412Coefficients.aCoefficients( 31, 9 ) = 0.865570777116694254343770343821098281832847401233011859346737;
    rungeKutta1412Coefficients.aCoefficients( 31, 10 ) = 3.32416449114093083106799552786572018336860092936986407160200;
    rungeKutta1412Coefficients.aCoefficients( 31, 13 ) = -12.0102223315977933882352385148661841260301942633996815127277;
    rungeKutta1412Coefficients.aCoefficients( 31, 14 ) = 0.476601466242493239430442776862061899602963782003580209476163;
    rungeKutta1412Coefficients.aCoefficients( 31, 15 ) = -29.0243011221036390525802623213654099596251221332470910692353;
    rungeKutta1412Coefficients.aCoefficients( 31, 23 ) = 29.0243011221036390525802623213654099596251221332470910692353;
    rungeKutta1412Coefficients.aCoefficients( 31, 24 ) = -0.476601466242493239430442776862061899602963782003580209476163;
    rungeKutta1412Coefficients.aCoefficients( 31, 25 ) = 12.0102223315977933882352385148661841260301942633996815127277;
    rungeKutta1412Coefficients.aCoefficients( 31, 27 ) = -3.32416449114093083106799552786572018336860092936986407160200;
    rungeKutta1412Coefficients.aCoefficients( 31, 28 ) = -0.865570777116694254343770343821098281832847401233011859346737;
    rungeKutta1412Coefficients.aCoefficients( 31, 29 ) = 0.254297048076270161384068506997153122141835626976703920846242;
    rungeKutta1412Coefficients.aCoefficients( 31, 30 ) = -0.109993425580724703919462404865068340845119058295846426463652;

    rungeKutta1412Coefficients.aCoefficients( 32, 0 ) = -0.833333333333333333333333333333333333333333333333333333333333;
    rungeKutta1412Coefficients.aCoefficients( 32, 1 ) = 1.38888888888888888888888888888888888888888888888888888888889;
    rungeKutta1412Coefficients.aCoefficients( 32, 4 ) = -0.750000000000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.aCoefficients( 32, 6 ) = -0.492529543718026304422682049114021320200214681580657784719074;
    rungeKutta1412Coefficients.aCoefficients( 32, 30 ) = 0.492529543718026304422682049114021320200214681580657784719074;
    rungeKutta1412Coefficients.aCoefficients( 32, 31 ) = 0.750000000000000000000000000000000000000000000000000000000000;

    rungeKutta1412Coefficients.aCoefficients( 33, 0 ) = 0.111111111111111111111111111111111111111111111111111111111111;
    rungeKutta1412Coefficients.aCoefficients( 33, 2 ) = -0.222222222222222222222222222222222222222222222222222222222222;
    rungeKutta1412Coefficients.aCoefficients( 33, 32 ) = 0.222222222222222222222222222222222222222222222222222222222222;

    rungeKutta1412Coefficients.aCoefficients( 34, 0 ) = 0.285835140388971558796088842163836414852927537894596466840753;
    rungeKutta1412Coefficients.aCoefficients( 34, 1 ) = 0.291666666666666666666666666666666666666666666666666666666667;
    rungeKutta1412Coefficients.aCoefficients( 34, 2 ) = 0.218750000000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.aCoefficients( 34, 4 ) = 0.164062500000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.aCoefficients( 34, 6 ) = 0.218194354945556658327188241581352107093288824322187941141516;
    rungeKutta1412Coefficients.aCoefficients( 34, 7 ) = 0.180392898478697766863635221946775437719620053641849228562435;
    rungeKutta1412Coefficients.aCoefficients( 34, 9 ) = 0.205713839404845018859120755122929542277570094982808905393991;
    rungeKutta1412Coefficients.aCoefficients( 34, 10 ) = 0.242715791581770239970282927959446515762745971386670541948576;
    rungeKutta1412Coefficients.aCoefficients( 34, 11 ) = 0.246465780813629305833609291181891407799228103869305705137021;
    rungeKutta1412Coefficients.aCoefficients( 34, 12 ) = -3.44991940790890824979834154601622662060370460614931644223924;
    rungeKutta1412Coefficients.aCoefficients( 34, 13 ) = 0.228875562160036081760729060738458584294220372552740218459295;
    rungeKutta1412Coefficients.aCoefficients( 34, 14 ) = 0.283290599702151415321527419056733335978436595493855789831434;
    rungeKutta1412Coefficients.aCoefficients( 34, 15 ) = 3.21085125837766640960131490544236787005557320332238705967955;
    rungeKutta1412Coefficients.aCoefficients( 34, 16 ) = -0.223538777364845699920233756214162507964125230083674032084065;
    rungeKutta1412Coefficients.aCoefficients( 34, 17 ) = -0.707121157204419073518727286207487212130091231955206160635271;
    rungeKutta1412Coefficients.aCoefficients( 34, 18 ) = 3.21123345150287080408174729202856500893260034443022374267639;
    rungeKutta1412Coefficients.aCoefficients( 34, 19 ) = 1.40954348309669766030414474301123175769045945573548986335553;
    rungeKutta1412Coefficients.aCoefficients( 34, 20 ) = -0.151362053443742613121602276742518111090963026203676055891793;
    rungeKutta1412Coefficients.aCoefficients( 34, 21 ) = 0.372350574527014276454724080214619984397121028202148298716575;
    rungeKutta1412Coefficients.aCoefficients( 34, 22 ) = 0.252978746406361336722199907762141285915775728129414319261111;
    rungeKutta1412Coefficients.aCoefficients( 34, 23 ) = -3.21085125837766640960131490544236787005557320332238705967955;
    rungeKutta1412Coefficients.aCoefficients( 34, 24 ) = -0.283290599702151415321527419056733335978436595493855789831434;
    rungeKutta1412Coefficients.aCoefficients( 34, 25 ) = -0.228875562160036081760729060738458584294220372552740218459295;
    rungeKutta1412Coefficients.aCoefficients( 34, 26 ) = -0.246465780813629305833609291181891407799228103869305705137021;
    rungeKutta1412Coefficients.aCoefficients( 34, 27 ) = -0.242715791581770239970282927959446515762745971386670541948576;
    rungeKutta1412Coefficients.aCoefficients( 34, 28 ) = -0.205713839404845018859120755122929542277570094982808905393991;
    rungeKutta1412Coefficients.aCoefficients( 34, 29 ) = -0.180392898478697766863635221946775437719620053641849228562435;
    rungeKutta1412Coefficients.aCoefficients( 34, 30 ) = -0.218194354945556658327188241581352107093288824322187941141516;
    rungeKutta1412Coefficients.aCoefficients( 34, 31 ) = -0.164062500000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.aCoefficients( 34, 32 ) = -0.218750000000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.aCoefficients( 34, 33 ) = -0.291666666666666666666666666666666666666666666666666666666667;
    
    // Define c-coefficients for the Runge-Kutta method.
    rungeKutta1412Coefficients.cCoefficients = Eigen::VectorXd::Zero( 35 );
    rungeKutta1412Coefficients.cCoefficients( 0 ) = 0.000000000000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.cCoefficients( 1 ) = 0.111111111111111111111111111111111111111111111111111111111111;
    rungeKutta1412Coefficients.cCoefficients( 2 ) = 0.555555555555555555555555555555555555555555555555555555555556;
    rungeKutta1412Coefficients.cCoefficients( 3 ) = 0.833333333333333333333333333333333333333333333333333333333333;
    rungeKutta1412Coefficients.cCoefficients( 4 ) = 0.333333333333333333333333333333333333333333333333333333333333;
    rungeKutta1412Coefficients.cCoefficients( 5 ) = 1.00000000000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.cCoefficients( 6 ) = 0.669986979272772921764683785505998513938845229638460353285142;
    rungeKutta1412Coefficients.cCoefficients( 7 ) = 0.297068384213818357389584716808219413223332094698915687379168;
    rungeKutta1412Coefficients.cCoefficients( 8 ) = 0.727272727272727272727272727272727272727272727272727272727273;
    rungeKutta1412Coefficients.cCoefficients( 9 ) = 0.140152799042188765276187487966946717629806463082532936287323;
    rungeKutta1412Coefficients.cCoefficients( 10 ) = 0.700701039770150737151099854830749337941407049265546408969222;
    rungeKutta1412Coefficients.cCoefficients( 11 ) = 0.363636363636363636363636363636363636363636363636363636363636;
    rungeKutta1412Coefficients.cCoefficients( 12 ) = 0.263157894736842105263157894736842105263157894736842105263158;
    rungeKutta1412Coefficients.cCoefficients( 13 ) = 0.0392172246650270859125196642501208648863714315266128052078483;
    rungeKutta1412Coefficients.cCoefficients( 14 ) = 0.812917502928376762983393159278036506189612372617238550774312;
    rungeKutta1412Coefficients.cCoefficients( 15 ) = 0.166666666666666666666666666666666666666666666666666666666667;
    rungeKutta1412Coefficients.cCoefficients( 16 ) = 0.900000000000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.cCoefficients( 17 ) = 0.0641299257451966923312771193896682809481096651615083225402924;
    rungeKutta1412Coefficients.cCoefficients( 18 ) = 0.204149909283428848927744634301023405027149505241333751628870;
    rungeKutta1412Coefficients.cCoefficients( 19 ) = 0.395350391048760565615671369827324372352227297456659450554577;
    rungeKutta1412Coefficients.cCoefficients( 20 ) = 0.604649608951239434384328630172675627647772702543340549445423;
    rungeKutta1412Coefficients.cCoefficients( 21 ) = 0.795850090716571151072255365698976594972850494758666248371130;
    rungeKutta1412Coefficients.cCoefficients( 22 ) = 0.935870074254803307668722880610331719051890334838491677459708;
    rungeKutta1412Coefficients.cCoefficients( 23 ) = 0.166666666666666666666666666666666666666666666666666666666667;
    rungeKutta1412Coefficients.cCoefficients( 24 ) = 0.812917502928376762983393159278036506189612372617238550774312;
    rungeKutta1412Coefficients.cCoefficients( 25 ) = 0.0392172246650270859125196642501208648863714315266128052078483;
    rungeKutta1412Coefficients.cCoefficients( 26 ) = 0.363636363636363636363636363636363636363636363636363636363636;
    rungeKutta1412Coefficients.cCoefficients( 27 ) = 0.700701039770150737151099854830749337941407049265546408969222;
    rungeKutta1412Coefficients.cCoefficients( 28 ) = 0.140152799042188765276187487966946717629806463082532936287323;
    rungeKutta1412Coefficients.cCoefficients( 29 ) = 0.297068384213818357389584716808219413223332094698915687379168;
    rungeKutta1412Coefficients.cCoefficients( 30 ) = 0.669986979272772921764683785505998513938845229638460353285142;
    rungeKutta1412Coefficients.cCoefficients( 31 ) = 0.333333333333333333333333333333333333333333333333333333333333;
    rungeKutta1412Coefficients.cCoefficients( 32 ) = 0.555555555555555555555555555555555555555555555555555555555556;
    rungeKutta1412Coefficients.cCoefficients( 33 ) = 0.111111111111111111111111111111111111111111111111111111111111;
    rungeKutta1412Coefficients.cCoefficients( 34 ) = 1.00000000000000000000000000000000000000000000000000000000000;
    
    // Define b-coefficients for the Runge-Kutta method.
    rungeKutta1412Coefficients.bCoefficients = Eigen::MatrixXd::Zero( 2, 35 );
    rungeKutta1412Coefficients.bCoefficients( 0, 0 ) = 0.0178571428571428571428571428571428571428571428571428571428571;
    rungeKutta1412Coefficients.bCoefficients( 0, 1 ) = 439.0 / 64000.0;
    rungeKutta1412Coefficients.bCoefficients( 0, 2 ) = 0.0117187500000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.bCoefficients( 0, 4 ) = 0.0175781250000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.bCoefficients( 0, 6 ) = 0.0234375000000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.bCoefficients( 0, 7 ) = 0.0292968750000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.bCoefficients( 0, 9 ) = 0.0351562500000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.bCoefficients( 0, 10 ) = 0.0410156250000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.bCoefficients( 0, 11 ) = 0.0468750000000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.bCoefficients( 0, 13 ) = 0.0527343750000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.bCoefficients( 0, 14 ) = 0.0585937500000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.bCoefficients( 0, 15 ) = 0.0644531250000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.bCoefficients( 0, 17 ) = 0.105352113571753019691496032887878162227673083080523884041670;
    rungeKutta1412Coefficients.bCoefficients( 0, 18 ) = 0.170561346241752182382120338553874085887555487802790804737501;
    rungeKutta1412Coefficients.bCoefficients( 0, 19 ) = 0.206229397329351940783526485701104894741914286259542454077972;
    rungeKutta1412Coefficients.bCoefficients( 0, 20 ) = 0.206229397329351940783526485701104894741914286259542454077972;
    rungeKutta1412Coefficients.bCoefficients( 0, 21 ) = 0.170561346241752182382120338553874085887555487802790804737501;
    rungeKutta1412Coefficients.bCoefficients( 0, 22 ) = 0.105352113571753019691496032887878162227673083080523884041670;
    rungeKutta1412Coefficients.bCoefficients( 0, 23 ) = -0.0644531250000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.bCoefficients( 0, 24 ) = -0.0585937500000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.bCoefficients( 0, 25 ) = -0.0527343750000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.bCoefficients( 0, 26 ) = -0.0468750000000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.bCoefficients( 0, 27 ) = -0.0410156250000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.bCoefficients( 0, 28 ) = -0.0351562500000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.bCoefficients( 0, 29 ) = -0.0292968750000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.bCoefficients( 0, 30 ) = -0.0234375000000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.bCoefficients( 0, 31 ) = -0.0175781250000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.bCoefficients( 0, 32 ) = -0.0117187500000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.bCoefficients( 0, 33 ) = -439.0 / 64000.0;
    rungeKutta1412Coefficients.bCoefficients( 0, 34 ) = 0.0178571428571428571428571428571428571428571428571428571428571;
    
    rungeKutta1412Coefficients.bCoefficients( 1, 0 ) = 0.0178571428571428571428571428571428571428571428571428571428571;
    rungeKutta1412Coefficients.bCoefficients( 1, 1 ) = 0.00585937500000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.bCoefficients( 1, 2 ) = 0.0117187500000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.bCoefficients( 1, 4 ) = 0.0175781250000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.bCoefficients( 1, 6 ) = 0.0234375000000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.bCoefficients( 1, 7 ) = 0.0292968750000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.bCoefficients( 1, 9 ) = 0.0351562500000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.bCoefficients( 1, 10 ) = 0.0410156250000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.bCoefficients( 1, 11 ) = 0.0468750000000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.bCoefficients( 1, 13 ) = 0.0527343750000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.bCoefficients( 1, 14 ) = 0.0585937500000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.bCoefficients( 1, 15 ) = 0.0644531250000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.bCoefficients( 1, 17 ) = 0.105352113571753019691496032887878162227673083080523884041670;
    rungeKutta1412Coefficients.bCoefficients( 1, 18 ) = 0.170561346241752182382120338553874085887555487802790804737501;
    rungeKutta1412Coefficients.bCoefficients( 1, 19 ) = 0.206229397329351940783526485701104894741914286259542454077972;
    rungeKutta1412Coefficients.bCoefficients( 1, 20 ) = 0.206229397329351940783526485701104894741914286259542454077972;
    rungeKutta1412Coefficients.bCoefficients( 1, 21 ) = 0.170561346241752182382120338553874085887555487802790804737501;
    rungeKutta1412Coefficients.bCoefficients( 1, 22 ) = 0.105352113571753019691496032887878162227673083080523884041670;
    rungeKutta1412Coefficients.bCoefficients( 1, 23 ) = -0.0644531250000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.bCoefficients( 1, 24 ) = -0.0585937500000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.bCoefficients( 1, 25 ) = -0.0527343750000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.bCoefficients( 1, 26 ) = -0.0468750000000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.bCoefficients( 1, 27 ) = -0.0410156250000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.bCoefficients( 1, 28 ) = -0.0351562500000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.bCoefficients( 1, 29 ) = -0.0292968750000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.bCoefficients( 1, 30 ) = -0.0234375000000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.bCoefficients( 1, 31 ) = -0.0175781250000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.bCoefficients( 1, 32 ) = -0.0117187500000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.bCoefficients( 1, 33 ) = -0.00585937500000000000000000000000000000000000000000000000000000;
    rungeKutta1412Coefficients.bCoefficients( 1, 34 ) = 0.0178571428571428571428571428571428571428571428571428571428571;

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
