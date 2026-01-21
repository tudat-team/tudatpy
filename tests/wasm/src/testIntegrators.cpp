/*    Copyright (c) 2010-2024, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Additional integrator WASM tests - RK78 and Bulirsch-Stoer.
 */

#include "wasmTestFramework.h"

// Mathematics
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/math/integrators/rungeKuttaVariableStepSizeIntegrator.h"
#include "tudat/math/integrators/bulirschStoerVariableStepsizeIntegrator.h"
#include "tudat/math/integrators/adamsBashforthMoultonIntegrator.h"

using namespace tudat;

void testRungeKutta78Integrator()
{
    std::cout << "\n=== Runge-Kutta 78 (Fehlberg) Integrator ===" << std::endl;

    using namespace numerical_integrators;

    // Test 1: Exponential growth (same as RK4 test but with adaptive stepping)
    // Reference: dy/dt = y, y(0) = 1, exact: y = e^t
    {
        auto stateDerivative = [](const double t, const Eigen::VectorXd& state) -> Eigen::VectorXd {
            return state;
        };

        Eigen::VectorXd initialState(1);
        initialState << 1.0;

        // Create RK78 integrator with adaptive stepping
        RungeKuttaCoefficients rkCoeffs =
            RungeKuttaCoefficients::get(CoefficientSets::rungeKuttaFehlberg78);

        // Constructor args: (coeffs, func, t0, state, minStep, maxStep, initialStep, relTol, absTol)
        RungeKuttaVariableStepSizeIntegrator<double, Eigen::VectorXd> integrator(
            rkCoeffs, stateDerivative, 0.0, initialState,
            std::numeric_limits<double>::epsilon(),  // min step
            std::numeric_limits<double>::infinity(), // max step
            1e-6,   // initial step
            1e-15, 1e-15);  // relative and absolute tolerance

        // Integrate to t = 1 using integrateTo() method
        Eigen::VectorXd finalState = integrator.integrateTo(1.0, 1e-6);

        double computedY = finalState(0);
        double exactY = std::exp(1.0);

        checkClose("RKF78 exponential growth", computedY, exactY, 1e-12);
    }

    // Test 2: Harmonic oscillator with adaptive stepping
    {
        auto harmonicDerivative = [](const double t, const Eigen::VectorXd& state) -> Eigen::VectorXd {
            Eigen::VectorXd derivative(2);
            derivative(0) = state(1);
            derivative(1) = -state(0);
            return derivative;
        };

        Eigen::VectorXd initialState(2);
        initialState << 1.0, 0.0;  // x=1, v=0

        RungeKuttaCoefficients rkCoeffs =
            RungeKuttaCoefficients::get(CoefficientSets::rungeKuttaFehlberg78);

        // Constructor args: (coeffs, func, t0, state, minStep, maxStep, initialStep, relTol, absTol)
        RungeKuttaVariableStepSizeIntegrator<double, Eigen::VectorXd> integrator(
            rkCoeffs, harmonicDerivative, 0.0, initialState,
            std::numeric_limits<double>::epsilon(),  // min step
            std::numeric_limits<double>::infinity(), // max step
            0.01,   // initial step
            1e-12, 1e-12);  // relative and absolute tolerance

        // Integrate to t = 2*pi (one full period) using integrateTo() method
        double tEnd = 2.0 * mathematical_constants::PI;
        Eigen::VectorXd finalState = integrator.integrateTo(tEnd, 0.01);

        // After one period, should return to initial state
        double computedX = finalState(0);
        double computedV = finalState(1);

        checkClose("RKF78 harmonic oscillator x after 2pi", computedX, 1.0, 1e-8);
        checkClose("RKF78 harmonic oscillator v after 2pi", computedV, 0.0, 1e-8);
    }
}

void testRungeKutta87DormandPrinceIntegrator()
{
    std::cout << "\n=== Runge-Kutta 87 (Dormand-Prince) Integrator ===" << std::endl;

    using namespace numerical_integrators;

    // Test 1: Exponential growth - dy/dt = y, y(0) = 1, exact: y = e^t
    {
        auto stateDerivative = [](const double t, const Eigen::VectorXd& state) -> Eigen::VectorXd {
            return state;
        };

        Eigen::VectorXd initialState(1);
        initialState << 1.0;

        RungeKuttaCoefficients rkCoeffs =
            RungeKuttaCoefficients::get(CoefficientSets::rungeKutta87DormandPrince);

        RungeKuttaVariableStepSizeIntegrator<double, Eigen::VectorXd> integrator(
            rkCoeffs, stateDerivative, 0.0, initialState,
            std::numeric_limits<double>::epsilon(),
            std::numeric_limits<double>::infinity(),
            1e-6, 1e-15, 1e-15);

        Eigen::VectorXd finalState = integrator.integrateTo(1.0, 1e-6);

        double computedY = finalState(0);
        double exactY = std::exp(1.0);

        checkClose("RKDP87 exponential growth", computedY, exactY, 1e-13);
    }

    // Test 2: Harmonic oscillator
    {
        auto harmonicDerivative = [](const double t, const Eigen::VectorXd& state) -> Eigen::VectorXd {
            Eigen::VectorXd derivative(2);
            derivative(0) = state(1);
            derivative(1) = -state(0);
            return derivative;
        };

        Eigen::VectorXd initialState(2);
        initialState << 1.0, 0.0;

        RungeKuttaCoefficients rkCoeffs =
            RungeKuttaCoefficients::get(CoefficientSets::rungeKutta87DormandPrince);

        RungeKuttaVariableStepSizeIntegrator<double, Eigen::VectorXd> integrator(
            rkCoeffs, harmonicDerivative, 0.0, initialState,
            std::numeric_limits<double>::epsilon(),
            std::numeric_limits<double>::infinity(),
            0.01, 1e-12, 1e-12);

        double tEnd = 2.0 * mathematical_constants::PI;
        Eigen::VectorXd finalState = integrator.integrateTo(tEnd, 0.01);

        checkClose("RKDP87 harmonic oscillator x after 2pi", finalState(0), 1.0, 1e-9);
        checkClose("RKDP87 harmonic oscillator v after 2pi", finalState(1), 0.0, 1e-9);
    }

    // Test 3: Non-autonomous ODE from native test
    // dy/dt = y - t^2 + 1 (used in Matlab benchmark data)
    {
        auto nonAutonomousDerivative = [](const double t, const Eigen::VectorXd& state) -> Eigen::VectorXd {
            Eigen::VectorXd derivative(1);
            derivative(0) = state(0) - t * t + 1.0;
            return derivative;
        };

        // Initial condition from native test: t0 = 0.0, y0 = 0.5
        Eigen::VectorXd initialState(1);
        initialState << 0.5;

        RungeKuttaCoefficients rkCoeffs =
            RungeKuttaCoefficients::get(CoefficientSets::rungeKutta87DormandPrince);

        RungeKuttaVariableStepSizeIntegrator<double, Eigen::VectorXd> integrator(
            rkCoeffs, nonAutonomousDerivative, 0.0, initialState,
            std::numeric_limits<double>::epsilon(),
            std::numeric_limits<double>::infinity(),
            0.1, 1e-14, 1e-14);

        // Integrate to t = 2.0
        Eigen::VectorXd finalState = integrator.integrateTo(2.0, 0.1);

        // Exact solution: y(t) = (t+1)^2 - 0.5*e^t
        // y(2) = 9 - 0.5*e^2 = 9 - 3.694528... = 5.305471...
        double exactY = 9.0 - 0.5 * std::exp(2.0);

        checkClose("RKDP87 non-autonomous ODE", finalState(0), exactY, 1e-12);
    }

    // Test 4: Compare RKDP87 vs RKF78 (both should give similar results)
    {
        auto stateDerivative = [](const double t, const Eigen::VectorXd& state) -> Eigen::VectorXd {
            Eigen::VectorXd derivative(2);
            derivative(0) = state(1);
            derivative(1) = -state(0) - 0.1 * state(1);  // Damped oscillator
            return derivative;
        };

        Eigen::VectorXd initialState(2);
        initialState << 1.0, 0.0;

        // RKDP87
        RungeKuttaCoefficients dpCoeffs =
            RungeKuttaCoefficients::get(CoefficientSets::rungeKutta87DormandPrince);
        RungeKuttaVariableStepSizeIntegrator<double, Eigen::VectorXd> dpIntegrator(
            dpCoeffs, stateDerivative, 0.0, initialState,
            std::numeric_limits<double>::epsilon(),
            std::numeric_limits<double>::infinity(),
            0.01, 1e-14, 1e-14);

        // RKF78
        RungeKuttaCoefficients fCoeffs =
            RungeKuttaCoefficients::get(CoefficientSets::rungeKuttaFehlberg78);
        RungeKuttaVariableStepSizeIntegrator<double, Eigen::VectorXd> fIntegrator(
            fCoeffs, stateDerivative, 0.0, initialState,
            std::numeric_limits<double>::epsilon(),
            std::numeric_limits<double>::infinity(),
            0.01, 1e-14, 1e-14);

        double tEnd = 10.0;
        Eigen::VectorXd dpFinal = dpIntegrator.integrateTo(tEnd, 0.01);
        Eigen::VectorXd fFinal = fIntegrator.integrateTo(tEnd, 0.01);

        checkClose("RKDP87 vs RKF78 x agreement", dpFinal(0), fFinal(0), 1e-10);
        checkClose("RKDP87 vs RKF78 v agreement", dpFinal(1), fFinal(1), 1e-10);
    }
}

void testRungeKuttaFehlberg45Integrator()
{
    std::cout << "\n=== Runge-Kutta-Fehlberg 45 Integrator ===" << std::endl;

    using namespace numerical_integrators;

    // Test 1: Exponential growth - dy/dt = y, y(0) = 1, exact: y = e^t
    {
        auto stateDerivative = [](const double t, const Eigen::VectorXd& state) -> Eigen::VectorXd {
            return state;
        };

        Eigen::VectorXd initialState(1);
        initialState << 1.0;

        RungeKuttaCoefficients rkCoeffs =
            RungeKuttaCoefficients::get(CoefficientSets::rungeKuttaFehlberg45);

        RungeKuttaVariableStepSizeIntegrator<double, Eigen::VectorXd> integrator(
            rkCoeffs, stateDerivative, 0.0, initialState,
            std::numeric_limits<double>::epsilon(),
            std::numeric_limits<double>::infinity(),
            1e-4, 1e-12, 1e-12);

        Eigen::VectorXd finalState = integrator.integrateTo(1.0, 1e-4);

        double computedY = finalState(0);
        double exactY = std::exp(1.0);

        checkClose("RKF45 exponential growth", computedY, exactY, 2e-10);
    }

    // Test 2: Harmonic oscillator
    {
        auto harmonicDerivative = [](const double t, const Eigen::VectorXd& state) -> Eigen::VectorXd {
            Eigen::VectorXd derivative(2);
            derivative(0) = state(1);
            derivative(1) = -state(0);
            return derivative;
        };

        Eigen::VectorXd initialState(2);
        initialState << 1.0, 0.0;

        RungeKuttaCoefficients rkCoeffs =
            RungeKuttaCoefficients::get(CoefficientSets::rungeKuttaFehlberg45);

        RungeKuttaVariableStepSizeIntegrator<double, Eigen::VectorXd> integrator(
            rkCoeffs, harmonicDerivative, 0.0, initialState,
            std::numeric_limits<double>::epsilon(),
            std::numeric_limits<double>::infinity(),
            0.01, 1e-10, 1e-10);

        double tEnd = 2.0 * mathematical_constants::PI;
        Eigen::VectorXd finalState = integrator.integrateTo(tEnd, 0.01);

        checkClose("RKF45 harmonic oscillator x after 2pi", finalState(0), 1.0, 1e-6);
        checkClose("RKF45 harmonic oscillator v after 2pi", finalState(1), 0.0, 1e-6);
    }

    // Test 3: Lotka-Volterra (predator-prey) equations
    // dx/dt = ax - bxy (prey)
    // dy/dt = dxy - cy (predator)
    {
        const double a = 1.1, b = 0.4, c = 0.4, d = 0.1;

        auto loktaVolterraDerivative = [a, b, c, d](const double t, const Eigen::VectorXd& state) -> Eigen::VectorXd {
            Eigen::VectorXd derivative(2);
            double x = state(0);  // prey
            double y = state(1);  // predator
            derivative(0) = a * x - b * x * y;
            derivative(1) = d * x * y - c * y;
            return derivative;
        };

        Eigen::VectorXd initialState(2);
        initialState << 10.0, 10.0;  // Initial populations

        RungeKuttaCoefficients rkCoeffs =
            RungeKuttaCoefficients::get(CoefficientSets::rungeKuttaFehlberg45);

        RungeKuttaVariableStepSizeIntegrator<double, Eigen::VectorXd> integrator(
            rkCoeffs, loktaVolterraDerivative, 0.0, initialState,
            std::numeric_limits<double>::epsilon(),
            std::numeric_limits<double>::infinity(),
            0.01, 1e-10, 1e-10);

        // Integrate for 50 time units
        Eigen::VectorXd finalState = integrator.integrateTo(50.0, 0.01);

        // Populations should remain positive and bounded
        checkTrue("RKF45 Lotka-Volterra prey > 0", finalState(0) > 0.0);
        checkTrue("RKF45 Lotka-Volterra predator > 0", finalState(1) > 0.0);
        checkTrue("RKF45 Lotka-Volterra prey < 100", finalState(0) < 100.0);
        checkTrue("RKF45 Lotka-Volterra predator < 100", finalState(1) < 100.0);
    }
}

void testBulirschStoerIntegrator()
{
    std::cout << "\n=== Bulirsch-Stoer Variable Step-Size Integrator ===" << std::endl;

    using namespace numerical_integrators;

    // Test 1: Fixed-step Bulirsch-Stoer for simple exponential ODE
    // The fixed-step constructor avoids the recursive step rejection that can cause
    // stack overflow in WASM due to limited stack size.
    // Reference: dy/dt = y, y(0) = 1, exact: y = e^t
    {
        auto stateDerivative = [](const double t, const Eigen::VectorXd& state) -> Eigen::VectorXd {
            return state;
        };

        Eigen::VectorXd initialState(1);
        initialState << 1.0;
        double t0 = 0.0;
        double tEnd = 1.0;
        double stepSize = 0.1;

        // Use the fixed-step constructor (2-argument step size)
        // This avoids the step size controller which can cause recursive calls
        BulirschStoerVariableStepSizeIntegrator<double, Eigen::VectorXd> bsIntegrator(
            getBulirschStoerStepSequence(bulirsch_stoer_sequence, 6),
            stateDerivative, t0, initialState, stepSize);

        Eigen::VectorXd bsFinalState = bsIntegrator.integrateTo(tEnd, stepSize);

        double computed = bsFinalState(0);
        double exact = std::exp(1.0);

        // BS with fixed steps should still be very accurate due to extrapolation
        checkClose("BS fixed-step exponential growth", computed, exact, 1e-8);
    }

    // Test 2: Fixed-step Bulirsch-Stoer for harmonic oscillator
    {
        auto harmonicDerivative = [](const double t, const Eigen::VectorXd& state) -> Eigen::VectorXd {
            Eigen::VectorXd deriv(2);
            deriv(0) = state(1);
            deriv(1) = -state(0);
            return deriv;
        };

        Eigen::VectorXd initialState(2);
        initialState << 1.0, 0.0;
        double t0 = 0.0;
        double tEnd = 2.0 * mathematical_constants::PI;
        double stepSize = 0.2;

        BulirschStoerVariableStepSizeIntegrator<double, Eigen::VectorXd> bsIntegrator(
            getBulirschStoerStepSequence(bulirsch_stoer_sequence, 6),
            harmonicDerivative, t0, initialState, stepSize);

        Eigen::VectorXd bsFinalState = bsIntegrator.integrateTo(tEnd, stepSize);

        // After one period, should return close to initial state
        checkClose("BS fixed-step harmonic x after 2pi", bsFinalState(0), 1.0, 1e-6);
        checkClose("BS fixed-step harmonic v after 2pi", bsFinalState(1), 0.0, 1e-6);
    }

    // Test 3: Compare BS with RK78 using relaxed tolerances to avoid deep recursion
    // Use larger tolerances (1e-8) to prevent repeated step rejections
    {
        auto stateDerivative = [](const double t, const Eigen::VectorXd& state) -> Eigen::VectorXd {
            Eigen::VectorXd deriv(1);
            deriv(0) = state(0) - t * t + 1.0;
            return deriv;
        };

        Eigen::VectorXd initialState(1);
        initialState << 0.5;
        double t0 = 0.5;
        double tEnd = 1.0;  // Shorter interval to reduce recursion risk

        // Use relaxed tolerances (1e-8) to minimize step rejections in WASM
        // Also use min/max step bounds that don't cause infinite recursion
        BulirschStoerVariableStepSizeIntegrator<double, Eigen::VectorXd> bsIntegrator(
            getBulirschStoerStepSequence(bulirsch_stoer_sequence, 4),
            stateDerivative, t0, initialState,
            0.001,   // min step
            0.5,     // max step
            0.05,    // initial step
            1e-8, 1e-8);  // relaxed tolerances

        Eigen::VectorXd bsFinalState = bsIntegrator.integrateTo(tEnd, 0.05);

        // Compare with RK78
        RungeKuttaCoefficients rkCoeffs =
            RungeKuttaCoefficients::get(CoefficientSets::rungeKuttaFehlberg78);

        RungeKuttaVariableStepSizeIntegrator<double, Eigen::VectorXd> rkIntegrator(
            rkCoeffs, stateDerivative, t0, initialState,
            0.001, 0.5, 0.05, 1e-8, 1e-8);

        Eigen::VectorXd rkFinalState = rkIntegrator.integrateTo(tEnd, 0.05);

        double bsResult = bsFinalState(0);
        double rkResult = rkFinalState(0);

        // Both should agree reasonably well
        checkClose("BS vs RK78 agreement (relaxed tol)", bsResult, rkResult, 1e-6);
    }
}

void testAdamsBashforthMoultonIntegrator()
{
    std::cout << "\n=== Adams-Bashforth-Moulton Integrator ===" << std::endl;

    using namespace numerical_integrators;

    // Test 1: Exponential growth - dy/dt = y, y(0) = 1, exact: y = e^t
    {
        auto stateDerivative = [](const double t, const Eigen::VectorXd& state) -> Eigen::VectorXd {
            return state;
        };

        Eigen::VectorXd initialState(1);
        initialState << 1.0;

        Eigen::VectorXd relTol(1), absTol(1);
        relTol << 1e-12;
        absTol << 1e-12;

        AdamsBashforthMoultonIntegrator<double, Eigen::VectorXd> integrator(
            stateDerivative, 0.0, initialState,
            1e-6,   // min step
            1.0,    // max step
            0.01,   // initial step
            relTol, absTol);

        Eigen::VectorXd finalState = integrator.integrateTo(1.0, 0.01);

        double computedY = finalState(0);
        double exactY = std::exp(1.0);

        checkClose("ABM exponential growth", computedY, exactY, 1e-10);
    }

    // Test 2: Harmonic oscillator
    {
        auto harmonicDerivative = [](const double t, const Eigen::VectorXd& state) -> Eigen::VectorXd {
            Eigen::VectorXd derivative(2);
            derivative(0) = state(1);
            derivative(1) = -state(0);
            return derivative;
        };

        Eigen::VectorXd initialState(2);
        initialState << 1.0, 0.0;

        Eigen::VectorXd relTol(2), absTol(2);
        relTol << 1e-12, 1e-12;
        absTol << 1e-12, 1e-12;

        AdamsBashforthMoultonIntegrator<double, Eigen::VectorXd> integrator(
            harmonicDerivative, 0.0, initialState,
            1e-6, 1.0, 0.01, relTol, absTol);

        double tEnd = 2.0 * mathematical_constants::PI;
        Eigen::VectorXd finalState = integrator.integrateTo(tEnd, 0.01);

        checkClose("ABM harmonic oscillator x after 2pi", finalState(0), 1.0, 1e-8);
        checkClose("ABM harmonic oscillator v after 2pi", finalState(1), 0.0, 1e-8);
    }

    // Test 3: Compare ABM vs RKF78 on damped oscillator
    {
        auto dampedDerivative = [](const double t, const Eigen::VectorXd& state) -> Eigen::VectorXd {
            Eigen::VectorXd derivative(2);
            derivative(0) = state(1);
            derivative(1) = -state(0) - 0.1 * state(1);
            return derivative;
        };

        Eigen::VectorXd initialState(2);
        initialState << 1.0, 0.0;

        // ABM
        Eigen::VectorXd relTol(2), absTol(2);
        relTol << 1e-12, 1e-12;
        absTol << 1e-12, 1e-12;

        AdamsBashforthMoultonIntegrator<double, Eigen::VectorXd> abmIntegrator(
            dampedDerivative, 0.0, initialState,
            1e-6, 1.0, 0.01, relTol, absTol);

        // RKF78
        RungeKuttaCoefficients rkCoeffs =
            RungeKuttaCoefficients::get(CoefficientSets::rungeKuttaFehlberg78);
        RungeKuttaVariableStepSizeIntegrator<double, Eigen::VectorXd> rkIntegrator(
            rkCoeffs, dampedDerivative, 0.0, initialState,
            std::numeric_limits<double>::epsilon(),
            std::numeric_limits<double>::infinity(),
            0.01, 1e-12, 1e-12);

        double tEnd = 10.0;
        Eigen::VectorXd abmFinal = abmIntegrator.integrateTo(tEnd, 0.01);
        Eigen::VectorXd rkFinal = rkIntegrator.integrateTo(tEnd, 0.01);

        checkClose("ABM vs RKF78 damped x", abmFinal(0), rkFinal(0), 1e-8);
        checkClose("ABM vs RKF78 damped v", abmFinal(1), rkFinal(1), 1e-8);
    }

    // Test 4: Non-autonomous ODE (same as RKDP87 test)
    {
        auto nonAutonomousDerivative = [](const double t, const Eigen::VectorXd& state) -> Eigen::VectorXd {
            Eigen::VectorXd derivative(1);
            derivative(0) = state(0) - t * t + 1.0;
            return derivative;
        };

        Eigen::VectorXd initialState(1);
        initialState << 0.5;

        Eigen::VectorXd relTol(1), absTol(1);
        relTol << 1e-12;
        absTol << 1e-12;

        AdamsBashforthMoultonIntegrator<double, Eigen::VectorXd> integrator(
            nonAutonomousDerivative, 0.0, initialState,
            1e-6, 1.0, 0.01, relTol, absTol);

        Eigen::VectorXd finalState = integrator.integrateTo(2.0, 0.01);

        // Exact solution: y(t) = (t+1)^2 - 0.5*e^t
        double exactY = 9.0 - 0.5 * std::exp(2.0);

        checkClose("ABM non-autonomous ODE", finalState(0), exactY, 1e-10);
    }
}
