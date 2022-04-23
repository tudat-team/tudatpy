/*    Copyright (c) 2010-2022, Delft University of Technology
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
 *
 */

#ifndef TUDAT_RUNGE_KUTTA_COEFFICIENTS_H
#define TUDAT_RUNGE_KUTTA_COEFFICIENTS_H

#include <memory>

#include <Eigen/Core>

#include <iostream>

namespace tudat
{
namespace numerical_integrators
{

// Enum of predefined coefficient sets.
//! @get_docstring(CoefficientSets.__docstring__)
enum CoefficientSets
{
    undefinedCoefficientSet = -1,
    forwardEuler,
    rungeKutta4Classic,
    explicitMidPoint,
    explicitTrapezoidRule,
    ralston,
    rungeKutta3,
    ralston3,
    SSPRK3,
    ralston4,
    threeEighthRuleRK4,
    heunEuler,
    rungeKuttaFehlberg12,
    rungeKuttaFehlberg45,
    rungeKuttaFehlberg56,
    rungeKuttaFehlberg78,
    rungeKutta87DormandPrince,
    rungeKuttaFehlberg89,
    rungeKuttaVerner89,
    rungeKuttaFeagin108,
    rungeKuttaFeagin1210,
    rungeKuttaFeagin1412
};

// Struct that defines the coefficients of a Runge-Kutta integrator
/*
 * Struct that defines coefficients of a Runge-Kutta integrators
 */
struct RungeKuttaCoefficients
{
    // Enum of order estimates that can be integrated.
    enum OrderEstimateToIntegrate { lower, higher };

    // Main table of the Butcher tableau.
    Eigen::MatrixXd aCoefficients;

    // Bottom rows of the Butcher tableau.
    Eigen::MatrixXd bCoefficients;

    // First column of the Butcher tableau.
    Eigen::VectorXd cCoefficients;

    // Order of the higher order estimate.
    unsigned int higherOrder;

    // Order of the lower order estimate.
    unsigned int lowerOrder;

    // Order estimate to integrate.
    OrderEstimateToIntegrate orderEstimateToIntegrate;

    // Save if the coefficient set corresponds to a fixed step size.
    bool isFixedStepSize;

    // Name of the coefficients.
    std::string name;

    // Default constructor.
    /*
     * Default constructor that initializes coefficients to 0.
     */
    RungeKuttaCoefficients( ) :
        aCoefficients( ),
        bCoefficients( ),
        cCoefficients( ),
        higherOrder( 0 ),
        lowerOrder( 0 ),
        orderEstimateToIntegrate( lower ),
        isFixedStepSize( false ),
        name( "Undefined" )        
    { }

    // Constructor.
    /*
     * Constructor that sets the coefficients.
     * \param aCoefficients_ Main table of the Butcher tableau.
     * \param bCoefficients_ Bottom rows of the Butcher tableau.
     * \param cCoefficients_ First column of the Butcher tableau.
     * \param higherOrder_ Order of the integrator.
     * \param lowerOrder_ Order of the embedded low-order integrator.
     * \param order Enum denoting whether to use the lower or higher order scheme for numerical
     * integration.
     * \param isFixedStepSize_ Boolean denoting whether the coefficient set is made for a fixed step size.
     * \param name_ Name of the coefficient set.
     */
    RungeKuttaCoefficients( const Eigen::MatrixXd& aCoefficients_,
                            const Eigen::MatrixXd& bCoefficients_,
                            const Eigen::MatrixXd& cCoefficients_,
                            const unsigned int higherOrder_,
                            const unsigned int lowerOrder_,
                            OrderEstimateToIntegrate order,
                            const bool isFixedStepSize_ = false,
                            std::string name_ = "Undefined" ) :
        aCoefficients( aCoefficients_ ),
        bCoefficients( bCoefficients_ ),
        cCoefficients( cCoefficients_ ),
        higherOrder( higherOrder_ ),
        lowerOrder( lowerOrder_ ),
        orderEstimateToIntegrate( order ),
        isFixedStepSize( isFixedStepSize_ ),
        name( name_ )
    { }

    // Get coefficients for a specified coefficient set.
    /*
     * Returns coefficients for a specified coefficient set.
     * \param coefficientSet The set to get the coefficients for.
     * \return The requested coefficient set.
     */
    static const RungeKuttaCoefficients& get( CoefficientSets coefficientSet );

};

// Typedef for shared-pointer to RungeKuttaCoefficients object.
typedef std::shared_ptr< RungeKuttaCoefficients > RungeKuttaCoefficientsPointer;

// Function to print the coefficients of a Runge-Kutta integrator.
void printButcherTableau( CoefficientSets coefficientSet );

} // namespace numerical_integrators
} // namespace tudat

#endif // TUDAT_RUNGE_KUTTA_COEFFICIENTS_H
