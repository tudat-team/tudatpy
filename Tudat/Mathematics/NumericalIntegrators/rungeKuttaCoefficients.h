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
 *      Burden, R.L., Faires, J.D. Numerical Analysis, 7th Edition, Books/Cole, 2001.
 *
 */

#ifndef TUDAT_RUNGE_KUTTA_COEFFICIENTS_H
#define TUDAT_RUNGE_KUTTA_COEFFICIENTS_H

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

namespace tudat
{
namespace numerical_integrators
{

//! Struct that defines the coefficients of a Runge-Kutta integrator
/*!
 * Struct that defines coefficients of a Runge-Kutta integrators
 */
struct RungeKuttaCoefficients
{
    //! Enum of order estimates that can be integrated.
    enum OrderEstimateToIntegrate { lower, higher };

    //! Main table of the Butcher tableau.
    Eigen::MatrixXd aCoefficients;

    //! Bottom rows of the Butcher tableau.
    Eigen::MatrixXd bCoefficients;

    //! First column of the Butcher tableau.
    Eigen::VectorXd cCoefficients;

    //! Order of the higher order estimate.
    unsigned int higherOrder;

    //! Order of the lower order estimate.
    unsigned int lowerOrder;

    //! Order estimate to integrate.
    OrderEstimateToIntegrate orderEstimateToIntegrate;

    //! Default constructor.
    /*!
     * Default constructor that initializes coefficients to 0.
     */
    RungeKuttaCoefficients( ) :
        aCoefficients( ),
        bCoefficients( ),
        cCoefficients( ),
        higherOrder( 0 ),
        lowerOrder( 0 ),
        orderEstimateToIntegrate( lower )
    { }

    //! Constructor.
    /*!
     * Constructor that sets the coefficients.
     * \param aCoefficients_ Main table of the Butcher tableau.
     * \param bCoefficients_ Bottom rows of the Butcher tableau.
     * \param cCoefficients_ First column of the Butcher tableau.
     * \param higherOrder_ Order of the integrator.
     * \param lowerOrder_ Order of the embedded low-order integrator.
     * \param order Enum denoting whether to use the lower or higher order scheme for numerical
     * integration.
     */
    RungeKuttaCoefficients( const Eigen::MatrixXd& aCoefficients_,
                            const Eigen::MatrixXd& bCoefficients_,
                            const Eigen::MatrixXd& cCoefficients_,
                            const unsigned int higherOrder_,
                            const unsigned int lowerOrder_,
                            OrderEstimateToIntegrate order ) :
        aCoefficients( aCoefficients_ ),
        bCoefficients( bCoefficients_ ),
        cCoefficients( cCoefficients_ ),
        higherOrder( higherOrder_ ),
        lowerOrder( lowerOrder_ ),
        orderEstimateToIntegrate( order )
    { }

    //! Enum of predefined coefficient sets.
    enum CoefficientSets
    {
        undefinedCoefficientSet = -1,
        rungeKuttaFehlberg45,
        rungeKuttaFehlberg56,
        rungeKuttaFehlberg78,
        rungeKutta87DormandPrince
    };

    //! Get coefficients for a specified coefficient set.
    /*!
     * Returns coefficients for a specified coefficient set.
     * \param coefficientSet The set to get the coefficients for.
     * \return The requested coefficient set.
     */
    static const RungeKuttaCoefficients& get( CoefficientSets coefficientSet );
};

//! Typedef for shared-pointer to RungeKuttaCoefficients object.
typedef boost::shared_ptr< RungeKuttaCoefficients > RungeKuttaCoefficientsPointer;

} // namespace numerical_integrators
} // namespace tudat

#endif // TUDAT_RUNGE_KUTTA_COEFFICIENTS_H
