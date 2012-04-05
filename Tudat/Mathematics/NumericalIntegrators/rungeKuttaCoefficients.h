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
 *      120203    B. Tong Minh      File created.
 *      120327    K. Kumar          Added Runge-Kutta 87 (Dormand and Prince) enum option; added
 *                                  lower-, higher-order, and order to integrate variables.
 *
 *    References
 *      Burden, R.L., Faires, J.D. Numerical Analysis, 7th Edition, Books/Cole, 2001.
 *
 */

#ifndef TUDAT_RUNGE_KUTTA_COEFFICIENTS_H
#define TUDAT_RUNGE_KUTTA_COEFFICIENTS_H

#include <Eigen/Core>

namespace tudat
{
namespace mathematics
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

typedef struct RungeKuttaCoefficients RungeKuttaCoefficients;

} // namespace integrators
} // namespace mathematics
} // namespace tudat

#endif // TUDAT_RUNGE_KUTTA_COEFFICIENTS_H
