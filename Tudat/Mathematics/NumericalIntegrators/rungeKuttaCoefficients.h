/*    Copyright (c) 2010-2013, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120203    B. Tong Minh      File created.
 *      120327    K. Kumar          Added Runge-Kutta 87 (Dormand and Prince) enum option; added
 *                                  lower-, higher-order, and order to integrate variables.
 *      130118    K. Kumar          Removed unused typedef.
 *      130121    K. Kumar          Added shared-ptr typedef.
 *      130916    K. Kumar          Removed RKF56.
 *
 *    References
 *      Burden, R.L., Faires, J.D. Numerical Analysis, 7th Edition, Books/Cole, 2001.
 *
 *    Notes
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
