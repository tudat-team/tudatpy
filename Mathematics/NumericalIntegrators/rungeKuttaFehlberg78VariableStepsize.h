/*! \file rungeKuttaFehlberg78VariableStepsize.h
 *    Header file that defines the 7(8)th-order, variable stepsize, Runge-Kutta
 *    integrator.
 *
 *    Path              : /Mathematics/NumericalIntegrators/
 *    Version           : 3
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Author            : E.A.G. Heeren
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : E.A.G.Heeren@student.tudelft.nl
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl

 *    Checker           : E.A.G. Heeren
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : E.A.G.Heeren@student.tudelft.nl
 *
 *    Date created      : 26 August, 2011
 *    Last modified     : 12 September, 2011
 *
 *    References
 *      Eagle, C.D. Celestial Computing with MATLAB, Science Software,
 *          http://www.cdeagle.com/ccmatlab/ccmatlab.pdf, last accessed:
 *          28 August, 2011.
 *
 *    Notes
 *
 *    Copyright (c) 2010-2011 Delft University of Technology.
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
 *      110826    K. Kumar          File created.
 *      110909    E.A.G. Heeren     Modified to be included in Tudat. Added if-statement
 *                                  to ensure a maximum stepsize change.
 *      110912    K. Kumar          Minor changes.
 */

#ifndef RUNGEKUTTAFEHLBERG78VARIABLESTEPSIZE_H
#define RUNGEKUTTAFEHLBERG78VARIABLESTEPSIZE_H

// Include statements.
#include <Eigen/Core>
#include "Mathematics/NumericalIntegrators/integrator.h"
#include "Astrodynamics/States/state.h"

//! Tudat library namespace.
/*!
 * The Tudat library namespace.
 */
namespace tudat
{

//! 7th/8th-order, variable stepsize, Runge-Kutta integrator class.
/*!
 * Implementation of 7th/8th-order, variable stepsize, Runge-Kutta integrator.
 */
class RungeKuttaFehlberg78VariableStepsize : public tudat::Integrator
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    RungeKuttaFehlberg78VariableStepsize( ) :
        minimumStepsize_( 1.0e-15 ), relativeErrorTolerance_( 1.0e-12 ), currentTime_( -0.0 ),
        errorInState_( -0.0 ), previousStepsize_( -0.0 ), truncationError_( -0.0 ),
        dampedAbsoluteTolerance_( -0.0 ), relativeTruncationError_( -0.0 ),
        aCoefficients_( Eigen::MatrixXd::Zero( 13, 13 ) ),
        bCoefficients_( Eigen::VectorXd::Zero( 13 ) ),
        cCoefficients_( Eigen::VectorXd::Zero( 13 ) ), fMatrix( ) { setCoefficients_( ); }

    //! Set relative error tolerance.
    /*!
     * Sets the relative error tolerance.
     * \param relativeErrorTolerance Relative error tolerance.
     */
    void setRelativeErrorTolerance( double relativeErrorTolerance )
    {
        // Set relative error tolerance.
        relativeErrorTolerance_ = relativeErrorTolerance;
    }

    //! Set minimum stepsize.
    /*!
     * Sets the minimum stepsize.
     * \param minimumStepsize Minimum stepsize.
     */
    void setMinimumStepsize( double minimumStepsize )
    {
        // Set minimumStepsize.
        minimumStepsize_ = minimumStepsize;
    }

    //! Integrate.
    /*!
     * This function executes an integration.
     */
    void integrate( );

    //! Overload ostream to print class information.
    /*!
     * Overloads ostream to print class information.
     * \param stream Stream object.
     * \param rungeKuttaFehlberg78VariableStepsize Runge-Kutta 7(8)th-order
     *        variable-stepsize integrator.
     * \return Stream object.
     */
    friend std::ostream& operator<<( std::ostream& stream,
                                     RungeKuttaFehlberg78VariableStepsize&
                                     rungeKuttaFehlberg78VariableStepsize );

protected:

private:

    //! Minimum stepsize during integration.
    /*!
     * Minimum stepsize during integration.
     */
    double minimumStepsize_;

    //! Relative error tolerance.
    /*!
    * The relative error tolerance of the 7(th)-order Runge-Kutta-Fehlberg Method.
    */
    double relativeErrorTolerance_;

    //! Current time.
    /*!
     * Current time.
     */
    double currentTime_;

    //! Error in state.
    /*!
     * Error in state.
     */
    double errorInState_;

    //! Previous stepsize.
    /*!
     * Previous stepsize.
     */
    double previousStepsize_;

    //! Truncation error.
    /*!
     * Truncation error.
     */
    double truncationError_;

    //! Damped absolute error tolerance.
    /*!
     * Damped absolute error tolerance.
     */
    double dampedAbsoluteTolerance_;

    //! Relative truncation error.
    /*!
     * Relative truncation error.
     */
    double relativeTruncationError_;

    //! a-Coefficients.
    /*!
    * MatrixXd containing the a-coefficients for the Runge-Kutta-Fehlberg
    * 7(8)th-order integration method.
    */
    Eigen::MatrixXd aCoefficients_;

    //! b-Coefficients.
    /*!
    * RowVectorXd containing the b-coefficients for the Runge-Kutta-Fehlberg
    * 7(8)th-order integration method.
    */
    Eigen::RowVectorXd bCoefficients_;

    //! c-Coefficients.
    /*!
    * VectorXd containing the c-coefficients for the Runge-Kutta-Fehlberg
    * 7(8)th-order integration method.
    */
    Eigen::VectorXd cCoefficients_;

    //! f-Matrix.
    /*!
    * f-Matrix for the Runge-Kutta-Fehlberg 7(8)th-order integration method.
    */
    Eigen::MatrixXd fMatrix;

    //! Current state.
    /*!
     * Current state.
     */
    tudat::State currentState_;

    //! Set coefficients of integration scheme.
    /*!
     * Sets the coefficients of the 7th-/8-th order, variable stepsize,
     * Runge-Kutta integration scheme.
     */
    void setCoefficients_( );
};

}

#endif // RUNGEKUTTAFEHLBERG78VARIABLESTEPSIZE_H

// End of file.
