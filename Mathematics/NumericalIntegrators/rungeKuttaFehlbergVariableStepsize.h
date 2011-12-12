/*! \file variableStepsizeRungeKuttaIntegrator.h
 *    Header file that defines the general, variable stepsize, Runge-Kutta
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
 *    Author            : F.M. Engelen
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : F.M.Engelen@student.tudelft.nl
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl

 *    Checker           : E.A.G. Heeren
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : E.A.G.Heeren@student.tudelft.nl
 *
 *    Date created      : 26 August, 2011
 *    Last modified     : 10 November, 2011
 *
 *    References
 *      Eagle, C.D. Celestial Computing with MATLAB, Science Software,
 *          http://www.cdeagle.com/ccmatlab/ccmatlab.pdf, last accessed: 28 August, 2011.
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
 *      110909    E.A.G. Heeren     Modified to be included in Tudat. Added if-statement to ensure
 *                                  a maximum stepsize change.
 *      110912    K. Kumar          Minor changes.
 *      111021    F.M. Engelen      Added interface for RKF45 and RKF56.
 *      111110    E.A.G Heeren      Minor changes.
 */

#ifndef RUNGEKUTTAFEHLBERGVARIABLESTEPSIZE_H
#define RUNGEKUTTAFEHLBERGVARIABLESTEPSIZE_H

// Include statements.
#include <Eigen/Core>
#include "Astrodynamics/States/state.h"
#include "Mathematics/NumericalIntegrators/integrator.h"

//! Tudat library namespace.
/*!
 * The Tudat library namespace.
 */
namespace tudat
{

//! The variable stepsize integrator class.
/*!
 * Implementation of the variable stepsize integrator.
 */
class RungeKuttaFehlBergVariableStepsizeIntegrator : public Integrator
{
public:

    //! Enum of the different Runge-Kutta-Fehlberg, variable stepsize, integrators available.
    /*!
     * Enum of the different integrators available.
     */
    enum RungeKuttaFehlbergVariableStepsizeIntegratorType
    {
        rungeKuttaFelhberg45VariableStepsize, rungeKuttaFelhberg56VariableStepsize,
        rungeKuttaFelhberg78VariableStepsize
    };

    //! Default constructor.
    /*!
     * Default constructor.
     */
    RungeKuttaFehlBergVariableStepsizeIntegrator(
            RungeKuttaFehlbergVariableStepsizeIntegratorType
            rungeKuttaFehlbergVariableStepsizeIntegratorType );

    //! Set relative error tolerance.
    /*!
     * Sets the relative error tolerance.
     * \param relativeErrorTolerance Relative error tolerance.
     */
    void setRelativeErrorTolerance( double relativeErrorTolerance )
    { relativeErrorTolerance_ = relativeErrorTolerance; }

    //! Set minimum stepsize.
    /*!
     * Sets the minimum stepsize.
     * \param minimumStepsize Minimum stepsize.
     */
    void setMinimumStepsize( double minimumStepsize ) { minimumStepsize_ = minimumStepsize; }

    //! Integrate.
    /*!
     * This function executes an integration.
     */
    void integrate( );

    //! Overload ostream to print class information.
    /*!
     * Overloads ostream to print class information.
     * \param stream Stream object.
     * \param rungeKuttaFehlbergVariableStepsize Runge-Kutta-Fehlberg variable-stepize integrator.
     * \return Stream object.
     */
    friend std::ostream& operator<<( std::ostream& stream,
                                     RungeKuttaFehlBergVariableStepsizeIntegrator&
                                     rungeKuttaFehlBergVariableStepsizeIntegrator );

protected:

private:

    //! Integrator order.
    /*!
     * The order of the intergrator, which is used to propagate.
     */
    unsigned int integratorOrder_;

    //! Number of stages of the integrator.
    /*!
     * The number of stages of the integrator.
     */
    unsigned int numberOfStages_;

    //! Current time.
    /*!
     * Current time.
     */
    double currentTime_;

    //! Damped absolute error tolerance.
    /*!
     * Damped absolute error tolerance.
     */
    double dampedAbsoluteTolerance_;

    //! Error in state.
    /*!
     * Error in state.
     */
    double errorInState_;

    //! Minimum stepsize during integration.
    /*!
     * Minimum stepsize during integration.
     */
    double minimumStepsize_;

    //! Previous stepsize.
    /*!
     * Previous stepsize.
     */
    double previousStepsize_;

    //! Relative error tolerance.
    /*!
    * The relative error tolerance of the integration Method.
    */
    double relativeErrorTolerance_;

    //! Relative truncation error.
    /*!
     * Relative truncation error.
     */
    double relativeTruncationError_;

    //! Truncation error.
    /*!
     * Truncation error.
     */
    double truncationError_;

    //! a-Coefficients.
    /*!
    * MatrixXd containing the a-coefficients for the Runge-Kutta integration method.
    */
    Eigen::MatrixXd aCoefficients_;

    //! b-Coefficients.
    /*!
    * MatrixXd containing the b-coefficients for the Runge-Kutta integration method.
    */
    Eigen::MatrixXd bCoefficients_;

    //! c-Coefficients.
    /*!
    * VectorXd containing the c-coefficients for the Runge-Kutta-Fehlberg
    * 7(8)th-order integration method.
    */
    Eigen::VectorXd cCoefficients_;

    //! f-Matrix.
    /*!
    * f-Matrix for the General Runge-Kutta integration method.
    */
    Eigen::MatrixXd fMatrix;

    //! Current state.
    /*!
     * Current state.
     */
    State currentState_;

    //! Set coefficients of integration scheme.
    /*!
     * Sets the coefficients of the variable stepsize,
     * Runge-Kutta integration scheme.
     */
    void setCoefficients_( RungeKuttaFehlbergVariableStepsizeIntegratorType coefficientSet );
};

}

#endif // RUNGEKUTTAFEHLBERGVARIABLESTEPSIZE_H

// End of file.
