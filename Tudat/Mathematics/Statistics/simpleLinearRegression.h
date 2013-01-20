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
 *      110701    K. Kumar          Creation of code.
 *      110726    K. Kumar          Changed filename and class name.
 *      110802    K. Kumar          Added computeFitErrors( ) function; removed
 *                                  DataModeling class inheritance; changed
 *                                  filename and class name.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *      120508    K. Kumar          Updated code so input data is set at construction.
 *      120521    A. Ronse          Updated namespaces.
 *
 *    References
 *      Press W.H., et al. Numerical Recipes in C++: The Art of Scientific Computing. Cambridge
 *          University Press, February 2002.
 *
 *    Notes
 *      This method does not include individual weights for the observations. The independent
 *      variable input data must be distinct for this method work. In addition, a full Linear Least
 *      Squares (LLS) method must be added in future that uses the general formulation of the
 *      normal equations.
 *
 */

#ifndef TUDAT_SIMPLE_LINEAR_REGRESSION_H
#define TUDAT_SIMPLE_LINEAR_REGRESSION_H

#include <map>

#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

namespace tudat
{
namespace statistics
{

//! Simple linear regression class.
/*!
 * Simple linear regression class.
 */
class SimpleLinearRegression
{
public:

    //! Type definition of map of input data.
    /*!
     * Type definition of map of input data.
     */
    typedef std::map< double, double > InputDataMap;

    //! Default constructor.
    /*!
     * Default constructor.
     */
    SimpleLinearRegression( const InputDataMap& inputDataToFit )
        : chiSquared_( TUDAT_NAN ),
          coefficientOfConstantTerm_( TUDAT_NAN ),
          coefficientOfLinearTerm_( TUDAT_NAN ),
          standardDeviationOfCoefficientOfConstantTerm_( TUDAT_NAN ),
          standardDeviationOfCoefficientOfLinearTerm_( TUDAT_NAN ),
          sumOfDependentVariableData_( TUDAT_NAN ),
          sumOfIndependentVariableData_( TUDAT_NAN ),
          sumOfTemporaryVariableSquared_( TUDAT_NAN ),
          inputDataToFit_( inputDataToFit )
    { }

    //! Get coefficient of constant term of fit.
    /*!
     * Returns coefficient of constant term of fit.
     * \return Coefficient of constant term.
     */
    double getCoefficientOfConstantTerm( ) const { return coefficientOfConstantTerm_; }

    //! Get coefficient of linear term of fit.
    /*!
     * Returns coefficient of linear term of fit.
     * \return Coefficient of linear term.
     */
    double getCoefficientOfLinearTerm( ) const { return coefficientOfLinearTerm_; }

    //! Get chi-squared value.
    /*!
     * Returns chi-squared value of fit.
     * \return Chi-squared value.
     */
    double getChiSquared( ) const { return chiSquared_; }

    //! Get standard deviation of coefficient of constant term of fit.
    /*!
     * Returns standard deviation of coefficient of constant term of fit.
     * \return Standard deviation of coefficient of constant term.
     */
    double getStandardDeviationOfCoefficientOfConstantTerm( ) const
    {
        return standardDeviationOfCoefficientOfConstantTerm_;
    }

    //! Get standard deviation of coefficient of linear term of fit.
    /*!
     * Returns standard deviation of coefficient of linear term of fit.
     * \return Standard deviation of coefficient of linear term.
     */
    double getStandardDeviationOfCoefficientOfLinearTerm( ) const
    {
        return standardDeviationOfCoefficientOfLinearTerm_;
    }

    //! Compute fit.
    /*!
     * Computes simple linear regression fitting with unit weights for all the
     * observations. The equation of a straight line used here is given as:
     * \f[
     *      y = a + b * x
     * \f]
     * where \f$ b \f$ is the coefficient of the linear term, and \f$ a \f$ is
     * the coefficient of the constant term. The algorithm implemented is based
     * on (Press W.H., et al., 2002).
     */
    void computeFit( );

    //! Compute fit errors.
    /*!
     * Computes the standard deviations and chi-squared value based on the
     * computed linear fit. This functin must be called after the computeFit( )
     * function. The algorithm implemented is based on
     * (Press W.H., et al., 2002).
     */
    void computeFitErrors( );

protected:

private:

    //! Chi-squared value.
    /*!
     * Chi-squared value of simple linear regression fit.
     */
    double chiSquared_;

    //! Coefficient of constant term.
    /*!
     * Coefficient of constant term of simple linear regression fit. This
     * corresponds to the coefficient \f$ a \f$.
     */
    double coefficientOfConstantTerm_;

    //! Coefficient of linear term.
    /*!
     * Coefficient of linear term of simple linear regression fit. This
     * corresponds to the coefficient \f$ b \f$.
     */
    double coefficientOfLinearTerm_;

    //! Standard deviation of coefficient of constant term.
    /*!
     * Standard deviation of coefficient of constant term of simple linear
     * regression fit.
     */
    double standardDeviationOfCoefficientOfConstantTerm_;

    //! Standard deviation of coefficient of linear term.
    /*!
     * Standard deviation of coefficient of linear term of simple linear
     * regression fit.
     */
    double standardDeviationOfCoefficientOfLinearTerm_;

    //! Sum of dependent variable data.
    /*!
     * Sum of dependent variable data.
     */
    double sumOfDependentVariableData_;

    //! Sum of independent variable data.
    /*!
     * Sum of independent variable data.
     */
    double sumOfIndependentVariableData_;

    //! Sum of temporary variables squared.
    /*!
     * Sum of temporary variables squared, used by the simple linear regression
     * fit algorithm.
     */
    double sumOfTemporaryVariableSquared_;

    //! Input data to fit.
    /*!
     * Input data to fit.
     */
    InputDataMap inputDataToFit_;

    //! Sum independent variable input data.
    /*!
     * Computes sum of independent variable input data. The independent
     * variable data is stored as the key in the map inputDataToFit_. This is
     * an auxilliary function for the computeFit( ) function.
     */
    void sumIndependentVariableData_( );

    //! Sum dependent variable input data.
    /*!
     * Computes sum of dependent variable input data. The dependent variable
     * data is stored as the value in the map inputDataToFit_. This is
     * an auxilliary function for the computeFit( ) function.
     */
    void sumDependentVariableData_( );
};

} // namespace statistics
} // namespace tudat

#endif // TUDAT_SIMPLE_LINEAR_REGRESSION_H
