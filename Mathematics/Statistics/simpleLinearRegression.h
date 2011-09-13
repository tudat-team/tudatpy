/*! \file simpleLinearRegression.h
 *    Header file that defines the simple linear regression method in Tudat.
 *
 *    Path              : /Mathematics/Statistics/
 *    Version           : 3
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : J.C.P Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 1 July, 2011
 *    Last modified     : 2 August, 2011
 *
 *    References
 *      Press W.H., et al. Numerical Recipes in C++: The Art of
 *          Scientific Computing. Cambridge University Press, February 2002.
 *
 *    Notes
 *      The independent variable input data must be distinct for this method to
 *      work. This method does not include individual weights for the
 *      observations. This must be added in future. In addition, a full Linear
 *      Least Squares method must be added in future that uses the general
 *      formulation of the normal equations.
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
 *      110701    K. Kumar          First creation of code.
 *      110726    K. Kumar          Changed filename and class name.
 *      110802    K. Kumar          Added computeFitErrors() function; removed
 *                                  DataModeling class inheritance; changed
 *                                  filename and class name.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 */

#ifndef SIMPLELINEARREGRESSION_H
#define SIMPLELINEARREGRESSION_H

// Include statements.
#include <map>

//! Simple linear regression class.
/*!
 * Simple linear regression class.
 */
class SimpleLinearRegression
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    SimpleLinearRegression( );

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~SimpleLinearRegression( ) { }

    //! Set input data.
    /*!
     * Sets input data to be fitted in a map.
     * \param inputDataToFit Input data to fit.
     */
    void setInputData( const std::map< double, double >& inputDataToFit ) 
    {
        inputDataToFit_ = inputDataToFit;
    }

    //! Get coefficient of constant term of fit.
    /*!
     * Returns coefficient of constant term of fit.
     * \return Coefficient of constant term.
     */
    double& getCoefficientOfConstantTerm( ) 
    {
        return coefficientOfConstantTerm_;
    }

    //! Get coefficient of linear term of fit.
    /*!
     * Returns coefficient of linear term of fit.
     * \return Coefficient of linear term.
     */
    double& getCoefficientOfLinearTerm( )
    {
        return coefficientOfLinearTerm_;
    }

    //! Get chi-squared value.
    /*!
     * Returns chi-squared value of fit.
     * \return Chi-squared value.
     */
    double& getChiSquared( )
    {
        return chiSquared_;
    }

    //! Get standard deviation of coefficient of constant term of fit.
    /*!
     * Returns standard deviation of coefficient of constant term of fit.
     * \return Standard deviation of coefficient of constant term.
     */
    double& getStandardDeviationOfCoefficientOfConstantTerm( )
    {
        return standardDeviationOfCoefficientOfConstantTerm_;
    }

    //! Get standard deviation of coefficient of linear term of fit.
    /*!
     * Returns standard deviation of coefficient of linear term of fit.
     * \return Standard deviation of coefficient of linear term.
     */
    double& getStandardDeviationOfCoefficientOfLinearTerm( )
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
     * computed linear fit. This functin must be called after the computeFit()
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
    std::map< double, double > inputDataToFit_;

    //! Iterator for input data to fit.
    /*!
     * Map iterator for input data to fit.
     */
    std::map< double, double >::iterator iteratorInputDataToFit_;

    //! Sum independent variable input data.
    /*!
     * Computes sum of independent variable input data. The independent
     * variable data is stored as the key in the map inputDataToFit_. This is
     * an auxilliary function for the computeFit() function.
     */
    void sumIndependentVariableData_( );

    //! Sum dependent variable input data.
    /*!
     * Computes sum of dependent variable input data. The dependent variable
     * data is stored as the value in the map inputDataToFit_. This is
     * an auxilliary function for the computeFit() function.
     */
    void sumDependentVariableData_( );
};

#endif // LINEARREGRESSION_H

// End of file.
