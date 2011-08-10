/*! \file simpleLinearRegression.cpp
 *    Source file that defines the simple linear regression method in Tudat.
 *
 *    Path              : /Mathematics/Statistics/
 *    Version           : 4
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
 *    Last modified     : 10 August, 2011
 *
 *    References
 *      Press W.H., et al. Numerical Recipes in C++: The Art of
 *          Scientific Computing. Cambridge University Press, February 2002.
 *
 *    Notes
 *
 *    Copyright (c) 2010 Delft University of Technology.
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
 *      110802    K. Kumar          Added computeFitErrors() function and
 *                                  split code from computeFit(); changed
 *                                  filename and class name.
 *      110810    J. Leloux         Corrected doxygen documentation.
 */

// Include statements.
#include "simpleLinearRegression.h"

//! Default constructor.
SimpleLinearRegression::SimpleLinearRegression( )
    : chiSquared_( -0.0 ), coefficientOfConstantTerm_( -0.0 ),
      coefficientOfLinearTerm_( -0.0 ),
      standardDeviationOfCoefficientOfConstantTerm_( -0.0 ),
      standardDeviationOfCoefficientOfLinearTerm_( -0.0 ),
      sumOfDependentVariableData_( -0.0 ),
      sumOfIndependentVariableData_( -0.0 )
{
}

//! Default destructor.
SimpleLinearRegression::~SimpleLinearRegression( )
{
}

//! Set input data.
void SimpleLinearRegression::
setInputData( const map< double, double >& inputDataToFit )
{
    inputDataToFit_ = inputDataToFit;
}

//! Get coefficient of constant term of fit.
double& SimpleLinearRegression::getCoefficientOfConstantTerm( )
{
    return coefficientOfConstantTerm_;
}

//! Get coefficient of linear term of fit.
double& SimpleLinearRegression::getCoefficientOfLinearTerm( )
{
    return coefficientOfLinearTerm_;
}

//! Get chi-squared value.
double& SimpleLinearRegression::getChiSquared( )
{
    return chiSquared_;
}

//! Get standard deviation of coefficient of constant term of fit.
double& SimpleLinearRegression::getStandardDeviationOfCoefficientOfConstantTerm( )
{
    return standardDeviationOfCoefficientOfConstantTerm_;
}

//! Get standard deviation of coefficient of linear term of fit.
double& SimpleLinearRegression::getStandardDeviationOfCoefficientOfLinearTerm( )
{
    return standardDeviationOfCoefficientOfLinearTerm_;
}

//! Compute fit.
void SimpleLinearRegression::computeFit( )
{
    // Algorithm given here is taken from pg. 780 - 785 of (Press at al, 2002).

    // Declare local variables.
    // Declare and initialize temporary variable.
    double temporaryVariable_ = 0.0;

    // Declare and initialize sum of squared temporary variable values.
    sumOfTemporaryVariableSquared_ = 0.0;

    // Declare and initialize number of data points.
    unsigned int numberOfDataPoints_ = inputDataToFit_.size( );

    // Compute sum of independent and dependent variable data.
    sumIndependentVariableData_( );
    sumDependentVariableData_( );

    // Compute weighted mean of independent variable data.
    double meanOfIndependentVariableData_
            = sumOfIndependentVariableData_ / numberOfDataPoints_;

    // Loop over input data to fit.
    for ( iteratorInputDataToFit_ = inputDataToFit_.begin( );
          iteratorInputDataToFit_ != inputDataToFit_.end( );
          iteratorInputDataToFit_++ )
    {
        // Update value of temporary variable with independent variable value.
        temporaryVariable_ = iteratorInputDataToFit_->first
                - meanOfIndependentVariableData_;

        // Update sum of squared temporary variable value.
        sumOfTemporaryVariableSquared_ += pow( temporaryVariable_, 2.0 );

        // Update value of coefficient of linear term (b).
        coefficientOfLinearTerm_ += temporaryVariable_
                * iteratorInputDataToFit_->second;
    }

    // Solve for coefficient of linear term (b).
    coefficientOfLinearTerm_ /= sumOfTemporaryVariableSquared_;

    // Compute coefficient of constant term (a).
    coefficientOfConstantTerm_ = ( sumOfDependentVariableData_
                                   - ( sumOfIndependentVariableData_
                                       * coefficientOfLinearTerm_ ) )
                                 / numberOfDataPoints_;
}

//! Compute fit errors.
void SimpleLinearRegression::computeFitErrors( )
{
    // Declare local variables.
    // Declare and initialize number of data points.
    unsigned int numberOfDataPoints_ = inputDataToFit_.size( );

    // Declare standard deviation obtained from chi-squared distribution.
    double standardDeviationFromChiSquared_ = 0.0;

    // Compute standard deviation of coefficient of constant term (sigma_a).
    standardDeviationOfCoefficientOfConstantTerm_
            = sqrt( ( 1.0
                      + ( pow( sumOfIndependentVariableData_, 2.0 )
                          / ( numberOfDataPoints_
                              * sumOfTemporaryVariableSquared_ ) ) )
                    / numberOfDataPoints_ );

    // Compute standard deviation of coefficient of linear term (sigma_b).
    standardDeviationOfCoefficientOfLinearTerm_
            = sqrt( 1.0 / sumOfTemporaryVariableSquared_ );

    // Compute chi-squared.
    for ( iteratorInputDataToFit_ = inputDataToFit_.begin( );
          iteratorInputDataToFit_ != inputDataToFit_.end( );
          iteratorInputDataToFit_++ )
    {
        chiSquared_ += pow( iteratorInputDataToFit_->second
                            - coefficientOfConstantTerm_
                            - ( coefficientOfLinearTerm_
                                * iteratorInputDataToFit_->first ), 2.0 );
    }

    // Adjust standard deviations based on typical standard deviation from
    // chi-squared distribution.
    if ( numberOfDataPoints_ > 2 )
    {
        standardDeviationFromChiSquared_
                = sqrt( chiSquared_
                        / static_cast< double >( numberOfDataPoints_ - 2 ) );
    }

    standardDeviationOfCoefficientOfConstantTerm_
            += standardDeviationFromChiSquared_;

    standardDeviationOfCoefficientOfLinearTerm_
            += standardDeviationFromChiSquared_;
}

//! Sum independent variable input data.
void SimpleLinearRegression::sumIndependentVariableData_( )
{
    for ( iteratorInputDataToFit_ = inputDataToFit_.begin( );
          iteratorInputDataToFit_ != inputDataToFit_.end( );
          iteratorInputDataToFit_++ )
    {
        sumOfIndependentVariableData_ += iteratorInputDataToFit_->first;
    }
}

//! Sum dependent variable input data.
void SimpleLinearRegression::sumDependentVariableData_( )
{
    for ( iteratorInputDataToFit_ = inputDataToFit_.begin( );
          iteratorInputDataToFit_ != inputDataToFit_.end( );
          iteratorInputDataToFit_++ )
    {
        sumOfDependentVariableData_ += iteratorInputDataToFit_->second;
    }
}

// End of file.
