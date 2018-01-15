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
 *      Press W.H., et al. Numerical Recipes in C++: The Art of Scientific Computing. Cambridge
 *          University Press, February 2002.
 *
 */

#include <cmath>

#include "Tudat/Mathematics/Statistics/simpleLinearRegression.h"

namespace tudat
{
namespace statistics
{

//! Compute fit.
void SimpleLinearRegression::computeFit( )
{
    // Algorithm given here is taken from pg. 780 - 785 of (Press at al, 2002).

    // Initialize sum of squared temporary variable values.
    sumOfTemporaryVariableSquared_ = 0.0;

    // Declare and initialize number of data points.
    std::size_t numberOfDataPoints_ = inputDataToFit_.size( );

    // Compute sum of independent and dependent variable data.
    sumIndependentVariableData_( );
    sumDependentVariableData_( );

    // Compute weighted mean of independent variable data.
    double meanOfIndependentVariableData_ = sumOfIndependentVariableData_ / numberOfDataPoints_;

    // Set coefficient to zero.
    coefficientOfLinearTerm_ = 0.0;
    coefficientOfConstantTerm_ = 0.0;

    // Loop over input data to fit.
    for ( InputDataMap::iterator iteratorInputDataToFit_ = inputDataToFit_.begin( );
          iteratorInputDataToFit_ != inputDataToFit_.end( ); iteratorInputDataToFit_++ )
    {
        // Update value of temporary variable with independent variable value.
        double temporaryVariable_ = iteratorInputDataToFit_->first - meanOfIndependentVariableData_;

        // Update sum of squared temporary variable value.
        sumOfTemporaryVariableSquared_ += temporaryVariable_ * temporaryVariable_;

        // Update value of coefficient of linear term (b).
        coefficientOfLinearTerm_ += temporaryVariable_ * iteratorInputDataToFit_->second;
    }

    // Solve for coefficient of linear term (b).
    coefficientOfLinearTerm_ /= sumOfTemporaryVariableSquared_;

    // Compute coefficient of constant term (a).
    coefficientOfConstantTerm_ = ( sumOfDependentVariableData_
                                   - ( sumOfIndependentVariableData_ * coefficientOfLinearTerm_ ) )
            / numberOfDataPoints_;
}

//! Compute fit errors.
void SimpleLinearRegression::computeFitErrors( )
{
    // Declare local variables.
    // Declare and initialize number of data points.
    std::size_t numberOfDataPoints_ = inputDataToFit_.size( );

    // Declare standard deviation obtained from chi-squared distribution.
    double standardDeviationFromChiSquared_ = 0.0;

    // Compute standard deviation of coefficient of constant term (sigma_a).
    standardDeviationOfCoefficientOfConstantTerm_
            = std::sqrt( ( 1.0 + ( sumOfIndependentVariableData_ * sumOfIndependentVariableData_
                                   / ( numberOfDataPoints_ * sumOfTemporaryVariableSquared_ ) ) )
                         / numberOfDataPoints_ );

    // Compute standard deviation of coefficient of linear term (sigma_b).
    standardDeviationOfCoefficientOfLinearTerm_
            = std::sqrt( 1.0 / sumOfTemporaryVariableSquared_ );

    // Compute chi-squared.
    chiSquared_ = 0.0;
    for ( InputDataMap::iterator iteratorInputDataToFit_ = inputDataToFit_.begin( );
          iteratorInputDataToFit_ != inputDataToFit_.end( ); iteratorInputDataToFit_++ )
    {
        chiSquared_ += std::pow( iteratorInputDataToFit_->second - coefficientOfConstantTerm_
                                 - ( coefficientOfLinearTerm_
                                     * iteratorInputDataToFit_->first ), 2.0 );
    }

    // Adjust standard deviations based on typical standard deviation from
    // chi-squared distribution.
    if ( numberOfDataPoints_ > 2 )
    {
        standardDeviationFromChiSquared_
                = std::sqrt( chiSquared_ / static_cast< double >( numberOfDataPoints_ - 2 ) );
    }

    standardDeviationOfCoefficientOfConstantTerm_ *= standardDeviationFromChiSquared_;

    standardDeviationOfCoefficientOfLinearTerm_ *= standardDeviationFromChiSquared_;
}

//! Sum independent variable input data.
void SimpleLinearRegression::sumIndependentVariableData_( )
{
    sumOfIndependentVariableData_ = 0.0;

    for ( InputDataMap::iterator iteratorInputDataToFit_ = inputDataToFit_.begin( );
          iteratorInputDataToFit_ != inputDataToFit_.end( ); iteratorInputDataToFit_++ )
    {
        sumOfIndependentVariableData_ += iteratorInputDataToFit_->first;
    }
}

//! Sum dependent variable input data.
void SimpleLinearRegression::sumDependentVariableData_( )
{
    sumOfDependentVariableData_ = 0.0;

    for ( InputDataMap::iterator iteratorInputDataToFit_ = inputDataToFit_.begin( );
          iteratorInputDataToFit_ != inputDataToFit_.end( ); iteratorInputDataToFit_++ )
    {
        sumOfDependentVariableData_ += iteratorInputDataToFit_->second;
    }
}

} // namespace statistics
} // namespace tudat
