/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */


#include "baseFunctionsHodographicShaping.h"
#include <math.h>


namespace tudat
{
namespace shape_based_methods
{


//! Function, derivative, and integral of the sine function.

double SineFunctionHodographicShaping::evaluateFunction( double independentVariable )
{
    return std::sin( frequency_ * independentVariable );
}

double SineFunctionHodographicShaping::evaluateDerivative( double independentVariable )
{
    return frequency_ * std::cos( frequency_ * independentVariable );
}

double SineFunctionHodographicShaping::evaluateIntegral( double independentVariable )
{
    return - std::cos( frequency_ * independentVariable ) / frequency_;
}


//! Function, derivative, and integral of the cosine function.

double CosineFunctionHodographicShaping::evaluateFunction( double independentVariable )
{
    return std::cos( frequency_ * independentVariable );
}

double CosineFunctionHodographicShaping::evaluateDerivative( double independentVariable )
{
    return - frequency_ * std::sin( frequency_ * independentVariable );
}

double CosineFunctionHodographicShaping::evaluateIntegral( double independentVariable )
{
    return std::sin( frequency_ * independentVariable ) / frequency_;
}


//! Function, derivative, and integral of the exponential function.

double ExponentialFunctionHodographicShaping::evaluateFunction( double independentVariable )
{
    return std::exp( exponent_ * independentVariable );
}

double ExponentialFunctionHodographicShaping::evaluateDerivative( double independentVariable )
{
    return exponent_ * std::exp( exponent_ * independentVariable );
}

double ExponentialFunctionHodographicShaping::evaluateIntegral( double independentVariable )
{
    return std::exp( exponent_ * independentVariable ) / exponent_;
}


//! Function, derivative, and integral of the scaled exponential function.
double ScaledExponentialFunctionHodographicShaping::evaluateFunction( double independentVariable )
{
    return std::exp( scaleFactor_ * exponent_ * independentVariable );
}

double ScaledExponentialFunctionHodographicShaping::evaluateDerivative( double independentVariable )
{
    return scaleFactor_ * exponent_ * std::exp( scaleFactor_ * exponent_ * independentVariable );
}

double ScaledExponentialFunctionHodographicShaping::evaluateIntegral( double independentVariable )
{
    return std::exp( scaleFactor_ * exponent_ * independentVariable ) / ( scaleFactor_ * exponent_ );
}


//! Function, derivative, and integral of exponential times sine function.

double ExponentialSineFunctionHodographicShaping::evaluateFunction( double independentVariable )
{
    return std::exp( exponentExponentialFunction_ * independentVariable )
           * std::sin( frequencySineFunction_ * independentVariable );
}

double ExponentialSineFunctionHodographicShaping::evaluateDerivative( double independentVariable )
{
    return std::exp( exponentExponentialFunction_ * independentVariable )
           * ( exponentExponentialFunction_ * std::sin( frequencySineFunction_ * independentVariable )
               + frequencySineFunction_ * std::cos( frequencySineFunction_ * independentVariable ) );
}

double ExponentialSineFunctionHodographicShaping::evaluateIntegral( double independentVariable )
{
    return std::exp( exponentExponentialFunction_ * independentVariable )
            / ( std::pow( exponentExponentialFunction_, 2.0 ) + std::pow( frequencySineFunction_, 2.0 ) )
            * ( - frequencySineFunction_ * std::cos( frequencySineFunction_ * independentVariable )
                + exponentExponentialFunction_ * std::sin( frequencySineFunction_ * independentVariable ) );
}


//! Function, derivative, and integral of scaled exponential times sine function.

double ScaledExponentialSineFunctionHodographicShaping::evaluateFunction( double independentVariable )
{
    return std::exp( scaleFactor_ * exponentExponentialFunction_ * independentVariable )
           * std::sin( frequencySineFunction_ * independentVariable );
}

double ScaledExponentialSineFunctionHodographicShaping::evaluateDerivative( double independentVariable )
{
    return std::exp( scaleFactor_ * exponentExponentialFunction_ * independentVariable )
            * ( scaleFactor_ * exponentExponentialFunction_ * std::sin( frequencySineFunction_ * independentVariable )
                + frequencySineFunction_ * std::cos( frequencySineFunction_ * independentVariable ) );
}

double ScaledExponentialSineFunctionHodographicShaping::evaluateIntegral( double independentVariable )
{
    return std::exp( scaleFactor_ * exponentExponentialFunction_ * independentVariable )
           / ( std::pow( scaleFactor_ * exponentExponentialFunction_, 2.0 ) + std::pow( frequencySineFunction_, 2.0 ) )
           * ( - frequencySineFunction_ * std::cos( frequencySineFunction_ * independentVariable )
               + scaleFactor_ * exponentExponentialFunction_ * std::sin( frequencySineFunction_ * independentVariable ) );
}


//! Function, derivative, and integral of exponential times cosine function.

double ExponentialCosineFunctionHodographicShaping::evaluateFunction( double independentVariable )
{
    return std::exp( exponentExponentialFunction_ * independentVariable )
           * std::cos( frequencyCosineFunction_ * independentVariable );
}

double ExponentialCosineFunctionHodographicShaping::evaluateDerivative( double independentVariable )
{
    return std::exp( exponentExponentialFunction_ * independentVariable )
           * ( exponentExponentialFunction_ * std::cos( frequencyCosineFunction_ * independentVariable )
               - frequencyCosineFunction_ * std::sin( frequencyCosineFunction_ * independentVariable ) );
}

double ExponentialCosineFunctionHodographicShaping::evaluateIntegral( double independentVariable )
{
    return std::exp( exponentExponentialFunction_ * independentVariable )
            / ( std::pow( exponentExponentialFunction_, 2.0 ) + std::pow( frequencyCosineFunction_, 2.0 ) )
            * ( frequencyCosineFunction_ * std::sin( frequencyCosineFunction_ * independentVariable )
                + exponentExponentialFunction_ * std::cos( frequencyCosineFunction_ * independentVariable ) );
}


//! Function, derivative, and integral of the scaled exponential times cosine function.

double ScaledExponentialCosineFunctionHodographicShaping::evaluateFunction( double independentVariable )
{
    return std::exp( scaleFactor_ * exponentExponentialFunction_ * independentVariable )
           * std::cos( frequencyCosineFunction_ * independentVariable );
}

double ScaledExponentialCosineFunctionHodographicShaping::evaluateDerivative( double independentVariable )
{
    return std::exp( scaleFactor_ * exponentExponentialFunction_ * independentVariable )
            * ( scaleFactor_ * exponentExponentialFunction_ * std::cos( frequencyCosineFunction_ * independentVariable )
                - frequencyCosineFunction_ * std::sin( frequencyCosineFunction_ * independentVariable ) );
}

double ScaledExponentialCosineFunctionHodographicShaping::evaluateIntegral( double independentVariable )
{
    return std::exp( scaleFactor_ * exponentExponentialFunction_ * independentVariable )
           / ( std::pow( scaleFactor_ * exponentExponentialFunction_, 2.0 ) + std::pow( frequencyCosineFunction_, 2.0 ) )
           * ( frequencyCosineFunction_ * std::sin( frequencyCosineFunction_ * independentVariable )
               + scaleFactor_ * exponentExponentialFunction_ * std::cos( frequencyCosineFunction_ * independentVariable ) );
}


//! Function, derivative, and integral of the power function.

double PowerFunctionHodographicShaping::evaluateFunction( double independentVariable )
{
    return std::pow( independentVariable , exponent_ );
}

double PowerFunctionHodographicShaping::evaluateDerivative( double independentVariable )
{
    return exponent_ * std::pow( independentVariable , ( exponent_ - 1.0 ) );
}

double PowerFunctionHodographicShaping::evaluateIntegral( double independentVariable )
{
    return std::pow( independentVariable , ( exponent_ + 1.0 ) ) / ( exponent_ + 1.0 );
}


//! Function, derivative, and integral of scaled power function.

double ScaledPowerFunctionHodographicShaping::evaluateFunction( double independentVariable )
{
    return scaleFactor_ * std::pow( independentVariable , exponent_ );
}

double ScaledPowerFunctionHodographicShaping::evaluateDerivative( double independentVariable )
{
    return scaleFactor_ * exponent_ * std::pow( independentVariable , ( exponent_ - 1.0 ) );
}

double ScaledPowerFunctionHodographicShaping::evaluateIntegral( double independentVariable )
{
    return scaleFactor_ * std::pow( independentVariable , ( exponent_ + 1.0 ) ) / ( exponent_ + 1.0 );
}


//! Function, derivative, and integral of the power times sine function.

double PowerSineFunctionHodographicShaping::evaluateFunction( double independentVariable )
{
    return std::pow( independentVariable, exponentPowerFunction_ )
           * std::sin( frequencySineFunction_ * independentVariable );
}

double PowerSineFunctionHodographicShaping::evaluateDerivative( double independentVariable )
{
    return exponentPowerFunction_ * std::pow( independentVariable, exponentPowerFunction_ - 1.0 )
           * std::sin( frequencySineFunction_ * independentVariable )
           + std::pow( independentVariable, exponentPowerFunction_ )
           * frequencySineFunction_ * std::cos( frequencySineFunction_ * independentVariable ) ;
}

double PowerSineFunctionHodographicShaping::evaluateIntegral( double independentVariable )
{

    double integralValue = 0.0;
    double integrationScaleFactor = 1.0;

    // Integration by parts.
    for ( int counterSuccessiveIntegrations = 1; counterSuccessiveIntegrations <= exponentPowerFunction_ + 1 ; counterSuccessiveIntegrations++ )
    {
        double doubleCounterSuccessiveIntegrations = static_cast<double>( counterSuccessiveIntegrations );

        if ( ( counterSuccessiveIntegrations % 2 ) == 0.0 )
        {
            integralValue += std::pow( - 1.0, doubleCounterSuccessiveIntegrations / 2.0 + 1.0 ) * integrationScaleFactor
                             * std::pow( independentVariable, exponentPowerFunction_ - doubleCounterSuccessiveIntegrations + 1.0 )
                             / std::pow( frequencySineFunction_, doubleCounterSuccessiveIntegrations )
                             * std::sin( frequencySineFunction_ * independentVariable );
        }
        else
        {
            integralValue += std::pow( - 1.0, ( doubleCounterSuccessiveIntegrations + 1.0 ) / 2.0 ) * integrationScaleFactor
                             * std::pow( independentVariable, exponentPowerFunction_ - doubleCounterSuccessiveIntegrations + 1.0 )
                             / std::pow( frequencySineFunction_, doubleCounterSuccessiveIntegrations )
                             * std::cos( frequencySineFunction_ * independentVariable );
        }

        integrationScaleFactor *= exponentPowerFunction_ - doubleCounterSuccessiveIntegrations + 1.0;
    }

    return integralValue;
}


//! Function, derivative, and integral of scaled power times sine function.

double ScaledPowerSineFunctionHodographicShaping::evaluateFunction( double independentVariable )
{
    return scaleFactor_ * std::pow( independentVariable, exponentPowerFunction_ )
           * std::sin( frequencySineFunction_ * independentVariable );
}

double ScaledPowerSineFunctionHodographicShaping::evaluateDerivative( double independentVariable )
{
    return scaleFactor_ *
           ( exponentPowerFunction_ * std::pow( independentVariable, exponentPowerFunction_ - 1.0 )
               * std::sin( frequencySineFunction_ * independentVariable )
             + std::pow( independentVariable, exponentPowerFunction_ )
               * frequencySineFunction_ * std::cos( frequencySineFunction_ * independentVariable ) ) ;
}

double ScaledPowerSineFunctionHodographicShaping::evaluateIntegral( double independentVariable )
{
    double integralValue = 0.0;
    double integrationScaleFactor = 1.0;

    // Integration by parts.
    for ( int counterSuccessiveIntegrations = 1 ; counterSuccessiveIntegrations <= exponentPowerFunction_ + 1 ; counterSuccessiveIntegrations++ )
    {
        double doubleCounterSuccessiveIntegrations = static_cast< double >( counterSuccessiveIntegrations );
        if ( ( counterSuccessiveIntegrations % 2 ) == 0.0 )
        {
            integralValue += std::pow( - 1.0, doubleCounterSuccessiveIntegrations / 2.0 + 1.0 ) * integrationScaleFactor
                             * std::pow( independentVariable, exponentPowerFunction_ - doubleCounterSuccessiveIntegrations + 1.0 )
                             / std::pow( frequencySineFunction_, doubleCounterSuccessiveIntegrations )
                             * std::sin( frequencySineFunction_ * independentVariable );
        }
        else
        {
            integralValue += std::pow( - 1.0, ( doubleCounterSuccessiveIntegrations + 1.0 ) / 2.0 ) * integrationScaleFactor
                             * std::pow( independentVariable, exponentPowerFunction_ - doubleCounterSuccessiveIntegrations + 1.0 )
                             / std::pow( frequencySineFunction_, doubleCounterSuccessiveIntegrations )
                             * std::cos( frequencySineFunction_ * independentVariable );
        }
        integrationScaleFactor *= exponentPowerFunction_ - doubleCounterSuccessiveIntegrations + 1.0;
    }
    return scaleFactor_ * integralValue;
}


//!  Function, derivative, and integral of the power times cosine function.

double PowerCosineFunctionHodographicShaping::evaluateFunction( double independentVariable )
{
    return std::pow( independentVariable, exponentPowerFunction_ )
           * std::cos( frequencyCosineFunction_ * independentVariable );
}

double PowerCosineFunctionHodographicShaping::evaluateDerivative( double independentVariable )
{
    return exponentPowerFunction_ * std::pow( independentVariable, exponentPowerFunction_ - 1.0 )
           * std::cos( frequencyCosineFunction_ * independentVariable )
           - std::pow( independentVariable, exponentPowerFunction_ )
           * frequencyCosineFunction_ * std::sin( frequencyCosineFunction_ * independentVariable ) ;
}

double PowerCosineFunctionHodographicShaping::evaluateIntegral( double independentVariable )
{

    double integralValue = 0.0;
    double integrationScaleFactor = 1.0;

    // Integration by parts.
    for ( int counterSuccessiveIntegrations = 1; counterSuccessiveIntegrations <= exponentPowerFunction_+1 ; counterSuccessiveIntegrations++ )
    {
        double doubleCounterSuccessiveIntegration = static_cast<double>( counterSuccessiveIntegrations );

        if ( ( counterSuccessiveIntegrations % 2 ) == 0.0 )
        {
            integralValue += std::pow( - 1.0, doubleCounterSuccessiveIntegration / 2.0 + 1.0 ) * integrationScaleFactor
                             * std::pow( independentVariable, exponentPowerFunction_ - doubleCounterSuccessiveIntegration + 1.0 )
                             / std::pow( frequencyCosineFunction_, doubleCounterSuccessiveIntegration )
                             * std::cos( frequencyCosineFunction_ * independentVariable );
        }
        else
        {
            integralValue += std::pow( - 1.0, ( doubleCounterSuccessiveIntegration - 1.0 ) / 2.0 ) * integrationScaleFactor
                             * std::pow( independentVariable, exponentPowerFunction_ - doubleCounterSuccessiveIntegration + 1.0 )
                             / std::pow( frequencyCosineFunction_, doubleCounterSuccessiveIntegration )
                             * std::sin( frequencyCosineFunction_ * independentVariable );
        }

        integrationScaleFactor *= exponentPowerFunction_ - doubleCounterSuccessiveIntegration + 1.0;
    }

    return integralValue;
}


//! Function, derivative, and integral of scaled power times cosine function.

double ScaledPowerCosineFunctionHodographicShaping::evaluateFunction( double independentVariable )
{
    return scaleFactor_ * std::pow( independentVariable, exponentPowerFunction_ )
           * std::cos( frequencyCosineFunction_ * independentVariable );
}

double ScaledPowerCosineFunctionHodographicShaping::evaluateDerivative( double independentVariable )
{
    return scaleFactor_ *
           ( exponentPowerFunction_ * std::pow( independentVariable, exponentPowerFunction_ - 1.0 )
               * std::cos( frequencyCosineFunction_ * independentVariable )
             - std::pow( independentVariable, exponentPowerFunction_ )
               * frequencyCosineFunction_ * std::sin( frequencyCosineFunction_ * independentVariable ) );
}

double ScaledPowerCosineFunctionHodographicShaping::evaluateIntegral( double independentVariable )
{
    double integralValue = 0.0;
    double integrationScaleFactor = 1.0;

    // Integration by parts.
    for ( int counterSuccessiveIntegrations = 1 ; counterSuccessiveIntegrations <= exponentPowerFunction_ + 1 ; counterSuccessiveIntegrations++ )
    {
        double doubleCounterSuccessiveIntegrations = static_cast< double >( counterSuccessiveIntegrations );
        if ( ( counterSuccessiveIntegrations % 2 ) == 0.0 )
        {
            integralValue += std::pow( - 1.0, doubleCounterSuccessiveIntegrations / 2.0 + 1.0 ) * integrationScaleFactor
                             * std::pow( independentVariable, exponentPowerFunction_ - doubleCounterSuccessiveIntegrations + 1.0 )
                             / std::pow( frequencyCosineFunction_, doubleCounterSuccessiveIntegrations )
                             * std::cos( frequencyCosineFunction_ * independentVariable );
        }
        else
        {
            integralValue += std::pow( - 1.0, ( doubleCounterSuccessiveIntegrations - 1.0 ) / 2.0 ) * integrationScaleFactor
                             * std::pow( independentVariable, exponentPowerFunction_ - doubleCounterSuccessiveIntegrations + 1.0 )
                             / std::pow( frequencyCosineFunction_, doubleCounterSuccessiveIntegrations )
                             * std::sin( frequencyCosineFunction_ * independentVariable );
        }
        integrationScaleFactor *= exponentPowerFunction_ - doubleCounterSuccessiveIntegrations + 1.0;
    }
    return scaleFactor_ * integralValue;
}


} // namespace shape_based_methods
} // namespace tudat
