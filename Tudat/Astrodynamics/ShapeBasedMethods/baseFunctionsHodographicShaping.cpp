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

double SineFunction::evaluateFunction( double independentVariable )
{
    return std::sin( frequency_ * independentVariable );
}

double SineFunction::evaluateDerivative( double independentVariable )
{
    return frequency_ * std::cos( frequency_ * independentVariable );
}

double SineFunction::evaluateIntegral( double independentVariable )
{
    return - std::cos( frequency_ * independentVariable ) / frequency_;
}


//! Function, derivative, and integral of the cosine function.

double CosineFunction::evaluateFunction( double independentVariable )
{
    return std::cos( frequency_ * independentVariable );
}

double CosineFunction::evaluateDerivative( double independentVariable )
{
    return - frequency_ * std::sin( frequency_ * independentVariable );
}

double CosineFunction::evaluateIntegral( double independentVariable )
{
    return std::sin( frequency_ * independentVariable ) / frequency_;
}


//! Function, derivative, and integral of the exponential function.

double ExponentialFunction::evaluateFunction( double independentVariable )
{
    return std::exp( exponent_ * independentVariable );
}

double ExponentialFunction::evaluateDerivative( double independentVariable )
{
    return exponent_ * std::exp( exponent_ * independentVariable );
}

double ExponentialFunction::evaluateIntegral( double independentVariable )
{
    return std::exp( exponent_ * independentVariable ) / exponent_;
}


//! Function, derivative, and integral of the scaled exponential function.
double ScaledExponentialFunction::evaluateFunction( double independentVariable )
{
    return std::exp( scaleFactor_ * exponent_ * independentVariable );
}

double ScaledExponentialFunction::evaluateDerivative( double independentVariable )
{
    return scaleFactor_ * exponent_ * std::exp( scaleFactor_ * exponent_ * independentVariable );
}

double ScaledExponentialFunction::evaluateIntegral( double independentVariable )
{
    return std::exp( scaleFactor_ * exponent_ * independentVariable ) / ( scaleFactor_ * exponent_ );
}


//! Function, derivative, and integral of exponential times sine function.

double ExponentialSineFunction::evaluateFunction( double independentVariable )
{
    return std::exp( exponentExponentialFunction_ * independentVariable )
           * std::sin( frequencySineFunction_ * independentVariable );
}

double ExponentialSineFunction::evaluateDerivative( double independentVariable )
{
    return std::exp( exponentExponentialFunction_ * independentVariable )
           * ( exponentExponentialFunction_ * std::sin( frequencySineFunction_ * independentVariable )
               + frequencySineFunction_ * std::cos( frequencySineFunction_ * independentVariable ) );
}

double ExponentialSineFunction::evaluateIntegral( double independentVariable )
{
    return std::exp( exponentExponentialFunction_ * independentVariable )
            / ( std::pow( exponentExponentialFunction_, 2.0 ) + std::pow( frequencySineFunction_, 2.0 ) )
            * ( - frequencySineFunction_ * std::cos( frequencySineFunction_ * independentVariable )
                + exponentExponentialFunction_ * std::sin( frequencySineFunction_ * independentVariable ) );
}


//! Function, derivative, and integral of scaled exponential times sine function.

double ScaledExponentialSineFunction::evaluateFunction( double independentVariable )
{
    return std::exp( scaleFactor_ * exponentExponentialFunction_ * independentVariable )
           * std::sin( frequencySineFunction_ * independentVariable );
}

double ScaledExponentialSineFunction::evaluateDerivative( double independentVariable )
{
    return std::exp( scaleFactor_ * exponentExponentialFunction_ * independentVariable )
            * ( scaleFactor_ * exponentExponentialFunction_ * std::sin( frequencySineFunction_ * independentVariable )
                + frequencySineFunction_ * std::cos( frequencySineFunction_ * independentVariable ) );
}

double ScaledExponentialSineFunction::evaluateIntegral( double independentVariable )
{
    return std::exp( scaleFactor_ * exponentExponentialFunction_ * independentVariable )
           / ( std::pow( scaleFactor_ * exponentExponentialFunction_, 2.0 ) + std::pow( frequencySineFunction_, 2.0 ) )
           * ( - frequencySineFunction_ * std::cos( frequencySineFunction_ * independentVariable )
               + scaleFactor_ * exponentExponentialFunction_ * std::sin( frequencySineFunction_ * independentVariable ) );
}


//! Function, derivative, and integral of exponential times cosine function.

double ExponentialCosineFunction::evaluateFunction( double independentVariable )
{
    return std::exp( exponentExponentialFunction_ * independentVariable )
           * std::cos( frequencyCosineFunction_ * independentVariable );
}

double ExponentialCosineFunction::evaluateDerivative( double independentVariable )
{
    return std::exp( exponentExponentialFunction_ * independentVariable )
           * ( exponentExponentialFunction_ * std::cos( frequencyCosineFunction_ * independentVariable )
               - frequencyCosineFunction_ * std::sin( frequencyCosineFunction_ * independentVariable ) );
}

double ExponentialCosineFunction::evaluateIntegral( double independentVariable )
{
    return std::exp( exponentExponentialFunction_ * independentVariable )
            / ( std::pow( exponentExponentialFunction_, 2.0 ) + std::pow( frequencyCosineFunction_, 2.0 ) )
            * ( frequencyCosineFunction_ * std::sin( frequencyCosineFunction_ * independentVariable )
                + exponentExponentialFunction_ * std::cos( frequencyCosineFunction_ * independentVariable ) );
}


//! Function, derivative, and integral of the scaled exponential times cosine function.

double ScaledExponentialCosineFunction::evaluateFunction( double independentVariable )
{
    return std::exp( scaleFactor_ * exponentExponentialFunction_ * independentVariable )
           * std::cos( frequencyCosineFunction_ * independentVariable );
}

double ScaledExponentialCosineFunction::evaluateDerivative( double independentVariable )
{
    return std::exp( scaleFactor_ * exponentExponentialFunction_ * independentVariable )
            * ( scaleFactor_ * exponentExponentialFunction_ * std::cos( frequencyCosineFunction_ * independentVariable )
                - frequencyCosineFunction_ * std::sin( frequencyCosineFunction_ * independentVariable ) );
}

double ScaledExponentialCosineFunction::evaluateIntegral( double independentVariable )
{
    return std::exp( scaleFactor_ * exponentExponentialFunction_ * independentVariable )
           / ( std::pow( scaleFactor_ * exponentExponentialFunction_, 2.0 ) + std::pow( frequencyCosineFunction_, 2.0 ) )
           * ( frequencyCosineFunction_ * std::sin( frequencyCosineFunction_ * independentVariable )
               + scaleFactor_ * exponentExponentialFunction_ * std::cos( frequencyCosineFunction_ * independentVariable ) );
}


//! Function, derivative, and integral of the power function.

double PowerFunction::evaluateFunction( double independentVariable )
{
    return std::pow( independentVariable , exponent_ );
}

double PowerFunction::evaluateDerivative( double independentVariable )
{
    return exponent_ * std::pow( independentVariable , ( exponent_ - 1.0 ) );
}

double PowerFunction::evaluateIntegral( double independentVariable )
{
    return std::pow( independentVariable , ( exponent_ + 1.0 ) ) / ( exponent_ + 1.0 );
}


//! Function, derivative, and integral of scaled power function.

double ScaledPowerFunction::evaluateFunction( double independentVariable )
{
    return scaleFactor_ * std::pow( independentVariable , exponent_ );
}

double ScaledPowerFunction::evaluateDerivative( double independentVariable )
{
    return scaleFactor_ * exponent_ * std::pow( independentVariable , ( exponent_ - 1.0 ) );
}

double ScaledPowerFunction::evaluateIntegral( double independentVariable )
{
    return scaleFactor_ * std::pow( independentVariable , ( exponent_ + 1.0 ) ) / ( exponent_ + 1.0 );
}


//! Function, derivative, and integral of the power times sine function.

double PowerSineFunction::evaluateFunction( double independentVariable )
{
    return std::pow( independentVariable, exponentPowerFunction_ )
           * std::sin( frequencySineFunction_ * independentVariable );
}

double PowerSineFunction::evaluateDerivative( double independentVariable )
{
    return exponentPowerFunction_ * std::pow( independentVariable, exponentPowerFunction_ - 1.0 )
           * std::sin( frequencySineFunction_ * independentVariable )
           + std::pow( independentVariable, exponentPowerFunction_ )
           * frequencySineFunction_ * std::cos( frequencySineFunction_ * independentVariable ) ;
}

double PowerSineFunction::evaluateIntegral( double independentVariable )
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

double ScaledPowerSineFunction::evaluateFunction( double independentVariable )
{
    return scaleFactor_ * std::pow( independentVariable, exponentPowerFunction_ )
           * std::sin( frequencySineFunction_ * independentVariable );
}

double ScaledPowerSineFunction::evaluateDerivative( double independentVariable )
{
    return scaleFactor_ *
           ( exponentPowerFunction_ * std::pow( independentVariable, exponentPowerFunction_ - 1.0 )
               * std::sin( frequencySineFunction_ * independentVariable )
             + std::pow( independentVariable, exponentPowerFunction_ )
               * frequencySineFunction_ * std::cos( frequencySineFunction_ * independentVariable ) ) ;
}

double ScaledPowerSineFunction::evaluateIntegral( double independentVariable )
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

double PowerCosineFunction::evaluateFunction( double independentVariable )
{
    return std::pow( independentVariable, exponentPowerFunction_ )
           * std::cos( frequencyCosineFunction_ * independentVariable );
}

double PowerCosineFunction::evaluateDerivative( double independentVariable )
{
    return exponentPowerFunction_ * std::pow( independentVariable, exponentPowerFunction_ - 1.0 )
           * std::cos( frequencyCosineFunction_ * independentVariable )
           - std::pow( independentVariable, exponentPowerFunction_ )
           * frequencyCosineFunction_ * std::sin( frequencyCosineFunction_ * independentVariable ) ;
}

double PowerCosineFunction::evaluateIntegral( double independentVariable )
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

double ScaledPowerCosineFunction::evaluateFunction( double independentVariable )
{
    return scaleFactor_ * std::pow( independentVariable, exponentPowerFunction_ )
           * std::cos( frequencyCosineFunction_ * independentVariable );
}

double ScaledPowerCosineFunction::evaluateDerivative( double independentVariable )
{
    return scaleFactor_ *
           ( exponentPowerFunction_ * std::pow( independentVariable, exponentPowerFunction_ - 1.0 )
               * std::cos( frequencyCosineFunction_ * independentVariable )
             - std::pow( independentVariable, exponentPowerFunction_ )
               * frequencyCosineFunction_ * std::sin( frequencyCosineFunction_ * independentVariable ) );
}

double ScaledPowerCosineFunction::evaluateIntegral( double independentVariable )
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
