/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include <stdexcept>
#include <cmath>
#include "tudat/astro/low_thrust/shape_based/createBaseFunctionHodographicShaping.h"


namespace tudat
{
namespace shape_based_methods
{

//! Function to create a base function for hodographic shaping.
std::shared_ptr< BaseFunctionHodographicShaping > createBaseFunctionHodographicShaping(
        const baseFunctionHodographicShapingType baseFunctionType,
        const std::shared_ptr< BaseFunctionHodographicShapingSettings > baseFunctionSettings )
{

    // Declare return object.
    std::shared_ptr< BaseFunctionHodographicShaping > baseFunctionHodographicShaping;

    // Check which type of base function is to be created.
    switch( baseFunctionType )
    {
    case constant:
    {
        baseFunctionHodographicShaping = std::make_shared< ConstantFunctionHodographicShaping >( );
        break;
    }
    case sine:
    {
        std::shared_ptr< TrigonometricFunctionHodographicShapingSettings > trigonometricFunctionSettings =
                std::dynamic_pointer_cast< TrigonometricFunctionHodographicShapingSettings >( baseFunctionSettings );

        if ( trigonometricFunctionSettings != nullptr )
        {
            baseFunctionHodographicShaping = std::make_shared< SineFunctionHodographicShaping >( trigonometricFunctionSettings->frequency_ );
        }
        else
        {
            throw std::runtime_error( "Error when creating sine base function for hodographic shaping, "
                                      "unconsistent settings were provided." );
        }
        break;
    }
    case cosine:
    {
        std::shared_ptr< TrigonometricFunctionHodographicShapingSettings > trigonometricFunctionSettings =
                std::dynamic_pointer_cast< TrigonometricFunctionHodographicShapingSettings >( baseFunctionSettings );

        if ( trigonometricFunctionSettings != nullptr )
        {
            baseFunctionHodographicShaping = std::make_shared< CosineFunctionHodographicShaping >( trigonometricFunctionSettings->frequency_ );
        }
        else
        {
            throw std::runtime_error( "Error when creating cosine base function for hodographic shapign, "
                                      "unconsistent settings were provided." );
        }
        break;
    }
    case exponential:
    {
        std::shared_ptr< ExponentialFunctionHodographicShapingSettings > exponentialFunctionSettings =
                std::dynamic_pointer_cast< ExponentialFunctionHodographicShapingSettings >( baseFunctionSettings );

        if ( exponentialFunctionSettings != nullptr )
        {
            baseFunctionHodographicShaping = std::make_shared< ExponentialFunctionHodographicShaping >( exponentialFunctionSettings->exponent_ );
        }
        else
        {
            throw std::runtime_error( "Error when creating exponential base function for hodographic shaping,"
                                      " unconsistent settings were provided." );
        }
        break;
    }
    case scaledExponential:
    {
        std::shared_ptr< ExponentialFunctionHodographicShapingSettings > exponentialFunctionSettings =
                std::dynamic_pointer_cast< ExponentialFunctionHodographicShapingSettings >( baseFunctionSettings );

        if ( exponentialFunctionSettings != nullptr )
        {
            baseFunctionHodographicShaping = std::make_shared< ScaledExponentialFunctionHodographicShaping >(
                        exponentialFunctionSettings->exponent_,
                        exponentialFunctionSettings->scaleFactor_ );
        }
        else
        {
            throw std::runtime_error( "Error when creating scaled exponential base function for hodographic shaping, "
                                      "unconsistent settings were provided." );
        }
        break;
    }
    case exponentialSine:
    {
        std::shared_ptr< ExponentialTimesTrigonometricFunctionHodographicShapingSettings > exponentialTimeSineFunctionSettings =
                std::dynamic_pointer_cast< ExponentialTimesTrigonometricFunctionHodographicShapingSettings >( baseFunctionSettings );

        if ( exponentialTimeSineFunctionSettings != nullptr )
        {
            baseFunctionHodographicShaping = std::make_shared< ExponentialSineFunctionHodographicShaping >(
                        exponentialTimeSineFunctionSettings->exponent_,
                        exponentialTimeSineFunctionSettings->frequency_ );
        }
        else
        {
            throw std::runtime_error( "Error when creating exponential time sine base function for hodographic shaping, "
                                      "unconsistent settings were provided." );
        }
        break;
    }
    case scaledExponentialSine:
    {
        std::shared_ptr< ExponentialTimesTrigonometricFunctionHodographicShapingSettings > exponentialTimesSineFunctionSettings =
                std::dynamic_pointer_cast< ExponentialTimesTrigonometricFunctionHodographicShapingSettings >( baseFunctionSettings );

        if ( exponentialTimesSineFunctionSettings != nullptr )
        {
            baseFunctionHodographicShaping = std::make_shared< ScaledExponentialSineFunctionHodographicShaping >(
                        exponentialTimesSineFunctionSettings->exponent_,
                        exponentialTimesSineFunctionSettings->frequency_,
                        exponentialTimesSineFunctionSettings->scaleFactor_ );
        }
        else
        {
            throw std::runtime_error( "Error when creating scaled exponential time sine base function for hodographic shaping, "
                                      "unconsistent settings were provided." );
        }
        break;
    }
    case exponentialCosine:
    {
        std::shared_ptr< ExponentialTimesTrigonometricFunctionHodographicShapingSettings > exponentialTimesCosineFunctionSettings =
                std::dynamic_pointer_cast< ExponentialTimesTrigonometricFunctionHodographicShapingSettings >( baseFunctionSettings );

        if ( exponentialTimesCosineFunctionSettings != nullptr )
        {
            baseFunctionHodographicShaping = std::make_shared< ExponentialCosineFunctionHodographicShaping >(
                        exponentialTimesCosineFunctionSettings->exponent_,
                        exponentialTimesCosineFunctionSettings->frequency_ );
        }
        else
        {
            throw std::runtime_error( "Error when creating exponential time cosine base function for hodographic shaping, "
                                      "unconsistent settings were provided." );
        }
        break;
    }
    case scaledExponentialCosine:
    {
        std::shared_ptr< ExponentialTimesTrigonometricFunctionHodographicShapingSettings > exponentialTimesCosineFunctionSettings =
                std::dynamic_pointer_cast< ExponentialTimesTrigonometricFunctionHodographicShapingSettings >( baseFunctionSettings );

        if ( exponentialTimesCosineFunctionSettings != nullptr )
        {
            baseFunctionHodographicShaping = std::make_shared< ScaledExponentialCosineFunctionHodographicShaping >(
                        exponentialTimesCosineFunctionSettings->exponent_,
                        exponentialTimesCosineFunctionSettings->frequency_,
                        exponentialTimesCosineFunctionSettings->scaleFactor_ );
        }
        else
        {
            throw std::runtime_error( "Error when creating scaled exponential time cosine base function for hodographic shaping, "
                                      "unconsistent settings were provided." );
        }
        break;
    }
    case power:
    {
        std::shared_ptr< PowerFunctionHodographicShapingSettings > powerFunctionSettings =
                std::dynamic_pointer_cast< PowerFunctionHodographicShapingSettings >( baseFunctionSettings );

        if ( powerFunctionSettings != nullptr )
        {
            baseFunctionHodographicShaping = std::make_shared< PowerFunctionHodographicShaping >( powerFunctionSettings->exponent_ );
        }
        else
        {
            throw std::runtime_error( "Error when creating power base function for hodographic shaping, "
                                      "unconsistent settings were provided." );
        }
        break;
    }
    case scaledPower:
    {
        std::shared_ptr< PowerFunctionHodographicShapingSettings > powerFunctionSettings =
                std::dynamic_pointer_cast< PowerFunctionHodographicShapingSettings >( baseFunctionSettings );

        if ( powerFunctionSettings != nullptr )
        {
            baseFunctionHodographicShaping = std::make_shared< ScaledPowerFunctionHodographicShaping >(
                        powerFunctionSettings->exponent_,
                        std::pow( powerFunctionSettings->scaleFactor_, powerFunctionSettings->exponent_ ) );
        }
        else
        {
            throw std::runtime_error( "Error when creating scaled power base function for hodographic shaping, "
                                      "unconsistent settings were provided." );
        }
        break;
    }
    case powerSine:
    {
        std::shared_ptr< PowerTimesTrigonometricFunctionHodographicShapingSettings > powerTimesSineFunctionSettings =
                std::dynamic_pointer_cast< PowerTimesTrigonometricFunctionHodographicShapingSettings >( baseFunctionSettings );

        if ( powerTimesSineFunctionSettings != nullptr )
        {
            baseFunctionHodographicShaping = std::make_shared< PowerSineFunctionHodographicShaping >( powerTimesSineFunctionSettings->exponent_,
                                                                                    powerTimesSineFunctionSettings->frequency_ );
        }
        else
        {
            throw std::runtime_error( "Error when creating power times sine base function for hodographic shaping, "
                                      "unconsistent settings were provided." );
        }
        break;
    }
    case scaledPowerSine:
    {
        std::shared_ptr< PowerTimesTrigonometricFunctionHodographicShapingSettings > powerTimesSineFunctionSettings =
                std::dynamic_pointer_cast< PowerTimesTrigonometricFunctionHodographicShapingSettings >( baseFunctionSettings );

        if ( powerTimesSineFunctionSettings != nullptr )
        {
            baseFunctionHodographicShaping = std::make_shared< ScaledPowerSineFunctionHodographicShaping >(
                        powerTimesSineFunctionSettings->exponent_,
                        powerTimesSineFunctionSettings->frequency_,
                        std::pow( powerTimesSineFunctionSettings->scaleFactor_, powerTimesSineFunctionSettings->exponent_ ) );
        }
        else
        {
            throw std::runtime_error( "Error when creating scaled power times sine base function for hodographic shaping, "
                                      "unconsistent settings were provided." );
        }
        break;
    }
    case powerCosine:
    {
        std::shared_ptr< PowerTimesTrigonometricFunctionHodographicShapingSettings > powerTimesCosineFunctionSettings =
                std::dynamic_pointer_cast< PowerTimesTrigonometricFunctionHodographicShapingSettings >( baseFunctionSettings );

        if ( powerTimesCosineFunctionSettings != nullptr )
        {
            baseFunctionHodographicShaping = std::make_shared< PowerCosineFunctionHodographicShaping >(
                        powerTimesCosineFunctionSettings->exponent_,
                        powerTimesCosineFunctionSettings->frequency_ );
        }
        else
        {
            throw std::runtime_error( "Error when creating power times cosine base function for hodographic shaping, "
                                      "unconsistent settings were provided." );
        }
        break;
    }
    case scaledPowerCosine:
    {
        std::shared_ptr< PowerTimesTrigonometricFunctionHodographicShapingSettings > powerTimesCosineSettings =
                std::dynamic_pointer_cast< PowerTimesTrigonometricFunctionHodographicShapingSettings >( baseFunctionSettings );

        if ( powerTimesCosineSettings != nullptr )
        {
            baseFunctionHodographicShaping = std::make_shared< ScaledPowerCosineFunctionHodographicShaping >(
                        powerTimesCosineSettings->exponent_,
                        powerTimesCosineSettings->frequency_,
                        std::pow( powerTimesCosineSettings->scaleFactor_, powerTimesCosineSettings->exponent_ ) );
        }
        else
        {
            throw std::runtime_error( "Error when creating scaled power times cosine base function for hodographic shaping, "
                                      "unconsistent settings were provided." );
        }
        break;
    }
    default:
    {
        throw std::runtime_error(
                    "Error, did not recognize base function type for hodographic shaping." );
    }

    }

    return baseFunctionHodographicShaping;

};

} // namespace shape_based_methods
} // namespace tudat
