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


#include <cmath>
#include <stdexcept>

#include "tudat/astro/low_thrust/shape_based/baseFunctionsSphericalShaping.h"


namespace tudat
{
namespace shape_based_methods
{

std::shared_ptr< BaseFunctionSphericalShaping > createBaseFunctionSphericalShaping(
        const baseFunctionSphericalShapingType baseFunctionType )
{
    // Declare return object.
    std::shared_ptr< BaseFunctionSphericalShaping > baseFunctionSphericalShaping;

    // Check which type of base function is to be created.
    switch( baseFunctionType )
    {
    case constantSphericalShaping:
    {
        baseFunctionSphericalShaping = std::make_shared< ConstantFunctionSphericalShaping >( );
        break;
    }
    case linearSphericalShaping:
    {
        baseFunctionSphericalShaping = std::make_shared< LinearFunctionSphericalShaping >( );
        break;
    }
    case squaredSphericalShaping:
    {
        baseFunctionSphericalShaping = std::make_shared< SquaredFunctionSphericalShaping >( );
        break;
    }
    case cosineSphericalShaping:
    {
        baseFunctionSphericalShaping = std::make_shared< CosineFunctionSphericalShaping >( );
        break;
    }
    case powerCosineSphericalShaping:
    {
        baseFunctionSphericalShaping = std::make_shared< PowerCosineFunctionSphericalShaping >( );
        break;
    }
    case sineSphericalShaping:
    {
        baseFunctionSphericalShaping = std::make_shared< SineFunctionSphericalShaping >( );
        break;
    }
    case powerSineSphericalShaping:
    {
        baseFunctionSphericalShaping = std::make_shared< PowerSineFunctionSphericalShaping >( );
        break;
    }
    default:
    {
        throw std::runtime_error(
                    "Error, did not recognize base function type for spherical shaping" );
    }

    }
    return baseFunctionSphericalShaping;
}

} // namespace shape_based_methods
} // namespace tudat
