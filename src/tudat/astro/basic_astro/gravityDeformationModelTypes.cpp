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

#include "tudat/astro/basic_astro/gravityDeformationModelTypes.h"

namespace tudat
{

namespace basic_astrodynamics
{

GravityDeformationType getGravityDeformationModelType(
        std::shared_ptr< basic_astrodynamics::GravityDeformationModel > gravityDeformationModel )
{
    if( std::dynamic_pointer_cast< basic_astrodynamics::MaxwellGravityDeformationModel >( gravityDeformationModel ) != nullptr )
    {
        return maxwell_deformation;
    }

    throw std::runtime_error( "Error, could not identify gravity deformation type" );
}

std::string getGravityDeformationModelName( const GravityDeformationType gravityDeformationType )
{
    switch( gravityDeformationType )
    {
        case maxwell_deformation:
            return "Maxwell deformation ";
        default:
            throw std::runtime_error( "Error, could not identify gravity deformation type" );
    }
}

}  // namespace basic_astrodynamics

}  // namespace tudat
