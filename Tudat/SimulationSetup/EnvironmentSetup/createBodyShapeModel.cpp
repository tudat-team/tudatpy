/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#if USE_CSPICE
#include "Tudat/External/SpiceInterface/spiceInterface.h"
#endif
#include "Tudat/Astrodynamics/BasicAstrodynamics/sphericalBodyShapeModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/oblateSpheroidBodyShapeModel.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createBodyShapeModel.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to create a body shape model.
std::shared_ptr< basic_astrodynamics::BodyShapeModel > createBodyShapeModel(
        const std::shared_ptr< BodyShapeSettings > shapeSettings,
        const std::string& body )
{
    using namespace tudat::basic_astrodynamics;

    std::shared_ptr< BodyShapeModel > shapeModel;

    // Check body shape type
    switch( shapeSettings->getBodyShapeType( ) )
    {
    case spherical:
    {
        // Check input consistency
        std::shared_ptr< SphericalBodyShapeSettings > sphericalShapeSettings =
                std::dynamic_pointer_cast< SphericalBodyShapeSettings >( shapeSettings );
        if( sphericalShapeSettings == nullptr )
        {
            throw std::runtime_error( "Error, expected spherical shape settings for body " + body );
        }
        else
        {
            // Creat spherical shape model
            shapeModel = std::make_shared< SphericalBodyShapeModel >(
                sphericalShapeSettings->getRadius( ) );
        }
        break;
    }
    case oblate_spheroid:
    {
        // Check input consistency
        std::shared_ptr< OblateSphericalBodyShapeSettings > oblateSpheroidShapeSettings =
                std::dynamic_pointer_cast< OblateSphericalBodyShapeSettings >( shapeSettings );
        if( oblateSpheroidShapeSettings == nullptr )
        {
           throw std::runtime_error( "Error, expected oblate spherical shape settings for body " + body );
        }
        else
        {
            // Creat oblate spheroid shape model
            shapeModel = std::make_shared< OblateSpheroidBodyShapeModel >(
                        oblateSpheroidShapeSettings->getEquatorialRadius( ),
                        oblateSpheroidShapeSettings->getFlattening( ) );
        }
        break;
    }
#if USE_CSPICE
    case spherical_spice:
    {
        // Retrieve radius from Spice and create spherical shape model.
        shapeModel = std::make_shared< SphericalBodyShapeModel >(
                    spice_interface::getAverageRadius( body ) );
        break;
    }
#endif
    default:
       throw std::runtime_error( "Error, did not recognize body shape settings for " + body );

    }
    return shapeModel;
}

} // namespace simulation_setup

} // namespace tudat
