#include <iostream>

#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/sphericalBodyShapeModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/oblateSpheroidBodyShapeModel.h"
#include "Tudat/SimulationSetup/createBodyShapeModel.h"

namespace tudat
{

namespace simulation_setup
{

boost::shared_ptr< basic_astrodynamics::BodyShapeModel > createBodyShapeModel(
        const boost::shared_ptr< BodyShapeSettings > shapeSettings,
        const std::string& body )
{
    using namespace tudat::basic_astrodynamics;

    boost::shared_ptr< BodyShapeModel > shapeModel;

    switch( shapeSettings->getBodyShapeType( ) )
    {
    case spherical:
    {
        boost::shared_ptr< SphericalBodyShapeSettings > sphericalShapeSettings =
                boost::dynamic_pointer_cast< SphericalBodyShapeSettings >( shapeSettings );

        if( sphericalShapeSettings == NULL )
        {
            std::cerr<<"Error, expected spherical shape settings for body "<<body<<std::endl;
        }
        else
        {
            shapeModel = boost::make_shared< SphericalBodyShapeModel >( sphericalShapeSettings->getRadius( ) );
        }
        break;
    }
    case oblate_spheroid:
    {
        boost::shared_ptr< OblateSphericalBodyShapeSettings > oblateSpheroidShapeSettings =
                boost::dynamic_pointer_cast< OblateSphericalBodyShapeSettings >( shapeSettings );

        if( oblateSpheroidShapeSettings == NULL )
        {
            std::cerr<<"Error, expected oblate spherical shape settings for body "<<body<<std::endl;
        }
        else
        {
            shapeModel = boost::make_shared< OblateSpheroidBodyShapeModel >(
                        oblateSpheroidShapeSettings->getEquatorialRadius( ),
                        oblateSpheroidShapeSettings->getFlattening( ) );
        }
        break;
    }
    case spherical_spice:
    {
        shapeModel = boost::make_shared< SphericalBodyShapeModel >(
                    spice_interface::getAverageRadius( body ) );
        break;
    }
    default:
        std::cerr<<"Error, did not recognize body shape settings for "<<body<<std::endl;

    }
    return shapeModel;
}

}

}
