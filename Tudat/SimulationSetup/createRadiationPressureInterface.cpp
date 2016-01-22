#include <boost/bind.hpp>

#include "Tudat/SimulationSetup/createRadiationPressureInterface.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/sphericalBodyShapeModel.h"
#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"

namespace tudat
{

namespace simulation_setup
{

void getOccultingBodiesInformation(
        const NamedBodyMap& bodyMap, const std::vector< std::string >& occultingBodies,
        std::vector< boost::function< Eigen::Vector3d( ) > >& occultingBodyPositions,
        std::vector< double >& occultingBodyRadii )
{
    for( unsigned int i = 0; i < occultingBodies.size( ); i++ )
    {
        if( bodyMap.count( occultingBodies[ i ] ) == 0 )
        {
            std::cerr<<"Error, could not find body "<<occultingBodies[ i ]<<" in body map when making occulting body settings"<<std::endl;
        }
        else
        {
            boost::shared_ptr< basic_astrodynamics::BodyShapeModel > shapeModel =
                    bodyMap.at( occultingBodies.at( i   ) )->getShapeModel( );
            occultingBodyPositions.push_back( boost::bind( &Body::getPosition, bodyMap.at( occultingBodies[ i ] ) ) );
            occultingBodyRadii.push_back( shapeModel->getAverageRadius( ) );
        }
    }
}

boost::shared_ptr< electro_magnetism::RadiationPressureInterface > createRadiationPressureInterface(
        const boost::shared_ptr< RadiationPressureInterfaceSettings > radiationPressureInterfaceSettings,
        const std::string& bodyName, const NamedBodyMap& bodyMap )
{
    boost::shared_ptr< electro_magnetism::RadiationPressureInterface > radiationPressureInterface;

    switch( radiationPressureInterfaceSettings->getRadiationPressureType( ) )
    {
    case cannon_ball:
    {
        boost::shared_ptr< CannonBallRadiationPressureInterfaceSettings > cannonBallSettings =
                boost::dynamic_pointer_cast< CannonBallRadiationPressureInterfaceSettings >(
                    radiationPressureInterfaceSettings );
        if( cannonBallSettings == NULL )
        {
            std::cerr<<"Error when making cannon ball radiation interface, type does not match object"<<std::endl;
        }
        boost::shared_ptr< Body > sourceBody =
                bodyMap.at( radiationPressureInterfaceSettings->getSourceBody( ) );

        if( sourceBody == NULL )
        {
            std::cerr<<"Error when making cannon ball radiation interface, source "<<
                       radiationPressureInterfaceSettings->getSourceBody( )<<" is not a celestial body"<<std::endl;
        }

        std::vector< std::string > occultingBodies = cannonBallSettings->getOccultingBodies( );
        std::vector< boost::function< Eigen::Vector3d( ) > > occultingBodyPositions;
        std::vector< double > occultingBodyRadii;

        getOccultingBodiesInformation( bodyMap, occultingBodies, occultingBodyPositions, occultingBodyRadii );

        double sourceRadius;

        if( occultingBodyPositions.size( ) > 0 )
        {
            boost::shared_ptr< basic_astrodynamics::BodyShapeModel > sourceShapeModel = sourceBody->getShapeModel( );

            if( sourceShapeModel == NULL )
            {
                std::cerr<<"Error when making occulted body, source body "<<radiationPressureInterfaceSettings->getSourceBody( )<<
                           " does not have a shape"<<std::endl;
            }
            else
            {
                sourceRadius = sourceShapeModel->getAverageRadius( );
            }
        }
        else
        {
            sourceRadius = 0.0;
        }

        boost::function< double( ) > radiatedPowerFunction;
        if( defaultRadiatedPowerValues.count(
                    radiationPressureInterfaceSettings->getSourceBody( ) ) == 0 )
        {
            throw( "" );
        }
        else
        {
            radiatedPowerFunction = boost::lambda::constant(
                        defaultRadiatedPowerValues.at(
                            radiationPressureInterfaceSettings->getSourceBody( ) ) );
        }

        radiationPressureInterface = boost::make_shared< electro_magnetism::RadiationPressureInterface >(
                    radiatedPowerFunction,
                    boost::bind( &Body::getPosition, sourceBody ),
                    boost::bind( &Body::getPosition, bodyMap.at( bodyName ) ),
                    cannonBallSettings->getRadiationPressureCoefficient( ),
                    cannonBallSettings->getArea( ), occultingBodyPositions, occultingBodyRadii,
                    sourceRadius );

    }
    default:
        throw( "" );
    }

    return radiationPressureInterface;
}

}

}
