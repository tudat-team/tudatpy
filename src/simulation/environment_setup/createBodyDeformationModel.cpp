#include "tudat/astro/ground_stations/basicTidalBodyDeformation.h"
#include "tudat/simulation/environment_setup/createBodyDeformationModel.h"

namespace tudat
{

namespace simulation_setup
{


std::shared_ptr< basic_astrodynamics::BodyDeformationModel > createBodyDeformationModel(
        const std::shared_ptr< BodyDeformationSettings > bodyDeformationSettings,
        const std::string body,
        const SystemOfBodies& bodyMap )
{
    using namespace tudat::basic_astrodynamics;
    using namespace tudat::gravitation;

    std::shared_ptr< BodyDeformationModel > bodyDeformationModel;

    switch( bodyDeformationSettings->getBodyDeformationType( ) )
    {
    case basic_solid_body:
    {
        std::shared_ptr< BasicSolidBodyDeformationSettings > basicSolidBodyDeformationSettings =
                std::dynamic_pointer_cast< BasicSolidBodyDeformationSettings >( bodyDeformationSettings );

        if( basicSolidBodyDeformationSettings == nullptr )
        {
            throw std::runtime_error( "Error when creating body deformation model, expected basic solid body settings for " + body );
        }
        else
        {
            std::vector< std::string > deformingBodies = basicSolidBodyDeformationSettings->getDeformingBodies( );
            std::vector< std::function< Eigen::Vector6d( const double ) > > deformingBodyEphemerides;
            std::vector< std::function< double( ) > > gravitionalParametersOfDeformingBodies;

            for( unsigned int i = 0; i < deformingBodies.size( ); i++ )
            {
                if( bodyMap.count( deformingBodies.at( i ) ) == 0 )
                {
                    throw std::runtime_error( "Error when making basic solid body deformation model, deforming body not found: " +
                                              deformingBodies.at( i ) );
                }

                deformingBodyEphemerides.push_back(
                            std::bind( &Body::getStateInBaseFrameFromEphemeris< double, double >, bodyMap.at( deformingBodies.at( i ) ), std::placeholders::_1 ) );

                std::shared_ptr< gravitation::GravityFieldModel > gravityFieldModel =
                        bodyMap.at( deformingBodies.at( i ) )->getGravityFieldModel( );
                if( gravityFieldModel == nullptr )
                {
                    throw std::runtime_error( "Error, no gravity field model of " + deformingBodies.at( i ) +
                                              " found when making basic body deformation of " + body );
                }
                gravitionalParametersOfDeformingBodies.push_back(
                            std::bind( &GravityFieldModel::getGravitationalParameter, gravityFieldModel ) );
            }

            std::shared_ptr< gravitation::GravityFieldModel > gravityFieldModel = bodyMap.at( body )->getGravityFieldModel( );
            if( gravityFieldModel == NULL )
            {
                throw std::runtime_error("Error, no gravity field model of " + body +
                                         " found when making basic body deformation of " + body );
            }

            std::function< double( ) > gravitionalParameterOfDeformedBody =
                    std::bind( &GravityFieldModel::getGravitationalParameter, gravityFieldModel );


            double deformationReferenceRadius = basicSolidBodyDeformationSettings->getBodyReferenceRadius( );
            if( deformationReferenceRadius != deformationReferenceRadius )
            {
                std::shared_ptr< BodyShapeModel > bodyShapeModel =
                        bodyMap.at( body )->getShapeModel( );
                if( bodyShapeModel == nullptr )
                {
                    throw std::runtime_error("Error, when making basic body deformation of " + body + ", no reference radius, and no shape model specified" );
                }
                else
                {
                    deformationReferenceRadius = bodyShapeModel->getAverageRadius( );
                }
            }

            bodyDeformationModel = std::make_shared< basic_astrodynamics::BasicTidalBodyDeformation >(
                        std::bind( &Body::getStateInBaseFrameFromEphemeris< double, double >, bodyMap.at( body ),
                                   std::placeholders::_1 ),
                        deformingBodyEphemerides,
                        std::bind( &ephemerides::RotationalEphemeris::getRotationToTargetFrame, bodyMap.at( body )->getRotationalEphemeris( ),
                                   std::placeholders::_1 ),
                        gravitionalParameterOfDeformedBody,
                        gravitionalParametersOfDeformingBodies,
                        deformationReferenceRadius,
                        basicSolidBodyDeformationSettings->getDisplacementLoveNumbers( ) );

        }
        break;
    }
    default:
        throw std::runtime_error( "Error, did not recognize body deformation settings type " +
                                  std::to_string( bodyDeformationSettings->getBodyDeformationType( ) ) +
                                  " of body " + body );
    }

    return bodyDeformationModel;
}

}

}
