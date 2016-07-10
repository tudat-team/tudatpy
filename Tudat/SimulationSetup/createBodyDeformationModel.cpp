#include <boost/lambda/lambda.hpp>

#include "Tudat/Mathematics/Interpolators/interpolator.h"

#include "Tudat/Astrodynamics/GroundStations/basicTidalBodyDeformation.h"
#include "Tudat/SimulationSetup/createBodyDeformationModel.h"

namespace tudat
{

namespace simulation_setup
{


boost::shared_ptr< site_displacements::BodyDeformationModel > createBodyDeformationModel(
        boost::shared_ptr< BodyDeformationSettings > bodyDeformationSettings,
        const std::string body,
        const NamedBodyMap bodyMap )
{
    using namespace tudat::site_displacements;
    using namespace tudat::gravitation;

    boost::shared_ptr< BodyDeformationModel > bodyDeformationModel;

    switch( bodyDeformationSettings->getBodyDeformationType( ) )
    {
    case basic_solid_body:
    {
        boost::shared_ptr< BasicSolidBodyDeformationSettings > basicSolidBodyDeformationSettings =
                boost::dynamic_pointer_cast< BasicSolidBodyDeformationSettings >( bodyDeformationSettings );
        if( basicSolidBodyDeformationSettings == NULL )
        {
            std::cerr<<"Error, expected basic solid body deformation settings for "<<body<<std::endl;
        }
        else
        {
            std::vector< std::string > deformingBodies = basicSolidBodyDeformationSettings->getDeformingBodies( );
            std::vector< boost::function< basic_mathematics::Vector6d( const double ) > > deformingBodyEphemerides;
            std::vector< boost::function< double( ) > > gravitionalParametersOfDeformingBodies;

            for( unsigned int i = 0; i < deformingBodies.size( ); i++ )
            {
                if( bodyMap.count( deformingBodies[ i ] ) == 0 )
                {
                    std::cerr<<"Error when making basic solid body deformation model, deforming body not found: "<<deformingBodies[ i ]<<std::endl;
                }

                deformingBodyEphemerides.push_back(
                            boost::bind( &Body::getStateInBaseFrameFromEphemeris, bodyMap.at( deformingBodies[ i ] ), _1 ) );

                boost::shared_ptr< gravitation::GravityFieldModel > gravityFieldModel =
                        bodyMap.at( deformingBodies[ i ] )->getGravityFieldModel( );
                if( gravityFieldModel == NULL )
                {
                    std::cerr<<"Error, no gravity field model of "<<deformingBodies[ i ]<<" found when making basic body deformation of "<<body<<std::endl;
                }
                gravitionalParametersOfDeformingBodies.push_back(
                            boost::bind( &GravityFieldModel::getGravitationalParameter, gravityFieldModel ) );
            }

            boost::shared_ptr< gravitation::GravityFieldModel > gravityFieldModel =
                    bodyMap.at( body )->getGravityFieldModel( );
            if( gravityFieldModel == NULL )
            {
                std::cerr<<"Error, no gravity field model of "<<body<<" found when making basic body deformation of "<<body<<std::endl;
            }

            boost::function< double( ) > gravitionalParameterOfDeformedBody =
                    boost::bind( &GravityFieldModel::getGravitationalParameter, gravityFieldModel );


            bodyDeformationModel = boost::make_shared< BasicTidalBodyDeformation >(
                        boost::bind( &Body::getStateInBaseFrameFromEphemeris, bodyMap.at( body ), _1 ),
                        deformingBodyEphemerides,
                        boost::bind( &ephemerides::RotationalEphemeris::getRotationToTargetFrame,
                                     bodyMap.at( body )->getRotationalEphemeris( ), _1, basic_astrodynamics::JULIAN_DAY_ON_J2000 ),
                        basicSolidBodyDeformationSettings->getMaximumDeformationOrders( ),
                        gravitionalParameterOfDeformedBody,
                        gravitionalParametersOfDeformingBodies,
                        boost::lambda::constant( basicSolidBodyDeformationSettings->getBodyReferenceRadius( ) ),
                        basicSolidBodyDeformationSettings->getDisplacementLoveNumbers( ),
                        basicSolidBodyDeformationSettings->getDisplacementShidaNumbers( ) );

        }
        break;
    }
    default:
        std::cerr<<"Error, body deformation type not recognized"<<std::endl;
    }

    return bodyDeformationModel;
}


}

}

