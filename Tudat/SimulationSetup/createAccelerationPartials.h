#ifndef TUDAT_CREATEACCELERATIONPARTIALS_H
#define TUDAT_CREATEACCELERATIONPARTIALS_H

#include <boost/make_shared.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"

#include "Tudat/SimulationSetup/body.h"
#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/accelerationPartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/centralGravityAccelerationPartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/radiationPressureAccelerationPartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/thirdBodyGravityPartial.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModelTypes.h"

namespace tudat
{

namespace orbit_determination
{

namespace partial_derivatives
{

template< typename InitialStateParameterType = double >
boost::shared_ptr< AccelerationPartial > createAnalyticalAccelerationPartial(
        boost::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > accelerationModel,
        const std::pair< std::string, boost::shared_ptr< simulation_setup::Body > > acceleratedBody,
        const std::pair< std::string, boost::shared_ptr< simulation_setup::Body > > acceleratingBody,
        const simulation_setup::NamedBodyMap bodyMap,
        const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > > parametersToEstimate =
        boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > >( ),
        const std::map< std::string, std::map< std::string, std::vector< boost::shared_ptr< AccelerationPartial > > > >& accelerationPartialsMap =
        std::map< std::string, std::map< std::string, std::vector< boost::shared_ptr< AccelerationPartial > > > >( ) )
{

    using namespace basic_astrodynamics;
    using namespace electro_magnetism;
    using namespace aerodynamics;

    AvailableAcceleration accelerationType = getAccelerationModelType( accelerationModel );
    boost::shared_ptr< AccelerationPartial > accelerationPartial;

    // Identify current acceleration model type
    switch( accelerationType )
    {
    case central_gravity:

        // Check if identifier is consistent with type.
        if( boost::dynamic_pointer_cast< CentralGravitationalAccelerationModel3d >( accelerationModel ) == NULL )
        {
            std::cerr<<"Acceleration class type does not match acceleration type enum set when making acceleration partial"<<std::endl;
        }
        else
        {
                // Create partial-calculating object.
                accelerationPartial = boost::make_shared< CentralGravitationPartial >
                        ( boost::dynamic_pointer_cast< CentralGravitationalAccelerationModel3d >( accelerationModel ),
                          acceleratedBody.first, acceleratingBody.first );
        }
        break;

    case third_body_central_gravity:
        // Check if identifier is consistent with type.
        if( boost::dynamic_pointer_cast< ThirdBodyCentralGravityAcceleration >( accelerationModel ) == NULL )
        {
            std::cerr<<"Acceleration class type does not match acceleration type enum set when making acceleration partial"<<std::endl;
        }
        else
        {
            boost::shared_ptr< ThirdBodyCentralGravityAcceleration > thirdBodyAccelerationModel  =
                    boost::dynamic_pointer_cast< ThirdBodyCentralGravityAcceleration >( accelerationModel );

            boost::shared_ptr< CentralGravitationPartial > accelerationPartialForBodyUndergoingAcceleration =
                    boost::dynamic_pointer_cast< CentralGravitationPartial >(
                        createAnalyticalAccelerationPartial(
                            thirdBodyAccelerationModel->getAccelerationModelForBodyUndergoingAcceleration( ),
                            acceleratedBody, acceleratingBody, bodyMap, parametersToEstimate, accelerationPartialsMap ) );
            boost::shared_ptr< CentralGravitationPartial > accelerationPartialForCentralBody =
                    boost::dynamic_pointer_cast< CentralGravitationPartial >(
                        createAnalyticalAccelerationPartial(
                            thirdBodyAccelerationModel->getAccelerationModelForCentralBody( ),
                            std::make_pair( thirdBodyAccelerationModel->getCentralBodyName( ),
                                            bodyMap.at( thirdBodyAccelerationModel->getCentralBodyName( ) ) ),
                            acceleratingBody, bodyMap, parametersToEstimate, accelerationPartialsMap  ) );
            accelerationPartial = boost::make_shared< ThirdBodyGravityPartial< CentralGravitationPartial > >(
                        accelerationPartialForBodyUndergoingAcceleration,
                        accelerationPartialForCentralBody, acceleratedBody.first, acceleratingBody.first,
                        thirdBodyAccelerationModel->getCentralBodyName( ) );

        }
        break;

    case cannon_ball_radiation_pressure:
    {
        // Check if identifier is consistent with type.
        boost::shared_ptr< CannonBallRadiationPressureAcceleration > radiationPressureAcceleration =
                boost::dynamic_pointer_cast< CannonBallRadiationPressureAcceleration >( accelerationModel );
        if( radiationPressureAcceleration == NULL )
        {
            std::cerr<<"Acceleration class type does not match acceleration type enum (cannon ball radiation pressure) set when making acceleration partial"<<std::endl;
        }
        else
        {
            std::map< std::string, boost::shared_ptr< RadiationPressureInterface > > radiationPressureInterfaces =
                    acceleratedBody.second->getRadiationPressureInterfaces( );

            if( radiationPressureInterfaces.count( acceleratingBody.first ) == 0 )
            {
                std::cerr<<"No radiation pressure coefficient interface found when making acceleration partial."<<std::endl;
            }
            else
            {
                boost::shared_ptr< RadiationPressureInterface > radiationPressureInterface =
                        radiationPressureInterfaces.at( acceleratingBody.first );

                // Create partial-calculating object.
                accelerationPartial = boost::make_shared< CannonBallRadiationPressurePartial >
                        ( radiationPressureInterface, radiationPressureAcceleration->getMassFunction( ),
                          acceleratedBody.first, acceleratingBody.first );
            }
        }
        break;
    }
    default:
        std::cerr<<"Acceleration model "<<accelerationType<<" not found when making acceleration partial"<<std::endl;
        break;
    }

    return accelerationPartial;
}

//! This function creates acceleration partial objects from acceleration models and list of bodies' states of which derivatives are needed.
/*!
 *  This function creates acceleration partial objects from acceleration models and list of bodies' states of
 *  which derivatives are needed. The return type is an AccelerationPartialsMap, a standardized type for communicating such lists
 *  of these objects.
 *  \param accelerationMap Map of maps containing list of acceleration models, identifying which acceleration acts on which body.
 *  \param bodyMap List of body objects constituting environment for calculations.
 *  \param parametersToEstimate List of parameters which are to be estimated.
 *  \return List of acceleration-partial-calculating objects in AccelerationPartialsMap type.
 */
template< typename InitialStateParameterType >
StateDerivativePartialsMap createAccelerationPartialsMap(
        const basic_astrodynamics::AccelerationMap& accelerationMap,
        const simulation_setup::NamedBodyMap& bodyMap,
        const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > > parametersToEstimate =
        boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > >( ) )

{
    // Declare return map.
    StateDerivativePartialsMap accelerationPartialsList;
    std::map< std::string, std::map< std::string, std::vector< boost::shared_ptr< AccelerationPartial > > > > accelerationPartialsMap;

    std::vector< boost::shared_ptr< estimatable_parameters::EstimatableParameter<
            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > initialDynamicalParameters =
            parametersToEstimate->getEstimatedInitialStateParameters( );
    accelerationPartialsList.resize( initialDynamicalParameters.size( ) );

    // Iterate over list of bodies of which the partials of the accelerations acting on them are required.
    for( basic_astrodynamics::AccelerationMap::const_iterator accelerationIterator = accelerationMap.begin( );
         accelerationIterator != accelerationMap.end( ); accelerationIterator++ )
    {
        for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
        {
            if( initialDynamicalParameters.at( i )->getParameterName( ).second.first == accelerationIterator->first )
            {
                if( ( initialDynamicalParameters.at( i )->getParameterName( ).first == estimatable_parameters::initial_body_state ) )
                {
                    const std::string acceleratedBody = accelerationIterator->first;

                    boost::shared_ptr< simulation_setup::Body > acceleratedBodyObject = bodyMap.at( acceleratedBody );

                    // Retrieve list of accelerations acting on current body.
                    basic_astrodynamics::SingleBodyAccelerationMap accelerationVector =
                            accelerationMap.at( acceleratedBody );

                    // Declare list of acceleration partials of current body.
                    std::vector< boost::shared_ptr< StateDerivativePartial > > accelerationPartialVector;

                    // Iterate over all acceleration models and generate their partial-calculating objects.
                    for(  basic_astrodynamics::SingleBodyAccelerationMap::iterator
                          innerAccelerationIterator = accelerationVector.begin( );
                          innerAccelerationIterator != accelerationVector.end( ); innerAccelerationIterator++ )
                    {
                        std::string acceleratingBody = innerAccelerationIterator->first;

                        boost::shared_ptr< simulation_setup::Body > acceleratingBodyObject;
                        if( acceleratingBody != "" )
                        {
                            acceleratingBodyObject = bodyMap.at( acceleratingBody );
                        }

                        for( unsigned int j = 0; j < innerAccelerationIterator->second.size( ); j++ )
                        {

                            boost::shared_ptr< AccelerationPartial > currentAccelerationPartial =
                                    createAnalyticalAccelerationPartial( innerAccelerationIterator->second[ j ],
                                                                         std::make_pair( acceleratedBody, acceleratedBodyObject ),
                                                                         std::make_pair( acceleratingBody, acceleratingBodyObject ),
                                                                         bodyMap, parametersToEstimate, accelerationPartialsMap );

                            accelerationPartialVector.push_back( currentAccelerationPartial );
                            accelerationPartialsMap[ acceleratedBody ][ acceleratingBody ].push_back( currentAccelerationPartial );
                        }
                    }

                    // Add partials of current body's accelerations to vector.
                    accelerationPartialsList[ i ] = accelerationPartialVector;

                }
                else
                {
                    std::cerr<<"Error when making acceleration partials map, could not identify parameter"<<std::endl;
                }
            }
        }
    }
    return accelerationPartialsList;

}

}

}

}
#endif // TUDAT_CREATEACCELERATIONPARTIALS_H
