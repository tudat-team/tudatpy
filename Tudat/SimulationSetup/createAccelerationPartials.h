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

//! Function to create a single acceleration partial derivative object.
/*!
 *  Function to create a single acceleration partial derivative object.
 *  \param accelerationModel Acceleration model for which a partial derivative is to be computed.
 *  \param acceleratedBody Pair of name and object of body undergoing acceleration
 *  \param acceleratingBody Pair of name and object of body exerting acceleration
 *  \param bodyMap List of all body objects
 *  \return Single acceleration partial derivative object.
 */
template< typename InitialStateParameterType = double >
boost::shared_ptr< AccelerationPartial > createAnalyticalAccelerationPartial(
        boost::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > accelerationModel,
        const std::pair< std::string, boost::shared_ptr< simulation_setup::Body > > acceleratedBody,
        const std::pair< std::string, boost::shared_ptr< simulation_setup::Body > > acceleratingBody,
        const simulation_setup::NamedBodyMap bodyMap )
{

    using namespace basic_astrodynamics;
    using namespace electro_magnetism;
    using namespace aerodynamics;

    boost::shared_ptr< AccelerationPartial > accelerationPartial;

    // Identify current acceleration model type
    AvailableAcceleration accelerationType = getAccelerationModelType( accelerationModel );
    switch( accelerationType )
    {
    case central_gravity:

        // Check if identifier is consistent with type.
        if( boost::dynamic_pointer_cast< CentralGravitationalAccelerationModel3d >( accelerationModel ) == NULL )
        {
            throw std::runtime_error( "Acceleration class type does not match acceleration type (central_gravity) when making acceleration partial" );
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
            throw std::runtime_error( "Acceleration class type does not match acceleration type (third_body_central_gravity) when making acceleration partial" );
        }
        else
        {
            boost::shared_ptr< ThirdBodyCentralGravityAcceleration > thirdBodyAccelerationModel  =
                    boost::dynamic_pointer_cast< ThirdBodyCentralGravityAcceleration >( accelerationModel );

            // Create partials for constituent central gravity accelerations
            boost::shared_ptr< CentralGravitationPartial > accelerationPartialForBodyUndergoingAcceleration =
                    boost::dynamic_pointer_cast< CentralGravitationPartial >(
                        createAnalyticalAccelerationPartial(
                            thirdBodyAccelerationModel->getAccelerationModelForBodyUndergoingAcceleration( ),
                            acceleratedBody, acceleratingBody, bodyMap ) );
            boost::shared_ptr< CentralGravitationPartial > accelerationPartialForCentralBody =
                    boost::dynamic_pointer_cast< CentralGravitationPartial >(
                        createAnalyticalAccelerationPartial(
                            thirdBodyAccelerationModel->getAccelerationModelForCentralBody( ),
                            std::make_pair( thirdBodyAccelerationModel->getCentralBodyName( ),
                                            bodyMap.at( thirdBodyAccelerationModel->getCentralBodyName( ) ) ),
                            acceleratingBody, bodyMap  ) );

            // Create partial-calculating object.
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
            throw std::runtime_error( "Acceleration class type does not match acceleration type (cannon_ball_radiation_pressure) when making acceleration partial" );
        }
        else
        {
            std::map< std::string, boost::shared_ptr< RadiationPressureInterface > > radiationPressureInterfaces =
                    acceleratedBody.second->getRadiationPressureInterfaces( );

            if( radiationPressureInterfaces.count( acceleratingBody.first ) == 0 )
            {
                throw std::runtime_error( "No radiation pressure coefficient interface found when making acceleration partial." );
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
        std::string errorMessage = "Acceleration model " + boost::lexical_cast< std::string >( accelerationType ) +
                " not found when making acceleration partial";
        throw std::runtime_error( errorMessage );
        break;
    }

    return accelerationPartial;
}

//! This function creates acceleration partial objects for translational dynamics
/*!
 *  This function creates acceleration partial objects for translational dynamics from acceleration models and
 *  list of bodies' states of which derivatives are needed. The return type is an StateDerivativePartialsMap,
 *  a standardized type for communicating such lists of these objects.
 *  \param accelerationMap Map of maps containing list of acceleration models, identifying which acceleration acts on which
 *   body.
 *  \param bodyMap List of body objects constituting environment for calculations.
 *  \param parametersToEstimate List of parameters which are to be estimated.
 *  \return List of acceleration-partial-calculating objects in StateDerivativePartialsMap type.
 */
template< typename InitialStateParameterType >
StateDerivativePartialsMap createAccelerationPartialsMap(
        const basic_astrodynamics::AccelerationMap& accelerationMap,
        const simulation_setup::NamedBodyMap& bodyMap,
        const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > >
        parametersToEstimate )
{
    // Declare return map.
    StateDerivativePartialsMap accelerationPartialsList;
    std::map< std::string, std::map< std::string, std::vector< boost::shared_ptr< AccelerationPartial > > > >
            accelerationPartialsMap;

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
                if( ( initialDynamicalParameters.at( i )->getParameterName( ).first ==
                      estimatable_parameters::initial_body_state ) )
                {
                    // Get object for body undergoing acceleration
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
                        // Get object for body exerting acceleration
                        std::string acceleratingBody = innerAccelerationIterator->first;
                        boost::shared_ptr< simulation_setup::Body > acceleratingBodyObject;
                        if( acceleratingBody != "" )
                        {
                            acceleratingBodyObject = bodyMap.at( acceleratingBody );
                        }

                        for( unsigned int j = 0; j < innerAccelerationIterator->second.size( ); j++ )
                        {
                            // Create single partial object
                            boost::shared_ptr< AccelerationPartial > currentAccelerationPartial =
                                    createAnalyticalAccelerationPartial(
                                        innerAccelerationIterator->second[ j ],
                                        std::make_pair( acceleratedBody, acceleratedBodyObject ),
                                        std::make_pair( acceleratingBody, acceleratingBodyObject ),
                                        bodyMap );

                            accelerationPartialVector.push_back( currentAccelerationPartial );
                            accelerationPartialsMap[ acceleratedBody ][ acceleratingBody ].push_back(
                                        currentAccelerationPartial );
                        }
                    }

                    // Add partials of current body's accelerations to vector.
                    accelerationPartialsList[ i ] = accelerationPartialVector;

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
