
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModelTypes.h"
#include "Tudat/Astrodynamics/Propagators/createEnvironmentUpdater.h"

namespace tudat
{

namespace propagators
{


void checkValidityOfRequiredEnvironmentUpdates(
        std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >&
        requestedUpdates,
        const simulation_setup::NamedBodyMap& bodyMap )
{
    using namespace propagators;

    std::map< propagators::EnvironmentModelsToUpdate, std::vector< int > > fullMapIndicesToRemove;

    for( std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >::iterator
         updateIterator = requestedUpdates.begin( );
         updateIterator != requestedUpdates.end( ); updateIterator++ )
    {
        std::vector< int > indicesToRemove;
        for( unsigned int i = 0; i < updateIterator->second.size( ); i++ )
        {
            if( updateIterator->second.at( i ) != "" )
            {
                if( bodyMap.count( updateIterator->second.at( i ) ) == 0 )
                {
                    throw std::runtime_error( "Error when making environment model update settings, could not find body " +
                                              boost::lexical_cast< std::string>( updateIterator->second.at( i ) ) );
                }

                switch( updateIterator->first )
                {
                case body_transational_state_update:
                {
                    if( bodyMap.at( updateIterator->second.at( i ) )->getEphemeris( ) == NULL )
                    {
                        throw std::runtime_error( "Error when making environment model update settings, could not find ephemeris of body " +
                                   boost::lexical_cast< std::string>( updateIterator->second.at( i ) ) );
                        indicesToRemove.push_back( i );
                    }
                    break;
                }
                case body_rotational_state_update:
                {
                    if( bodyMap.at( updateIterator->second.at( i ) )->getRotationalEphemeris( ) == NULL )
                    {
                        throw std::runtime_error( "Error when making environment model update settings, could not find rotational ephemeris of body " +
                                   boost::lexical_cast< std::string>( updateIterator->second.at( i ) ) );
                        indicesToRemove.push_back( i );
                    }
                    break;
                }
                case vehicle_flight_conditions_update:
                {
                    boost::shared_ptr< aerodynamics::FlightConditions > flightConditions = bodyMap.at( updateIterator->second.at( i ) )->getFlightConditions( );
                    if( flightConditions == NULL )
                    {
                        throw std::runtime_error( "Error when making environment model update settings, could not find flight conditions of body " +
                                   boost::lexical_cast< std::string>( updateIterator->second.at( i ) ) );
                        indicesToRemove.push_back( i );
                    }
                    break;
                }
                case radiation_pressure_interface_update:
                {
                    std::map< std::string, boost::shared_ptr< electro_magnetism::RadiationPressureInterface > >
                     radiationPressureInterfaces = bodyMap.at( updateIterator->second.at( i ) )->getRadiationPressureInterfaces( );
                    if( radiationPressureInterfaces.size( ) == 0 )
                    {
                        throw std::runtime_error( "Error when making environment model update settings, could not find radiation pressure interface of body " +
                                   boost::lexical_cast< std::string>( updateIterator->second.at( i ) ) );
                        indicesToRemove.push_back( i );
                    }
                    break;
                }
                case body_mass_update:
                    break;
                }
            }
            fullMapIndicesToRemove[ updateIterator->first ] = indicesToRemove;
        }
    }

    for( std::map< propagators::EnvironmentModelsToUpdate, std::vector< int > >::iterator removalIterator =
         fullMapIndicesToRemove.begin( ); removalIterator != fullMapIndicesToRemove.end( ); removalIterator++ )
    {
        for( int i = removalIterator->second.size( ) - 1; i >= 0; i-- )
        {
            requestedUpdates[ removalIterator->first ].erase( requestedUpdates[ removalIterator->first ].begin( ) + removalIterator->second.at( i ) );
        }
    }
}

void addEnvironmentUpdates(
        std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >& environmentUpdateList,
        const std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > updatesToAdd )
{
    for( std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >::const_iterator environmentUpdateIterator =
         updatesToAdd.begin( ); environmentUpdateIterator != updatesToAdd.end( ); environmentUpdateIterator++ )
    {
        bool addCurrentUpdate = 0;
        for( unsigned int i = 0; i < environmentUpdateIterator->second.size( ); i++ )
        {
            addCurrentUpdate = 0;
            if( environmentUpdateList.count( environmentUpdateIterator->first ) == 0 )
            {
                addCurrentUpdate = 1;
            }
            else if( std::find( environmentUpdateList.at( environmentUpdateIterator->first ).begin( ),
                                environmentUpdateList.at( environmentUpdateIterator->first ).end( ),
                                environmentUpdateIterator->second.at( i ) ) ==
                     environmentUpdateList.at( environmentUpdateIterator->first ).end( ) )
            {
                addCurrentUpdate = 1;
            }

            if( addCurrentUpdate )
            {
                environmentUpdateList[ environmentUpdateIterator->first ].push_back( environmentUpdateIterator->second.at( i ) );
            }
        }
    }
}


std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > createTranslationalEquationsOfMotionEnvironmentUpdaterSettings(
        const basic_astrodynamics::AccelerationMap& translationalAccelerationModels, const simulation_setup::NamedBodyMap& bodyMap )
{
    using namespace basic_astrodynamics;
    using namespace propagators;

    std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > environmentModelsToUpdate;
    std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > singleAccelerationUpdateNeeds;

    for( AccelerationMap::const_iterator acceleratedBodyIterator = translationalAccelerationModels.begin( );
         acceleratedBodyIterator != translationalAccelerationModels.end( ); acceleratedBodyIterator++ )
    {
        for( SingleBodyAccelerationMap::const_iterator accelerationModelIterator = acceleratedBodyIterator->second.begin( );
             accelerationModelIterator != acceleratedBodyIterator->second.end( ); accelerationModelIterator++ )
        {
            singleAccelerationUpdateNeeds.clear( );
            for( unsigned int i = 0; i < accelerationModelIterator->second.size( ); i++ )
            {
                AvailableAcceleration currentAccelerationModelType =
                        getAccelerationModelType( accelerationModelIterator->second.at( i ) );

                singleAccelerationUpdateNeeds[ body_transational_state_update ].push_back( accelerationModelIterator->first );
                singleAccelerationUpdateNeeds[ body_transational_state_update ].push_back( acceleratedBodyIterator->first );

                switch( currentAccelerationModelType )
                {
                case aerodynamic:
                    singleAccelerationUpdateNeeds[ body_rotational_state_update ].push_back( accelerationModelIterator->first );
                    singleAccelerationUpdateNeeds[ vehicle_flight_conditions_update ].push_back( acceleratedBodyIterator->first );
                    singleAccelerationUpdateNeeds[ body_mass_update ].push_back( acceleratedBodyIterator->first );
                    break;
                case cannon_ball_radiation_pressure:
                    singleAccelerationUpdateNeeds[ radiation_pressure_interface_update ].push_back( acceleratedBodyIterator->first );
                    singleAccelerationUpdateNeeds[ body_mass_update ].push_back( acceleratedBodyIterator->first );
                    break;
                case spherical_harmonic_gravity:
                    singleAccelerationUpdateNeeds[ body_rotational_state_update ].push_back( accelerationModelIterator->first );
                    break;
                case third_body_spherical_harmonic_gravity:
                    singleAccelerationUpdateNeeds[ body_rotational_state_update ].push_back( accelerationModelIterator->first );
                    break;
                default:
                    break;
                }

            }

            checkValidityOfRequiredEnvironmentUpdates( singleAccelerationUpdateNeeds, bodyMap );
            addEnvironmentUpdates( environmentModelsToUpdate, singleAccelerationUpdateNeeds );
        }
    }

    return environmentModelsToUpdate;
}


std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > createFullEnvironmentUpdaterSettings(
        const simulation_setup::NamedBodyMap& bodyMap )
{
    using namespace basic_astrodynamics;
    using namespace electro_magnetism;
    using namespace gravitation;
    using namespace propagators;

    std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > environmentModelsToUpdate;
    std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > singleAccelerationUpdateNeeds;

    // Iterate over all bodies.
    for( simulation_setup::NamedBodyMap::const_iterator bodyIterator = bodyMap.begin( ); bodyIterator != bodyMap.end( ); bodyIterator++ )
    {
        singleAccelerationUpdateNeeds.clear( );

        // Check if current body is a vehicle.

        // Check if current body has flight conditions set.
        if( bodyIterator->second ->getFlightConditions( ) != NULL )
        {
            // If vehicle has flight conditions, add flight conditions update function to update list.
            singleAccelerationUpdateNeeds[ vehicle_flight_conditions_update ].push_back( bodyIterator->first );
        }

        // Get body radiation pressure interface(s) (one per source)
        std::map< std::string, boost::shared_ptr< RadiationPressureInterface > > radiationPressureInterfaces =
                bodyIterator->second->getRadiationPressureInterfaces( );

        // Add each interface update function to update list.
        for( std::map< std::string, boost::shared_ptr< RadiationPressureInterface > > ::iterator iterator = radiationPressureInterfaces.begin( );
             iterator != radiationPressureInterfaces.end( ); iterator++ )
        {
            singleAccelerationUpdateNeeds[ radiation_pressure_interface_update ].push_back( bodyIterator->first );
        }

        // If body has rotation model, update rotational state in each time step.
        boost::shared_ptr< ephemerides::RotationalEphemeris > rotationalEphemeris = bodyIterator->second->getRotationalEphemeris( );
        if( rotationalEphemeris != NULL )
        {
            singleAccelerationUpdateNeeds[ body_rotational_state_update ].push_back( bodyIterator->first );
        }

        singleAccelerationUpdateNeeds[ body_mass_update ].push_back( bodyIterator->first );

        checkValidityOfRequiredEnvironmentUpdates( singleAccelerationUpdateNeeds, bodyMap );
        addEnvironmentUpdates( environmentModelsToUpdate, singleAccelerationUpdateNeeds );
    }
    return environmentModelsToUpdate;
}

}

}


