
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModelTypes.h"
#include "Tudat/Astrodynamics/Propagators/createEnvironmentUpdater.h"

namespace tudat
{

namespace propagators
{


void checkValidityOfRequiredEnvironmentUpdates(
        std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >& requestedUpdates,
        const NamedBodyMap& bodyMap )
{
    using namespace propagators;

    std::map< propagators::EnvironmentModelsToUpdate, std::vector< int > > fullMapIndicesToRemove;

    for( std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >::iterator updateIterator =
         requestedUpdates.begin( ); updateIterator != requestedUpdates.end( ); updateIterator++ )
    {
        std::vector< int > indicesToRemove;
        for( unsigned int i = 0; i < updateIterator->second.size( ); i++ )
        {
            if( updateIterator->second.at( i ) != "" )
            {
                if( bodyMap.count( updateIterator->second.at( i ) ) == 0 )
                {
                    std::cerr<<"Error when making environment model update settings, could not find body "<<updateIterator->second.at( i )<<std::endl;
                }

                switch( updateIterator->first )
                {
                case body_transational_state_update:
                {
                    if( bodyMap.at( updateIterator->second.at( i ) )->getEphemeris( ) == NULL )
                    {
                        std::cerr<<"Error when making environment model update settings, could not find ephemeris of body "<<
                                   updateIterator->second.at( i )<<std::endl;
                        indicesToRemove.push_back( i );
                    }
                    break;
                }
                case body_rotational_state_update:
                {
                    if( bodyMap.at( updateIterator->second.at( i ) )->getRotationalEphemeris( ) == NULL )
                    {
                        std::cerr<<"Error when making environment model update settings, could not find rotational ephemeris of body "<<
                                   updateIterator->second.at( i )<<std::endl;
                        indicesToRemove.push_back( i );
                    }
                    break;
                }
                case spherical_harmonic_gravity_field_update:
                {
                    boost::shared_ptr< gravitation::SphericalHarmonicsGravityField > gravityFieldModel =
                            boost::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >(
                                boost::dynamic_pointer_cast< bodies::CelestialBody >( bodyMap.at( updateIterator->second.at( i ) ) )
                                ->getGravityFieldModel( ) );
                    if( gravityFieldModel == NULL )
                    {
                        std::cerr<<"Error when making environment model update settings, could not find spherical harmonic gravity field of body "<<
                                   updateIterator->second.at( i )<<std::endl;
                        indicesToRemove.push_back( i );
                    }
                    else if( boost::dynamic_pointer_cast< gravitation::TimeDependentSphericalHarmonicsGravityField >( gravityFieldModel ) == NULL )
                    {
                        indicesToRemove.push_back( i );
                    }
                    break;
                }
                case vehicle_flight_conditions_update:
                {
                    boost::shared_ptr< aerodynamics::FlightConditions > flightConditions =
                            boost::dynamic_pointer_cast< bodies::Vehicle >( bodyMap.at( updateIterator->second.at( i ) ) )->getFlightConditions( );
                    if( flightConditions == NULL )
                    {
                        std::cerr<<"Error when making environment model update settings, could not find flight conditions of body "<<
                                   updateIterator->second.at( i )<<std::endl;
                        indicesToRemove.push_back( i );
                    }
                    break;
                }
                case radiation_pressure_interface_update:
                    break;
                case vehicle_part_orientation_update:
                {
                    boost::shared_ptr< bodies::CurrentVehiclePartOrientations > vehiclePartOrientationModel =
                            boost::dynamic_pointer_cast< bodies::Vehicle >( bodyMap.at( updateIterator->second.at( i ) ) )
                            ->getVehiclePartOrientationModel( );
                    if( vehiclePartOrientationModel == NULL )
                    {
                        indicesToRemove.push_back( i );
                    }
                    break;
                }
                case body_mass_update:
                {
                    if( boost::dynamic_pointer_cast< bodies::Vehicle >( bodyMap.at( updateIterator->second.at( i ) ) ) == NULL )
                    {
                        std::cerr<<"Error when making environment model update settings "<<updateIterator->second.at( i )<<
                                   " is not a vehicle and does not have a mass model "<<std::endl;
                        indicesToRemove.push_back( i );
                    }
                    break;
                }
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

std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > createRotationalEquationsOfMotionEnvironmentUpdaterSettings(
        const basic_astrodynamics::TorqueModelMap& torqueModels, const NamedBodyMap& bodyMap )
{
    using namespace basic_astrodynamics;

    std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > environmentModelsToUpdate;
    std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > singleTorqueUpdateNeeds;

    for( TorqueModelMap::const_iterator acceleratedBodyIterator = torqueModels.begin( );
         acceleratedBodyIterator != torqueModels.end( ); acceleratedBodyIterator++ )
    {
        for( SingleBodyTorqueModelMap::const_iterator torqueModelIterator = acceleratedBodyIterator->second.begin( );
             torqueModelIterator != acceleratedBodyIterator->second.end( ); torqueModelIterator++ )
        {
            singleTorqueUpdateNeeds.clear( );
            for( unsigned int i = 0; i < torqueModelIterator->second.size( ); i++ )
            {
                AvailableTorque currentAccelerationModelType =
                        getTorqueModelType( torqueModelIterator->second.at( i ) );

                switch( currentAccelerationModelType )
                {
                case second_order_gravitational_torque:
                    singleTorqueUpdateNeeds[ body_transational_state_update ].push_back( acceleratedBodyIterator->first );
                    singleTorqueUpdateNeeds[ body_transational_state_update ].push_back( torqueModelIterator->first );

                    break;
                case dissipative_torque_model:
                    break;
                default:
                    std::cerr<<"Error, update information not found for torque model "<<currentAccelerationModelType<<std::endl;
                    break;
                }
            }

            checkValidityOfRequiredEnvironmentUpdates( singleTorqueUpdateNeeds, bodyMap );
            addEnvironmentUpdates( environmentModelsToUpdate, singleTorqueUpdateNeeds );
        }
    }

    return environmentModelsToUpdate;
}

std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > createTranslationalEquationsOfMotionEnvironmentUpdaterSettings(
        const AccelerationMap& translationalAccelerationModels, const NamedBodyMap& bodyMap )
{
    using namespace astrodynamics::acceleration_models;
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
                case constant_drag_aerodynamic:
                    singleAccelerationUpdateNeeds[ body_rotational_state_update ].push_back( accelerationModelIterator->first );
                    singleAccelerationUpdateNeeds[ vehicle_flight_conditions_update ].push_back( acceleratedBodyIterator->first );
                    singleAccelerationUpdateNeeds[ body_mass_update ].push_back( acceleratedBodyIterator->first );
                    break;
                case aerodynamic:
                    singleAccelerationUpdateNeeds[ body_rotational_state_update ].push_back( accelerationModelIterator->first );
                    singleAccelerationUpdateNeeds[ vehicle_flight_conditions_update ].push_back( acceleratedBodyIterator->first );
                    singleAccelerationUpdateNeeds[ vehicle_part_orientation_update ].push_back( acceleratedBodyIterator->first  );
                    singleAccelerationUpdateNeeds[ body_mass_update ].push_back( acceleratedBodyIterator->first );
                    break;
                case cannon_ball_radiation_pressure:
                    singleAccelerationUpdateNeeds[ radiation_pressure_interface_update ].push_back( acceleratedBodyIterator->first );
                    singleAccelerationUpdateNeeds[ body_mass_update ].push_back( acceleratedBodyIterator->first );
                    break;
                case spherical_harmonic_gravity:
                    singleAccelerationUpdateNeeds[ spherical_harmonic_gravity_field_update ].push_back( accelerationModelIterator->first );
                    singleAccelerationUpdateNeeds[ body_rotational_state_update ].push_back( accelerationModelIterator->first );
                    break;
                case third_body_spherical_harmonic_gravity:
                    singleAccelerationUpdateNeeds[ spherical_harmonic_gravity_field_update ].push_back( accelerationModelIterator->first );
                    singleAccelerationUpdateNeeds[ body_rotational_state_update ].push_back( accelerationModelIterator->first );
                    break;
                case mutual_spherical_harmonic_gravity:
                    singleAccelerationUpdateNeeds[ spherical_harmonic_gravity_field_update ].push_back( accelerationModelIterator->first );
                    singleAccelerationUpdateNeeds[ body_rotational_state_update ].push_back( accelerationModelIterator->first );
                    singleAccelerationUpdateNeeds[ spherical_harmonic_gravity_field_update ].push_back(  acceleratedBodyIterator->first );
                    singleAccelerationUpdateNeeds[ body_rotational_state_update ].push_back(  acceleratedBodyIterator->first );
                    break;
                case third_body_mutual_spherical_harmonic_gravity:
                    singleAccelerationUpdateNeeds[ spherical_harmonic_gravity_field_update ].push_back( accelerationModelIterator->first );
                    singleAccelerationUpdateNeeds[ body_rotational_state_update ].push_back( accelerationModelIterator->first );
                    singleAccelerationUpdateNeeds[ spherical_harmonic_gravity_field_update ].push_back(  acceleratedBodyIterator->first );
                    singleAccelerationUpdateNeeds[ body_rotational_state_update ].push_back(  acceleratedBodyIterator->first );
                    break;
                case panelled_radiation_pressure:
                    singleAccelerationUpdateNeeds[ radiation_pressure_interface_update ].push_back( acceleratedBodyIterator->first );
                    singleAccelerationUpdateNeeds[ vehicle_part_orientation_update ].push_back( acceleratedBodyIterator->first  );
                    singleAccelerationUpdateNeeds[ body_mass_update ].push_back( acceleratedBodyIterator->first );
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

std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > createProperTimeEquationEnvironmentUpdaterSettings(
        const boost::shared_ptr< RelativisticTimeStatePropagatorSettings > stateDerivativeModel, const NamedBodyMap& bodyMap )
{
    using namespace relativity;

    std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > environmentModelsToUpdate;
    environmentModelsToUpdate[ body_transational_state_update ].push_back( stateDerivativeModel->getReferencePointId( ).first );

    switch( stateDerivativeModel->getRelativisticStateDerivativeType( ) )
    {
    case first_order_barycentric_to_bodycentric:
    {
        boost::shared_ptr< FirstOrderBodycentricRelativisticTimePropagatorSettings > firstOrderSettings =
                boost::dynamic_pointer_cast< FirstOrderBodycentricRelativisticTimePropagatorSettings >( stateDerivativeModel );
        for( unsigned int i = 0; i < firstOrderSettings->getExternalBodyList( ).size( ); i++ )
        {
            environmentModelsToUpdate[ body_transational_state_update ].push_back( firstOrderSettings->getExternalBodyList( ).at( i ) );
        }

        std::map< std::string, std::pair< int, int > > shOrders = firstOrderSettings->getSphericalHarmonicGravityExpansions( );
        for( std::map< std::string, std::pair< int, int > >::const_iterator shIterator = shOrders.begin( );
             shIterator != shOrders.end( ); shIterator++ )
        {
            environmentModelsToUpdate[ spherical_harmonic_gravity_field_update ].push_back( shIterator->first );
        }

        break;
    }
    case second_order_barycentric_to_bodycentric:
    {
        boost::shared_ptr< SecondOrderBodyCenteredRelativisticTimeConverterSettings > secondOrderSettings =
                boost::dynamic_pointer_cast< SecondOrderBodyCenteredRelativisticTimeConverterSettings >( stateDerivativeModel );
        for( unsigned int i = 0; i < secondOrderSettings->getExternalBodyList( ).size( ); i++ )
        {
            environmentModelsToUpdate[ body_transational_state_update ].push_back( secondOrderSettings->getExternalBodyList( ).at( i ) );
        }

        std::map< std::string, std::pair< int, int > > shOrders = secondOrderSettings->getSphericalHarmonicGravityExpansions( );
        for( std::map< std::string, std::pair< int, int > >::const_iterator shIterator = shOrders.begin( ); shIterator != shOrders.end( );
             shIterator++ )
        {
            environmentModelsToUpdate[ spherical_harmonic_gravity_field_update ].push_back( shIterator->first );
        }

        break;
    }
    case first_order_bodycentric_to_topocentric:
    {
        boost::shared_ptr< BodycenteredToTopocentricTimePropagatorSettings > firstOrderSettings =
                boost::dynamic_pointer_cast< BodycenteredToTopocentricTimePropagatorSettings >( stateDerivativeModel );
        for( unsigned int i = 0; i < firstOrderSettings->getTopocentricExternalBodies( ).size( ); i++ )
        {
            environmentModelsToUpdate[ body_transational_state_update ].push_back( firstOrderSettings->getTopocentricExternalBodies( ).at( i ) );
        }

        if( firstOrderSettings->getMaximumSphericalHarmonicDegree( ) > 0 )
        {
            environmentModelsToUpdate[ spherical_harmonic_gravity_field_update ].push_back( firstOrderSettings->getReferencePointId( ).first );
        }

        break;
    }
    case direct_from_metric:
    {
        if( evaluatedMetricObjects.count( stateDerivativeModel->getReferencePointId( ) ) == 0 )
        {
            evaluatedMetricObjects[ stateDerivativeModel->getReferencePointId( ) ] = baseMetric->Clone( );
        }
        boost::shared_ptr< Metric > metricToUse = evaluatedMetricObjects[ stateDerivativeModel->getReferencePointId( ) ];

        std::vector< std::string > bodyList;
        std::vector< std::string > shList;

        if( boost::dynamic_pointer_cast< SolarSystemMetric >( metricToUse ) != NULL )
        {
            boost::shared_ptr< SolarSystemMetric > solarSystemMetric =
                    boost::dynamic_pointer_cast< SolarSystemMetric >( metricToUse );
            bodyList = solarSystemMetric->getBodyList( );
            std::map< int, boost::shared_ptr< SphericalHarmonicWrapper > > shFunctionList =
                    solarSystemMetric->getBodySphericalHarmonicGravityWrappers( );
            for( std::map< int, boost::shared_ptr< SphericalHarmonicWrapper > >::const_iterator shIterator = shFunctionList.begin( );
                 shIterator != shFunctionList.end( ); shIterator++ )
            {
                if( shIterator->first < static_cast< int >( bodyList.size( ) ) )
                {
                    shList.push_back( bodyList.at( shIterator->first ) );
                }
                else
                {
                    std::cerr<<"Error when making proper time updater direct from metric, sh index is too high"<<std::endl;
                }
            }
        }
        else if( boost::dynamic_pointer_cast< HarmonicSchwarzschildMetric >( metricToUse ) != NULL )
        {

        }
        else
        {
            std::cerr<<"Error when making proper time updater direct from metric, did not recognize metric type"<<std::endl;
        }

        for( unsigned int i = 0; i < bodyList.size( ); i++ )
        {
            environmentModelsToUpdate[ body_transational_state_update ].push_back( bodyList.at( i ) );
        }

        for( unsigned int i = 0; i < shList.size( ); i++ )
        {
            environmentModelsToUpdate[ spherical_harmonic_gravity_field_update ].push_back( shList.at( i ) );
        }

        break;
    }
    default:
        std::cerr<<"Error when making proper time updater, did not recognize equation type: "<<
                   stateDerivativeModel->getRelativisticStateDerivativeType( )<<std::endl;
    }
    checkValidityOfRequiredEnvironmentUpdates( environmentModelsToUpdate, bodyMap );

    return environmentModelsToUpdate;
}

std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > createFullEnvironmentUpdaterSettings(
        const NamedBodyMap& bodyMap )
{
    using namespace astrodynamics::acceleration_models;
    using namespace electro_magnetism;
    using namespace gravitation;
    using namespace propagators;

    std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > environmentModelsToUpdate;
    std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > singleAccelerationUpdateNeeds;

    // Iterate over all bodies.
    for( NamedBodyMap::const_iterator bodyIterator = bodyMap.begin( ); bodyIterator != bodyMap.end( ); bodyIterator++ )
    {
        singleAccelerationUpdateNeeds.clear( );

        // Check if current body is a vehicle.
        boost::shared_ptr< bodies::Vehicle > vehicleToEvaluate =
                boost::dynamic_pointer_cast< bodies::Vehicle >( bodyIterator->second );
        if( vehicleToEvaluate != NULL )
        {
            // Check if current body has flight conditions set.
            if( vehicleToEvaluate->getFlightConditions( ) != NULL )
            {
                // If vehicle has flight conditions, add flight conditions update function to update list.
                singleAccelerationUpdateNeeds[ vehicle_flight_conditions_update ].push_back( bodyIterator->first );
            }

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

        // If body has time dependent gravity model, update
        boost::shared_ptr< bodies::CelestialBody > currentCelestialBody = boost::dynamic_pointer_cast< bodies::CelestialBody >
                ( bodyIterator->second );
        if( currentCelestialBody != NULL )
        {
            boost::shared_ptr< TimeDependentSphericalHarmonicsGravityField > gravityField =
                    boost::dynamic_pointer_cast< TimeDependentSphericalHarmonicsGravityField >
                    ( currentCelestialBody->getGravityFieldModel( ) );
            if( gravityField != NULL )
            {
                singleAccelerationUpdateNeeds[ spherical_harmonic_gravity_field_update ].push_back( bodyIterator->first );

            }
        }

        if( vehicleToEvaluate != NULL )
        {
            singleAccelerationUpdateNeeds[ body_mass_update ].push_back( bodyIterator->first );

            // Check if body has part orientation model that requires updating
            if( vehicleToEvaluate->getVehiclePartOrientationModel( ) != NULL )
            {
                singleAccelerationUpdateNeeds[ vehicle_part_orientation_update ].push_back( bodyIterator->first );

            }
        }

        checkValidityOfRequiredEnvironmentUpdates( singleAccelerationUpdateNeeds, bodyMap );
        addEnvironmentUpdates( environmentModelsToUpdate, singleAccelerationUpdateNeeds );
    }
    return environmentModelsToUpdate;
}

}

}


