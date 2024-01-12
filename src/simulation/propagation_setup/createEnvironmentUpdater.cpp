/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/basic_astro/accelerationModelTypes.h"
#include "tudat/astro/basic_astro/torqueModelTypes.h"
#include "tudat/simulation/propagation_setup/createEnvironmentUpdater.h"
#include "tudat/simulation/environment_setup/createFlightConditions.h"

namespace tudat
{

namespace propagators
{

//! Function to check whether the requested environment updates are
//! possible with the existing environment.
void checkValidityOfRequiredEnvironmentUpdates(
        const std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >&
        requestedUpdates, const simulation_setup::SystemOfBodies& bodies )
{
    using namespace propagators;

    // Iterate over all environment update types.
    for( std::map< propagators::EnvironmentModelsToUpdate,
         std::vector< std::string > >::const_iterator
         updateIterator = requestedUpdates.begin( );
         updateIterator != requestedUpdates.end( ); updateIterator++ )
    {
        // Iterate over all updated bodies.
        for( unsigned int i = 0; i < updateIterator->second.size( ); i++ )
        {
            // Ignore global required updates.
            if( updateIterator->second.at( i ) != "" )
            {
                // Check if body exists.
                if( bodies.count( updateIterator->second.at( i ) ) == 0 )
                {
                    throw std::runtime_error(
                                "Error when making environment model update settings, could not find body "
                                + updateIterator->second.at( i ) );
                }

                // Check if requested environment model exists.
                switch( updateIterator->first )
                {
                case body_translational_state_update:
                {
                    if( bodies.at( updateIterator->second.at( i ) )->getEphemeris( ) == nullptr )
                    {
                        throw std::runtime_error(
                                    "Error when making environment model update settings, could not find ephemeris of body "
                                    + updateIterator->second.at( i ) );
                    }
                    break;
                }
                case body_rotational_state_update:
                {
                    if( ( bodies.at( updateIterator->second.at( i ) )->
                          getRotationalEphemeris( ) == nullptr )
//                            &&
//                            ( bodies.at( updateIterator->second.at( i ) )->getDependentOrientationCalculator( ) == nullptr )
                            )
                    {
                        throw std::runtime_error(
                                    "Error when making environment model update settings, could not find rotational ephemeris of body "
                                    + updateIterator->second.at( i ) );
                    }
                    break;
                }
                case spherical_harmonic_gravity_field_update:
                {
                    std::shared_ptr< gravitation::SphericalHarmonicsGravityField >
                            gravityFieldModel =
                            std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >(
                                bodies.at( updateIterator->second.at( i ) )->getGravityFieldModel( ) );
                    if( gravityFieldModel == nullptr )
                    {
                        throw std::runtime_error(
                                    "Error when making environment model update settings, could not find spherical harmonic gravity field of body "
                                    + updateIterator->second.at( i ) );
                    }
                    break;
                }
                case vehicle_flight_conditions_update:
                {
                    std::shared_ptr< aerodynamics::FlightConditions > flightConditions = bodies.at(
                                updateIterator->second.at( i ) )->getFlightConditions( );
                    if( flightConditions == nullptr )
                    {
                        throw std::runtime_error(
                                    "Error when making environment model update settings, could not find flight conditions of body "
                                    + updateIterator->second.at( i ) );
                    }
                    break;
                }
                case body_mass_update:
                {
                    if( bodies.at( updateIterator->second.at( i ) )->getMassProperties( ) == nullptr )
                    {
                        throw std::runtime_error(
                                    "Error when making environment model update settings, no mass properties of body "
                                    + updateIterator->second.at( i ) );
                    }

                    break;
                }
                case body_mass_distribution_update:
                {
                    if( bodies.at( updateIterator->second.at( i ) )->getMassProperties( ) == nullptr )
                    {
                        throw std::runtime_error(
                            "Error when making environment model update settings, no mass properties of body "
                            + updateIterator->second.at( i ) );
                    }

                    break;
                }
                case radiation_source_model_update:
                {
                    std::shared_ptr< electromagnetism::RadiationSourceModel > radiationSourceModel = bodies.at(
                            updateIterator->second.at( i ) )->getRadiationSourceModel();
                    if( radiationSourceModel == nullptr )
                    {
                        throw std::runtime_error(
                                "Error when making environment model update settings, could not find radiation source model of body "
                                + updateIterator->second.at( i ) );
                    }
                    break;
                }
                case radiation_pressure_target_model_update:
                {
                    std::shared_ptr< electromagnetism::RadiationPressureTargetModel > radiationPressureTargetModel = bodies.at(
                            updateIterator->second.at( i ) )->getRadiationPressureTargetModel();
                    if( radiationPressureTargetModel == nullptr )
                    {
                        throw std::runtime_error(
                                "Error when making environment model update settings, could not find radiation pressure target model of body "
                                + updateIterator->second.at( i ) );
                    }
                    break;
                }
                case body_segment_orientation_update:
                {
                    std::shared_ptr< system_models::VehicleSystems > systemModels = bodies.at(
                        updateIterator->second.at( i ) )->getVehicleSystems( );
                    if( systemModels == nullptr )
                    {
                        throw std::runtime_error(
                            "Error when making environment model update settings, could not find vehicle systems of body "
                            + updateIterator->second.at( i ) );
                    }
                    break;
                }
                }
            }
        }
    }
}

//! Function that removes propagated states from the updated environment variables
void removePropagatedStatesFomEnvironmentUpdates(
        std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >& environmentModelsToUpdate,
        const std::map< IntegratedStateType, std::vector< std::tuple< std::string, std::string, PropagatorType > > >& integratedStateList )
{
    // Iterate over all propagated states
    for( auto it = integratedStateList.begin( ); it != integratedStateList.end( ); it++ )
    {
        // Iterate over all propagated bodies for current state
        for( unsigned int i = 0; i < it->second.size( ); i++ )
        {
            switch( it->first )
            {
            // Check for propagated translational states in update list, and remove if necessary
            case translational_state:
                if( environmentModelsToUpdate.count( body_translational_state_update ) > 0 )
                {
                    std::vector< std::string > bodiesToUpdate = environmentModelsToUpdate.at( body_translational_state_update );
                    std::vector< std::string >::iterator findIterator =
                            std::find( bodiesToUpdate.begin( ), bodiesToUpdate.end( ), std::get< 0 >( it->second.at( i ) ) );

                    if( findIterator != bodiesToUpdate.end( ) )
                    {
                        bodiesToUpdate.erase( findIterator );
                        environmentModelsToUpdate[ body_translational_state_update ] = bodiesToUpdate;

                    }
                }
                break;
                // Check for propagated rotational states in update list, and remove if necessary
            case rotational_state:
                if( environmentModelsToUpdate.count( body_rotational_state_update ) > 0 )
                {
                    std::vector< std::string > bodiesToUpdate = environmentModelsToUpdate.at( body_rotational_state_update );
                    std::vector< std::string >::iterator findIterator =
                            std::find( bodiesToUpdate.begin( ), bodiesToUpdate.end( ), std::get< 0 >( it->second.at( i ) ) );

                    if( findIterator != bodiesToUpdate.end( ) )
                    {
                        bodiesToUpdate.erase( findIterator );
                        environmentModelsToUpdate[ body_rotational_state_update ] = bodiesToUpdate;

                    }
                }
                break;
                // Check for propagated mass states in update list, and remove if necessary
            case body_mass_state:
//                if( environmentModelsToUpdate.count( body_mass_update ) > 0 )
//                {
//                    std::vector< std::string > bodiesToUpdate = environmentModelsToUpdate.at( body_mass_update );
//                    std::vector< std::string >::iterator findIterator =
//                            std::find( bodiesToUpdate.begin( ), bodiesToUpdate.end( ), std::get< 0 >( it->second.at( i ) ) );
//
//                    if( findIterator != bodiesToUpdate.end( ) )
//                    {
//                        bodiesToUpdate.erase( findIterator );
//                        environmentModelsToUpdate[ body_mass_update ] = bodiesToUpdate;
//
//                    }
//                }
                break;
            case custom_state:
                break;
            default:
                throw std::runtime_error( "Error when removing propagated states from environment updates, state type " +
                                          std::to_string( it->first ) + " not recognized." );
            }
        }
    }
}

//! Get list of required environment model update settings from torque models.
std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >
createRotationalEquationsOfMotionEnvironmentUpdaterSettings(
        const basic_astrodynamics::TorqueModelMap& torqueModels, const simulation_setup::SystemOfBodies& bodies,
        const std::vector< std::string > bodiesToIntegrate )
{
    using namespace basic_astrodynamics;

    std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > environmentModelsToUpdate;
    std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > singleTorqueUpdateNeeds;

    for( unsigned  int i = 0; i < bodiesToIntegrate.size( ); i++ )
    {
        singleTorqueUpdateNeeds[ body_mass_distribution_update ].push_back(
                bodiesToIntegrate.at( i ) );
    }
    addEnvironmentUpdates( environmentModelsToUpdate, singleTorqueUpdateNeeds );


    // Iterate over all bodies on which torques are being exerting
    for( TorqueModelMap::const_iterator acceleratedBodyIterator = torqueModels.begin( );
         acceleratedBodyIterator != torqueModels.end( ); acceleratedBodyIterator++ )
    {
        // Iterate over all bodies exerting on current body
        for( SingleBodyTorqueModelMap::const_iterator torqueModelIterator = acceleratedBodyIterator->second.begin( );
             torqueModelIterator != acceleratedBodyIterator->second.end( ); torqueModelIterator++ )
        {
            singleTorqueUpdateNeeds.clear( );
            for( unsigned int i = 0; i < torqueModelIterator->second.size( ); i++ )
            {
                AvailableTorque currentTorqueModelType =
                        getTorqueModelType( torqueModelIterator->second.at( i ) );

                switch( currentTorqueModelType )
                {
                case second_order_gravitational_torque:
                    singleTorqueUpdateNeeds[ body_translational_state_update ].push_back(
                                torqueModelIterator->first );
                    singleTorqueUpdateNeeds[ body_translational_state_update ].push_back(
                                acceleratedBodyIterator->first );
                    break;
                case spherical_harmonic_gravitational_torque:
                    singleTorqueUpdateNeeds[ body_translational_state_update ].push_back(
                                torqueModelIterator->first );
                    singleTorqueUpdateNeeds[ body_translational_state_update ].push_back(
                                acceleratedBodyIterator->first );
                    singleTorqueUpdateNeeds[ spherical_harmonic_gravity_field_update ].push_back(
                                acceleratedBodyIterator->first );
                    break;

                case aerodynamic_torque:
                    singleTorqueUpdateNeeds[ body_translational_state_update ].push_back(
                                torqueModelIterator->first );
                    singleTorqueUpdateNeeds[ body_translational_state_update ].push_back(
                                acceleratedBodyIterator->first );
                    singleTorqueUpdateNeeds[ body_rotational_state_update ].push_back(
                                torqueModelIterator->first );
                    singleTorqueUpdateNeeds[ vehicle_flight_conditions_update ].push_back(
                                acceleratedBodyIterator->first );
                    singleTorqueUpdateNeeds[ body_mass_update ].push_back(
                        acceleratedBodyIterator->first );
                    singleTorqueUpdateNeeds[ body_mass_distribution_update ].push_back(
                        acceleratedBodyIterator->first );
                    break;
                case inertial_torque:
                    break;
                case dissipative_torque:
                    break;
                case custom_torque:
                    break;
                default:
                    std::cerr << "Error, update information not found for torque model " << currentTorqueModelType << std::endl;
                    break;
                }
            }

            addEnvironmentUpdates( environmentModelsToUpdate, singleTorqueUpdateNeeds );
        }
    }

    return environmentModelsToUpdate;
}


//! Get list of required environment model update settings from translational acceleration models.
std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >
createTranslationalEquationsOfMotionEnvironmentUpdaterSettings(
        const basic_astrodynamics::AccelerationMap& translationalAccelerationModels,
        const simulation_setup::SystemOfBodies& bodies )
{
    using namespace basic_astrodynamics;
    using namespace propagators;

    std::map< propagators::EnvironmentModelsToUpdate,
            std::vector< std::string > > environmentModelsToUpdate;
    std::map< propagators::EnvironmentModelsToUpdate,
            std::vector< std::string > > singleAccelerationUpdateNeeds;

    // Iterate over all bodies being accelerated
    for( AccelerationMap::const_iterator acceleratedBodyIterator
         = translationalAccelerationModels.begin( );
         acceleratedBodyIterator != translationalAccelerationModels.end( );
         acceleratedBodyIterator++ )
    {
        // Iterate over all bodies exerting acceleration.
        for( SingleBodyAccelerationMap::const_iterator accelerationModelIterator
             = acceleratedBodyIterator->second.begin( );
             accelerationModelIterator != acceleratedBodyIterator->second.end( );
             accelerationModelIterator++ )
        {
            singleAccelerationUpdateNeeds.clear( );
            for( unsigned int i = 0; i < accelerationModelIterator->second.size( ); i++ )
            {
                AvailableAcceleration currentAccelerationModelType =
                        getAccelerationModelType( accelerationModelIterator->second.at( i ) );

                // Add translational state of both bodies to update list for current acceleration model.
                if( translationalAccelerationModels.count( accelerationModelIterator->first ) == 0 )
                {
                    singleAccelerationUpdateNeeds[ body_translational_state_update ].push_back(
                                accelerationModelIterator->first );
                }

                // Check acceleration model type and change environment update list accordingly.
                switch( currentAccelerationModelType )
                {
                case point_mass_gravity:
                    break;
                case third_body_point_mass_gravity:
                {
                    std::shared_ptr< gravitation::ThirdBodyCentralGravityAcceleration >
                            thirdBodyAcceleration = std::dynamic_pointer_cast<
                            gravitation::ThirdBodyCentralGravityAcceleration >(
                                accelerationModelIterator->second.at( i ) );
                    if( thirdBodyAcceleration != nullptr && translationalAccelerationModels.count(
                                thirdBodyAcceleration->getCentralBodyName( ) ) == 0 )
                    {
                        if( translationalAccelerationModels.count( thirdBodyAcceleration->getCentralBodyName( ) ) == 0 )
                        {
                            singleAccelerationUpdateNeeds[ body_translational_state_update ].push_back(
                                        thirdBodyAcceleration->getCentralBodyName( ) );
                        }
                    }
                    else if( thirdBodyAcceleration == nullptr )
                    {
                        throw std::runtime_error(
                                    std::string( "Error, incompatible input (ThirdBodyCentralGravityAcceleration) to" )
                                    + std::string(  "createTranslationalEquationsOfMotionEnvironmentUpdaterSettings" ) );
                    }
                    break;
                }
                case aerodynamic:
                    singleAccelerationUpdateNeeds[ body_rotational_state_update ].push_back(
                                accelerationModelIterator->first );
                    singleAccelerationUpdateNeeds[ vehicle_flight_conditions_update ].push_back(
                                acceleratedBodyIterator->first );
                    singleAccelerationUpdateNeeds[ body_mass_update ].push_back(
                                acceleratedBodyIterator->first );
                    break;
                case radiation_pressure:
                {
                    const auto sourceName = accelerationModelIterator->first;
                    const auto targetName = acceleratedBodyIterator->first;

                    singleAccelerationUpdateNeeds[ body_mass_update ].push_back(targetName);
                    singleAccelerationUpdateNeeds[ radiation_source_model_update ].push_back(sourceName);
                    singleAccelerationUpdateNeeds[ radiation_pressure_target_model_update ].push_back(targetName);

                    auto radiationPressureAcceleration =
                            std::dynamic_pointer_cast<electromagnetism::RadiationPressureAcceleration>(
                                    accelerationModelIterator->second.at(i));

                    // Update original source in case of paneled source
                    auto paneledRadiationSourceModel =
                            std::dynamic_pointer_cast<electromagnetism::PaneledRadiationSourceModel>(
                                    radiationPressureAcceleration->getSourceModel());
                    if (paneledRadiationSourceModel != nullptr) {
                        // Only paneled source, not point source needs rotational state and has original source
                        singleAccelerationUpdateNeeds[ body_rotational_state_update ].push_back(sourceName);

                        for (const auto& originalSourceName :
                                paneledRadiationSourceModel->getSourcePanelRadiosityModelUpdater()->getOriginalSourceBodyNames())
                        {
                            singleAccelerationUpdateNeeds[ radiation_source_model_update ].push_back(originalSourceName);
                            singleAccelerationUpdateNeeds[ body_translational_state_update ].push_back(originalSourceName);
                            // No original source rotational state update necessary because only isotropic point sources
                            // are supported, which are rotation-invariant
                        }

                        for (const auto& bodyName :
                                paneledRadiationSourceModel->getSourcePanelRadiosityModelUpdater()->getOriginalSourceToSourceOccultingBodyNames())
                        {
                            singleAccelerationUpdateNeeds[ body_translational_state_update ].push_back(bodyName);
                        }
                    }

                    // Update target rotation in case of paneled target
                    auto paneledRadiationPressureTargetModel =
                            std::dynamic_pointer_cast<electromagnetism::PaneledRadiationPressureTargetModel>(
                                    radiationPressureAcceleration->getTargetModel());
                    if (paneledRadiationPressureTargetModel != nullptr)
                    {
                        // Only paneled target, not cannonball target needs rotational state
                        singleAccelerationUpdateNeeds[ body_rotational_state_update ].push_back(targetName );

                        // Positions of bodies tracked by panel normal need to be updated
                        for( unsigned int j = 0; j < paneledRadiationPressureTargetModel->getBodyFixedPanels( ).size( ); j++ )
                        {
                            std::string bodyTrackedByPanel = paneledRadiationPressureTargetModel->getBodyFixedPanels( ).at( j )->getTrackedBody( );
                            if( bodyTrackedByPanel != bodyTrackedByPanel )
                            {
                                singleAccelerationUpdateNeeds[ body_translational_state_update ].push_back( bodyTrackedByPanel );
                            }
                        }

                        if( paneledRadiationPressureTargetModel->getSegmentFixedPanels( ).size( ) > 0 )
                        {
                            singleAccelerationUpdateNeeds[ body_segment_orientation_update ].push_back(targetName);
                        }
                        for( auto it : paneledRadiationPressureTargetModel->getSegmentFixedPanels( ) )
                        {
                            for( unsigned int j = 0; j < it.second.size( ); j++ )
                            {
                                std::string bodyTrackedByPanel = it.second.at( j )->getTrackedBody( );
                                if( bodyTrackedByPanel != bodyTrackedByPanel )
                                {
                                    singleAccelerationUpdateNeeds[ body_translational_state_update ].push_back(bodyTrackedByPanel);
                                }
                            }
                        }
                    }

                    // Update occulting body positions
                    auto sourceToTargetOccultingBodyNames =
                            radiationPressureAcceleration->getSourceToTargetOccultationModel()->getOccultingBodyNames();
                    for (auto bodyName : sourceToTargetOccultingBodyNames)
                    {
                        singleAccelerationUpdateNeeds[ body_translational_state_update ].push_back(bodyName);
                    }

                    break;
                }
                case spherical_harmonic_gravity:
                    singleAccelerationUpdateNeeds[ body_rotational_state_update ].push_back(
                                accelerationModelIterator->first );
                    singleAccelerationUpdateNeeds[ spherical_harmonic_gravity_field_update ].
                            push_back( accelerationModelIterator->first );
                    break;
                case mutual_spherical_harmonic_gravity:
                    singleAccelerationUpdateNeeds[ body_rotational_state_update ].push_back(
                                accelerationModelIterator->first );
                    singleAccelerationUpdateNeeds[ spherical_harmonic_gravity_field_update ].push_back(
                                accelerationModelIterator->first );
                    singleAccelerationUpdateNeeds[ body_rotational_state_update ].push_back(
                                acceleratedBodyIterator->first );
                    singleAccelerationUpdateNeeds[ spherical_harmonic_gravity_field_update ].push_back(
                                acceleratedBodyIterator->first );
                    break;
                case polyhedron_gravity:
                    singleAccelerationUpdateNeeds[ body_rotational_state_update ].push_back(
                            accelerationModelIterator->first );
                    break;
                case ring_gravity:
                    singleAccelerationUpdateNeeds[ body_rotational_state_update ].push_back(
                            accelerationModelIterator->first );
                    break;
                case third_body_spherical_harmonic_gravity:
                {
                    singleAccelerationUpdateNeeds[ body_rotational_state_update ].push_back(
                                accelerationModelIterator->first );
                    singleAccelerationUpdateNeeds[ spherical_harmonic_gravity_field_update ].
                            push_back( accelerationModelIterator->first );

                    std::shared_ptr<
                            gravitation::ThirdBodySphericalHarmonicsGravitationalAccelerationModel >
                            thirdBodyAcceleration = std::dynamic_pointer_cast<
                            gravitation::ThirdBodySphericalHarmonicsGravitationalAccelerationModel >(
                                accelerationModelIterator->second.at( i ) );;
                    if( thirdBodyAcceleration != nullptr && translationalAccelerationModels.count(
                                thirdBodyAcceleration->getCentralBodyName( ) ) == 0  )
                    {
                        singleAccelerationUpdateNeeds[ body_translational_state_update ].push_back(
                                    thirdBodyAcceleration->getCentralBodyName( ) );
                    }
                    else if( thirdBodyAcceleration == nullptr )
                    {
                        throw std::runtime_error(
                                    std::string( "Error, incompatible input (ThirdBodySphericalHarmonicsGravitational" )
                                    + std::string( "AccelerationModel) to createTranslationalEquationsOfMotion ")
                                    + std::string( "EnvironmentUpdaterSettings" ) );
                    }
                    break;
                }
                case third_body_mutual_spherical_harmonic_gravity:
                {
                    singleAccelerationUpdateNeeds[ body_rotational_state_update ].push_back(
                                accelerationModelIterator->first );
                    singleAccelerationUpdateNeeds[ spherical_harmonic_gravity_field_update ].push_back(
                                accelerationModelIterator->first );
                    singleAccelerationUpdateNeeds[ body_rotational_state_update ].push_back(
                                acceleratedBodyIterator->first );
                    singleAccelerationUpdateNeeds[ spherical_harmonic_gravity_field_update ].push_back(
                                acceleratedBodyIterator->first );

                    std::shared_ptr< gravitation::ThirdBodyMutualSphericalHarmonicsGravitationalAccelerationModel >
                            thirdBodyAcceleration = std::dynamic_pointer_cast<
                            gravitation::ThirdBodyMutualSphericalHarmonicsGravitationalAccelerationModel >(
                                accelerationModelIterator->second.at( i ) );
                    if( thirdBodyAcceleration != nullptr && translationalAccelerationModels.count(
                                thirdBodyAcceleration->getCentralBodyName( ) ) == 0  )
                    {
                        singleAccelerationUpdateNeeds[ body_translational_state_update ].push_back(
                                    thirdBodyAcceleration->getCentralBodyName( ) );
                        singleAccelerationUpdateNeeds[ body_rotational_state_update ].push_back(
                                    thirdBodyAcceleration->getCentralBodyName( ) );
                        singleAccelerationUpdateNeeds[ spherical_harmonic_gravity_field_update ].push_back(
                                    thirdBodyAcceleration->getCentralBodyName( ) );
                    }
                    else if( thirdBodyAcceleration == nullptr )
                    {
                        throw std::runtime_error(
                                    std::string( "Error, incompatible input (ThirdBodyMutualSphericalHarmonicsGravitational" ) +
                                    std::string( "AccelerationModel) to createTranslationalEquationsOfMotion ") +
                                    std::string( "EnvironmentUpdaterSettings" ) );
                    }
                    break;
                }
                case third_body_polyhedron_gravity:
                {
                    singleAccelerationUpdateNeeds[ body_rotational_state_update ].push_back(
                        accelerationModelIterator->first );

                    std::shared_ptr< gravitation::ThirdBodyPolyhedronGravitationalAccelerationModel >
                        thirdBodyAcceleration = std::dynamic_pointer_cast<
                            gravitation::ThirdBodyPolyhedronGravitationalAccelerationModel >(
                                accelerationModelIterator->second.at( i ) );

                    if( thirdBodyAcceleration != nullptr && translationalAccelerationModels.count(
                                thirdBodyAcceleration->getCentralBodyName( ) ) == 0  )
                    {
                        singleAccelerationUpdateNeeds[ body_translational_state_update ].push_back(
                                    thirdBodyAcceleration->getCentralBodyName( ) );
                    }
                    else if( thirdBodyAcceleration == nullptr )
                    {
                        throw std::runtime_error(
                                    std::string( "Error, incompatible input (ThirdBodyPolyhedronGravitational" )
                                    + std::string( "AccelerationModel) to createTranslationalEquationsOfMotion ")
                                    + std::string( "EnvironmentUpdaterSettings" ) );
                    }
                    break;
                }
                case third_body_ring_gravity:
                {
                    singleAccelerationUpdateNeeds[ body_rotational_state_update ].push_back(
                        accelerationModelIterator->first );

                    std::shared_ptr< gravitation::ThirdBodyRingGravitationalAccelerationModel >
                        thirdBodyAcceleration = std::dynamic_pointer_cast<
                            gravitation::ThirdBodyRingGravitationalAccelerationModel >(
                                accelerationModelIterator->second.at( i ) );

                    if( thirdBodyAcceleration != nullptr && translationalAccelerationModels.count(
                                thirdBodyAcceleration->getCentralBodyName( ) ) == 0 )
                    {
                        singleAccelerationUpdateNeeds[ body_translational_state_update ].push_back(
                                    thirdBodyAcceleration->getCentralBodyName( ) );
                    }
                    else if( thirdBodyAcceleration == nullptr )
                    {
                        throw std::runtime_error(
                                "Error, incompatible input (ThirdBodyRingGravitationalAccelerationModel "
                                "to createTranslationalEquationsOfMotion EnvironmentUpdaterSettings" );
                    }
                    break;
                }
                case thrust_acceleration:
                {
                    std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > thrustModelUpdates =
                            std::dynamic_pointer_cast< propulsion::ThrustAcceleration >(
                                accelerationModelIterator->second.at( i ) )->getRequiredModelUpdates( );
                    addEnvironmentUpdates( singleAccelerationUpdateNeeds, thrustModelUpdates );

                    singleAccelerationUpdateNeeds[ body_mass_update ].push_back(
                                acceleratedBodyIterator->first );

                    break;
                }
                case relativistic_correction_acceleration:
                {
                    std::shared_ptr< relativity::RelativisticAccelerationCorrection >
                            accelerationCorrection = std::dynamic_pointer_cast< relativity::RelativisticAccelerationCorrection >(
                                accelerationModelIterator->second.at( i ) );
                    if( accelerationCorrection->getCalculateDeSitterCorrection( ) )
                    {
                        std::string primaryBody = accelerationCorrection->getPrimaryBodyName( );
                        if( translationalAccelerationModels.count(primaryBody ) == 0 )
                        {
                            singleAccelerationUpdateNeeds[ body_translational_state_update ].push_back(
                                        primaryBody );
                        }
                    }
                    break;
                }
                case direct_tidal_dissipation_in_central_body_acceleration:
                    singleAccelerationUpdateNeeds[ body_rotational_state_update ].push_back(
                                accelerationModelIterator->first );
                    singleAccelerationUpdateNeeds[ spherical_harmonic_gravity_field_update ].
                            push_back( accelerationModelIterator->first );
                    break;
                case direct_tidal_dissipation_in_orbiting_body_acceleration:
                    singleAccelerationUpdateNeeds[ body_rotational_state_update ].push_back(
                                accelerationModelIterator->first );
                    singleAccelerationUpdateNeeds[ spherical_harmonic_gravity_field_update ].
                            push_back( accelerationModelIterator->first );
                    break;
                case empirical_acceleration:
                    break;
                case momentum_wheel_desaturation_acceleration:
                    break;
                case yarkovsky_acceleration:
                        break;
                case custom_acceleration:
                    break;
                case einstein_infeld_hoffmann_acceleration:
                {
                    std::shared_ptr< relativity::EinsteinInfeldHoffmannAcceleration > eihAcceleration =
                        std::dynamic_pointer_cast< relativity::EinsteinInfeldHoffmannAcceleration >(
                            accelerationModelIterator->second.at( i ) );
                    if( eihAcceleration == nullptr )
                    {
                        throw std::runtime_error( "Error when getting environment updates for EIH acceleration, acceleration object is incompatible" );
                    }
                    std::vector< std::string > bodiesExertingEihAcceleration = eihAcceleration->getBodiesExertingAcceleration( );
                    for( unsigned int j = 0; j < bodiesExertingEihAcceleration.size( ); j++ )
                    {
                        singleAccelerationUpdateNeeds[ body_translational_state_update ].push_back( bodiesExertingEihAcceleration.at( j ) );
                    }
                    break;
                }
                default:
                    throw std::runtime_error( std::string( "Error when setting acceleration model update needs, model type not recognized: " ) +
                                              std::to_string( currentAccelerationModelType ) );
                    break;
                }

            }

            // Add requested updates of current acceleration model to
            // full list of environment updates.
            addEnvironmentUpdates( environmentModelsToUpdate, singleAccelerationUpdateNeeds );
        }
    }

    return environmentModelsToUpdate;
}

//! Get list of required environment model update settings from mass rate models.
std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >
createMassPropagationEnvironmentUpdaterSettings(
        const std::map< std::string, std::vector< std::shared_ptr< basic_astrodynamics::MassRateModel > > > massRateModels,
        const simulation_setup::SystemOfBodies& bodies )
{
    using namespace basic_astrodynamics;
    using namespace propagators;

    std::map< propagators::EnvironmentModelsToUpdate,
            std::vector< std::string > > environmentModelsToUpdate;
    std::map< propagators::EnvironmentModelsToUpdate,
            std::vector< std::string > > singleRateModelUpdateNeeds;

    // Iterate over all bodies with mass rate model.
    for( std::map< std::string, std::vector< std::shared_ptr< MassRateModel > > >::const_iterator massRateModelIterator =
         massRateModels.begin( ); massRateModelIterator != massRateModels.end( ); massRateModelIterator++ )
    {
        for( unsigned int i = 0; i < massRateModelIterator->second.size( ); i++ )
        {
            singleRateModelUpdateNeeds.clear( );

            // Identify mass rate type and set required environment update settings.
            AvailableMassRateModels currentAccelerationModelType =
                    getMassRateModelType( massRateModelIterator->second.at( i ) );
            switch( currentAccelerationModelType )
            {
            case custom_mass_rate_model:
                break;
            case from_thrust_mass_rate_model:
                break;
            default:
                throw std::runtime_error( std::string( "Error when setting mass rate model update needs, model type not recognized: " ) +
                                          std::to_string( currentAccelerationModelType ) );

            }

            // Add requested updates of current acceleration model to
            // full list of environment updates.
            addEnvironmentUpdates( environmentModelsToUpdate, singleRateModelUpdateNeeds );
        }
    }

    return environmentModelsToUpdate;
}

//! Function to update environment to allow all required updates to be made
void checkAndModifyEnvironmentForDependentVariableSaving(
        const EnvironmentModelsToUpdate updateType,
        const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSaveSettings,
        const simulation_setup::SystemOfBodies& bodies )
{
    if( bodies.count( dependentVariableSaveSettings->associatedBody_ ) == 0 )
    {
        throw std::runtime_error( "Error when setting update for computation of dependent variable: "+
                                  getDependentVariableId( dependentVariableSaveSettings ) +
                                  ". Body <" + dependentVariableSaveSettings->associatedBody_ + "> does not exist" );
    }


    switch( updateType )
    {
    case vehicle_flight_conditions_update:
        if( bodies.at( dependentVariableSaveSettings->associatedBody_ )->getFlightConditions( ) == nullptr )
        {
            if( ( dependentVariableSaveSettings->secondaryBody_ != "" ) &&
                  ( dependentVariableSaveSettings->secondaryBody_ != "SSB" ) &&
                    bodies.count( dependentVariableSaveSettings->secondaryBody_ ) == 0 )
            {
                throw std::runtime_error( "Error when setting update for computation of dependent variable "+
                                          getDependentVariableId( dependentVariableSaveSettings ) +
                                          ". PROBLEM: Body <" + dependentVariableSaveSettings->secondaryBody_ + "> does not exist" );
            }
            simulation_setup::addFlightConditions(
                        bodies, dependentVariableSaveSettings->associatedBody_, dependentVariableSaveSettings->secondaryBody_ );
        }
        break;
    default:
        break;
    }
}

//! Function to create environment update settings for a single dependent variable
std::map< propagators::EnvironmentModelsToUpdate,
std::vector< std::string > > createEnvironmentUpdaterSettingsForDependentVariables(
        const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSaveSettings,
        const simulation_setup::SystemOfBodies& bodies )
{
    std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >  variablesToUpdate;
    switch( dependentVariableSaveSettings->dependentVariableType_ )
    {
    case mach_number_dependent_variable:
        variablesToUpdate[ vehicle_flight_conditions_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_rotational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        break;
    case altitude_dependent_variable:
        variablesToUpdate[ vehicle_flight_conditions_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_rotational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        break;
    case airspeed_dependent_variable:
        variablesToUpdate[ vehicle_flight_conditions_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_rotational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        break;
    case local_density_dependent_variable:
        variablesToUpdate[ vehicle_flight_conditions_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_rotational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        break;
    case relative_speed_dependent_variable:
        if( dependentVariableSaveSettings->associatedBody_ != "SSB" )
        {
            variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        }
        if( dependentVariableSaveSettings->secondaryBody_ != "SSB" )
        {
            variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        }
        break;
    case relative_position_dependent_variable:
        if( dependentVariableSaveSettings->associatedBody_ != "SSB" )
        {
            variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        }
        if( dependentVariableSaveSettings->secondaryBody_ != "SSB" )
        {
            variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        }
        break;
    case relative_distance_dependent_variable:
        if( dependentVariableSaveSettings->associatedBody_ != "SSB" )
        {
            variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        }
        if( dependentVariableSaveSettings->secondaryBody_ != "SSB" )
        {
            variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        }
        break;
    case relative_velocity_dependent_variable:
        if( dependentVariableSaveSettings->associatedBody_ != "SSB" )
        {
            variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        }
        if( dependentVariableSaveSettings->secondaryBody_ != "SSB" )
        {
            variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        }
        break;
    case total_acceleration_norm_dependent_variable:
        break;
    case single_acceleration_norm_dependent_variable:
        break;
    case total_acceleration_dependent_variable:
        break;
    case single_acceleration_dependent_variable:
        break;
    case aerodynamic_force_coefficients_dependent_variable:
        variablesToUpdate[ vehicle_flight_conditions_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_rotational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        break;
    case aerodynamic_moment_coefficients_dependent_variable:
        variablesToUpdate[ vehicle_flight_conditions_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_rotational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        break;
    case aerodynamic_control_surface_free_force_coefficients_dependent_variable:
        variablesToUpdate[ vehicle_flight_conditions_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_rotational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        break;
    case aerodynamic_control_surface_free_moment_coefficients_dependent_variable:
        variablesToUpdate[ vehicle_flight_conditions_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_rotational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        break;
    case aerodynamic_control_surface_force_coefficients_increment_dependent_variable:
        variablesToUpdate[ vehicle_flight_conditions_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_rotational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        break;
    case aerodynamic_control_surface_moment_coefficients_increment_dependent_variable:
        variablesToUpdate[ vehicle_flight_conditions_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_rotational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        break;
   case inertial_to_body_fixed_rotation_matrix_variable:
        variablesToUpdate[ body_rotational_state_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        break;
    case intermediate_aerodynamic_rotation_matrix_variable:
        variablesToUpdate[ vehicle_flight_conditions_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_rotational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        break;
    case relative_body_aerodynamic_orientation_angle_variable:
        variablesToUpdate[ vehicle_flight_conditions_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_rotational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        break;
    case body_fixed_airspeed_based_velocity_variable:
        variablesToUpdate[ vehicle_flight_conditions_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_rotational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        break;
    case total_aerodynamic_g_load_variable:
        variablesToUpdate[ vehicle_flight_conditions_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_rotational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        break;
    case stagnation_point_heat_flux_dependent_variable:
        variablesToUpdate[ vehicle_flight_conditions_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_rotational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        break;
    case local_temperature_dependent_variable:
        variablesToUpdate[ vehicle_flight_conditions_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_rotational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        break;
    case local_dynamic_pressure_dependent_variable:
        variablesToUpdate[ vehicle_flight_conditions_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_rotational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        break;
    case local_aerodynamic_heat_rate_dependent_variable:
        variablesToUpdate[ vehicle_flight_conditions_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_rotational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        break;
    case geodetic_latitude_dependent_variable:
        variablesToUpdate[ vehicle_flight_conditions_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_rotational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        break;
    case body_fixed_groundspeed_based_velocity_variable:
        variablesToUpdate[ vehicle_flight_conditions_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_rotational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        break;
    case total_mass_rate_dependent_variables:
        break;
    case total_torque_norm_dependent_variable:
        break;
    case single_torque_norm_dependent_variable:
        break;
    case total_torque_dependent_variable:
        break;
    case single_torque_dependent_variable:
        break;
    case keplerian_state_dependent_variable:
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        break;
    case modified_equinocial_state_dependent_variable:
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        break;
    case spherical_harmonic_acceleration_terms_dependent_variable:
        break;
    case spherical_harmonic_acceleration_norm_terms_dependent_variable:
        break;
    case body_fixed_relative_cartesian_position:
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        variablesToUpdate[ body_rotational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        break;
    case body_fixed_relative_spherical_position:
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        variablesToUpdate[ body_rotational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        break;
    case euler_angles_to_body_fixed_313:
        variablesToUpdate[ body_rotational_state_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        break;
    case tnw_to_inertial_frame_rotation_dependent_variable:
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        break;
    case rsw_to_inertial_frame_rotation_dependent_variable:
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        break;
    case control_surface_deflection_dependent_variable:
        variablesToUpdate[ vehicle_flight_conditions_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        break;
//    case radiation_pressure_dependent_variable:
//        variablesToUpdate[ radiation_pressure_interface_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
//        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
//        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
//        break;
    case periapsis_altitude_dependent_variable:
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        break;
    case apoapsis_altitude_dependent_variable:
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        break;
    case total_gravity_field_variation_acceleration:
        break;
    case single_gravity_field_variation_acceleration:
        break;
    case single_gravity_field_variation_acceleration_terms:
        break;
    case acceleration_partial_wrt_body_translational_state:
        break;
    case total_acceleration_partial_wrt_body_translational_state:
        break;
    case total_spherical_harmonic_cosine_coefficient_variation:
        variablesToUpdate[ spherical_harmonic_gravity_field_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        break;
    case total_spherical_harmonic_sine_coefficient_variation:
        variablesToUpdate[ spherical_harmonic_gravity_field_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        break;
    case current_body_mass_dependent_variable:
        variablesToUpdate[ body_mass_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        break;
    case radiation_pressure_coefficient_dependent_variable:
        variablesToUpdate[ radiation_pressure_target_model_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        break;
    case custom_dependent_variable:
        break;
    case gravity_field_potential_dependent_variable:
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        break;
    case gravity_field_laplacian_of_potential_dependent_variable:
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->secondaryBody_ );
        break;
    case minimum_constellation_distance:
    {
        std::shared_ptr< MinimumConstellationDistanceDependentVariableSaveSettings > minimumDistanceDependentVariable =
                std::dynamic_pointer_cast< MinimumConstellationDistanceDependentVariableSaveSettings >( dependentVariableSaveSettings );
        if( minimumDistanceDependentVariable == nullptr )
        {
            throw std::runtime_error( "Error when getting environment updates for dependent variables, type for minimum_constellation_distance not consistent" );
        }
        else
        {
            variablesToUpdate[ body_translational_state_update ].push_back( minimumDistanceDependentVariable->associatedBody_ );
            for( unsigned int i = 0; i < minimumDistanceDependentVariable->bodiesToCheck_.size( ); i++ )
            {
                variablesToUpdate[ body_translational_state_update ].push_back( minimumDistanceDependentVariable->bodiesToCheck_.at( i ) );
            }
        }
        break;
    }
    case minimum_constellation_ground_station_distance:
    {
        std::shared_ptr< MinimumConstellationStationDistanceDependentVariableSaveSettings > minimumDistanceDependentVariable =
                std::dynamic_pointer_cast< MinimumConstellationStationDistanceDependentVariableSaveSettings >( dependentVariableSaveSettings );
        if( minimumDistanceDependentVariable == nullptr )
        {
            throw std::runtime_error( "Error when getting environment updates for dependent variables, type for minimum_constellation_ground_station_distance not consistent" );
        }
        else
        {
            variablesToUpdate[ body_translational_state_update ].push_back( minimumDistanceDependentVariable->associatedBody_ );
            variablesToUpdate[ body_rotational_state_update ].push_back( minimumDistanceDependentVariable->associatedBody_ );

            for( unsigned int i = 0; i < minimumDistanceDependentVariable->bodiesToCheck_.size( ); i++ )
            {
                variablesToUpdate[ body_translational_state_update ].push_back( minimumDistanceDependentVariable->bodiesToCheck_.at( i ) );
            }
        }
        break;
    }
    case received_irradiance:
    case received_fraction:
    case visible_and_emitting_source_panel_count:
    case visible_source_area:
    {
        // Update target position, this should also update the source position
        variablesToUpdate[ body_translational_state_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        break;
    }
    case body_center_of_mass:
        variablesToUpdate[ body_mass_distribution_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        break;
    case body_inertia_tensor:
        variablesToUpdate[ body_mass_distribution_update ].push_back( dependentVariableSaveSettings->associatedBody_ );
        break;
    default:
        throw std::runtime_error( "Error when getting environment updates for dependent variables, parameter " +
                                  std::to_string( dependentVariableSaveSettings->dependentVariableType_ ) + " not found." );
    }

    if( variablesToUpdate.count( vehicle_flight_conditions_update ) > 0 )
    {
        checkAndModifyEnvironmentForDependentVariableSaving(
                    vehicle_flight_conditions_update, dependentVariableSaveSettings, bodies );
    }

    return variablesToUpdate;
}

//! Create environment update settings for dependent variables
std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > createEnvironmentUpdaterSettings(
        const std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > >& dependentVariableList,
        const simulation_setup::SystemOfBodies& bodies )
{
    std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > environmentModelsToUpdate;

    if( dependentVariableList.size( ) > 0 )
    {
        for( unsigned int i = 0; i < dependentVariableList.size( ); i++ )
        {
            std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > currentEnvironmentModelsToUpdate =
                    createEnvironmentUpdaterSettingsForDependentVariables( dependentVariableList.at( i ), bodies );
            addEnvironmentUpdates( environmentModelsToUpdate, currentEnvironmentModelsToUpdate );
        }
    }
    return environmentModelsToUpdate;
}

//! Create environment update settings for termination settings
std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > createEnvironmentUpdaterSettings(
        const std::shared_ptr< PropagationTerminationSettings > terminationSettings,
        const simulation_setup::SystemOfBodies& bodies )
{
    std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > environmentModelsToUpdate;
    switch( terminationSettings->terminationType_ )
    {
    case time_stopping_condition:
        break;
    case cpu_time_stopping_condition:
        break;
    case dependent_variable_stopping_condition:
    {
        std::shared_ptr< PropagationDependentVariableTerminationSettings > dependentVariableTerminationSettings =
                std::dynamic_pointer_cast< PropagationDependentVariableTerminationSettings >( terminationSettings );
        std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >
                environmentModelsToUpdateForSingleTerminationSetting =
                createEnvironmentUpdaterSettingsForDependentVariables(
                    dependentVariableTerminationSettings->dependentVariableSettings_, bodies );
        addEnvironmentUpdates( environmentModelsToUpdate, environmentModelsToUpdateForSingleTerminationSetting );
        break;
    }
    case hybrid_stopping_condition:
    {
        std::shared_ptr< PropagationHybridTerminationSettings > hybridTerminationSettings =
                std::dynamic_pointer_cast< PropagationHybridTerminationSettings >( terminationSettings );
        for( unsigned int i = 0; i < hybridTerminationSettings->terminationSettings_.size( ); i++ )
        {
            std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >
                    environmentModelsToUpdateForSingleTerminationSetting = createEnvironmentUpdaterSettings(
                        hybridTerminationSettings->terminationSettings_.at( i ), bodies );
            addEnvironmentUpdates( environmentModelsToUpdate, environmentModelsToUpdateForSingleTerminationSetting );
        }
        break;
    }
    case custom_stopping_condition:
        break;
    case non_sequential_stopping_condition:
        break;
    default:
        throw std::runtime_error( "Error when creating environment updater settings for termination conditions, type not found" );
    }

    return environmentModelsToUpdate;
}

//! Function to create 'brute-force' update settings, in which each environment model is updated.
std::map< propagators::EnvironmentModelsToUpdate,
std::vector< std::string > > createFullEnvironmentUpdaterSettings(
        const simulation_setup::SystemOfBodies& bodies )
{
    using namespace basic_astrodynamics;
    using namespace electromagnetism;
    using namespace gravitation;
    using namespace propagators;

    std::map< propagators::EnvironmentModelsToUpdate,
            std::vector< std::string > > environmentModelsToUpdate;
    std::map< propagators::EnvironmentModelsToUpdate,
            std::vector< std::string > > singleAccelerationUpdateNeeds;

    // Iterate over all bodies.
    for( auto bodyIterator : bodies.getMap( )  )
    {
        singleAccelerationUpdateNeeds.clear( );

        // Check if current body is a vehicle.

        // Check if current body has flight conditions set.
        if( bodyIterator.second ->getFlightConditions( ) != nullptr )
        {
            // If vehicle has flight conditions, add flight conditions update function to update list.
            singleAccelerationUpdateNeeds[ vehicle_flight_conditions_update ].
                    push_back( bodyIterator.first );
        }


        // If body has radiation source model, add its update to update list
        if( bodyIterator.second->getRadiationSourceModel() != nullptr )
        {
            singleAccelerationUpdateNeeds[ radiation_source_model_update ].
                    push_back( bodyIterator.first );
        }

        // If body has radiation pressure target model, add its update to update list
        if( bodyIterator.second->getRadiationPressureTargetModel() != nullptr )
        {
            singleAccelerationUpdateNeeds[ radiation_pressure_target_model_update ].
                    push_back( bodyIterator.first );
        }

        // If body has rotation model, update rotational state in each time step.;
        if( ( bodyIterator.second->getRotationalEphemeris( ) != nullptr )
//                ||
//                ( bodyIterator.second->getDependentOrientationCalculator( ) != nullptr )
                )
        {
            singleAccelerationUpdateNeeds[ body_rotational_state_update ].
                    push_back( bodyIterator.first );
        }


        std::shared_ptr< TimeDependentSphericalHarmonicsGravityField > gravityField =
                std::dynamic_pointer_cast< TimeDependentSphericalHarmonicsGravityField >
                ( bodyIterator.second->getGravityFieldModel( ) );
        if( gravityField != nullptr )
        {
            singleAccelerationUpdateNeeds[ spherical_harmonic_gravity_field_update ].
                    push_back( bodyIterator.first );
        }

        singleAccelerationUpdateNeeds[ body_mass_update ].push_back( bodyIterator.first );

        // Check whether requested updates are possible.
        checkValidityOfRequiredEnvironmentUpdates( singleAccelerationUpdateNeeds, bodies );

        // Add requested updates of current acceleration model to full list of environment updates.
        addEnvironmentUpdates( environmentModelsToUpdate, singleAccelerationUpdateNeeds );
    }
    return environmentModelsToUpdate;
}

} // namespace propagators

} // namespace tudat
