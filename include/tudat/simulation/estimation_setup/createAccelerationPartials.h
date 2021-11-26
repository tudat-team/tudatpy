/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEACCELERATIONPARTIALS_H
#define TUDAT_CREATEACCELERATIONPARTIALS_H

#include <memory.h>

#include "tudat/astro/basic_astro/accelerationModel.h"

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/astro/orbit_determination/acceleration_partials/accelerationPartial.h"
#include "tudat/astro/orbit_determination/acceleration_partials/centralGravityAccelerationPartial.h"
#include "tudat/astro/orbit_determination/acceleration_partials/radiationPressureAccelerationPartial.h"
#include "tudat/astro/orbit_determination/acceleration_partials/thirdBodyGravityPartial.h"
#include "tudat/astro/orbit_determination/acceleration_partials/relativisticAccelerationPartial.h"
#include "tudat/astro/orbit_determination/acceleration_partials/sphericalHarmonicAccelerationPartial.h"
#include "tudat/astro/orbit_determination/acceleration_partials/aerodynamicAccelerationPartial.h"
#include "tudat/astro/orbit_determination/acceleration_partials/mutualSphericalHarmonicGravityPartial.h"
#include "tudat/astro/orbit_determination/acceleration_partials/empiricalAccelerationPartial.h"
#include "tudat/astro/orbit_determination/acceleration_partials/directTidalDissipationAccelerationPartial.h"
#include "tudat/astro/orbit_determination/acceleration_partials/panelledRadiationPressureAccelerationPartial.h"
#include "tudat/astro/orbit_determination/acceleration_partials/thrustAccelerationPartial.h"
#include "tudat/astro/orbit_determination/observation_partials/rotationMatrixPartial.h"
#include "tudat/simulation/estimation_setup/createCartesianStatePartials.h"
#include "tudat/astro/basic_astro/accelerationModelTypes.h"
#include "tudat/astro/orbit_determination/acceleration_partials/tidalLoveNumberPartialInterface.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to create a list of objects that can be used to compute partials of tidal gravity field variations
/*!
 * Function to create a list of objects that can be used to compute partials of tidal gravity field variations
 * \param bodies List of all body objects
 * \param acceleratingBodyName Name of body for which tidal gravity field variation objects are to be created
 * \return List of tidal gravity field variation objects, one for each such field variation object of bodyacceleratingBodyName
 */
std::vector< std::shared_ptr< orbit_determination::TidalLoveNumberPartialInterface > > createTidalLoveNumberInterfaces(
        const SystemOfBodies& bodies,
        const std::string& acceleratingBodyName );

//! Function to create a single acceleration partial derivative object.
/*!
 *  Function to create a single acceleration partial derivative object.
 *  \param accelerationModel Acceleration model for which a partial derivative is to be computed.
 *  \param acceleratedBody Pair of name and object of body undergoing acceleration
 *  \param acceleratingBody Pair of name and object of body exerting acceleration
 *  \param bodies List of all body objects
 *  \param parametersToEstimate List of parameters that are to be estimated. Empty by default, only required for selected
 *  types of partials (e.g. spherical harmonic acceleration w.r.t. rotational parameters).
 *  \return Single acceleration partial derivative object.
 */
template< typename InitialStateParameterType = double >
std::shared_ptr< acceleration_partials::AccelerationPartial > createAnalyticalAccelerationPartial(
        std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > accelerationModel,
        const std::pair< std::string, std::shared_ptr< simulation_setup::Body > > acceleratedBody,
        const std::pair< std::string, std::shared_ptr< simulation_setup::Body > > acceleratingBody,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > >
        parametersToEstimate =
        std::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > >( ) )
{
    using namespace gravitation;
    using namespace basic_astrodynamics;
    using namespace electromagnetism;
    using namespace aerodynamics;
    using namespace acceleration_partials;

    std::shared_ptr< acceleration_partials::AccelerationPartial > accelerationPartial = nullptr;

    // Identify current acceleration model type
    AvailableAcceleration accelerationType = getAccelerationModelType( accelerationModel );
    switch( accelerationType )
    {
    case point_mass_gravity:

        // Check if identifier is consistent with type.
        if( std::dynamic_pointer_cast< CentralGravitationalAccelerationModel3d >( accelerationModel ) == nullptr )
        {
            throw std::runtime_error( "Acceleration class type does not match acceleration type (point_mass_gravity) "
                                      "when making acceleration partial." );
        }
        else
        {
            // Create partial-calculating object.
            accelerationPartial = std::make_shared< CentralGravitationPartial >
                    ( std::dynamic_pointer_cast< CentralGravitationalAccelerationModel3d >( accelerationModel ),
                      acceleratedBody.first, acceleratingBody.first );
        }
        break;
    case relativistic_correction_acceleration:

        // Check if identifier is consistent with type.
        if( std::dynamic_pointer_cast< relativity::RelativisticAccelerationCorrection >( accelerationModel ) == nullptr )
        {
            throw std::runtime_error( "Acceleration class type does not match acceleration type "
                                      "(relativistic_correction_acceleration) when making acceleration partial." );
        }
        else
        {
            // Create partial-calculating object.
            accelerationPartial = std::make_shared< RelativisticAccelerationPartial  >
                    ( std::dynamic_pointer_cast< relativity::RelativisticAccelerationCorrection >( accelerationModel ),
                      acceleratedBody.first, acceleratingBody.first );
        }
        break;
    case direct_tidal_dissipation_in_central_body_acceleration:
    {
        // Check if identifier is consistent with type.
        if( std::dynamic_pointer_cast< gravitation::DirectTidalDissipationAcceleration >( accelerationModel ) == nullptr )
        {
            throw std::runtime_error( "Acceleration class type does not match acceleration type "
                                      "(direct_tidal_dissipation_acceleration) when making acceleration partial." );
        }
        else
        {
            // Create partial-calculating object.
            accelerationPartial = std::make_shared< DirectTidalDissipationAccelerationPartial  >
                    ( std::dynamic_pointer_cast< gravitation::DirectTidalDissipationAcceleration >( accelerationModel ),
                      acceleratedBody.first, acceleratingBody.first );
        }
        break;
    }
    case direct_tidal_dissipation_in_orbiting_body_acceleration:
    {
        // Check if identifier is consistent with type.
        if( std::dynamic_pointer_cast< gravitation::DirectTidalDissipationAcceleration >( accelerationModel ) == nullptr )
        {
            throw std::runtime_error( "Acceleration class type does not match acceleration type "
                                      "(direct_tidal_dissipation_acceleration) when making acceleration partial." );
        }
        else
        {
            // Create partial-calculating object.
            accelerationPartial = std::make_shared< DirectTidalDissipationAccelerationPartial  >
                    ( std::dynamic_pointer_cast< gravitation::DirectTidalDissipationAcceleration >( accelerationModel ),
                      acceleratedBody.first, acceleratingBody.first );
        }
        break;
    }
    case third_body_point_mass_gravity:
        // Check if identifier is consistent with type.
        if( std::dynamic_pointer_cast< ThirdBodyCentralGravityAcceleration >( accelerationModel ) == nullptr )
        {
            throw std::runtime_error( "Acceleration class type does not match acceleration type (third_body_point_mass_gravity) "
                                      "when making acceleration partial." );
        }
        else
        {
            std::shared_ptr< ThirdBodyCentralGravityAcceleration > thirdBodyAccelerationModel  =
                    std::dynamic_pointer_cast< ThirdBodyCentralGravityAcceleration >( accelerationModel );

            // Create partials for constituent central gravity accelerations
            std::shared_ptr< CentralGravitationPartial > accelerationPartialForBodyUndergoingAcceleration =
                    std::dynamic_pointer_cast< CentralGravitationPartial >(
                        createAnalyticalAccelerationPartial(
                            thirdBodyAccelerationModel->getAccelerationModelForBodyUndergoingAcceleration( ),
                            acceleratedBody, acceleratingBody, bodies, parametersToEstimate ) );
            std::shared_ptr< CentralGravitationPartial > accelerationPartialForCentralBody =
                    std::dynamic_pointer_cast< CentralGravitationPartial >(
                        createAnalyticalAccelerationPartial(
                            thirdBodyAccelerationModel->getAccelerationModelForCentralBody( ),
                            std::make_pair( thirdBodyAccelerationModel->getCentralBodyName( ),
                                            bodies.at( thirdBodyAccelerationModel->getCentralBodyName( ) ) ),
                            acceleratingBody, bodies, parametersToEstimate ) );

            // Create partial-calculating object.
            accelerationPartial = std::make_shared< ThirdBodyGravityPartial< CentralGravitationPartial > >(
                        accelerationPartialForBodyUndergoingAcceleration,
                        accelerationPartialForCentralBody, acceleratedBody.first, acceleratingBody.first,
                        thirdBodyAccelerationModel->getCentralBodyName( ) );

        }
        break;
    case spherical_harmonic_gravity:
    {
        // Check if identifier is consistent with type.
        std::shared_ptr< SphericalHarmonicsGravitationalAccelerationModel > sphericalHarmonicAcceleration =
                std::dynamic_pointer_cast< SphericalHarmonicsGravitationalAccelerationModel >( accelerationModel );
        if( sphericalHarmonicAcceleration == nullptr )
        {
            throw std::runtime_error(
                        "Acceleration class type does not match acceleration type enum (spher. harm. grav.) set when making "
                        "acceleration partial." );
        }
        else
        {
            std::map< std::pair< estimatable_parameters::EstimatebleParametersEnum, std::string >,
                    std::shared_ptr< observation_partials::RotationMatrixPartial > >
                    rotationMatrixPartials = observation_partials::createRotationMatrixPartials(
                        parametersToEstimate, acceleratingBody.first, bodies );

            // If body has gravity field variations, create partial objects
            std::vector< std::shared_ptr< orbit_determination::TidalLoveNumberPartialInterface > >
                    currentBodyLoveNumberPartialInterfaces;
            if( acceleratingBody.second->getGravityFieldVariationSet( ) != nullptr )
            {
                currentBodyLoveNumberPartialInterfaces = createTidalLoveNumberInterfaces(
                            bodies, acceleratingBody.first );
            }

            // Create partial-calculating object.
            accelerationPartial = std::make_shared< SphericalHarmonicsGravityPartial >
                    ( acceleratedBody.first, acceleratingBody.first,
                      sphericalHarmonicAcceleration, rotationMatrixPartials, currentBodyLoveNumberPartialInterfaces );
        }
        break;
    }
    case third_body_spherical_harmonic_gravity:
        // Check if identifier is consistent with type.
        if( std::dynamic_pointer_cast< ThirdBodySphericalHarmonicsGravitationalAccelerationModel >( accelerationModel ) == nullptr )
        {
            throw std::runtime_error( "Acceleration class type does not match acceleration type "
                                      "(third_body_spherical_harmonic_gravity) when making acceleration partia.l" );
        }
        else
        {
            std::shared_ptr< ThirdBodySphericalHarmonicsGravitationalAccelerationModel > thirdBodyAccelerationModel  =
                    std::dynamic_pointer_cast< ThirdBodySphericalHarmonicsGravitationalAccelerationModel >(
                        accelerationModel );

            // Create partials for constituent central gravity accelerations
            std::shared_ptr< SphericalHarmonicsGravityPartial > accelerationPartialForBodyUndergoingAcceleration =
                    std::dynamic_pointer_cast< SphericalHarmonicsGravityPartial >(
                        createAnalyticalAccelerationPartial(
                            thirdBodyAccelerationModel->getAccelerationModelForBodyUndergoingAcceleration( ),
                            acceleratedBody, acceleratingBody, bodies, parametersToEstimate ) );
            std::shared_ptr< SphericalHarmonicsGravityPartial > accelerationPartialForCentralBody =
                    std::dynamic_pointer_cast< SphericalHarmonicsGravityPartial >(
                        createAnalyticalAccelerationPartial(
                            thirdBodyAccelerationModel->getAccelerationModelForCentralBody( ),
                            std::make_pair( thirdBodyAccelerationModel->getCentralBodyName( ),
                                            bodies.at( thirdBodyAccelerationModel->getCentralBodyName( ) ) ),
                            acceleratingBody, bodies, parametersToEstimate  ) );

            // Create partial-calculating object.
            accelerationPartial = std::make_shared< ThirdBodyGravityPartial< SphericalHarmonicsGravityPartial > >(
                        accelerationPartialForBodyUndergoingAcceleration,
                        accelerationPartialForCentralBody, acceleratedBody.first, acceleratingBody.first,
                        thirdBodyAccelerationModel->getCentralBodyName( ) );

        }
        break;
    case mutual_spherical_harmonic_gravity:
    {
        // Check if identifier is consistent with type.
        std::shared_ptr< MutualSphericalHarmonicsGravitationalAccelerationModel > mutualSphericalHarmonicAcceleration =
                std::dynamic_pointer_cast< MutualSphericalHarmonicsGravitationalAccelerationModel >( accelerationModel );
        if( mutualSphericalHarmonicAcceleration == nullptr )
        {
            throw std::runtime_error( "Acceleration class type does not match acceleration type enum (mut. spher. harm. grav.) "
                                      "set when making acceleration partial." );
        }
        else
        {
            std::shared_ptr< SphericalHarmonicsGravityPartial > accelerationPartialOfShExpansionOfBodyExertingAcceleration =
                    std::dynamic_pointer_cast< SphericalHarmonicsGravityPartial >(
                        createAnalyticalAccelerationPartial(
                            mutualSphericalHarmonicAcceleration->getAccelerationModelFromShExpansionOfBodyExertingAcceleration( ),
                            acceleratedBody, acceleratingBody,bodies, parametersToEstimate ) );
            std::shared_ptr< SphericalHarmonicsGravityPartial > accelerationPartialOfShExpansionOfBodyUndergoingAcceleration =
                    std::dynamic_pointer_cast< SphericalHarmonicsGravityPartial >(
                        createAnalyticalAccelerationPartial(
                            mutualSphericalHarmonicAcceleration->getAccelerationModelFromShExpansionOfBodyUndergoingAcceleration( ),
                            acceleratingBody, acceleratedBody, bodies, parametersToEstimate ) );
            accelerationPartial = std::make_shared< MutualSphericalHarmonicsGravityPartial >(
                        accelerationPartialOfShExpansionOfBodyExertingAcceleration,
                        accelerationPartialOfShExpansionOfBodyUndergoingAcceleration, acceleratedBody.first, acceleratingBody.first,
                        mutualSphericalHarmonicAcceleration->getUseCentralBodyFixedFrame( ) );
        }
        break;
    }
    case third_body_mutual_spherical_harmonic_gravity:
    {
        // Check if identifier is consistent with type.
        if( std::dynamic_pointer_cast< ThirdBodyMutualSphericalHarmonicsGravitationalAccelerationModel >( accelerationModel ) == nullptr )
        {
            throw std::runtime_error( "Acceleration class type does not match acceleration type "
                                      "(third_body_mutual_spherical_harmonic_gravity) enum set when making acceleration partial." );
        }
        else
        {
            std::shared_ptr< ThirdBodyMutualSphericalHarmonicsGravitationalAccelerationModel > thirdBodyAccelerationModel  =
                    std::dynamic_pointer_cast< ThirdBodyMutualSphericalHarmonicsGravitationalAccelerationModel >( accelerationModel );

            std::shared_ptr< MutualSphericalHarmonicsGravityPartial > accelerationPartialForBodyUndergoingAcceleration =
                    std::dynamic_pointer_cast< MutualSphericalHarmonicsGravityPartial >(
                        createAnalyticalAccelerationPartial(
                            thirdBodyAccelerationModel->getAccelerationModelForBodyUndergoingAcceleration( ),
                            acceleratedBody, acceleratingBody, bodies, parametersToEstimate  ) );
            std::shared_ptr< MutualSphericalHarmonicsGravityPartial > accelerationPartialForCentralBody =
                    std::dynamic_pointer_cast< MutualSphericalHarmonicsGravityPartial >(
                        createAnalyticalAccelerationPartial(
                            thirdBodyAccelerationModel->getAccelerationModelForCentralBody( ),
                            std::make_pair( thirdBodyAccelerationModel->getCentralBodyName( ),
                                            bodies.at( thirdBodyAccelerationModel->getCentralBodyName( ) ) ),
                            acceleratingBody, bodies, parametersToEstimate ) );
            accelerationPartial = std::make_shared< ThirdBodyGravityPartial< MutualSphericalHarmonicsGravityPartial > >(
                        accelerationPartialForBodyUndergoingAcceleration,
                        accelerationPartialForCentralBody, acceleratedBody.first, acceleratingBody.first,
                        thirdBodyAccelerationModel->getCentralBodyName( ) );

        }
        break;
    }
    case cannon_ball_radiation_pressure:
    {
        // Check if identifier is consistent with type.
        std::shared_ptr< CannonBallRadiationPressureAcceleration > radiationPressureAcceleration =
                std::dynamic_pointer_cast< CannonBallRadiationPressureAcceleration >( accelerationModel );
        if( radiationPressureAcceleration == nullptr )
        {
            throw std::runtime_error( "Acceleration class type does not match acceleration type (cannon_ball_radiation_pressure) "
                                      "when making acceleration partial." );
        }
        else
        {
            std::map< std::string, std::shared_ptr< RadiationPressureInterface > > radiationPressureInterfaces =
                    acceleratedBody.second->getRadiationPressureInterfaces( );

            if( radiationPressureInterfaces.count( acceleratingBody.first ) == 0 )
            {
                throw std::runtime_error( "No radiation pressure coefficient interface found when making acceleration partial." );
            }
            else
            {
                std::shared_ptr< RadiationPressureInterface > radiationPressureInterface =
                        radiationPressureInterfaces.at( acceleratingBody.first );

                // Create partial-calculating object.
                accelerationPartial = std::make_shared< CannonBallRadiationPressurePartial >
                        ( radiationPressureInterface, radiationPressureAcceleration->getMassFunction( ),
                          acceleratedBody.first, acceleratingBody.first );
            }
        }
        break;
    }
    case aerodynamic:
    {
        // Check if identifier is consistent with type.
        std::shared_ptr< AerodynamicAcceleration > aerodynamicAcceleration =
                std::dynamic_pointer_cast< AerodynamicAcceleration >( accelerationModel );
        if( aerodynamicAcceleration == nullptr )
        {
            throw std::runtime_error( "Acceleration class type does not match acceleration type (aerodynamic) when making "
                                      "acceleration partial." );
        }
        else
        {
            std::shared_ptr< AtmosphericFlightConditions > flightConditions =
                    std::dynamic_pointer_cast< AtmosphericFlightConditions >(
                        acceleratedBody.second->getFlightConditions( ) );

            if( flightConditions == nullptr )
            {
                throw std::runtime_error( "No flight conditions found when making acceleration partial." );
            }
            else
            {
                // Create partial-calculating object.
                accelerationPartial = std::make_shared< AerodynamicAccelerationPartial >
                        ( aerodynamicAcceleration,
                          flightConditions,
                          std::bind( &Body::getState, acceleratedBody.second ),
                          std::bind( &Body::setState, acceleratedBody.second, std::placeholders::_1 ),
                          acceleratedBody.first, acceleratingBody.first );
            }
        }
        break;
    }
    case empirical_acceleration:
    {
        // Check if identifier is consistent with type.
        std::shared_ptr< EmpiricalAcceleration > empiricalAcceleration =
                std::dynamic_pointer_cast< EmpiricalAcceleration >( accelerationModel );
        if( empiricalAcceleration == nullptr )
        {
            std::cerr << "Acceleration class type does not match acceleration type enum (rel. corr.) "
                         "set when making acceleration partial." << std::endl;

        }
        else
        {
            // Create partial-calculating object.
            accelerationPartial = std::make_shared< EmpiricalAccelerationPartial >( empiricalAcceleration,
                                                                                    acceleratedBody.first, acceleratingBody.first );
        }
        break;
    }
    case momentum_wheel_desaturation_acceleration:
    {
        // Check if identifier is consistent with type.
        std::shared_ptr< propulsion::MomentumWheelDesaturationThrustAcceleration > thrustAcceleration =
                std::dynamic_pointer_cast< propulsion::MomentumWheelDesaturationThrustAcceleration >( accelerationModel );
        if( thrustAcceleration == nullptr )
        {
            std::cerr << "Acceleration class type does not match acceleration type enum (mom. wheel desat.) "
                         "set when making acceleration partial." << std::endl;

        }
        else
        {
            // Create partial-calculating object.
            accelerationPartial = std::make_shared< MomentumWheelDesaturationPartial >(
                        thrustAcceleration, acceleratedBody.first );
        }
        break;
    }
    case custom_acceleration:
        std::cerr<<"Warning, custom acceleration partials implicitly set to zero - depending on thrust guidance model, this may provide biased results for variational equations"<<std::endl;
        break;
    case thrust_acceleration:
        std::cerr<<"Warning, thrust acceleration partials implicitly set to zero - depending on thrust guidance model, this may provide biased results for variational equations"<<std::endl;
        break;
    default:
        std::string errorMessage = "Acceleration model " + std::to_string( accelerationType ) +
                " not found when making acceleration partial";
        throw std::runtime_error( errorMessage );
        break;
    }

    return accelerationPartial;
}

//extern template std::shared_ptr< acceleration_partials::AccelerationPartial > createAnalyticalAccelerationPartial< double >(
//        std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > accelerationModel,
//        const std::pair< std::string, std::shared_ptr< simulation_setup::Body > > acceleratedBody,
//        const std::pair< std::string, std::shared_ptr< simulation_setup::Body > > acceleratingBody,
//        const simulation_setup::SystemOfBodies& bodies,
//        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > >
//        parametersToEstimate );
//#if( TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS )
//extern template std::shared_ptr< acceleration_partials::AccelerationPartial > createAnalyticalAccelerationPartial< long double >(
//        std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > accelerationModel,
//        const std::pair< std::string, std::shared_ptr< simulation_setup::Body > > acceleratedBody,
//        const std::pair< std::string, std::shared_ptr< simulation_setup::Body > > acceleratingBody,
//        const simulation_setup::SystemOfBodies& bodies,
//        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< long double > >
//        parametersToEstimate );
//#endif

//! This function creates acceleration partial objects for translational dynamics
/*!
 *  This function creates acceleration partial objects for translational dynamics from acceleration models and
 *  list of bodies' states of which derivatives are needed. The return type is an StateDerivativePartialsMap,
 *  a standardized type for communicating such lists of these objects.
 *  \param accelerationMap Map of maps containing list of acceleration models, identifying which acceleration acts on which
 *   body.
 *  \param bodies List of body objects constituting environment for calculations.
 *  \param parametersToEstimate List of parameters which are to be estimated.
 *  \return List of acceleration-partial-calculating objects in StateDerivativePartialsMap type.
 */
template< typename InitialStateParameterType >
orbit_determination::StateDerivativePartialsMap createAccelerationPartialsMap(
        const basic_astrodynamics::AccelerationMap& accelerationMap,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > >
        parametersToEstimate )
{
    // Declare return map.
    orbit_determination::StateDerivativePartialsMap accelerationPartialsList;
    std::map< std::string, std::map< std::string,
            std::vector< std::shared_ptr< acceleration_partials::AccelerationPartial > > > >
            accelerationPartialsMap;

    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter<
            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > initialDynamicalParameters =
            getListOfTranslationalStateParametersToEstimate( parametersToEstimate );
    accelerationPartialsList.resize( initialDynamicalParameters.size( ) );

    // Iterate over list of bodies of which the partials of the accelerations acting on them are required.
    int bodyCounter = 0;
    for( basic_astrodynamics::AccelerationMap::const_iterator accelerationIterator = accelerationMap.begin( );
         accelerationIterator != accelerationMap.end( ); accelerationIterator++ )
    {
        for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
        {
            if( initialDynamicalParameters.at( i )->getParameterName( ).second.first == accelerationIterator->first )
            {
                if( ( initialDynamicalParameters.at( i )->getParameterName( ).first == estimatable_parameters::initial_body_state ) ||
                        ( initialDynamicalParameters.at( i )->getParameterName( ).first == estimatable_parameters::arc_wise_initial_body_state ) )
                {
                    // Get object for body undergoing acceleration
                    const std::string acceleratedBody = accelerationIterator->first;
                    std::shared_ptr< simulation_setup::Body > acceleratedBodyObject = bodies.at( acceleratedBody );

                    // Retrieve list of accelerations acting on current body.
                    basic_astrodynamics::SingleBodyAccelerationMap accelerationVector =
                            accelerationMap.at( acceleratedBody );

                    // Declare list of acceleration partials of current body.
                    std::vector< std::shared_ptr< orbit_determination::StateDerivativePartial > > accelerationPartialVector;

                    // Iterate over all acceleration models and generate their partial-calculating objects.
                    for(  basic_astrodynamics::SingleBodyAccelerationMap::iterator
                          innerAccelerationIterator = accelerationVector.begin( );
                          innerAccelerationIterator != accelerationVector.end( ); innerAccelerationIterator++ )
                    {
                        // Get object for body exerting acceleration
                        std::string acceleratingBody = innerAccelerationIterator->first;
                        std::shared_ptr< simulation_setup::Body > acceleratingBodyObject;
                        if( acceleratingBody != "" )
                        {
                            acceleratingBodyObject = bodies.at( acceleratingBody );
                        }

                        for( unsigned int j = 0; j < innerAccelerationIterator->second.size( ); j++ )
                        {
                            // Create single partial object
                            std::shared_ptr< acceleration_partials::AccelerationPartial > currentAccelerationPartial =
                                    createAnalyticalAccelerationPartial(
                                        innerAccelerationIterator->second[ j ],
                                        std::make_pair( acceleratedBody, acceleratedBodyObject ),
                                        std::make_pair( acceleratingBody, acceleratingBodyObject ),
                                        bodies, parametersToEstimate );

                            if( currentAccelerationPartial != nullptr )
                            {
                                accelerationPartialVector.push_back( currentAccelerationPartial );
                                accelerationPartialsMap[ acceleratedBody ][ acceleratingBody ].push_back(
                                            currentAccelerationPartial );
                            }
                        }
                    }

                    // Add partials of current body's accelerations to vector.
                    accelerationPartialsList[ i ] = accelerationPartialVector;

                    bodyCounter++;
                }
            }
        }
    }
    return accelerationPartialsList;
}

//extern template orbit_determination::StateDerivativePartialsMap createAccelerationPartialsMap< double >(
//const basic_astrodynamics::AccelerationMap& accelerationMap,
//const simulation_setup::SystemOfBodies& bodies,
//const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > >
//parametersToEstimate );

//#if( TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS )
//extern template orbit_determination::StateDerivativePartialsMap createAccelerationPartialsMap< long double >(
//const basic_astrodynamics::AccelerationMap& accelerationMap,
//const simulation_setup::SystemOfBodies& bodies,
//const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< long double > >
//parametersToEstimate );
//#endif

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_CREATEACCELERATIONPARTIALS_H
