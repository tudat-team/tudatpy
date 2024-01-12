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
#include "tudat/astro/orbit_determination/acceleration_partials/polyhedronAccelerationPartial.h"
#include "tudat/astro/orbit_determination/acceleration_partials/aerodynamicAccelerationPartial.h"
#include "tudat/astro/orbit_determination/acceleration_partials/mutualSphericalHarmonicGravityPartial.h"
#include "tudat/astro/orbit_determination/acceleration_partials/empiricalAccelerationPartial.h"
#include "tudat/astro/orbit_determination/acceleration_partials/einsteinInfeldHoffmannPartials.h"
#include "tudat/astro/orbit_determination/acceleration_partials/directTidalDissipationAccelerationPartial.h"
#include "tudat/astro/orbit_determination/acceleration_partials/panelledRadiationPressureAccelerationPartial.h"
#include "tudat/astro/orbit_determination/acceleration_partials/thrustAccelerationPartial.h"
#include "tudat/astro/orbit_determination/acceleration_partials/yarkovskyAccelerationPartial.h"
#include "tudat/astro/orbit_determination/acceleration_partials/customAccelerationPartial.h"
#include "tudat/astro/orbit_determination/acceleration_partials/fullRadiationPressureAccelerationPartial.h"
#include "tudat/astro/orbit_determination/observation_partials/rotationMatrixPartial.h"
#include "tudat/simulation/estimation_setup/createCartesianStatePartials.h"
#include "tudat/astro/basic_astro/accelerationModelTypes.h"
#include "tudat/astro/orbit_determination/acceleration_partials/tidalLoveNumberPartialInterface.h"
#include "tudat/simulation/propagation_setup/environmentUpdater.h"

namespace tudat
{

namespace simulation_setup
{

template< typename ParameterScalarType >
std::shared_ptr< estimatable_parameters::CustomAccelerationPartialCalculator > createCustomAccelerationPartial(
    const std::shared_ptr< estimatable_parameters::CustomAccelerationPartialSettings > customPartialSettings,
    const std::shared_ptr< estimatable_parameters::EstimatableParameter< ParameterScalarType > > parameter,
    const SystemOfBodies& bodies )
{
    std::shared_ptr< estimatable_parameters::CustomAccelerationPartialCalculator > customPartialCalculator = nullptr;
    std::shared_ptr< estimatable_parameters::NumericalAccelerationPartialSettings > numericalCustomPartialSettings =
        std::dynamic_pointer_cast< estimatable_parameters::NumericalAccelerationPartialSettings >( customPartialSettings );
    std::shared_ptr< estimatable_parameters::AnalyticalAccelerationPartialSettings > analyticalCustomPartialSettings =
        std::dynamic_pointer_cast< estimatable_parameters::AnalyticalAccelerationPartialSettings >( customPartialSettings );

    if( numericalCustomPartialSettings!= nullptr )
    {
        if( parameter->getParameterName( ).first != estimatable_parameters::initial_body_state )
        {
            throw std::runtime_error( "Error, only initial cartesian state supported for custom numerical acceleration partial" );
        }
        std::string bodyName = parameter->getParameterName( ).second.first;
        if( bodies.count( bodyName ) == 0 )
        {
            throw std::runtime_error( "Error when creating custom numerical acceleration partial. Could not find body " + bodyName );
        }

        std::function< Eigen::Vector6d( ) > bodyStateGetFunction = std::bind( &Body::getState, bodies.at( bodyName ) );
        std::function< void( const Eigen::Vector6d& ) > bodyStateSetFunction = std::bind( &Body::setState, bodies.at( bodyName ), std::placeholders::_1 );

        std::function< void( const double ) > environmentUpdateFunction = std::bind(
            &propagators::EnvironmentUpdater< double, double >::updateEnvironment,
                std::make_shared< propagators::EnvironmentUpdater< double, double > >(
                    bodies, numericalCustomPartialSettings->environmentUpdateSettings_ ),
                std::placeholders::_1,
                std::unordered_map< propagators::IntegratedStateType, Eigen::VectorXd >( ),
                std::vector< propagators::IntegratedStateType >( ) );

        customPartialCalculator = std::make_shared< estimatable_parameters::NumericalAccelerationPartialWrtStateCalculator >(
            numericalCustomPartialSettings->parameterPerturbation_,
            bodyStateGetFunction,
            bodyStateSetFunction,
            environmentUpdateFunction );
    }
    else if( analyticalCustomPartialSettings!= nullptr )
    {
//        if( parameter->getParameterName( ).first != estimatable_parameters::initial_body_state )
//        {
//            throw std::runtime_error( "Error, only initial cartesian state supported for custom numerical acceleration partial" );
//        }

        customPartialCalculator = std::make_shared< estimatable_parameters::AnalyticalAccelerationPartialCalculator >(
            analyticalCustomPartialSettings->accelerationPartialFunction_,
            parameter->getParameterSize( ),
            parameter->getParameterDescription( ) );
    }
    else
    {
        throw std::runtime_error( "Error when creating custom acceleration partial, only numerical partial is supported" );
    }
    return customPartialCalculator;
}

template< typename ParameterScalarType >
std::shared_ptr< estimatable_parameters::CustomSingleAccelerationPartialCalculatorSet > createCustomAccelerationPartial(
    const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterScalarType > > parameterSet,
    const std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > accelerationModel,
    const std::pair< std::string, std::shared_ptr< simulation_setup::Body > > acceleratedBody,
    const std::pair< std::string, std::shared_ptr< simulation_setup::Body > > acceleratingBody,
    const SystemOfBodies& bodies )
{
    basic_astrodynamics::AvailableAcceleration accelerationType = getAccelerationModelType( accelerationModel );

    std::shared_ptr< estimatable_parameters::CustomSingleAccelerationPartialCalculatorSet > partialCalculatorSet =
        std::make_shared< estimatable_parameters::CustomSingleAccelerationPartialCalculatorSet >( );
    for( unsigned int i = 0; i < parameterSet->getEstimatedInitialStateParameters( ).size( ); i++ )
    {
        std::vector< std::shared_ptr< estimatable_parameters::CustomAccelerationPartialSettings > > customPartialSettings =
            parameterSet->getEstimatedInitialStateParameters( ).at( i )->getCustomPartialSettings( );
        for( unsigned int j = 0; j < customPartialSettings.size( ); j++ )
        {
            if( customPartialSettings.at( j )->accelerationMatches( acceleratedBody.first, acceleratingBody.first, accelerationType ) )
            {
                switch( parameterSet->getEstimatedInitialStateParameters( ).at( i )->getParameterName( ).first )
                {
                case estimatable_parameters::initial_body_state:
                {
                    partialCalculatorSet->customInitialStatePartials_[ parameterSet->getEstimatedInitialStateParameters( ).at( i )->getParameterName( ) ] =
                        createCustomAccelerationPartial(
                            customPartialSettings.at( j ),
                            parameterSet->getEstimatedInitialStateParameters( ).at( i ),
                            bodies );

                    break;
                }

                default:
                    throw std::runtime_error( "Error when creating custom partial calculator for " +
                                              parameterSet->getEstimatedInitialStateParameters( ).at( i )->getParameterDescription( ) +
                                              ", parameter not supported for this option" );
                }
            }
        }
    }

    for( unsigned int i = 0; i < parameterSet->getEstimatedDoubleParameters( ).size( ); i++ )
    {
        if( parameterSet->getEstimatedDoubleParameters( ).at( i )->getCustomPartialSettings( ).size( ) > 0 )
        {
            throw std::runtime_error( "Error when creating custom partial calculator for " +
                                      parameterSet->getEstimatedDoubleParameters( ).at( i )->getParameterDescription( ) +
                                      ", parameter not supported for this option" );
        }
    }

    for( unsigned int i = 0; i < parameterSet->getEstimatedVectorParameters( ).size( ); i++ )
    {

        std::vector< std::shared_ptr< estimatable_parameters::CustomAccelerationPartialSettings > > customPartialSettings =
            parameterSet->getEstimatedVectorParameters( ).at( i )->getCustomPartialSettings( );
        for( unsigned int j = 0; j < customPartialSettings.size( ); j++ )
        {
            if( customPartialSettings.at( j )->accelerationMatches( acceleratedBody.first, acceleratingBody.first, accelerationType ) )
            {
                switch( parameterSet->getEstimatedVectorParameters( ).at( i )->getParameterName( ).first )
                {
                case estimatable_parameters::custom_estimated_parameter:
                {
                    partialCalculatorSet->customVectorParameterPartials_[ parameterSet->getEstimatedVectorParameters( ).at( i )->getParameterName( ) ] =
                        createCustomAccelerationPartial(
                            customPartialSettings.at( j ),
                            parameterSet->getEstimatedVectorParameters( ).at( i ),
                            bodies );

                    break;
                }

                default:
                    throw std::runtime_error( "Error when creating custom partial calculator for " +
                                              parameterSet->getEstimatedVectorParameters( ).at( i )->getParameterDescription( ) +
                                              ", parameter not supported for this option" );
                }
            }
        }
    }
    return partialCalculatorSet;
}

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


std::shared_ptr< estimatable_parameters::NumericalAccelerationPartialSettings > getDefaultPanelledSurfaceRadiationPressurePartialSettings(
    const std::string bodyUndergoingAcceleration,
    const std::string bodyExertingAcceleration );

template< typename InitialStateParameterType = double >
std::shared_ptr< acceleration_partials::AccelerationPartial > createRadiationPressureAccelerationPartial(
    std::shared_ptr< electromagnetism::RadiationPressureAcceleration > radiationPressureAccelerationModel,
    const std::pair< std::string, std::shared_ptr< simulation_setup::Body > > acceleratedBody,
    const std::pair< std::string, std::shared_ptr< simulation_setup::Body > > acceleratingBody,
    const simulation_setup::SystemOfBodies& bodies,
    const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > >
    parametersToEstimate )
{
    using namespace electromagnetism;
    std::shared_ptr< acceleration_partials::AccelerationPartial > accelerationPartial = nullptr;

    if( std::dynamic_pointer_cast< CannonballRadiationPressureTargetModel >( radiationPressureAccelerationModel->getTargetModel( ) ) != nullptr &&
        std::dynamic_pointer_cast< IsotropicPointSourceRadiationPressureAcceleration >( radiationPressureAccelerationModel ) != nullptr )
    {
        accelerationPartial = std::make_shared< acceleration_partials::CannonBallRadiationPressurePartial >
            ( std::dynamic_pointer_cast< CannonballRadiationPressureTargetModel >( radiationPressureAccelerationModel->getTargetModel( ) ),
              std::dynamic_pointer_cast< IsotropicPointSourceRadiationPressureAcceleration >( radiationPressureAccelerationModel ),
              acceleratedBody.first, acceleratingBody.first );
    }
    else if( std::dynamic_pointer_cast< PaneledRadiationPressureTargetModel >( radiationPressureAccelerationModel->getTargetModel( ) ) != nullptr &&
             std::dynamic_pointer_cast< IsotropicPointSourceRadiationPressureAcceleration >( radiationPressureAccelerationModel ) != nullptr )
    {
        accelerationPartial = std::make_shared< acceleration_partials::PanelledRadiationPressurePartial >(
            std::dynamic_pointer_cast< IsotropicPointSourceRadiationPressureAcceleration >( radiationPressureAccelerationModel ),
            std::dynamic_pointer_cast< PaneledRadiationPressureTargetModel >( radiationPressureAccelerationModel->getTargetModel( ) ),
            acceleratedBody.first, acceleratingBody.first );
    }
    else
    {
        std::shared_ptr< estimatable_parameters::CustomSingleAccelerationPartialCalculatorSet >
            customPartialCalculator;
        if( parametersToEstimate == nullptr )
        {
            throw std::runtime_error( "Error when creating custom partials for non-isotropic source radiation pressure; no estimated parameters found" );
        }

        if( customPartialCalculator == nullptr )
        {
            for( unsigned  int i = 0; i < parametersToEstimate->getEstimatedInitialStateParameters( ).size( ); i++ )
            {
                if( !parametersToEstimate->getEstimatedInitialStateParameters( ).at( i )->hasCustomAccelerationPartialSettings(
                    acceleratedBody.first, acceleratingBody.first, basic_astrodynamics::radiation_pressure ) )
                {
                    parametersToEstimate->getEstimatedInitialStateParameters( ).at( i )->addCustomPartialSettings(
                        getDefaultPanelledSurfaceRadiationPressurePartialSettings( acceleratedBody.first, acceleratingBody.first ) );
                }

            }
            customPartialCalculator = createCustomAccelerationPartial< InitialStateParameterType >(
                parametersToEstimate, radiationPressureAccelerationModel, acceleratedBody, acceleratingBody, bodies );

        }

        // Create partial-calculating object.
        accelerationPartial = std::make_shared< acceleration_partials::RadiationPressureAccelerationPartial >
            ( acceleratedBody.first, acceleratingBody.first,
              std::dynamic_pointer_cast< PaneledSourceRadiationPressureAcceleration >( radiationPressureAccelerationModel ), customPartialCalculator );

    }
    return accelerationPartial;
}

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

    // Identify current acceleration model type
    AvailableAcceleration accelerationType = getAccelerationModelType( accelerationModel );

    std::shared_ptr< acceleration_partials::AccelerationPartial > accelerationPartial = nullptr;

    std::shared_ptr< estimatable_parameters::CustomSingleAccelerationPartialCalculatorSet >
        customPartialCalculator;
    if( parametersToEstimate != nullptr )
    {
        customPartialCalculator = createCustomAccelerationPartial(
            parametersToEstimate, accelerationModel, acceleratedBody, acceleratingBody, bodies );

        if( ( customPartialCalculator->getNumberOfCustomPartials( ) > 0 ) &&
            ( accelerationType != custom_acceleration ) &&
            ( accelerationType != radiation_pressure ) )
        {
            throw std::runtime_error( "Error, custom acceleration partials only supported for custom and radiation pressure acceleration (at present)" );
        }
    }
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
                                      "(third_body_spherical_harmonic_gravity) when making acceleration partial." );
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
    case polyhedron_gravity:
    {
        // Check if identifier is consistent with type.
        std::shared_ptr< PolyhedronGravitationalAccelerationModel > polyhedronAcceleration =
                std::dynamic_pointer_cast< PolyhedronGravitationalAccelerationModel >( accelerationModel );
        if( polyhedronAcceleration == nullptr )
        {
            throw std::runtime_error(
                        "Acceleration class type does not match acceleration type enum (polyhedron_gravity) set when making "
                        "acceleration partial." );
        }
        else
        {
            std::map< std::pair< estimatable_parameters::EstimatebleParametersEnum, std::string >,
                    std::shared_ptr< observation_partials::RotationMatrixPartial > >
                    rotationMatrixPartials = observation_partials::createRotationMatrixPartials(
                        parametersToEstimate, acceleratingBody.first, bodies );

            // Create partial-calculating object.
            accelerationPartial = std::make_shared< PolyhedronGravityPartial >
                    ( acceleratedBody.first, acceleratingBody.first, polyhedronAcceleration, rotationMatrixPartials );

        }
        break;
    }
    case third_body_polyhedron_gravity:
    {
        // Check if identifier is consistent with type.
        if( std::dynamic_pointer_cast< ThirdBodyPolyhedronGravitationalAccelerationModel >( accelerationModel ) == nullptr )
        {
            throw std::runtime_error( "Acceleration class type does not match acceleration type "
                                      "(third_body_polyhedron_gravity) when making acceleration partial." );
        }
        else
        {
            std::shared_ptr< ThirdBodyPolyhedronGravitationalAccelerationModel > thirdBodyAccelerationModel  =
                    std::dynamic_pointer_cast< ThirdBodyPolyhedronGravitationalAccelerationModel >(
                        accelerationModel );

            // Create partials for constituent central gravity accelerations
            std::shared_ptr< PolyhedronGravityPartial > accelerationPartialForBodyUndergoingAcceleration =
                    std::dynamic_pointer_cast< PolyhedronGravityPartial >(
                        createAnalyticalAccelerationPartial(
                            thirdBodyAccelerationModel->getAccelerationModelForBodyUndergoingAcceleration( ),
                            acceleratedBody, acceleratingBody, bodies, parametersToEstimate ) );
            std::shared_ptr< PolyhedronGravityPartial > accelerationPartialForCentralBody =
                    std::dynamic_pointer_cast< PolyhedronGravityPartial >(
                        createAnalyticalAccelerationPartial(
                            thirdBodyAccelerationModel->getAccelerationModelForCentralBody( ),
                            std::make_pair( thirdBodyAccelerationModel->getCentralBodyName( ),
                                            bodies.at( thirdBodyAccelerationModel->getCentralBodyName( ) ) ),
                            acceleratingBody, bodies, parametersToEstimate  ) );

            // Create partial-calculating object.
            accelerationPartial = std::make_shared< ThirdBodyGravityPartial< PolyhedronGravityPartial > >(
                        accelerationPartialForBodyUndergoingAcceleration,
                        accelerationPartialForCentralBody, acceleratedBody.first, acceleratingBody.first,
                        thirdBodyAccelerationModel->getCentralBodyName( ) );

        }
        break;
    }
    case ring_gravity:
    {
        // Check if identifier is consistent with type.
        std::shared_ptr< RingGravitationalAccelerationModel > ringAcceleration =
                std::dynamic_pointer_cast< RingGravitationalAccelerationModel >( accelerationModel );
        if( ringAcceleration == nullptr )
        {
            throw std::runtime_error(
                        "Acceleration class type does not match acceleration type enum (ring_gravity) set when making "
                        "acceleration partial." );
        }
        else
        {
            std::map< std::pair< estimatable_parameters::EstimatebleParametersEnum, std::string >,
                    std::shared_ptr< observation_partials::RotationMatrixPartial > >
                    rotationMatrixPartials = observation_partials::createRotationMatrixPartials(
                        parametersToEstimate, acceleratingBody.first, bodies );

            // Create partial-calculating object.
            accelerationPartial = std::make_shared< RingGravityPartial >
                    ( acceleratedBody.first, acceleratingBody.first, ringAcceleration, rotationMatrixPartials );

        }
        break;
    }
    case third_body_ring_gravity:
    {
        // Check if identifier is consistent with type.
        if( std::dynamic_pointer_cast< ThirdBodyRingGravitationalAccelerationModel >( accelerationModel ) == nullptr )
        {
            throw std::runtime_error( "Acceleration class type does not match acceleration type "
                                      "(third_body_ring_gravity) when making acceleration partial." );
        }
        else
        {
            std::shared_ptr< ThirdBodyRingGravitationalAccelerationModel > thirdBodyAccelerationModel  =
                    std::dynamic_pointer_cast< ThirdBodyRingGravitationalAccelerationModel >(
                        accelerationModel );

            // Create partials for constituent central gravity accelerations
            std::shared_ptr< RingGravityPartial > accelerationPartialForBodyUndergoingAcceleration =
                    std::dynamic_pointer_cast< RingGravityPartial >(
                        createAnalyticalAccelerationPartial(
                            thirdBodyAccelerationModel->getAccelerationModelForBodyUndergoingAcceleration( ),
                            acceleratedBody, acceleratingBody, bodies, parametersToEstimate ) );
            std::shared_ptr< RingGravityPartial > accelerationPartialForCentralBody =
                    std::dynamic_pointer_cast< RingGravityPartial >(
                        createAnalyticalAccelerationPartial(
                            thirdBodyAccelerationModel->getAccelerationModelForCentralBody( ),
                            std::make_pair( thirdBodyAccelerationModel->getCentralBodyName( ),
                                            bodies.at( thirdBodyAccelerationModel->getCentralBodyName( ) ) ),
                            acceleratingBody, bodies, parametersToEstimate  ) );

            // Create partial-calculating object.
            accelerationPartial = std::make_shared< ThirdBodyGravityPartial< RingGravityPartial > >(
                        accelerationPartialForBodyUndergoingAcceleration,
                        accelerationPartialForCentralBody, acceleratedBody.first, acceleratingBody.first,
                        thirdBodyAccelerationModel->getCentralBodyName( ) );

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
            throw std::runtime_error(
                        "Acceleration class type does not match acceleration type enum (empirical) set when making acceleration partial." );

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
           throw std::runtime_error(
                        "Acceleration class type does not match acceleration type enum (mom. wheel desat.) set when making acceleration partial." );

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
    {
        // Check if identifier is consistent with type.
        std::shared_ptr< CustomAccelerationModel > customAccelerationModel =
            std::dynamic_pointer_cast< CustomAccelerationModel >( accelerationModel );
        if( customAccelerationModel == nullptr )
        {
            throw std::runtime_error(
                "Acceleration class type does not match acceleration type enum (custom_acceleration) set when making "
                "acceleration partial." );
        }
        else
        {
            // Create partial-calculating object.
            accelerationPartial = std::make_shared< CustomAccelerationPartial >
                ( acceleratedBody.first, acceleratingBody.first, customAccelerationModel, customPartialCalculator );

        }
        break;
    }
    case radiation_pressure:
    {
        // Check if identifier is consistent with type.
        std::shared_ptr< RadiationPressureAcceleration > radiationPressureAccelerationModel =
            std::dynamic_pointer_cast< RadiationPressureAcceleration >( accelerationModel );
        if( radiationPressureAccelerationModel == nullptr )
        {
            throw std::runtime_error(
                "Acceleration class type does not match acceleration type enum (radiation_pressure) set when making "
                "acceleration partial." );
        }
        else
        {
            // Create partial-calculating object.
            accelerationPartial = createRadiationPressureAccelerationPartial(
                radiationPressureAccelerationModel, acceleratedBody, acceleratingBody, bodies, parametersToEstimate );

        }
        break;
    }
    case thrust_acceleration:
    {
        // Check if identifier is consistent with type.
        std::shared_ptr< propulsion::ThrustAcceleration > thrustAcceleration =
                std::dynamic_pointer_cast< propulsion::ThrustAcceleration >( accelerationModel );
        if( thrustAcceleration == nullptr )
        {
            throw std::runtime_error(
                        "Acceleration class type does not match acceleration type enum (thrust) set when making acceleration partial." );

        }
        else
        {
            std::map< std::pair< estimatable_parameters::EstimatebleParametersEnum, std::string >,
                    std::shared_ptr< observation_partials::RotationMatrixPartial > >
                    rotationMatrixPartials;
            if( parametersToEstimate != nullptr )
            {
                    rotationMatrixPartials = observation_partials::createRotationMatrixPartials(
                        parametersToEstimate, acceleratingBody.first, bodies );
            }


            // Create partial-calculating object.
            accelerationPartial = std::make_shared< ThrustAccelerationPartial >(
                        thrustAcceleration, acceleratedBody.first, rotationMatrixPartials );
        }
        break;
    }
    case yarkovsky_acceleration:
    {
        // Check if identifier is consistent with type.
        std::shared_ptr< YarkovskyAcceleration > yarkovskyAcceleration =
            std::dynamic_pointer_cast< YarkovskyAcceleration >( accelerationModel );
        if( yarkovskyAcceleration == nullptr )
        {
            throw std::runtime_error(
                "Acceleration class type does not match acceleration type enum (yarkovsky_acceleration) set when making "
                "acceleration partial." );
        }
        else
        {
            // Create partial-calculating object.
            accelerationPartial = std::make_shared< YarkovskyAccelerationPartial >
                ( yarkovskyAcceleration, acceleratedBody.first, acceleratingBody.first );

        }
        break;
    }
    default:
        std::string errorMessage = "Acceleration model " + std::to_string( accelerationType ) +
                " not found when making acceleration partial";
        throw std::runtime_error( errorMessage );
        break;
    }

    return accelerationPartial;
}

std::map< std::string, std::shared_ptr< acceleration_partials::AccelerationPartial > > createEihAccelerationPartials(
    const std::map< std::string, std::shared_ptr< relativity::EinsteinInfeldHoffmannAcceleration > > eihAccelerations );

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
    std::map< std::string, std::shared_ptr< relativity::EinsteinInfeldHoffmannAcceleration > > eihAccelerations;
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
                            if( basic_astrodynamics::getAccelerationModelType( innerAccelerationIterator->second[ j ] ) !=
                                basic_astrodynamics::einstein_infeld_hoffmann_acceleration )
                            {
                                // Create single partial object
                                std::shared_ptr<acceleration_partials::AccelerationPartial> currentAccelerationPartial =
                                    createAnalyticalAccelerationPartial(
                                        innerAccelerationIterator->second[ j ],
                                        std::make_pair( acceleratedBody, acceleratedBodyObject ),
                                        std::make_pair( acceleratingBody, acceleratingBodyObject ),
                                        bodies, parametersToEstimate );

                                if ( currentAccelerationPartial != nullptr )
                                {
                                    accelerationPartialVector.push_back( currentAccelerationPartial );
                                    accelerationPartialsMap[ acceleratedBody ][ acceleratingBody ].push_back(
                                        currentAccelerationPartial );
                                }
                            }
                            else
                            {
                                std::shared_ptr< relativity::EinsteinInfeldHoffmannAcceleration > eihAcceleration = std::dynamic_pointer_cast<
                                    relativity::EinsteinInfeldHoffmannAcceleration >( innerAccelerationIterator->second[ j ] );
                                if( eihAcceleration == nullptr )
                                {
                                    throw std::runtime_error( "Error when making acceleration partials, eih acceleration has incomaptible type" );
                                }
                                else if( eihAccelerations.count( acceleratedBody ) != 0 )
                                {
                                    throw std::runtime_error( "Error when making acceleration partials, multiple eih accelerations found for single body " + acceleratedBody );
                                }
                                else
                                {
                                    eihAccelerations[ acceleratedBody ] = eihAcceleration;
                                }

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

    if( eihAccelerations.size( ) > 0 )
    {
        std::map< std::string, std::shared_ptr< acceleration_partials::AccelerationPartial > >  eihPartials = createEihAccelerationPartials(
            eihAccelerations );
        unsigned int additionCounter =  0;
        for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
        {
            if ( eihPartials.count( initialDynamicalParameters.at( i )->getParameterName( ).second.first ) > 0 )
            {
                accelerationPartialsList[ i ].push_back( eihPartials.at( initialDynamicalParameters.at( i )->getParameterName( ).second.first ) );
                additionCounter++;
            }
        }

        if( !( additionCounter == eihPartials.size() ) )
        {
            throw std::runtime_error( "Error when making acceleration partials, not all eih partials have been used" );
        }
    }
    return accelerationPartialsList;
}


} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_CREATEACCELERATIONPARTIALS_H
