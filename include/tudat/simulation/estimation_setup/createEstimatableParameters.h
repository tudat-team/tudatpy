/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEESTIMATABLEPARAMETERS_H
#define TUDAT_CREATEESTIMATABLEPARAMETERS_H

#include "tudat/astro/basic_astro/accelerationModel.h"

#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/initialTranslationalState.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/initialRotationalState.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/initialMassState.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/constantDragCoefficient.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/constantRotationRate.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/constantRotationalOrientation.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/empiricalAccelerationCoefficients.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/gravitationalParameter.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/observationBiasParameter.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/groundStationPosition.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/sphericalHarmonicCosineCoefficients.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/sphericalHarmonicSineCoefficients.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/radiationPressureCoefficient.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/ppnParameters.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/equivalencePrincipleViolationParameter.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/tidalLoveNumber.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/directTidalTimeLag.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/inverseTidalQualityFactor.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/meanMomentOfInertiaParameter.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/desaturationDeltaV.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/periodicSpinVariation.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/polarMotionAmplitude.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/coreFactor.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/freeCoreNutationRate.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/desaturationDeltaV.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/longitudeLibrationAmplitude.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/constantThrust.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/yarkovskyParameter.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/referencePointPosition.h"
#include "tudat/astro/relativity/metric.h"
#include "tudat/astro/basic_astro/accelerationModelTypes.h"
#include "tudat/simulation/estimation_setup/estimatableParameterSettings.h"
#include "tudat/simulation/propagation_setup/dynamicsSimulator.h"
#include "tudat/simulation/environment_setup/body.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to get a list of acceleration models that is to be linked to the given parameter
/*!
 *  Function to get a list of acceleration models that is to be linked to the given parameter, from single-arc propagator settings.
 *  For selected parameter types, this function finds the acceleration models to which they have to be linked to fully create
 *  the parameter objects. If  parameter type needs no acceleration, or no compatibel acceleration is found, an empty list is
 *  returned.
 *  \param propagatorSettings Single-arc propagator settings, from which acceleration models are to be extracted
 *  \param parameterSettings Settings for parameter settings for which acceleration models are to be found
 *  \return List of acceleration models that is to be linked to parameter defined by parameterSettings
 */
template< typename StateScalarType, typename TimeType >
std::vector< std::shared_ptr< basic_astrodynamics::AccelerationModel3d > > getAccelerationModelsListForParameters(
        const std::shared_ptr< propagators::SingleArcPropagatorSettings< StateScalarType, TimeType > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > parameterSettings )
{
    using namespace estimatable_parameters;
    std::vector< std::shared_ptr< basic_astrodynamics::AccelerationModel3d > > accelerationModelList;

    // Retrieve acceleration models
    basic_astrodynamics::AccelerationMap accelerationModelMap = getAccelerationMapFromPropagatorSettings(
                propagatorSettings );

    // Check parameter type
    switch( parameterSettings->parameterType_.first )
    {
    //  Empirical acceleration coefficeints need to be linked to empirical acceleration object
    case empirical_acceleration_coefficients:
    {
        std::shared_ptr< EmpiricalAccelerationEstimatableParameterSettings > empiricalAccelerationSettings =
                std::dynamic_pointer_cast< EmpiricalAccelerationEstimatableParameterSettings >( parameterSettings );

        // Check if acceleration model with required bodies undergoing/exerting accelerations exist
        if( accelerationModelMap.count( empiricalAccelerationSettings->parameterType_.second.first ) != 0 )
        {
            if( accelerationModelMap.at( empiricalAccelerationSettings->parameterType_.second.first ).count(
                        empiricalAccelerationSettings->parameterType_.second.second ) != 0 )

            {
                // Retrieve acceleration model.
                std::vector< std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > >
                        accelerationModelListToCheck =
                        accelerationModelMap.at( empiricalAccelerationSettings->parameterType_.second.first ).at(
                            empiricalAccelerationSettings->parameterType_.second.second );
                for( unsigned int i = 0; i < accelerationModelListToCheck.size( ); i++ )
                {
                    if( basic_astrodynamics::getAccelerationModelType( accelerationModelListToCheck[ i ] ) ==
                            basic_astrodynamics::empirical_acceleration )
                    {
                        accelerationModelList.push_back( accelerationModelListToCheck[ i ] );
                    }
                }
            }
        }
        break;
    }
        // Arc-wise empirical acceleration coefficeints need to be linked to empirical acceleration object
    case arc_wise_empirical_acceleration_coefficients:
    {
        std::shared_ptr< ArcWiseEmpiricalAccelerationEstimatableParameterSettings > empiricalAccelerationSettings =
                std::dynamic_pointer_cast< ArcWiseEmpiricalAccelerationEstimatableParameterSettings >( parameterSettings );

        // Check if acceleration model with required bodies undergoing/exerting accelerations exist
        if( accelerationModelMap.count( empiricalAccelerationSettings->parameterType_.second.first ) != 0 )
        {
            if( accelerationModelMap.at( empiricalAccelerationSettings->parameterType_.second.first ).count(
                        empiricalAccelerationSettings->parameterType_.second.second ) != 0 )

            {
                // Retrieve acceleration model.
                std::vector< std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > >
                        accelerationModelListToCheck =
                        accelerationModelMap.at( empiricalAccelerationSettings->parameterType_.second.first ).at(
                            empiricalAccelerationSettings->parameterType_.second.second );
                for( unsigned int i = 0; i < accelerationModelListToCheck.size( ); i++ )
                {
                    if( basic_astrodynamics::getAccelerationModelType( accelerationModelListToCheck[ i ] ) ==
                            basic_astrodynamics::empirical_acceleration )
                    {
                        accelerationModelList.push_back( accelerationModelListToCheck[ i ] );
                    }
                }
            }
        }
        break;
    }
    // Direct tidal time lags need to be linked to direct tidal acceleration
    case direct_dissipation_tidal_time_lag:
    {
        std::shared_ptr< DirectTidalTimeLagEstimatableParameterSettings > dissipationTimeLagSettings =
                std::dynamic_pointer_cast< DirectTidalTimeLagEstimatableParameterSettings >( parameterSettings );
        std::string currentBodyName =  parameterSettings ->parameterType_.second.first;
        if( dissipationTimeLagSettings == nullptr )
        {
            throw std::runtime_error( "Error, expected dissipation time lag parameter settings." );
        }
        else
        {
            std::vector< std::shared_ptr< gravitation::DirectTidalDissipationAcceleration > > tidalAccelerationModelList =
                    gravitation::getTidalDissipationAccelerationModels(
                        accelerationModelMap, currentBodyName, dissipationTimeLagSettings->deformingBodies_ );
            for( unsigned int i = 0; i < tidalAccelerationModelList.size( ); i++ )
            {
                accelerationModelList.push_back( tidalAccelerationModelList.at( i ) );
            }

        }
        break;
    }
        // Desaturation Delta V needs to be linked to destauration acceleration
    case desaturation_delta_v_values:
    {
        // Check if acceleration model with required bodies undergoing/exerting accelerations exist
        if( accelerationModelMap.count( parameterSettings->parameterType_.second.first ) != 0 )
        {
            if( accelerationModelMap.at( parameterSettings->parameterType_.second.first ).count(
                        parameterSettings->parameterType_.second.first ) != 0 )

            {
                // Retrieve acceleration model.
                std::vector< std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > >
                        accelerationModelListToCheck =
                        accelerationModelMap.at( parameterSettings->parameterType_.second.first ).at(
                            parameterSettings->parameterType_.second.first );
                for( unsigned int i = 0; i < accelerationModelListToCheck.size( ); i++ )
                {
                    if( basic_astrodynamics::getAccelerationModelType( accelerationModelListToCheck[ i ] ) ==
                            basic_astrodynamics::momentum_wheel_desaturation_acceleration )
                    {
                        accelerationModelList.push_back( accelerationModelListToCheck[ i ] );
                    }
                }
            }
        }
        break;
    }
    // Inverse tidal quality factor to be linked to direct tidal acceleration
    case inverse_tidal_quality_factor:
    {
        std::shared_ptr< InverseTidalQualityFactorEstimatableParameterSettings > qualityFactorSettings =
                std::dynamic_pointer_cast< InverseTidalQualityFactorEstimatableParameterSettings >( parameterSettings );
        std::string currentBodyName =  parameterSettings ->parameterType_.second.first;
        if( qualityFactorSettings == nullptr )
        {
            throw std::runtime_error( "Error, expected inverse tidal quality factor parameter settings." );
        }
        else
        {
            std::vector< std::shared_ptr< gravitation::DirectTidalDissipationAcceleration > > tidalAccelerationModelList =
                    gravitation::getTidalDissipationAccelerationModels( accelerationModelMap, currentBodyName, qualityFactorSettings->deformingBodies_ );
            for( unsigned int i = 0; i < tidalAccelerationModelList.size( ); i++ )
            {
                accelerationModelList.push_back( tidalAccelerationModelList.at( i ) );
            }

        }
        break;
    }
    case yarkovsky_parameter:
    {
        if( parameterSettings == nullptr )
        {
            throw std::runtime_error( "Error, expected Yarkovsky parameter settings." );
        }
        else
        {
            if( accelerationModelMap.at( parameterSettings->parameterType_.second.first ).count(
                parameterSettings->parameterType_.second.second ) != 0 )

            {
                // Retrieve acceleration model.
                std::vector< std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > >
                    accelerationModelListToCheck =
                    accelerationModelMap.at( parameterSettings->parameterType_.second.first ).at(
                        parameterSettings->parameterType_.second.second );
                for( unsigned int i = 0; i < accelerationModelListToCheck.size( ); i++ )
                {
                    if( basic_astrodynamics::getAccelerationModelType( accelerationModelListToCheck[ i ] ) ==
                        basic_astrodynamics::yarkovsky_acceleration )
                    {
                        accelerationModelList.push_back( accelerationModelListToCheck[ i ] );
                    }
                }
            }

        }
        break;
    }
    default:
        break;
    }
    return accelerationModelList;
}

//! Function to get a list of acceleration models that is to be linked to the given parameter
/*!
 *  Function to get a list of acceleration models that is to be linked to the given parameter, from multi-arc propagator settings.
 *  For selected parameter types, this function finds the acceleration models to which they have to be linked to fully create
 *  the parameter objects. If  parameter type needs no acceleration, or no compatible acceleration is found, an empty list is
 *  returned.
 *  \param propagatorSettings Single-arc propagator settings, from which acceleration models are to be extracted
 *  \param parameterSettings Settings for parameter settings for which acceleration models are to be found
 *  \return List of acceleration models (from all arcs) that is to be linked to parameter defined by parameterSettings
 */
template< typename StateScalarType, typename TimeType >
std::vector< std::shared_ptr< basic_astrodynamics::AccelerationModel3d > > getAccelerationModelsListForParameters(
        const std::shared_ptr< propagators::MultiArcPropagatorSettings< StateScalarType, TimeType > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > parameterSettings )
{
    std::vector< std::shared_ptr< basic_astrodynamics::AccelerationModel3d > > accelerationModelList;
    for( unsigned int i = 0; i < propagatorSettings->getSingleArcSettings( ).size( ); i++ )
    {
        std::vector< std::shared_ptr< basic_astrodynamics::AccelerationModel3d > > singleArcAccelerationModelList =
                getAccelerationModelsListForParameters(
                    propagatorSettings->getSingleArcSettings( ).at( i ), parameterSettings );
        accelerationModelList.insert(
                    accelerationModelList.end( ), singleArcAccelerationModelList.begin( ), singleArcAccelerationModelList.end( ) );
    }
    return accelerationModelList;
}

//! Function to get a list of acceleration models that is to be linked to the given parameter
/*!
 *  Function to get a list of acceleration models that is to be linked to the given parameter, from hybrid-arc propagator settings.
 *  For selected parameter types, this function finds the acceleration models to which they have to be linked to fully create
 *  the parameter objects. If  parameter type needs no acceleration, or no compatible acceleration is found, an empty list is
 *  returned.
 *  \param propagatorSettings Single-arc propagator settings, from which acceleration models are to be extracted
 *  \param parameterSettings Settings for parameter settings for which acceleration models are to be found
 *  \return List of acceleration models (from all arcs) that is to be linked to parameter defined by parameterSettings
 */
template< typename StateScalarType, typename TimeType >
std::vector< std::shared_ptr< basic_astrodynamics::AccelerationModel3d > > getAccelerationModelsListForParameters(
        const std::shared_ptr< propagators::HybridArcPropagatorSettings< StateScalarType, TimeType > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > parameterSettings )
{
    std::vector< std::shared_ptr< basic_astrodynamics::AccelerationModel3d > > multiArcAccelerationModelList;
    for( unsigned int i = 0; i < propagatorSettings->getMultiArcPropagatorSettings( )->getSingleArcSettings( ).size( ); i++ )
    {
        std::vector< std::shared_ptr< basic_astrodynamics::AccelerationModel3d > > singleArcAccelerationModelList =
                getAccelerationModelsListForParameters(
                    propagatorSettings->getMultiArcPropagatorSettings( )->getSingleArcSettings( ).at( i ), parameterSettings );
        multiArcAccelerationModelList.insert(
                    multiArcAccelerationModelList.end( ), singleArcAccelerationModelList.begin( ), singleArcAccelerationModelList.end( ) );
    }

    std::vector< std::shared_ptr< basic_astrodynamics::AccelerationModel3d > > singleArcAccelerationModelList =
            getAccelerationModelsListForParameters(
                propagatorSettings->getSingleArcPropagatorSettings( ), parameterSettings );

    if( singleArcAccelerationModelList.size( ) != 0 && multiArcAccelerationModelList.size( ) != 0 )
    {
        std::cerr<<"Warning when linking parameter to acceleration model in hybrid arc propagation. Dependencies found in both single- and multi-arc segments."<<std::endl;
    }

    std::vector< std::shared_ptr< basic_astrodynamics::AccelerationModel3d > > accelerationModelList =
            multiArcAccelerationModelList;
    accelerationModelList.insert(
                accelerationModelList.end( ), singleArcAccelerationModelList.begin( ), singleArcAccelerationModelList.end( ) );

    return accelerationModelList;
}

//! Function to get a list of acceleration models that is to be linked to the given parameter
/*!
 *  Function to get a list of acceleration models that is to be linked to the given parameter, from any propagator settings.
 *  For selected parameter types, this function finds the acceleration models to which they have to be linked to fully create
 *  the parameter objects. If  parameter type needs no acceleration, or no compatible acceleration is found, an empty list is
 *  returned.
 *  \param propagatorSettings Single-arc propagator settings, from which acceleration models are to be extracted
 *  \param parameterSettings Settings for parameter settings for which acceleration models are to be found
 *  \return List of acceleration models (from all arcs if applicable) that is to be linked to parameter defined by
 *  parameterSettings
 */
template< typename StateScalarType, typename TimeType >
std::vector< std::shared_ptr< basic_astrodynamics::AccelerationModel3d > > getAccelerationModelsListForParametersFromBase(
        const std::shared_ptr< propagators::PropagatorSettings< StateScalarType > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > parameterSettings )
{
    std::vector< std::shared_ptr< basic_astrodynamics::AccelerationModel3d > > accelerationModelList;

    if( std::dynamic_pointer_cast< propagators::SingleArcPropagatorSettings< StateScalarType, TimeType > >( propagatorSettings ) != nullptr )
    {
        accelerationModelList = getAccelerationModelsListForParameters(
                    std::dynamic_pointer_cast< propagators::SingleArcPropagatorSettings< StateScalarType, TimeType > >( propagatorSettings ),
                    parameterSettings );
    }
    else if( std::dynamic_pointer_cast< propagators::MultiArcPropagatorSettings< StateScalarType, TimeType > >( propagatorSettings ) != nullptr )
    {
        accelerationModelList = getAccelerationModelsListForParameters(
                    std::dynamic_pointer_cast< propagators::MultiArcPropagatorSettings< StateScalarType, TimeType > >( propagatorSettings ),
                    parameterSettings );
    }
    else if( std::dynamic_pointer_cast< propagators::HybridArcPropagatorSettings< StateScalarType, TimeType > >( propagatorSettings ) != nullptr )
    {
        accelerationModelList = getAccelerationModelsListForParameters(
                    std::dynamic_pointer_cast< propagators::HybridArcPropagatorSettings< StateScalarType, TimeType > >( propagatorSettings ),
                    parameterSettings );
    }

    if( accelerationModelList.size( ) == 0 )
    {
        throw std::runtime_error( "Error when getting acceleration model for parameter " +
                                  std::to_string( parameterSettings->parameterType_.first ) + ", no acceleration model found." );
    }

    return accelerationModelList;
}

template< typename InitialStateParameterType = double, typename TimeType = double >
std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > > getInitialStateParameterSettings(
        const std::shared_ptr< propagators::PropagatorSettings< InitialStateParameterType > > propagatorSettings,
        const SystemOfBodies& bodies,
        const std::vector< double > arcStartTimes = std::vector< double >( ) );

template< typename InitialStateParameterType = double, typename TimeType = double >
std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > > getInitialMultiArcParameterSettings(
        const std::shared_ptr< propagators::MultiArcPropagatorSettings< InitialStateParameterType, TimeType > > propagatorSettings,
        const SystemOfBodies& bodies,
        const std::vector< double > arcStartTimes )
{
    using namespace estimatable_parameters;
    using namespace propagators;

    std::vector< std::shared_ptr< SingleArcPropagatorSettings< InitialStateParameterType, TimeType > > > singleArcSettings =
            propagatorSettings->getSingleArcSettings( );
    std::vector< std::shared_ptr< TranslationalStatePropagatorSettings< InitialStateParameterType, TimeType > > >
            singleArcTranslationalSettings;

    std::vector< std::string > propagatedBodies;
    std::vector< std::vector< std::string > > centralBodiesPerArc;
    std::vector< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > initialStates;

    for( unsigned int i = 0; i < singleArcSettings.size( ); i++ )
    {
        singleArcTranslationalSettings.push_back(
                    std::dynamic_pointer_cast< TranslationalStatePropagatorSettings< InitialStateParameterType, TimeType > >(
                        singleArcSettings.at( i ) ) );
        if( singleArcTranslationalSettings.at( i ) == nullptr )
        {
            throw std::runtime_error( "Only translational state supported when auto-creating multi-arc initial state settings" );
        }
        else
        {

            initialStates.push_back( singleArcTranslationalSettings.at( i )->getInitialStates( ) );
            centralBodiesPerArc.push_back( singleArcTranslationalSettings.at( i )->centralBodies_ );
            if( i == 0 )
            {
                propagatedBodies = singleArcTranslationalSettings.at( i )->bodiesToIntegrate_;
            }
            else
            {
                if( !( propagatedBodies == ( singleArcTranslationalSettings.at( i )->bodiesToIntegrate_ ) ) )
                {
                    throw std::runtime_error( "Only equal bodies per arc supported when auto-creating multi-arc initial state settings" );
                }
            }
        }
    }

    std::vector< std::vector< std::string > > centralBodiesPerBody;
    centralBodiesPerBody.resize( centralBodiesPerArc.at( 0 ).size( ) );

    for( unsigned int i = 0; i < centralBodiesPerArc.size( ); i++ )
    {
        for( unsigned int j = 0; j < centralBodiesPerArc.at( i ).size( ); j++ )
        {
            if( i == 0 )
            {
                centralBodiesPerBody.at( j ).resize( centralBodiesPerArc.size( ) );
            }
            centralBodiesPerBody.at( j ).at( i ) = centralBodiesPerArc.at( i ).at( j );
        }
    }

    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > > arcwiseInitialStates;
    for( unsigned int i = 0; i < propagatedBodies.size( ); i++ )
    {
        Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > multiArcInitialStateValue =
                Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 >::Zero( 6 * initialStates.size( ) );
        for( unsigned int j = 0; j < initialStates.size( ); j++ )
        {
            multiArcInitialStateValue.segment( j * 6, 6 ) = initialStates.at( j ).segment( i * 6, 6 );
        }
        arcwiseInitialStates.push_back(
                    std::make_shared<
                    ArcWiseInitialTranslationalStateEstimatableParameterSettings< InitialStateParameterType > >(
                        propagatedBodies.at( i ),
                        multiArcInitialStateValue,
                        arcStartTimes,
                        centralBodiesPerBody.at( i ),
                        bodies.getFrameOrientation( ) ) );
    }

    return arcwiseInitialStates;
}

template< typename InitialStateParameterType = double, typename TimeType = double >
std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > > getInitialHybridArcParameterSettings(
        const std::shared_ptr< propagators::HybridArcPropagatorSettings< InitialStateParameterType, TimeType > > propagatorSettings,
        const SystemOfBodies& bodies,
        const std::vector< double > arcStartTimes )
{
    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > > multiArcParameters =
            getInitialMultiArcParameterSettings< InitialStateParameterType, TimeType >(
                propagatorSettings->getMultiArcPropagatorSettings( ), bodies, arcStartTimes );
    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > > singleArcParameters =
            getInitialStateParameterSettings< InitialStateParameterType, TimeType >(
                propagatorSettings->getSingleArcPropagatorSettings( ), bodies );
    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > > hybirdArcParameters = multiArcParameters;

    hybirdArcParameters.insert( hybirdArcParameters.end( ), singleArcParameters.begin( ), singleArcParameters.end( ) );
    return hybirdArcParameters;
}


template< typename InitialStateParameterType, typename TimeType >
std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > > getInitialStateParameterSettings(
        const std::shared_ptr< propagators::PropagatorSettings< InitialStateParameterType > > propagatorSettings,
        const SystemOfBodies& bodies,
        const std::vector< double > arcStartTimes )
{
    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > > initialStateParameterSettings;

    using namespace propagators;

    // Process single-arc settings
    if( std::dynamic_pointer_cast< SingleArcPropagatorSettings< InitialStateParameterType, TimeType > >( propagatorSettings ) != nullptr )
    {
        std::shared_ptr< SingleArcPropagatorSettings< InitialStateParameterType, TimeType > > singleArcSettings =
                std::dynamic_pointer_cast< SingleArcPropagatorSettings< InitialStateParameterType, TimeType > >( propagatorSettings );
        switch( singleArcSettings->getStateType( ) )
        {
        case hybrid:
        {
            std::shared_ptr< MultiTypePropagatorSettings< InitialStateParameterType, TimeType > > multiTypePropagatorSettings =
                    std::dynamic_pointer_cast< MultiTypePropagatorSettings< InitialStateParameterType, TimeType > >( propagatorSettings );


            std::map< IntegratedStateType, std::vector< std::shared_ptr< SingleArcPropagatorSettings< InitialStateParameterType, TimeType > > > >
                    propagatorSettingsMap = multiTypePropagatorSettings->propagatorSettingsMap_;
            for( auto propIterator : propagatorSettingsMap )
            {
                for( unsigned int i = 0; i < propIterator.second.size( ); i++ )
                {
                    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > >
                            singleTypeinitialStateParameterSettings =  getInitialStateParameterSettings< InitialStateParameterType, TimeType >(
                                propIterator.second.at( i ), bodies );
                    initialStateParameterSettings.insert(
                                initialStateParameterSettings.end( ),
                                singleTypeinitialStateParameterSettings.begin( ),
                                singleTypeinitialStateParameterSettings.end( ) );
                }
            }
            break;
        }
        case translational_state:
        {
            std::shared_ptr< TranslationalStatePropagatorSettings< InitialStateParameterType, TimeType > > translationalPropagatorSettings =
                    std::dynamic_pointer_cast< TranslationalStatePropagatorSettings< InitialStateParameterType, TimeType > >( propagatorSettings );

            // Retrieve estimated and propagated translational states, and check equality.
            std::vector< std::string > propagatedBodies = translationalPropagatorSettings->bodiesToIntegrate_;
            std::vector< std::string > centralBodies = translationalPropagatorSettings->centralBodies_;

            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > initialStates =  translationalPropagatorSettings->getInitialStates( );
            for( unsigned int i = 0; i < propagatedBodies.size( ); i++ )
            {
                initialStateParameterSettings.push_back(
                            std::make_shared< estimatable_parameters::InitialTranslationalStateEstimatableParameterSettings<
                            InitialStateParameterType > >(
                                propagatedBodies.at( i ), initialStates.segment( i * 6, 6 ), centralBodies.at( i ),
                                bodies.getFrameOrientation( ) ) );
            }
            break;

        }
        case rotational_state:
        {
            std::shared_ptr< RotationalStatePropagatorSettings< InitialStateParameterType, TimeType > > rotationalPropagatorSettings =
                    std::dynamic_pointer_cast< RotationalStatePropagatorSettings< InitialStateParameterType, TimeType > >( propagatorSettings );

            // Retrieve estimated and propagated translational states, and check equality.
            std::vector< std::string > propagatedBodies = rotationalPropagatorSettings->bodiesToIntegrate_;

            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > initialStates =  rotationalPropagatorSettings->getInitialStates( );
            for( unsigned int i = 0; i < propagatedBodies.size( ); i++ )
            {
                initialStateParameterSettings.push_back(
                            std::make_shared< estimatable_parameters::InitialRotationalStateEstimatableParameterSettings<
                            InitialStateParameterType > >(
                                propagatedBodies.at( i ), initialStates.segment( i * 7, 7 ).template cast< InitialStateParameterType >( ), bodies.getFrameOrientation( ) ) );
            }
            break;
        }
        case body_mass_state:
        {
            std::shared_ptr< MassPropagatorSettings< InitialStateParameterType, TimeType > > massPropagatorSettings =
                    std::dynamic_pointer_cast< MassPropagatorSettings< InitialStateParameterType, TimeType > >( propagatorSettings );

            std::vector< std::string > propagatedBodies = massPropagatorSettings->bodiesWithMassToPropagate_;
            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > initialStates =  massPropagatorSettings->getInitialStates( );
            for( unsigned int i = 0; i < propagatedBodies.size( ); i++ )
            {
                initialStateParameterSettings.push_back(
                            std::make_shared< estimatable_parameters::InitialMassEstimatableParameterSettings<
                            InitialStateParameterType > >(
                                propagatedBodies.at( i ), initialStates( i ) ) );
            }
            break;
        }
        case custom_state:
        {
            throw std::runtime_error( "Error, cannot estimate initial custom state" );
        }
        default:
            throw std::runtime_error( "Error, did not recognize single-arc state type when identifying propagator settings for estimatable parameter settings." );
        }
    }
    else if( std::dynamic_pointer_cast< MultiArcPropagatorSettings< InitialStateParameterType, TimeType > >( propagatorSettings ) != nullptr )
    {
        std::shared_ptr< MultiArcPropagatorSettings< InitialStateParameterType, TimeType > > multiArcSettings =
                std::dynamic_pointer_cast< MultiArcPropagatorSettings< InitialStateParameterType, TimeType > >( propagatorSettings ) ;
        if( arcStartTimes.size( ) == 0 )
        {
            throw std::runtime_error( "Error when parsing propagator settings for estimatable parameter settings; multi-arc settings found, but no arc times" );
        }
        initialStateParameterSettings = getInitialMultiArcParameterSettings(
                    multiArcSettings, bodies, arcStartTimes );
    }
    else if( std::dynamic_pointer_cast< HybridArcPropagatorSettings< InitialStateParameterType, TimeType > >( propagatorSettings ) != nullptr )
    {
        std::shared_ptr< HybridArcPropagatorSettings< InitialStateParameterType, TimeType > > hybridArcSettings =
                std::dynamic_pointer_cast< HybridArcPropagatorSettings< InitialStateParameterType, TimeType > >( propagatorSettings );
        if( arcStartTimes.size( ) == 0 )
        {
            throw std::runtime_error( "Error when parsing propagator settings for estimatable parameter settings; hybric-arc settings found, but no arc times" );
        }
        initialStateParameterSettings = getInitialHybridArcParameterSettings(
                    hybridArcSettings, bodies, arcStartTimes );
    }

    return initialStateParameterSettings;
}


//! Function to create interface object for estimating parameters representing an initial dynamical state.
/*!
 *  Function to create interface object for estimating parameters representing an initial dynamical state.
 *  \param bodies Map of body objects containing the fll simulation environment.
 *  \param parameterSettings Object defining the parameter interface that is to be created.
 *  \return Interface object for estimating an initial state.
 */
template< typename InitialStateParameterType = double >
std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::Matrix
< InitialStateParameterType, Eigen::Dynamic, 1 > > > createInitialDynamicalStateParameterToEstimate(
        const SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSettings >& parameterSettings )
{
    using namespace tudat::estimatable_parameters;

    std::shared_ptr< EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > >
            initialStateParameterToEstimate;

    // Check consistency of input.
    if( !isParameterDynamicalPropertyInitialState( parameterSettings->parameterType_.first ) )
    {
        std::string errorMessage = "Error when requesting to make initial state parameter " +
                std::to_string( parameterSettings->parameterType_.first ) + " of " +
                parameterSettings->parameterType_.second.first +
                ", parameter is not an initial state parameter ";
        throw std::runtime_error( errorMessage );
    }
    else
    {
        // Identify state that is to be estimation
        switch( parameterSettings->parameterType_.first )
        {
        case initial_body_state:

            // Check consistency of input.
            if( std::dynamic_pointer_cast<
                    InitialTranslationalStateEstimatableParameterSettings< InitialStateParameterType > >(
                        parameterSettings ) == nullptr )
            {
                throw std::runtime_error( "Error when making body initial state parameter, settings type is incompatible" );
            }
            else
            {
                std::shared_ptr< InitialTranslationalStateEstimatableParameterSettings< InitialStateParameterType > >
                        initialStateSettings =
                        std::dynamic_pointer_cast<
                        InitialTranslationalStateEstimatableParameterSettings< InitialStateParameterType > >(
                            parameterSettings );

                Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > initialTranslationalState;

                // If initial time is not defined, use preset initial state
                if( ! ( initialStateSettings->initialTime_ == initialStateSettings->initialTime_  ) )
                {
                    initialTranslationalState = initialStateSettings->initialStateValue_;


                }
                // Compute initial state from environment
                else
                {
                    initialTranslationalState = propagators::getInitialStateOfBody
                            < double, InitialStateParameterType >(
                                initialStateSettings->parameterType_.second.first, initialStateSettings->centralBody_,
                                bodies, initialStateSettings->initialTime_ );

                }

                // Create translational state estimation interface object
                initialStateParameterToEstimate =
                        std::make_shared< InitialTranslationalStateParameter< InitialStateParameterType > >(
                            initialStateSettings->parameterType_.second.first, initialTranslationalState,
                            initialStateSettings->centralBody_,
                            initialStateSettings->frameOrientation_ );
            }
            break;
        case arc_wise_initial_body_state:
            if( std::dynamic_pointer_cast< ArcWiseInitialTranslationalStateEstimatableParameterSettings<
                    InitialStateParameterType > >( parameterSettings ) == nullptr )
            {
                throw std::runtime_error(
                            "Error when making body initial state parameter, settings type is incompatible" );
            }
            else
            {
                std::shared_ptr< ArcWiseInitialTranslationalStateEstimatableParameterSettings< InitialStateParameterType > >
                        initialStateSettings =  std::dynamic_pointer_cast<
                        ArcWiseInitialTranslationalStateEstimatableParameterSettings< InitialStateParameterType > >(
                            parameterSettings );

                if( initialStateSettings->isStateSet_ )
                {
                    initialStateParameterToEstimate = std::make_shared< ArcWiseInitialTranslationalStateParameter<
                            InitialStateParameterType > >(
                                initialStateSettings->parameterType_.second.first,
                                initialStateSettings->arcStartTimes_,
                                initialStateSettings->initialStateValue_,
                                initialStateSettings->centralBodies_,
                                initialStateSettings->frameOrientation_ );
                }
                else
                {
                    initialStateParameterToEstimate = std::make_shared< ArcWiseInitialTranslationalStateParameter<
                            InitialStateParameterType > >(
                                initialStateSettings->parameterType_.second.first, initialStateSettings->arcStartTimes_,
                                propagators::getInitialArcWiseStateOfBody< double, InitialStateParameterType >(
                                    initialStateSettings->parameterType_.second.first,
                                    initialStateSettings->centralBodies_, bodies,
                                    initialStateSettings->arcStartTimes_ ),
                                initialStateSettings->centralBodies_, initialStateSettings->frameOrientation_ );
                }
            }
            break;
        case initial_rotational_body_state:

            // Check consistency of input.
            if( std::dynamic_pointer_cast<
                    InitialRotationalStateEstimatableParameterSettings< InitialStateParameterType > >(
                        parameterSettings ) == nullptr )
            {
                throw std::runtime_error( "Error when making body initial state parameter, settings type is incompatible" );
            }
            else
            {
                std::shared_ptr< InitialRotationalStateEstimatableParameterSettings< InitialStateParameterType > >
                        initialStateSettings = std::dynamic_pointer_cast<
                        InitialRotationalStateEstimatableParameterSettings< InitialStateParameterType > >(
                            parameterSettings );

                Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > initialRotationalState;

                // If initial time is not defined, use preset initial state
                if( ! ( initialStateSettings->initialTime_ == initialStateSettings->initialTime_  ) )
                {
                    initialRotationalState = initialStateSettings->initialStateValue_;


                }
                // Compute initial state from environment
                else
                {
                    initialRotationalState = propagators::getInitialRotationalStateOfBody
                            < double, InitialStateParameterType >(
                                initialStateSettings->parameterType_.second.first, initialStateSettings->baseOrientation_,
                                bodies, initialStateSettings->initialTime_ );

                }

                // Create rotational state estimation interface object
                initialStateParameterToEstimate =
                        std::make_shared< InitialRotationalStateParameter< InitialStateParameterType > >(
                            initialStateSettings->parameterType_.second.first, initialRotationalState,
                            std::bind( &Body::getBodyInertiaTensor,
                                       bodies.at( initialStateSettings->parameterType_.second.first ) ),
                            initialStateSettings->baseOrientation_ );
            }
            break;
        case initial_mass_state:
        {
            // Check consistency of input.
            if( std::dynamic_pointer_cast<
                    InitialMassEstimatableParameterSettings< InitialStateParameterType > >(
                        parameterSettings ) == nullptr )
            {
                throw std::runtime_error( "Error when making body initial mass state parameter, settings type is incompatible" );
            }
            else
            {
                std::shared_ptr< InitialMassEstimatableParameterSettings< InitialStateParameterType > >
                        initialStateSettings = std::dynamic_pointer_cast<
                        InitialMassEstimatableParameterSettings< InitialStateParameterType > >(
                            parameterSettings );

                InitialStateParameterType initialMass = initialStateSettings->initialStateValue_;
                initialStateParameterToEstimate =
                        std::make_shared< InitialMassStateParameter< InitialStateParameterType > >(
                            initialStateSettings->parameterType_.second.first,
                            ( Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 >( 1 ) << initialMass ).finished( ) );
            }
            break;
        }

        default:
            std::string errorMessage = "Error, could not create parameter for initial state of type " +
                    std::to_string( parameterSettings->parameterType_.first );
            throw std::runtime_error( errorMessage );
        }
    }

    if( parameterSettings->customPartialSettings_.size( ) != 0 )
    {
        initialStateParameterToEstimate->setCustomPartialSettings( parameterSettings->customPartialSettings_ );
    }

    return initialStateParameterToEstimate;
}

//! Function to create an interface object for estimating a parameter defined by a single double value
/*!
 * Function to create an interface object for estimating a parameter defined by a single double value
 * \param doubleParameterName Object defining the parameter interface that is to be created.
 * \param bodies Map of body objects containing the fll simulation environment.
 * \param propagatorSettings Object defining all settigns for the propagator; empty by default (only required for
 * selected parameters).
 * \return Interface object for estimating parameter.
 */
template< typename InitialStateParameterType, typename TimeType >
std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > createDoubleParameterToEstimate(
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSettings >& doubleParameterName,
        const SystemOfBodies& bodies, const std::shared_ptr< propagators::PropagatorSettings< InitialStateParameterType > > propagatorSettings =
        std::shared_ptr< propagators::PropagatorSettings< InitialStateParameterType > >( ) )
{
    using namespace simulation_setup;
    using namespace ephemerides;
    using namespace gravitation;
    using namespace estimatable_parameters;


    std::shared_ptr< EstimatableParameter< double > > doubleParameterToEstimate;

    // Check input consistency.
    if( isDoubleParameter( doubleParameterName->parameterType_.first ) != true )
    {
        std::string errorMessage = "Error when requesting to make double parameter " +
                std::to_string( doubleParameterName->parameterType_.first ) + " of " +
                doubleParameterName->parameterType_.second.first +
                ", parameter is not a double parameter ";
        throw std::runtime_error( errorMessage );
    }
    else
    {
        // Check if body associated with parameter exists.
        std::string currentBodyName = doubleParameterName->parameterType_.second.first;
        std::shared_ptr< Body > currentBody;

        if( ( currentBodyName != "global_metric" ) && ( currentBodyName != "" ) && ( bodies.count( currentBodyName ) == 0 ) )
        {
            std::string errorMessage = "Error when creating parameters to estimate, body " +
                    currentBodyName + "  not in system of bodies " +
                    std::to_string( doubleParameterName->parameterType_.first );
            throw std::runtime_error( errorMessage );
        }
        else if( ( currentBodyName != "" ) && ( currentBodyName != "global_metric" ) )
        {
            currentBody = bodies.at( currentBodyName );
        }

        // Identify parameter type.
        switch( doubleParameterName->parameterType_.first )
        {
        case gravitational_parameter:
        {
            if( currentBody->getGravityFieldModel( )== nullptr )
            {
                std::string errorMessage = "Error, body " +
                        currentBodyName + " has no gravity field, cannot estimate gravitational parameter.";
                throw std::runtime_error( errorMessage );
            }
            else
            {
                std::shared_ptr< GravityFieldModel > gravityFieldModel = currentBody->getGravityFieldModel( );
                doubleParameterToEstimate = std::make_shared< GravitationalParameter >
                        ( gravityFieldModel, currentBodyName );
            }
            break;
        }
        case radiation_pressure_coefficient:
        {
            if( currentBody->getRadiationPressureTargetModel( ) == nullptr )
            {
                std::string errorMessage = "Error, no radiation pressure target model found in body " +
                        currentBodyName + " when making Cr parameter.";
                throw std::runtime_error( errorMessage );
            }
            else if( std::dynamic_pointer_cast< electromagnetism::CannonballRadiationPressureTargetModel >( currentBody->getRadiationPressureTargetModel( ) ) == nullptr )
            {
                std::string errorMessage = "Error, no cannonball radiation pressure target model found in body " +
                                           currentBodyName + " when making Cr parameter.";
                throw std::runtime_error( errorMessage );
            }
            else
            {
                doubleParameterToEstimate = std::make_shared< RadiationPressureCoefficient >(
                    std::dynamic_pointer_cast< electromagnetism::CannonballRadiationPressureTargetModel >( currentBody->getRadiationPressureTargetModel( ) ),
                            currentBodyName );
            }
            break;
        }
        case constant_rotation_rate:
        {
            if( std::dynamic_pointer_cast< SimpleRotationalEphemeris >( currentBody->getRotationalEphemeris( ) ) == nullptr )
            {
                std::string errorMessage = "Warning, no simple rotational ephemeris present in body " + currentBodyName +
                        " when making constant rotation rate parameter";
                throw std::runtime_error( errorMessage );
            }
            else
            {
                doubleParameterToEstimate = std::make_shared< RotationRate >(
                            std::dynamic_pointer_cast< ephemerides::SimpleRotationalEphemeris >
                            ( currentBody->getRotationalEphemeris( ) ), currentBodyName );
            }
            break;
        }
        case constant_drag_coefficient:
        {
            if( currentBody->getAerodynamicCoefficientInterface( ) == nullptr )
            {
                std::string errorMessage = "Error, body " +
                        currentBodyName + " has no coefficient interface, cannot estimate constant drag coefficient.";
                throw std::runtime_error( errorMessage );
            }
            else if( std::dynamic_pointer_cast< aerodynamics::CustomAerodynamicCoefficientInterface >(
                         currentBody->getAerodynamicCoefficientInterface( ) ) == nullptr )
            {
                std::string errorMessage = "Error, body " +
                        currentBodyName + " has no custom coefficient interface, cannot estimate constant drag coefficient.";
                throw std::runtime_error( errorMessage );
            }
            else
            {
                doubleParameterToEstimate = std::make_shared< ConstantDragCoefficient >
                        ( std::dynamic_pointer_cast< aerodynamics::CustomAerodynamicCoefficientInterface >(
                              currentBody->getAerodynamicCoefficientInterface( ) ), currentBodyName );
            }
            break;
        }
        case ppn_parameter_gamma:
        {
            doubleParameterToEstimate = std::make_shared< PPNParameterGamma >( relativity::ppnParameterSet );
            break;
        }
        case ppn_parameter_beta:
        {
            doubleParameterToEstimate = std::make_shared< PPNParameterBeta >( relativity::ppnParameterSet );
            break;
        }
        case equivalence_principle_lpi_violation_parameter:
        {
            doubleParameterToEstimate = std::make_shared< EquivalencePrincipleLpiViolationParameter >( );
            break;
        }
        case direct_dissipation_tidal_time_lag:
        {
            if( propagatorSettings == nullptr )
            {
                throw std::runtime_error( "Error when creating direct_dissipation_tidal_time_lag parameter, no propagatorSettings provided." );
            }

            // Check input consistency
            std::shared_ptr< DirectTidalTimeLagEstimatableParameterSettings > dissipationTimeLagSettings =
                    std::dynamic_pointer_cast< DirectTidalTimeLagEstimatableParameterSettings >( doubleParameterName );
            if( dissipationTimeLagSettings == nullptr )
            {
                throw std::runtime_error( "Error, expected dissipation time lag parameter settings." );
            }
            else
            {
                std::vector< std::shared_ptr< basic_astrodynamics::AccelerationModel3d > > associatedAccelerationModels =
                        getAccelerationModelsListForParametersFromBase< InitialStateParameterType, TimeType >( propagatorSettings, doubleParameterName );
                std::vector< std::shared_ptr< DirectTidalDissipationAcceleration > > associatedTidalAccelerationModels;
                for( unsigned int i = 0; i < associatedAccelerationModels.size( ); i++ )
                {
                    // Create parameter object
                    if( std::dynamic_pointer_cast< DirectTidalDissipationAcceleration >( associatedAccelerationModels.at( i ) )
                            != nullptr )
                    {
                        associatedTidalAccelerationModels.push_back(
                                    std::dynamic_pointer_cast< DirectTidalDissipationAcceleration >( associatedAccelerationModels.at( i ) ) );
                    }
                    else
                    {
                        throw std::runtime_error(
                                    "Error, expected DirectTidalDissipationAcceleration in list when creating direct_dissipation_tidal_time_lag parameter" );
                    }
                }
                doubleParameterToEstimate = std::make_shared< DirectTidalTimeLag >(
                            associatedTidalAccelerationModels, currentBodyName, dissipationTimeLagSettings->deformingBodies_ );
            }
            break;
        }
        case mean_moment_of_inertia:
        {
            if( currentBody == nullptr )
            {
                std::string errorMessage = "Error, body is nullptr when making mean moment of inertia parameter.";
                throw std::runtime_error( errorMessage );
            }
            else if( std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( currentBody->getGravityFieldModel( ) ) == nullptr )
            {
                std::string errorMessage = "Error, body gravity field is not spherical harmonic when making mean moment of inertia parameter.";
                throw std::runtime_error( errorMessage );
            }
            else
            {
                auto gravityFieldModel =
                        std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( currentBody->getGravityFieldModel( ) );
                doubleParameterToEstimate = std::make_shared< MeanMomentOfInertiaParameter >
                        ( std::bind( &SphericalHarmonicsGravityField::getScaledMeanMomentOfInertia, gravityFieldModel ),
                          std::bind( &SphericalHarmonicsGravityField::setScaledMeanMomentOfInertia, gravityFieldModel, std::placeholders::_1 ),
                          currentBodyName );
            }
            break;
        }
        case core_factor:
        {
            if( std::dynamic_pointer_cast< PlanetaryRotationModel >( currentBody->getRotationalEphemeris( ) ) == nullptr )
            {
                std::string errorMessage = "Warning, no full planetary rotational ephemeris" + currentBodyName +
                        " when making free core parameter";
                throw std::runtime_error( errorMessage );
            }
            else
            {
                doubleParameterToEstimate = std::make_shared< CoreFactor >
                        ( std::dynamic_pointer_cast< PlanetaryRotationModel >( currentBody->getRotationalEphemeris( ) ),
                          currentBodyName );

            }
            break;
        }
        case free_core_nutation_rate:
        {
            if( std::dynamic_pointer_cast< PlanetaryRotationModel >( currentBody->getRotationalEphemeris( ) ) == nullptr )
            {
                std::string errorMessage = "Warning, no full planetary rotational ephemeris" + currentBodyName +
                        " when making free core nutation rate parameter";
                throw std::runtime_error( errorMessage );
            }
            else
            {
                doubleParameterToEstimate = std::make_shared< FreeCoreNutationRate >
                        ( std::dynamic_pointer_cast< PlanetaryRotationModel > ( currentBody->getRotationalEphemeris( ) ), currentBodyName);

            }
            break;
        }
        case scaled_longitude_libration_amplitude:
        {
            if( std::dynamic_pointer_cast< SynchronousRotationalEphemeris >( currentBody->getRotationalEphemeris( ) ) == nullptr )
            {
                std::string errorMessage = "Warning, no synchronous rotation model present in body " + currentBodyName +
                        " when making longitude libration parameter";
                throw std::runtime_error( errorMessage );
            }
            else
            {
                std::shared_ptr< LongitudeLibrationCalculator > longitudeLibrationCalculator =
                        std::dynamic_pointer_cast< SynchronousRotationalEphemeris >( currentBody->getRotationalEphemeris( ) )->
                        getLongitudeLibrationCalculator( );

                if( std::dynamic_pointer_cast< DirectLongitudeLibrationCalculator >( longitudeLibrationCalculator ) == nullptr )
                {
                    std::string errorMessage = "Warning, no direct libration model " + currentBodyName +
                            " when making scaled longitude libration parameter";
                    throw std::runtime_error( errorMessage );
                }
                else
                {

                    doubleParameterToEstimate = std::make_shared< ScaledLongitudeLibrationAmplitude >
                            ( std::dynamic_pointer_cast< DirectLongitudeLibrationCalculator >( longitudeLibrationCalculator ),
                              currentBodyName );
                }

            }
            break;
        }
        case constant_thrust_magnitude_parameter:
        {
            if( currentBody->getVehicleSystems( ) == nullptr )
            {
                throw std::runtime_error( "Error when creating constant thrust magnitude for body " + currentBodyName +
                                          ", body has no vehicle systems" );
            }
            else
            {
                if( currentBody->getVehicleSystems( )->getEngineModels( ).count(
                            doubleParameterName->parameterType_.second.second ) == 0 )
                {
                    throw std::runtime_error( "Error when creating constant thrust magnitude for engine " +
                                              doubleParameterName->parameterType_.second.second + " on body " +
                                              currentBodyName + ", engine does not exist" );
                }
                else
                {
                    std::shared_ptr< propulsion::ThrustMagnitudeWrapper > thrustWrapper =
                            currentBody->getVehicleSystems( )->getEngineModels( ).at(
                                                       doubleParameterName->parameterType_.second.second )->getThrustMagnitudeWrapper( );
                    if( thrustWrapper == nullptr )
                    {
                        throw std::runtime_error( "Error when creating constant thrust magnitude for engine " +
                                                  doubleParameterName->parameterType_.second.second + " on body " +
                                                  currentBodyName + ", engine does not have thrust magnitude model." );
                    }
                    else
                    {
                        if( std::dynamic_pointer_cast< propulsion::ConstantThrustMagnitudeWrapper >( thrustWrapper ) == nullptr )
                        {
                            throw std::runtime_error( "Error when creating constant thrust magnitude for engine " +
                                                      doubleParameterName->parameterType_.second.second + " on body " +
                                                      currentBodyName + ", engine thrust magnitude model does not support constant thrust." );
                        }
                        else
                        {
                            doubleParameterToEstimate = std::make_shared< ConstantThrustMagnitudeParameter >
                                    ( std::dynamic_pointer_cast< propulsion::ConstantThrustMagnitudeWrapper >( thrustWrapper ),
                                      currentBodyName, doubleParameterName->parameterType_.second.second );
                        }
                    }
                }
            }

            break;
        }
        case constant_specific_impulse:
        {
            if( currentBody->getVehicleSystems( ) == nullptr )
            {
                throw std::runtime_error( "Error when creating constant specific impulse for body " + currentBodyName +
                                          ", body has no vehicle systems" );
            }
            else
            {
                if( currentBody->getVehicleSystems( )->getEngineModels( ).count(
                            doubleParameterName->parameterType_.second.second ) == 0 )
                {
                    throw std::runtime_error( "Error when creating constant specific impulse for engine " +
                                              doubleParameterName->parameterType_.second.second + " on body " +
                                              currentBodyName + ", engine does not exist" );
                }
                else
                {
                    std::shared_ptr< propulsion::ThrustMagnitudeWrapper > thrustWrapper =
                            currentBody->getVehicleSystems( )->getEngineModels( ).at(
                                                       doubleParameterName->parameterType_.second.second )->getThrustMagnitudeWrapper( );
                    if( thrustWrapper == nullptr )
                    {
                        throw std::runtime_error( "Error when creating constant specific impulse for engine " +
                                                  doubleParameterName->parameterType_.second.second + " on body " +
                                                  currentBodyName + ", engine does not have thrust magnitude model." );
                    }
                    else
                    {
                        if( std::dynamic_pointer_cast< propulsion::ConstantThrustMagnitudeWrapper >( thrustWrapper ) == nullptr )
                        {
                            throw std::runtime_error( "Error when creating constant specific impulse for engine " +
                                                      doubleParameterName->parameterType_.second.second + " on body " +
                                                      currentBodyName + ", engine thrust magnitude model does not support constant thrust." );
                        }
                        else
                        {
                            doubleParameterToEstimate = std::make_shared< ConstantSpecificImpulseParameter< propulsion::ConstantThrustMagnitudeWrapper  > >
                                    ( std::dynamic_pointer_cast< propulsion::ConstantThrustMagnitudeWrapper >( thrustWrapper ),
                                      currentBodyName, doubleParameterName->parameterType_.second.second );
                        }
                    }
                }
            }

            break;
        }
        case inverse_tidal_quality_factor:
        {
            if( propagatorSettings == nullptr )
            {
                throw std::runtime_error( "Error when creating inverse_tidal_quality_factor parameter, no propagatorSettings provided." );
            }

            // Check input consistency
            std::shared_ptr< InverseTidalQualityFactorEstimatableParameterSettings > qualityFactorSettings =
                    std::dynamic_pointer_cast< InverseTidalQualityFactorEstimatableParameterSettings >( doubleParameterName );
            if( qualityFactorSettings == nullptr )
            {
                throw std::runtime_error( "Error, expected inverse tidal quality factor parameter settings." );
            }
            else
            {
                std::vector< std::shared_ptr< basic_astrodynamics::AccelerationModel3d > > associatedAccelerationModels =
                        getAccelerationModelsListForParametersFromBase< InitialStateParameterType, TimeType >( propagatorSettings, doubleParameterName );
                std::vector< std::shared_ptr< DirectTidalDissipationAcceleration > > associatedTidalAccelerationModels;
                for( unsigned int i = 0; i < associatedAccelerationModels.size( ); i++ )
                {
                    // Create parameter object
                    if( std::dynamic_pointer_cast< DirectTidalDissipationAcceleration >( associatedAccelerationModels.at( i ) )
                        != nullptr )
                    {
                        associatedTidalAccelerationModels.push_back(
                                std::dynamic_pointer_cast< DirectTidalDissipationAcceleration >( associatedAccelerationModels.at( i ) ) );
                    }
                    else
                    {
                        throw std::runtime_error(
                                "Error, expected DirectTidalDissipationAcceleration in list when creating inverse_tidal_quality_factor parameter" );
                    }
                }
                doubleParameterToEstimate = std::make_shared< InverseTidalQualityFactor >(
                        associatedTidalAccelerationModels, currentBodyName, qualityFactorSettings->deformingBodies_ );
            }
            break;
        }
        case yarkovsky_parameter:
        {
            if( propagatorSettings == nullptr )
            {
                throw std::runtime_error( "Error when creating yarkovsky_parameter parameter, no propagatorSettings provided." );
            }

            std::vector< std::shared_ptr< basic_astrodynamics::AccelerationModel3d > > associatedAccelerationModels =
                getAccelerationModelsListForParametersFromBase< InitialStateParameterType, TimeType >( propagatorSettings, doubleParameterName );
            std::vector< std::shared_ptr< electromagnetism::YarkovskyAcceleration > > associatedYarkovskyAccelerationModels;
            for( unsigned int i = 0; i < associatedAccelerationModels.size( ); i++ )
            {
                // Create parameter object
                if( std::dynamic_pointer_cast< electromagnetism::YarkovskyAcceleration >( associatedAccelerationModels.at( i ) )
                    != nullptr )
                {
                    associatedYarkovskyAccelerationModels.push_back(
                        std::dynamic_pointer_cast< electromagnetism::YarkovskyAcceleration >( associatedAccelerationModels.at( i ) ) );
                }
                else
                {
                    throw std::runtime_error(
                        "Error, expected YarkovskyAcceleration in list when creating yarkovsky_parameter parameter" );
                }
            }
            if( associatedYarkovskyAccelerationModels.size( ) != 1 )
            {
                throw std::runtime_error(
                    "Error, expected single YarkovskyAcceleration in list when creating yarkovsky_parameter parameter, found " +
                    std::to_string( associatedYarkovskyAccelerationModels.size( ) ) );
            }
            doubleParameterToEstimate = std::make_shared< YarkovskyParameter >(
                associatedYarkovskyAccelerationModels.at( 0 ), currentBodyName, doubleParameterName->parameterType_.second.second );
            break;
        }
        default:
            throw std::runtime_error( "Warning, this double parameter has not yet been implemented when making parameters" );
            break;
        }
    }

    if( doubleParameterName->customPartialSettings_.size( ) != 0 )
    {
        doubleParameterToEstimate->setCustomPartialSettings( doubleParameterName->customPartialSettings_ );
    }

    return doubleParameterToEstimate;
}

//! Function to create an interface object for estimating a parameter defined by a list of double values
/*!
 * Function to create an interface object for estimating a parameter defined by a list of single double values
 * \param vectorParameterName Object defining the parameter interface that is to be created.
 * \param bodies Map of body objects containing the fll simulation environment.
 * \param propagatorSettings Object defining all settigns for the propagator; empty by default (only required for
 * selected parameters).
 * \return Interface object for estimating parameter.
 */
template< typename InitialStateParameterType, typename TimeType >
std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > createVectorParameterToEstimate(
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSettings >& vectorParameterName,
        const SystemOfBodies& bodies,
        const std::shared_ptr< propagators::PropagatorSettings< InitialStateParameterType > > propagatorSettings =
        std::shared_ptr< propagators::PropagatorSettings< InitialStateParameterType > >( ) )
{
    using namespace simulation_setup;
    using namespace ephemerides;
    using namespace gravitation;
    using namespace estimatable_parameters;

    std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > vectorParameterToEstimate;

    // Check input consistency.
    if( isDoubleParameter( vectorParameterName->parameterType_.first ) != false )
    {
        std::string errorMessage = "Error when requesting to make vector parameter " +
                std::to_string( vectorParameterName->parameterType_.first ) +
                " of  " + std::string( vectorParameterName->parameterType_.second.first ) +
                ", parameter is not a vector parameter ";
        throw std::runtime_error( errorMessage );
    }
    else
    {
        // Check if body associated with parameter exists.
        std::string currentBodyName = vectorParameterName->parameterType_.second.first;
        std::shared_ptr< Body > currentBody;
        if( ( currentBodyName != "" ) && ( bodies.count( currentBodyName ) == 0 ) )
        {
            std::string errorMessage = "Warning when creating parameters to estimate, body " +
                    currentBodyName  + "not in system of bodies " + std::to_string( vectorParameterName->parameterType_.first );
            throw std::runtime_error( errorMessage );
        }
        else if( ( currentBodyName != "" ) )
        {
            currentBody = bodies.at( currentBodyName );
        }

        // Identify parameter type.
        switch( vectorParameterName->parameterType_.first )
        {
        case constant_additive_observation_bias:
        {
            std::shared_ptr< ConstantObservationBiasEstimatableParameterSettings > biasSettings =
                    std::dynamic_pointer_cast< ConstantObservationBiasEstimatableParameterSettings >( vectorParameterName );
            if( biasSettings == nullptr )
            {
                throw std::runtime_error( "Error when creating constant observation bias, input is inconsistent" );
            }
            else
            {
                vectorParameterToEstimate = std::make_shared< ConstantObservationBiasParameter >(
                            std::function< Eigen::VectorXd( ) >( ),
                            std::function< void( const Eigen::VectorXd& ) >( ),
                            biasSettings->linkEnds_.linkEnds_, biasSettings->observableType_, true );
            }
            break;
        }
        case constant_relative_observation_bias:
        {
            std::shared_ptr< ConstantObservationBiasEstimatableParameterSettings > biasSettings =
                    std::dynamic_pointer_cast< ConstantObservationBiasEstimatableParameterSettings >( vectorParameterName );
            if( biasSettings == nullptr )
            {
                throw std::runtime_error( "Error when creating constant observation bias, input is inconsistent" );
            }
            else
            {
                vectorParameterToEstimate = std::make_shared< ConstantObservationBiasParameter >(
                            std::function< Eigen::VectorXd( ) >( ),
                            std::function< void( const Eigen::VectorXd& ) >( ),
                            biasSettings->linkEnds_.linkEnds_, biasSettings->observableType_, false );
            }
            break;
        }
        case arcwise_constant_additive_observation_bias:
        {
            std::shared_ptr< ArcWiseConstantObservationBiasEstimatableParameterSettings > biasSettings =
                    std::dynamic_pointer_cast< ArcWiseConstantObservationBiasEstimatableParameterSettings >( vectorParameterName );
            if( biasSettings == nullptr )
            {
                throw std::runtime_error( "Error when creating arcwise constant observation bias, input is inconsistent" );
            }
            else
            {
                vectorParameterToEstimate = std::make_shared< ArcWiseObservationBiasParameter >(
                            biasSettings->arcStartTimes_,
                            std::function< std::vector< Eigen::VectorXd >( ) >( ),
                            std::function< void( const std::vector< Eigen::VectorXd >& ) >( ),
                            observation_models::getLinkEndIndicesForLinkEndTypeAtObservable(
                                biasSettings->observableType_, biasSettings->linkEndForTime_, biasSettings->linkEnds_.size( ) ).at( 0 ),
                            biasSettings->linkEnds_.linkEnds_, biasSettings->observableType_, true );
            }
            break;
        }
        case arcwise_constant_relative_observation_bias:
        {
            std::shared_ptr< ArcWiseConstantObservationBiasEstimatableParameterSettings > biasSettings =
                    std::dynamic_pointer_cast< ArcWiseConstantObservationBiasEstimatableParameterSettings >( vectorParameterName );
            if( biasSettings == nullptr )
            {
                throw std::runtime_error( "Error when creating arcwise constant relative observation bias, input is inconsistent" );
            }
            else
            {
                vectorParameterToEstimate = std::make_shared< ArcWiseObservationBiasParameter >(
                            biasSettings->arcStartTimes_,
                            std::function< std::vector< Eigen::VectorXd >( ) >( ),
                            std::function< void( const std::vector< Eigen::VectorXd >& ) >( ),
                            observation_models::getLinkEndIndicesForLinkEndTypeAtObservable(
                                biasSettings->observableType_, biasSettings->linkEndForTime_, biasSettings->linkEnds_.size( ) ).at( 0 ),
                            biasSettings->linkEnds_.linkEnds_, biasSettings->observableType_, false );
            }
            break;
        }
        case constant_time_drift_observation_bias:
        {
            std::shared_ptr< ConstantTimeDriftBiasEstimatableParameterSettings > biasSettings =
                    std::dynamic_pointer_cast< ConstantTimeDriftBiasEstimatableParameterSettings >( vectorParameterName );
            if( biasSettings == nullptr )
            {
                throw std::runtime_error( "Error when creating constant time drift bias, input is inconsistent" );
            }
            else
            {
                vectorParameterToEstimate = std::make_shared< ConstantTimeDriftBiasParameter >(
                        std::function< Eigen::VectorXd( ) >( ),
                        std::function< void( const Eigen::VectorXd& ) >( ),
                        observation_models::getLinkEndIndicesForLinkEndTypeAtObservable(
                                biasSettings->observableType_, biasSettings->linkEndForTime_, biasSettings->linkEnds_.size( ) ).at( 0 ),
                        biasSettings->linkEnds_, biasSettings->observableType_, biasSettings->referenceEpoch_ );
            }
            break;
        }
        case arc_wise_time_drift_observation_bias:
        {
            std::shared_ptr< ArcWiseTimeDriftBiasEstimatableParameterSettings > timeBiasSettings =
                    std::dynamic_pointer_cast< ArcWiseTimeDriftBiasEstimatableParameterSettings >( vectorParameterName );
            if( timeBiasSettings == nullptr )
            {
                throw std::runtime_error( "Error when creating arcwise time drift bias, input is inconsistent" );
            }
            else
            {
                vectorParameterToEstimate = std::make_shared< ArcWiseTimeDriftBiasParameter >(
                        timeBiasSettings->arcStartTimes_,
                        std::function< std::vector< Eigen::VectorXd >( ) >( ),
                        std::function< void( const std::vector< Eigen::VectorXd >& ) >( ),
                        observation_models::getLinkEndIndicesForLinkEndTypeAtObservable(
                                timeBiasSettings->observableType_, timeBiasSettings->linkEndForTime_, timeBiasSettings->linkEnds_.size( ) ).at( 0 ),
                        timeBiasSettings->linkEnds_, timeBiasSettings->observableType_, timeBiasSettings->referenceEpochs_ );
            }
            break;
        }
        case constant_time_observation_bias:
        {
            std::shared_ptr< ConstantTimeBiasEstimatableParameterSettings > biasSettings =
                    std::dynamic_pointer_cast< ConstantTimeBiasEstimatableParameterSettings >( vectorParameterName );
            if( biasSettings == nullptr )
            {
                throw std::runtime_error( "Error when creating constant time bias, input is inconsistent" );
            }
            else
            {
                vectorParameterToEstimate = std::make_shared< ConstantTimeBiasParameter >(
                        std::function< Eigen::VectorXd( ) >( ),
                        std::function< void( const Eigen::VectorXd& ) >( ),
                        biasSettings->linkEndForTime_,
                        biasSettings->linkEnds_,
                        biasSettings->observableType_ );
            }
            break;
        }
        case arc_wise_time_observation_bias:
        {
            std::shared_ptr< ArcWiseTimeBiasEstimatableParameterSettings > timeBiasSettings =
                    std::dynamic_pointer_cast< ArcWiseTimeBiasEstimatableParameterSettings >( vectorParameterName );
            if( timeBiasSettings == nullptr )
            {
                throw std::runtime_error( "Error when creating arcwise time bias, input is inconsistent" );
            }
            else
            {
                vectorParameterToEstimate = std::make_shared< ArcWiseTimeBiasParameter >(
                        timeBiasSettings->arcStartTimes_,
                        std::function< std::vector< Eigen::VectorXd >( ) >( ),
                        std::function< void( const std::vector< Eigen::VectorXd >& ) >( ),
                        timeBiasSettings->linkEndForTime_,
                        timeBiasSettings->linkEnds_, timeBiasSettings->observableType_ );
            }
            break;
        }
        case rotation_pole_position:
            if( std::dynamic_pointer_cast< SimpleRotationalEphemeris >( currentBody->getRotationalEphemeris( ) ) == nullptr )
            {
                std::string errorMessage = "Warning, no simple rotational ephemeris present in body " + currentBodyName +
                        " when making constant rotation orientation parameter";
                throw std::runtime_error( errorMessage );
            }
            else
            {
                vectorParameterToEstimate = std::make_shared< ConstantRotationalOrientation >
                        ( std::dynamic_pointer_cast< ephemerides::SimpleRotationalEphemeris >
                          ( currentBody->getRotationalEphemeris( ) ), currentBodyName );

            }
            break;

        case spherical_harmonics_cosine_coefficient_block:
        {
            std::shared_ptr< GravityFieldModel > gravityField = currentBody->getGravityFieldModel( );
            std::shared_ptr< SphericalHarmonicsGravityField > shGravityField =
                    std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( gravityField );
            if( shGravityField == nullptr )
            {
                std::string errorMessage = "Error, requested spherical harmonic cosine coefficient block parameter of " +
                        std::string( vectorParameterName->parameterType_.second.first ) +
                        ", but body does not have a spherical harmonic gravity field.";
                throw std::runtime_error( errorMessage );
            }
            else
            {
                // Check if spherical harmonic gravity field is static or time-dependent; set associated
                // functions accordingly
                std::shared_ptr< TimeDependentSphericalHarmonicsGravityField > timeDependentShField =
                        std::dynamic_pointer_cast< TimeDependentSphericalHarmonicsGravityField >( shGravityField );

                std::function< Eigen::MatrixXd( ) > getCosineCoefficientsFunction;
                std::function< void( Eigen::MatrixXd ) > setCosineCoefficientsFunction;

                if( timeDependentShField == nullptr )
                {
                    getCosineCoefficientsFunction = std::bind( &SphericalHarmonicsGravityField::getCosineCoefficients,
                                                               shGravityField );
                    setCosineCoefficientsFunction = std::bind( &SphericalHarmonicsGravityField::setCosineCoefficients,
                                                               shGravityField, std::placeholders::_1 );
                }
                else
                {
                    getCosineCoefficientsFunction = std::bind(
                                &TimeDependentSphericalHarmonicsGravityField::getNominalCosineCoefficients,
                                timeDependentShField );
                    setCosineCoefficientsFunction = std::bind(
                                &TimeDependentSphericalHarmonicsGravityField::setNominalCosineCoefficients,
                                timeDependentShField, std::placeholders::_1 );
                }

                // Create cosine coefficients estimation object.
                std::shared_ptr< SphericalHarmonicEstimatableParameterSettings > blockParameterSettings =
                        std::dynamic_pointer_cast< SphericalHarmonicEstimatableParameterSettings >( vectorParameterName );
                if( blockParameterSettings != nullptr )
                {
                    vectorParameterToEstimate = std::make_shared< SphericalHarmonicsCosineCoefficients >(
                                getCosineCoefficientsFunction,
                                setCosineCoefficientsFunction,
                                blockParameterSettings->blockIndices_,
                                vectorParameterName->parameterType_.second.first );
                }
                else
                {
                    throw std::runtime_error( "Error, expected SphericalHarmonicEstimatableParameterSettings for cosine coefficients" );
                }
            }
            break;
        }
        case spherical_harmonics_sine_coefficient_block:
        {
            std::shared_ptr< GravityFieldModel > gravityField = currentBody->getGravityFieldModel( );
            std::shared_ptr< SphericalHarmonicsGravityField > shGravityField =
                    std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( gravityField );
            if( shGravityField == nullptr )
            {
                std::string errorMessage = "Error, requested spherical harmonic sine coefficient block parameter of " +
                        std::string( vectorParameterName->parameterType_.second.first ) +
                        ", but body does not have a spherical harmonic gravity field.";
                throw std::runtime_error( errorMessage );

            }
            else
            {
                std::shared_ptr< SphericalHarmonicEstimatableParameterSettings > blockParameterSettings =
                        std::dynamic_pointer_cast< SphericalHarmonicEstimatableParameterSettings >( vectorParameterName );

                // Check if spherical harmonic gravity field is static or time-dependent; set associated
                // functions accordingly
                std::function< Eigen::MatrixXd( ) > getSineCoefficientsFunction;
                std::function< void( Eigen::MatrixXd ) > setSineCoefficientsFunction;
                std::shared_ptr< TimeDependentSphericalHarmonicsGravityField > timeDependentShField =
                        std::dynamic_pointer_cast< TimeDependentSphericalHarmonicsGravityField >( shGravityField );

                if( timeDependentShField == nullptr )
                {
                    getSineCoefficientsFunction = std::bind( &SphericalHarmonicsGravityField::getSineCoefficients,
                                                             shGravityField );
                    setSineCoefficientsFunction = std::bind( &SphericalHarmonicsGravityField::setSineCoefficients,
                                                             shGravityField, std::placeholders::_1 );
                }
                else
                {
                    getSineCoefficientsFunction = std::bind(
                                &TimeDependentSphericalHarmonicsGravityField::getNominalSineCoefficients,
                                timeDependentShField );
                    setSineCoefficientsFunction = std::bind(
                                &TimeDependentSphericalHarmonicsGravityField::setNominalSineCoefficients,
                                timeDependentShField, std::placeholders::_1 );
                }

                // Create sine coefficients estimation object.
                if( blockParameterSettings != nullptr )
                {
                    vectorParameterToEstimate = std::make_shared< SphericalHarmonicsSineCoefficients >(
                                getSineCoefficientsFunction,
                                setSineCoefficientsFunction,
                                blockParameterSettings->blockIndices_,
                                vectorParameterName->parameterType_.second.first );
                }
                else
                {
                    throw std::runtime_error( "Error, expected SphericalHarmonicEstimatableParameterSettings for sine coefficients" );
                }
            }

            break;
        }
        case ground_station_position:
        {

            if( currentBody->getGroundStationMap( ).count( vectorParameterName->parameterType_.second.second ) == 0 )
            {
                std::string errorMessage =
                        "Error, requested ground station position parameter of "
                        + vectorParameterName->parameterType_.second.first + " "
                        + vectorParameterName->parameterType_.second.second + " , but ground station was not found";
                throw std::runtime_error( errorMessage );
            }
            else
            {
                std::shared_ptr< ground_stations::GroundStationState > groundStationState =
                        currentBody->getGroundStation( vectorParameterName->parameterType_.second.second )->
                        getNominalStationState( );
                if( groundStationState == nullptr )
                {
                    std::string errorMessage =
                            "Error, requested ground station position parameter of " +
                            vectorParameterName->parameterType_.second.first + " " +
                            vectorParameterName->parameterType_.second.second +
                            "  but nominal ground station state is nullptr";
                    throw std::runtime_error( errorMessage );
                }
                else
                {
                    vectorParameterToEstimate = std::make_shared< GroundStationPosition  >(
                                groundStationState, vectorParameterName->parameterType_.second.first,
                                vectorParameterName->parameterType_.second.second );
                }
            }
            break;
        }
        case reference_point_position:
        {
            if( currentBody->getVehicleSystems( ) == nullptr )
            {
                std::string errorMessage =
                    "Error, requested reference point position parameter of "
                    + vectorParameterName->parameterType_.second.first + " "
                    + vectorParameterName->parameterType_.second.second + " , but no system models found";
                throw std::runtime_error( errorMessage );
            }
            else
            {
                vectorParameterToEstimate = std::make_shared< ReferencePointPosition  >(
                    currentBody->getVehicleSystems( ),
                    vectorParameterName->parameterType_.second.first,
                    vectorParameterName->parameterType_.second.second );
            }
            break;
        }
        case empirical_acceleration_coefficients:
        {
            if( propagatorSettings == nullptr )
            {
                throw std::runtime_error( "Error when creating empirical_acceleration_coefficients parameter, no propagatorSettings provided." );
            }

            // Check input consistency
            std::shared_ptr< EmpiricalAccelerationEstimatableParameterSettings > empiricalAccelerationSettings =
                    std::dynamic_pointer_cast< EmpiricalAccelerationEstimatableParameterSettings >( vectorParameterName );
            if( empiricalAccelerationSettings == nullptr )
            {
                throw std::runtime_error(
                            "Error when trying to make constant empirical acceleration coefficients parameter, settings type inconsistent" );
            }
            else
            {

                std::vector< std::shared_ptr< basic_astrodynamics::AccelerationModel3d > > associatedAccelerationModels =
                        getAccelerationModelsListForParametersFromBase< InitialStateParameterType, TimeType >( propagatorSettings, vectorParameterName );
                std::vector< std::shared_ptr< basic_astrodynamics::EmpiricalAcceleration > > empiricalAccelerations;
                for( unsigned int i = 0; i < associatedAccelerationModels.size( ); i++ )
                {
                    // Create parameter object
                    if( std::dynamic_pointer_cast< basic_astrodynamics::EmpiricalAcceleration >( associatedAccelerationModels.at( i ) )
                            != nullptr )
                    {
                        empiricalAccelerations.push_back(
                                    std::dynamic_pointer_cast< basic_astrodynamics::EmpiricalAcceleration >( associatedAccelerationModels.at( i ) ) );
                    }
                    else
                    {
                        throw std::runtime_error(
                                    "Error, expected EmpiricalAcceleration in list when creating empirical_acceleration_coefficients parameter" );
                    }
                }

                // Create empirical acceleration parameter
                vectorParameterToEstimate = std::make_shared< EmpiricalAccelerationCoefficientsParameter >(
                            empiricalAccelerations,
                            empiricalAccelerationSettings->parameterType_.second.first,
                            empiricalAccelerationSettings->parameterType_.second.second,
                            empiricalAccelerationSettings->componentsToEstimate_ );


            }
            break;
        }
        case arc_wise_radiation_pressure_coefficient:
        {
            // Check input consistency
            std::shared_ptr< ArcWiseRadiationPressureCoefficientEstimatableParameterSettings > radiationPressureCoefficientSettings =
                    std::dynamic_pointer_cast< ArcWiseRadiationPressureCoefficientEstimatableParameterSettings >( vectorParameterName );
            if( radiationPressureCoefficientSettings == nullptr )
            {
                throw std::runtime_error(
                            "Error when trying to make arc-wise radiation pressure coefficients parameter, settings type inconsistent" );
            }
            else
            {
                if( currentBody->getRadiationPressureTargetModel( ) == nullptr )
                {
                    std::string errorMessage = "Error, no radiation pressure interfaces found in body " +
                            currentBodyName + " when making arcwise Cr parameter.";
                    throw std::runtime_error( errorMessage );
                }
                else if( std::dynamic_pointer_cast< electromagnetism::CannonballRadiationPressureTargetModel >( currentBody->getRadiationPressureTargetModel( ) ) == nullptr )
                {
                    std::string errorMessage = "Error, no cannonball radiation pressure target model found in body " +
                                               currentBodyName + " when making arcwise Cr parameter.";
                    throw std::runtime_error( errorMessage );
                }
                else
                {
                    vectorParameterToEstimate = std::make_shared< ArcWiseRadiationPressureCoefficient >(
                        std::dynamic_pointer_cast< electromagnetism::CannonballRadiationPressureTargetModel >( currentBody->getRadiationPressureTargetModel( ) ),
                        radiationPressureCoefficientSettings->arcStartTimeList_,
                        currentBodyName );
                }
                break;
            }
            break;
        }
        case arc_wise_constant_drag_coefficient:
        {
            // Check input consistency
            std::shared_ptr< ArcWiseDragCoefficientEstimatableParameterSettings > dragCoefficientSettings =
                    std::dynamic_pointer_cast< ArcWiseDragCoefficientEstimatableParameterSettings >( vectorParameterName );
            if( dragCoefficientSettings == nullptr )
            {
                throw std::runtime_error(
                            "Error when trying to make arc-wise radiation pressure coefficients parameter, settings type inconsistent" );
            }
            else
            {
                if( currentBody->getAerodynamicCoefficientInterface( ) == nullptr )
                {
                    std::string errorMessage = "Error, no aerodynamic coefficient interfaces found in body " +
                            currentBodyName + " when making arcwise Cd parameter.";
                    throw std::runtime_error( errorMessage );
                }
                else if( std::dynamic_pointer_cast< aerodynamics::CustomAerodynamicCoefficientInterface >(
                             currentBody->getAerodynamicCoefficientInterface( ) ) == nullptr )
                {
                    std::string errorMessage = "Error, incompatible aerodynamic coefficient interfaces found in body " +
                            currentBodyName + " when making arcwise Cd parameter.";
                    throw std::runtime_error( errorMessage );
                }
                else
                {
                    vectorParameterToEstimate = std::make_shared< ArcWiseConstantDragCoefficient >(
                                std::dynamic_pointer_cast< aerodynamics::CustomAerodynamicCoefficientInterface >(
                                    currentBody->getAerodynamicCoefficientInterface( ) ),
                                dragCoefficientSettings->arcStartTimeList_,
                                currentBodyName );
                }
                break;
            }
            break;

        }
        case arc_wise_empirical_acceleration_coefficients:
        {
            if( propagatorSettings == nullptr )
            {
                throw std::runtime_error( "Error when creating arc_wise_empirical_acceleration_coefficients parameter, no propagatorSettings provided." );
            }

            // Check input consistency
            std::shared_ptr< ArcWiseEmpiricalAccelerationEstimatableParameterSettings > empiricalAccelerationSettings =
                    std::dynamic_pointer_cast< ArcWiseEmpiricalAccelerationEstimatableParameterSettings >( vectorParameterName );
            if( empiricalAccelerationSettings == nullptr )
            {
                throw std::runtime_error(
                            "Error when trying to make constant empirical acceleration coefficients parameter, settings type inconsistent" );
            }
            else
            {

                std::vector< std::shared_ptr< basic_astrodynamics::AccelerationModel3d > > associatedAccelerationModels =
                        getAccelerationModelsListForParametersFromBase< InitialStateParameterType, TimeType >( propagatorSettings, vectorParameterName );
                std::vector< std::shared_ptr< basic_astrodynamics::EmpiricalAcceleration > > empiricalAccelerations;
                for( unsigned int i = 0; i < associatedAccelerationModels.size( ); i++ )
                {
                    // Create parameter object
                    if( std::dynamic_pointer_cast< basic_astrodynamics::EmpiricalAcceleration >( associatedAccelerationModels.at( i ) )
                            != nullptr )
                    {
                        empiricalAccelerations.push_back(
                                    std::dynamic_pointer_cast< basic_astrodynamics::EmpiricalAcceleration >( associatedAccelerationModels.at( i ) ) );
                    }
                    else
                    {
                        throw std::runtime_error(
                                    "Error, expected EmpiricalAcceleration in list when creating arc_wise_empirical_acceleration_coefficients parameter" );
                    }
                }
                // Create arcwise empirical acceleration parameter
                vectorParameterToEstimate = std::make_shared< ArcWiseEmpiricalAccelerationCoefficientsParameter >(
                            empiricalAccelerations, empiricalAccelerationSettings->parameterType_.second.first,
                            empiricalAccelerationSettings->parameterType_.second.second,
                            empiricalAccelerationSettings->componentsToEstimate_, empiricalAccelerationSettings->arcStartTimeList_ );

            }

            break;
        }
        case full_degree_tidal_love_number:
        {
            // Check input consistency
            std::shared_ptr< FullDegreeTidalLoveNumberEstimatableParameterSettings > tidalLoveNumberSettings =
                    std::dynamic_pointer_cast< FullDegreeTidalLoveNumberEstimatableParameterSettings >( vectorParameterName );
            if( tidalLoveNumberSettings == nullptr )
            {
                throw std::runtime_error( "Error, expected tidal love number parameter settings." );
            }
            else
            {
                // Check consistency of body gravity field
                std::shared_ptr< GravityFieldModel > gravityField = currentBody->getGravityFieldModel( );
                std::shared_ptr< TimeDependentSphericalHarmonicsGravityField > timeDepGravityField =
                        std::dynamic_pointer_cast< TimeDependentSphericalHarmonicsGravityField >( gravityField );
                if( timeDepGravityField == nullptr )
                {
                    throw std::runtime_error(
                                "Error, requested tidal love number parameter of " +
                                vectorParameterName->parameterType_.second.first +
                                ", but body does not have a time dependent spherical harmonic gravity field." );
                }
                else if( currentBody->getGravityFieldVariationSet( ) == nullptr )
                {
                    throw std::runtime_error( "Error, requested tidal love number parameter of " +
                                              vectorParameterName->parameterType_.second.first +
                                              ", but body does not have gravity field variations" );
                }
                else
                {

                    // Get associated gravity field variation
                    std::shared_ptr< gravitation::BasicSolidBodyTideGravityFieldVariations > gravityFieldVariation =
                            std::dynamic_pointer_cast< gravitation::BasicSolidBodyTideGravityFieldVariations >(
                                currentBody->getGravityFieldVariationSet( )->getDirectTidalGravityFieldVariation(
                                    tidalLoveNumberSettings->deformingBodies_ ) );

                    // Create parameter object
                    if( gravityFieldVariation != nullptr )
                    {
                        vectorParameterToEstimate = std::make_shared< FullDegreeTidalLoveNumber >(
                                    gravityFieldVariation, currentBodyName, tidalLoveNumberSettings->degree_,
                                    tidalLoveNumberSettings->useComplexValue_ );
                    }
                    else
                    {
                        throw std::runtime_error(
                                    "Error, expected BasicSolidBodyTideGravityFieldVariations for tidal love number" );
                    }
                }
            }
            break;
        }
        case single_degree_variable_tidal_love_number:
        {
            // Check input consistency
            std::shared_ptr< SingleDegreeVariableTidalLoveNumberEstimatableParameterSettings > tidalLoveNumberSettings =
                    std::dynamic_pointer_cast< SingleDegreeVariableTidalLoveNumberEstimatableParameterSettings >( vectorParameterName );
            if( tidalLoveNumberSettings == nullptr )
            {
                throw std::runtime_error( "Error, expected variable tidal love number parameter settings " );
            }
            else
            {
                // Check consistency of body gravity field
                std::shared_ptr< TimeDependentSphericalHarmonicsGravityField > timeDepGravityField =
                        std::dynamic_pointer_cast< TimeDependentSphericalHarmonicsGravityField >(
                            currentBody->getGravityFieldModel( ) );
                if( timeDepGravityField == nullptr )
                {
                    throw std::runtime_error(
                                "Error, requested variable tidal love number parameter of " +
                                vectorParameterName->parameterType_.second.first +
                                ", but body does not have a time dependent spherical harmonic gravity field." );
                }
                else if( currentBody->getGravityFieldVariationSet( ) == nullptr )
                {
                    throw std::runtime_error( "Error, requested variable tidal love number parameter of " +
                                              vectorParameterName->parameterType_.second.first +
                                              ", but body does not have gravity field variations" );
                }
                else
                {
                    // Get associated gravity field variation
                    std::shared_ptr< gravitation::BasicSolidBodyTideGravityFieldVariations > gravityFieldVariation =
                            std::dynamic_pointer_cast< gravitation::BasicSolidBodyTideGravityFieldVariations >(
                                currentBody->getGravityFieldVariationSet( )->getDirectTidalGravityFieldVariation(
                                    tidalLoveNumberSettings->deformingBodies_ ) );

                    // Create parameter object
                    if( gravityFieldVariation != nullptr )
                    {
                        std::vector< int > orders = tidalLoveNumberSettings->orders_;
                        if( std::find( orders.begin( ), orders.end( ), 0 ) != orders.end( ) &&
                                tidalLoveNumberSettings->useComplexValue_ )
                        {
                            std::cerr << "Warning, creating parameter to estimate complex Love number at order 0, but imaginary part has no influence on dynamcis" << std::endl;
                        }
                        vectorParameterToEstimate = std::make_shared< SingleDegreeVariableTidalLoveNumber >(
                                    gravityFieldVariation, currentBodyName, tidalLoveNumberSettings->degree_,
                                    tidalLoveNumberSettings->orders_, tidalLoveNumberSettings->useComplexValue_ );
                    }
                    else
                    {
                        throw std::runtime_error(
                                    "Error, expected BasicSolidBodyTideGravityFieldVariations for variable tidal love number" );
                    }
                }
            }
            break;
        }
        case desaturation_delta_v_values:
        {
            if( propagatorSettings == nullptr )
            {
                throw std::runtime_error( "Error when creating desaturation_delta_v_values parameter, no propagatorSettings provided." );
            }
            // Check input consistency.
            std::string acceleratedBody = vectorParameterName->parameterType_.second.first;

            // Retrieve acceleration model.
            std::vector< std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > > desaturationAccelerationModels =
                    getAccelerationModelsListForParametersFromBase< InitialStateParameterType, TimeType >( propagatorSettings, vectorParameterName );

            if( desaturationAccelerationModels.size( ) == 0 )
            {
                throw std::runtime_error( "Error when making desaturation Delta V parameter, no acceleration models found in list" );

            }
            else if( desaturationAccelerationModels.size( ) > 1 )
            {
                throw std::runtime_error( "Error when making desaturation Delta V parameter, multiple acceleration models found in list" );

            }
            else
            {
                // Create desaturation deltaV values parameter.
                vectorParameterToEstimate = std::make_shared< DesaturationDeltaV >(
                            std::dynamic_pointer_cast< propulsion::MomentumWheelDesaturationThrustAcceleration >(
                                desaturationAccelerationModels.at( 0 ) ), acceleratedBody );
            }

            break;
        }
        case periodic_spin_variation:
        {
            if( std::dynamic_pointer_cast< PlanetaryRotationModel >( currentBody->getRotationalEphemeris( ) ) == nullptr )
            {
                std::string errorMessage = "Warning, no full planetary rotational ephemeris" + currentBodyName +
                        " when making periodic spin variation parameter";
                throw std::runtime_error( errorMessage );
            }
            else
            {

                vectorParameterToEstimate = std::make_shared< PeriodicSpinVariation >
                        ( std::dynamic_pointer_cast< PlanetaryRotationModel > ( currentBody->getRotationalEphemeris( ) ), currentBodyName);

            }
            break;
        }
        case polar_motion_amplitude:
        {
            if( std::dynamic_pointer_cast< PlanetaryRotationModel >( currentBody->getRotationalEphemeris( ) ) == nullptr )
            {
                std::string errorMessage = "Warning, no full planetary rotational ephemeris" + currentBodyName +
                        " when making polar motion amplitude parameter";
                throw std::runtime_error( errorMessage );
            }
            else
            {
                vectorParameterToEstimate = std::make_shared< PolarMotionAmplitude >
                        ( std::dynamic_pointer_cast< PlanetaryRotationModel > ( currentBody->getRotationalEphemeris( ) ), currentBodyName);
            }
            break;
        }
        case custom_estimated_parameter:
        {
            std::shared_ptr< CustomEstimatableParameterSettings > customParameterSettings =
                std::dynamic_pointer_cast< CustomEstimatableParameterSettings >( vectorParameterName );
            if( customParameterSettings == nullptr )
            {
                throw std::runtime_error( "Error, expected variable custom parameter settings " );
            }
            else
            {
                vectorParameterToEstimate = std::make_shared< CustomEstimatableParameter >
                    ( customParameterSettings->parameterType_.second.second,
                      customParameterSettings->parameterSize_,
                      customParameterSettings->getParameterFunction_,
                      customParameterSettings->setParameterFunction_ );
            }
            break;
        }
        default:
            std::string errorMessage = "Warning, this vector parameter (" +
                    std::to_string( vectorParameterName->parameterType_.first ) +
                    ") has not yet been implemented when making parameters";
            throw std::runtime_error( errorMessage );

            break;
        }
    }

    if( vectorParameterName->customPartialSettings_.size( ) != 0 )
    {
        vectorParameterToEstimate->setCustomPartialSettings( vectorParameterName->customPartialSettings_ );
    }

    return vectorParameterToEstimate;
}


//! Function checking whether the direct tidal parameters to be estimated are not incompatible
//! Tidal time lag and inverse of tidal quality factor cannot be simultaneously estimated for the same body and deforming bodies.
template< typename InitialStateParameterType = double >
bool checkCompatibilityDirectTidalParameters(
        const std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > >& parameterNames )
{
    using namespace tudat::estimatable_parameters;

    bool compatibleParameters = true;
    for ( unsigned int i = 0; i < parameterNames.size( ); i++ )
    {
        if ( parameterNames[ i ]->parameterType_.first == direct_dissipation_tidal_time_lag )
        {
            std::string associatedBody = parameterNames[ i ]->parameterType_.second.first;
            std::vector< std::string > deformingBodies = std::dynamic_pointer_cast< DirectTidalTimeLagEstimatableParameterSettings >( parameterNames[ i ] )->deformingBodies_;

            for ( unsigned int j = 0 ; j < parameterNames.size( ) ; j++ )
            {
                if ( i != j )
                {
                    if ( parameterNames[ j ]->parameterType_.first == inverse_tidal_quality_factor )
                    {
                        std::string associatedBody2 = parameterNames[ j ]->parameterType_.second.first;
                        std::vector< std::string > deformingBodies2 =
                                std::dynamic_pointer_cast< InverseTidalQualityFactorEstimatableParameterSettings >( parameterNames[ j ] )->deformingBodies_;

                        if ( ( associatedBody == associatedBody2 ) && ( deformingBodies == deformingBodies2 ) )
                        {
                            compatibleParameters = false;
                        }
                    }
                }
            }
        }
    }
    return compatibleParameters;
}

//! Function to create the interface object for estimating any number/type of parameters.
/*!
 *  Function to create the interface object for estimating any number/type of parameters. This can include both
 *  environmental parameters and initial dynamical states. The types of parameters are defined by the parameterNames m
 *  input variables
 *  \param parameterNames List of objects defining the parameters that are to be estimated.
 *  \param bodies Map of body objects containing the fll simulation environment.
 * \param propagatorSettings Object defining all settigns for the propagator; empty by default (only required for
 * selected parameters).
 *  \return Interface object for estimating a set of parameters.
 */
template< typename InitialStateParameterType = double, typename TimeType = double >
std::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > > createParametersToEstimate(
        const std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > >& parameterNames,
        const SystemOfBodies& bodies,
        const std::shared_ptr< propagators::PropagatorSettings< InitialStateParameterType > > propagatorSettings =
        std::shared_ptr< propagators::PropagatorSettings< InitialStateParameterType > >( ),
        const std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > >& considerParameterNames =
                std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > >( ) )

{
    using namespace tudat::estimatable_parameters;

    std::vector< std::shared_ptr< EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > >
            initialDynamicalParametersToEstimate;
    std::vector< std::shared_ptr< EstimatableParameter< double > > > doubleParametersToEstimate;
    std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > > vectorParametersToEstimate;

    // Check that tidal parameters are not inconsistent
    if ( !checkCompatibilityDirectTidalParameters( parameterNames ) )
    {
        throw std::runtime_error( "Error, tidal time lag and inverse tidal quality factor cannot be simultaneously estimated"
                                  " for the same bodies." );
    }

    // Iterate over all parameters.
    bool vectorParameterIsFound = 0;
    bool parameterOrderWarningPrinted = 0;

    for( unsigned int i = 0; i < parameterNames.size( ); i++ )
    {
        // Create initial dynamical parameters.
        if( isParameterDynamicalPropertyInitialState( parameterNames.at( i )->parameterType_.first ) )
        {
            initialDynamicalParametersToEstimate.push_back(
                        createInitialDynamicalStateParameterToEstimate< InitialStateParameterType >(
                            bodies, parameterNames.at( i ) ) );
        }
        // Create parameters defined by single double value
        else if( isDoubleParameter( parameterNames[ i ]->parameterType_.first ) == true )
        {
            doubleParametersToEstimate.push_back( createDoubleParameterToEstimate< InitialStateParameterType, TimeType >(
                                                      parameterNames[ i ], bodies, propagatorSettings ) );
            if( vectorParameterIsFound == true && parameterOrderWarningPrinted == false )
            {
                std::cerr<<"Warning when creating estimated parameters. The parameters will be ordered such that all parameters (excluding initial states) "<<
                           "defined by a single variable will be stored before those represented by a list of variables. "<<
                           "The parameter order will be different than those in your parameter settings. It is recommended that you "<<
                           "check the parameter order by calling the print_parameter_names(Python)/printEstimatableParameterEntries(C++) function"<<std::endl;
                parameterOrderWarningPrinted = true;
            }
        }
        // Create parameters defined by list of double values
        else if( isDoubleParameter( parameterNames[ i ]->parameterType_.first ) == false )
        {
            vectorParametersToEstimate.push_back( createVectorParameterToEstimate< InitialStateParameterType, TimeType >(
                                                      parameterNames[ i ], bodies, propagatorSettings ) );
            vectorParameterIsFound = true;
        }
        else
        {
            std::string errorMessage = "Error, parameter type of  " +
                    std::string( parameterNames[ i ]->parameterType_.second.first ) + "of " +
                    std::to_string( parameterNames[ i ]->parameterType_.first ) +
                    "not recognized when making estimatable parameter set.";

            throw std::runtime_error( errorMessage );
        }
    }

    std::shared_ptr< EstimatableParameterSet< InitialStateParameterType > > considerParameters;
    if ( !considerParameterNames.empty( ) )
    {
        considerParameters = createParametersToEstimate( considerParameterNames, bodies, propagatorSettings );
    }

    return std::make_shared< EstimatableParameterSet< InitialStateParameterType > >(
                doubleParametersToEstimate, vectorParametersToEstimate, initialDynamicalParametersToEstimate, considerParameters );
}

//! Function to get the multi-arc parameter equivalent of a single-arc initial state parameter
/*!
 *  Function to get the multi-arc parameter equivalent of a single-arc initial state parameter. The initial state arcs are
 *  provided as input to this function.
 *  \param singleArcParameter Single-arc parameter object for which the multi-arc equivalent is to be created.
 *  \param arcStartTimes Vector of start times for separate arcs.
 *  \return Multi-arc parameter equivalent of single-arc initial state parameter input
 */
template< typename StateScalarType >
std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > >
getAssociatedMultiArcParameter(
        const std::shared_ptr< estimatable_parameters::EstimatableParameter<
        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > singleArcParameter,
        const std::vector< double >& arcStartTimes )
{
    std::shared_ptr< estimatable_parameters::EstimatableParameter<
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > >  multiArcParameter;

    // Check state type
    switch( singleArcParameter->getParameterName( ).first )
    {
    case estimatable_parameters::initial_body_state:
    {
        // Check input consistency
        std::shared_ptr< estimatable_parameters::InitialTranslationalStateParameter< StateScalarType > >
                singleArcTranslationalStateParameter =
                std::dynamic_pointer_cast< estimatable_parameters::InitialTranslationalStateParameter< StateScalarType > >(
                    singleArcParameter );
        if( singleArcTranslationalStateParameter == nullptr )
        {
            throw std::runtime_error(
                        "Error when getting multi-arc parameter from single-arc equivalent, single-arc translational state is inconsistent " );
        }

        // Retrieve single-arc initial state
        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > singleArcInitialState =
                singleArcTranslationalStateParameter->getParameterValue( );

        // Create multi-arc initial states. First arc initial state is taken from single-arc, other initial states set to zero.
        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > multiArcInitialStates =
                Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( 6 * arcStartTimes.size( ) );
        multiArcInitialStates.segment( 0, 6 ) = singleArcInitialState;

        // Creater multi-arc parameter
        std::vector< std::string > centralBodyList;
        for( unsigned int i = 0; i < arcStartTimes.size( ); i++ )
        {
            centralBodyList.push_back( singleArcTranslationalStateParameter->getCentralBody( ) );
        }
        multiArcParameter = std::make_shared< estimatable_parameters::ArcWiseInitialTranslationalStateParameter<
                StateScalarType > >(
                    singleArcTranslationalStateParameter->getParameterName( ).second.first,
                    arcStartTimes,
                    multiArcInitialStates,
                    centralBodyList,
                    singleArcTranslationalStateParameter->getFrameOrientation( ) );
        break;
    }
    default:
        throw std::runtime_error( "Error when getting multi-arc parameter from single-arc equivalent, parameter type " +
                                  boost::lexical_cast< std::string >( singleArcParameter->getParameterName( ).first ) +
                                  " not recognized." );
    }
    return multiArcParameter;
}


//! Function to get initial state vector of estimated dynamical states.
/*!
 *  Function to get initial state vector of estimated dynamical states (i.e. presently estimated state at propagation
 *  start time.
 *  \param estimatableParameters Object containing all parameters that are to be estimated.
 *  \param propagatorSettings Object containing propagation settings to be used
 *  \return State vector of estimated dynamics at propagation start time.
 */
template< typename InitialStateParameterType = double, typename TimeType = double >
void setInitialStateVectorFromParameterSet(
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > > estimatableParameters,
        const std::shared_ptr< propagators::PropagatorSettings< InitialStateParameterType > > propagatorSettings )
{
    typedef Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > VectorType;

    // Retrieve initial dynamical parameters.
    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter<
            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > initialDynamicalParameters =
            estimatableParameters->getEstimatedInitialStateParameters( );

    // Initialize state vector.
    Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > initialStateVector =
            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 >::Zero(
                estimatableParameters->getInitialDynamicalStateParameterSize( ), 1 );

    if( std::dynamic_pointer_cast< propagators::SingleArcPropagatorSettings< InitialStateParameterType, TimeType  > >( propagatorSettings ) != nullptr )
    {
        int vectorSize = 0;
        // Iterate over list of bodies of which the partials of the accelerations acting on them are required.
        for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
        {
            if( isParameterDynamicalPropertyInitialState( initialDynamicalParameters.at( i )->getParameterName( ).first ) )
            {
                int currentParameterSize = initialDynamicalParameters.at( i )->getParameterSize( );
                initialStateVector.block( vectorSize, 0, currentParameterSize, 1 ) = initialDynamicalParameters.at( i )->getParameterValue( );

                vectorSize += currentParameterSize;
            }
        }

        propagatorSettings->resetInitialStates( initialStateVector.block( 0, 0, vectorSize, 1 ) );
    }
    else if( std::dynamic_pointer_cast< propagators::MultiArcPropagatorSettings< InitialStateParameterType, TimeType > >( propagatorSettings ) )
    {
        std::shared_ptr< propagators::MultiArcPropagatorSettings< InitialStateParameterType, TimeType > > multiArcSettings =
                std::dynamic_pointer_cast< propagators::MultiArcPropagatorSettings< InitialStateParameterType, TimeType > >(
                    propagatorSettings );

        std::vector< std::shared_ptr< propagators::SingleArcPropagatorSettings< InitialStateParameterType, TimeType > > >
                singleArcSettings = multiArcSettings->getSingleArcSettings( );
        int numberOfArcs = singleArcSettings.size( );

        // Counting in how many arcs each body has already been propagated/estimated.
        std::map< propagators::IntegratedStateType, std::map< std::string, unsigned int > > initialStatesBodiesCounter;

        for( int i = 0; i < numberOfArcs; i++ )
        {
            std::map< propagators::IntegratedStateType, std::map< std::pair< std::string, std::string >, VectorType > > currentArcInitialStates;

            std::shared_ptr< propagators::TranslationalStatePropagatorSettings< InitialStateParameterType, TimeType > > singleArcTranslationalStatePropagatorSettings
            = std::dynamic_pointer_cast< propagators::TranslationalStatePropagatorSettings< InitialStateParameterType, TimeType > >( singleArcSettings.at( i ) );

            if ( singleArcTranslationalStatePropagatorSettings == nullptr )
            {
                throw std::runtime_error( "Propagator settings for arc " + std::to_string( i ) + " are not translational state propagator settings. "
                                                                            " Other propagator settings not supported (yet) for multi-arc propagation." );
            }
            else
            {

                for ( unsigned int j = 0; j < initialDynamicalParameters.size( ); j++ ) {
                    VectorType currentParameterValue = initialDynamicalParameters.at(j)->getParameterValue();
                    int currentParameterSize = initialDynamicalParameters.at(j)->getParameterSize();
                    std::pair<std::string, std::string> bodyIdentifier = initialDynamicalParameters.at(j)->getParameterName().second;
//                    std::cout << "body identifier: " << bodyIdentifier.first << " & " << bodyIdentifier.second << "\n\n";

                    auto itr = std::find( singleArcTranslationalStatePropagatorSettings->bodiesToIntegrate_.begin( ),
                                                                          singleArcTranslationalStatePropagatorSettings->bodiesToIntegrate_.end( ), bodyIdentifier.first );
                    if ( itr != singleArcTranslationalStatePropagatorSettings->bodiesToIntegrate_.cend( ) )
                    {
//                        std::cout << "body " << bodyIdentifier.first << " - detected for arc " << i << "\n\n";


                        switch ( initialDynamicalParameters.at( j )->getParameterName( ).first )
                        {
                            case estimatable_parameters::arc_wise_initial_body_state:
                            {
//                                std::cout << "current parameter size: " << currentParameterSize << "\n\n";

                                if ( ( initialStatesBodiesCounter.count( propagators::translational_state )  == 0 )
                                || initialStatesBodiesCounter.at( propagators::translational_state ).count( bodyIdentifier.first ) == 0 )
                                {
                                    initialStatesBodiesCounter[ propagators::translational_state ][ bodyIdentifier.first ] = 0;
                                }
                                int index = initialStatesBodiesCounter.at( propagators::translational_state ).at( bodyIdentifier.first );
//                                std::cout << "index: " << index << "\n\n";


//                                if ( currentParameterSize / numberOfArcs != 6 )
//                                {
//                                    throw std::runtime_error("Error when moving initial states from parameters to propagator settings. Incompatible multi-arc translational state size found");
//                                }
//                                std::cout << "full parameters values: " << currentParameterValue.transpose( ) << "\n\n";
//                                std::cout << "subset parameters values: " << currentParameterValue.segment( index * 6, 6 ).transpose( ) << "\n\n";
                                currentArcInitialStates[ propagators::translational_state ][ bodyIdentifier ] = currentParameterValue.segment( index * 6, 6 );

                                // update counter of propagated/estimated bodies
                                initialStatesBodiesCounter.at( propagators::translational_state ).at( bodyIdentifier.first ) ++;

                                break;
                            }
                            default:
                                throw std::runtime_error("Error when moving initial states from parameters to propagator settings. Multi-arc parameter type not recognized");
                        }
                    }
                }
            }
            propagators::resetSingleArcInitialStates( singleArcSettings.at( i ), currentArcInitialStates );
            multiArcSettings->updateInitialStateFromConsituentSettings( );
        }
    }
    else if( std::dynamic_pointer_cast< propagators::HybridArcPropagatorSettings< InitialStateParameterType, TimeType > >( propagatorSettings ) )
    {
        std::shared_ptr< propagators::HybridArcPropagatorSettings< InitialStateParameterType, TimeType > > hybridArcSettings =
                std::dynamic_pointer_cast< propagators::HybridArcPropagatorSettings< InitialStateParameterType, TimeType > >( propagatorSettings );

        std::shared_ptr< propagators::SingleArcPropagatorSettings< InitialStateParameterType, TimeType > > singleArcSettings =
                hybridArcSettings->getSingleArcPropagatorSettings( );
        std::shared_ptr< propagators::MultiArcPropagatorSettings< InitialStateParameterType, TimeType > > multiArcSettings =
                hybridArcSettings->getMultiArcPropagatorSettings( );

        setInitialStateVectorFromParameterSet< InitialStateParameterType, TimeType >(
                    std::make_shared< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > >(
                        std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > >( ),
                        std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > >( ),
                        estimatableParameters->getEstimatedSingleArcInitialStateParameters( ) ), singleArcSettings );
        setInitialStateVectorFromParameterSet< InitialStateParameterType, TimeType >(
                    std::make_shared< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > >(
                        std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > >( ),
                        std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > >( ),
                        estimatableParameters->getEstimatedMultiArcInitialStateParameters( ) ), multiArcSettings );

        hybridArcSettings->setInitialStatesFromConstituents( );

    }
    else
    {
        throw std::runtime_error( "Error when setting initial state vector from parameter set; type not identified" );
    }

}

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_CREATEESTIMATABLEPARAMETERS_H
