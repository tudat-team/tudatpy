/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/Gravitation/basicSolidBodyTideGravityFieldVariations.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/constantDragCoefficient.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/constantRotationRate.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/constantRotationalOrientation.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/empiricalAccelerationCoefficients.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/gravitationalParameter.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/observationBiasParameter.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/groundStationPosition.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/sphericalHarmonicCosineCoefficients.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/sphericalHarmonicSineCoefficients.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/radiationPressureCoefficient.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/ppnParameters.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/equivalencePrincipleViolationParameter.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/tidalLoveNumber.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/directTidalTimeLag.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/meanMomentOfInertiaParameter.h"
#include "Tudat/Astrodynamics/Relativity/metric.h"
#include "Tudat/SimulationSetup/EstimationSetup/createEstimatableParameters.h"

namespace tudat
{

namespace simulation_setup
{


using namespace simulation_setup;
using namespace ephemerides;
using namespace gravitation;
using namespace estimatable_parameters;

//! Function to create an interface object for estimating a parameter defined by a single double value
std::shared_ptr< EstimatableParameter< double > > createDoubleParameterToEstimate(
        const std::shared_ptr< EstimatableParameterSettings >& doubleParameterName,
        const NamedBodyMap& bodyMap, const basic_astrodynamics::AccelerationMap& accelerationModelMap )
{
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

        if( ( currentBodyName != "global_metric" ) && ( currentBodyName != "" ) && ( bodyMap.count( currentBodyName ) == 0 ) )
        {
            std::string errorMessage = "Error when creating parameters to estimate, body " +
                    currentBodyName + "  not in body map " +
                    std::to_string( doubleParameterName->parameterType_.first );
            throw std::runtime_error( errorMessage );
        }
        else if( ( currentBodyName != "" ) && ( currentBodyName != "global_metric" ) )
        {
            currentBody = bodyMap.at( currentBodyName );
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
            if( currentBody->getRadiationPressureInterfaces( ).size( ) == 0 )
            {
                std::string errorMessage = "Error, no radiation pressure interfaces found in body " +
                        currentBodyName + " when making Cr parameter.";
                throw std::runtime_error( errorMessage );
            }
            else if( currentBody->getRadiationPressureInterfaces( ).size( ) > 1 )
            {
                std::string errorMessage = "Error, multiple radiation pressure interfaces found in body " +
                        currentBodyName + " when making Cr parameter.";
                throw std::runtime_error( errorMessage );
            }
            else
            {
                doubleParameterToEstimate = std::make_shared< RadiationPressureCoefficient >(
                            currentBody->getRadiationPressureInterfaces( ).begin( )->second,
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
            // Check input consistency
            std::shared_ptr< DirectTidalTimeLagEstimatableParameterSettings > dissipationTimeLagSettings =
                    std::dynamic_pointer_cast< DirectTidalTimeLagEstimatableParameterSettings >( doubleParameterName );
            if( dissipationTimeLagSettings == nullptr )
            {
                throw std::runtime_error( "Error, expected dissipation time lag parameter settings." );
            }
            else
            {
                std::vector< std::shared_ptr< DirectTidalDissipationAcceleration > > assiciatedTidalAccelerationModels =
                        getTidalDissipationAccelerationModels(
                            accelerationModelMap, currentBodyName, dissipationTimeLagSettings->deformingBodies_ );

                // Create parameter object
                if( assiciatedTidalAccelerationModels.size( ) != 0 )
                {
                    doubleParameterToEstimate = std::make_shared< DirectTidalTimeLag >(
                                assiciatedTidalAccelerationModels, currentBodyName, dissipationTimeLagSettings->deformingBodies_ );
                }
                else
                {
                    throw std::runtime_error(
                                "Error, expected DirectTidalDissipationAcceleration list for tidal time lag" );
                }
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
            else
            {
                doubleParameterToEstimate = std::make_shared< MeanMomentOfInertiaParameter >
                        ( std::bind( &simulation_setup::Body::getScaledMeanMomentOfInertia, currentBody ),
                          std::bind( &simulation_setup::Body::setScaledMeanMomentOfInertia, currentBody, std::placeholders::_1 ),
                          currentBodyName );
            }
            break;
        }
        default:
            throw std::runtime_error( "Warning, this double parameter has not yet been implemented when making parameters" );
            break;
        }
    }

    return doubleParameterToEstimate;
}

//! Function to create an interface object for estimating a parameter defined by a list of double values
std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > createVectorParameterToEstimate(
        const std::shared_ptr< EstimatableParameterSettings >& vectorParameterName,
        const NamedBodyMap& bodyMap, const basic_astrodynamics::AccelerationMap& accelerationModelMap )
{
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
        if( ( currentBodyName != "" ) && ( bodyMap.count( currentBodyName ) == 0 ) )
        {
            std::string errorMessage = "Warning when creating parameters to estimate, body " +
                    currentBodyName  + "not in body map " + std::to_string( vectorParameterName->parameterType_.first );
            throw std::runtime_error( errorMessage );
        }
        else if( ( currentBodyName != "" ) )
        {
            currentBody = bodyMap.at( currentBodyName );
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
                            biasSettings->linkEnds_, biasSettings->observableType_, true );
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
                            biasSettings->linkEnds_, biasSettings->observableType_, false );
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
                            biasSettings->linkEnds_, biasSettings->observableType_, true );
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
                            biasSettings->linkEnds_, biasSettings->observableType_, false );
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
        case empirical_acceleration_coefficients:
        {
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
                // Check if required acceleration model exists
                if( accelerationModelMap.count( empiricalAccelerationSettings->parameterType_.second.first ) == 0 )
                {
                    std::string errorMessage =
                            "Error, did not find accelerations on body " + 
                            std::string( empiricalAccelerationSettings->parameterType_.second.first ) +
                            " when making constant empirical acceleration coefficients parameter";
                    throw std::runtime_error( errorMessage );
                }
                else if( accelerationModelMap.at( empiricalAccelerationSettings->parameterType_.second.first ).count(
                             empiricalAccelerationSettings->centralBody_ ) == 0 )
                {
                    std::string errorMessage =
                            "Error, did not find accelerations on body " +
                            empiricalAccelerationSettings->parameterType_.second.first +
                            " due to body " + empiricalAccelerationSettings->centralBody_ +
                            " when making constant empirical acceleration coefficients parameter";
                    throw std::runtime_error( errorMessage );
                }
                else
                {
                    // Retrieve acceleration model.
                    std::shared_ptr< basic_astrodynamics::EmpiricalAcceleration > empiricalAcceleration;
                    std::vector< std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > >
                            accelerationModelList =
                            accelerationModelMap.at( empiricalAccelerationSettings->parameterType_.second.first ).at(
                                empiricalAccelerationSettings->centralBody_ );
                    for( unsigned int i = 0; i < accelerationModelList.size( ); i++ )
                    {
                        if( basic_astrodynamics::getAccelerationModelType( accelerationModelList[ i ] ) ==
                                basic_astrodynamics::empirical_acceleration )
                        {
                            empiricalAcceleration = std::dynamic_pointer_cast< basic_astrodynamics::EmpiricalAcceleration >(
                                        accelerationModelList[ i ] );
                        }
                    }

                    if( empiricalAcceleration == nullptr )
                    {
                        throw std::runtime_error(
                                    "Error when making constant empirical acceleration coefficients parameter, could not find acceleration model" );
                    }
                    else
                    {
                        // Create empirical acceleration parameter
                        vectorParameterToEstimate = std::make_shared< EmpiricalAccelerationCoefficientsParameter >(
                                    empiricalAcceleration, empiricalAccelerationSettings->parameterType_.second.first,
                                    empiricalAccelerationSettings->componentsToEstimate_ );
                    }
                }
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
                if( currentBody->getRadiationPressureInterfaces( ).size( ) == 0 )
                {
                    std::string errorMessage = "Error, no radiation pressure interfaces found in body " +
                            currentBodyName + " when making Cr parameter.";
                    throw std::runtime_error( errorMessage );
                }
                else if( currentBody->getRadiationPressureInterfaces( ).size( ) > 1 )
                {
                    std::string errorMessage = "Error, multiple radiation pressure interfaces found in body " +
                            currentBodyName + " when making Cr parameter.";
                    throw std::runtime_error( errorMessage );
                }
                else
                {
                    vectorParameterToEstimate = std::make_shared< ArcWiseRadiationPressureCoefficient >(
                                currentBody->getRadiationPressureInterfaces( ).begin( )->second,
                                radiationPressureCoefficientSettings->arcStartTimeList_,
                                currentBodyName );
                }
                break;
            }
            break;
        }
        case arc_wise_empirical_acceleration_coefficients:
        {
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
                // Check if required acceleration model exists
                if( accelerationModelMap.count( empiricalAccelerationSettings->parameterType_.second.first ) == 0 )
                {
                    std::string errorMessage =
                            "Error, did not find accelerations on body " +  empiricalAccelerationSettings->parameterType_.second.first +
                            " when making arcwise empirical acceleration coefficients parameter";
                    throw std::runtime_error( errorMessage );

                }
                else if( accelerationModelMap.at( empiricalAccelerationSettings->parameterType_.second.first ).count(
                             empiricalAccelerationSettings->centralBody_ ) == 0 )
                {
                    std::string errorMessage =
                            "Error, did not find accelerations on body " +
                            empiricalAccelerationSettings->parameterType_.second.first +
                            " due to body " +  empiricalAccelerationSettings->centralBody_ +
                            " when making arcwise empirical acceleration coefficients parameter";
                    throw std::runtime_error( errorMessage );
                }
                else
                {
                    // Retrieve acceleration model.
                    std::shared_ptr< basic_astrodynamics::EmpiricalAcceleration > empiricalAcceleration;
                    std::vector< std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > > accelerationModelList =
                            accelerationModelMap.at( empiricalAccelerationSettings->parameterType_.second.first ).at(
                                empiricalAccelerationSettings->centralBody_ );
                    for( unsigned int i = 0; i < accelerationModelList.size( ); i++ )
                    {
                        if( basic_astrodynamics::getAccelerationModelType( accelerationModelList[ i ] ) ==
                                basic_astrodynamics::empirical_acceleration )
                        {
                            empiricalAcceleration = std::dynamic_pointer_cast< basic_astrodynamics::EmpiricalAcceleration >(
                                        accelerationModelList[ i ] );
                        }
                    }

                    if( empiricalAcceleration == nullptr )
                    {
                        throw std::runtime_error(
                                    "Error when making constant empirical acceleration coefficients parameter, could not find acceleration model" );
                    }
                    else
                    {
                        // Create arcwise empirical acceleration parameter
                        vectorParameterToEstimate = std::make_shared< ArcWiseEmpiricalAccelerationCoefficientsParameter >(
                                    empiricalAcceleration, empiricalAccelerationSettings->parameterType_.second.first,
                                    empiricalAccelerationSettings->componentsToEstimate_, empiricalAccelerationSettings->arcStartTimeList_ );
                    }
                }
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
        default:
            std::string errorMessage = "Warning, this vector parameter (" +
                    std::to_string( vectorParameterName->parameterType_.first ) +
                    ") has not yet been implemented when making parameters";
            throw std::runtime_error( errorMessage );

            break;
        }
    }

    return vectorParameterToEstimate;
}

}

}
