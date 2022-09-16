/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATESTATEDERIVATIVEMODEL_H
#define TUDAT_CREATESTATEDERIVATIVEMODEL_H

#include <string>
#include <boost/bind/bind.hpp>
using namespace boost::placeholders;


#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/astro/propagators/singleStateTypeDerivative.h"
#include "tudat/simulation/propagation_setup/propagationSettings.h"
#include "tudat/astro/reference_frames/aerodynamicAngleCalculator.h"
#include "tudat/astro/propagators/nBodyCowellStateDerivative.h"
#include "tudat/astro/propagators/nBodyEnckeStateDerivative.h"
#include "tudat/astro/propagators/nBodyGaussKeplerStateDerivative.h"
#include "tudat/astro/propagators/nBodyGaussModifiedEquinoctialStateDerivative.h"
#include "tudat/astro/propagators/nBodyUnifiedStateModelQuaternionsStateDerivative.h"
#include "tudat/astro/propagators/nBodyUnifiedStateModelModifiedRodriguesParametersStateDerivative.h"
#include "tudat/astro/propagators/nBodyUnifiedStateModelExponentialMapStateDerivative.h"
#include "tudat/astro/propagators/rotationalMotionStateDerivative.h"
#include "tudat/astro/propagators/rotationalMotionQuaternionsStateDerivative.h"
#include "tudat/astro/propagators/rotationalMotionModifiedRodriguesParametersStateDerivative.h"
#include "tudat/astro/propagators/rotationalMotionExponentialMapStateDerivative.h"
#include "tudat/astro/propagators/bodyMassStateDerivative.h"
#include "tudat/astro/propagators/customStateDerivative.h"
#include "tudat/astro/propagators/stateDerivativeCircularRestrictedThreeBodyProblem.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/math/integrators/createNumericalIntegrator.h"

namespace tudat
{


namespace reference_frames
{


class FromBodyAerodynamicAngleInterface: public BodyFixedAerodynamicAngleInterface
{
public:
    FromBodyAerodynamicAngleInterface(
            const std::shared_ptr< simulation_setup::Body > body ):
    BodyFixedAerodynamicAngleInterface( body_fixed_angles_from_body ),
    body_( body ){ }

    virtual ~FromBodyAerodynamicAngleInterface( ){ }

    Eigen::Vector3d getAngles( const double time,
                               const Eigen::Matrix3d& trajectoryToInertialFrame )
    {
        return computeBodyFixedAeroAngles(
                    body_->getCurrentRotationMatrixToLocalFrame( ), trajectoryToInertialFrame );
    }

private:

    std::shared_ptr< simulation_setup::Body > body_;

};

}

namespace propagators
{

//! Function to create object handling frame origin transformations during numerical integration
/*!
 * Function to create object handling frame origin transformations during numerical integration
 * (needed by NBodyStateDerivative).
 * \param centralBodies Names of central bodies, belonging to the
 * entries in the bodiesToIntegrate vector of same index.
 * \param bodiesToIntegrate Names of bodies that are to be integrated numerically.
 * \param bodies List of body objects used in simulation.
 * \return Object handling frame origin transformations during numerical integration
 */
template< typename StateScalarType, typename TimeType >
std::shared_ptr< CentralBodyData< StateScalarType, TimeType > > createCentralBodyData(
        const std::vector< std::string >& centralBodies,
        const std::vector< std::string >& bodiesToIntegrate,
        const simulation_setup::SystemOfBodies& bodies )
{

    // Check whether the bodies that are to e integrated exist in bodies
    for( unsigned int i = 0; i < bodiesToIntegrate.size( ); i++ )
    {
        if( bodies.count( bodiesToIntegrate.at( i ) ) == 0 )
        {
            throw std::runtime_error(
                        "Warning when creating CentralBodyData, body " + bodiesToIntegrate.at( i )
                        + " not present in provided system of bodies." );
        }
    }

    std::vector< std::string > centralBodiesToUse;
    // Generate central body data for requested settings; if no
    // central bodies provided, get inertial SSB as central body
    if( centralBodies.size( ) == 0 )
    {
        for( unsigned int i = 0; i < bodiesToIntegrate.size( ); i++ )
        {
            centralBodiesToUse.push_back( "SSB" );
        }
    }
    else
    {
        centralBodiesToUse = centralBodies;
    }

    std::map< std::string, std::function< Eigen::Matrix< StateScalarType, 6, 1 >
            ( const TimeType ) > > bodyStateFunctions;

    // Retrieve frame origin state functions
    for( unsigned int i = 0; i < centralBodiesToUse.size( ); i++ )
    {
        if( centralBodiesToUse.at( i )  != "SSB" )
        {
            bodyStateFunctions[ centralBodiesToUse.at( i ) ] =
                    std::bind( &simulation_setup::Body::getStateInBaseFrameFromEphemeris
                               < StateScalarType, TimeType >,
                               bodies.at( centralBodiesToUse.at( i ) ), std::placeholders::_1 );
        }
        else
        {
            bodyStateFunctions[ centralBodiesToUse.at( i ) ] = [ ]( const TimeType ){ return
                        Eigen::Matrix< StateScalarType, 6, 1 >::Zero( ); };
        }
    }

    // Get state function of global frame origin w.r.t. barycenter
    std::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > globalFrameOriginBarycentricFunction;
    std::string globalFrameOrigin = simulation_setup::getGlobalFrameOrigin( bodies );
    if( globalFrameOrigin == "SSB" )
    {
        globalFrameOriginBarycentricFunction = [ ]( const TimeType ){ return
                    Eigen::Matrix< StateScalarType, 6, 1 >::Zero( ); };
    }
    else
    {
        globalFrameOriginBarycentricFunction =
                std::bind( &simulation_setup::Body::getGlobalFrameOriginBarycentricStateFromEphemeris< StateScalarType, TimeType >,
                           bodies.at( globalFrameOrigin ), std::placeholders::_1 );
    }

    return std::make_shared< CentralBodyData< StateScalarType, TimeType > >(
                centralBodiesToUse, bodiesToIntegrate, bodyStateFunctions,
                globalFrameOriginBarycentricFunction, globalFrameOrigin );
}

//! Function to create a translational state derivative model.
/*!
 *  Function to create a translational state derivative model from
 *  propagation settings and environment.
 *  \param translationPropagatorSettings Settings for the translational dynamics model.
 *  \param bodies List of body objects in the environment
 *  \param propagationStartTime Time from which numerical propagation starts.
 *  \return Translational state derivative model (instance of derived class of NBodyStateDerivative)
 */
template< typename StateScalarType = double, typename TimeType = double >
std::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > >
createTranslationalStateDerivativeModel(
        const std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > >
        translationPropagatorSettings,
        const simulation_setup::SystemOfBodies& bodies,
        const TimeType propagationStartTime )
{

    // Create object for frame origin transformations.
    std::shared_ptr< CentralBodyData< StateScalarType, TimeType > > centralBodyData =
            createCentralBodyData< StateScalarType, TimeType >(
                translationPropagatorSettings->centralBodies_,
                translationPropagatorSettings->bodiesToIntegrate_,
                bodies );

    std::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > stateDerivativeModel;

    // Check propagator type and create corresponding state derivative object.
    switch( translationPropagatorSettings->propagator_ )
    {
    case cowell:
    {
        stateDerivativeModel = std::make_shared<
                NBodyCowellStateDerivative< StateScalarType, TimeType > >
                ( translationPropagatorSettings->getAccelerationsMap( ), centralBodyData,
                  translationPropagatorSettings->bodiesToIntegrate_ );
        break;
    }
    case encke:
    {
        // Calculate initial Kepler elements for Encke propagator
        std::vector< Eigen::Matrix< StateScalarType, 6, 1 > > initialKeplerElements;
        initialKeplerElements.resize( translationPropagatorSettings->bodiesToIntegrate_.size( ) );
        std::vector< std::string > centralBodies = translationPropagatorSettings->centralBodies_;

        for( unsigned int i = 0; i < translationPropagatorSettings->bodiesToIntegrate_.size( ); i++ )
        {
            if( bodies.count( centralBodies[ i ] ) == 0 )
            {
                std::string errorMessage = "Error when creating Encke propagator, did not find central body " + centralBodies[ i ];
                throw std::runtime_error( errorMessage );
            }
            initialKeplerElements[ i ] = orbital_element_conversions::convertCartesianToKeplerianElements< StateScalarType >(
                        translationPropagatorSettings->getInitialStates( ).segment( i * 6, 6 ), static_cast< StateScalarType >(
                            bodies.at( centralBodies[ i ] )->getGravityFieldModel( )->getGravitationalParameter( ) ) );
        }

        // Create Encke state derivative object.
        stateDerivativeModel = std::make_shared< NBodyEnckeStateDerivative< StateScalarType, TimeType > >
                ( translationPropagatorSettings->getAccelerationsMap( ), centralBodyData, translationPropagatorSettings->bodiesToIntegrate_,
                  initialKeplerElements, propagationStartTime );

        break;
    }
    case gauss_keplerian:
    {
        // Create Keplerian state derivative object.
        stateDerivativeModel = std::make_shared< NBodyGaussKeplerStateDerivative< StateScalarType, TimeType > >
                ( translationPropagatorSettings->getAccelerationsMap( ), centralBodyData,
                  translationPropagatorSettings->bodiesToIntegrate_ );
        break;
    }
    case gauss_modified_equinoctial:
    {
        std::vector< Eigen::Matrix< StateScalarType, 6, 1 > > initialKeplerElements;
        std::vector< std::string > centralBodies = translationPropagatorSettings->centralBodies_;

        for( unsigned int i = 0; i < translationPropagatorSettings->bodiesToIntegrate_.size( ); i++ )
        {
            if( bodies.count( centralBodies[ i ] ) == 0 )
            {
                std::string errorMessage = "Error when creating modified equinoctial propagator, did not find central body " + centralBodies[ i ];
                throw std::runtime_error( errorMessage );
            }
            initialKeplerElements.push_back( orbital_element_conversions::convertCartesianToKeplerianElements< StateScalarType >(
                                                 translationPropagatorSettings->getInitialStates( ).segment( i * 6, 6 ), static_cast< StateScalarType >(
                                                     bodies.at( centralBodies[ i ] )->getGravityFieldModel( )->getGravitationalParameter( ) ) ) );
        }

        // Create modified equinoctial state derivative object.:
        stateDerivativeModel = std::make_shared< NBodyGaussModifiedEquinictialStateDerivative< StateScalarType, TimeType > >
                ( translationPropagatorSettings->getAccelerationsMap( ), centralBodyData,
                  translationPropagatorSettings->bodiesToIntegrate_, initialKeplerElements );

        break;
    }
    case unified_state_model_quaternions:
    {
        // Create USM7 state derivative object.
        stateDerivativeModel = std::make_shared< NBodyUnifiedStateModelQuaternionsStateDerivative< StateScalarType, TimeType > >
                ( translationPropagatorSettings->getAccelerationsMap( ), centralBodyData,
                  translationPropagatorSettings->bodiesToIntegrate_ );
        break;
    }
    case unified_state_model_modified_rodrigues_parameters:
    {
        // Create USM6 state derivative object.
        stateDerivativeModel = std::make_shared< NBodyUnifiedStateModelModifiedRodriguesParametersStateDerivative< StateScalarType, TimeType > >
                ( translationPropagatorSettings->getAccelerationsMap( ), centralBodyData,
                  translationPropagatorSettings->bodiesToIntegrate_ );
        break;
    }
    case unified_state_model_exponential_map:
    {
        // Create USMEM state derivative object.
        stateDerivativeModel = std::make_shared< NBodyUnifiedStateModelExponentialMapStateDerivative< StateScalarType, TimeType > >
                ( translationPropagatorSettings->getAccelerationsMap( ), centralBodyData,
                  translationPropagatorSettings->bodiesToIntegrate_ );
        break;
    }
    default:
        throw std::runtime_error( "Error, did not recognize translational state propagation type: " +
                                  std::to_string( translationPropagatorSettings->propagator_ ) );
    }
    return stateDerivativeModel;
}

//! Function to create a rotational dynamics state derivative model.
/*!
 *  Function to create a rotational dynamics state derivative model from propagation settings and environment.
 *  \param rotationPropagatorSettings Settings for the rotational dynamics model.
 *  \param bodies List of body objects in the environment
 *  \param startTime propagation start time
 *  \return Rotational dynamics state derivative model.
 */
template< typename StateScalarType = double, typename TimeType = double >
std::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > createRotationalStateDerivativeModel(
        const std::shared_ptr< RotationalStatePropagatorSettings< StateScalarType > > rotationPropagatorSettings,
        const simulation_setup::SystemOfBodies& bodies, const TimeType startTime )
{
    std::vector< std::function< Eigen::Matrix3d( ) > > momentOfInertiaFunctions;
    for( unsigned int i = 0; i < rotationPropagatorSettings->bodiesToIntegrate_.size( ); i++ )
    {
        momentOfInertiaFunctions.push_back(
                    std::bind( &simulation_setup::Body::getBodyInertiaTensor,
                               bodies.at( rotationPropagatorSettings->bodiesToIntegrate_.at( i ) ) ) );
    }

    // Check propagator type and create corresponding state derivative object.
    std::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > stateDerivativeModel;
    switch( rotationPropagatorSettings->propagator_ )
    {
    case quaternions:
    {
        stateDerivativeModel = std::make_shared< RotationalMotionQuaternionsStateDerivative< StateScalarType, TimeType > >(
                    rotationPropagatorSettings->getTorqueModelsMap( ), rotationPropagatorSettings->bodiesToIntegrate_,
                    momentOfInertiaFunctions );
        break;
    }
    case modified_rodrigues_parameters:
    {
        stateDerivativeModel = std::make_shared< RotationalMotionModifiedRodriguesParametersStateDerivative< StateScalarType, TimeType > >(
                    rotationPropagatorSettings->getTorqueModelsMap( ), rotationPropagatorSettings->bodiesToIntegrate_,
                    momentOfInertiaFunctions );
        break;
    }
    case exponential_map:
    {
        stateDerivativeModel = std::make_shared< RotationalMotionExponentialMapStateDerivative< StateScalarType, TimeType > >(
                    rotationPropagatorSettings->getTorqueModelsMap( ), rotationPropagatorSettings->bodiesToIntegrate_,
                    momentOfInertiaFunctions );
        break;
    }
    default:
        throw std::runtime_error( "Error, did not recognize rotational state propagation type: " +
                                  std::to_string( rotationPropagatorSettings->propagator_ ) );
    }

    return stateDerivativeModel;
}

//! Function to create a mass state derivative model.
/*!
 *  Function to create a mass state derivative model from propagation settings and environment.
 *  \param massPropagatorSettings Settings for the mass dynamics model.
 *  \param bodies List of body objects in the environment
 *  \return Mass state derivative model.
 */
template< typename StateScalarType = double, typename TimeType = double >
std::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > createBodyMassStateDerivativeModel(
        const std::shared_ptr< MassPropagatorSettings< StateScalarType > > massPropagatorSettings,
        const  simulation_setup::SystemOfBodies& bodies )
{
    return std::make_shared< propagators::BodyMassStateDerivative< StateScalarType, TimeType > >(
                massPropagatorSettings->getMassRateModelsMap( ),
                massPropagatorSettings->bodiesWithMassToPropagate_ );
}

//! Function to create a state derivative model.
/*!
 *  Function to create a state derivative model from propagation settings and the environment.
 *  \param propagatorSettings Settings for the dynamical model.
 *  \param bodies List of body objects in the environment
 *  \param propagationStartTime Time from which numerical propagation starts.
 *  \return State derivative model (instance of required derived class of SingleStateTypeDerivative)
 */
template< typename StateScalarType = double, typename TimeType = double >
std::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > >
createStateDerivativeModel(
        const std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > propagatorSettings,
        const simulation_setup::SystemOfBodies& bodies,
        const TimeType propagationStartTime )
{
    std::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > stateDerivativeModel;

    // Check dynamics type and call associated function to create
    // specific type of state derivative model.
    switch( propagatorSettings->getStateType( ) )
    {
    case translational_state:
    {
        // Check input consistency.
        std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > >
                translationPropagatorSettings =
                std::dynamic_pointer_cast<
                TranslationalStatePropagatorSettings< StateScalarType > >( propagatorSettings );
        if( translationPropagatorSettings == nullptr )
        {
            throw std::runtime_error(
                        "Error, expected translational state propagation settings when making state derivative model" );
        }
        else
        {
            stateDerivativeModel = createTranslationalStateDerivativeModel< StateScalarType, TimeType >(
                        translationPropagatorSettings, bodies, propagationStartTime );
        }
        break;
    }
    case rotational_state:
    {
        std::shared_ptr< RotationalStatePropagatorSettings< StateScalarType > > rotationPropagatorSettings =
                std::dynamic_pointer_cast< RotationalStatePropagatorSettings< StateScalarType > >( propagatorSettings );
        if( rotationPropagatorSettings == nullptr )
        {
            throw std::runtime_error(
                        "Error, expected rotation state propagation settings when making state derivative model" );
        }
        else
        {
            stateDerivativeModel = createRotationalStateDerivativeModel< StateScalarType, TimeType >(
                        rotationPropagatorSettings, bodies, propagationStartTime );
        }
        break;
    }
    case body_mass_state:
    {
        // Check input consistency.
        std::shared_ptr< MassPropagatorSettings< StateScalarType > > massPropagatorSettings =
                std::dynamic_pointer_cast< MassPropagatorSettings< StateScalarType > >( propagatorSettings );
        if( massPropagatorSettings == nullptr )
        {
            throw std::runtime_error(
                        "Error, expected mass propagation settings when making state derivative model" );
        }
        else
        {
            stateDerivativeModel = createBodyMassStateDerivativeModel< StateScalarType, TimeType >(
                        massPropagatorSettings, bodies );
        }
        break;
    }
    case custom_state:
    {
        // Check input consistency.
        std::shared_ptr< CustomStatePropagatorSettings< StateScalarType, TimeType > > customPropagatorSettings =
                std::dynamic_pointer_cast< CustomStatePropagatorSettings< StateScalarType, TimeType > >(
                    propagatorSettings );
        if( customPropagatorSettings == nullptr )
        {
            throw std::runtime_error(
                        "Error, expected custom propagation settings when making state derivative model" );
        }
        else
        {
            stateDerivativeModel = std::make_shared< CustomStateDerivative< StateScalarType, TimeType > >(
                        customPropagatorSettings->stateDerivativeFunction_, customPropagatorSettings->stateSize_ );
        }
        break;
    }
    default:
        throw std::runtime_error(
                    "Error, could not process state type "
                    + std::to_string( propagatorSettings->getStateType( ) )
                    + " when making state derivative model" );
    }
    return stateDerivativeModel;
}

//! Function that finalized multi-type propagator creation by ensuring that any mutual dependencies are correctly set
/*!
 *  Function that finalized multi-type propagator creation by ensuring that any mutual dependencies are correctly set
 *  \param propagatorSettings Settings for the numerical propagation
 *  \param bodies List of body objects that comprises the environment
 */
template< typename StateScalarType = double >
void setMultiTypePropagationClosure(
        const std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > propagatorSettings,
        const simulation_setup::SystemOfBodies& bodies )
{
    // Cast to multi-type settings, and perform closure if
    std::shared_ptr< MultiTypePropagatorSettings< StateScalarType > > multiTypePropagatorSettings =
            std::dynamic_pointer_cast< MultiTypePropagatorSettings< StateScalarType > >( propagatorSettings );

    if( multiTypePropagatorSettings != nullptr )
    {
        // Perform closure for the case where both translational and rotational states are propagated
        if( multiTypePropagatorSettings->propagatorSettingsMap_.count( translational_state ) > 0 &&
                multiTypePropagatorSettings->propagatorSettingsMap_.count( rotational_state ) > 0 )
        {
            std::vector< std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > > translationalStateSettings =
                    multiTypePropagatorSettings->propagatorSettingsMap_.at( translational_state );
            std::vector< std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > > rotationalStateSettings =
                    multiTypePropagatorSettings->propagatorSettingsMap_.at( rotational_state );

            // Iterate over all accelerations, and identify those bodies for which an aerodynamic acceleration is applied
            std::vector< std::string > bodiesWithAerodynamicAcceleration;
            for( unsigned int i = 0; i < translationalStateSettings.size( ); i++ )
            {
                std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > > currentTranslationalState =
                        std::dynamic_pointer_cast< TranslationalStatePropagatorSettings< StateScalarType > >(
                            translationalStateSettings.at( i ) );
                basic_astrodynamics::AccelerationMap currentAccelerationsMap = currentTranslationalState->getAccelerationsMap( );

                for( basic_astrodynamics::AccelerationMap::const_iterator accelerationIterator = currentAccelerationsMap.begin( );
                     accelerationIterator != currentAccelerationsMap.end( ); accelerationIterator++ )
                {
                    for( basic_astrodynamics::SingleBodyAccelerationMap::const_iterator singleBodyIterator =
                         accelerationIterator->second.begin( ); singleBodyIterator != accelerationIterator->second.end( );
                         singleBodyIterator++ )
                    {
                        for( unsigned int j = 0; j < singleBodyIterator->second.size( ); j++ )
                        {
                            if( std::find( bodiesWithAerodynamicAcceleration.begin( ), bodiesWithAerodynamicAcceleration.end( ),
                                           accelerationIterator->first ) == bodiesWithAerodynamicAcceleration.end( ) )
                            {
                                if( basic_astrodynamics::getAccelerationModelType( singleBodyIterator->second.at( j ) ) ==
                                        basic_astrodynamics::aerodynamic )
                                {
                                    bodiesWithAerodynamicAcceleration.push_back( accelerationIterator->first );
                                }
                            }
                        }
                    }
                }
            }

            // Iterate over all settings and identify those bodies for which rotational dynamics is propagated.
            std::vector< std::string > bodiesWithPropagatedRotation;
            for( unsigned int i = 0; i < rotationalStateSettings.size( ); i++ )
            {
                std::shared_ptr< RotationalStatePropagatorSettings< StateScalarType > > currentTranslationalState =
                        std::dynamic_pointer_cast< RotationalStatePropagatorSettings< StateScalarType > >(
                            rotationalStateSettings.at( i ) );
                std::vector< std::string > currentBodiesWithPropagatedRotation = currentTranslationalState->bodiesToIntegrate_;
                for( unsigned int j = 0; j < currentBodiesWithPropagatedRotation.size( ); j++ )
                {
                    if( std::find( bodiesWithPropagatedRotation.begin( ), bodiesWithPropagatedRotation.end( ),
                                   currentBodiesWithPropagatedRotation.at( j ) ) == bodiesWithPropagatedRotation.end( ) )
                    {
                        bodiesWithPropagatedRotation.push_back( currentBodiesWithPropagatedRotation.at( j ) );
                    }
                }
            }

            // Find bodies for which both aerodynamic acceleration is used and rotational propagation is performed.
            std::vector< std::string > bodiesWithAerodynamicRotationalClosure;
            for( unsigned int i = 0; i < bodiesWithPropagatedRotation.size( ); i++ )
            {
                if( std::find( bodiesWithAerodynamicAcceleration.begin( ), bodiesWithAerodynamicAcceleration.end( ),
                               bodiesWithPropagatedRotation.at( i ) ) != bodiesWithAerodynamicAcceleration.end( ) )
                {
                    bodiesWithAerodynamicRotationalClosure.push_back( bodiesWithPropagatedRotation.at( i ) );
                }
            }

            // Ensure that vehicle orientation is correctly set for aerodynamic acceleration/torque
            for( unsigned int i = 0; i < bodiesWithAerodynamicRotationalClosure.size( ); i++ )
            {
                std::shared_ptr< aerodynamics::FlightConditions > currentFlightConditions =
                        bodies.at( bodiesWithAerodynamicRotationalClosure.at( i ) )->getFlightConditions( );
                std::shared_ptr< reference_frames::FromBodyAerodynamicAngleInterface > rotationInterface =
                        std::make_shared< reference_frames::FromBodyAerodynamicAngleInterface >(
                            bodies.at( bodiesWithAerodynamicRotationalClosure.at( i ) ) );
                currentFlightConditions->getAerodynamicAngleCalculator( )->setBodyFixedAngleInterface( rotationInterface );
            }
        }
    }
}


//! Function to create a list of state derivative models.
/*!
 *  Function to create a list of state derivative models from
 *  propagation settings and the environment.
 *  \param propagatorSettings Settings for the dynamical model.
 *  \param bodies List of body objects in the environment.
 *  \param propagationStartTime Time from which numerical propagation starts.
 *  \return List of state derivative models (instances of required
 *  derived class of SingleStateTypeDerivative)
 */
template< typename StateScalarType = double, typename TimeType = double >
std::vector< std::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > >
createStateDerivativeModels(
        const std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > propagatorSettings,
        const simulation_setup::SystemOfBodies& bodies,
        const TimeType propagationStartTime )
{
    std::vector< std::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > >
            stateDerivativeModels;

    // Check type of state derivative model and call associated create function.
    switch( propagatorSettings->getStateType( ) )
    {
    // If hybrid, call create function separately for each entry.
    case hybrid:
    {
        std::shared_ptr< MultiTypePropagatorSettings< StateScalarType > > multiTypePropagatorSettings =
                std::dynamic_pointer_cast< MultiTypePropagatorSettings< StateScalarType > >( propagatorSettings );

        // Iterate over all propagation settings
        for( typename std::map< IntegratedStateType,
             std::vector< std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > > >::iterator
             propagatorIterator = multiTypePropagatorSettings->propagatorSettingsMap_.begin( );
             propagatorIterator != multiTypePropagatorSettings->propagatorSettingsMap_.end( ); propagatorIterator++ )
        {
            for( unsigned int i = 0; i < propagatorIterator->second.size( ); i++ )
            {
                // Call create function for current propagation settings.
                if( propagatorIterator->first != hybrid )
                {
                    stateDerivativeModels.push_back( createStateDerivativeModel< StateScalarType, TimeType >(
                                                         propagatorIterator->second.at( i ), bodies, propagationStartTime ) );
                }
                else
                {
                    throw std::runtime_error(
                                "Error when making state derivative model, cannot process nested hybrid propagators" );
                }
            }
        }

        setMultiTypePropagationClosure( propagatorSettings, bodies );

        break;
    }
        // If not hybrid, call create function for single object directly.
    default:
        stateDerivativeModels.push_back( createStateDerivativeModel< StateScalarType, TimeType >(
                                             propagatorSettings, bodies, propagationStartTime ) );
    }

    return stateDerivativeModels;
}

//extern template std::vector< std::shared_ptr< SingleStateTypeDerivative< double, double > > > createStateDerivativeModels< double, double >(
//        const std::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings,
//        const simulation_setup::SystemOfBodies& bodies,
//        const double propagationStartTime );
//extern template std::shared_ptr< SingleStateTypeDerivative< double, double > > createStateDerivativeModel< double, double >(
//        const std::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings,
//        const simulation_setup::SystemOfBodies& bodies,
//        const double propagationStartTime );

//#if( TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS )
//extern template std::vector< std::shared_ptr< SingleStateTypeDerivative< long double, double > > > createStateDerivativeModels< long double, double >(
//        const std::shared_ptr< SingleArcPropagatorSettings< long double > > propagatorSettings,
//        const simulation_setup::SystemOfBodies& bodies,
//        const double propagationStartTime );
//extern template std::vector< std::shared_ptr< SingleStateTypeDerivative< double, Time > > > createStateDerivativeModels< double, Time >(
//        const std::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings,
//        const simulation_setup::SystemOfBodies& bodies,
//        const Time propagationStartTime );
//extern template std::vector< std::shared_ptr< SingleStateTypeDerivative< long double, Time > > > createStateDerivativeModels< long double, Time >(
//        const std::shared_ptr< SingleArcPropagatorSettings< long double > > propagatorSettings,
//        const simulation_setup::SystemOfBodies& bodies,
//        const Time propagationStartTime );
//extern template std::shared_ptr< SingleStateTypeDerivative< long double, double > > createStateDerivativeModel< long double, double >(
//        const std::shared_ptr< SingleArcPropagatorSettings< long double > > propagatorSettings,
//        const simulation_setup::SystemOfBodies& bodies,
//        const double propagationStartTime );
//extern template std::shared_ptr< SingleStateTypeDerivative< double, Time > > createStateDerivativeModel< double, Time >(
//        const std::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings,
//        const simulation_setup::SystemOfBodies& bodies,
//        const Time propagationStartTime );
//extern template std::shared_ptr< SingleStateTypeDerivative< long double, Time > > createStateDerivativeModel< long double, Time >(
//        const std::shared_ptr< SingleArcPropagatorSettings< long double > > propagatorSettings,
//        const simulation_setup::SystemOfBodies& bodies,
//        const Time propagationStartTime );
//#endif

//! Function to convert a list of state derivative models to a map sorted by state type
/*!
 *  Function to convert a list of state derivative models to a map sorted by state type
 *  \param stateDerivativeModelList List of state derivative models
 *  \return Map of state derivative models
 */
template< typename StateScalarType = double, typename TimeType = double >
std::unordered_map< IntegratedStateType, std::vector< std::shared_ptr<
SingleStateTypeDerivative< StateScalarType, TimeType > > > > getStateDerivativeModelMapFromVector(
        const std::vector< std::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > >& stateDerivativeModelList )
{
    std::unordered_map< IntegratedStateType, std::vector< std::shared_ptr<
            SingleStateTypeDerivative< StateScalarType, TimeType > > > > stateDerivativeModelsMap;
    for( unsigned int i = 0; i < stateDerivativeModelList.size( ); i++ )
    {
        stateDerivativeModelsMap[ stateDerivativeModelList.at( i )->getIntegratedStateType( ) ].push_back(
                    stateDerivativeModelList.at( i ) );
    }
    return stateDerivativeModelsMap;
}

//! Function to create a map of state derivative models.
/*!
 *  Function to create a map of state derivative models from
 *  propagation settings and the environment.
 *  \param propagatorSettings Settings for the dynamical model.
 *  \param bodies List of body objects in the environment.
 *  \param propagationStartTime Time from which numerical propagation starts.
 *  \return Map of state derivative models (instances of required
 *  derived class of SingleStateTypeDerivative)
 */
template< typename StateScalarType = double, typename TimeType = double >
std::unordered_map< IntegratedStateType, std::vector< std::shared_ptr<
SingleStateTypeDerivative< StateScalarType, TimeType > > > >
createStateDerivativeModelMap(
        const std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > propagatorSettings,
        const simulation_setup::SystemOfBodies& bodies,
        const TimeType propagationStartTime )
{
    return getStateDerivativeModelMapFromVector( createStateDerivativeModels(
                                                     propagatorSettings, bodies, propagationStartTime ) );
}

//! Function to create an integrator to propagate the dynamics (in normalized units) in CR3BP
/*!
 *  Function to create an integrator to propagate the dynamics (in normalized units) in Circularly Restricted Three-Body Problem.
 * \param integratorSettings Settings for the numerical integration
 * \param massParameter Normalized mass parameter
 * \param initialState Initial normalized state
 * \return Integrator object for propagation of CR3BP with requested settings
 */
std::shared_ptr< numerical_integrators::NumericalIntegrator< double, Eigen::Vector6d > > createCR3BPIntegrator(
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const double massParameter,
        const Eigen::Vector6d& initialState );

//! Function to propagate the dynamics (in normalized units) in CR3BP
/*!
 *  Function to propagate the dynamics (in normalized units) in Circularly Restricted Three-Body Problem.
 * \param integratorSettings Settings for the numerical integration
 * \param massParameter Normalized mass parameter
 * \param initialState Initial normalized state
 * \param finalTime End time for the numerical integration
 * \param propagateToExactFinalTime Boolean denoting whether to terminate exactly on the final time (if true), or on the first
 * step that exceeds the final time
 * \return Propagated state history of normalized dynamics in CR3BP.
 */
std::map< double, Eigen::Vector6d > performCR3BPIntegration(
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const double massParameter,
        const Eigen::Vector6d& initialState,
        const double finalTime,
        const bool propagateToExactFinalTime = false );


} // namespace propagators

} // namespace tudat

#endif // TUDAT_CREATESTATEDERIVATIVEMODEL_H
