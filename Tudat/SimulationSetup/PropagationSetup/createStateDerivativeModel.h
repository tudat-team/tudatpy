/*    Copyright (c) 2010-2018, Delft University of Technology
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
#include <boost/bind.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Astrodynamics/Propagators/singleStateTypeDerivative.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationSettings.h"
#include "Tudat/Astrodynamics/ReferenceFrames/aerodynamicAngleCalculator.h"
#include "Tudat/Astrodynamics/Propagators/nBodyCowellStateDerivative.h"
#include "Tudat/Astrodynamics/Propagators/nBodyEnckeStateDerivative.h"
#include "Tudat/Astrodynamics/Propagators/nBodyGaussKeplerStateDerivative.h"
#include "Tudat/Astrodynamics/Propagators/nBodyGaussModifiedEquinoctialStateDerivative.h"
#include "Tudat/Astrodynamics/Propagators/rotationalMotionStateDerivative.h"
#include "Tudat/Astrodynamics/Propagators/bodyMassStateDerivative.h"
#include "Tudat/Astrodynamics/Propagators/customStateDerivative.h"
#include "Tudat/Astrodynamics/Propagators/stateDerivativeCircularRestrictedThreeBodyProblem.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h"

namespace tudat
{

namespace propagators
{

//! Function to create object handling frame origin transformations during numerical integration
/*!
 * Function to create object handling frame origin transformations during numerical integration
 * (needed by NBodyStateDerivative).
 * \param centralBodies Names of central bodies, belonging to the
 * entries in the bodiesToIntegrate vector of same index.
 * \param bodiesToIntegrate Names of bodies that are to be integrated numerically.
 * \param bodyMap List of body objects used in simulation.
 * \return Object handling frame origin transformations during numerical integration
 */
template< typename StateScalarType, typename TimeType >
boost::shared_ptr< CentralBodyData< StateScalarType, TimeType > > createCentralBodyData(
        const std::vector< std::string >& centralBodies,
        const std::vector< std::string >& bodiesToIntegrate,
        const simulation_setup::NamedBodyMap& bodyMap )
{

    // Check whether the bodies that are to e integrated exist in bodyMap
    for( unsigned int i = 0; i < bodiesToIntegrate.size( ); i++ )
    {
        if( bodyMap.count( bodiesToIntegrate.at( i ) ) == 0 )
        {
            throw std::runtime_error(
                        "Warning when creating CentralBodyData, body " + bodiesToIntegrate.at( i )
                        + " not present in provided body map." );
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

    std::map< std::string, boost::function< Eigen::Matrix< StateScalarType, 6, 1 >
            ( const TimeType ) > > bodyStateFunctions;

    // Retrieve frame origin state functions
    for( unsigned int i = 0; i < centralBodiesToUse.size( ); i++ )
    {
        if( centralBodiesToUse.at( i )  != "SSB" )
        {
            bodyStateFunctions[ centralBodiesToUse.at( i ) ] =
                    boost::bind( &simulation_setup::Body::getStateInBaseFrameFromEphemeris
                                 < StateScalarType, TimeType >,
                                 bodyMap.at( centralBodiesToUse.at( i ) ), _1 );
        }
        else
        {
            bodyStateFunctions[ centralBodiesToUse.at( i ) ] = boost::lambda::constant(
                        Eigen::Matrix< StateScalarType, 6, 1 >::Zero( ) );
        }
    }

    // Get state function of global frame origin w.r.t. barycenter
    boost::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > globalFrameOriginBarycentricFunction;
    std::string globalFrameOrigin = simulation_setup::getGlobalFrameOrigin( bodyMap );
    if( globalFrameOrigin == "SSB" )
    {
        globalFrameOriginBarycentricFunction = boost::lambda::constant(
                    Eigen::Matrix< StateScalarType, 6, 1 >::Zero( ) );
    }
    else
    {
        globalFrameOriginBarycentricFunction =
                boost::bind( &simulation_setup::Body::getGlobalFrameOriginBarycentricStateFromEphemeris< StateScalarType, TimeType >,
                             bodyMap.at( globalFrameOrigin ), _1 );
    }

    return boost::make_shared< CentralBodyData< StateScalarType, TimeType > >(
                centralBodiesToUse, bodiesToIntegrate, bodyStateFunctions,
                globalFrameOriginBarycentricFunction, globalFrameOrigin );
}

//! Function to create a translational state derivative model.
/*!
 *  Function to create a translational state derivative model from
 *  propagation settings and environment.
 *  \param translationPropagatorSettings Settings for the translational dynamics model.
 *  \param bodyMap List of body objects in the environment
 *  \param propagationStartTime Time from which numerical propagation starts.
 *  \return Translational state derivative model (instance of derived class of NBodyStateDerivative)
 */
template< typename StateScalarType = double, typename TimeType = double >
boost::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > >
createTranslationalStateDerivativeModel(
        const boost::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > >
        translationPropagatorSettings,
        const simulation_setup::NamedBodyMap& bodyMap,
        const TimeType propagationStartTime )
{

    // Create object for frame origin transformations.
    boost::shared_ptr< CentralBodyData< StateScalarType, TimeType > > centralBodyData =
            createCentralBodyData< StateScalarType, TimeType >(
                translationPropagatorSettings->centralBodies_,
                translationPropagatorSettings->bodiesToIntegrate_,
                bodyMap );

    boost::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > stateDerivativeModel;

    // Check propagator type and create corresponding state derivative object.
    switch( translationPropagatorSettings->propagator_ )
    {
    case cowell:
    {
        stateDerivativeModel = boost::make_shared<
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
            if( bodyMap.count( centralBodies[ i ] ) == 0 )
            {
                std::string errorMessage = "Error when creating Encke propagator, did not find central body " + centralBodies[ i ];
                throw std::runtime_error( errorMessage );
            }
            initialKeplerElements[ i ] = orbital_element_conversions::convertCartesianToKeplerianElements< StateScalarType >(
                        translationPropagatorSettings->getInitialStates( ).segment( i * 6, 6 ), static_cast< StateScalarType >(
                            bodyMap.at( centralBodies[ i ] )->getGravityFieldModel( )->getGravitationalParameter( ) ) );
        }

        // Create Encke state derivative object.
        stateDerivativeModel = boost::make_shared< NBodyEnckeStateDerivative< StateScalarType, TimeType > >
                ( translationPropagatorSettings->getAccelerationsMap( ), centralBodyData, translationPropagatorSettings->bodiesToIntegrate_,
                  initialKeplerElements, propagationStartTime );

        break;
    }
    case gauss_keplerian:
    {
        // Create Encke state derivative object.
        stateDerivativeModel = boost::make_shared< NBodyGaussKeplerStateDerivative< StateScalarType, TimeType > >
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
            if( bodyMap.count( centralBodies[ i ] ) == 0 )
            {
                std::string errorMessage = "Error when creating Encke propagator, did not find central body " + centralBodies[ i ];
                throw std::runtime_error( errorMessage );
            }
            initialKeplerElements.push_back( orbital_element_conversions::convertCartesianToKeplerianElements< StateScalarType >(
                                                 translationPropagatorSettings->getInitialStates( ).segment( i * 6, 6 ), static_cast< StateScalarType >(
                                                     bodyMap.at( centralBodies[ i ] )->getGravityFieldModel( )->getGravitationalParameter( ) ) ) );
        }

        // Create Encke state derivative object.:
        stateDerivativeModel = boost::make_shared< NBodyGaussModifiedEquinictialStateDerivative< StateScalarType, TimeType > >
                ( translationPropagatorSettings->getAccelerationsMap( ), centralBodyData,
                  translationPropagatorSettings->bodiesToIntegrate_, initialKeplerElements );

        break;
    }
    default:
        throw std::runtime_error(
                    "Error, did not recognize translational state propagation type: " +
                    std::to_string( translationPropagatorSettings->propagator_ ) );
    }
    return stateDerivativeModel;
}

//! Function to create a rotational dynamics state derivative model.
/*!
 *  Function to create a rotational dynamics state derivative model from propagation settings and environment.
 *  \param rotationPropagatorSettings Settings for the rotational dynamics model.
 *  \param bodyMap List of body objects in the environment
 *  \param startTime Propagation start time
 *  \return Rotational dynamics state derivative model.
 */
template< typename StateScalarType = double, typename TimeType = double >
boost::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > createRotationalStateDerivativeModel(
        const boost::shared_ptr< RotationalStatePropagatorSettings< StateScalarType > > rotationPropagatorSettings,
        const simulation_setup::NamedBodyMap& bodyMap, const TimeType startTime )
{
    std::vector< boost::function< Eigen::Matrix3d( ) > > momentOfInertiaFunctions;
    for( unsigned int i = 0; i < rotationPropagatorSettings->bodiesToIntegrate_.size( ); i++ )
    {
        momentOfInertiaFunctions.push_back(
                    boost::bind( &simulation_setup::Body::getBodyInertiaTensor,
                                 bodyMap.at( rotationPropagatorSettings->bodiesToIntegrate_.at( i ) ) ) );
    }
    return boost::make_shared< RotationalMotionStateDerivative< StateScalarType, TimeType > >(
                rotationPropagatorSettings->getTorqueModelsMap( ), rotationPropagatorSettings->bodiesToIntegrate_,
                momentOfInertiaFunctions );
}



//! Function to create a mass state derivative model.
/*!
 *  Function to create a mass state derivative model from propagation settings and environment.
 *  \param massPropagatorSettings Settings for the mass dynamics model.
 *  \param bodyMap List of body objects in the environment
 *  \return Mass state derivative model.
 */
template< typename StateScalarType = double, typename TimeType = double >
boost::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > createBodyMassStateDerivativeModel(
        const boost::shared_ptr< MassPropagatorSettings< StateScalarType > > massPropagatorSettings,
        const  simulation_setup::NamedBodyMap& bodyMap )
{
    return boost::make_shared< propagators::BodyMassStateDerivative< StateScalarType, TimeType > >(
                massPropagatorSettings->getMassRateModelsMap( ),
                massPropagatorSettings->bodiesWithMassToPropagate_ );
}

//! Function to create a state derivative model.
/*!
 *  Function to create a state derivative model from propagation settings and the environment.
 *  \param propagatorSettings Settings for the dynamical model.
 *  \param bodyMap List of body objects in the environment
 *  \param propagationStartTime Time from which numerical propagation starts.
 *  \return State derivative model (instance of required derived class of SingleStateTypeDerivative)
 */
template< typename StateScalarType = double, typename TimeType = double >
boost::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > >
createStateDerivativeModel(
        const boost::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > propagatorSettings,
        const simulation_setup::NamedBodyMap& bodyMap,
        const TimeType propagationStartTime )
{
    boost::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > stateDerivativeModel;

    // Check dynamics type and call associated function to create
    // specific type of state derivative model.
    switch( propagatorSettings->getStateType( ) )
    {
    case transational_state:
    {
        // Check input consistency.
        boost::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > >
                translationPropagatorSettings =
                boost::dynamic_pointer_cast<
                TranslationalStatePropagatorSettings< StateScalarType > >( propagatorSettings );
        if( translationPropagatorSettings == NULL )
        {
            throw std::runtime_error(
                        "Error, expected translational state propagation settings when making state derivative model" );
        }
        else
        {
            stateDerivativeModel = createTranslationalStateDerivativeModel< StateScalarType, TimeType >(
                        translationPropagatorSettings, bodyMap, propagationStartTime );
        }
        break;
    }
    case rotational_state:
    {
        boost::shared_ptr< RotationalStatePropagatorSettings< StateScalarType > > rotationPropagatorSettings =
                boost::dynamic_pointer_cast< RotationalStatePropagatorSettings< StateScalarType > >( propagatorSettings );
        if( rotationPropagatorSettings == NULL )
        {
            throw std::runtime_error(
                        "Error, expected rotation state propagation settings when making state derivative model" );
        }
        else
        {
            stateDerivativeModel = createRotationalStateDerivativeModel< StateScalarType, TimeType >(
                        rotationPropagatorSettings, bodyMap, propagationStartTime );
        }
        break;
    }
    case body_mass_state:
    {
        // Check input consistency.
        boost::shared_ptr< MassPropagatorSettings< StateScalarType > > massPropagatorSettings =
                boost::dynamic_pointer_cast< MassPropagatorSettings< StateScalarType > >( propagatorSettings );
        if( massPropagatorSettings == NULL )
        {
            throw std::runtime_error(
                        "Error, expected mass propagation settings when making state derivative model" );
        }
        else
        {
            stateDerivativeModel = createBodyMassStateDerivativeModel< StateScalarType, TimeType >(
                        massPropagatorSettings, bodyMap );
        }
        break;
    }
    case custom_state:
    {
        // Check input consistency.
        boost::shared_ptr< CustomStatePropagatorSettings< StateScalarType, TimeType > > customPropagatorSettings =
                boost::dynamic_pointer_cast< CustomStatePropagatorSettings< StateScalarType, TimeType > >(
                    propagatorSettings );
        if( customPropagatorSettings == NULL )
        {
            throw std::runtime_error(
                        "Error, expected custom propagation settings when making state derivative model" );
        }
        else
        {
            stateDerivativeModel = boost::make_shared< CustomStateDerivative< StateScalarType, TimeType > >(
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
 *  \param propagatorSettings Settings for teh numerical propagation
 *  \param bodyMap List of body objects that comprises the environment
 */
template< typename StateScalarType = double >
void setMultiTypePropagationClosure(
        const boost::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > propagatorSettings,
        const simulation_setup::NamedBodyMap& bodyMap )
{
    // Cast to multi-type settings, and perform closure if
    boost::shared_ptr< MultiTypePropagatorSettings< StateScalarType > > multiTypePropagatorSettings =
            boost::dynamic_pointer_cast< MultiTypePropagatorSettings< StateScalarType > >( propagatorSettings );
    if( multiTypePropagatorSettings != NULL )
    {
        // Perform closure for the case where both translational and rotational states are propagated
        if( multiTypePropagatorSettings->propagatorSettingsMap_.count( transational_state ) > 0 &&
                multiTypePropagatorSettings->propagatorSettingsMap_.count( rotational_state ) > 0 )
        {
            std::vector< boost::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > > translationalStateSettings =
                    multiTypePropagatorSettings->propagatorSettingsMap_.at( transational_state );
            std::vector< boost::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > > rotationalStateSettings =
                    multiTypePropagatorSettings->propagatorSettingsMap_.at( rotational_state );

            // Iterate over all accelerations, and identify those bodies for which an aerodynamic acceleration is applied
            std::vector< std::string > bodiesWithAerodynamicAcceleration;
            for( unsigned int i = 0; i < translationalStateSettings.size( ); i++ )
            {
                boost::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > > currentTranslationalState =
                        boost::dynamic_pointer_cast< TranslationalStatePropagatorSettings< StateScalarType > >(
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
                boost::shared_ptr< RotationalStatePropagatorSettings< StateScalarType > > currentTranslationalState =
                        boost::dynamic_pointer_cast< RotationalStatePropagatorSettings< StateScalarType > >(
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
            for( unsigned int i = 0; i < bodiesWithPropagatedRotation.size( ); i++ )
            {
                boost::shared_ptr< aerodynamics::FlightConditions > currentFlightConditions =
                        bodyMap.at( bodiesWithPropagatedRotation.at( i ) )->getFlightConditions( );
                reference_frames::setAerodynamicDependentOrientationCalculatorClosure(
                            boost::bind( &simulation_setup::Body::getCurrentRotationToLocalFrame,
                                         bodyMap.at( bodiesWithPropagatedRotation.at( i ) ) ),
                            currentFlightConditions->getAerodynamicAngleCalculator( ) );
            }
        }
    }
}


//! Function to create a list of state derivative models.
/*!
 *  Function to create a list of state derivative models from
 *  propagation settings and the environment.
 *  \param propagatorSettings Settings for the dynamical model.
 *  \param bodyMap List of body objects in the environment.
 *  \param propagationStartTime Time from which numerical propagation starts.
 *  \return List of state derivative models (instances of required
 *  derived class of SingleStateTypeDerivative)
 */
template< typename StateScalarType = double, typename TimeType = double >
std::vector< boost::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > >
createStateDerivativeModels(
        const boost::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > propagatorSettings,
        const simulation_setup::NamedBodyMap& bodyMap,
        const TimeType propagationStartTime )
{
    std::vector< boost::shared_ptr< SingleStateTypeDerivative< StateScalarType, TimeType > > >
            stateDerivativeModels;

    // Check type of state derivative model and call associated create function.
    switch( propagatorSettings->getStateType( ) )
    {
    // If hybrid, call create function separately for each entry.
    case hybrid:
    {
        boost::shared_ptr< MultiTypePropagatorSettings< StateScalarType > > multiTypePropagatorSettings =
                boost::dynamic_pointer_cast< MultiTypePropagatorSettings< StateScalarType > >( propagatorSettings );

        // Iterate over all propagation settings
        for( typename std::map< IntegratedStateType,
             std::vector< boost::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > > >::iterator
             propagatorIterator = multiTypePropagatorSettings->propagatorSettingsMap_.begin( );
             propagatorIterator != multiTypePropagatorSettings->propagatorSettingsMap_.end( ); propagatorIterator++ )
        {
            for( unsigned int i = 0; i < propagatorIterator->second.size( ); i++ )
            {
                // Call create function for current propagation settings.
                if( propagatorIterator->first != hybrid )
                {
                    stateDerivativeModels.push_back( createStateDerivativeModel< StateScalarType, TimeType >(
                                                         propagatorIterator->second.at( i ), bodyMap, propagationStartTime ) );
                }
                else
                {
                    throw std::runtime_error(
                                "Error when making state derivative model, cannot process nested hybrid propagators" );
                }
            }
        }

        setMultiTypePropagationClosure( propagatorSettings, bodyMap );

        break;
    }
        // If not hybrid, call create function for single object directly.
    default:
        stateDerivativeModels.push_back( createStateDerivativeModel< StateScalarType, TimeType >(
                                             propagatorSettings, bodyMap, propagationStartTime ) );
    }

    return stateDerivativeModels;
}

//! Function to create an integrator to propagate the dynamics (in normalized units) in CR3BP
/*!
 *  Function to create an integrator to propagate the dynamics (in normalized units) in Circularly Restricted Three-Body Problem.
 * \param integratorSettings Settings for the numerical integration
 * \param massParameter Normalized mass parameter
 * \param initialState Initial normalized state
 * \return Integrator object for propagation of CR3BP with requested settings
 */
boost::shared_ptr< numerical_integrators::NumericalIntegrator< double, Eigen::Vector6d > > createCR3BPIntegrator(
        const boost::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const double massParameter,
        const Eigen::Vector6d& initialState );

//! Function to propagate the dynamics (in normalized units) in CR3BP
/*!
 *  Function to propagate the dynamics (in normalized units) in Circularly Restricted Three-Body Problem.
 * \param integratorSettings Settings for the numerical integration
 * \param massParameter Normalized mass parameter
 * \param initialState Initial normalized state
 * \param finalTime End time for the numerical integration
 * \return Propagated state history of normalized dynamics in CR3BP.
 */
std::map< double, Eigen::Vector6d > performCR3BPIntegration(
        const boost::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const double massParameter,
        const Eigen::Vector6d& initialState,
        const double finalTime  );

} // namespace propagators

} // namespace tudat

#endif // TUDAT_CREATESTATEDERIVATIVEMODEL_H
