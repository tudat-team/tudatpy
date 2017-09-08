/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_PROPAGATOR_H
#define TUDAT_JSONINTERFACE_PROPAGATOR_H

#include <Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h>
#include <Tudat/SimulationSetup/PropagationSetup/propagationSettings.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/sphericalBodyShapeModel.h>

#include "Tudat/External/JsonInterface/Support/valueAccess.h"
#include "Tudat/External/JsonInterface/Support/valueConversions.h"

#include "termination.h"
#include "variable.h"
#include "acceleration.h"
#include "massRateModel.h"
#include "torque.h"
#include "export.h"

namespace tudat
{

namespace json_interface
{

// StateType

//! Possible ways of providing initial translational states.
enum StateType
{
    cartesianComponents,
    keplerianComponents
};

//! Map of `StateType` string representations.
static std::map< StateType, std::string > stateTypes =
{
    { cartesianComponents, "cartesian" },
    { keplerianComponents, "keplerian" }
};

//! `StateType` not supported by `json_interface`.
static std::vector< StateType > unsupportedStateTypes = { };

//! Convert `StateType` to `json`.
inline void to_json( json& jsonObject, const StateType& stateType )
{
    jsonObject = json_interface::stringFromEnum( stateType, stateTypes );
}

//! Convert `json` to `StateType`.
inline void from_json( const json& jsonObject, StateType& stateType )
{
    stateType = json_interface::enumFromString( jsonObject.get< std::string >( ), stateTypes );
}


//! Get the associated key (defined in a body JSON object) for an integrated state.
/*!
 * @copybrief getAssociatedKey
 * \param integratedStateType The integrated state type for which the associated key is requested.
 * \return The associated key for the integrated state.
 */
std::string getAssociatedKey( const propagators::IntegratedStateType integratedStateType )
{
    switch ( integratedStateType )
    {
    case propagators::translational_state:
        return Keys::Body::initialState;
    case propagators::body_mass_state:
        return Keys::Body::mass;
    case propagators::rotational_state:
        return Keys::Body::rotationalState;
    default:
        std::cerr << "Unknown associated JSON key for state type " << integratedStateType << std::endl;
        throw;
    }
}

//! Get the system initial state for a \p jsonObject containing the settings for a propagator.
/*!
 * @copybrief getInitialStates If the property "initialStates" is defined, this value will be returned.
 * If it is not defined, the function will try to look for the states in the properties defined in the root object at
 * "bodies". If this process fails, an empty matrix of size 0x1 will be returned.
 * \param jsonObject The `json` object containing the settings for a propagator (and optionally the root object).
 * \return Initial cobined state for the bodies to be propagated.
 */
template< typename StateScalarType >
Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > getInitialStates( const json& jsonObject )
{
    using namespace propagators;
    using namespace json_interface;
    using K = Keys::Propagator;

    if ( defined( jsonObject, K::initialStates ) )
    {
        return getValue< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >( jsonObject, K::initialStates );
    }
    else
    {
        // Integrated state type
        const IntegratedStateType integratedStateType =
                getValue( jsonObject, K::integratedStateType, translational_state );

        // State size and associated stateKey
        const unsigned int stateSize = getSingleIntegrationSize( integratedStateType );
        const std::string stateKey = getAssociatedKey( integratedStateType );

        // Bodies to propagate
        const std::vector< std::string > bodiesToPropagate =
                getValue< std::vector< std::string > >( jsonObject, K::bodiesToPropagate );

        // System initial state
        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialStates =
                Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero(
                    stateSize * bodiesToPropagate.size( ), 1 );

        // Get state for each body
        for ( unsigned int i = 0; i < bodiesToPropagate.size( ); ++i )
        {
            initialStates.segment( i * stateSize, stateSize ) =
                    getValue< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >(
                        jsonObject, SpecialKeys::root / Keys::bodies / bodiesToPropagate.at( i ) / stateKey );
        }

        return initialStates;
    }
}

//! Determine initial states for the propagator object contained in \p jsonObject (if not provided).
/*!
 * Determine initial states for the propagator object contained in \p jsonObject (if not provided).
 * The initial states can be inferred either from the state properties of the body settings (e.g. body.initialState,
 * body.mass, etc.) or from the ephemeris of the body objects in \p bodyMap at the initial time determined from
 * \p integratorSettings. If the initial states cannot be inferred, the initialStates of the propagators in \p
 * jsonObject won't be updated.
 * \param jsonObject The root `json` object to be updated with the inferred initial states (returned by reference).
 * \param bodyMap Body map containing only bodies to be propagated with valid ephemeris.
 * \param integratorSettings Integrator settings containing the initial epoch for the ephemeris to be used.
 */
template< typename TimeType, typename StateScalarType >
void determineInitialStates(
        json& jsonObject,
        const simulation_setup::NamedBodyMap& bodyMap,
        const boost::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > >& integratorSettings )
{
    using namespace propagators;
    using namespace basic_astrodynamics;
    using namespace orbital_element_conversions;
    using namespace simulation_setup;
    using KP = Keys::Propagator;
    using KS = Keys::Body::State;

    // Get propagators as an array of json (even if only one object is provided)
    json jsonPropagators = jsonObject.at( Keys::propagator );
    if ( jsonPropagators.is_object( ) )
    {
        jsonPropagators = json( );
        jsonPropagators[ 0 ] = jsonObject.at( Keys::propagator );
    }

    // Update propagators at jsonObject with initial states retrieved from bodies ephemeris
    bool usedEphemeris = false;
    if ( jsonPropagators.size( ) == 1 )
    {
        json& jsonPropagator = jsonPropagators.front( );
        if ( ! defined( jsonPropagator, KP::initialStates ) )
        {
            const IntegratedStateType integratedStateType =
                    getValue( jsonPropagator, KP::integratedStateType, translational_state );

            if ( integratedStateType == translational_state )
            {
                try
                {
                    const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialStates =
                            getInitialStatesOfBodies< TimeType, StateScalarType >(
                                getValue< std::vector< std::string > >( jsonPropagator, KP::bodiesToPropagate ),
                                getValue< std::vector< std::string > >( jsonPropagator, KP::centralBodies ),
                                bodyMap,
                                integratorSettings->initialTime_ );
                    jsonPropagator[ KP::initialStates ] = initialStates;
                    usedEphemeris = true;
                } catch ( ... ) { }
            }
        }
    }

    if ( ! usedEphemeris )
    {
        // Update propagators at jsonObject with initial states stored at JSON bodies settings
        for ( json& jsonPropagator : jsonPropagators )
        {
            if ( ! defined( jsonPropagator, KP::initialStates ) )
            {
                // Integrated state type
                const IntegratedStateType integratedStateType =
                        getValue( jsonPropagator, KP::integratedStateType, translational_state );

                // State size and associated stateKey
                const unsigned int stateSize = getSingleIntegrationSize( integratedStateType );
                const std::string stateKey = getAssociatedKey( integratedStateType );

                // Bodies to propagate
                const std::vector< std::string > bodiesToPropagate =
                        getValue< std::vector< std::string > >( jsonPropagator, KP::bodiesToPropagate );

                // System initial state
                Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialStates =
                        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero(
                            stateSize * bodiesToPropagate.size( ), 1 );

                // Get state for each body
                for ( unsigned int i = 0; i < bodiesToPropagate.size( ); ++i )
                {
                    const std::string bodyName = bodiesToPropagate.at( i );
                    const KeyPath stateKeyPath = Keys::bodies / bodyName / stateKey;
                    const json jsonState = getValue< json >( jsonObject, stateKeyPath );

                    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > bodyState( 0 );

                    // Instead of a vector, an object can be used to provide initial translational state
                    if ( integratedStateType == translational_state && ! isConvertibleToArray( jsonState ) )
                    {
                        const StateType stateType = getValue< StateType >( jsonState, KS::type );
                        switch ( stateType ) {
                        case cartesianComponents:
                        {
                            bodyState = Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( stateSize );
                            updateFromJSONIfDefined( bodyState( xCartesianPositionIndex ), jsonState, KS::x );
                            updateFromJSONIfDefined( bodyState( yCartesianPositionIndex ), jsonState, KS::y );
                            updateFromJSONIfDefined( bodyState( zCartesianPositionIndex ), jsonState, KS::z );
                            updateFromJSONIfDefined( bodyState( xCartesianVelocityIndex ), jsonState, KS::vx );
                            updateFromJSONIfDefined( bodyState( yCartesianVelocityIndex ), jsonState, KS::vy );
                            updateFromJSONIfDefined( bodyState( zCartesianVelocityIndex ), jsonState, KS::vz );
                            break;
                        }
                        case keplerianComponents:
                        {
                            // Get central body
                            std::string centralBodyName;
                            if ( defined( jsonState, KS::centralBody ) )
                            {
                                centralBodyName = getValue< std::string >( jsonState, KS::centralBody );
                            }
                            else
                            {
                                centralBodyName = getValue< std::string >( jsonPropagator, KP::centralBodies / i );
                            }
                            const boost::shared_ptr< Body > centralBody = bodyMap.at( centralBodyName );
                            const double mu = centralBody->getGravityFieldModel( )->getGravitationalParameter( );
                            const double R = centralBody->getShapeModel( )->getAverageRadius( );
                            bool usedAverageRadius = false;


                            // Detemrine semiMajorAxis, eccentricity, argumentOfPeriapsis

                            StateScalarType semiMajorAxis = TUDAT_NAN;
                            StateScalarType eccentricity = TUDAT_NAN;
                            StateScalarType argumentOfPeriapsis = TUDAT_NAN;

                            if ( defined( jsonState, KS::radius ) || defined( jsonState, KS::altitude ) )  // circular
                            {
                                if ( defined( jsonState, KS::altitude ) )
                                {
                                    semiMajorAxis = R + getValue< StateScalarType >( jsonState, KS::altitude );
                                    usedAverageRadius = true;
                                }
                                else
                                {
                                    semiMajorAxis = getValue< StateScalarType >( jsonState, KS::radius );
                                }
                                eccentricity = 0.0;
                                argumentOfPeriapsis = 0.0;
                            }
                            else  // generic
                            {
                                argumentOfPeriapsis = getValue< StateScalarType >( jsonState, KS::argumentOfPeriapsis, 0.0 );

                                if ( defined( jsonState, KS::apoapsisDistance ) || defined( jsonState, KS::apoapsisAltitude ) ||
                                     defined( jsonState, KS::periapsisDistance ) || defined( jsonState, KS::periapsisAltitude ) )
                                {
                                    StateScalarType apoapsisDistance = TUDAT_NAN;
                                    if ( defined( jsonState, KS::apoapsisAltitude ) )
                                    {
                                        apoapsisDistance = R + getValue< StateScalarType >( jsonState, KS::apoapsisAltitude );
                                        usedAverageRadius = true;
                                    }
                                    else if ( defined( jsonState, KS::apoapsisDistance ) )
                                    {
                                        apoapsisDistance = getValue< StateScalarType >( jsonState, KS::apoapsisDistance );
                                    }

                                    StateScalarType periapsisDistance = TUDAT_NAN;
                                    if ( defined( jsonState, KS::periapsisAltitude ) )
                                    {
                                        periapsisDistance = R + getValue< StateScalarType >( jsonState, KS::periapsisAltitude );
                                        usedAverageRadius = true;
                                    }
                                    else if ( defined( jsonState, KS::periapsisDistance ) )
                                    {
                                        periapsisDistance = getValue< StateScalarType >( jsonState, KS::periapsisDistance );
                                    }

                                    if ( ! isNaN( apoapsisDistance) && ! isNaN( periapsisDistance ) )
                                    {
                                        // r_a, r_p -> a, e
                                        semiMajorAxis = ( apoapsisDistance + periapsisDistance ) / 2.0;
                                        eccentricity = ( apoapsisDistance - periapsisDistance ) / ( 2.0 * semiMajorAxis );
                                    }
                                    else if ( ! isNaN( apoapsisDistance) )
                                    {
                                        if ( defined( jsonState, KS::semiMajorAxis ) )
                                        {
                                            semiMajorAxis = getValue< StateScalarType >( jsonState, KS::semiMajorAxis );
                                            // r_a, a -> e
                                            eccentricity = apoapsisDistance / semiMajorAxis - 1.0;
                                        }
                                        else
                                        {
                                            eccentricity = getValue< StateScalarType >( jsonState, KS::eccentricity );
                                            // r_a, e -> a
                                            semiMajorAxis = apoapsisDistance / ( 1.0 + eccentricity );
                                        }
                                    }
                                    else
                                    {
                                        if ( defined( jsonState, KS::semiMajorAxis ) )
                                        {
                                            semiMajorAxis = getValue< StateScalarType >( jsonState, KS::semiMajorAxis );
                                            // r_p, a -> e
                                            eccentricity = 1.0 - periapsisDistance / semiMajorAxis;
                                        }
                                        else
                                        {
                                            eccentricity = getValue< StateScalarType >( jsonState, KS::eccentricity );
                                            // r_p, e -> a
                                            semiMajorAxis = periapsisDistance / ( 1.0 - eccentricity );
                                        }
                                    }
                                }
                                else
                                {
                                    eccentricity = getValue< StateScalarType >( jsonState, KS::eccentricity, 0.0 );

                                    if ( defined( jsonState, KS::semiLatusRectum ) )
                                    {
                                        semiMajorAxis = getValue< StateScalarType >( jsonState, KS::semiLatusRectum ) /
                                                ( 1.0 - std::pow( eccentricity, 2.0 ) );
                                    }
                                    else if ( defined( jsonState, KS::meanMotion ) || defined( jsonState, KS::period ) )
                                    {
                                        StateScalarType meanMotion;
                                        if ( defined( jsonState, KS::meanMotion ) )
                                        {
                                            meanMotion = getValue< StateScalarType >( jsonState, KS::meanMotion );
                                        }
                                        else
                                        {
                                            meanMotion = 2.0 * M_PI / getValue< StateScalarType >( jsonState, KS::period );
                                        }
                                        semiMajorAxis = std::pow( ( mu / std::pow( meanMotion, 2.0 ) ), 1.0/3.0 );
                                    }
                                    else
                                    {
                                        semiMajorAxis = getValue< StateScalarType >( jsonState, KS::semiMajorAxis );
                                    }
                                }
                            }


                            // Determine inclination, longitudeOfAscendingNode

                            StateScalarType inclination = getValue< StateScalarType >( jsonState, KS::inclination, 0.0 );
                            StateScalarType longitudeOfAscendingNode = getValue< StateScalarType >( jsonState, KS::longitudeOfAscendingNode, 0.0 );


                            // Determine trueAnomaly

                            StateScalarType trueAnomaly;
                            if ( defined( jsonState, KS::meanAnomaly ) || defined( jsonState, KS::eccentricAnomaly ) )
                            {
                                StateScalarType eccentricAnomaly;
                                if ( defined( jsonState, KS::meanAnomaly ) )
                                {
                                    eccentricAnomaly = convertMeanAnomalyToEccentricAnomaly(
                                                eccentricity,
                                                getValue< StateScalarType >( jsonState, KS::meanAnomaly ) );
                                }
                                else
                                {
                                    eccentricAnomaly = getValue< StateScalarType >( jsonState, KS::eccentricAnomaly );
                                }
                                trueAnomaly = convertEccentricAnomalyToTrueAnomaly( eccentricAnomaly, eccentricity );
                            }
                            else
                            {
                                trueAnomaly = getValue< StateScalarType >( jsonState, KS::trueAnomaly, 0.0 );
                            }


                            // Full state

                            Eigen::Matrix< StateScalarType, 6, 1 > keplerianElements;
                            keplerianElements( semiMajorAxisIndex ) = semiMajorAxis;
                            keplerianElements( eccentricityIndex ) = eccentricity;
                            keplerianElements( inclinationIndex ) = inclination;
                            keplerianElements( argumentOfPeriapsisIndex ) = argumentOfPeriapsis;
                            keplerianElements( longitudeOfAscendingNodeIndex ) = longitudeOfAscendingNode;
                            keplerianElements( trueAnomalyIndex ) = trueAnomaly;

                            if ( usedAverageRadius && ! boost::dynamic_pointer_cast< SphericalBodyShapeModel >( centralBody->getShapeModel( ) ) )
                            {
                                std::cout << "Using average radius of a non-spherical body (" << centralBodyName
                                          << ") to determine the initial state of " << bodyName << "." << std::endl;
                            }

                            // Convert to Cartesian elements
                            bodyState = convertKeplerianToCartesianElements( keplerianElements, mu );
                            break;
                        }
                        default:
                            handleUnimplementedEnumValue( stateType, stateTypes, unsupportedStateTypes );
                        }
                    }

                    // Could not get the state as Cartesian, Keplerian components... then try to convert directly
                    if ( bodyState.rows( ) == 0 )
                    {
                        bodyState = getAs< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >( jsonState );
                    }

                    // Update initial states segment for current body
                    initialStates.segment( i * stateSize, stateSize ) = bodyState;
                }

                jsonPropagator[ KP::initialStates ] = initialStates;
            }
        }
    }

    jsonObject[ Keys::propagator ] = jsonPropagators;
}

} // namespace json_interface


namespace propagators
{

// IntegratedStateType

//! Map of `IntegratedStateType`s string representations.
static std::map< IntegratedStateType, std::string > integratedStateTypes =
{
    { hybrid, "hybrid" },
    { translational_state, "translational" },
    { rotational_state, "rotational" },
    { body_mass_state, "mass" },
    { relativistic_time_rate, "relativisticTimeRate" },
    { custom_state, "custom" }
};

//! `IntegratedStateType`s not supported by `json_interface`.
static std::vector< IntegratedStateType > unsupportedIntegratedStateTypes =
{
    relativistic_time_rate,
    custom_state
};

//! Convert `IntegratedStateType` to `json`.
inline void to_json( json& jsonObject, const IntegratedStateType& integratedStateType )
{
    jsonObject = json_interface::stringFromEnum( integratedStateType, integratedStateTypes );
}

//! Convert `json` to `IntegratedStateType`.
inline void from_json( const json& jsonObject, IntegratedStateType& integratedStateType )
{
    integratedStateType = json_interface::enumFromString( jsonObject.get< std::string >( ), integratedStateTypes );
}


// TranslationalPropagatorType

//! Map of `TranslationalPropagatorType`s string representations.
static std::map< TranslationalPropagatorType, std::string > translationalPropagatorTypes =
{
    { cowell, "cowell" },
    { encke, "encke" },
    { gauss_keplerian, "gaussKeplerian" },
    { gauss_modified_equinoctial, "gaussModifiedEquinoctial" }
};

//! `TranslationalPropagatorType`s not supported by `json_interface`.
static std::vector< TranslationalPropagatorType > unsupportedTranslationalPropagatorTypes = { };

//! Convert `TranslationalPropagatorType` to `json`.
inline void to_json( json& jsonObject, const TranslationalPropagatorType& translationalPropagatorType )
{
    jsonObject = json_interface::stringFromEnum( translationalPropagatorType, translationalPropagatorTypes );
}

//! Convert `json` to `TranslationalPropagatorType`.
inline void from_json( const json& jsonObject, TranslationalPropagatorType& translationalPropagatorType )
{
    translationalPropagatorType =
            json_interface::enumFromString( jsonObject.get< std::string >( ), translationalPropagatorTypes );
}

/*
// PropagatorSettings

//! Create a `json` object from a shared pointer to a `PropagatorSettings` object.
template< typename StateScalarType >
void to_json( json& jsonObject, const boost::shared_ptr< PropagatorSettings< StateScalarType > >& propagatorSettings )
{
    if ( ! propagatorSettings )
    {
        return;
    }
    using namespace json_interface;
    using K = Keys::Propagator;

    const boost::shared_ptr< MultiArcPropagatorSettings< StateScalarType > > multiArcPropagatorSettings =
            boost::dynamic_pointer_cast< MultiArcPropagatorSettings< StateScalarType > >( propagatorSettings );
    if ( multiArcPropagatorSettings )  // return an array, each element is a SingleArc
    {
        jsonObject = multiArcPropagatorSettings->getSingleArcSettings( );
    }
    else  // assume SingleArc: return an object
    {
        const boost::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > singleArcPropagatorSettings =
                boost::dynamic_pointer_cast< SingleArcPropagatorSettings< StateScalarType > >( propagatorSettings );
        jsonObject = singleArcPropagatorSettings;
    }
}

//! Create a shared pointer to a `PropagatorSettings` object from a `json` object.
template< typename StateScalarType >
void from_json( const json& jsonObject, boost::shared_ptr< PropagatorSettings< StateScalarType > >& propagatorSettings )
{
    using namespace simulation_setup;
    using namespace json_interface;
    using K = Keys::Propagator;

    if ( jsonObject.is_array( ) )  // MultiArc: each element is a SingleArc
    {
        propagatorSettings = boost::make_shared< MultiArcPropagatorSettings< StateScalarType > >(
                    getAs< std::vector< boost::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > > >(
                        jsonObject ) );
    }
    else  // assume is_object -> SingleArc
    {
        propagatorSettings = getAs< boost::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > >( jsonObject );
    }
}
*/

// MultiTypePropagatorSettings

//! Create a `json` object from a shared pointer to a `MultiTypePropagatorSettings` object.
template< typename StateScalarType >
void to_json( json& jsonObject,
              const boost::shared_ptr< MultiTypePropagatorSettings< StateScalarType > >& multiTypePropagatorSettings )
{
    if ( ! multiTypePropagatorSettings )
    {
        return;
    }
    using namespace simulation_setup;
    using namespace json_interface;

    jsonObject[ Keys::propagator ] = getFlattenedMapValues< std::map, IntegratedStateType,
            boost::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > >(
                multiTypePropagatorSettings->propagatorSettingsMap_ );
    jsonObject[ Keys::termination ] = multiTypePropagatorSettings->getTerminationSettings( );
    if ( ! isNaN( multiTypePropagatorSettings->getPrintInterval( ) ) )
    {
        jsonObject[ Keys::options ][ Keys::Options::printInterval ] = multiTypePropagatorSettings->getPrintInterval( );
    }
}

//! Create a shared pointer to a `MultiTypePropagatorSettings` object from a `json` object.
template< typename StateScalarType >
void from_json( const json& jsonObject,
                boost::shared_ptr< MultiTypePropagatorSettings< StateScalarType > >& multiTypePropagatorSettings )
{
    using namespace simulation_setup;
    using namespace json_interface;

    // Termination settings. If not provided, stop when epoch > simulation.finalEpoch
    std::vector< boost::shared_ptr< PropagationTerminationSettings > > terminationConditions;

    // Find user-defined conditions (and determine if time condition is missing)
    bool timeConditionMissing = true;
    if ( defined( jsonObject, Keys::termination ) )
    {
        boost::shared_ptr< PropagationHybridTerminationSettings > userConditions =
                getValue< boost::shared_ptr< PropagationHybridTerminationSettings > >( jsonObject, Keys::termination );
        terminationConditions.push_back( userConditions );

        for ( boost::shared_ptr< PropagationTerminationSettings > condition : userConditions->terminationSettings_ )
        {
            boost::shared_ptr< PropagationTimeTerminationSettings > timeCondition =
                    boost::dynamic_pointer_cast< PropagationTimeTerminationSettings >( condition );
            if ( timeCondition )
            {
                timeConditionMissing = false;
                break;
            }
        }
    }

    // If user did not provide conditions, or if finalEpoch is defined but the time condition is missing, create it
    if ( ! defined( jsonObject, Keys::termination ) ||
         ( defined( jsonObject, Keys::finalEpoch ) && timeConditionMissing ) )
    {
        terminationConditions.push_back( boost::make_shared< PropagationTimeTerminationSettings >(
                                             getValue< double >( jsonObject, Keys::finalEpoch ) ) );
    }

    // If there's only one condition in total (either user-provided time, user-provided dependent or created time)
    // use it as termination settings.
    boost::shared_ptr< PropagationTerminationSettings > terminationSettings;
    if ( terminationConditions.size( ) == 1 )
    {
        terminationSettings = terminationConditions.at( 0 );
    }
    // If there are two conditions (those provided by the user, and the created time condition), combine them into
    // hybrid termination settings satisfying any of the two conditions.
    else
    {
        terminationSettings = boost::make_shared< PropagationHybridTerminationSettings >( terminationConditions, true );
    }

    // Determine save settings from variables to be exported
    std::vector< std::string > addedVariablesIDs;
    std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
    std::vector< boost::shared_ptr< ExportSettings > > exportSettingsVector;
    updateFromJSONIfDefined( exportSettingsVector, jsonObject, Keys::xport );
    for ( const boost::shared_ptr< ExportSettings > exportSettings : exportSettingsVector )
    {
        for ( const boost::shared_ptr< propagators::VariableSettings > variable : exportSettings->variables_ )
        {
            boost::shared_ptr< SingleDependentVariableSaveSettings > dependentVariable =
                    boost::dynamic_pointer_cast< SingleDependentVariableSaveSettings >( variable );
            if ( dependentVariable )
            {
                const std::string variableID = getVariableId( dependentVariable );
                if ( ! contains( addedVariablesIDs, variableID ) )
                {
                    addedVariablesIDs.push_back( variableID );
                    dependentVariables.push_back( dependentVariable );
                }
            }
        }
    }

    // Print interval
    const double printInterval =
            getValue< double >( jsonObject, Keys::options / Keys::Options::printInterval, TUDAT_NAN );

    multiTypePropagatorSettings = boost::make_shared< MultiTypePropagatorSettings< StateScalarType > >(
                getValue< std::vector< boost::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > > >(
                    jsonObject, Keys::propagator ),
                terminationSettings,
                boost::make_shared< DependentVariableSaveSettings >( dependentVariables, false ),
                printInterval );
}


// SingleArcPropagatorSettings

//! Create a `json` object from a shared pointer to a `SingleArcPropagatorSettings` object.
template< typename StateScalarType >
void to_json( json& jsonObject,
              const boost::shared_ptr< SingleArcPropagatorSettings< StateScalarType > >& singleArcPropagatorSettings )
{
    if ( ! singleArcPropagatorSettings )
    {
        return;
    }
    using namespace simulation_setup;
    using namespace json_interface;
    using K = Keys::Propagator;

    // Common keys
    const IntegratedStateType integratedStateType = singleArcPropagatorSettings->getStateType( );
    jsonObject[ K::integratedStateType ] = integratedStateType;
    if ( singleArcPropagatorSettings->getInitialStates( ).rows( ) > 0 )
    {
        jsonObject[ K::initialStates ] = singleArcPropagatorSettings->getInitialStates( );
    }

    switch ( integratedStateType )
    {
    case translational_state:
    {
        boost::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > >
                translationalStatePropagatorSettings = boost::dynamic_pointer_cast<
                TranslationalStatePropagatorSettings< StateScalarType > >( singleArcPropagatorSettings );
        enforceNonNullPointer( translationalStatePropagatorSettings );
        jsonObject[ K::type ] = translationalStatePropagatorSettings->propagator_;
        jsonObject[ K::centralBodies ] = translationalStatePropagatorSettings->centralBodies_;
        jsonObject[ K::bodiesToPropagate ] = translationalStatePropagatorSettings->bodiesToIntegrate_;
        jsonObject[ K::accelerations ] = translationalStatePropagatorSettings->accelerationSettingsMap_;
        return;
    }
    case body_mass_state:
    {
        boost::shared_ptr< MassPropagatorSettings< StateScalarType > > massPropagatorSettings =
                boost::dynamic_pointer_cast< MassPropagatorSettings< StateScalarType > >( singleArcPropagatorSettings );
        enforceNonNullPointer( massPropagatorSettings );
        jsonObject[ K::bodiesToPropagate ] = massPropagatorSettings->bodiesWithMassToPropagate_;
        jsonObject[ K::massRateModels ] = massPropagatorSettings->massRateSettingsMap_;
        return;
    }
    case rotational_state:
    {
        boost::shared_ptr< RotationalStatePropagatorSettings< StateScalarType > > rotationalStatePropagatorSettings =
                boost::dynamic_pointer_cast< RotationalStatePropagatorSettings< StateScalarType > >(
                    singleArcPropagatorSettings );
        enforceNonNullPointer( rotationalStatePropagatorSettings );
        jsonObject[ K::bodiesToPropagate ] = rotationalStatePropagatorSettings->bodiesToIntegrate_;
        jsonObject[ K::torques ] = rotationalStatePropagatorSettings->torqueSettingsMap_;
        return;
    }
    case hybrid:
    {
        std::cerr << "Multitype (hybrid) propagation is implicitly supported by providing a list of propagators, "
                  << "but multitype propagators cannot be nested inside multitype propagators." << std::endl;
        throw;
    }
    default:
        handleUnimplementedEnumValue( integratedStateType, integratedStateTypes, unsupportedIntegratedStateTypes );
    }
}

//! Create a shared pointer to a `SingleArcPropagatorSettings` object from a `json` object.
template< typename StateScalarType >
void from_json( const json& jsonObject,
                boost::shared_ptr< SingleArcPropagatorSettings< StateScalarType > >& singleArcPropagatorSettings )
{
    using namespace simulation_setup;
    using namespace json_interface;
    using K = Keys::Propagator;

    // Integrated state type
    const IntegratedStateType integratedStateType =
            getValue( jsonObject, K::integratedStateType, translational_state );

    // No termination settings ( epoch > TUDAT_NAN will always be false )
    boost::shared_ptr< PropagationTerminationSettings > terminationSettings =
            boost::make_shared< PropagationTimeTerminationSettings >( TUDAT_NAN );

    // Named of bodies to propagate
    const std::vector< std::string > bodiesToPropagate =
            getValue< std::vector< std::string > >( jsonObject, K::bodiesToPropagate );

    // Initial states
    const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialStates =
            getValue< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >( jsonObject, K::initialStates );

    switch ( integratedStateType )
    {
    case translational_state:
    {
        TranslationalStatePropagatorSettings< StateScalarType > defaults(
        { }, SelectedAccelerationMap( ), { }, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >( ),
                    boost::shared_ptr< PropagationTerminationSettings >( ) );
        singleArcPropagatorSettings = boost::make_shared< TranslationalStatePropagatorSettings< StateScalarType > >(
                    getValue< std::vector< std::string > >( jsonObject, K::centralBodies ),
                    getValue< SelectedAccelerationMap >( jsonObject, K::accelerations ),
                    bodiesToPropagate,
                    initialStates,
                    terminationSettings,
                    getValue( jsonObject, K::type, defaults.propagator_ ) );
        return;
    }
    case body_mass_state:
    {
        MassPropagatorSettings< StateScalarType > defaults(
        { }, SelectedMassRateModelMap( ), Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >( ), NULL );
        singleArcPropagatorSettings = boost::make_shared< MassPropagatorSettings< StateScalarType > >(
                    bodiesToPropagate,
                    getValue< SelectedMassRateModelMap >( jsonObject, K::massRateModels ),
                    initialStates,
                    terminationSettings );
        return;
    }
    case rotational_state:
    {
        RotationalStatePropagatorSettings< StateScalarType > defaults(
                    SelectedTorqueMap( ), { }, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >( ), NULL );
        singleArcPropagatorSettings = boost::make_shared< RotationalStatePropagatorSettings< StateScalarType > >(
                    getValue< SelectedTorqueMap >( jsonObject, K::torques ),
                    bodiesToPropagate,
                    initialStates,
                    terminationSettings );
        return;
    }
    case hybrid:
    {
        std::cerr << "Multitype (hybrid) propagation is implicitly supported by providing a list of propagators, "
                  << "but multitype propagators cannot be nested inside multitype propagators." << std::endl;
        throw;
    }
    default:
        handleUnimplementedEnumValue( integratedStateType, integratedStateTypes, unsupportedIntegratedStateTypes );
    }
}

} // namespace propagators

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_PROPAGATOR_H
