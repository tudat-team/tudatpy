/*    Copyright (c) 2010-2018, Delft University of Technology
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

#include "Tudat/SimulationSetup/PropagationSetup/propagationSettings.h"
#include "Tudat/JsonInterface/Support/valueAccess.h"
#include "Tudat/JsonInterface/Support/valueConversions.h"

#include "Tudat/JsonInterface/Propagation/termination.h"
#include "Tudat/JsonInterface/Propagation/state.h"
#include "Tudat/JsonInterface/Propagation/variable.h"
#include "Tudat/JsonInterface/Propagation/acceleration.h"
#include "Tudat/JsonInterface/Propagation/massRateModel.h"
#include "Tudat/JsonInterface/Propagation/torque.h"
#include "Tudat/JsonInterface/Propagation/export.h"

namespace tudat
{

namespace json_interface
{

//! Get the associated key (defined in a body JSON object) for an integrated state.
/*!
 * @copybrief getAssociatedKey
 * \param integratedStateType The integrated state type for which the associated key is requested.
 * \return The associated key for the integrated state.
 */
inline std::string getAssociatedKey( const propagators::IntegratedStateType integratedStateType )
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
        nlohmann::json& jsonObject,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > >& integratorSettings )
{
    using namespace propagators;
    using namespace basic_astrodynamics;
    using namespace simulation_setup;
    using K = Keys::Propagator;

    // Get propagators
    nlohmann::json jsonPropagators = jsonObject.at( Keys::propagators );

    // Update propagators at jsonObject with initial states stored at JSON bodies settings
    for ( nlohmann::json& jsonPropagator : jsonPropagators )
    {
        if ( ! isDefined( jsonPropagator, K::initialStates ) )
        {
            // Integrated state type
            const IntegratedStateType integratedStateType =
                    getValue( jsonPropagator, K::integratedStateType, translational_state );

            // State size and associated stateKey
            const unsigned int stateSize = getSingleIntegrationSize( integratedStateType );
            const std::string stateKey = getAssociatedKey( integratedStateType );

            // Bodies to propagate
            const std::vector< std::string > bodiesToPropagate =
                    getValue< std::vector< std::string > >( jsonPropagator, K::bodiesToPropagate );

            // System initial state
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialStates =
                    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero(
                        stateSize * bodiesToPropagate.size( ), 1 );

            // Get state for each body
            for ( unsigned int i = 0; i < bodiesToPropagate.size( ); ++i )
            {
                const KeyPath stateKeyPath = Keys::bodies / bodiesToPropagate.at( i ) / stateKey;
                Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > bodyState;
                try
                {
                    if ( integratedStateType == translational_state )
                    {
                        Eigen::Matrix< StateScalarType, 6, 1 > stateToAdd;
                        std::string initialStateOrigin;

                        stateToAdd.setZero( );

                        const std::string centralBodyName = getValue< std::vector< std::string > >(
                                    jsonPropagator, K::centralBodies ).at( i );
                        try
                        {
                            const nlohmann::json jsonState = getValue< nlohmann::json >( jsonObject, stateKeyPath );
                            initialStateOrigin = getValue< std::string >(
                                        jsonState, Keys::Body::initialStateOrigin );

                            if( initialStateOrigin != centralBodyName )
                            {
                                if( centralBodyName == "SSB" )
                                {
                                    std::cerr<<"Error, found SSB as propagation origin, but not as initial state origin when reading JSON file. This is currently unsupported"<<std::endl;
                                }
                                stateToAdd = getInitialStateOfBody< TimeType, StateScalarType >(
                                            initialStateOrigin, centralBodyName,
                                            bodyMap, integratorSettings->initialTime_ );
                            }
                        }
                        catch( ... )
                        {
                            initialStateOrigin = centralBodyName;
                        }

                        std::shared_ptr< simulation_setup::Body > stateOriginBody;
                        if( bodyMap.count( initialStateOrigin ) == 0 )
                        {
                            stateOriginBody = nullptr;
                        }
                        else
                        {
                            stateOriginBody = bodyMap.at( initialStateOrigin );
                        }
                        bodyState = getCartesianState< StateScalarType >(
                                    jsonObject, stateKeyPath, stateOriginBody,
                                    integratorSettings->initialTime_ ) + stateToAdd;
                    }
                    else
                    {
                        bodyState = getValue< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >(
                                    jsonObject, stateKeyPath );
                    }
                }
                catch( ... )
                {
                    if ( integratedStateType == translational_state )
                    {
                        bodyState =
                                getInitialStateOfBody< TimeType, StateScalarType >(
                                    bodiesToPropagate.at( i ),
                                    getValue< std::vector< std::string > >( jsonPropagator, K::centralBodies ).at( i ),
                                    bodyMap,
                                    integratorSettings->initialTime_ );
                        jsonPropagator[ K::initialStates ] = initialStates;
                    }
                    else
                    {
                        throw std::runtime_error( "Error when getting initial state from JSON file, state type not recognized" );
                    }
                }
                initialStates.segment( i * stateSize, stateSize ) = bodyState;
            }

            // Update system initial states
            jsonPropagator[ K::initialStates ] = initialStates;
        }
    }

    jsonObject[ Keys::propagators ] = jsonPropagators;
}

//! Update dependent variable save settings of propagator from export settings object.
/*!
 * Update dependent variable save settings of propagator from export settings object.
 * \param propagatorSettings The propagator settings object to be updated (passed by reference).
 * \param exportSettingsVector The export settings settings object containing the variables to be saved.
 */
template< typename StateScalarType >
void resetDependentVariableSaveSettings(
        std::shared_ptr< propagators::MultiTypePropagatorSettings< StateScalarType > >& propagatorSettings,
        const std::vector< std::shared_ptr< ExportSettings > >& exportSettingsVector )
{
    using namespace propagators;

    // Determine save settings from variables to be exported
    std::vector< std::string > addedVariablesIDs;
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
    for ( const std::shared_ptr< ExportSettings > exportSettings : exportSettingsVector )
    {
        for ( const std::shared_ptr< propagators::VariableSettings > variable : exportSettings->variables_ )
        {
            std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariable =
                    std::dynamic_pointer_cast< SingleDependentVariableSaveSettings >( variable );
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

    propagatorSettings->resetDependentVariablesToSave(
                std::make_shared< DependentVariableSaveSettings >( dependentVariables, false ) );
}

//! Get end epoch for propagator. Returns `TUDAT_NAN` if there is no time termination condition.
template< typename TimeType >
TimeType getTerminationEpoch(
        const std::shared_ptr< propagators::PropagationTerminationSettings >& terminationSettings )
{
    using namespace propagators;

    std::shared_ptr< PropagationTimeTerminationSettings > timeTerminationSettings =
            std::dynamic_pointer_cast< PropagationTimeTerminationSettings >( terminationSettings );
    if ( timeTerminationSettings )
    {
        return timeTerminationSettings->terminationTime_;
    }

    std::shared_ptr< PropagationHybridTerminationSettings > hybridTerminationSettings =
            std::dynamic_pointer_cast< PropagationHybridTerminationSettings >( terminationSettings );
    if ( hybridTerminationSettings )
    {
        for ( unsigned int i = 0; i < hybridTerminationSettings->terminationSettings_.size( ); ++i )
        {
            const TimeType endEpoch =
                    getTerminationEpoch< TimeType >( hybridTerminationSettings->terminationSettings_.at( i ) );
            if ( ! isNaN( endEpoch ) )
            {
                return endEpoch;
            }
        }
    }

    return TUDAT_NAN;
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
    { custom_state, "custom" }
};

//! `IntegratedStateType`s not supported by `json_interface`.
static std::vector< IntegratedStateType > unsupportedIntegratedStateTypes =
{
    hybrid,         // propagators contained in a multi-type propagator cannot be hybrid
    custom_state
};

//! Convert `IntegratedStateType` to `json`.
inline void to_json( nlohmann::json& jsonObject, const IntegratedStateType& integratedStateType )
{
    jsonObject = json_interface::stringFromEnum( integratedStateType, integratedStateTypes );
}

//! Convert `json` to `IntegratedStateType`.
inline void from_json( const nlohmann::json& jsonObject, IntegratedStateType& integratedStateType )
{
    integratedStateType = json_interface::enumFromString( jsonObject, integratedStateTypes );
}


// TranslationalPropagatorType

//! Map of `TranslationalPropagatorType`s string representations.
static std::map< TranslationalPropagatorType, std::string > translationalPropagatorTypes =
{
    { cowell, "cowell" },
    { encke, "encke" },
    { gauss_keplerian, "gaussKeplerian" },
    { gauss_modified_equinoctial, "gaussModifiedEquinoctial" },
    { unified_state_model_quaternions, "unifiedStateModelQuaternions" },
    { unified_state_model_modified_rodrigues_parameters, "unifiedStateModelModifiedRodriguesParameters" },
    { unified_state_model_exponential_map, "unifiedStateModelExponentialMap" }
};

//! `TranslationalPropagatorType`s not supported by `json_interface`.
static std::vector< TranslationalPropagatorType > unsupportedTranslationalPropagatorTypes = { };

//! Convert `TranslationalPropagatorType` to `json`.
inline void to_json( nlohmann::json& jsonObject, const TranslationalPropagatorType& translationalPropagatorType )
{
    jsonObject = json_interface::stringFromEnum( translationalPropagatorType, translationalPropagatorTypes );
}

//! Convert `json` to `TranslationalPropagatorType`.
inline void from_json( const nlohmann::json& jsonObject, TranslationalPropagatorType& translationalPropagatorType )
{
    translationalPropagatorType = json_interface::enumFromString( jsonObject, translationalPropagatorTypes );
}


// RotationalPropagatorType

//! Map of `RotationalPropagatorType`s string representations.
static std::map< RotationalPropagatorType, std::string > rotationalPropagatorTypes =
{
    { quaternions, "quaternions" },
    { modified_rodrigues_parameters, "modifiedRodriguesParameters" },
    { exponential_map, "exponentialMap" }
};

//! `RotationalPropagatorType`s not supported by `json_interface`.
static std::vector< RotationalPropagatorType > unsupportedRotationalPropagatorTypes = { };

//! Convert `RotationalPropagatorType` to `json`.
inline void to_json( nlohmann::json& jsonObject, const RotationalPropagatorType& rotationalPropagatorType )
{
    jsonObject = json_interface::stringFromEnum( rotationalPropagatorType, rotationalPropagatorTypes );
}

//! Convert `json` to `RotationalPropagatorType`.
inline void from_json( const nlohmann::json& jsonObject, RotationalPropagatorType& rotationalPropagatorType )
{
    rotationalPropagatorType = json_interface::enumFromString( jsonObject, rotationalPropagatorTypes );
}


// MultiTypePropagatorSettings

//! Create a `json` object from a shared pointer to a `MultiTypePropagatorSettings` object.
template< typename StateScalarType >
void to_json( nlohmann::json& jsonObject,
              const std::shared_ptr< MultiTypePropagatorSettings< StateScalarType > >& multiTypePropagatorSettings )
{
    if ( ! multiTypePropagatorSettings )
    {
        return;
    }
    using namespace simulation_setup;
    using namespace json_interface;

    jsonObject[ Keys::propagators ] = getFlattenedMapValues< std::map, IntegratedStateType,
            std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > >(
                multiTypePropagatorSettings->propagatorSettingsMap_ );
    jsonObject[ Keys::termination ] = multiTypePropagatorSettings->getTerminationSettings( );
    if ( ! isNaN( multiTypePropagatorSettings->getPrintInterval( ) ) )
    {
        jsonObject[ Keys::options ][ Keys::Options::printInterval ] = multiTypePropagatorSettings->getPrintInterval( );
    }
}


//! Create a shared pointer to a `MultiTypePropagatorSettings` object from a `json` object.
template< typename StateScalarType >
void from_json( const nlohmann::json& jsonObject,
                std::shared_ptr< MultiTypePropagatorSettings< StateScalarType > >& multiTypePropagatorSettings )
{
    using namespace simulation_setup;
    using namespace json_interface;

    // Termination settings. If not provided, stop when epoch > simulation.finalEpoch
    std::vector< std::shared_ptr< PropagationTerminationSettings > > terminationConditions;

    // Find user-defined conditions (and determine if time condition is missing)
    bool timeConditionMissing = true;
    if ( isDefined( jsonObject, Keys::termination ) )
    {
        std::shared_ptr< PropagationHybridTerminationSettings > userConditions =
                getValue< std::shared_ptr< PropagationHybridTerminationSettings > >( jsonObject, Keys::termination );
        terminationConditions.push_back( userConditions );

        for ( std::shared_ptr< PropagationTerminationSettings > condition : userConditions->terminationSettings_ )
        {
            std::shared_ptr< PropagationTimeTerminationSettings > timeCondition =
                    std::dynamic_pointer_cast< PropagationTimeTerminationSettings >( condition );
            if ( timeCondition )
            {
                timeConditionMissing = false;
                break;
            }
        }
    }

    // If user did not provide conditions, or if finalEpoch is defined but the time condition is missing, create it
    if ( ! isDefined( jsonObject, Keys::termination ) ||
         ( isDefined( jsonObject, Keys::finalEpoch ) && timeConditionMissing ) )
    {
        terminationConditions.push_back( std::make_shared< PropagationTimeTerminationSettings >(
                                             getValue< double >( jsonObject, Keys::finalEpoch ) ) );
    }

    // If there's only one condition in total (either user-provided time, user-provided dependent or created time)
    // use it as termination settings.
    std::shared_ptr< PropagationTerminationSettings > terminationSettings;
    if ( terminationConditions.size( ) == 1 )
    {
        terminationSettings = terminationConditions.at( 0 );
    }
    // If there are two conditions (those provided by the user, and the created time condition), combine them into
    // hybrid termination settings satisfying any of the two conditions.
    else
    {
        terminationSettings = std::make_shared< PropagationHybridTerminationSettings >( terminationConditions, true );
    }

    multiTypePropagatorSettings = std::make_shared< MultiTypePropagatorSettings< StateScalarType > >(
                getValue< std::vector< std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > > >(
                    jsonObject, Keys::propagators ),
                terminationSettings,
                std::shared_ptr< DependentVariableSaveSettings >( ),
                getValue< double >( jsonObject, Keys::options / Keys::Options::printInterval, TUDAT_NAN ) );
}



// SingleArcPropagatorSettings

////! Create a `json` object from a shared pointer to a `MultiTypePropagatorSettings` object.
//template< typename StateScalarType >
//void to_json( nlohmann::json& jsonObject,
//              const std::shared_ptr< MultiTypePropagatorSettings< StateScalarType > >& multiTypePropagatorSettings )
//{

//}


//! Create a `json` object from a shared pointer to a `SingleArcPropagatorSettings` object.
template< typename StateScalarType >
void to_json( nlohmann::json& jsonObject,
              const std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > >& singleArcPropagatorSettings )
{
    if ( ! singleArcPropagatorSettings )
    {
        return;
    }
    using namespace simulation_setup;
    using namespace json_interface;
    using K = Keys::Propagator;

    if( singleArcPropagatorSettings->getStateType( ) == hybrid )
    {
        std::shared_ptr< MultiTypePropagatorSettings< StateScalarType > > multiTypePropagatorSettings =
                std::dynamic_pointer_cast< MultiTypePropagatorSettings< StateScalarType > >( singleArcPropagatorSettings );
        if ( ! multiTypePropagatorSettings )
        {
            return;
        }
        using namespace simulation_setup;
        using namespace json_interface;

        jsonObject[ Keys::propagators ] = getFlattenedMapValues< std::map, IntegratedStateType,
                std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > >(
                    multiTypePropagatorSettings->propagatorSettingsMap_ );
        jsonObject[ Keys::termination ] = multiTypePropagatorSettings->getTerminationSettings( );
        if ( ! isNaN( multiTypePropagatorSettings->getPrintInterval( ) ) )
        {
            jsonObject[ Keys::options ][ Keys::Options::printInterval ] = multiTypePropagatorSettings->getPrintInterval( );
        }
    }
    else
    {

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
            std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > >
                    translationalStatePropagatorSettings = std::dynamic_pointer_cast<
                    TranslationalStatePropagatorSettings< StateScalarType > >( singleArcPropagatorSettings );
            assertNonnullptrPointer( translationalStatePropagatorSettings );
            jsonObject[ K::type ] = translationalStatePropagatorSettings->propagator_;
            jsonObject[ K::centralBodies ] = translationalStatePropagatorSettings->centralBodies_;
            jsonObject[ K::bodiesToPropagate ] = translationalStatePropagatorSettings->bodiesToIntegrate_;
            jsonObject[ K::accelerations ] = translationalStatePropagatorSettings->getAccelerationSettingsMap( );
            return;
        }
        case body_mass_state:
        {
            std::shared_ptr< MassPropagatorSettings< StateScalarType > > massPropagatorSettings =
                    std::dynamic_pointer_cast< MassPropagatorSettings< StateScalarType > >( singleArcPropagatorSettings );
            assertNonnullptrPointer( massPropagatorSettings );
            jsonObject[ K::bodiesToPropagate ] = massPropagatorSettings->bodiesWithMassToPropagate_;
            jsonObject[ K::massRateModels ] = massPropagatorSettings->getMassRateSettingsMap( );
            return;
        }
        case rotational_state:
        {
            std::shared_ptr< RotationalStatePropagatorSettings< StateScalarType > > rotationalStatePropagatorSettings =
                    std::dynamic_pointer_cast< RotationalStatePropagatorSettings< StateScalarType > >(
                        singleArcPropagatorSettings );

            assertNonnullptrPointer( rotationalStatePropagatorSettings );
            jsonObject[ K::type ] = rotationalStatePropagatorSettings->propagator_;
            jsonObject[ K::bodiesToPropagate ] = rotationalStatePropagatorSettings->bodiesToIntegrate_;
            jsonObject[ K::torques ] = rotationalStatePropagatorSettings->getTorqueSettingsMap( );
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
}

//! Create a shared pointer to a `SingleArcPropagatorSettings` object from a `json` object.
template< typename StateScalarType >
void from_json( const nlohmann::json& jsonObject,
                std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > >& singleArcPropagatorSettings )
{
    using namespace simulation_setup;
    using namespace json_interface;
    using K = Keys::Propagator;

    // Integrated state type
    const IntegratedStateType integratedStateType =
            getValue( jsonObject, K::integratedStateType, translational_state );

    // Named of bodies to propagate
    const std::vector< std::string > bodiesToPropagate =
            getValue< std::vector< std::string > >( jsonObject, K::bodiesToPropagate );

    // Initial states
    const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialStates =
            getValue< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >( jsonObject, K::initialStates );

    // No termination settings ( epoch > TUDAT_NAN will always be false )
    std::shared_ptr< PropagationTerminationSettings > terminationSettings =
            std::make_shared< PropagationTimeTerminationSettings >( TUDAT_NAN );

    switch ( integratedStateType )
    {
    case translational_state:
    {
        TranslationalStatePropagatorSettings< StateScalarType > defaults(
        { }, SelectedAccelerationMap( ), { }, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >( ),
                    std::shared_ptr< PropagationTerminationSettings >( ) );
        singleArcPropagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< StateScalarType > >(
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
        { }, SelectedMassRateModelMap( ), Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >( ), nullptr );
        singleArcPropagatorSettings = std::make_shared< MassPropagatorSettings< StateScalarType > >(
                    bodiesToPropagate,
                    getValue< SelectedMassRateModelMap >( jsonObject, K::massRateModels ),
                    initialStates,
                    terminationSettings );
        return;
    }
    case rotational_state:
    {
        RotationalStatePropagatorSettings< StateScalarType > defaults(
                    SelectedTorqueMap( ), { }, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >( ), nullptr );
        singleArcPropagatorSettings = std::make_shared< RotationalStatePropagatorSettings< StateScalarType > >(
                    getValue< SelectedTorqueMap >( jsonObject, K::torques ),
                    bodiesToPropagate,
                    initialStates,
                    terminationSettings,
                    getValue( jsonObject, K::type, defaults.propagator_ ) );
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
