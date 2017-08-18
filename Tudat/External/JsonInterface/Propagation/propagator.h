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

    // Termination settings. If not provided, stop when epoch > simulation.endEpoch
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

    // If user did not provide conditions, or if endEpoch is defined but the time condition is missing, create it
    if ( ! defined( jsonObject, Keys::termination ) ||
         ( defined( jsonObject, Keys::endEpoch ) && timeConditionMissing ) )
    {
        terminationConditions.push_back( boost::make_shared< PropagationTimeTerminationSettings >(
                                             getEpoch< double >( jsonObject, Keys::endEpoch ) ) );
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
        for ( const boost::shared_ptr< propagators::VariableSettings > variable : exportSettings->variables )
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
            getNumeric( jsonObject, Keys::options / Keys::Options::printInterval, TUDAT_NAN, true );

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

    switch ( integratedStateType )
    {
    case translational_state:
    {
        TranslationalStatePropagatorSettings< StateScalarType > defaults(
        { }, SelectedAccelerationMap( ), { }, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >( 0 ), NULL );
        singleArcPropagatorSettings = boost::make_shared< TranslationalStatePropagatorSettings< StateScalarType > >(
                    getValue< std::vector< std::string > >( jsonObject, K::centralBodies ),
                    getValue< SelectedAccelerationMap >( jsonObject, K::accelerations ),
                    getValue< std::vector< std::string > >( jsonObject, K::bodiesToPropagate ),
                    getValue( jsonObject, K::initialStates, defaults.getInitialStates( ) ),
                    terminationSettings,
                    getValue( jsonObject, K::type, defaults.propagator_ ) );
        return;
    }
    case body_mass_state:
    {
        MassPropagatorSettings< StateScalarType > defaults(
        { }, SelectedMassRateModelMap( ), Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >( ), NULL );
        singleArcPropagatorSettings = boost::make_shared< MassPropagatorSettings< StateScalarType > >(
                    getValue< std::vector< std::string > >( jsonObject, K::bodiesToPropagate ),
                    getValue< SelectedMassRateModelMap >( jsonObject, K::massRateModels ),
                    getValue( jsonObject, K::initialStates, defaults.getInitialStates( ) ),
                    terminationSettings );
        return;
    }
    case rotational_state:
    {
        RotationalStatePropagatorSettings< StateScalarType > defaults(
                    SelectedTorqueMap( ), { }, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >( ), NULL );
        singleArcPropagatorSettings = boost::make_shared< RotationalStatePropagatorSettings< StateScalarType > >(
                    getValue< SelectedTorqueMap >( jsonObject, K::torques ),
                    getValue< std::vector< std::string > >( jsonObject, K::bodiesToPropagate ),
                    getValue( jsonObject, K::initialStates, defaults.getInitialStates( ) ),
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


namespace json_interface
{

//! Infer initial states for \p multiTypePropagatorSettings (if not provided) from ephemeris of bodies in \p bodyMap
//! at the initial time of \p integratorSettings.
/*!
 * Infer initial states for \p multiTypePropagatorSettings (if not provided) from ephemeris of bodies in \p bodyMap
 * at the initial time of \p integratorSettings.
 * \param multiTypePropagatorSettings Propagator settings whose initial states may be updated (returned by reference).
 * \param bodyMap Body map containing only bodies to be propagated with valid ephemeris.
 * \param integratorSettings Integrator settings containing the initial epoch for the ephemeris to be used.
 * \throws std::exception If initial states were not provided and could not be inferred. This happens when \p
 * multiTypePropagatorSettings contains more than one propagator, or when it contains no translational propagator,
 * or when some of the bodies to integrate contains no (valid) ephemeris (loaded by Spice).
 */
template< typename TimeType, typename StateScalarType >
void inferInitialStatesIfNecessary(
        boost::shared_ptr< propagators::MultiTypePropagatorSettings< StateScalarType > >& multiTypePropagatorSettings,
        const simulation_setup::NamedBodyMap& bodyMap,
        const boost::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > >& integratorSettings )
{
    using namespace propagators;

    // Try to get initial states from bodies' epehemeris if not provided
    if ( multiTypePropagatorSettings->getInitialStates( ).rows( ) == 0 )
    {
        if ( multiTypePropagatorSettings->propagatorSettingsMap_.size( ) == 1 )
        {
            try
            {
                const std::vector< boost::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > >
                        translationalPropagators =
                        multiTypePropagatorSettings->propagatorSettingsMap_.at( transational_state );
                if ( translationalPropagators.size( ) == 1 )
                {
                    const boost::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > >
                            translationalPropagator = boost::dynamic_pointer_cast<
                            TranslationalStatePropagatorSettings< StateScalarType > >(
                                translationalPropagators.front( ) );
                    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > systemInitialState =
                            getInitialStatesOfBodies( translationalPropagator->bodiesToIntegrate_,
                                                      translationalPropagator->centralBodies_,
                                                      bodyMap, integratorSettings->initialTime_ );
                    translationalPropagator->resetInitialStates( systemInitialState );
                    multiTypePropagatorSettings->resetInitialStates( systemInitialState );
                }
            }
            catch ( ... ) { }
        }
    }

    // If initial states could not be inferred, propagation cannot be run
    if ( multiTypePropagatorSettings->getInitialStates( ).rows( ) == 0 )
    {
        std::cerr << "Could not determine initial integration state. "
                  << "Please provide it to the propagator settings using the key \""
                  << Keys::Propagator::initialStates << "\"." << std::endl;
        throw;
    }
}

} // namespace json_interface

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_PROPAGATOR_H
