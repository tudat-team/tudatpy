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

#include <Tudat/SimulationSetup/PropagationSetup/propagationSettings.h>

#include "Tudat/External/JsonInterface/Support/valueAccess.h"
#include "Tudat/External/JsonInterface/Support/valueConversions.h"

#include "termination.h"
#include "variable.h"
#include "acceleration.h"
// #include "massRate.h"
// #include "torque.h"

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


// PropagatorSettings

//! Create a `json` object from a shared pointer to a `PropagatorSettings` object.
template< typename StateScalarType = double >
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
template< typename StateScalarType = double >
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


// SingleArcPropagatorSettings

//! Create a `json` object from a shared pointer to a `SingleArcPropagatorSettings` object.
template< typename StateScalarType = double >
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
    jsonObject[ K::initialStates ] = singleArcPropagatorSettings->getInitialStates( );
    jsonObject[ K::termination ] = singleArcPropagatorSettings->getTerminationSettings( );
    if( singleArcPropagatorSettings->getDependentVariablesToSave( ) )
    {
        assignIfNotEmpty( jsonObject, K::computeVariables,
                          singleArcPropagatorSettings->getDependentVariablesToSave( )->dependentVariables_ );
    }
    assignIfNotNaN( jsonObject, K::printInterval, singleArcPropagatorSettings->getPrintInterval( ) );

    switch ( integratedStateType )
    {
    case hybrid:
    {
        boost::shared_ptr< MultiTypePropagatorSettings< StateScalarType > > multiTypePropagatorSettings =
                boost::dynamic_pointer_cast< MultiTypePropagatorSettings< StateScalarType > >(
                    singleArcPropagatorSettings );
        enforceNonNullPointer( multiTypePropagatorSettings );
        jsonObject[ K::propagators ] = getFlattenedMapValues< std::map, IntegratedStateType,
                boost::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > >(
                    multiTypePropagatorSettings->propagatorSettingsMap_ );
        return;
    }
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
    default:
        handleUnimplementedEnumValue( integratedStateType, integratedStateTypes, unsupportedIntegratedStateTypes );
    }
}

//! Create a shared pointer to a `SingleArcPropagatorSettings` object from a `json` object.
template< typename StateScalarType = double >
void from_json( const json& jsonObject,
                boost::shared_ptr< SingleArcPropagatorSettings< StateScalarType > >& singleArcPropagatorSettings )
{
    using namespace simulation_setup;
    using namespace json_interface;
    using K = Keys::Propagator;

    // Integrated state type
    const IntegratedStateType integratedStateType =
            getValue( jsonObject, K::integratedStateType, translational_state );

    // Termination settings. If not provided, stop when epoch > simulation.endEpoch
    const boost::shared_ptr< PropagationTerminationSettings > terminationSettings =
            getValue< boost::shared_ptr< PropagationTerminationSettings > >(
                jsonObject, K::termination, boost::make_shared< PropagationTimeTerminationSettings >(
                    getEpoch< double >( jsonObject, SpecialKeys::root / Keys::endEpoch ) ) );

    // Save settings
    boost::shared_ptr< DependentVariableSaveSettings > saveSettings;
    if ( defined( jsonObject, K::computeVariables) )
    {
        saveSettings = boost::make_shared< DependentVariableSaveSettings >(
                    getValue< std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > >(
                        jsonObject, K::computeVariables  ), false );
    }

    switch ( integratedStateType )
    {
    case hybrid:
    {
        MultiTypePropagatorSettings< StateScalarType > defaults(
                    std::vector< boost::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > >( ), NULL );
        singleArcPropagatorSettings = boost::make_shared< MultiTypePropagatorSettings< StateScalarType > >(
                    getValue< std::vector< boost::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > > >(
                        jsonObject, K::propagators ),
                    terminationSettings,
                    saveSettings,
                    getNumeric( jsonObject, K::printInterval, defaults.getPrintInterval( ), true ) );
        return;
    }
    case translational_state:
    {
        TranslationalStatePropagatorSettings< StateScalarType > defaults(
        { }, SelectedAccelerationMap( ), { }, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >( ), NULL );
        singleArcPropagatorSettings = boost::make_shared< TranslationalStatePropagatorSettings< StateScalarType > >(
                    getValue< std::vector< std::string > >( jsonObject, K::centralBodies ),
                    getValue< SelectedAccelerationMap >( jsonObject, K::accelerations ),
                    getValue< std::vector< std::string > >( jsonObject, K::bodiesToPropagate ),
                    getValue< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >( jsonObject, K::initialStates ),
                    terminationSettings,
                    getValue( jsonObject, K::type, defaults.propagator_ ),
                    saveSettings,
                    getNumeric( jsonObject, K::printInterval, defaults.getPrintInterval( ), true ) );
        return;
    }
    default:
        handleUnimplementedEnumValue( integratedStateType, integratedStateTypes, unsupportedIntegratedStateTypes );
    }
}

} // namespace propagators

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_PROPAGATOR_H
