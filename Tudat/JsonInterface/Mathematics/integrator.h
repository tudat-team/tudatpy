/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_INTEGRATOR_H
#define TUDAT_JSONINTERFACE_INTEGRATOR_H

#include "Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h"
#include "Tudat/JsonInterface/Support/valueAccess.h"
#include "Tudat/JsonInterface/Support/valueConversions.h"

namespace tudat
{

namespace numerical_integrators
{

//! Map of `AvailableIntegrators` string representations.
static std::map< AvailableIntegrators, std::string > integratorTypes =
{
    { rungeKutta4, "rungeKutta4" },
    { euler, "euler" },
    { rungeKuttaVariableStepSize, "rungeKuttaVariableStepSize" },
    { adamsBashforthMoulton, "adamsBashforthMoulton" },
    { bulirschStoer, "bulirschStoer" },
};

//! `AvailableIntegrators` not supported by `json_interface`.
static std::vector< AvailableIntegrators > unsupportedIntegratorTypes = { };

//! Convert `AvailableIntegrators` to `json`.
inline void to_json( nlohmann::json& jsonObject, const AvailableIntegrators& availableIntegrator )
{
    jsonObject = json_interface::stringFromEnum( availableIntegrator, integratorTypes );
}

//! Convert `json` to `AvailableIntegrators`.
inline void from_json( const nlohmann::json& jsonObject, AvailableIntegrators& availableIntegrator )
{
    availableIntegrator = json_interface::enumFromString( jsonObject, integratorTypes );
}


//! Map of `RungeKuttaCoefficients::CoefficientSets` string representations.
static std::map< RungeKuttaCoefficients::CoefficientSets, std::string > rungeKuttaCoefficientSets =
{
    { RungeKuttaCoefficients::rungeKuttaFehlberg45, "rungeKuttaFehlberg45" },
    { RungeKuttaCoefficients::rungeKuttaFehlberg56, "rungeKuttaFehlberg56" },
    { RungeKuttaCoefficients::rungeKuttaFehlberg78, "rungeKuttaFehlberg78" },
    { RungeKuttaCoefficients::rungeKutta87DormandPrince, "rungeKutta87DormandPrince" }
};

//! `RungeKuttaCoefficients::CoefficientSets` not supported by `json_interface`.
static std::vector< RungeKuttaCoefficients::CoefficientSets > unsupportedRungeKuttaCoefficientSets = { };

//! Convert `RungeKuttaCoefficients::CoefficientSets` to `json`.
inline void to_json( nlohmann::json& jsonObject, const RungeKuttaCoefficients::CoefficientSets& rungeKuttaCoefficientSet )
{
    jsonObject = json_interface::stringFromEnum( rungeKuttaCoefficientSet, rungeKuttaCoefficientSets );
}

//! Convert `json` to `RungeKuttaCoefficients::CoefficientSets`.
inline void from_json( const nlohmann::json& jsonObject, RungeKuttaCoefficients::CoefficientSets& rungeKuttaCoefficientSet )
{
    rungeKuttaCoefficientSet =
            json_interface::enumFromString( jsonObject, rungeKuttaCoefficientSets );
}


//! Create a `json` object from a shared pointer to an `IntegratorSettings` object.
template< typename TimeType >
void to_json( nlohmann::json& jsonObject, const boost::shared_ptr< IntegratorSettings< TimeType > >& integratorSettings )
{
    if ( ! integratorSettings )
    {
        return;
    }
    using namespace json_interface;
    using K = Keys::Integrator;

    // Common keys
    const AvailableIntegrators integratorType = integratorSettings->integratorType_;
    jsonObject[ K::type ] = integratorType;
    jsonObject[ K::initialTime ] = integratorSettings->initialTime_;
    jsonObject[ K::saveFrequency ] = integratorSettings->saveFrequency_;
    jsonObject[ K::assessPropagationTerminationConditionDuringIntegrationSubsteps ] =
            integratorSettings->assessPropagationTerminationConditionDuringIntegrationSubsteps_;

    switch ( integratorType )
    {
    case rungeKutta4:
    case euler:
        jsonObject[ K::stepSize ] = integratorSettings->initialTimeStep_;
        return;
    case rungeKuttaVariableStepSize:
    {
        boost::shared_ptr< RungeKuttaVariableStepSizeSettings< TimeType > > rungeKuttaVariableStepSizeSettings =
                boost::dynamic_pointer_cast< RungeKuttaVariableStepSizeSettings< TimeType > >( integratorSettings );
        assertNonNullPointer( rungeKuttaVariableStepSizeSettings );
        jsonObject[ K::rungeKuttaCoefficientSet ] =
                stringFromEnum( rungeKuttaVariableStepSizeSettings->coefficientSet_, rungeKuttaCoefficientSets );
        jsonObject[ K::initialStepSize ] = rungeKuttaVariableStepSizeSettings->initialTimeStep_;
        jsonObject[ K::minimumStepSize ] = rungeKuttaVariableStepSizeSettings->minimumStepSize_;
        jsonObject[ K::maximumStepSize ] = rungeKuttaVariableStepSizeSettings->maximumStepSize_;
        jsonObject[ K::relativeErrorTolerance ] = rungeKuttaVariableStepSizeSettings->relativeErrorTolerance_;
        jsonObject[ K::absoluteErrorTolerance ] = rungeKuttaVariableStepSizeSettings->absoluteErrorTolerance_;
        jsonObject[ K::safetyFactorForNextStepSize ] =
                rungeKuttaVariableStepSizeSettings->safetyFactorForNextStepSize_;
        jsonObject[ K::maximumFactorIncreaseForNextStepSize ] =
                rungeKuttaVariableStepSizeSettings->maximumFactorIncreaseForNextStepSize_;
        jsonObject[ K::minimumFactorDecreaseForNextStepSize ] =
                rungeKuttaVariableStepSizeSettings->minimumFactorDecreaseForNextStepSize_;
        return;
    }
    case adamsBashforthMoulton:
    {
        boost::shared_ptr< AdamsBashforthMoultonSettings< TimeType > > adamsBashforthMoultonSettings =
                boost::dynamic_pointer_cast< AdamsBashforthMoultonSettings< TimeType > >( integratorSettings );
        assertNonNullPointer( adamsBashforthMoultonSettings );
        jsonObject[ K::initialStepSize ] = adamsBashforthMoultonSettings->initialTimeStep_;
        jsonObject[ K::minimumStepSize ] = adamsBashforthMoultonSettings->minimumStepSize_;
        jsonObject[ K::maximumStepSize ] = adamsBashforthMoultonSettings->maximumStepSize_;
        jsonObject[ K::relativeErrorTolerance ] = adamsBashforthMoultonSettings->relativeErrorTolerance_;
        jsonObject[ K::absoluteErrorTolerance ] = adamsBashforthMoultonSettings->absoluteErrorTolerance_;
        jsonObject[ K::minimumOrder ] = adamsBashforthMoultonSettings->minimumOrder_;
        jsonObject[ K::maximumOrder ] = adamsBashforthMoultonSettings->maximumOrder_;
        jsonObject[ K::bandwidth ] = adamsBashforthMoultonSettings->bandwidth_;
        return;
    }
    case bulirschStoer:
    {
        boost::shared_ptr< BulirschStoerIntegratorSettings< TimeType > > bulirschStoerSettings =
                boost::dynamic_pointer_cast< BulirschStoerIntegratorSettings< TimeType > >( integratorSettings );
        assertNonNullPointer( bulirschStoerSettings );
        jsonObject[ K::initialStepSize ] = bulirschStoerSettings->initialTimeStep_;
        jsonObject[ K::minimumStepSize ] = bulirschStoerSettings->minimumStepSize_;
        jsonObject[ K::maximumStepSize ] = bulirschStoerSettings->maximumStepSize_;
        jsonObject[ K::relativeErrorTolerance ] = bulirschStoerSettings->relativeErrorTolerance_;
        jsonObject[ K::absoluteErrorTolerance ] = bulirschStoerSettings->absoluteErrorTolerance_;
        jsonObject[ K::extrapolationSequence ] =  bulirschStoerSettings->extrapolationSequence_;
        jsonObject[ K::maximumNumberOfSteps ] =  bulirschStoerSettings->maximumNumberOfSteps_;
        jsonObject[ K::safetyFactorForNextStepSize ] =
                bulirschStoerSettings->safetyFactorForNextStepSize_;
        jsonObject[ K::maximumFactorIncreaseForNextStepSize ] =
                bulirschStoerSettings->maximumFactorIncreaseForNextStepSize_;
        jsonObject[ K::minimumFactorDecreaseForNextStepSize ] =
                bulirschStoerSettings->minimumFactorDecreaseForNextStepSize_;

        return;
    }
    default:
        handleUnimplementedEnumValue( integratorType, integratorTypes, unsupportedIntegratorTypes );
    }
}

//! Create a `json` object from a shared pointer to an `IntegratorSettings` object.
template< typename TimeType >
void from_json( const nlohmann::json& jsonObject, boost::shared_ptr< IntegratorSettings< TimeType > >& integratorSettings )
{
    using namespace json_interface;
    using RungeKuttaCoefficientSet = RungeKuttaCoefficients::CoefficientSets;
    using K = Keys::Integrator;

    // Read JSON settings shared by all supported integrators
    const AvailableIntegrators integratorType = getValue( jsonObject, K::type, rungeKutta4 );
    const TimeType initialTime =
            getValue< TimeType >( jsonObject, { K::initialTime, SpecialKeys::root / Keys::initialEpoch } );

    // Create IntegratorSettings pointer from JSON settings
    switch ( integratorType )
    {
    case euler:
    case rungeKutta4:
    {
        IntegratorSettings< TimeType > defaults( integratorType, 0.0, 0.0 );
        integratorSettings = boost::make_shared< IntegratorSettings< TimeType > >(
                    integratorType,
                    initialTime,
                    getValue< TimeType >( jsonObject, K::stepSize ),
                    getValue( jsonObject, K::saveFrequency, defaults.saveFrequency_ ),
                    getValue( jsonObject, K::assessPropagationTerminationConditionDuringIntegrationSubsteps,
                              defaults.assessPropagationTerminationConditionDuringIntegrationSubsteps_ ) );
        return;
    }
    case rungeKuttaVariableStepSize:
    {
        RungeKuttaVariableStepSizeSettings< TimeType > defaults(
                    integratorType, 0.0, 0.0, RungeKuttaCoefficientSet::rungeKuttaFehlberg45, 0.0, 0.0 );

        integratorSettings = boost::make_shared< RungeKuttaVariableStepSizeSettings< TimeType > >(
                    integratorType,
                    initialTime,
                    getValue< TimeType >( jsonObject, K::initialStepSize ),
                    getValue< RungeKuttaCoefficientSet >( jsonObject, K::rungeKuttaCoefficientSet ),
                    getValue< TimeType >( jsonObject, K::minimumStepSize ),
                    getValue< TimeType >( jsonObject, K::maximumStepSize ),
                    getValue( jsonObject, K::relativeErrorTolerance, defaults.relativeErrorTolerance_ ),
                    getValue( jsonObject, K::absoluteErrorTolerance, defaults.absoluteErrorTolerance_ ),
                    getValue( jsonObject, K::saveFrequency, defaults.saveFrequency_ ),
                    getValue( jsonObject, K::assessPropagationTerminationConditionDuringIntegrationSubsteps,
                              defaults.assessPropagationTerminationConditionDuringIntegrationSubsteps_ ),
                    getValue( jsonObject, K::safetyFactorForNextStepSize,
                              defaults.safetyFactorForNextStepSize_ ),
                    getValue( jsonObject, K::maximumFactorIncreaseForNextStepSize,
                              defaults.maximumFactorIncreaseForNextStepSize_ ),
                    getValue( jsonObject, K::minimumFactorDecreaseForNextStepSize,
                              defaults.minimumFactorDecreaseForNextStepSize_ ) );
        return;
    }
    case adamsBashforthMoulton:
    {
        AdamsBashforthMoultonSettings< TimeType > defaults(
                    0.0, 0.0, 0.0, 0.0 );

        integratorSettings = boost::make_shared< AdamsBashforthMoultonSettings< TimeType > >(
                    initialTime,
                    getValue< TimeType >( jsonObject, K::initialStepSize ),
                    getValue< TimeType >( jsonObject, K::minimumStepSize ),
                    getValue< TimeType >( jsonObject, K::maximumStepSize ),
                    getValue( jsonObject, K::relativeErrorTolerance, defaults.relativeErrorTolerance_ ),
                    getValue( jsonObject, K::absoluteErrorTolerance, defaults.absoluteErrorTolerance_ ),
                    getValue( jsonObject, K::minimumOrder, defaults.minimumOrder_ ),
                    getValue( jsonObject, K::maximumOrder, defaults.maximumOrder_ ),
                    getValue( jsonObject, K::saveFrequency, defaults.saveFrequency_ ),
                    getValue( jsonObject, K::assessPropagationTerminationConditionDuringIntegrationSubsteps,
                              defaults.assessPropagationTerminationConditionDuringIntegrationSubsteps_ ),
                    getValue( jsonObject, K::bandwidth,
                              defaults.bandwidth_ ) );
        return;
    }
    case bulirschStoer:
    {
        BulirschStoerIntegratorSettings< TimeType > defaults(
                    0.0, 0.0, bulirsch_stoer_sequence, 6, std::numeric_limits< double >::epsilon( ),
                    std::numeric_limits< double >::infinity( ) );

        integratorSettings = boost::make_shared< BulirschStoerIntegratorSettings< TimeType > >(
                    initialTime,
                    getValue< TimeType >( jsonObject, K::initialStepSize ),
                    getValue( jsonObject, K::extrapolationSequence, defaults.extrapolationSequence_ ),
                    getValue( jsonObject, K::maximumNumberOfSteps, defaults.maximumNumberOfSteps_ ),
                    getValue< TimeType >( jsonObject, K::minimumStepSize ),
                    getValue< TimeType >( jsonObject, K::maximumStepSize ),
                    getValue( jsonObject, K::relativeErrorTolerance, defaults.relativeErrorTolerance_ ),
                    getValue( jsonObject, K::absoluteErrorTolerance, defaults.absoluteErrorTolerance_ ),
                    getValue( jsonObject, K::saveFrequency, defaults.saveFrequency_ ),
                    getValue( jsonObject, K::assessPropagationTerminationConditionDuringIntegrationSubsteps,
                              defaults.assessPropagationTerminationConditionDuringIntegrationSubsteps_ ),
                    getValue( jsonObject, K::safetyFactorForNextStepSize,
                              defaults.safetyFactorForNextStepSize_ ),
                    getValue( jsonObject, K::maximumFactorIncreaseForNextStepSize,
                              defaults.maximumFactorIncreaseForNextStepSize_ ),
                    getValue( jsonObject, K::minimumFactorDecreaseForNextStepSize,
                              defaults.minimumFactorDecreaseForNextStepSize_ ) );
        return;
    }
    default:
        handleUnimplementedEnumValue( integratorType, integratorTypes, unsupportedIntegratorTypes );
    }
}

} // namespace numerical_integrators

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_INTEGRATOR_H
