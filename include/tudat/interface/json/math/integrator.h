/*    Copyright (c) 2010-2019, Delft University of Technology
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

#include "tudat/math/integrators/createNumericalIntegrator.h"
#include "tudat/interface/json/support/valueAccess.h"
#include "tudat/interface/json/support/valueConversions.h"

namespace tudat
{

namespace numerical_integrators
{

//! Map of `AvailableIntegrators` string representations.
static std::map< AvailableIntegrators, std::string > integratorTypes =
{
    { rungeKutta4, "rungeKutta4" },
    { euler, "euler" },
    { rungeKuttaFixedStepSize, "rungeKuttaFixedStepSize" },
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


//! Map of `CoefficientSets` string representations.
static std::map< CoefficientSets, std::string > rungeKuttaCoefficientSets =
{
    { forwardEuler, "forwardEuler" },
    { rungeKutta4Classic, "rungeKutta4" },
    { rungeKuttaFehlberg45, "rungeKuttaFehlberg45" },
    { rungeKuttaFehlberg56, "rungeKuttaFehlberg56" },
    { rungeKuttaFehlberg78, "rungeKuttaFehlberg78" },
    { rungeKutta87DormandPrince, "rungeKutta87DormandPrince" }
};

//! `CoefficientSets` not supported by `json_interface`.
static std::vector< CoefficientSets > unsupportedRungeKuttaCoefficientSets = { };

//! Convert `CoefficientSets` to `json`.
inline void to_json( nlohmann::json& jsonObject, const CoefficientSets& rungeKuttaCoefficientSet )
{
    jsonObject = json_interface::stringFromEnum( rungeKuttaCoefficientSet, rungeKuttaCoefficientSets );
}

//! Convert `json` to `CoefficientSets`.
inline void from_json( const nlohmann::json& jsonObject, CoefficientSets& rungeKuttaCoefficientSet )
{
    rungeKuttaCoefficientSet =
            json_interface::enumFromString( jsonObject, rungeKuttaCoefficientSets );
}


//! Create a `json` object from a shared pointer to an `IntegratorSettings` object.
template< typename TimeType >
void to_json( nlohmann::json& jsonObject, const std::shared_ptr< IntegratorSettings< TimeType > >& integratorSettings )
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
    jsonObject[ K::assessTerminationOnMinorSteps ] =
            integratorSettings->assessTerminationOnMinorSteps_;

    switch ( integratorType )
    {
    case rungeKutta4:
    case euler:
        jsonObject[ K::stepSize ] = integratorSettings->initialTimeStep_;
        return;
    case rungeKuttaVariableStepSize:
    {
        // Create Runge-Kutta base object
        std::shared_ptr< RungeKuttaVariableStepSizeBaseSettings< TimeType > > rungeKuttaVariableStepSizeSettings =
                std::dynamic_pointer_cast< RungeKuttaVariableStepSizeBaseSettings< TimeType > >( integratorSettings );
        assertNonnullptrPointer( rungeKuttaVariableStepSizeSettings );

        // Check which integrator settings is requested
        if ( rungeKuttaVariableStepSizeSettings->areTolerancesDefinedAsScalar_ )
        {
            // Integrator with scalar tolerances
            std::shared_ptr< RungeKuttaVariableStepSizeSettingsScalarTolerances< TimeType > > scalarTolerancesIntegratorSettings =
                    std::dynamic_pointer_cast< RungeKuttaVariableStepSizeSettingsScalarTolerances< TimeType > >( integratorSettings );

            jsonObject[ K::rungeKuttaCoefficientSet ] =
                    stringFromEnum( scalarTolerancesIntegratorSettings->coefficientSet_, rungeKuttaCoefficientSets );
            jsonObject[ K::initialStepSize ] = scalarTolerancesIntegratorSettings->initialTimeStep_;
            jsonObject[ K::minimumStepSize ] = scalarTolerancesIntegratorSettings->minimumStepSize_;
            jsonObject[ K::maximumStepSize ] = scalarTolerancesIntegratorSettings->maximumStepSize_;
            jsonObject[ K::relativeErrorTolerance ] = scalarTolerancesIntegratorSettings->relativeErrorTolerance_;
            jsonObject[ K::absoluteErrorTolerance ] = scalarTolerancesIntegratorSettings->absoluteErrorTolerance_;
            jsonObject[ K::areTolerancesDefinedAsScalar ] = scalarTolerancesIntegratorSettings->areTolerancesDefinedAsScalar_;
            jsonObject[ K::safetyFactorForNextStepSize ] = scalarTolerancesIntegratorSettings->safetyFactorForNextStepSize_;
            jsonObject[ K::maximumFactorIncreaseForNextStepSize ] =
                    scalarTolerancesIntegratorSettings->maximumFactorIncreaseForNextStepSize_;
            jsonObject[ K::minimumFactorDecreaseForNextStepSize ] =
                    scalarTolerancesIntegratorSettings->minimumFactorDecreaseForNextStepSize_;
        }
        else
        {
            throw std::runtime_error( "Error while creating Runge-Kutta variable step-size integrator via JSON interface. RK "
                                      "integrators with vector tolerances are not yet supported via JSON." );
        }
        return;
    }
    case adamsBashforthMoulton:
    {
        std::shared_ptr< AdamsBashforthMoultonSettings< TimeType > > adamsBashforthMoultonSettings =
                std::dynamic_pointer_cast< AdamsBashforthMoultonSettings< TimeType > >( integratorSettings );
        assertNonnullptrPointer( adamsBashforthMoultonSettings );
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
        std::shared_ptr< BulirschStoerIntegratorSettings< TimeType > > bulirschStoerSettings =
                std::dynamic_pointer_cast< BulirschStoerIntegratorSettings< TimeType > >( integratorSettings );
        assertNonnullptrPointer( bulirschStoerSettings );
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
void from_json( const nlohmann::json& jsonObject, std::shared_ptr< IntegratorSettings< TimeType > >& integratorSettings )
{
    using namespace json_interface;
    using RungeKuttaCoefficientSet = CoefficientSets;
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
        integratorSettings = std::make_shared< IntegratorSettings< TimeType > >(
                    integratorType,
                    initialTime,
                    getValue< TimeType >( jsonObject, K::stepSize ),
                    getValue( jsonObject, K::saveFrequency, defaults.saveFrequency_ ),
                    getValue( jsonObject, K::assessTerminationOnMinorSteps,
                              defaults.assessTerminationOnMinorSteps_ ) );
        return;
    }
    case rungeKuttaVariableStepSize:
    {
        // Check which constructor to use
        if ( getValue< bool >( jsonObject, K::areTolerancesDefinedAsScalar, true ) )
        {
            // Scalar tolerances
            RungeKuttaVariableStepSizeSettingsScalarTolerances< TimeType > defaults(
                        0.0, 0.0, RungeKuttaCoefficientSet::rungeKuttaFehlberg45, 0.0, 0.0 );

            integratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettingsScalarTolerances< TimeType > >(
                        initialTime,
                        getValue< TimeType >( jsonObject, K::initialStepSize ),
                        getValue< RungeKuttaCoefficientSet >( jsonObject, K::rungeKuttaCoefficientSet ),
                        getValue< TimeType >( jsonObject, K::minimumStepSize ),
                        getValue< TimeType >( jsonObject, K::maximumStepSize ),
                        getValue( jsonObject, K::relativeErrorTolerance, defaults.relativeErrorTolerance_ ),
                        getValue( jsonObject, K::absoluteErrorTolerance, defaults.absoluteErrorTolerance_ ),
                        getValue( jsonObject, K::saveFrequency, defaults.saveFrequency_ ),
                        getValue( jsonObject, K::assessTerminationOnMinorSteps,
                                  defaults.assessTerminationOnMinorSteps_ ),
                        getValue( jsonObject, K::safetyFactorForNextStepSize,
                                  defaults.safetyFactorForNextStepSize_ ),
                        getValue( jsonObject, K::maximumFactorIncreaseForNextStepSize,
                                  defaults.maximumFactorIncreaseForNextStepSize_ ),
                        getValue( jsonObject, K::minimumFactorDecreaseForNextStepSize,
                                  defaults.minimumFactorDecreaseForNextStepSize_ ) );
        }
        else
        {
            throw std::runtime_error( "Error while creating Runge-Kutta variable step-size integrator from JSON object. RK "
                                      "integrators with vector tolerances are not yet supported via JSON." );
        }
        return;
    }
    case adamsBashforthMoulton:
    {
        AdamsBashforthMoultonSettings< TimeType > defaults(
                    0.0, 0.0, 0.0, 0.0 );

        integratorSettings = std::make_shared< AdamsBashforthMoultonSettings< TimeType > >(
                    initialTime,
                    getValue< TimeType >( jsonObject, K::initialStepSize ),
                    getValue< TimeType >( jsonObject, K::minimumStepSize ),
                    getValue< TimeType >( jsonObject, K::maximumStepSize ),
                    getValue( jsonObject, K::relativeErrorTolerance, defaults.relativeErrorTolerance_ ),
                    getValue( jsonObject, K::absoluteErrorTolerance, defaults.absoluteErrorTolerance_ ),
                    getValue( jsonObject, K::minimumOrder, defaults.minimumOrder_ ),
                    getValue( jsonObject, K::maximumOrder, defaults.maximumOrder_ ),
                    getValue( jsonObject, K::saveFrequency, defaults.saveFrequency_ ),
                    getValue( jsonObject, K::assessTerminationOnMinorSteps,
                              defaults.assessTerminationOnMinorSteps_ ),
                    getValue( jsonObject, K::bandwidth,
                              defaults.bandwidth_ ) );
        return;
    }
    case bulirschStoer:
    {
        BulirschStoerIntegratorSettings< TimeType > defaults(
                    0.0, 0.0, bulirsch_stoer_sequence, 6, std::numeric_limits< double >::epsilon( ),
                    std::numeric_limits< double >::infinity( ) );

        integratorSettings = std::make_shared< BulirschStoerIntegratorSettings< TimeType > >(
                    initialTime,
                    getValue< TimeType >( jsonObject, K::initialStepSize ),
                    getValue( jsonObject, K::extrapolationSequence, defaults.extrapolationSequence_ ),
                    getValue( jsonObject, K::maximumNumberOfSteps, defaults.maximumNumberOfSteps_ ),
                    getValue< TimeType >( jsonObject, K::minimumStepSize ),
                    getValue< TimeType >( jsonObject, K::maximumStepSize ),
                    getValue( jsonObject, K::relativeErrorTolerance, defaults.relativeErrorTolerance_ ),
                    getValue( jsonObject, K::absoluteErrorTolerance, defaults.absoluteErrorTolerance_ ),
                    getValue( jsonObject, K::saveFrequency, defaults.saveFrequency_ ),
                    getValue( jsonObject, K::assessTerminationOnMinorSteps,
                              defaults.assessTerminationOnMinorSteps_ ),
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
