/*    Copyright (c) 2010-2017, Delft University of Technology
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

#include <Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h>

#include "jsonInterface.h"

namespace tudat
{

namespace numerical_integrators
{

//! Map of `AvailableIntegrators` supported by `json_interface`.
static std::map< std::string, AvailableIntegrators > availableIntegrators =
{
    { "rungeKutta4",                rungeKutta4 },
    { "euler",                      euler },
    { "rungeKuttaVariableStepSize", rungeKuttaVariableStepSize }
};

//! Convert `AvailableIntegrators` to `json`.
void to_json( json& jsonObject, const AvailableIntegrators& availableIntegrator );

//! Convert `json` to `AvailableIntegrators`.
void from_json( const json& jsonObject, AvailableIntegrators& availableIntegrator );

//! Map of `RungeKuttaCoefficients::CoefficientSets` supported by `json_interface`.
static std::map< std::string, RungeKuttaCoefficients::CoefficientSets > rungeKuttaCoefficientSets =
{
    { "rungeKuttaFehlberg45",      RungeKuttaCoefficients::rungeKuttaFehlberg45 },
    { "rungeKuttaFehlberg56",      RungeKuttaCoefficients::rungeKuttaFehlberg56 },
    { "rungeKuttaFehlberg78",      RungeKuttaCoefficients::rungeKuttaFehlberg78 },
    { "rungeKutta87DormandPrince", RungeKuttaCoefficients::rungeKutta87DormandPrince }
};
void to_json( json& jsonObject, const RungeKuttaCoefficients::CoefficientSets& rungeKuttaCoefficientSet );
void from_json( const json& jsonObject, RungeKuttaCoefficients::CoefficientSets& rungeKuttaCoefficientSet );

//! Create a `json` object from a shared pointer to an `IntegratorSettings` object.
//! Called automatically by `nlohmann::json` when using a constructor such as `json( integratorSettings )`.
template< typename TimeType = double >
void to_json( json& jsonObject, const boost::shared_ptr< IntegratorSettings< TimeType > >& integratorSettings )
{
    using namespace json_interface;
    using Keys = Keys::Integrator;

    // Initialise
    jsonObject = json( );

    // Common keys
    jsonObject[ Keys::type ] = stringFromEnum( integratorSettings->integratorType_, availableIntegrators );
    jsonObject[ Keys::initialTime ] = integratorSettings->initialTime_;
    jsonObject[ Keys::initialTimeStep ] = integratorSettings->initialTimeStep_;
    jsonObject[ Keys::saveFrequency ] = integratorSettings->saveFrequency_;

    // Integrator-specific keys
    boost::shared_ptr< RungeKuttaVariableStepSizeSettings< TimeType > > rungeKuttaVariableStepSizeSettings =
            boost::dynamic_pointer_cast< RungeKuttaVariableStepSizeSettings< TimeType > >( integratorSettings );
    if ( rungeKuttaVariableStepSizeSettings )
    {
        jsonObject[ Keys::rungeKuttaCoefficientSet ] =
                stringFromEnum( rungeKuttaVariableStepSizeSettings->coefficientSet_, rungeKuttaCoefficientSets );
        jsonObject[ Keys::minimumStepSize ] = rungeKuttaVariableStepSizeSettings->minimumStepSize_;
        jsonObject[ Keys::maximumStepSize ] = rungeKuttaVariableStepSizeSettings->maximumStepSize_;
        jsonObject[ Keys::relativeErrorTolerance ] = rungeKuttaVariableStepSizeSettings->relativeErrorTolerance_;
        jsonObject[ Keys::absoluteErrorTolerance ] = rungeKuttaVariableStepSizeSettings->absoluteErrorTolerance_;
        jsonObject[ Keys::safetyFactorForNextStepSize ] =
                rungeKuttaVariableStepSizeSettings->safetyFactorForNextStepSize_;
        jsonObject[ Keys::maximumFactorIncreaseForNextStepSize ] =
                rungeKuttaVariableStepSizeSettings->maximumFactorIncreaseForNextStepSize_;
        jsonObject[ Keys::minimumFactorDecreaseForNextStepSize ] =
                rungeKuttaVariableStepSizeSettings->minimumFactorDecreaseForNextStepSize_;
        return;
    }
}

} // namespace numerical_integrators


namespace json_interface
{

//! Create a shared pointer to an `IntegratorSettings` object from a `json` object.
/*!
 * Create a shared pointer to an `IntegratorSettings` object from a `json` object.
 * \param settings `json` object containing the settings for one integrator.
 * \param keyTree Key tree at which the object containing the integrator settings can be accessed.
 * Empty if `settings` contains ONLY the integrator settings.
 * \param fallbakInitialTime Initial time to be used if not defined in `settings`.
 * \return Shared pointer to an `IntegratorSettings` object.
 */
template< typename TimeType = double >
boost::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > createIntegratorSettings(
        const json& settings, const KeyTree& keyTree = { }, const TimeType fallbakInitialTime = TUDAT_NAN )
{
    using namespace numerical_integrators;
    using Keys = Keys::Integrator;

    // Read JSON settings shared by all supported integrators
    const auto integratorType = getValue< AvailableIntegrators >( settings, keyTree + Keys::type );
    const auto initialTime = getNumeric( settings, keyTree + Keys::initialTime, fallbakInitialTime );
    const auto initialTimeStep = getNumeric< TimeType >( settings, keyTree + Keys::initialTimeStep );
    const auto saveFrequency = getValuePointer< int >( settings, keyTree + Keys::saveFrequency );

    // Create IntegratorSettings pointer from JSON settings
    switch ( integratorType )
    {
    case euler:
    case rungeKutta4:
    {
        // Construct with mandatory arguments, and optional arguments as default
        IntegratorSettings< TimeType > integratorSettings( integratorType, initialTime, initialTimeStep );

        // Update optional arguments if provided
        if ( saveFrequency )
        {
            integratorSettings.saveFrequency_ = *saveFrequency;
        }

        // Return shared pointer
        return boost::make_shared< IntegratorSettings< TimeType > >( integratorSettings );
    }
    case rungeKuttaVariableStepSize:
    {
        // Read additional settings specific for this integrator
        const auto coefficientSet = getValue< RungeKuttaCoefficients::CoefficientSets >(
                    settings, keyTree + Keys::rungeKuttaCoefficientSet );
        const auto minimumStepSize = getNumeric< TimeType >( settings, keyTree + Keys::minimumStepSize );
        const auto maximumStepSize = getNumeric< TimeType >( settings, keyTree + Keys::maximumStepSize );
        const auto relativeErrorTolerance = getValuePointer< TimeType >(
                    settings, keyTree + Keys::relativeErrorTolerance );
        const auto absoluteErrorTolerance = getValuePointer< TimeType >(
                    settings, keyTree + Keys::absoluteErrorTolerance );
        const auto safetyFactorForNextStepSize = getValuePointer< TimeType >(
                    settings, keyTree + Keys::safetyFactorForNextStepSize );
        const auto maximumFactorIncreaseForNextStepSize = getValuePointer< TimeType >(
                    settings, keyTree + Keys::maximumFactorIncreaseForNextStepSize );
        const auto minimumFactorDecreaseForNextStepSize = getValuePointer< TimeType >(
                    settings, keyTree + Keys::minimumFactorDecreaseForNextStepSize );

        // Construct with mandatory arguments, and optional arguments as default
        RungeKuttaVariableStepSizeSettings< TimeType > integratorSettings(
                    integratorType, initialTime, initialTimeStep, coefficientSet, minimumStepSize, maximumStepSize );

        // Update optional arguments if provided
        if ( relativeErrorTolerance )
        {
            integratorSettings.relativeErrorTolerance_ = *relativeErrorTolerance;
        }
        if ( absoluteErrorTolerance )
        {
            integratorSettings.absoluteErrorTolerance_ = *absoluteErrorTolerance;
        }
        if ( saveFrequency )
        {
            integratorSettings.saveFrequency_ = *saveFrequency;
        }
        if ( safetyFactorForNextStepSize )
        {
            integratorSettings.safetyFactorForNextStepSize_ = *safetyFactorForNextStepSize;
        }
        if ( safetyFactorForNextStepSize )
        {
            integratorSettings.safetyFactorForNextStepSize_ = *safetyFactorForNextStepSize;
        }
        if ( maximumFactorIncreaseForNextStepSize )
        {
            integratorSettings.maximumFactorIncreaseForNextStepSize_ = *maximumFactorIncreaseForNextStepSize;
        }
        if ( minimumFactorDecreaseForNextStepSize )
        {
            integratorSettings.minimumFactorDecreaseForNextStepSize_ = *minimumFactorDecreaseForNextStepSize;
        }

        // Return shared pointer
        return boost::make_shared< RungeKuttaVariableStepSizeSettings< TimeType > >( integratorSettings );
    }
    default:
        throw std::runtime_error( stringFromEnum( integratorType, availableIntegrators )
                                  + " not supported by json_interface." );
    }
}

} // namespace json_interface

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_INTEGRATOR_H
