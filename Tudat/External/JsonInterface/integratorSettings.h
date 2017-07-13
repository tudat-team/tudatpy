/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_INTEGRATORSETTINGS_H
#define TUDAT_JSONINTERFACE_INTEGRATORSETTINGS_H

#include "jsonInterface.h"

#include <Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h>

namespace tudat
{

namespace numerical_integrators
{

//! Map of `AvailableIntegrators` supported by `json_interface`.
std::map< std::string, AvailableIntegrators > availableIntegrators =
{
    { "rungeKutta4",                rungeKutta4 },
    { "euler",                      euler },
    { "rungeKuttaVariableStepSize", rungeKuttaVariableStepSize }
};


//! Map of `RungeKuttaCoefficients::CoefficientSets` supported by `json_interface`.
std::map< std::string, RungeKuttaCoefficients::CoefficientSets > rungeKuttaCoefficientSets =
{
    { "rungeKuttaFehlberg45",      RungeKuttaCoefficients::rungeKuttaFehlberg45 },
    { "rungeKuttaFehlberg56",      RungeKuttaCoefficients::rungeKuttaFehlberg56 },
    { "rungeKuttaFehlberg78",      RungeKuttaCoefficients::rungeKuttaFehlberg78 },
    { "rungeKutta87DormandPrince", RungeKuttaCoefficients::rungeKutta87DormandPrince }
};


//! Create a `json` object from a shared pointer to an `IntegratorSettings` object.
//! Called automatically by `nlohmann::json` when using a constructor such as `json( integratorSettings )`.
template< typename TimeType = double >
void to_json( json& j, const boost::shared_ptr< IntegratorSettings< TimeType > >& integratorSettings ) {
    using namespace json_interface;
    // Initialise
    j = json();

    // Common keys
    j[ "type" ] = stringFromEnum( integratorSettings->integratorType_, availableIntegrators );
    j[ "initialTime" ] = integratorSettings->initialTime_;
    j[ "initialTimeStep" ] = integratorSettings->initialTimeStep_;
    j[ "saveFrequency" ] = integratorSettings->saveFrequency_;

    // Integrator-specific keys
    boost::shared_ptr< RungeKuttaVariableStepSizeSettings< TimeType > > rungeKuttaVariableStepSizeSettings =
            boost::dynamic_pointer_cast< RungeKuttaVariableStepSizeSettings< TimeType > >( integratorSettings );
    if ( rungeKuttaVariableStepSizeSettings != NULL )
    {
        j[ "rungeKuttaCoefficientSet" ] =
                stringFromEnum( rungeKuttaVariableStepSizeSettings->coefficientSet_, rungeKuttaCoefficientSets );
        j[ "minimumStepSize" ] = rungeKuttaVariableStepSizeSettings->minimumStepSize_;
        j[ "maximumStepSize" ] = rungeKuttaVariableStepSizeSettings->maximumStepSize_;
        j[ "relativeErrorTolerance" ] = rungeKuttaVariableStepSizeSettings->relativeErrorTolerance_;
        j[ "absoluteErrorTolerance" ] = rungeKuttaVariableStepSizeSettings->absoluteErrorTolerance_;
        j[ "safetyFactorForNextStepSize" ] = rungeKuttaVariableStepSizeSettings->safetyFactorForNextStepSize_;
        j[ "maximumFactorIncreaseForNextStepSize" ] =
                rungeKuttaVariableStepSizeSettings->maximumFactorIncreaseForNextStepSize_;
        j[ "minimumFactorDecreaseForNextStepSize" ] =
                rungeKuttaVariableStepSizeSettings->minimumFactorDecreaseForNextStepSize_;
    }
}

} // namespace numerical_integrators


namespace json_interface
{

//! Create a shared pointer to an `IntegratorSettings` object from a `json` object.
/*!
 * Create a shared pointer to an `IntegratorSettings` object from a `json` object.
 * \param settings `json` object containing only the settings for one integrator.
 * \return Shared pointer to an `IntegratorSettings` object.
 */
template< typename TimeType = double >
boost::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > createIntegratorSettings(
        const json &settings )
{
    using namespace numerical_integrators;

    // Read JSON settings shared by all supported integrators
    const AvailableIntegrators integratorType =
            enumFromString( getValue< std::string >( settings, "type" ), availableIntegrators );
    const TimeType initialTime = getNumber< TimeType >( settings, "initialTime" );
    const TimeType initialTimeStep = getNumber< TimeType >( settings, "initialTimeStep" );
    boost::shared_ptr< int > saveFrequency = getValuePointer< int >( settings, "saveFrequency" );

    // Create IntegratorSettings pointer from JSON settings
    if ( integratorType == euler || integratorType == rungeKutta4 )
    {
        // Construct with mandatory arguments, and optional arguments as default
        boost::shared_ptr< IntegratorSettings< TimeType > > integratorSettings =
                boost::make_shared< IntegratorSettings< TimeType > >( integratorType, initialTime, initialTimeStep );

        // Update optional arguments if provided
        if ( saveFrequency != NULL )
        {
            integratorSettings->saveFrequency_ = *saveFrequency;
        }
        return integratorSettings;
    }
    else
    {
        // Read additional settings specific for this integrator
        const RungeKuttaCoefficients::CoefficientSets coefficientSet = enumFromString(
                    getValue< std::string >( settings, "rungeKuttaCoefficientSet" ), rungeKuttaCoefficientSets );
        const TimeType minimumStepSize = getNumber< TimeType >( settings, "minimumStepSize" );
        const TimeType maximumStepSize = getNumber< TimeType >( settings, "maximumStepSize" );
        boost::shared_ptr< TimeType > relativeErrorTolerance =
                getValuePointer< TimeType >( settings, "relativeErrorTolerance" );
        boost::shared_ptr< TimeType > absoluteErrorTolerance =
                getValuePointer< TimeType >( settings, "absoluteErrorTolerance" );
        boost::shared_ptr< TimeType > safetyFactorForNextStepSize =
                getValuePointer< TimeType >( settings, "safetyFactorForNextStepSize" );
        boost::shared_ptr< TimeType > maximumFactorIncreaseForNextStepSize =
                getValuePointer< TimeType >( settings, "maximumFactorIncreaseForNextStepSize" );
        boost::shared_ptr< TimeType > minimumFactorDecreaseForNextStepSize =
                getValuePointer< TimeType >( settings, "minimumFactorDecreaseForNextStepSize" );

        // Construct with mandatory arguments, and optional arguments as default
        boost::shared_ptr< RungeKuttaVariableStepSizeSettings< TimeType > > integratorSettings =
                boost::make_shared< RungeKuttaVariableStepSizeSettings< TimeType > >(
                    integratorType, initialTime, initialTimeStep, coefficientSet, minimumStepSize, maximumStepSize );

        // Update optional arguments if provided
        if ( relativeErrorTolerance != NULL )
        {
            integratorSettings->relativeErrorTolerance_ = *relativeErrorTolerance;
        }
        if ( absoluteErrorTolerance != NULL )
        {
            integratorSettings->absoluteErrorTolerance_ = *absoluteErrorTolerance;
        }
        if ( saveFrequency != NULL )
        {
            integratorSettings->saveFrequency_ = *saveFrequency;
        }
        if ( safetyFactorForNextStepSize != NULL )
        {
            integratorSettings->safetyFactorForNextStepSize_ = *safetyFactorForNextStepSize;
        }
        if ( safetyFactorForNextStepSize != NULL )
        {
            integratorSettings->safetyFactorForNextStepSize_ = *safetyFactorForNextStepSize;
        }
        if ( maximumFactorIncreaseForNextStepSize != NULL )
        {
            integratorSettings->maximumFactorIncreaseForNextStepSize_ = *maximumFactorIncreaseForNextStepSize;
        }
        if ( minimumFactorDecreaseForNextStepSize != NULL )
        {
            integratorSettings->minimumFactorDecreaseForNextStepSize_ = *minimumFactorDecreaseForNextStepSize;
        }
        return integratorSettings;
    }
}

} // namespace json_interfaces


} // namespace tudat

#endif // TUDAT_JSONINTERFACE_INTEGRATORSETTINGS_H
