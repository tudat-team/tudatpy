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

#include "utilities.h"

#include <Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h>

namespace tudat
{

namespace json_interface
{

//! -DOC
numerical_integrators::AvailableIntegrators getIntegratorType( const json &settings )
{
    using namespace numerical_integrators;
    std::string stringValue = getValue< std::string >( settings, "type" );
    if ( stringValue == "euler" )
    {
        return euler;
    }
    else if ( stringValue == "rungeKutta4" )
    {
        return rungeKutta4;
    }
    else if ( stringValue == "rungeKuttaVariableStepSize" )
    {
        return rungeKuttaVariableStepSize;
    }
    else
    {
        throw std::runtime_error( "Integrator type \"" + stringValue + "\" not supported." );
    }
}

//! -DOC
numerical_integrators::RungeKuttaCoefficients::CoefficientSets getRungeKuttaCoefficientSet( const json &settings )
{
    using namespace numerical_integrators;
    std::string stringValue = getValue< std::string >( settings, "rungeKuttaCoefficientSet" );
    if ( stringValue == "rungeKuttaFehlberg45" )
    {
        return RungeKuttaCoefficients::rungeKuttaFehlberg45;
    }
    else if ( stringValue == "rungeKuttaFehlberg56" )
    {
        return RungeKuttaCoefficients::rungeKuttaFehlberg56;
    }
    else if ( stringValue == "rungeKuttaFehlberg78" )
    {
        return RungeKuttaCoefficients::rungeKuttaFehlberg78;
    }
    else if ( stringValue == "rungeKutta87DormandPrince" )
    {
        return RungeKuttaCoefficients::rungeKutta87DormandPrince;
    }
    else
    {
        throw std::runtime_error( "Runge Kutta coefficient set \"" + stringValue + "\" not supported." );
    }
}

//! -DOC
template< typename TimeType = double >
boost::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > createIntegratorSettings(
        const json &settings )
{
    using namespace numerical_integrators;

    // Read JSON settings shared by all supported integrators
    const AvailableIntegrators integratorType = getIntegratorType( settings );
    const TimeType initialTime = getValue< TimeType >( settings, "initialTime" );
    const TimeType initialTimeStep = getValue< TimeType >( settings, "initialTimeStep" );
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
        const RungeKuttaCoefficients::CoefficientSets coefficientSet = getRungeKuttaCoefficientSet( settings );
        const TimeType minimumStepSize = getValue< TimeType >( settings, "minimumStepSize" );
        const TimeType maximumStepSize = getValue< TimeType >( settings, "maximumStepSize" );
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


namespace numerical_integrators
{
//! Function to create a `json` object from an `IntegrationSettings` pointer.
//! Called automatically by `nlohmann::json` when writing `json( simulation )`.
template< typename TimeType = double >
void to_json( json& j,
              const boost::shared_ptr< IntegratorSettings< TimeType > >& integratorSettings ) {
    // Initialise
    j = json();

    // Common keys
    j[ "type" ] = integratorSettings->integratorType_;
    j[ "initialTime" ] = integratorSettings->initialTime_;
    j[ "initialTimeStep" ] = integratorSettings->initialTimeStep_;
    j[ "saveFrequency" ] = integratorSettings->saveFrequency_;

    // Integrator-specific keys
    boost::shared_ptr< RungeKuttaVariableStepSizeSettings< TimeType > > rungeKuttaVariableStepSizeSettings =
            boost::dynamic_pointer_cast< RungeKuttaVariableStepSizeSettings< TimeType > >( integratorSettings );
    if ( rungeKuttaVariableStepSizeSettings != NULL )
    {
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


} // namespace tudat

#endif // TUDAT_JSONINTERFACE_INTEGRATORSETTINGS_H
