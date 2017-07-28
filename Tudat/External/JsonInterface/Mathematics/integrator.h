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

#include "Tudat/External/JsonInterface/Support/valueAccess.h"
#include "Tudat/External/JsonInterface/Support/valueConversions.h"

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

//! Convert `RungeKuttaCoefficients::CoefficientSets` to `json`.
void to_json( json& jsonObject, const RungeKuttaCoefficients::CoefficientSets& rungeKuttaCoefficientSet );

//! Convert `json` to `RungeKuttaCoefficients::CoefficientSets`.
void from_json( const json& jsonObject, RungeKuttaCoefficients::CoefficientSets& rungeKuttaCoefficientSet );

//! Create a `json` object from a shared pointer to an `IntegratorSettings` object.
template< typename TimeType = double >
void to_json( json& jsonObject, const boost::shared_ptr< IntegratorSettings< TimeType > >& integratorSettings )
{
    if ( integratorSettings )
    {
        using namespace json_interface;
        using K = Keys::Integrator;

        // Common keys
        jsonObject[ K::type ] = stringFromEnum( integratorSettings->integratorType_, availableIntegrators );
        jsonObject[ K::initialTime ] = integratorSettings->initialTime_;
        jsonObject[ K::initialTimeStep ] = integratorSettings->initialTimeStep_;
        jsonObject[ K::saveFrequency ] = integratorSettings->saveFrequency_;

        // Integrator-specific keys
        boost::shared_ptr< RungeKuttaVariableStepSizeSettings< TimeType > > rungeKuttaVariableStepSizeSettings =
                boost::dynamic_pointer_cast< RungeKuttaVariableStepSizeSettings< TimeType > >( integratorSettings );
        if ( rungeKuttaVariableStepSizeSettings )
        {
            jsonObject[ K::rungeKuttaCoefficientSet ] =
                    stringFromEnum( rungeKuttaVariableStepSizeSettings->coefficientSet_, rungeKuttaCoefficientSets );
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
    }
}

//! Create a `json` object from a shared pointer to an `IntegratorSettings` object.
template< typename TimeType = double >
void from_json( const json& jsonObject, boost::shared_ptr< IntegratorSettings< TimeType > >& integratorSettings )
{
    using namespace json_interface;
    using RungeKuttaCoefficientSet = RungeKuttaCoefficients::CoefficientSets;
    using K = Keys::Integrator;

    // Fallback initial time (retrieved from simulation.startEpoch), to be used if not specified in integrator settings
    const TimeType fallbackInitialTime = 0.0; /*getNumeric< TimeType >(
                jsonObject, SpecialKeys::root / KeyPaths::Simulation::startEpoch, TUDAT_NAN, true );*/

    // Read JSON settings shared by all supported integrators
    const AvailableIntegrators integratorType = getValue< AvailableIntegrators >( jsonObject, K::type );
    const TimeType initialTime = getNumeric( jsonObject, K::initialTime, fallbackInitialTime );
    const TimeType initialTimeStep = getNumeric< TimeType >( jsonObject, K::initialTimeStep );

    // Create IntegratorSettings pointer from JSON settings
    switch ( integratorType )
    {
    case euler:
    case rungeKutta4:
    {
        IntegratorSettings< TimeType > defaults( integratorType, 0.0, 0.0 );
        integratorSettings = boost::make_shared< IntegratorSettings< TimeType > >(
                    integratorType, initialTime, initialTimeStep,
                    getNumeric( jsonObject, K::saveFrequency, defaults.saveFrequency_ ) );
        return;
    }
    case rungeKuttaVariableStepSize:
    {
        RungeKuttaVariableStepSizeSettings< TimeType > defaults(
                    integratorType, 0.0, 0.0, RungeKuttaCoefficientSet::rungeKuttaFehlberg45, 0.0, 0.0 );

        RungeKuttaVariableStepSizeSettings< TimeType > rkSettings(
                    integratorType, initialTime, initialTimeStep,
                    getValue< RungeKuttaCoefficientSet >( jsonObject, K::rungeKuttaCoefficientSet ),
                    getNumeric< TimeType >( jsonObject, K::minimumStepSize ),
                    getNumeric< TimeType >( jsonObject, K::maximumStepSize ),
                    getNumeric( jsonObject, K::relativeErrorTolerance, defaults.relativeErrorTolerance_ ),
                    getNumeric( jsonObject, K::absoluteErrorTolerance, defaults.absoluteErrorTolerance_ ),
                    getNumeric( jsonObject, K::saveFrequency, defaults.saveFrequency_ ),
                    getNumeric( jsonObject, K::safetyFactorForNextStepSize,
                                defaults.safetyFactorForNextStepSize_ ),
                    getNumeric( jsonObject, K::maximumFactorIncreaseForNextStepSize,
                                defaults.maximumFactorIncreaseForNextStepSize_ ),
                    getNumeric( jsonObject, K::minimumFactorDecreaseForNextStepSize,
                                defaults.minimumFactorDecreaseForNextStepSize_ ) );

        integratorSettings = boost::make_shared< RungeKuttaVariableStepSizeSettings< TimeType > >( rkSettings );
        return;
    }
    default:
        throw std::runtime_error( stringFromEnum( integratorType, availableIntegrators )
                                  + " not supported by json_interface." );
    }
}

} // namespace numerical_integrators

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_INTEGRATOR_H
