/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/math/interpolators/multiLinearInterpolator.h"
#include "tudat/simulation/environment_setup/thrustSettings.h"
#include "tudat/io/multiDimensionalArrayReader.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to create a list of functions that (compute and) return independent variables for thrust
std::vector< std::function< double( ) > > getPropulsionInputVariables(
        const std::shared_ptr< Body > bodyWithGuidance,
        const std::vector< propulsion::ThrustIndependentVariables > independentVariables,
        const std::vector< std::function< double( ) > > guidanceInputFunctions )
{
    std::vector< std::function< double( ) > > inputFunctions;
    std::shared_ptr< aerodynamics::FlightConditions > vehicleFlightConditions  =
            bodyWithGuidance->getFlightConditions( );

    std::shared_ptr< aerodynamics::AtmosphericFlightConditions > vehicleAtmosphericFlightConditions  =
             std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                bodyWithGuidance->getFlightConditions( ) );

    // Iterate over all dependent variables and create requested function.
    unsigned int numberOfCustomInputs = 0;
    for( unsigned int i = 0; i < independentVariables.size( ); i++ )
    {
        switch( independentVariables.at( i ) )
        {
        case propulsion::altitude_dependent_thrust:
            inputFunctions.push_back(
                        std::bind( &aerodynamics::FlightConditions::getCurrentAltitude, vehicleFlightConditions ) );
            break;
        case propulsion::density_dependent_thrust:
            if( vehicleAtmosphericFlightConditions == nullptr )
            {
                throw std::runtime_error( "Error when making propulsion input variables, atmospheric flight conditions required, but not found" );
            }
            inputFunctions.push_back(
                        std::bind( &aerodynamics::AtmosphericFlightConditions::getCurrentDensity, vehicleAtmosphericFlightConditions ) );
            break;
        case propulsion::dynamic_pressure_dependent_thrust:
            if( vehicleAtmosphericFlightConditions == nullptr )
            {
                throw std::runtime_error( "Error when making propulsion input variables, atmospheric flight conditions required, but not found" );
            }
            inputFunctions.push_back(
                        std::bind( &aerodynamics::AtmosphericFlightConditions::getCurrentDynamicPressure, vehicleAtmosphericFlightConditions ) );
            break;
        case propulsion::mach_number_dependent_thrust:
            if( vehicleAtmosphericFlightConditions == nullptr )
            {
                throw std::runtime_error( "Error when making propulsion input variables, atmospheric flight conditions required, but not found" );
            }
            inputFunctions.push_back(
                        std::bind( &aerodynamics::AtmosphericFlightConditions::getCurrentMachNumber, vehicleAtmosphericFlightConditions ) );
            break;
        case propulsion::pressure_dependent_thrust:
            if( vehicleAtmosphericFlightConditions == nullptr )
            {
                throw std::runtime_error( "Error when making propulsion input variables, atmospheric flight conditions required, but not found" );
            }
            inputFunctions.push_back(
                        std::bind( &aerodynamics::AtmosphericFlightConditions::getCurrentPressure, vehicleAtmosphericFlightConditions ) );
            break;
        case propulsion::guidance_input_dependent_thrust:
            inputFunctions.push_back( guidanceInputFunctions.at( numberOfCustomInputs ) );
            numberOfCustomInputs++;
            break;
        case propulsion::throttle_dependent_thrust:
            if( guidanceInputFunctions.size( ) >= numberOfCustomInputs )
            {
                throw std::runtime_error( "Error when creating propulsion indput dependent variables, insufficient user-defined inputs found" );
            }
            inputFunctions.push_back( guidanceInputFunctions.at( numberOfCustomInputs ) );
            numberOfCustomInputs++;
            break;
        default:
            throw std::runtime_error( "Error when getting parameterized thrust input variables, variable " +
                                      std::to_string( independentVariables.at( i ) ) + "not found" );
        }
    }

    // Check input consistency
    if( numberOfCustomInputs != guidanceInputFunctions.size( ) )
    {
        std::cerr << "Warning when creating propulsion indput dependent variables, not all user-defined inputs have been parsed" << std::endl;
    }

    return inputFunctions;
}


//! Interface function to multiply a maximum thrust by a multiplier to obtain the actual thrust
double multiplyMaximumThrustByScalingFactor(
        const std::function< double( const std::vector< double >& ) > maximumThrustFunction,
        const std::function< double( ) > maximumThrustMultiplier,
        const std::vector< double >& maximumThrustIndependentVariables )
{
    return maximumThrustMultiplier( ) * maximumThrustFunction( maximumThrustIndependentVariables );
}

//! Function to check the validity of the input data, and process the maximum thrust multiplier if provided
void ParameterizedThrustMagnitudeSettings::parseInputDataAndCheckConsistency(
        const std::shared_ptr< interpolators::Interpolator< double, double > > thrustMagnitudeInterpolator,
        const std::shared_ptr< interpolators::Interpolator< double, double > > specificImpulseInterpolator )
{
    int numberOfUserSpecifiedThrustInputs =
            std::count( thrustIndependentVariables_.begin( ), thrustIndependentVariables_.end( ),
                        propulsion::guidance_input_dependent_thrust );
    int numberOfMaximumThrustMultipliers =
            std::count( thrustIndependentVariables_.begin( ), thrustIndependentVariables_.end( ),
                        propulsion::throttle_dependent_thrust );

    // Check consistency of user-defined thrust input.
    if( ( numberOfUserSpecifiedThrustInputs + numberOfMaximumThrustMultipliers ) !=
            static_cast< int >( thrustGuidanceInputVariables_.size( ) ) )
    {
        throw std::runtime_error( "Error in parameterized thrust settings, inconsistent number of user-defined input variables for thrust" );
    }

    if( thrustMagnitudeInterpolator->getNumberOfDimensions( ) !=
            ( static_cast< int >( thrustIndependentVariables_.size( ) ) - numberOfMaximumThrustMultipliers ) )
    {
        throw std::runtime_error( "Error in parameterized thrust settings, thrust interpolator size has inconsistent size" );
    }

    if( numberOfMaximumThrustMultipliers > 1 )
    {
        throw std::runtime_error( "Error in parameterized thrust settings, only 1 maximum thrust multiplier may be defined." );
    }

    // Parse maximum thrust multiplier function
    if( numberOfMaximumThrustMultipliers == 1 )
    {
        std::vector< propulsion::ThrustIndependentVariables >::iterator findIterator =
                std::find( thrustIndependentVariables_.begin( ), thrustIndependentVariables_.end( ),
                           propulsion::throttle_dependent_thrust );
        int entryForMaximumThrustMultiplierInIndependentVariables =
                std::distance( thrustIndependentVariables_.begin( ), findIterator );
        int entryForMaximumThrustMultiplierInGuidanceFunctions = 0;

        for( unsigned int i = 0; i < thrustIndependentVariables_.size( ); i++ )
        {
            if( thrustIndependentVariables_.at( i ) == propulsion::guidance_input_dependent_thrust )
            {
                entryForMaximumThrustMultiplierInGuidanceFunctions++;
            }
            else if( thrustIndependentVariables_.at( i ) == propulsion::throttle_dependent_thrust)
            {
                break;
            }
        }

        std::function< double( const std::vector< double >& ) > maximumThrustMagnitudeFunction =
                thrustMagnitudeFunction_;
        thrustMagnitudeFunction_ = std::bind(
                    &multiplyMaximumThrustByScalingFactor, maximumThrustMagnitudeFunction,
                    thrustGuidanceInputVariables_.at( entryForMaximumThrustMultiplierInGuidanceFunctions ), std::placeholders::_1 );
        thrustIndependentVariables_.erase( thrustIndependentVariables_.begin( ) +
                                         entryForMaximumThrustMultiplierInIndependentVariables);
        thrustGuidanceInputVariables_.erase( thrustGuidanceInputVariables_.begin( ) +
                                             entryForMaximumThrustMultiplierInGuidanceFunctions );

    }

    if( specificImpulseInterpolator != nullptr )
    {
        // Check consistency of user-defined specific impulse input.
        int numberOfUserSpecifiedSpecificImpulseInputs =
                std::count( specificImpulseDependentVariables_.begin( ), specificImpulseDependentVariables_.end( ),
                            propulsion::guidance_input_dependent_thrust );
        int numberOfMaximumThrustMultipliersForSpecificImpulse =
                std::count( specificImpulseDependentVariables_.begin( ), specificImpulseDependentVariables_.end( ),
                            propulsion::throttle_dependent_thrust );

        if( numberOfUserSpecifiedSpecificImpulseInputs !=
                static_cast< int >( specificImpulseGuidanceInputVariables_.size( ) ) )
        {
            throw std::runtime_error( "Error in parameterized thrust settings, inconsistent number of user-defined input variables for specific impulse " );
        }
        if( specificImpulseInterpolator->getNumberOfDimensions( ) !=  static_cast< int >(
                    specificImpulseDependentVariables_.size( ) ) )
        {
            throw std::runtime_error( "Error in parameterized thrust settings, specific impulse interpolator size has inconsistent size" );
        }
        if( numberOfMaximumThrustMultipliersForSpecificImpulse != 0 )
        {
            throw std::runtime_error( "Error in parameterized thrust settings, no maximum thrust multiplier may be defined for specific impulse." );
        }
    }
    else
    {
        if( specificImpulseDependentVariables_.size( ) > 0 || specificImpulseGuidanceInputVariables_.size( ) > 0 )
        {
            throw std::runtime_error( "Error in parameterized thrust settings, no specific impulse interpolator define, but inptu data has been provided." );
        }
    }
}

//! Function to read a thrust or specific impulse interpolator from a file.
std::shared_ptr< interpolators::Interpolator< double, double > > readCoefficientInterpolatorFromFile(
        const std::string coefficientFile )
{
    std::shared_ptr< interpolators::Interpolator< double, double > > coefficientInterpolator;
    int numberOfThrustIndependentVariables = input_output::getNumberOfIndependentVariablesInCoefficientFile(
                coefficientFile );
    if( numberOfThrustIndependentVariables == 1 )
    {
        std::pair< boost::multi_array< double, 1 >, std::vector< std::vector< double > > > coefficientData =
                input_output::MultiArrayFileReader< 1 >::readMultiArrayAndIndependentVariables(
                    coefficientFile );
        coefficientInterpolator = std::make_shared< interpolators::MultiLinearInterpolator< double, double, 1 > >(
                    coefficientData.second, coefficientData.first );
    }
    else if( numberOfThrustIndependentVariables == 2 )
    {
        std::pair< boost::multi_array< double, 2 >, std::vector< std::vector< double > > > coefficientData =
                input_output::MultiArrayFileReader< 2 >::readMultiArrayAndIndependentVariables(
                    coefficientFile );
        coefficientInterpolator = std::make_shared< interpolators::MultiLinearInterpolator< double, double, 2 > >(
                    coefficientData.second, coefficientData.first );
    }
    else if( numberOfThrustIndependentVariables == 3 )
    {
        std::pair< boost::multi_array< double, 3 >, std::vector< std::vector< double > > > coefficientData =
                input_output::MultiArrayFileReader< 3 >::readMultiArrayAndIndependentVariables(
                    coefficientFile );
        coefficientInterpolator = std::make_shared< interpolators::MultiLinearInterpolator< double, double, 3 > >(
                    coefficientData.second, coefficientData.first );
    }
    else
    {
        throw std::runtime_error( "Error, reading of thrust magnitude files of mroe than 3 independent variables not yet implemented" );
    }
    return coefficientInterpolator;
}

//! Function to create thrust magnitude settings from guidance input and tables
std::shared_ptr< ParameterizedThrustMagnitudeSettings > createParameterizedThrustMagnitudeSettings(
        const std::shared_ptr< ThrustInputParameterGuidance > thrustInputParameterGuidance,
        const std::shared_ptr< interpolators::Interpolator< double, double > > thrustMagnitudeInterpolator,
        const std::vector< propulsion::ThrustIndependentVariables > thrustIndependentVariables,
        const std::shared_ptr< interpolators::Interpolator< double, double > > specificImpulseInterpolator,
        const std::vector< propulsion::ThrustIndependentVariables > specificImpulseDependentVariables )
{
    // Create guidance input functions.
    std::vector< std::function< double( ) > > thrustGuidanceInputVariables;
    for( int i = 0; i < thrustInputParameterGuidance->getNumberOfThrustInputParameters( ); i++ )
    {
        thrustGuidanceInputVariables.push_back(
                    std::bind( &ThrustInputParameterGuidance::getThrustInputGuidanceParameter,
                                 thrustInputParameterGuidance, i ) );
    }

    // Create specific impulse input functions.
    std::vector< std::function< double( ) > > specificImpulseGuidanceInputVariables;
    for( int i = 0; i < thrustInputParameterGuidance->getNumberOfSpecificImpulseInputParameters( ); i++ )
    {
        specificImpulseGuidanceInputVariables.push_back(
                    std::bind( &ThrustInputParameterGuidance::getSpecificImpulseInputGuidanceParameter,
                                 thrustInputParameterGuidance, i ) );
    }

    // Create thrust magnitude settings object.
    return std::make_shared< ParameterizedThrustMagnitudeSettings >(
                thrustMagnitudeInterpolator,
                thrustIndependentVariables,
                specificImpulseInterpolator,
                specificImpulseDependentVariables,
                thrustGuidanceInputVariables,
                specificImpulseGuidanceInputVariables,
                std::bind( &ThrustInputParameterGuidance::update, thrustInputParameterGuidance, std::placeholders::_1 ) );
}

//! Function to create thrust magnitude settings from guidance input and tables
std::shared_ptr< ParameterizedThrustMagnitudeSettings > createParameterizedThrustMagnitudeSettings(
        const std::shared_ptr< ThrustInputParameterGuidance > thrustInputParameterGuidance,
        const std::string thrustMagnitudeDataFile,
        const std::vector< propulsion::ThrustIndependentVariables > thrustIndependentVariables,
        const std::string specificImpulseDataFile,
        const std::vector< propulsion::ThrustIndependentVariables > specificImpulseDependentVariables )
{
    // Create thrust magnitude settings object.
    return createParameterizedThrustMagnitudeSettings(
                thrustInputParameterGuidance, readCoefficientInterpolatorFromFile( thrustMagnitudeDataFile) ,
                thrustIndependentVariables, readCoefficientInterpolatorFromFile( specificImpulseDataFile ),
                specificImpulseDependentVariables );
}

//! Function to create thrust magnitude settings from guidance input and tables, with constant specific impulse
std::shared_ptr< ParameterizedThrustMagnitudeSettings > createParameterizedThrustMagnitudeSettings(
        const std::shared_ptr< ThrustInputParameterGuidance > thrustInputParameterGuidance,
        const std::shared_ptr< interpolators::Interpolator< double, double > > thrustMagnitudeInterpolator,
        const std::vector< propulsion::ThrustIndependentVariables > thrustIndependentVariables,
        const double constantSpecificImpulse )
{
    // Create guidance input functions.
    std::vector< std::function< double( ) > > thrustGuidanceInputVariables;
    for( int i = 0; i < thrustInputParameterGuidance->getNumberOfThrustInputParameters( ); i++ )
    {
        thrustGuidanceInputVariables.push_back(
                    std::bind( &ThrustInputParameterGuidance::getThrustInputGuidanceParameter,
                                 thrustInputParameterGuidance, i ) );
    }

    // Create thrust magnitude settings object.
    return std::make_shared< ParameterizedThrustMagnitudeSettings >(
                thrustMagnitudeInterpolator,
                thrustIndependentVariables,
                constantSpecificImpulse,
                thrustGuidanceInputVariables,
                std::bind( &ThrustInputParameterGuidance::update, thrustInputParameterGuidance, std::placeholders::_1 ) );
}

//! Function to create thrust magnitude settings from guidance input and tables, with constant specific impulse
std::shared_ptr< ParameterizedThrustMagnitudeSettings > createParameterizedThrustMagnitudeSettings(
        const std::shared_ptr< ThrustInputParameterGuidance > thrustInputParameterGuidance,
        const std::string thrustMagnitudeDataFile,
        const std::vector< propulsion::ThrustIndependentVariables > thrustIndependentVariables,
        const double constantSpecificImpulse )
{
    // Create thrust magnitude settings object.
    return createParameterizedThrustMagnitudeSettings(
                thrustInputParameterGuidance, readCoefficientInterpolatorFromFile( thrustMagnitudeDataFile ),
                thrustIndependentVariables,
                constantSpecificImpulse );
}

//! Function to create a thrust magnitude settings based on interpolated maximum thrust, with throttle determined by
//! maximum allowed axial acceleration (constant specific impulse).
std::shared_ptr< ParameterizedThrustMagnitudeSettings > createAccelerationLimitedParameterizedThrustMagnitudeSettings(
        const SystemOfBodies& bodies,
        const std::string nameOfBodyWithGuidance,
        const double maximumAcceleration,
        const std::shared_ptr< interpolators::Interpolator< double, double > > thrustMagnitudeInterpolator,
        const std::vector< propulsion::ThrustIndependentVariables > thrustIndependentVariables,
        const double specificImpulse,
        const std::string nameOfCentralBody )
{

    std::shared_ptr< ThrustInputParameterGuidance > thrustGuidance =
            std::make_shared< AccelerationLimitedThrottleGuidance >(
                bodies, nameOfBodyWithGuidance, nameOfCentralBody, thrustIndependentVariables,
               thrustMagnitudeInterpolator, maximumAcceleration );
    return createParameterizedThrustMagnitudeSettings(
                thrustGuidance, thrustMagnitudeInterpolator, thrustIndependentVariables, specificImpulse );
}

//! Function to create a thrust magnitude settings based on interpolated maximum thrust, with throttle determined by
//! maximum allowed axial acceleration (constant specific impulse).
std::shared_ptr< ParameterizedThrustMagnitudeSettings > createAccelerationLimitedParameterizedThrustMagnitudeSettings(
        const SystemOfBodies& bodies,
        const std::string nameOfBodyWithGuidance,
        const double maximumAcceleration,
        const std::string thrustMagnitudeDataFile,
        const std::vector< propulsion::ThrustIndependentVariables > thrustIndependentVariables,
        const double specificImpulse,
        const std::string nameOfCentralBody )
{
    return createAccelerationLimitedParameterizedThrustMagnitudeSettings(
                bodies, nameOfBodyWithGuidance, maximumAcceleration,
                readCoefficientInterpolatorFromFile( thrustMagnitudeDataFile ),
                thrustIndependentVariables, specificImpulse, nameOfCentralBody );
}

//! Function to create a thrust magnitude settings based on interpolated maximum thrust, with throttle determined by
//! maximum allowed axial acceleration.
std::shared_ptr< ParameterizedThrustMagnitudeSettings > createAccelerationLimitedParameterizedThrustMagnitudeSettings(
        const SystemOfBodies& bodies,
        const std::string nameOfBodyWithGuidance,
        const double maximumAcceleration,
        const std::shared_ptr< interpolators::Interpolator< double, double > > thrustMagnitudeInterpolator,
        const std::vector< propulsion::ThrustIndependentVariables > thrustIndependentVariables,
        const std::shared_ptr< interpolators::Interpolator< double, double > > specificImpulseInterpolator,
        const std::vector< propulsion::ThrustIndependentVariables > specificImpulseDependentVariables,
        const std::string nameOfCentralBody )
{

    std::shared_ptr< ThrustInputParameterGuidance > thrustGuidance =
            std::make_shared< AccelerationLimitedThrottleGuidance >(
                bodies, nameOfBodyWithGuidance, nameOfCentralBody, thrustIndependentVariables,
               thrustMagnitudeInterpolator, maximumAcceleration );
    return createParameterizedThrustMagnitudeSettings(
                thrustGuidance, thrustMagnitudeInterpolator, thrustIndependentVariables,
                specificImpulseInterpolator, specificImpulseDependentVariables );
}

//! Function to create a thrust magnitude settings based on interpolated maximum thrust, with throttle determined by
//! maximum allowed axial acceleration.
std::shared_ptr< ParameterizedThrustMagnitudeSettings > createAccelerationLimitedParameterizedThrustMagnitudeSettings(
        const SystemOfBodies& bodies,
        const std::string nameOfBodyWithGuidance,
        const double maximumAcceleration,
        const std::string thrustMagnitudeDataFile,
        const std::vector< propulsion::ThrustIndependentVariables > thrustIndependentVariables,
        const std::string specificImpulseDataFile,
        const std::vector< propulsion::ThrustIndependentVariables > specificImpulseDependentVariables,
        const std::string nameOfCentralBody  )
{
    return createAccelerationLimitedParameterizedThrustMagnitudeSettings(
                bodies, nameOfBodyWithGuidance, maximumAcceleration,
                readCoefficientInterpolatorFromFile( thrustMagnitudeDataFile ),
                thrustIndependentVariables,
                readCoefficientInterpolatorFromFile( specificImpulseDataFile ),
                specificImpulseDependentVariables, nameOfCentralBody );
}



} // namespace simulation_setup

} // namespace tudat
