/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Mathematics/Interpolators/multiLinearInterpolator.h"
#include "Tudat/SimulationSetup/PropagationSetup/thrustSettings.h"
#include "Tudat/InputOutput/multiDimensionalArrayReader.h"

namespace tudat
{

namespace simulation_setup
{

//! Interface function to multiply a maximum thrust by a multiplier to obtain the actual thrust
double multiplyMaximumThrustByScalingFactor(
        const boost::function< double( const std::vector< double >& ) > maximumThrustFunction,
        const boost::function< double( ) > maximumThrustMultiplier,
        const std::vector< double >& maximumThrustIndependentVariables )
{
    return maximumThrustMultiplier( ) * maximumThrustFunction( maximumThrustIndependentVariables );
}

//! Function to check the validity of the input data, and process the maximum thrust multiplier if provided
void ParameterizedThrustMagnitudeSettings::parseInputDataAndCheckConsistency(
        const boost::shared_ptr< interpolators::Interpolator< double, double > > thrustMagnitudeInterpolator,
        const boost::shared_ptr< interpolators::Interpolator< double, double > > specificImpulseInterpolator )
{
    int numberOfUserSpecifiedThrustInputs =
            std::count( thrustDependentVariables_.begin( ), thrustDependentVariables_.end( ),
                        propulsion::guidance_input_dependent_thrust );
    int numberOfMaximumThrustMultipliers =
            std::count( thrustDependentVariables_.begin( ), thrustDependentVariables_.end( ),
                        propulsion::throttle_dependent_thrust );
;
    // Check consistency of user-defined thrust input.
    if( ( numberOfUserSpecifiedThrustInputs + numberOfMaximumThrustMultipliers ) !=
            static_cast< int >( thrustGuidanceInputVariables_.size( ) ) )
    {
        throw std::runtime_error( "Error in parameterized thrust settings, inconsistent number of user-defined input variables for thrust" );
    }
    if( thrustMagnitudeInterpolator->getNumberOfDimensions( ) !=
            ( static_cast< int >( thrustDependentVariables_.size( ) ) - numberOfMaximumThrustMultipliers ) )
    {
        throw std::runtime_error( "Error in parameterized thrust settings, thrust interpolator size has inconsistent size" );
    }
    if( numberOfMaximumThrustMultipliers > 1 )
    {
        throw std::runtime_error( "Error  in parameterized thrust settings, only 1 maximum thrust multiplier may be defined." );
    }

    // Parse maximum thrust multiplier function
    if( numberOfMaximumThrustMultipliers == 1 )
    {
        std::vector< propulsion::ThrustDependentVariables >::iterator findIterator =
                std::find( thrustDependentVariables_.begin( ), thrustDependentVariables_.end( ),
                           propulsion::throttle_dependent_thrust );
        int entryForMaximumThrustMultiplierInIndependentVariables =
                std::distance( thrustDependentVariables_.begin( ), findIterator );
        int entryForMaximumThrustMultiplierInGuidanceFunctions = 0;

        for( unsigned int i = 0; i < thrustDependentVariables_.size( ); i++ )
        {
            if( thrustDependentVariables_.at( i ) == propulsion::guidance_input_dependent_thrust )
            {
                entryForMaximumThrustMultiplierInGuidanceFunctions++;
            }
            else if( thrustDependentVariables_.at( i ) == propulsion::throttle_dependent_thrust)
            {
                break;
            }
        }

        boost::function< double( const std::vector< double >& ) > maximumThrustMagnitudeFunction =
                thrustMagnitudeFunction_;
        thrustMagnitudeFunction_ = boost::bind(
                    &multiplyMaximumThrustByScalingFactor, maximumThrustMagnitudeFunction,
                    thrustGuidanceInputVariables_.at( entryForMaximumThrustMultiplierInGuidanceFunctions ), _1 );
        thrustDependentVariables_.erase( thrustDependentVariables_.begin( ) +
                                         entryForMaximumThrustMultiplierInIndependentVariables);
        thrustGuidanceInputVariables_.erase( thrustGuidanceInputVariables_.begin( ) +
                                             entryForMaximumThrustMultiplierInGuidanceFunctions );

    }

    if( specificImpulseInterpolator != NULL )
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

//! Function to create thrust magnitude settings from guidance input and tables
boost::shared_ptr< ParameterizedThrustMagnitudeSettings > createParameterizedThrustMagnitudeSettings(
        const boost::shared_ptr< ThrustInputParameterGuidance > thrustInputParameterGuidance,
        const boost::shared_ptr< interpolators::Interpolator< double, double > > thrustMagnitudeInterpolator,
        const std::vector< propulsion::ThrustDependentVariables > thrustDependentVariables,
        const boost::shared_ptr< interpolators::Interpolator< double, double > > specificImpulseInterpolator,
        const std::vector< propulsion::ThrustDependentVariables > specificImpulseDependentVariables )
{
    // Create guidance input functions.
    std::vector< boost::function< double( ) > > thrustGuidanceInputVariables;
    for( int i = 0; i < thrustInputParameterGuidance->getNumberOfThrustInputParameters( ); i++ )
    {
        thrustGuidanceInputVariables.push_back(
                    boost::bind( &ThrustInputParameterGuidance::getThrustInputGuidanceParameter,
                                 thrustInputParameterGuidance, i ) );
    }

    // Create specific impulse input functions.
    std::vector< boost::function< double( ) > > specificImpulseGuidanceInputVariables;
    for( int i = 0; i < thrustInputParameterGuidance->getNumberOfSpecificImpulseInputParameters( ); i++ )
    {
        specificImpulseGuidanceInputVariables.push_back(
                    boost::bind( &ThrustInputParameterGuidance::getSpecificImpulseInputGuidanceParameter,
                                 thrustInputParameterGuidance, i ) );
    }

    // Create thrust magnitude settings object.
    return boost::make_shared< ParameterizedThrustMagnitudeSettings >(
                thrustMagnitudeInterpolator,
                thrustDependentVariables,
                specificImpulseInterpolator,
                specificImpulseDependentVariables,
                thrustGuidanceInputVariables,
                specificImpulseGuidanceInputVariables,
                boost::bind( &ThrustInputParameterGuidance::update, thrustInputParameterGuidance, _1 ) );
}

//! Function to create thrust magnitude settings from guidance input and tables
boost::shared_ptr< ParameterizedThrustMagnitudeSettings > createParameterizedThrustMagnitudeSettings(
        const boost::shared_ptr< ThrustInputParameterGuidance > thrustInputParameterGuidance,
        const std::string thrustMagnitudeDataFile,
        const std::vector< propulsion::ThrustDependentVariables > thrustDependentVariables,
        const std::string specificImpulseDataFile,
        const std::vector< propulsion::ThrustDependentVariables > specificImpulseDependentVariables )
{
    // Read thrust files.
    boost::shared_ptr< interpolators::Interpolator< double, double > > thrustMagnitudeInterpolator;
    int numberOfThrustIndependentVariables = input_output::getNumberOfIndependentVariablesInCoefficientFile(
                thrustMagnitudeDataFile );
    if( numberOfThrustIndependentVariables == 1 )
    {
        std::pair< boost::multi_array< double, 1 >, std::vector< std::vector< double > > > thrustData =
                input_output::MultiArrayFileReader< 1 >::readMultiArrayAndIndependentVariables(
                    thrustMagnitudeDataFile );
        thrustMagnitudeInterpolator = boost::make_shared< interpolators::MultiLinearInterpolator< double, double, 1 > >(
                    thrustData.second, thrustData.first );
    }
    else if( numberOfThrustIndependentVariables == 2 )
    {
        std::pair< boost::multi_array< double, 2 >, std::vector< std::vector< double > > > thrustData =
                input_output::MultiArrayFileReader< 2 >::readMultiArrayAndIndependentVariables(
                    thrustMagnitudeDataFile );
        thrustMagnitudeInterpolator = boost::make_shared< interpolators::MultiLinearInterpolator< double, double, 2 > >(
                    thrustData.second, thrustData.first );
    }
    else if( numberOfThrustIndependentVariables == 3 )
    {
        std::pair< boost::multi_array< double, 3 >, std::vector< std::vector< double > > > thrustData =
                input_output::MultiArrayFileReader< 3 >::readMultiArrayAndIndependentVariables(
                    thrustMagnitudeDataFile );
        thrustMagnitudeInterpolator = boost::make_shared< interpolators::MultiLinearInterpolator< double, double, 3 > >(
                    thrustData.second, thrustData.first );
    }
    else
    {
        throw std::runtime_error( "Error, reading of thrust magnitude files of mroe than 3 independent variables not yet implemented" );
    }


    // Read specific impulse files.
    boost::shared_ptr< interpolators::Interpolator< double, double > > specificImpulseInterpolator;
    int numberOfSpecificImpulseIndependentVariables = input_output::getNumberOfIndependentVariablesInCoefficientFile(
                specificImpulseDataFile );
    if( numberOfSpecificImpulseIndependentVariables == 1 )
    {
        std::pair< boost::multi_array< double, 1 >, std::vector< std::vector< double > > > specificImpulseData =
                input_output::MultiArrayFileReader< 1 >::readMultiArrayAndIndependentVariables(
                    specificImpulseDataFile );
        specificImpulseInterpolator = boost::make_shared< interpolators::MultiLinearInterpolator< double, double, 1 > >(
                    specificImpulseData.second, specificImpulseData.first );
    }
    else if( numberOfSpecificImpulseIndependentVariables == 2 )
    {
        std::pair< boost::multi_array< double, 2 >, std::vector< std::vector< double > > > specificImpulseData =
                input_output::MultiArrayFileReader< 2 >::readMultiArrayAndIndependentVariables(
                    specificImpulseDataFile );
        specificImpulseInterpolator = boost::make_shared< interpolators::MultiLinearInterpolator< double, double, 2 > >(
                    specificImpulseData.second, specificImpulseData.first );
    }
    else if( numberOfSpecificImpulseIndependentVariables == 3 )
    {
        std::pair< boost::multi_array< double, 3 >, std::vector< std::vector< double > > > specificImpulseData =
                input_output::MultiArrayFileReader< 3 >::readMultiArrayAndIndependentVariables(
                    specificImpulseDataFile );
        specificImpulseInterpolator = boost::make_shared< interpolators::MultiLinearInterpolator< double, double, 3 > >(
                    specificImpulseData.second, specificImpulseData.first );
    }
    else
    {
        throw std::runtime_error( "Error, reading of specific impulse files of mroe than 3 independent variables not yet implemented" );
    }

    // Create thrust magnitude settings object.
    return createParameterizedThrustMagnitudeSettings(
                thrustInputParameterGuidance, thrustMagnitudeInterpolator, thrustDependentVariables,
                specificImpulseInterpolator, specificImpulseDependentVariables );
}

//! Function to create thrust magnitude settings from guidance input and tables, with constant specific impulse
boost::shared_ptr< ParameterizedThrustMagnitudeSettings > createParameterizedThrustMagnitudeSettings(
        const boost::shared_ptr< ThrustInputParameterGuidance > thrustInputParameterGuidance,
        const boost::shared_ptr< interpolators::Interpolator< double, double > > thrustMagnitudeInterpolator,
        const std::vector< propulsion::ThrustDependentVariables > thrustDependentVariables,
        const double constantSpecificImpulse )
{
    // Create guidance input functions.
    std::vector< boost::function< double( ) > > thrustGuidanceInputVariables;
    for( int i = 0; i < thrustInputParameterGuidance->getNumberOfThrustInputParameters( ); i++ )
    {
        thrustGuidanceInputVariables.push_back(
                    boost::bind( &ThrustInputParameterGuidance::getThrustInputGuidanceParameter,
                                 thrustInputParameterGuidance, i ) );
    }

    // Create thrust magnitude settings object.
    return boost::make_shared< ParameterizedThrustMagnitudeSettings >(
                thrustMagnitudeInterpolator,
                thrustDependentVariables,
                constantSpecificImpulse,
                thrustGuidanceInputVariables,
                boost::bind( &ThrustInputParameterGuidance::update, thrustInputParameterGuidance, _1 ) );
}

//! Function to create thrust magnitude settings from guidance input and tables, with constant specific impulse
boost::shared_ptr< ParameterizedThrustMagnitudeSettings > createParameterizedThrustMagnitudeSettings(
        const boost::shared_ptr< ThrustInputParameterGuidance > thrustInputParameterGuidance,
        const std::string thrustMagnitudeDataFile,
        const std::vector< propulsion::ThrustDependentVariables > thrustDependentVariables,
        const double constantSpecificImpulse )
{
    boost::shared_ptr< interpolators::Interpolator< double, double > > thrustMagnitudeInterpolator;
    int numberOfThrustIndependentVariables = input_output::getNumberOfIndependentVariablesInCoefficientFile(
                thrustMagnitudeDataFile );

    // Read thrust files.
    if( numberOfThrustIndependentVariables == 1 )
    {
        std::pair< boost::multi_array< double, 1 >, std::vector< std::vector< double > > > thrustData =
                input_output::MultiArrayFileReader< 1 >::readMultiArrayAndIndependentVariables(
                    thrustMagnitudeDataFile );
        thrustMagnitudeInterpolator = boost::make_shared< interpolators::MultiLinearInterpolator< double, double, 1 > >(
                    thrustData.second, thrustData.first );
    }
    else if( numberOfThrustIndependentVariables == 2 )
    {
        std::pair< boost::multi_array< double, 2 >, std::vector< std::vector< double > > > thrustData =
                input_output::MultiArrayFileReader< 2 >::readMultiArrayAndIndependentVariables(
                    thrustMagnitudeDataFile );
        thrustMagnitudeInterpolator = boost::make_shared< interpolators::MultiLinearInterpolator< double, double, 2 > >(
                    thrustData.second, thrustData.first );
    }
    else if( numberOfThrustIndependentVariables == 3 )
    {
        std::pair< boost::multi_array< double, 3 >, std::vector< std::vector< double > > > thrustData =
                input_output::MultiArrayFileReader< 3 >::readMultiArrayAndIndependentVariables(
                    thrustMagnitudeDataFile );
        thrustMagnitudeInterpolator = boost::make_shared< interpolators::MultiLinearInterpolator< double, double, 3 > >(
                    thrustData.second, thrustData.first );
    }
    else
    {
        throw std::runtime_error( "Error, reading of thrust magnitude files of mroe than 3 independent variables not yet implemented" );
    }

    // Create thrust magnitude settings object.
    return createParameterizedThrustMagnitudeSettings(
                thrustInputParameterGuidance, thrustMagnitudeInterpolator, thrustDependentVariables,
                constantSpecificImpulse );
}



} // namespace simulation_setup

} // namespace tudat
