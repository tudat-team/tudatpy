/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/SimulationSetup/PropagationSetup/thrustSettings.h"

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
                        propulsion::maximum_thrust_multiplier );

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
                           propulsion::maximum_thrust_multiplier );
        int entryForMaximumThrustMultiplierInIndependentVariables =
                std::distance( thrustDependentVariables_.begin( ), findIterator );
        int entryForMaximumThrustMultiplierInGuidanceFunctions = 0;

        for( unsigned int i = 0; i < thrustDependentVariables_.size( ); i++ )
        {
            if( thrustDependentVariables_.at( i ) == propulsion::guidance_input_dependent_thrust )
            {
                entryForMaximumThrustMultiplierInGuidanceFunctions++;
            }
            else if( thrustDependentVariables_.at( i ) == propulsion::maximum_thrust_multiplier)
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
                            propulsion::maximum_thrust_multiplier );

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

} // namespace simulation_setup

} // namespace tudat
