/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/propulsion/thrustMagnitudeWrapper.h"

namespace tudat
{

namespace propulsion
{

    void CustomThrustMagnitudeWrapper::update( const double time )
    {
        if( !( currentTime_ == time ) )
        {
            currentThrustMagnitude_ = thrustMagnitudeFunction_( time );
            currentSpecificImpulse_ = specificImpulseFunction_( time );
            currentTime_ = time;
        }
    }

    void CustomThrustAccelerationMagnitudeWrapper::update( const double time )
    {
        if( !( currentTime_ == time ) )
        {
            currentThrustAccelerationMagnitude_ = thrustAccelerationMagnitudeFunction_( time );
            currentSpecificImpulse_ = specificImpulseFunction_( time );
            currentTime_ = time;
        }
    }


    void MeeCostatesBangBangThrustMagnitudeWrapper::update( const double time )
    {
        if( !( currentTime_ == time ) )
        {

            Eigen::VectorXd costates_ = costateFunction_( time );

            // Get the current state in cartesian coordinates and keplerian elements, and some convenient parameters
            Eigen::Vector6d currentState = thrustingBodyStateFunction_( ) - centralBodyStateFunction_( );
            double centralBodyGravitationalParameter = centralBodyGravitationalParameterFunction_( );

            // Obtain ModifiedEquinoctial elements, flag of 0 indicates that singularity occurs at 180 deg inclination.
            Eigen::Vector6d modifiedEquinoctialElements =
                    orbital_element_conversions::convertCartesianToModifiedEquinoctialElements(
                        currentState, centralBodyGravitationalParameter, 0 );

            // Retrieve modified equinoctial elements.
            double p = modifiedEquinoctialElements[ orbital_element_conversions::semiParameterIndex ];
            double f = modifiedEquinoctialElements[ orbital_element_conversions::fElementIndex ];
            double g = modifiedEquinoctialElements[ orbital_element_conversions::gElementIndex ];
            double h = modifiedEquinoctialElements[ orbital_element_conversions::hElementIndex ];
            double k = modifiedEquinoctialElements[ orbital_element_conversions::kElementIndex ];
            double L = modifiedEquinoctialElements[ orbital_element_conversions::trueLongitudeIndex ];

            double w1 = 1.0 + f * std::cos( L ) + g * std::sin( L );
            double w2 = 1.0 + h * h + k * k;

            // Compute all required auxiliary variables to compute optimal angle alpha.
            double lambdap = costates_[ orbital_element_conversions::semiParameterIndex ] * ( 2.0 * p ) / w1;
            double lambdaf1 = costates_[ orbital_element_conversions::fElementIndex ] * std::sin( L );
            double lambdag1 = costates_[ orbital_element_conversions::gElementIndex ] * std::cos( L );
            double lambdaf2 = costates_[ orbital_element_conversions::fElementIndex ] / w1 *
                    ( ( w1 + 1.0 ) * std::cos( L ) + f );
                    //costates_[ orbital_element_conversions::fElementIndex ] * ( ( w1 + 1.0 ) * std::cos( L ) + f ) / w1;
            double lambdag2 = costates_[ orbital_element_conversions::gElementIndex ] / w1 * ( ( w1 + 1.0 ) * std::sin( L ) + g );

            // Compute sinus of the optimal value of angle alpha.
            double sinOptimalAlpha = - ( lambdaf1 - lambdag1 ) /
                    std::sqrt( ( lambdaf1 - lambdag1 ) * ( lambdaf1 - lambdag1 ) +
                               ( lambdap + lambdaf2 + lambdag2 ) * ( lambdap + lambdaf2 + lambdag2 ) );

            // Compute cosinus of the optimal value of angle alpha.
            double cosOptimalAlpha = - ( lambdap + lambdaf2 + lambdag2 ) /
                    std::sqrt( ( lambdaf1 - lambdag1 ) * ( lambdaf1 - lambdag1 ) +
                               ( lambdap + lambdaf2 + lambdag2 ) * ( lambdap + lambdaf2 + lambdag2 ) );

            // Compute all required auxiliary variables to compute optimal angle beta.
            lambdap = costates_[ orbital_element_conversions::semiParameterIndex ] * ( 2.0 * p ) / w1 * cosOptimalAlpha;
            lambdaf1 = costates_[ orbital_element_conversions::fElementIndex ] * std::sin( L ) * sinOptimalAlpha;
            lambdag1 = costates_[ orbital_element_conversions::gElementIndex ] * std::cos( L ) * sinOptimalAlpha;
            lambdaf2 = costates_[ orbital_element_conversions::fElementIndex ] * ( ( 1.0 + w1 ) * std::cos( L ) + f ) / w1 * cosOptimalAlpha;
            lambdag2 = costates_[ orbital_element_conversions::gElementIndex ] * ( ( 1.0 + w1 ) * std::sin( L ) + g ) / w1 * cosOptimalAlpha;
            double lambdaf3 = costates_[ orbital_element_conversions::fElementIndex ] * ( g / w1 ) * ( h * std::sin( L ) - k * std::cos( L ) );
            double lambdag3 = costates_[ orbital_element_conversions::gElementIndex ] * ( f / w1 ) * ( h * std::sin( L ) - k * std::cos( L ) );
            double lambdah = costates_[ orbital_element_conversions::hElementIndex ] * ( w2 * std::cos( L ) ) / ( 2.0 * w1 );
            double lambdak = costates_[ orbital_element_conversions::kElementIndex ] * ( w2 * std::sin( L ) ) / ( 2.0 * w1 );

            // Compute sinus of optimal thrust angle beta.
            double sinOptimalBeta = - ( - lambdaf3 + lambdag3 + lambdah + lambdak ) /
                    std::sqrt( ( - lambdaf3 + lambdag3 + lambdah + lambdak ) * ( - lambdaf3 + lambdag3 + lambdah + lambdak )
                               + ( lambdap + lambdaf1 - lambdag1 + lambdaf2 + lambdag2 )
                               * ( lambdap + lambdaf1 - lambdag1 + lambdaf2 + lambdag2 ) );

            // Compute cosinus of optimal thrust angle beta.
            double cosOptimalBeta = - ( lambdap + lambdaf1 - lambdag1 + lambdaf2 + lambdag2 ) /
                    std::sqrt( ( - lambdaf3 + lambdag3 + lambdah + lambdak ) * ( - lambdaf3 + lambdag3 + lambdah + lambdak )
                               + ( lambdap + lambdaf1 - lambdag1 + lambdaf2 + lambdag2 )
                               * ( lambdap + lambdaf1 - lambdag1 + lambdaf2 + lambdag2 ) );


            // Switching function for the thrust magnitude.
            double thrustMagnitudeSwitchingCondition = ( 1.0 / thrustingBodyMassFunction_( ) ) *
                    ( lambdap * cosOptimalBeta + lambdah * sinOptimalBeta + lambdak * sinOptimalBeta
                    + lambdaf1 * cosOptimalBeta + lambdaf2 * cosOptimalBeta - lambdaf3 * sinOptimalBeta
                    - lambdag1 * cosOptimalBeta + lambdag2 * cosOptimalBeta + lambdag3 * sinOptimalBeta );


            // Compute current thrust magnitude and specific impulse.
            if ( thrustMagnitudeSwitchingCondition <= 0.0 )
            {
//                std::cout << "INSIDE THRUST MAGNITUDE FUNCTION, THRUST ON. " << "\n\n";
                currentThrustMagnitude_ = maximumThrustMagnitude_;
                currentSpecificImpulse_ = specificImpulseFunction_( time );
            }
            else
            {
//                std::cout << "INSIDE THRUST MAGNITUDE FUNCTION, THRUST OFF. " << "\n\n";
                currentThrustMagnitude_ = 0.0;
                currentSpecificImpulse_ = specificImpulseFunction_( time );
            }

            currentTime_ = time;
        }
    }

    void ParameterizedThrustMagnitudeWrapper::update( const double time )
    {
        if( !( currentTime_ == time ) )
        {
            if( !( inputUpdateFunction_ == nullptr ) )
            {
                inputUpdateFunction_( time );
            }
            // Retrieve thrust independent variables
            for( unsigned int i = 0; i < thrustInputVariableFunctions_.size( ); i++ )
            {
                currentThrustInputVariables_[ i ] = thrustInputVariableFunctions_.at( i )( );
            }

            // Compute thrust
            currentThrustMagnitude_ = thrustMagnitudeFunction_( currentThrustInputVariables_ );

            // Retrieve specific impulse independent variables
            for( unsigned int i = 0; i < specificImpulseInputVariableFunctions_.size( ); i++ )
            {
                currentSpecificImpulseInputVariables_[ i ] = specificImpulseInputVariableFunctions_.at( i )( );
            }

            // Compute specific impulse
            currentSpecificImpulse_ = specificImpulseFunction_( currentSpecificImpulseInputVariables_ );
        }
    }



    void CustomThrustVectorWrapper::update( const double time )
    {
        if( !( currentTime_ == time ) )
        {
            currentThrustVector_ = thrustVectorFunction_( time );
            currentThrustMagnitude_ = currentThrustVector_.norm( );
            currentSpecificImpulse_ = specificImpulseFunction_( time );
            currentTime_ = time;
        }
    }




} // namespace propulsion

} // namespace tudat

