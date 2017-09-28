/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_STATE_H
#define TUDAT_JSONINTERFACE_STATE_H

#include <Tudat/Astrodynamics/BasicAstrodynamics/sphericalBodyShapeModel.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>

#include "Tudat/JsonInterface/Support/valueAccess.h"
#include "Tudat/JsonInterface/Support/valueConversions.h"

namespace tudat
{

namespace json_interface
{

// StateType

//! Possible ways of providing initial translational states.
enum StateType
{
    cartesianComponents,
    keplerianComponents
};

//! Map of `StateType` string representations.
static std::map< StateType, std::string > stateTypes =
{
    { cartesianComponents, "cartesian" },
    { keplerianComponents, "keplerian" }
};

//! `StateType` not supported by `json_interface`.
static std::vector< StateType > unsupportedStateTypes = { };

//! Convert `StateType` to `json`.
inline void to_json( nlohmann::json& jsonObject, const StateType& stateType )
{
    jsonObject = json_interface::stringFromEnum( stateType, stateTypes );
}

//! Convert `json` to `StateType`.
inline void from_json( const nlohmann::json& jsonObject, StateType& stateType )
{
    stateType = json_interface::enumFromString( jsonObject, stateTypes );
}

//! -DOC
template< typename StateScalarType >
Eigen::Matrix< StateScalarType, 6, 1 > getCartesianState(
        const nlohmann::json& jsonObject,
        const KeyPath& keyPath = KeyPath( ),
        double centralBodyGravitationalParameter = TUDAT_NAN,
        double centralBodyAverageRadius = TUDAT_NAN )
{
    using namespace orbital_element_conversions;
    using K = Keys::Body::State;

    const nlohmann::json jsonState = getValue< nlohmann::json >( jsonObject, keyPath );
    Eigen::Matrix< StateScalarType, 6, 1 > bodyState = Eigen::Matrix< StateScalarType, 6, 1 >::Zero( );

    // Instead of a vector, an object can be used to provide initial translational state
    if ( jsonState.is_object( ) )
    {
        const StateType stateType = getValue< StateType >( jsonState, K::type );
        switch ( stateType ) {
        case cartesianComponents:
        {
            updateFromJSONIfDefined( bodyState( xCartesianPositionIndex ), jsonState, K::x );
            updateFromJSONIfDefined( bodyState( yCartesianPositionIndex ), jsonState, K::y );
            updateFromJSONIfDefined( bodyState( zCartesianPositionIndex ), jsonState, K::z );
            updateFromJSONIfDefined( bodyState( xCartesianVelocityIndex ), jsonState, K::vx );
            updateFromJSONIfDefined( bodyState( yCartesianVelocityIndex ), jsonState, K::vy );
            updateFromJSONIfDefined( bodyState( zCartesianVelocityIndex ), jsonState, K::vz );
            break;
        }
        case keplerianComponents:
        {
            if ( isNaN( centralBodyGravitationalParameter ) )
            {
                centralBodyGravitationalParameter =
                        getValue< double >( jsonState, K::centralBodyGravitationalParameter );
            }

            // Detemrine semiMajorAxis, eccentricity, argumentOfPeriapsis

            StateScalarType semiMajorAxis = TUDAT_NAN;
            StateScalarType eccentricity = TUDAT_NAN;
            StateScalarType argumentOfPeriapsis = TUDAT_NAN;

            if ( isDefined( jsonState, K::radius ) || isDefined( jsonState, K::altitude ) )  // circular
            {
                if ( isDefined( jsonState, K::altitude ) )
                {
                    if ( isNaN( centralBodyAverageRadius ) )
                    {
                        centralBodyAverageRadius = getValue< double >( jsonState, K::centralBodyAverageRadius );
                    }
                    semiMajorAxis = centralBodyAverageRadius + getValue< StateScalarType >( jsonState, K::altitude );
                }
                else
                {
                    semiMajorAxis = getValue< StateScalarType >( jsonState, K::radius );
                }
                eccentricity = 0.0;
                argumentOfPeriapsis = 0.0;
            }
            else  // generic
            {
                argumentOfPeriapsis = getValue< StateScalarType >( jsonState, K::argumentOfPeriapsis, 0.0 );

                if ( isDefined( jsonState, K::apoapsisDistance ) || isDefined( jsonState, K::apoapsisAltitude ) ||
                     isDefined( jsonState, K::periapsisDistance ) || isDefined( jsonState, K::periapsisAltitude ) )
                {
                    StateScalarType apoapsisDistance = TUDAT_NAN;
                    if ( isDefined( jsonState, K::apoapsisAltitude ) )
                    {
                        if ( isNaN( centralBodyAverageRadius ) )
                        {
                            centralBodyAverageRadius = getValue< double >( jsonState, K::centralBodyAverageRadius );
                        }
                        apoapsisDistance = centralBodyAverageRadius +
                                getValue< StateScalarType >( jsonState, K::apoapsisAltitude );
                    }
                    else if ( isDefined( jsonState, K::apoapsisDistance ) )
                    {
                        apoapsisDistance = getValue< StateScalarType >( jsonState, K::apoapsisDistance );
                    }

                    StateScalarType periapsisDistance = TUDAT_NAN;
                    if ( isDefined( jsonState, K::periapsisAltitude ) )
                    {
                        if ( isNaN( centralBodyAverageRadius ) )
                        {
                            centralBodyAverageRadius = getValue< double >( jsonState, K::centralBodyAverageRadius );
                        }
                        periapsisDistance = centralBodyAverageRadius +
                                getValue< StateScalarType >( jsonState, K::periapsisAltitude );
                    }
                    else if ( isDefined( jsonState, K::periapsisDistance ) )
                    {
                        periapsisDistance = getValue< StateScalarType >( jsonState, K::periapsisDistance );
                    }

                    if ( ! isNaN( apoapsisDistance) && ! isNaN( periapsisDistance ) )
                    {
                        // r_a, r_p -> a, e
                        semiMajorAxis = ( apoapsisDistance + periapsisDistance ) / 2.0;
                        eccentricity = ( apoapsisDistance - periapsisDistance ) / ( 2.0 * semiMajorAxis );
                    }
                    else if ( ! isNaN( apoapsisDistance) )
                    {
                        if ( isDefined( jsonState, K::semiMajorAxis ) )
                        {
                            semiMajorAxis = getValue< StateScalarType >( jsonState, K::semiMajorAxis );
                            // r_a, a -> e
                            eccentricity = apoapsisDistance / semiMajorAxis - 1.0;
                        }
                        else
                        {
                            eccentricity = getValue< StateScalarType >( jsonState, K::eccentricity );
                            // r_a, e -> a
                            semiMajorAxis = apoapsisDistance / ( 1.0 + eccentricity );
                        }
                    }
                    else
                    {
                        if ( isDefined( jsonState, K::semiMajorAxis ) )
                        {
                            semiMajorAxis = getValue< StateScalarType >( jsonState, K::semiMajorAxis );
                            // r_p, a -> e
                            eccentricity = 1.0 - periapsisDistance / semiMajorAxis;
                        }
                        else
                        {
                            eccentricity = getValue< StateScalarType >( jsonState, K::eccentricity );
                            // r_p, e -> a
                            semiMajorAxis = periapsisDistance / ( 1.0 - eccentricity );
                        }
                    }
                }
                else
                {
                    eccentricity = getValue< StateScalarType >( jsonState, K::eccentricity, 0.0 );

                    if ( isDefined( jsonState, K::semiLatusRectum ) )
                    {
                        semiMajorAxis = getValue< StateScalarType >( jsonState, K::semiLatusRectum ) /
                                ( 1.0 - std::pow( eccentricity, 2.0 ) );
                    }
                    else if ( isDefined( jsonState, K::meanMotion ) || isDefined( jsonState, K::period ) )
                    {
                        StateScalarType meanMotion;
                        if ( isDefined( jsonState, K::meanMotion ) )
                        {
                            meanMotion = getValue< StateScalarType >( jsonState, K::meanMotion );
                        }
                        else
                        {
                            meanMotion = 2.0 * mathematical_constants::PI /
                                    getValue< StateScalarType >( jsonState, K::period );
                        }
                        semiMajorAxis = std::pow( ( centralBodyGravitationalParameter /
                                                    std::pow( meanMotion, 2.0 ) ), 1.0/3.0 );
                    }
                    else
                    {
                        semiMajorAxis = getValue< StateScalarType >( jsonState, K::semiMajorAxis );
                    }
                }
            }


            // Determine inclination, longitudeOfAscendingNode

            StateScalarType inclination = getValue< StateScalarType >( jsonState, K::inclination, 0.0 );
            StateScalarType longitudeOfAscendingNode =
                    getValue< StateScalarType >( jsonState, K::longitudeOfAscendingNode, 0.0 );


            // Determine trueAnomaly

            StateScalarType trueAnomaly;
            if ( isDefined( jsonState, K::meanAnomaly ) || isDefined( jsonState, K::eccentricAnomaly ) )
            {
                StateScalarType eccentricAnomaly;
                if ( isDefined( jsonState, K::meanAnomaly ) )
                {
                    eccentricAnomaly = convertMeanAnomalyToEccentricAnomaly(
                                eccentricity,
                                getValue< StateScalarType >( jsonState, K::meanAnomaly ) );
                }
                else
                {
                    eccentricAnomaly = getValue< StateScalarType >( jsonState, K::eccentricAnomaly );
                }
                trueAnomaly = convertEccentricAnomalyToTrueAnomaly( eccentricAnomaly, eccentricity );
            }
            else
            {
                trueAnomaly = getValue< StateScalarType >( jsonState, K::trueAnomaly, 0.0 );
            }


            // Full state

            Eigen::Matrix< StateScalarType, 6, 1 > keplerianElements;
            keplerianElements( semiMajorAxisIndex ) = semiMajorAxis;
            keplerianElements( eccentricityIndex ) = eccentricity;
            keplerianElements( inclinationIndex ) = inclination;
            keplerianElements( argumentOfPeriapsisIndex ) = argumentOfPeriapsis;
            keplerianElements( longitudeOfAscendingNodeIndex ) = longitudeOfAscendingNode;
            keplerianElements( trueAnomalyIndex ) = trueAnomaly;


            // Convert to Cartesian elements

            bodyState = convertKeplerianToCartesianElements( keplerianElements, centralBodyGravitationalParameter );

            break;
        }
        default:
            handleUnimplementedEnumValue( stateType, stateTypes, unsupportedStateTypes );
        }
    }
    else  // could not get the state as Cartesian, Keplerian components... then try to convert directly
    {
        bodyState = getValue< Eigen::Matrix< StateScalarType, 6, 1 > >( jsonObject, keyPath );
    }

    return bodyState;
}

} // namespace json_interface

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_STATE_H
