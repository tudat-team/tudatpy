/*    Copyright (c) 2010-2018, Delft University of Technology
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

#include "Tudat/SimulationSetup/EnvironmentSetup/createBodies.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/sphericalBodyShapeModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/sphericalStateConversions.h"
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
    cartesianElements,
    keplerianElements,
    sphericalElements
};

//! Map of `StateType` string representations.
static std::map< StateType, std::string > stateTypes =
{
    { cartesianElements, "cartesian" },
    { keplerianElements, "keplerian" },
    { sphericalElements, "spherical" }
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

//! Get a Cartesian state from a JSON object.
/*!
 * Get a Cartesian state from a JSON object. The \p jsonObject can be either an array of 6 numbers (in which case it
 * is assumed to be directly the Cartesian state) or an object with different keys (such as vx, eccentricity, etc.)
 * and different possible types (Cartesian, Keplerian, spherical).
 * \param jsonObject The JSON object.
 * \param keyPath The key path at which the state array/object can be retrieved from \p jsonObject.
 * Empty if \p jsonObject is already the array/object (default value).
 * \param centralBody Associated central body from which the central gravitational parameter, average radius, ephemeris,
 * etc. are retrieved when converting certain elements types to Cartesian. This information can also be provided in the
 * \p jsonObject, in which case this can be set to `nullptr` (default value).
 * \param epoch The associated epoch for converting from spherical to Cartesian elements. This information can also be
 * provided in the \p jsonObject, in which case this can be set to `TUDAT_NAN` (default value).
 */
template< typename StateScalarType >
Eigen::Matrix< StateScalarType, 6, 1 > getCartesianState(
        const nlohmann::json& jsonObject,
        const KeyPath& keyPath = KeyPath( ),
        const std::shared_ptr< simulation_setup::Body >& centralBody = nullptr,
        double epoch = TUDAT_NAN )
{
    using namespace basic_astrodynamics;
    using namespace orbital_element_conversions;
    using K = Keys::Body::State;

    const nlohmann::json jsonState = getValue< nlohmann::json >( jsonObject, keyPath );
    Eigen::Matrix< StateScalarType, 6, 1 > bodyState = Eigen::Matrix< StateScalarType, 6, 1 >::Zero( );

    // Instead of a vector, an object can be used to provide initial translational state
    if ( jsonState.is_object( ) )
    {
        const double centralBodyAverageRadius = ! isDefined( jsonState, K::centralBodyAverageRadius ) && centralBody ?
                    centralBody->getShapeModel( )->getAverageRadius( ) :
                    getValue< double >( jsonState, K::centralBodyAverageRadius, TUDAT_NAN );
        bool usedCentralBodyAverageRadius = false;

        const StateType stateType = getValue< StateType >( jsonState, K::type );
        switch ( stateType ) {
        case cartesianElements:
        {
            updateFromJSONIfDefined( bodyState( xCartesianPositionIndex ), jsonState, K::x );
            updateFromJSONIfDefined( bodyState( yCartesianPositionIndex ), jsonState, K::y );
            updateFromJSONIfDefined( bodyState( zCartesianPositionIndex ), jsonState, K::z );
            updateFromJSONIfDefined( bodyState( xCartesianVelocityIndex ), jsonState, K::vx );
            updateFromJSONIfDefined( bodyState( yCartesianVelocityIndex ), jsonState, K::vy );
            updateFromJSONIfDefined( bodyState( zCartesianVelocityIndex ), jsonState, K::vz );
            break;
        }
        case keplerianElements:
        {
            const double centralBodyGravitationalParameter =
                    ! isDefined( jsonState, K::centralBodyGravitationalParameter ) && centralBody ?
                        centralBody->getGravityFieldModel( )->getGravitationalParameter( ) :
                        getValue< double >( jsonState, K::centralBodyGravitationalParameter );

            // Detemrine semiMajorAxis, eccentricity, argumentOfPeriapsis

            StateScalarType semiMajorAxis = TUDAT_NAN;
            StateScalarType eccentricity = TUDAT_NAN;
            StateScalarType argumentOfPeriapsis = TUDAT_NAN;

            if ( isDefined( jsonState, K::radius ) || isDefined( jsonState, K::altitude ) )  // circular
            {
                if ( isDefined( jsonState, K::altitude ) )
                {
                    semiMajorAxis = centralBodyAverageRadius + getValue< StateScalarType >( jsonState, K::altitude );
                    usedCentralBodyAverageRadius = true;
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
                        apoapsisDistance = centralBodyAverageRadius +
                                getValue< StateScalarType >( jsonState, K::apoapsisAltitude );
                        usedCentralBodyAverageRadius = true;
                    }
                    else if ( isDefined( jsonState, K::apoapsisDistance ) )
                    {
                        apoapsisDistance = getValue< StateScalarType >( jsonState, K::apoapsisDistance );
                    }

                    StateScalarType periapsisDistance = TUDAT_NAN;
                    if ( isDefined( jsonState, K::periapsisAltitude ) )
                    {
                        periapsisDistance = centralBodyAverageRadius +
                                getValue< StateScalarType >( jsonState, K::periapsisAltitude );
                        usedCentralBodyAverageRadius = true;
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
            bodyState = convertKeplerianToCartesianElements< StateScalarType >( keplerianElements, centralBodyGravitationalParameter );

            break;
        }
        case sphericalElements:
        {
            assertNonnullptrPointer( centralBody );

            double radius;
            if ( isDefined( jsonState, K::altitude ) )
            {
                radius = centralBodyAverageRadius + getValue< StateScalarType >( jsonState, K::altitude );
                usedCentralBodyAverageRadius = true;
            }
            else
            {
                radius = getValue< StateScalarType >( jsonState, K::radius );
            }

            Eigen::Matrix< StateScalarType, 6, 1 > sphericalElements;
            sphericalElements( radiusIndex ) = radius;
            sphericalElements( latitudeIndex ) = getValue< StateScalarType >( jsonState, K::latitude, 0.0 );
            sphericalElements( longitudeIndex ) = getValue< StateScalarType >( jsonState, K::longitude, 0.0 );
            sphericalElements( speedIndex ) = getValue< StateScalarType >( jsonState, K::speed );
            sphericalElements( flightPathIndex ) = getValue< StateScalarType >( jsonState, K::flightPathAngle, 0.0 );
            sphericalElements( headingAngleIndex ) = getValue< StateScalarType >( jsonState, K::headingAngle, 0.0 );

            // Convert to Cartesian elements
            if ( isDefined( jsonState, K::epoch ) || isNaN( epoch ) )
            {
                epoch = getValue< StateScalarType >( jsonState, K::epoch );
            }
            bodyState = tudat::ephemerides::transformStateToGlobalFrame< StateScalarType >(
                        convertSphericalOrbitalToCartesianState( sphericalElements ),
                        epoch, centralBody->getRotationalEphemeris( ) );

            break;
        }
        default:
            handleUnimplementedEnumValue( stateType, stateTypes, unsupportedStateTypes );
        }

        if ( usedCentralBodyAverageRadius )
        {
            if ( isNaN( centralBodyAverageRadius ) )
            {
                // Print error saying that key centralBodyAverageRadius is missing
                getValue< double >( jsonState, K::centralBodyAverageRadius );
            }
            else if ( centralBody )
            {
                if ( ! std::dynamic_pointer_cast< SphericalBodyShapeModel >( centralBody->getShapeModel( ) ) )
                {
                    std::cerr << "Using average radius of a non-spherical body "
                                 "to determine the initial Cartesian state of body." << std::endl;
                }
            }
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
