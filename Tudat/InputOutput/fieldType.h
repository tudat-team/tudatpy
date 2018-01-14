/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Notes
 *      This file contains the definitation for fieldtype, and lists several common field types.
 *
 */

#ifndef TUDAT_FIELD_TYPE_H
#define TUDAT_FIELD_TYPE_H

#include <cstring>

#include <boost/functional/hash.hpp>

namespace tudat
{
namespace input_output
{

//! FieldType is a hash number of a description.
/*!
 * FieldType is a hash number (see note at the end of this comment) created from a description for
 * that FieldType. This means that it is important to maintain different descriptions of each field,
 * as the same description leads to the same hash. It is strongly advised to use the following
 * convention:
 *
 * A description is built up from a category path. Each category starts with a capital letter and
 * must not contain any spaces. Following the category name is a colon followed by a space and the
 * description of the next sub-category or the actual string identifying the specific field. The
 * same conventions as for categories are used on the string identifier.
 *
 * Examples can be seen in FieldType.h
 *
 * Note: A hash number is a fixed length key that uniquely indexes a variable length key. Hash
 * numbers are mostly used to accelerate table lookup or data comparison tasks such as finding items
 * in a database (Wikipedia). Here, each field type (variable length) is associated with a hash
 * number of fixed length.
 */
typedef std::size_t FieldType;

namespace field_types
{

//! Convert string into hash number.
static FieldType hash_constructor( const char* text )
{
    // Convert text string into hash number.
    boost::hash< std::string > contructor;
    return contructor( text );
}

// -----------------------------
// General fields
// -----------------------------
namespace general
{
    //! Object name [-].
    static const FieldType name = hash_constructor( "General: Object_name" );

    //! Object identification code [-].
    static const FieldType id   = hash_constructor( "General: Object_id" );

    //! Parameter name [-].
    static const FieldType parameterName = hash_constructor( "General: Parameter_name" );

    //! Parameter value [-].
    static const FieldType parameterValue = hash_constructor( "General: Parameter_value" );
}

// -----------------------------
// Time related fields
// -----------------------------
namespace time
{
    //! Epoch Julian Date [JD, UT1].
    static const FieldType epoch = hash_constructor( "Time: Julian_date_epoch[UT1]" );
}

// -----------------------------
// State related fields
// -----------------------------
namespace state
{
    //! Cartesian X coordinate [m].
    static const FieldType cartesianXCoordinate
        = hash_constructor( "State: Cartesian: X_coordinate" );

    //! Cartesian Y coordinate [m].
    static const FieldType cartesianYCoordinate
        = hash_constructor( "State: Cartesian: Y_coordinate" );

    //! Cartesian Z coordinate [m].
    static const FieldType cartesianZCoordinate
    = hash_constructor( "State: Cartesian: Z_coordinate" );

    //! Cartesian X velocity [m/s].
    static const FieldType cartesianXVelocity
        = hash_constructor( "State: Cartesian: X_velocity" );

    //! Cartesian Y velocity [m/s].
    static const FieldType cartesianYVelocity
        = hash_constructor( "State: Cartesian: Y_velocity" );

    //! Cartesian Z velocity [m/s].
    static const FieldType cartesianZVelocity
        = hash_constructor( "State: Cartesian: Z_velocity" );

    //! Cartesian X acceleration [m/s].
    static const FieldType cartesianXAcceleration
        = hash_constructor( "State: Cartesian: X_acceleration" );

    //! Cartesian Y acceleration [m/s].
    static const FieldType cartesianYAcceleration
        = hash_constructor( "State: Cartesian: Y_acceleration" );

    //! Cartesian Z acceleration [m/s].
    static const FieldType cartesianZAcceleration
    = hash_constructor( "State: Cartesian: Z_acceleration" );

    //! Eccentricity, e [-].
    static const FieldType eccentricity
        = hash_constructor( "State: Kepler_element: Eccentricity" );

    //! Inclination with respect to xy-plane, i [rad].
    static const FieldType inclination
        = hash_constructor( "State: Kepler_element: Inclination" );

    //! Longitude of Ascending Node, OMEGA [rad].
    static const FieldType longitudeOfAscendingNode
        = hash_constructor( "State: Kepler_element: Longitude_of_the_ascending_node" );

    //! Argument of Perifocus, w [rad].
    static const FieldType argumentOfPeriapsis
        = hash_constructor( "State: Kepler_element: Argument_of_periapsis" );

    //! Time of periapsis passage [JD, UT1].
    static const FieldType timeOfPeriapsisPassage
        = hash_constructor( "State: Kepler_element: Time_of_periapsis" );

    //! Mean anomaly, M [rad].
    static const FieldType meanAnomaly
        = hash_constructor( "State: Kepler_element: Mean_anomaly" );

    //! True anomaly, theta [rad].
    static const FieldType trueAnomaly
        = hash_constructor( "State: Kepler_element: True_anomaly" );

    //! Semi-major axis, a [m].
    static const FieldType semiMajorAxis
        = hash_constructor( "State: Kepler_element: Semi_major_axis" );

    //! Apoapsis distance [m].
    static const FieldType apoapsisDistance
    = hash_constructor( "State: Kepler_element: Apoapsis_distance" );

    //! Periapsis distance [m]
    static const FieldType periapsisDistance
        = hash_constructor( "State: Kepler_element: Periapsis_distance" );

    //! Mean motion, n [rad/s].
    static const FieldType meanMotion = hash_constructor( "State: Mean_motion" );

    //! Orbital period [s].
    static const FieldType orbitalPeriod = hash_constructor( "State: Orbital_period" );
}

} // namespace field_types
} // namespace input_output
} // namespace tudat

#endif // TUDAT_FIELD_TYPE_H
