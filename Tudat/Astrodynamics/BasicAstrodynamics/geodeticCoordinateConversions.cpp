/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Montebruck O, Gill E. Satellite Orbits, Springer, 2000.
 *
 */

#include <vector>
#include <cmath>

#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/geodeticCoordinateConversions.h"

namespace tudat
{

namespace coordinate_conversions
{

//! Calculate the ellipticity of an ellipsoid.
double calculateEllipticity( const double flattening )
{
    // Montenbruck and Gill (2000), Eq. (5.86).
    return std::sqrt( 1.0 - ( 1.0 - flattening ) * ( 1.0 - flattening ) );
}

//! Calculate auxiliary quantities for geodetic coordinate conversions.
std::pair< double, double > calculateGeodeticCoordinatesAuxiliaryQuantities(
        const Eigen::Vector3d cartesianState,
        const double equatorialRadius,
        const double ellipticity,
        const double tolerance )
{
    // Precompute square of ellipticity.
    const double ellipticitySquared = ellipticity * ellipticity;

    // Initialize z-intercept value.
    double zInterceptOffset0 = ellipticitySquared * cartesianState.z( );

    // Declare variables to use in iteration.
    double oldZInterceptOffset0 = 0.0, interceptToSurfaceDistance = 0.0,
            sineOfGeodeticLatitude = 0.0, totalZ = 0.0, totalDistance = 0.0;

    // Perform iterative algorithm until tolerance is reached, Montenbruck and Gill (2000),
    // Eq. (5.87).
    do
    {
        // Set old z-intercept offset to current value.
        oldZInterceptOffset0 = zInterceptOffset0;

        // Calculate total z-distance from intercept and total distance from intercept point.
        totalZ = cartesianState.z( ) + zInterceptOffset0;
        totalDistance = std::sqrt( cartesianState.x( ) * cartesianState.x( ) +
                                   cartesianState.y( ) * cartesianState.y( ) +
                                   totalZ * totalZ );

        // Calculate sine of geodetic longitude for current values.
        sineOfGeodeticLatitude = totalZ / totalDistance;

        // Calculate auxiliary variables for current iteration.
        interceptToSurfaceDistance = equatorialRadius /
                std::sqrt( 1.0 - ellipticitySquared *
                           sineOfGeodeticLatitude * sineOfGeodeticLatitude );
        zInterceptOffset0 = interceptToSurfaceDistance * ellipticitySquared *
                sineOfGeodeticLatitude;

    } while ( std::fabs( oldZInterceptOffset0 - zInterceptOffset0 ) > tolerance );

    // Set and return final auxiliary variables.
    return std::make_pair( interceptToSurfaceDistance, zInterceptOffset0 );
}

//! Calculate the Cartesian position from geodetic coordinates.
Eigen::Vector3d convertGeodeticToCartesianCoordinates( const Eigen::Vector3d geodeticCoordinates,
                                                       const double equatorialRadius,
                                                       const double flattening )
{
    // Pre-compute values for efficiency.
    const double sineLatitude = std::sin( geodeticCoordinates.y( ) );
    const double centerOffset = equatorialRadius /
            std::sqrt( 1.0 - flattening * ( 2.0 - flattening ) * sineLatitude * sineLatitude );

    // Calculate Cartesian coordinates, Montenbruck and Gill (2000), Eq. (5.82).
    Eigen::Vector3d cartesianCoordinates = Eigen::Vector3d::Zero( );
    cartesianCoordinates.x( ) =
            ( centerOffset + geodeticCoordinates.x( ) ) *
            std::cos( geodeticCoordinates.y( ) ) * std::cos( geodeticCoordinates.z( ) );
    cartesianCoordinates.y( ) = cartesianCoordinates.x( ) * std::tan( geodeticCoordinates.z( ) );
    cartesianCoordinates.z( ) = ( ( 1.0 - flattening ) * ( 1.0 - flattening ) * centerOffset +
                                  geodeticCoordinates.x( ) ) * sineLatitude;

    // Return Cartesian coordinates.
    return cartesianCoordinates;
}

//! Calculate the altitude over an oblate spheroid of a position vector from auxiliary variables.
double calculateAltitudeOverOblateSpheroid( const Eigen::Vector3d cartesianState,
                                            const double zInterceptOffset,
                                            const double interceptToSurfaceDistance )
{
    // Calculate Cartesian coordinates, Montenbruck and Gill (2000), Eq. (5.88c).
    return std::sqrt( cartesianState.x( ) * cartesianState.x( ) +
                      cartesianState.y( ) * cartesianState.y( ) +
                      ( cartesianState.z( ) + zInterceptOffset ) *
                      ( cartesianState.z( ) + zInterceptOffset ) ) - interceptToSurfaceDistance;
}

//! Calculate the altitude over an oblate spheroid of a position vector.
double calculateAltitudeOverOblateSpheroid( const Eigen::Vector3d cartesianState,
                                            const double equatorialRadius,
                                            const double flattening,
                                            const double tolerance )
{
    // Calculate auxiliary variables.
    std::pair< double, double > auxiliaryVariables =
            calculateGeodeticCoordinatesAuxiliaryQuantities(
                cartesianState, equatorialRadius,  calculateEllipticity( flattening ), tolerance );

    // Calculate and return geodetic altitude.
    return calculateAltitudeOverOblateSpheroid(
                cartesianState, auxiliaryVariables.second, auxiliaryVariables.first );
}

//! Calculate the geodetic latitude from Cartesian position and offset of z-intercept.
double calculateGeodeticLatitude( const Eigen::Vector3d cartesianState,
                                  const double zInterceptOffset )
{
    // Calculate Cartesian coordinates, Montenbruck and Gill (2000), Eq. (5.88b).
    return std::atan2( cartesianState.z( ) +
                       zInterceptOffset, std::sqrt(
                           cartesianState.x( ) * cartesianState.x( ) +
                           cartesianState.y( ) * cartesianState.y( ) ) );
}

//! Calculate the geodetic latitude of a position vector.
double calculateGeodeticLatitude( const Eigen::Vector3d cartesianState,
                                  const double equatorialRadius,
                                  const double flattening,
                                  const double tolerance )
{
    // Calculate auxiliary variables.
    std::pair< double, double > auxiliaryVariables =
            calculateGeodeticCoordinatesAuxiliaryQuantities(
                cartesianState, equatorialRadius, calculateEllipticity( flattening ), tolerance );

    // Calculate and return geodetic latitude.
    return calculateGeodeticLatitude( cartesianState, auxiliaryVariables.second );

}

//! Calculate geodetic coordinates ( altitude, geodetic latitude, longitude ) of a position vector.
Eigen::Vector3d convertCartesianToGeodeticCoordinates( const Eigen::Vector3d cartesianCoordinates,
                                                       const double equatorialRadius,
                                                       const double flattening,
                                                       const double tolerance )
{
    using coordinate_conversions::convertCartesianToSpherical;

    // Determine spherical coordinates of position vector.
    const Eigen::Vector3d sphericalCoordinates
            = convertCartesianToSpherical( cartesianCoordinates );
    Eigen::Vector3d geodeticCoordinates = Eigen::Vector3d::Zero( );

    // Calculate auxiliary variables of geodetic coordinates.
    std::pair< double, double > auxiliaryVariables =
            calculateGeodeticCoordinatesAuxiliaryQuantities(
                cartesianCoordinates, equatorialRadius,
                calculateEllipticity( flattening ), tolerance );

    // Calculate altitude.
    geodeticCoordinates.x( ) = calculateAltitudeOverOblateSpheroid(
                cartesianCoordinates, auxiliaryVariables.second, auxiliaryVariables.first );

    // Calculate geodetic latitude.
    geodeticCoordinates.y( ) = calculateGeodeticLatitude(
                cartesianCoordinates, auxiliaryVariables.second );

    // Set longitude.
    geodeticCoordinates.z( ) = sphericalCoordinates.z( );
    return geodeticCoordinates;
}

} // namespace tudat

} // namespace coordinate_conversions
