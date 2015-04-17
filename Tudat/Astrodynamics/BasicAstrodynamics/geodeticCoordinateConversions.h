/*    Copyright (c) 2010-2015, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      130301    D. Dirkx          Migrated from personal code.
 *      130308    E.D. Brandon      Minor changes.
 *
 *    References
 *      Montebruck O, Gill E. Satellite Orbits, Springer, 2000.
 *
 *    Notes
 *
 */

#ifndef TUDAT_GEODETIC_COORDINATE_CONVERSIONS_H
#define TUDAT_GEODETIC_COORDINATE_CONVERSIONS_H

#include <utility>

#include <Eigen/Core>

namespace tudat
{

namespace coordinate_conversions
{

//! Calculate the ellipticity of an ellipsoid.
/*!
 * Calculates the ellipticity of an ellipsoid from its flattening. From Montenbruck & Gill (2000).
 * Note: The ellipticity is not to be confused with the eccentricity.
 * \param flattening Flattening of the ellipsoid.
 * \return Ellipticity of the ellipsoid.
 */
double calculateEllipticity( const double flattening );

//! Calculate auxiliary quantities for geodetic coordinate conversions.
/*!
 * Calculates auxiliary quantities for geodetic coordinate conversions.
 * (Montenbruck & Gill, 2010, Fig 5.12).
 * The auxiliary quantities are determined by creating a line through the cartesianPosition,
 * perpendicular to the surface and finding its intercept with the z-axis. The distance between
 * the body surface and the z-axis intercept is the first auxiliary quantity, the offset in
 * z-direction of the intersect point from the origin is the second auxiliary variable.
 * \param cartesianPosition Cartesian position in body-fixed frame where altitude is to
 *          be determined.
 * \param equatorialRadius Equatorial radius of oblate spheroid.
 * \param ellipticity Ellipticity of oblate spheroid.
 * \param tolerance Convergence criterion for iterative algorithm that is employed. Represents the
 *          required change of position (in m) between two iterations.
 * \return Auxiliary parameters for geodetic coordinate conversions.
 */
std::pair< double, double > calculateGeodeticCoordinatesAuxiliaryQuantities(
        const Eigen::Vector3d cartesianPosition,
        const double equatorialRadius,
        const double ellipticity,
        const double tolerance );

//! Calculate the Cartesian position from geodetic coordinates.
/*!
 * Calculates the Cartesian position from geodetic coordinates
 * (altitude, geodetic latitude, longitude).
 * \param geodeticCoordinates Geodetic coordinates w.r.t. given body.
 * \param equatorialRadius Equatorial radius of oblate spheroid.
 * \param flattening Flattening of oblate spheroid.
 * \return Cartesian position in body-fixed frame.
 */
Eigen::Vector3d convertGeodeticToCartesianCoordinates( const Eigen::Vector3d geodeticCoordinates,
                                                       const double equatorialRadius,
                                                       const double flattening );

//! Calculate the altitude over an oblate spheroid of a position vector from auxiliary variables.
/*!
 * Calculates the altitude over an oblate spheroid of a position vector.
 * This function gets the auxiliary variables(see calculateGeodeticCoordinatesAuxiliaryQuantities)
 * as input. These values are determined by drawing the line L from cartesianPosition
 * perpendicular to the ellipsoid and calculating the intercept with the z-axis.
 * \param cartesianPosition Cartesian position in body-fixed frame where altitude is to
 * be determined.
 * \param zInterceptOffset Offset of intercept of line L with z-axis from origin in z-direction.
 * \param interceptToSurfaceDistance Distance from intercept of line L with z-axis to body surface
 * \return Altitude above specified oblate spheroid at requested point.
 * \sa calculateGeodeticCoordinatesAuxiliaryQuantities
 */
double calculateAltitudeOverOblateSpheroid( const Eigen::Vector3d cartesianPosition,
                                            const double zInterceptOffset,
                                            const double interceptToSurfaceDistance );

//! Calculate the altitude over an oblate spheroid of a position vector.
/*!
 * Calculates the altitude over an oblate spheroid of a position vector.
 * The algorithm that is used is iterative, so that it requires a tolerance (in m) for the
 * difference of associated geodetic position between two iterations.
 * \param cartesianPosition Cartesian position in body-fixed frame where altitude is to
 * be determined.
 * \param equatorialRadius Equatorial radius of oblate spheroid.
 * \param flattening Flattening of oblate spheroid.
 * \param tolerance Convergence criterion for iterative algorithm that is employed. Represents the
 * required change of position (in m) between two iterations.
 * \return Altitude above specified oblate spheroid at requested point.
 */
double calculateAltitudeOverOblateSpheroid( const Eigen::Vector3d cartesianPosition,
                                            const double equatorialRadius,
                                            const double flattening,
                                            const double tolerance );

//! Calculate the geodetic latitude from Cartesian position and offset of z-intercept.
/*!
 * Calculates the geodetic latitude from Cartesian position and offset of z-intercept.
 * This intercept is determined by drawing the line from cartesianPosition perpendicular to the
 * ellipsoid and calculating the offset from the origin where it intercepts the z-axis.
 * \param cartesianPosition Cartesian position in body-fixed frame where geodetic latitude is to
 * be determined.
 * \param zInterceptOffset Offset from origin of intersection with z-axis of line perpendicular to
 * surface from cartesianPosition.
 * \return Geodetic latitude above specified oblate spheroid at requested point.
 */
double calculateGeodeticLatitude( const Eigen::Vector3d cartesianPosition,
                                  const double zInterceptOffset );

//! Calculate the geodetic latitude of a position vector.
/*!
 * Calculates the geodetic latitude of a position vector on an oblate spheroid.
 * The algorithm that is used is iterative, so that it requires a tolerance (in m) for the
 * difference of associated geodetic position between two iterations.
 * \param cartesianPosition Cartesian position in body-fixed frame where geodetic latitude is to
 * be determined.
 * \param equatorialRadius Equatorial radius of oblate spheroid.
 * \param flattening Flattening of oblate spheroid.
 * \param tolerance Convergence criterion for iterative algorithm that is employed. Represents the
 * required change of position (in m) between two iterations.
 * \return Geodetic latitude above specified oblate spheroid at requested point.
 */
double calculateGeodeticLatitude( const Eigen::Vector3d cartesianPosition,
                                  const double equatorialRadius,
                                  const double flattening,
                                  const double tolerance );

//! Calculate geodetic coordinates (altitude, geodetic latitude, longitude) of a position vector.
/*!
 * Calculates the geodetic coordinates (altitude, geodetic latitude, longitude)
 * of a position vector. The algorithm that is used is iterative, so that it requires a tolerance
 * (in m) for the difference of associated geodetic position between two iterations.
 * \param cartesianCoordinates Cartesian position in body-fixed frame where geodetic coordinates 
 *          are to be determined.
 * \param equatorialRadius Equatorial radius of oblate spheroid.
 * \param flattening Flattening of oblate spheroid.
 * \param tolerance Convergence criterion for iterative algorithm that is employed. Represents the
 *          required change of position (in m) between two iterations.
 * \return Geodetic coordinates at requested point.
 */
Eigen::Vector3d convertCartesianToGeodeticCoordinates( const Eigen::Vector3d cartesianCoordinates,
                                                       const double equatorialRadius,
                                                       const double flattening,
                                                       const double tolerance );

} // namespace coordinate_conversions

} // namespace tudat

#endif // TUDAT_GEODETIC_COORDINATE_CONVERSIONS_H
