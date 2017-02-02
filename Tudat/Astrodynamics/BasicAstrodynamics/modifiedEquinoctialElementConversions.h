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
 *      121128    R.C.A. Boon       File created (in progress).
 *      130124    R.C.A. Boon       Removed retrograde factor calculation
 *      130124    R.C.A. Boon       Added Cartesian conversions, improved doxygen comments
 *      130225    D. Dirkx          Added overloaded function for Kepler to MEE that determines
 *                                  retrogradeness based on Kepler state
 *      130301    R.C.A. Boon       Minor textual changes
 *      130305    R.C.A. Boon       Replaced Eigen::VectorXd by Eigen::Vector6d.
 *      140221    H.P. Gijsen       Moved the Modified Equinotial Elements indeces to stateVectorIndices.h
 *
 *    References
 *      Verified Interval Propagation, Bart R�mgens; Delft (2012?). Code archive.
 *      Modified Equinoctial Orbital Elements, author unknown;
 *          http://www.cdeagle.com/pdf/mee.pdf (2010?).
 *      Survey of Orbital Element Sets, Gerald R. Hintz; Journal of Guidance, Control and
 *          Dynamics (2008, Vol. 31 - Nr. 3).
 *      Code archive, E. Heeren (fellow Tudat developer).
 *
 *    Notes
 *
 */

#ifndef TUDAT_MODIFIED_EQUINOCTIAL_ELEMENT_CONVERSIONS_H
#define TUDAT_MODIFIED_EQUINOCTIAL_ELEMENT_CONVERSIONS_H

#include "Tudat/Basics/basicTypedefs.h"

namespace tudat
{

namespace orbital_element_conversions
{

//! Convert Keplerian to modified equinoctial orbital elements using implicit MEE equation set.
/*!
 * Converts Keplerian to modified equinoctial elements using the prograde/retrograde equation
 * determined by the values of the input Kepler elements. If input exceeds allowable ranges,
 * an error is thrown.
 * \param keplerianElements Vector containing Keplerian elements. Order of elements is important!
 *         keplerianElements( 0 ) = semi-major axis,                                            [m]
 *         keplerianElements( 1 ) = eccentricity,                                               [-]
 *         keplerianElements( 2 ) = inclination (in range [0,PI]),                            [rad]
 *         keplerianElements( 3 ) = argument of periapsis,                                    [rad]
 *         keplerianElements( 4 ) = longitude of ascending node,                              [rad]
 *         keplerianElements( 5 ) = true anomaly.                                             [rad]
 * \return Converted state in modified equinoctial elements. The order of elements is fixed!
 *         modifiedEquinoctialElements( 0 ) = semi-latus rectum,                                [m]
 *         modifiedEquinoctialElements( 1 ) = f-element,                                        [-]
 *         modifiedEquinoctialElements( 2 ) = g-element,                                        [-]
 *         modifiedEquinoctialElements( 3 ) = h-element,                                        [-]
 *         modifiedEquinoctialElements( 4 ) = k-element,                                        [-]
 *         modifiedEquinoctialElements( 5 ) = true longitude.                                 [rad]
 */
Eigen::Vector6d convertKeplerianToModifiedEquinoctialElements(
        const Eigen::Vector6d& keplerianElements );

//! Convert Keplerian to modified equinoctial orbital elements using MEE explicit equation set.
/*!
 * Converts Keplerian to modified equinoctial elements using one of two sets of equations specified
 * by the user. If input exceeds allowable ranges, an error is thrown.
 * \param keplerianElements Vector containing Keplerian elements. Order of elements is important!
 *         keplerianElements( 0 ) = semi-major axis,                                            [m]
 *         keplerianElements( 1 ) = eccentricity,                                               [-]
 *         keplerianElements( 2 ) = inclination (in range [0,PI]),                            [rad]
 *         keplerianElements( 3 ) = argument of periapsis,                                    [rad]
 *         keplerianElements( 4 ) = longitude of ascending node,                              [rad]
 *         keplerianElements( 5 ) = true anomaly.                                             [rad]
 * \param avoidSingularityAtPiInclination Boolean flag to indicate whether the set of equations for
 *          the inclination = PI singular case are to be used. Take note: the same set of equations
 *          is required for conversion back to Keplerian elements to retrieve original state!
 * \return Converted state in modified equinoctial elements. The order of elements is fixed!
 *         modifiedEquinoctialElements( 0 ) = semi-latus rectum,                                [m]
 *         modifiedEquinoctialElements( 1 ) = f-element,                                        [-]
 *         modifiedEquinoctialElements( 2 ) = g-element,                                        [-]
 *         modifiedEquinoctialElements( 3 ) = h-element,                                        [-]
 *         modifiedEquinoctialElements( 4 ) = k-element,                                        [-]
 *         modifiedEquinoctialElements( 5 ) = true longitude.                                 [rad]
 */
Eigen::Vector6d convertKeplerianToModifiedEquinoctialElements(
        const Eigen::Vector6d& keplerianElements,
        const bool avoidSingularityAtPiInclination );

//! Convert modified equinoctial to Keplerian orbital elements.
/*!
 * Converts modified equinoctial elements to Keplerian using one of two sets of equations specified
 * by the user.
 * \param modifiedEquinoctialElements Vector containing modified equinoctial elements. Order of
 *          elements is important!
 *         modifiedEquinoctialElements( 0 ) = semi-latus rectum,                                [m]
 *         modifiedEquinoctialElements( 1 ) = f-element,                                        [-]
 *         modifiedEquinoctialElements( 2 ) = g-element,                                        [-]
 *         modifiedEquinoctialElements( 3 ) = h-element,                                        [-]
 *         modifiedEquinoctialElements( 4 ) = k-element,                                        [-]
 *         modifiedEquinoctialElements( 5 ) = true longitude.                                 [rad]
 * \param avoidSingularityAtPiInclination Boolean flag to indicate whether the set of equations for
 *          the inclination = PI singular case are to be used. Take note: the same set of equations
 *          is required for conversion back to modified equinoctial elements to retrieve original
 *          state!
 * \return Converted state in Keplerian elements. The order of elements is fixed!
 *         keplerianElements( 0 ) = semi-major axis,                                            [m]
 *         keplerianElements( 1 ) = eccentricity,                                               [-]
 *         keplerianElements( 2 ) = inclination [0, 180],                                     [rad]
 *         keplerianElements( 3 ) = argument of periapsis,                                    [rad]
 *         keplerianElements( 4 ) = longitude of ascending node,                              [rad]
 *         keplerianElements( 5 ) = true anomaly.                                             [rad]
 */
Eigen::Vector6d convertModifiedEquinoctialToKeplerianElements(
        const Eigen::Vector6d& modifiedEquinoctialElements,
        const bool avoidSingularityAtPiInclination );

//! Convert Cartesian to modified equinoctial orbital elements using implicit MEE equation set.
/*!
 * Converts Cartesian to modified equinoctial elements using one of two sets of equations implicitly
 * determined from intermediate inclination.
 * \param cartesianElements Vector containing Cartesian elements. Order of elements is important!
 *         cartesianElements( 0 ) = x-position coordinate,                                      [m]
 *         cartesianElements( 1 ) = y-position coordinate,                                      [m]
 *         cartesianElements( 2 ) = z-position coordinate,                                      [m]
 *         cartesianElements( 3 ) = x-velocity coordinate,                                    [m/s]
 *         cartesianElements( 4 ) = y-velocity coordinate,                                    [m/s]
 *         cartesianElements( 5 ) = z-velocity coordinate.                                    [m/s]
 * \param centralBodyGravitationalParameter Gravitational parameter of central body.
 * \return Converted state in modified equinoctial elements. The order of elements is fixed!
 *         modifiedEquinoctialElements( 0 ) = semi-latus rectum,                                [m]
 *         modifiedEquinoctialElements( 1 ) = f-element,                                        [-]
 *         modifiedEquinoctialElements( 2 ) = g-element,                                        [-]
 *         modifiedEquinoctialElements( 3 ) = h-element,                                        [-]
 *         modifiedEquinoctialElements( 4 ) = k-element,                                        [-]
 *         modifiedEquinoctialElements( 5 ) = true longitude.                                 [rad]
 */
Eigen::Vector6d convertCartesianToModifiedEquinoctialElements(
        const Eigen::Vector6d& cartesianElements,
        const double centralBodyGravitationalParameter );

//! Convert Cartesian to modified equinoctial orbital elements using explicit MEE equation set.
/*!
 * Converts Cartesian to modified equinoctial elements using one of two sets of equations specified
 * by the user.
 * \param cartesianElements Vector containing Cartesian elements. Order of elements is important!
 *         cartesianElements( 0 ) = x-position coordinate,                                      [m]
 *         cartesianElements( 1 ) = y-position coordinate,                                      [m]
 *         cartesianElements( 2 ) = z-position coordinate,                                      [m]
 *         cartesianElements( 3 ) = x-velocity coordinate,                                    [m/s]
 *         cartesianElements( 4 ) = y-velocity coordinate,                                    [m/s]
 *         cartesianElements( 5 ) = z-velocity coordinate.                                    [m/s]
 * \param centralBodyGravitationalParameter Gravitational parameter of central body.
 * \param avoidSingularityAtPiInclination Boolean flag to indicate whether the set of equations for
 *          the inclination = PI singular case are to be used. Take note: the same set of equations
 *          is required for conversion back to Cartesian elements to retrieve original state!
 * \return Converted state in modified equinoctial elements. The order of elements is fixed!
 *         modifiedEquinoctialElements( 0 ) = semi-latus rectum,                                [m]
 *         modifiedEquinoctialElements( 1 ) = f-element,                                        [-]
 *         modifiedEquinoctialElements( 2 ) = g-element,                                        [-]
 *         modifiedEquinoctialElements( 3 ) = h-element,                                        [-]
 *         modifiedEquinoctialElements( 4 ) = k-element,                                        [-]
 *         modifiedEquinoctialElements( 5 ) = true longitude.                                 [rad]
 */
Eigen::Vector6d convertCartesianToModifiedEquinoctialElements(
        const Eigen::Vector6d& cartesianElements,
        const double centralBodyGravitationalParameter,
        const bool avoidSingularityAtPiInclination );

//! Convert modified equinoctial elements to Cartesian orbital elements.
/*!
 * Converts modified equinoctial elements to Cartesian orbital elements using one of two sets of
 * equations specified by the user.
 * \param modifiedEquinoctialElements Vector containing modified equinoctial elements. Order of
 * elements is important!
 *         modifiedEquinoctialElements( 0 ) = semi-latus rectum,                                [m]
 *         modifiedEquinoctialElements( 1 ) = f-element,                                        [-]
 *         modifiedEquinoctialElements( 2 ) = g-element,                                        [-]
 *         modifiedEquinoctialElements( 3 ) = h-element,                                        [-]
 *         modifiedEquinoctialElements( 4 ) = k-element,                                        [-]
 *         modifiedEquinoctialElements( 5 ) = true longitude.                                 [rad]
 * \param centralBodyGravitationalParameter Gravitational parameter of central body.
 * \param avoidSingularityAtPiInclination Boolean flag to indicate whether the set of equations for
 *          the inclination = PI singular case are to be used. Take note: the same set of equations
 *          is required for conversion back to modified equinoctial elements to retrieve original
 *          state!
 * \return Converted state in Cartesian elements. The order of elements is fixed!
 *         cartesianElements( 0 ) = x-position coordinate,                                      [m]
 *         cartesianElements( 1 ) = y-position coordinate,                                      [m]
 *         cartesianElements( 2 ) = z-position coordinate,                                      [m]
 *         cartesianElements( 3 ) = x-velocity coordinate,                                    [m/s]
 *         cartesianElements( 4 ) = y-velocity coordinate,                                    [m/s]
 *         cartesianElements( 5 ) = z-velocity coordinate.                                    [m/s]
 */
Eigen::Vector6d convertModifiedEquinoctialToCartesianElements(
        const Eigen::Vector6d& modifiedEquinoctialElements,
        const double centralBodyGravitationalParameter,
        const bool avoidSingularityAtPiInclination );

} // namespace orbital_element_conversions

} // namespace tudat

#endif // TUDAT_MODIFIED_EQUINOCTIAL_ELEMENT_CONVERSIONS_H
