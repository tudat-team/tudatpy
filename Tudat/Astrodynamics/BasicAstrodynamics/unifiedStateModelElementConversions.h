/*    Copyright (c) 2010-2016, Delft University of Technology
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
 *      160406    M. Van den Broeck File created.
 *
 *    References
 *      The Unified State Model: Derivation and Applications in Astrodynamics and Navigation,
 *          Vivek Vittaldevs; M.Sc. Thesis TU Delft (2010). Available on repository.tudelft.nl
 *
 *    Notes
 *
 */

#ifndef TUDAT_UNIFIED_STATE_MODEL_ELEMENT_CONVERSIONS_H
#define TUDAT_UNIFIED_STATE_MODEL_ELEMENT_CONVERSIONS_H

#include "Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h"

namespace tudat
{

namespace orbital_element_conversions
{

//! Convert Keplerian to Unified State Model elements.
/*!
 * Converts Keplerian to Unified State Model elements.
 * \param keplerianElements Vector containing Keplerian elements. Order of elements is important!
 *         keplerianElements( 0 ) = semi-major axis,                                            [m]
 *         keplerianElements( 1 ) = eccentricity,                                               [-]
 *         keplerianElements( 2 ) = inclination (in range [0,PI]),                            [rad]
 *         keplerianElements( 3 ) = argument of periapsis,                                    [rad]
 *         keplerianElements( 4 ) = longitude of ascending node,                              [rad]
 *         keplerianElements( 5 ) = true anomaly.                                             [rad]
 * \param centralBodyGravitationalParameter Gravitational parameter of central body.      [m^3/s^2]
 * \return convertedUnifiedStateModelElements Converted state in Unified State Model elements. The order of elements is fixed
 *         convertedUnifiedStateModelElements( 0 ) = C hodograph element,                     [m/s]
 *         convertedUnifiedStateModelElements( 1 ) = Rf1 hodograph element,                   [m/s]
 *         convertedUnifiedStateModelElements( 2 ) = Rf1 hodograph element,                   [m/s]
 *         convertedUnifiedStateModelElements( 3 ) = epsilon1 quaternion element,               [-]
 *         convertedUnifiedStateModelElements( 4 ) = epsilon2 quaternion element,               [-]
 *         convertedUnifiedStateModelElements( 5 ) = epsilon3 quaternion element,               [-]
 *         convertedUnifiedStateModelElements( 6 ) = eta quaternion element.                    [-]
 */
Eigen::Matrix< double, 7, 1 > convertKeplerianToUnifiedStateModelElements(
        const basic_mathematics::Vector6d& keplerianElements,
        const double centralBodyGravitationalParameter );

//! Convert Unified State Model elements to Keplerian elements.
/*!
 * Converts Unified State Model elements to Keplerian elements.
 * \param unifiedStateModelElements Vector containing Unified State Model elements. Order of elements is important!
 *         unifiedStateModelElements( 0 ) = C hodograph element,                              [m/s]
 *         unifiedStateModelElements( 1 ) = Rf1 hodograph element,                            [m/s]
 *         unifiedStateModelElements( 2 ) = Rf1 hodograph element,                            [m/s]
 *         unifiedStateModelElements( 3 ) = epsilon1 quaternion element,                        [-]
 *         unifiedStateModelElements( 4 ) = epsilon2 quaternion element,                        [-]
 *         unifiedStateModelElements( 5 ) = epsilon3 quaternion element,                        [-]
 *         unifiedStateModelElements( 6 ) = eta quaternion element.                             [-]
 * \param centralBodyGravitationalParameter Gravitational parameter of central body.      [m^3/s^2]
 * \return convertedKeplerianElements Converted state in Keplerian elements. The order of elements is fixed!
 *         convertedKeplerianElements( 0 ) = semi-major axis,                                   [m]
 *         convertedKeplerianElements( 1 ) = eccentricity,                                      [-]
 *         convertedKeplerianElements( 2 ) = inclination (in range [0,PI]),                   [rad]
 *         convertedKeplerianElements( 3 ) = argument of periapsis,                           [rad]
 *         convertedKeplerianElements( 4 ) = longitude of ascending node,                     [rad]
 *         convertedKeplerianElements( 5 ) = true anomaly.                                    [rad]
 */
basic_mathematics::Vector6d convertUnifiedStateModelToKeplerianElements(
        const Eigen::Matrix< double, 7, 1 >& unifiedStateModelElements,
        const double centralBodyGravitationalParameter );

} // close namespace orbital_element_conversion

} // close tudat

#endif // TUDAT_UNIFIED_STATE_MODEL_ELEMENT_CONVERSIONS_H
