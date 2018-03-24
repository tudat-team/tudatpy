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
 *      The Unified State Model: Derivation and Applications in Astrodynamics and Navigation,
 *          Vivek Vittaldevs; M.Sc. Thesis TU Delft (2010). Available on repository.tudelft.nl
 *
 */

#ifndef TUDAT_UNIFIED_STATE_MODEL_ELEMENT_CONVERSIONS_H
#define TUDAT_UNIFIED_STATE_MODEL_ELEMENT_CONVERSIONS_H

#include "Tudat/Basics/basicTypedefs.h"

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
        const Eigen::Vector6d& keplerianElements,
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
Eigen::Vector6d convertUnifiedStateModelToKeplerianElements(
        const Eigen::Matrix< double, 7, 1 >& unifiedStateModelElements,
        const double centralBodyGravitationalParameter );

} // namespace orbital_element_conversions

} // close tudat

#endif // TUDAT_UNIFIED_STATE_MODEL_ELEMENT_CONVERSIONS_H
