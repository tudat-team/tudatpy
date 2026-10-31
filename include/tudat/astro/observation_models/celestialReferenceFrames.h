/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CELESTIALREFERENCEFRAMES_H
#define TUDAT_CELESTIALREFERENCEFRAMES_H

#include <memory>

#include <Eigen/Core>

#include "tudat/astro/observation_models/observationAncillarySettings.h"

namespace tudat
{

namespace observation_models
{

//! Return the rotation from ICRF/J2000 to the celestial frame selected by ancillary settings (not defined for a custom pole).
/*!
 * The mean- and true-of-date frames use SOFA: iauPmat06/iauPnm06a (IAU 2006 and IAU 2006/2000A, which start from the ICRF), and
 * iauPmat76/iauPnm80 (IAU 1976 and IAU 1976/1980, which start from FK5 J2000) preceded by the ICRF to FK5 J2000 frame bias (iauBp00).
 * The B1950 frame is NAIF's FK4 frame: a single constant rotation, which does not include the E-terms of aberration, nor the
 * epoch-dependent FK4 equinox and proper-motion errors (see convertFk4B1950ToIcrfRightAscensionDeclination for observed positions).
 */
Eigen::Matrix3d getPositionAngleReferenceFrameRotationFromJ2000(
        const double observationTime,
        const std::shared_ptr< ObservationAncillarySimulationSettings >& ancillarySettings = nullptr );

//! Return the celestial north-pole direction selected by ancillary settings, expressed in ICRF/J2000.
/*!
 * \param observationTime Observation reception time in TDB seconds since J2000. This epoch is used for a time-dependent frame unless
 * an explicit reference epoch is present in the ancillary settings.
 * \param ancillarySettings Settings selecting a predefined celestial frame or custom pole. A null pointer selects J2000.
 * \return Unit north-pole direction expressed in ICRF/J2000.
 */
Eigen::Vector3d getPositionAngleReferencePoleInJ2000(
        const double observationTime,
        const std::shared_ptr< ObservationAncillarySimulationSettings >& ancillarySettings = nullptr );

//! Return whether the ancillary settings select directions that include stellar aberration (false by default).
bool getPositionAngleUseStellarAberration( const std::shared_ptr< ObservationAncillarySimulationSettings >& ancillarySettings = nullptr );

//! Convert an observed FK4 B1950.0 right ascension and declination to the ICRF (Tudat's J2000 frame).
/*!
 * Uses the standard SOFA transformations iauFk45z (FK4 B1950.0 to FK5 J2000.0, Aoki et al. 1983, which removes the E-terms of
 * aberration and accounts for the FK4 equinox and proper-motion system, assuming zero proper motion in FK5) followed by iauFk5hz
 * (FK5 to the Hipparcos/ICRS system, assuming zero proper motion in Hipparcos). Both steps depend on the epoch of observation.
 * \param rightAscension FK4 B1950.0 right ascension at the epoch of observation [rad]
 * \param declination FK4 B1950.0 declination at the epoch of observation [rad]
 * \param observationTime Epoch of observation, in TT seconds since J2000 (TDB may be used, as the difference is negligible)
 * \return ICRF right ascension and declination [rad]
 */
Eigen::Vector2d convertFk4B1950ToIcrfRightAscensionDeclination( const double rightAscension,
                                                                const double declination,
                                                                const double observationTime );

//! Print a warning (once per run) if stellar aberration is computed with a global frame origin other than the barycentre.
void warnIfStellarAberrationOriginIsNotBarycentric( const bool isGlobalFrameOriginBarycentric );

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_CELESTIALREFERENCEFRAMES_H
