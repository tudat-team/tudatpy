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

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_CELESTIALREFERENCEFRAMES_H
