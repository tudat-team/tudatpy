/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#ifndef TUDAT_CREATE_TIME22_GRAVITY_FIELD_VARIATIONS_H
#define TUDAT_CREATE_TIME22_GRAVITY_FIELD_VARIATIONS_H

#include "tudat/simulation/environment_setup/createEarthGravityFieldVariations.h"

namespace tudat
{
namespace simulation_setup
{
//! Read TiME22 revision 1.1 paired ICGEM fields from directory/ATM or directory/OCN at settings construction.
//! The published loading and Doodson-Warburg phase corrections are retained; no extra correction is applied.
//! Optional admittance replaces the 57 ocean majors by 435 targets; the atmospheric 16-by-16 map is the identity.
//! Only the selected coefficients (degrees two through 360) are retained in the returned canonical settings.
std::shared_ptr< GravityFieldVariationSettings > time22GravityFieldVariationSettings(
        const std::string& directory,
        gravitation::HarmonicGravityComponent component,
        int maximumDegree,
        int maximumOrder,
        bool applyAdmittance = true,
        const std::shared_ptr< ModelInterpolationSettings >& interpolation = nullptr );

}  // namespace simulation_setup
}  // namespace tudat
#endif
