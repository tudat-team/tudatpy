/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#ifndef TUDAT_AOD1B_GRAVITY_FIELD_VARIATIONS_H
#define TUDAT_AOD1B_GRAVITY_FIELD_VARIATIONS_H

#include <memory>
#include "tudat/io/readAod1bGravityField.h"
#include "tudat/astro/gravitation/gravityFieldVariations.h"
#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"

namespace tudat
{
namespace gravitation
{

//! Build additive RL07 anomalies, linearly interpolated in GPS seconds with no extrapolation.
//! The selected Earth rotation supplies the converter; the output block starts at degree two/order zero.
std::shared_ptr< GravityFieldVariations > createAod1bGravityFieldVariations(
        const std::vector< input_output::Aod1bGravityFieldData >& data,
        input_output::Aod1bComponent component,
        int maximumDegree,
        int maximumOrder,
        double outputGravitationalParameter,
        double outputReferenceRadius,
        const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter >& converter );

}  // namespace gravitation
}  // namespace tudat
#endif
