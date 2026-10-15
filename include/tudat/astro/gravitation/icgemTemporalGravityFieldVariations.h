/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#ifndef TUDAT_ICGEM_TEMPORAL_GRAVITY_FIELD_VARIATIONS_H
#define TUDAT_ICGEM_TEMPORAL_GRAVITY_FIELD_VARIATIONS_H

#include "tudat/io/readIcgemGravityField.h"
#include "tudat/astro/gravitation/gravityFieldVariations.h"
#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"

namespace tudat
{
namespace gravitation
{

//! Explicit interpretation of published trend and periodic units; no convention is inferred from a header.
enum class IcgemYearConvention { julianYear, fractionalGregorianYear };

//! Construct piecewise temporal anomalies relative to the supplied nominal construction snapshot.
//! Julian years use elapsed TAI SI seconds/31557600; fractional years use the declared source calendar.
//! Format 1.0 date-only tokens mean midnight; format 2.0 requires HHMM or explicit source-scale overrides.
//! Format 2.0 uses half-open gfct validity intervals; format 1.0 has no artificial validity bounds.
//! An exact converted interval start selects the interval beginning there.
//! Distinct starts or ends that round to the same double TDB cannot be distinguished and are rejected.
//! Later nominal-field updates add an ordinary offset to this adopted baseline; they do not rebase the model.
std::shared_ptr< GravityFieldVariations > createIcgemTemporalGravityFieldVariations(
        const input_output::IcgemGravityFieldData& data,
        basic_astrodynamics::TimeScales sourceTimeScale,
        IcgemYearConvention yearConvention,
        int maximumDegree,
        int maximumOrder,
        const std::map< std::string, double >& dateOverrides,
        const Eigen::MatrixXd& nominalCosine,
        const Eigen::MatrixXd& nominalSine,
        double outputGravitationalParameter,
        double outputReferenceRadius,
        const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter >& converter );

}  // namespace gravitation
}  // namespace tudat
#endif
