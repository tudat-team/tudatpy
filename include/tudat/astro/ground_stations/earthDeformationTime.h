#ifndef TUDAT_EARTH_DEFORMATION_TIME_H
#define TUDAT_EARTH_DEFORMATION_TIME_H

#include <memory>
#include "tudat/basics/timeType.h"

namespace tudat
{
namespace earth_orientation
{
class TerrestrialTimeScaleConverter;
}
namespace basic_astrodynamics
{
//! Convenience default for geocentric TAI/TT/TDB conversion without loading EOP
//! or short-period correction files. Sharing this converter is not required for precision.
std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > getGeocentricTimeConverter( );

//! For a double station query equal to the rounded TDB time of the first or last
//! EOP observation, return that observation's precise TDB time. This prevents
//! rounding from placing an otherwise valid query just outside EOP coverage.
//! Reject a double that cannot distinguish the two observation times.
//! Queries supplied as precise Time never call this function.
Time normalizeScalarTdb( double tdb, const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter >& converter );

//! Validate a closed TT interval for the current Earth models' scalar-UTC phase convention.
/*!
 * Reject ambiguous scalar UTC samples and UTC offset steps, including steps at either
 * endpoint. Smooth pre-1972 UTC drift is permitted. This is a limitation of the current
 * UTC phase representation, not a physical tide singularity. Custom argument callbacks
 * with a different time basis are responsible for their own interval policy.
 */
void validateEarthDeformationUtcInterval( const double lowerTt, const double upperTt );

}  // namespace basic_astrodynamics
}  // namespace tudat

#endif
