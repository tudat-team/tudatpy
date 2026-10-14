#ifndef TUDAT_EARTH_DEFORMATION_TIME_H
#define TUDAT_EARTH_DEFORMATION_TIME_H

namespace tudat
{
namespace basic_astrodynamics
{
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
