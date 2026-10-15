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
//! Time scale converter used by deformation and station motion models that are created without one.
/*!
 * This converter has no Earth orientation data, so it converts only between TAI, TT and TDB (at the geocentre), which is
 * all these models need when no Earth rotation model is used. The same converter is returned on each call, so the
 * data files of a full converter are not loaded for each model.
 */
std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > getGeocentricTimeConverter( );

//! Check that the velocity of an Earth deformation model can be computed numerically over the given TT interval.
/*!
 * The solid Earth tide and ocean tide loading models compute part of their tidal arguments from the time of day in UTC. If a
 * leap second (or another jump in TAI - UTC) occurs within the interval, or at its start or end, UTC jumps within the
 * interval, and a numerical derivative over this interval is wrong. In that case, this function throws a std::domain_error.
 * The gradual change of TAI - UTC before 1972 is allowed. The interval is also rejected if it contains an epoch during a leap
 * second, which cannot be expressed as a number of UTC seconds since J2000.
 * \param lowerTt Start of the interval (TT seconds since J2000)
 * \param upperTt End of the interval (TT seconds since J2000)
 */
void validateEarthDeformationUtcInterval( const double lowerTt, const double upperTt );

}  // namespace basic_astrodynamics
}  // namespace tudat

#endif
