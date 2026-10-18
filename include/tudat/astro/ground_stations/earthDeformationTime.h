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

}  // namespace basic_astrodynamics
}  // namespace tudat

#endif
