#include "tudat/astro/ground_stations/earthDeformationTime.h"

#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"

namespace tudat
{
namespace basic_astrodynamics
{

std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > getGeocentricTimeConverter( )
{
    static const auto converter = std::make_shared< earth_orientation::TerrestrialTimeScaleConverter >( nullptr, nullptr );
    return converter;
}

}  // namespace basic_astrodynamics
}  // namespace tudat
