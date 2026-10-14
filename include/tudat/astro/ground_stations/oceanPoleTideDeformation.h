#ifndef TUDAT_OCEAN_POLE_TIDE_DEFORMATION_H
#define TUDAT_OCEAN_POLE_TIDE_DEFORMATION_H

#include "tudat/io/readOceanPoleTideData.h"
#include "tudat/astro/ground_stations/bodyDeformationModel.h"
#include "tudat/astro/earth_orientation/polarMotionCalculator.h"

namespace tudat
{
namespace basic_astrodynamics
{
//! The official half-degree grid stops at +/-89.75 degrees: it does not define polar caps.
enum class OceanPoleTideLatitudePolicy { requireGridCoverage, clampToNearestGridLatitudeApproximation };

//! IERS 2010 Eq. (29), with radians as wobble input and metres as output (radial/north/east).
//! Defaults use Chapter 1 constants, not the older G in the published 2003-mean-pole test.
Eigen::Vector3d calculateOceanPoleTideLoading( const Eigen::Matrix< double, 3, 2 >& coefficients,
                                               const Eigen::Vector2d& wobbleRadians,
                                               const double scaleMetresPerRadian = 5340.428562229677 );

//! Bilinear interpolation in the model layer; longitude wraps only for grids spanning a complete period.
//! A valid coordinate outside the supplied coverage raises std::out_of_range under the strict policy.
Eigen::Matrix< double, 3, 2 > interpolateOceanPoleTideCoefficients(
        const input_output::OceanPoleTideGrid& grid,
        const double latitudeDegrees,
        const double longitudeDegrees,
        const OceanPoleTideLatitudePolicy latitudePolicy = OceanPoleTideLatitudePolicy::requireGridCoverage );

//! Conventional ocean pole-tide loading, IERS Section 7.1.5, using the 2018 secular pole.
//! Grid origins are retained without conversion. The official CM grid already includes its geocentre correction.
class OceanPoleTideDeformation : public BodyDeformationModel
{
public:
    OceanPoleTideDeformation( const input_output::OceanPoleTideGrid& grid,
                              const std::shared_ptr< earth_orientation::PolarMotionCalculator > polarMotionCalculator,
                              const OceanPoleTideLatitudePolicy latitudePolicy = OceanPoleTideLatitudePolicy::requireGridCoverage,
                              const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > converter = nullptr );
    OceanPoleTideDeformation( const std::string& gridFile,
                              const std::shared_ptr< earth_orientation::PolarMotionCalculator > polarMotionCalculator,
                              const OceanPoleTideLatitudePolicy latitudePolicy = OceanPoleTideLatitudePolicy::requireGridCoverage,
                              const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > converter = nullptr );

    const input_output::OceanPoleTideGrid& getGrid( ) const
    {
        return grid_;
    }
    OceanPoleTideLatitudePolicy getLatitudePolicy( ) const
    {
        return latitudePolicy_;
    }

protected:
    Eigen::Vector3d evaluateDisplacement( const Time& tdb, const std::shared_ptr< ground_stations::GroundStationState > station ) override;
    Eigen::Vector6d evaluateDisplacementState( const Time& tdb,
                                               const std::shared_ptr< ground_stations::GroundStationState > station ) override;

private:
    Eigen::Vector3d calculateDisplacementAtTt( const Time& tt, const std::shared_ptr< ground_stations::GroundStationState > station );
    input_output::OceanPoleTideGrid grid_;
    std::shared_ptr< earth_orientation::PolarMotionCalculator > polarMotionCalculator_;
    OceanPoleTideLatitudePolicy latitudePolicy_;
};
}  // namespace basic_astrodynamics
}  // namespace tudat
#endif
