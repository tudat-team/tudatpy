#ifndef TUDAT_OCEAN_POLE_TIDE_DEFORMATION_H
#define TUDAT_OCEAN_POLE_TIDE_DEFORMATION_H

#include "tudat/io/readOceanPoleTideData.h"
#include "tudat/astro/ground_stations/bodyDeformationModel.h"
#include "tudat/astro/earth_orientation/polarMotionCalculator.h"

namespace tudat
{
namespace basic_astrodynamics
{
//! IERS 2010 Eq. (29), with radians as wobble input and metres as output (radial/north/east).
//! Defaults use Chapter 1 constants, not the older G in the published 2003-mean-pole test.
Eigen::Vector3d calculateOceanPoleTideLoading( const Eigen::Matrix< double, 3, 2 >& coefficients,
                                               const Eigen::Vector2d& wobbleRadians,
                                               const double scaleMetresPerRadian = 5340.428562229677 );

//! Bilinear interpolation in the model layer; longitude wraps only for grids spanning a complete period.
//! The official half-degree grid stops at +/-89.75 degrees; latitudes up to one grid step beyond the first/last row
//! use that row. Coordinates further outside the coverage (for example of a regional grid) raise std::out_of_range.
Eigen::Matrix< double, 3, 2 > interpolateOceanPoleTideCoefficients( const input_output::OceanPoleTideGrid& grid,
                                                                    const double latitudeDegrees,
                                                                    const double longitudeDegrees );

//! Conventional ocean pole-tide loading, IERS Section 7.1.5, using the 2018 secular pole.
//! Grid origins are retained without conversion. The official CM grid already includes its geocentre correction.
class OceanPoleTideDeformation : public BodyDeformationModel
{
public:
    OceanPoleTideDeformation( const input_output::OceanPoleTideGrid& grid,
                              const std::shared_ptr< earth_orientation::PolarMotionCalculator > polarMotionCalculator,
                              const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > converter = nullptr );

    const input_output::OceanPoleTideGrid& getGrid( ) const
    {
        return grid_;
    }

protected:
    Eigen::Vector3d evaluateDisplacement( const Time& tdb, const std::shared_ptr< ground_stations::GroundStationState > station ) override;
    Eigen::Vector6d evaluateDisplacementState( const Time& tdb,
                                               const std::shared_ptr< ground_stations::GroundStationState > station ) override;

private:
    Eigen::Vector3d calculateDisplacementAtTt( const Time& tt, const std::shared_ptr< ground_stations::GroundStationState > station );

    //! Grid of ocean pole tide coefficients.
    input_output::OceanPoleTideGrid grid_;

    //! Object computing the position of the pole, from which the pole tide is computed.
    std::shared_ptr< earth_orientation::PolarMotionCalculator > polarMotionCalculator_;
};
}  // namespace basic_astrodynamics
}  // namespace tudat
#endif
