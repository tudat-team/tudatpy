#ifndef TUDAT_HIGH_FREQUENCY_EARTH_ORIENTATION_MODEL_H
#define TUDAT_HIGH_FREQUENCY_EARTH_ORIENTATION_MODEL_H

#include "tudat/astro/earth_orientation/shortPeriodEarthOrientationCorrectionCalculator.h"

namespace tudat
{
namespace earth_orientation
{

//! Select only the ocean-tide contribution; independent IERS libration is always retained.
enum class OceanTideEopModel { iers_2010, desai_sibois_2016 };

struct HighFrequencyEopCorrectionCalculators {
    std::shared_ptr< ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > > polarMotion;
    std::shared_ptr< ShortPeriodEarthOrientationCorrectionCalculator< double > > ut1;
};

//! Compose one ocean model and the conventional quasi-diurnal PM/semidiurnal UT1 libration.
//! Cutoffs use the original amplitude units (microarcseconds for PM, microseconds for UT1).
//! Daily EOP, CPO, precession/nutation and long-period zonal tides are not included here.
HighFrequencyEopCorrectionCalculators createHighFrequencyEopCorrectionCalculators(
        OceanTideEopModel oceanModel = OceanTideEopModel::iers_2010,
        double minimumPolarMotionAmplitude = 0.0,
        double minimumUt1Amplitude = 0.0 );

}  // namespace earth_orientation
}  // namespace tudat
#endif
