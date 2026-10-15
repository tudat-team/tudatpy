#ifndef TUDAT_EARTH_ORIENTATION_DATA_EVALUATOR_H
#define TUDAT_EARTH_ORIENTATION_DATA_EVALUATOR_H

#include <memory>
#include <vector>
#include <Eigen/Core>
#include "tudat/interface/sofa/sofaTimeConversions.h"
#include "tudat/basics/timeType.h"
#include "tudat/io/readEarthOrientationData.h"

namespace tudat
{
namespace earth_orientation
{

//! Method used to interpolate the daily Earth orientation parameters (EOP) between two data points.
enum class EarthOrientationInterpolation {
    //! Linear interpolation between the two data points around the requested epoch.
    linear,

    //! Lagrange interpolation through four data points: two on each side of the requested epoch, or the first or last four
    //! points of the data near its start or end.
    lagrange_four_point
};

//! What to do when celestial pole offsets (dX and dY) are needed at an epoch for which the EOP data does not contain them.
enum class MissingCelestialPolePolicy {
    //! Throw an exception.
    reject,

    //! Use zero offsets (and zero offset rates).
    zero
};

//! What to do when UT1 is needed at an epoch before the start of the EOP data.
enum class HistoricalUt1Policy {
    //! Throw an exception.
    reject,

    //! Compute UT1 from the table of historical values of Delta T (TT - UT1). This is only done before 1962, when the
    //! EOP data starts. Polar motion and celestial pole offsets are not available at these epochs.
    delta_t
};

//! Settings for the interpolation of the daily EOP data by an EarthOrientationDataEvaluator. The file readers only read
//! the data; these settings determine how it is used.
struct EarthOrientationEvaluationSettings {
    //! Largest time interval (in seconds) between two consecutive data points over which the data is interpolated. An epoch
    //! for which interpolation would use a longer interval throws an exception. Default: 36 hours, so a single missing day in
    //! daily data is rejected.
    double maximumGapSeconds = 36.0 * 3600.0;

    //! What to do when the data has no celestial pole offsets at the requested epoch.
    MissingCelestialPolePolicy missingCelestialPole = MissingCelestialPolePolicy::reject;

    //! What to do when UT1 is requested before the start of the data.
    HistoricalUt1Policy historicalUt1 = HistoricalUt1Policy::delta_t;

    //! Interpolation method.
    EarthOrientationInterpolation interpolation = EarthOrientationInterpolation::linear;
};

//! Daily EOP values interpolated to a single epoch, as returned by EarthOrientationDataEvaluator::evaluate.
/*!
 *  These are the values from the data file only, interpolated to the requested epoch. They do not include the short-period
 *  (tidal and libration) corrections, or the IAU precession-nutation model; those are added by the PolarMotionCalculator,
 *  PrecessionNutationCalculator and TerrestrialTimeScaleConverter to give the rotation angles in EarthOrientationEvaluatedAnglesAndRates.
 *  Each rate is the time derivative of the interpolated value, per TAI second. A rate is empty when it cannot be computed,
 *  which is the case at a data point next to a gap that is longer than the maximum interpolation gap.
 */
struct EarthOrientationInputDataValue {
    //! Polar motion x_p and y_p (radians).
    Eigen::Vector2d polarMotion;

    //! Rate of change of polar motion (radians per second).
    std::optional< Eigen::Vector2d > polarMotionRate;

    //! UT1 - TAI (seconds). Unlike UT1 - UTC, this does not jump at leap seconds.
    double ut1MinusTai;

    //! Rate of change of UT1 - TAI (seconds per second).
    std::optional< double > ut1MinusTaiRate;

    //! Celestial pole offsets dX and dY (radians); empty if the data does not contain them at this epoch.
    std::optional< Eigen::Vector2d > celestialPoleOffsets;

    //! Rate of change of the celestial pole offsets (radians per second).
    std::optional< Eigen::Vector2d > celestialPoleOffsetRate;

    //! Excess length of day (seconds); empty if the data does not contain it at this epoch.
    std::optional< double > lengthOfDay;
};

//! Interpolates the daily EOP data (polar motion, UT1, celestial pole offsets and length of day) to a given epoch.
/*!
 *  The EOP files give their values at midnight UTC on each day. This class converts these epochs to TAI, and converts
 *  UT1 - UTC to UT1 - TAI. The interpolation is then done in TAI, so it is not affected by the jumps of UTC at leap seconds.
 *  The values are interpolated linearly or with a four-point Lagrange interpolation (see EarthOrientationEvaluationSettings),
 *  and the time derivative of the interpolating function is returned with each value.
 *
 *  Epochs outside the data, and epochs in an interval between data points that is longer than the maximum gap, throw an
 *  exception. An epoch within fractional-hour clock roundoff of a data point is evaluated at that point. Its value is taken
 *  from the file and its rate from the interval after that point (or before it, at the last data point). The rate is not
 *  available if that interval is longer than the maximum gap. The recognition margin depends
 *  on long-double precision, not on the number of seconds since J2000.
 */
class EarthOrientationDataEvaluator
{
public:
    //! Constructor; stores a copy of the data and converts its epochs to TAI.
    EarthOrientationDataEvaluator( const input_output::EarthOrientationData& data,
                                   const EarthOrientationEvaluationSettings& settings = {} );

    //! Margin (in seconds) by which a requested TAI epoch may lie outside the first or last observation.
    //! Such an epoch is evaluated at that observation. The margin covers the rounding of an observation
    //! epoch to double (about 1e-7 s), and changes UT1 by less than 1e-13 s.
    static constexpr double coverageToleranceSeconds = 1.0E-6;

    //! Interpolated EOP values and their rates at the given TAI epoch (seconds since J2000). If requireCelestialPole is true
    //! and the data has no celestial pole offsets at this epoch, the missingCelestialPole setting is applied.
    //! Epochs within fractional-hour rounding error of a data point use that point before selecting the interval.
    EarthOrientationInputDataValue evaluate( const Time& tai, bool requireCelestialPole = false ) const;

    //! Celestial pole offsets dX and dY (radians) at the given TAI epoch, applying the missingCelestialPole setting.
    Eigen::Vector2d getCelestialPoleOffsets( const Time& tai ) const;

    //! Epochs of the data points in TAI (seconds since J2000).
    const std::vector< Time >& getEpochsTai( ) const
    {
        return epochs_;
    }
    const Time& getStartTai( ) const
    {
        return epochs_.front( );
    }
    const Time& getEndTai( ) const
    {
        return epochs_.back( );
    }
    const input_output::EarthOrientationData& getRawData( ) const
    {
        return data_;
    }
    const EarthOrientationEvaluationSettings& getSettings( ) const
    {
        return settings_;
    }

    //! Convert TAI to UTC (both in seconds since J2000). Throws during a leap second (23:59:60), which cannot be expressed
    //! as a number of UTC seconds since J2000.
    static Time convertTaiToScalarUtc( const Time& tai );

private:
    //! The EOP data, as read from the file.
    input_output::EarthOrientationData data_;

    //! Settings for the interpolation of the data.
    EarthOrientationEvaluationSettings settings_;

    //! Epochs of the data points in TAI (seconds since J2000).
    std::vector< Time > epochs_;

    //! UTC modified Julian day of each data point, used to look up its record in data_.
    std::vector< double > recordMjd_;

    //! UT1 - TAI (seconds) at each data point.
    std::vector< double > ut1MinusTai_;

    //! Indices of the data points used to interpolate between data points left and left + 1 (two for linear, four for
    //! Lagrange interpolation).
    std::vector< std::size_t > stencil( std::size_t left ) const;

    //! Whether all intervals between the given consecutive data points are at most the maximum gap.
    bool validStencil( const std::vector< std::size_t >& indices ) const;

    //! Throw if four-point interpolation of UT1 - TAI would make UT1 decrease with increasing TAI somewhere.
    void validateCubicUt1Monotonicity( ) const;
};

}  // namespace earth_orientation
}  // namespace tudat
#endif
