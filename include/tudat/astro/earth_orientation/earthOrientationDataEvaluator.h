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

enum class EarthOrientationInterpolation { linear, lagrange_four_point };
enum class MissingCelestialPolePolicy { reject, zero };
enum class HistoricalUt1Policy { reject, delta_t };

//! Model policies; raw readers do not interpolate or impose these choices.
struct EarthOrientationEvaluationSettings {
    double maximumGapSeconds = 36.0 * 3600.0;
    MissingCelestialPolePolicy missingCelestialPole = MissingCelestialPolePolicy::reject;
    HistoricalUt1Policy historicalUt1 = HistoricalUt1Policy::reject;
    EarthOrientationInterpolation interpolation = EarthOrientationInterpolation::linear;
};

struct EarthOrientationDataValue {
    Eigen::Vector2d polarMotion;
    std::optional< Eigen::Vector2d > polarMotionRate;
    double ut1MinusTai;
    std::optional< double > ut1MinusTaiRate;
    std::optional< Eigen::Vector2d > celestialPoleOffsets;
    std::optional< Eigen::Vector2d > celestialPoleOffsetRate;
    std::optional< double > lengthOfDay;
};

//! Interpolate daily EOP values and their time derivatives using TAI seconds.
//! Convert UTC sample times to TAI and UT1-UTC values to continuous UT1-TAI.
//! Use linear or four-point interpolation within the supplied observations only.
//! At a sample, take the derivative of the interval on its right, except at the
//! last sample. Near the beginning/end, use the first/last four observations.
//! An excessive gap prevents interpolation between samples; a raw sample value
//! remains available there, but its derivative is unavailable.
class EarthOrientationDataEvaluator
{
public:
    //! Own the supplied raw observations and prepare their continuous TAI times.
    EarthOrientationDataEvaluator( const input_output::EarthOrientationData& data,
                                   const EarthOrientationEvaluationSettings& settings = {} );

    //! Return observed values and the selected interpolant's rates with respect to TAI seconds.
    EarthOrientationDataValue evaluate( const Time& tai, bool requireCelestialPole = false ) const;

    Eigen::Vector2d getCelestialPoleOffsets( const Time& tai ) const;
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
    //! Preserve precise TT; only exact rounded double endpoint images denote raw boundaries.
    Time convertTtToTai( const double tt ) const;
    Time convertTtToTai( const Time& tt ) const;

    //! Scalar UTC cannot uniquely represent an inserted second; throw for it.
    static Time convertTaiToScalarUtc( const Time& tai );

private:
    //! The sole owned set of raw observations and their source metadata.
    input_output::EarthOrientationData data_;

    //! Interpolation, gap and missing-observation policies for these data.
    EarthOrientationEvaluationSettings settings_;

    //! Continuous TAI observation times and their matching raw UTC MJD keys.
    std::vector< Time > epochs_;

    std::vector< double > recordMjd_;

    //! Continuous UT1-TAI seconds at each observation, excluding UTC jumps.
    std::vector< double > ut1MinusTai_;

    std::vector< std::size_t > stencil( std::size_t left ) const;
    bool validStencil( const std::vector< std::size_t >& indices ) const;
    void validateCubicUt1Monotonicity( ) const;
};

}  // namespace earth_orientation
}  // namespace tudat
#endif
