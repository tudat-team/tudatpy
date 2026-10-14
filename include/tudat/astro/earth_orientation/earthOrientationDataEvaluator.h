#ifndef TUDAT_EARTH_ORIENTATION_DATA_EVALUATOR_H
#define TUDAT_EARTH_ORIENTATION_DATA_EVALUATOR_H

#include <algorithm>
#include <memory>
#include <vector>
#include <type_traits>
#include "tudat/math/interpolators/oneDimensionalInterpolator.h"
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

//! Strict, continuous TAI-axis evaluation of raw UTC-indexed EOP data.
//! Linear or explicitly selected four-point Lagrange interpolation.
//! No artificial buffer, extrapolation, order reduction or zero fallback.
//! At knots, use the interval to the right (left at the final sample). The
//! four-point window centers that interval, clamped to the first/last four samples.
//! Rates are absent if any selected-stencil gap is too long; raw values remain valid.
class EarthOrientationDataEvaluator
{
public:
    EarthOrientationDataEvaluator( const input_output::EarthOrientationData& data,
                                   const EarthOrientationEvaluationSettings& settings = {} );
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
    input_output::EarthOrientationData data_;
    EarthOrientationEvaluationSettings settings_;
    std::vector< Time > epochs_;
    std::vector< input_output::EarthOrientationRecord > records_;
    std::vector< double > ut1MinusTai_;
    std::vector< std::size_t > stencil( std::size_t left ) const;
    bool validStencil( const std::vector< std::size_t >& indices ) const;
    void validateCubicUt1Monotonicity( ) const;
};

//! Guarded UTC compatibility view. The production consumers use TAI directly.
enum class EarthOrientationDataField { polar_motion, celestial_pole, ut1_minus_utc };
template< typename Value >
class EarthOrientationDataInterpolator : public interpolators::OneDimensionalInterpolator< double, Value >
{
public:
    EarthOrientationDataInterpolator( const std::shared_ptr< EarthOrientationDataEvaluator >& evaluator,
                                      const EarthOrientationDataField field ):
        interpolators::OneDimensionalInterpolator< double, Value >( interpolators::throw_exception_at_boundary ), evaluator_( evaluator ),
        field_( field )
    {
        for( const auto& entry : evaluator->getRawData( ).records )
        {
            const double utc = ( entry.first - 51544.5 ) * 86400.0;
            this->independentValues_.push_back( utc );
            if constexpr( std::is_same< Value, double >::value )
            {
                this->dependentValues_.push_back( entry.second.ut1MinusUtc );
                atomicOffsets_.push_back( static_cast< double >( sofa_interface::convertUTCtoTAI( Time( utc ) ) - Time( utc ) ) );
            }
            else
            {
                if( field == EarthOrientationDataField::celestial_pole )
                    this->dependentValues_.push_back( entry.second.celestialPoleOffsets
                                                              ? Eigen::Vector2d( ( *entry.second.celestialPoleOffsets )[ 0 ],
                                                                                 ( *entry.second.celestialPoleOffsets )[ 1 ] )
                                                              : Eigen::Vector2d::Constant( TUDAT_NAN ).eval( ) );
                else
                    this->dependentValues_.push_back( Eigen::Vector2d( entry.second.polarMotion[ 0 ], entry.second.polarMotion[ 1 ] ) );
            }
        }
    }
    Value interpolate( const double utc ) override
    {
        const Time tai = sofa_interface::convertUTCtoTAI( Time( utc ) );
        if constexpr( std::is_same< Value, double >::value )
        {
            // Preserve the exact raw value at a sample; subtracting and adding
            // the atomic offset would otherwise lose bits in UT1-UTC.
            const auto node = std::lower_bound( this->independentValues_.begin( ), this->independentValues_.end( ), utc );
            if( node != this->independentValues_.end( ) && *node == utc )
            {
                return this->dependentValues_[ node - this->independentValues_.begin( ) ];
            }
            evaluator_->evaluate( tai );  // Apply the same strict coverage and gap checks.
            const auto& epochs = evaluator_->getEpochsTai( );
            const auto upper = std::upper_bound( epochs.begin( ), epochs.end( ), tai );
            const std::size_t right = upper == epochs.end( ) ? epochs.size( ) - 1 : upper - epochs.begin( );
            const std::size_t left = right - 1;
            const double fraction =
                    static_cast< double >( tai - epochs[ left ] ) / static_cast< double >( epochs[ right ] - epochs[ left ] );
            // Difference the small quantities first to avoid cancelling a large
            // DAT against UT1-TAI. This is algebraically the native TAI interpolation.
            if( evaluator_->getSettings( ).interpolation == EarthOrientationInterpolation::lagrange_four_point )
                return evaluator_->evaluate( tai ).ut1MinusTai + static_cast< double >( tai - Time( utc ) );
            return this->dependentValues_[ left ] +
                    fraction *
                    ( this->dependentValues_[ right ] - this->dependentValues_[ left ] -
                      ( atomicOffsets_[ right ] - atomicOffsets_[ left ] ) ) +
                    ( static_cast< double >( tai - Time( utc ) ) - atomicOffsets_[ left ] );
        }
        else
        {
            if( field_ == EarthOrientationDataField::celestial_pole ) return evaluator_->getCelestialPoleOffsets( tai );
            return evaluator_->evaluate( tai ).polarMotion;
        }
    }
    interpolators::InterpolatorTypes getInterpolatorType( ) override
    {
        return evaluator_->getSettings( ).interpolation == EarthOrientationInterpolation::linear ? interpolators::linear_interpolator
                                                                                                 : interpolators::lagrange_interpolator;
    }

private:
    std::shared_ptr< EarthOrientationDataEvaluator > evaluator_;
    EarthOrientationDataField field_;
    std::vector< double > atomicOffsets_;
};

}  // namespace earth_orientation
}  // namespace tudat
#endif
