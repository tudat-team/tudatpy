#include "tudat/astro/earth_orientation/earthOrientationDataEvaluator.h"
#include "tudat/interface/sofa/sofaTimeConversions.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <type_traits>

namespace tudat
{
namespace earth_orientation
{

// Own the raw observations once and prepare scalar indexes without borrowed map iterators.
EarthOrientationDataEvaluator::EarthOrientationDataEvaluator( const input_output::EarthOrientationData& data,
                                                              const EarthOrientationEvaluationSettings& settings ):
    data_( data ), settings_( settings )
{
    input_output::validateEarthOrientationData( data );
    if( data.records.size( ) < 2 )
    {
        throw std::runtime_error( "EOP evaluation requires at least two UTC epochs" );
    }
    if( !std::isfinite( settings.maximumGapSeconds ) || settings.maximumGapSeconds <= 0 )
    {
        throw std::runtime_error( "EOP maximum interpolation gap must be positive and finite" );
    }
    if( settings.interpolation != EarthOrientationInterpolation::linear &&
        settings.interpolation != EarthOrientationInterpolation::lagrange_four_point )
    {
        throw std::runtime_error( "Unsupported EOP interpolation policy" );
    }
    if( settings.interpolation == EarthOrientationInterpolation::lagrange_four_point && data.records.size( ) < 4 )
    {
        throw std::runtime_error( "Four-point EOP interpolation requires at least four raw epochs" );
    }
    for( const auto& entry : data.records )
    {
        // Convert each UTC MJD sample to seconds since J2000 on the continuous TAI clock.
        const Time utc = ( entry.first - 51544.5 ) * 86400.0;

        const Time tai = sofa_interface::convertUTCtoTAI( utc );

        if( !epochs_.empty( ) && tai <= epochs_.back( ) )
        {
            throw std::runtime_error( data.source + ": TAI sample epochs are not strictly increasing" );
        }
        epochs_.push_back( tai );
        recordMjd_.push_back( entry.first );
        ut1MinusTai_.push_back( entry.second.ut1MinusUtc - static_cast< double >( tai - utc ) );
        if( epochs_.size( ) > 1 && tai + ut1MinusTai_.back( ) <= epochs_[ epochs_.size( ) - 2 ] + ut1MinusTai_[ ut1MinusTai_.size( ) - 2 ] )
        {
            throw std::runtime_error( data.source + ": UT1 sample epochs are not monotonic" );
        }
    }
    if( settings.interpolation == EarthOrientationInterpolation::lagrange_four_point )
    {
        validateCubicUt1Monotonicity( );
    }
}

std::vector< std::size_t > EarthOrientationDataEvaluator::stencil( const std::size_t left ) const
{
    if( settings_.interpolation == EarthOrientationInterpolation::linear ) return { left, left + 1 };
    const std::size_t first = std::min( left == 0 ? 0 : left - 1, epochs_.size( ) - 4 );
    return { first, first + 1, first + 2, first + 3 };
}

bool EarthOrientationDataEvaluator::validStencil( const std::vector< std::size_t >& indices ) const
{
    for( std::size_t i = 1; i < indices.size( ); ++i )
        if( static_cast< double >( epochs_[ indices[ i ] ] - epochs_[ indices[ i - 1 ] ] ) > settings_.maximumGapSeconds ) return false;
    return true;
}

void EarthOrientationDataEvaluator::validateCubicUt1Monotonicity( ) const
{
    for( std::size_t left = 0; left + 1 < epochs_.size( ); ++left )
    {
        const auto indices = stencil( left );
        if( !validStencil( indices ) ) continue;  // Interior evaluation is rejected for this complete stencil.
        const double duration = static_cast< double >( epochs_[ left + 1 ] - epochs_[ left ] );
        std::array< double, 4 > coefficient{ 0.0, 0.0, 0.0, 0.0 };
        for( const std::size_t node : indices )
        {
            const double x = static_cast< double >( epochs_[ node ] - epochs_[ left ] ) / duration;
            std::array< double, 4 > basis{ 1.0, 0.0, 0.0, 0.0 };
            int degree = 0;
            for( const std::size_t other : indices )
            {
                if( other == node ) continue;
                const double y = static_cast< double >( epochs_[ other ] - epochs_[ left ] ) / duration;
                std::array< double, 4 > next{ 0.0, 0.0, 0.0, 0.0 };
                for( int k = 0; k <= degree; ++k )
                {
                    next[ k ] -= y * basis[ k ] / ( x - y );
                    next[ k + 1 ] += basis[ k ] / ( x - y );
                }
                basis = next;
                ++degree;
            }
            for( int k = 0; k < 4; ++k ) coefficient[ k ] += ( ut1MinusTai_[ node ] - ut1MinusTai_[ left ] ) * basis[ k ];
        }
        const auto derivative = [ &coefficient, duration ]( const double u ) {
            return 1.0 + ( coefficient[ 1 ] + 2.0 * coefficient[ 2 ] * u + 3.0 * coefficient[ 3 ] * u * u ) / duration;
        };
        double minimum = std::min( derivative( 0.0 ), derivative( 1.0 ) );
        if( coefficient[ 3 ] != 0.0 )
        {
            const double extremum = -coefficient[ 2 ] / ( 3.0 * coefficient[ 3 ] );
            if( extremum > 0.0 && extremum < 1.0 ) minimum = std::min( minimum, derivative( extremum ) );
        }
        if( !std::isfinite( minimum ) || minimum <= 0.0 )
            throw std::runtime_error( data_.source + ": four-point interpolation gives non-positive dUT1/dTAI inside an interval" );
    }
}

// Evaluate values and selected-interpolant rates from the same two or four observations.
EarthOrientationDataValue EarthOrientationDataEvaluator::evaluate( const Time& requestedTai, const bool requireCelestialPole ) const
{
    if( !std::isfinite( static_cast< double >( requestedTai ) ) || requestedTai < epochs_.front( ) - coverageToleranceSeconds ||
        requestedTai > epochs_.back( ) + coverageToleranceSeconds )
    {
        throw std::runtime_error( data_.source + ": requested epoch outside raw EOP coverage" );
    }

    // An epoch within the coverage tolerance outside the first or last observation (for instance the
    // rounded double value of that observation's epoch) is evaluated at that observation.
    const Time tai = std::min( std::max( requestedTai, epochs_.front( ) ), epochs_.back( ) );

    // Select the right-hand interval at an exact sample, except at the final sample.
    const auto upper = std::upper_bound( epochs_.begin( ), epochs_.end( ), tai );

    const std::size_t right = upper == epochs_.end( ) ? epochs_.size( ) - 1 : upper - epochs_.begin( );

    const std::size_t left = right - 1;

    const std::optional< std::size_t > exactNode = tai == epochs_[ left ]
            ? std::optional< std::size_t >( left )
            : ( tai == epochs_[ right ] ? std::optional< std::size_t >( right ) : std::nullopt );

    const auto indices = stencil( left );

    const bool valid = validStencil( indices );

    if( !valid && !exactNode )
    {
        throw std::runtime_error( data_.source + ": requested epoch has a gap in its complete EOP stencil" );
    }

    // Exact raw values remain available across a gap, but no interpolated rates are assigned.
    EarthOrientationDataValue value;

    if( valid )
    {
        // Borrow only the selected records for this call; the owned map is never modified.
        std::array< const input_output::EarthOrientationRecord*, 4 > selectedRecords{};

        for( std::size_t j = 0; j < indices.size( ); ++j )
        {
            selectedRecords[ j ] = &data_.records.at( recordMjd_[ indices[ j ] ] );
        }

        // Scale interpolation abscissae by this interval's duration in TAI seconds.
        const double duration = static_cast< double >( epochs_[ right ] - epochs_[ left ] );

        const double u = static_cast< double >( tai - epochs_[ left ] ) / duration;

        // Node abscissae, value weights and derivative weights in inverse TAI seconds.
        std::array< double, 4 > x{}, weight{}, rate{};

        for( std::size_t j = 0; j < indices.size( ); ++j )
        {
            x[ j ] = static_cast< double >( epochs_[ indices[ j ] ] - epochs_[ left ] ) / duration;
        }
        for( std::size_t j = 0; j < indices.size( ); ++j )
        {
            double w = 1.0, derivative = 0.0;

            for( std::size_t k = 0; k < indices.size( ); ++k )
            {
                if( k != j )
                {
                    w *= ( u - x[ k ] ) / ( x[ j ] - x[ k ] );
                }
            }
            for( std::size_t k = 0; k < indices.size( ); ++k )
            {
                if( k == j )
                {
                    continue;
                }
                double term = 1.0 / ( x[ j ] - x[ k ] );

                for( std::size_t m = 0; m < indices.size( ); ++m )
                {
                    if( m != j && m != k )
                    {
                        term *= ( u - x[ m ] ) / ( x[ j ] - x[ m ] );
                    }
                }
                derivative += term;
            }
            weight[ j ] = w;
            rate[ j ] = derivative / duration;
        }
        // Sum differences from the first observation to avoid cancellation of large UTC offsets.
        // The same sum is used for scalar (double) and vector (Eigen::Vector2d) quantities.
        const auto interpolate = [ &indices, &weight, &rate ]( const auto& get ) {
            using Value = std::decay_t< decltype( get( 0 ) ) >;

            const Value anchor = get( 0 );

            Value result = anchor;

            Value derivative = anchor * 0.0;

            for( std::size_t j = 1; j < indices.size( ); ++j )
            {
                const Value difference = get( j ) - anchor;

                result += weight[ j ] * difference;
                derivative += rate[ j ] * difference;
            }
            return std::make_pair( result, derivative );
        };

        const auto ut1 = interpolate( [ this, &indices ]( std::size_t sample ) { return ut1MinusTai_[ indices[ sample ] ]; } );

        value.ut1MinusTai = ut1.first;
        value.ut1MinusTaiRate = ut1.second;
        // Polar-motion values are radians and rates are radians per TAI second.
        const auto pm = interpolate( [ &selectedRecords ]( std::size_t sample ) {
            return Eigen::Vector2d( selectedRecords[ sample ]->polarMotion[ 0 ], selectedRecords[ sample ]->polarMotion[ 1 ] );
        } );

        value.polarMotion = pm.first;
        value.polarMotionRate = pm.second;
        // CPO interpolation requires observed corrections at every selected sample.
        const bool allCpo = std::all_of( selectedRecords.begin( ), selectedRecords.begin( ) + indices.size( ), []( const auto* record ) {
            return record->celestialPoleOffsets.has_value( );
        } );

        if( allCpo )
        {
            // Interpolate the measured celestial-pole offsets and their TAI rates together.
            const auto cpo = interpolate( [ &selectedRecords ]( std::size_t sample ) {
                return Eigen::Vector2d( ( *selectedRecords[ sample ]->celestialPoleOffsets )[ 0 ],
                                        ( *selectedRecords[ sample ]->celestialPoleOffsets )[ 1 ] );
            } );

            value.celestialPoleOffsets = cpo.first;
            value.celestialPoleOffsetRate = cpo.second;
        }
        // LOD remains optional and is reported only when the complete stencil supplies it.
        const bool allLod = std::all_of( selectedRecords.begin( ), selectedRecords.begin( ) + indices.size( ), []( const auto* record ) {
            return record->lengthOfDay.has_value( );
        } );

        if( allLod )
        {
            value.lengthOfDay =
                    interpolate( [ &selectedRecords ]( std::size_t sample ) { return *selectedRecords[ sample ]->lengthOfDay; } ).first;
        }
    }
    if( exactNode )
    {
        // Recover exact published sample values independently of interpolated rates.
        const auto& record = data_.records.at( recordMjd_[ *exactNode ] );

        value.polarMotion = Eigen::Vector2d( record.polarMotion[ 0 ], record.polarMotion[ 1 ] );
        value.ut1MinusTai = ut1MinusTai_[ *exactNode ];
        value.celestialPoleOffsets = record.celestialPoleOffsets
                ? std::optional< Eigen::Vector2d >(
                          Eigen::Vector2d( ( *record.celestialPoleOffsets )[ 0 ], ( *record.celestialPoleOffsets )[ 1 ] ) )
                : std::nullopt;
        value.lengthOfDay = record.lengthOfDay;
    }
    if( requireCelestialPole && !value.celestialPoleOffsets )
    {
        if( settings_.missingCelestialPole == MissingCelestialPolePolicy::zero )
        {
            value.celestialPoleOffsets = Eigen::Vector2d::Zero( );
            value.celestialPoleOffsetRate = Eigen::Vector2d::Zero( );
        }
        else
        {
            throw std::runtime_error( data_.source + ": celestial pole offsets missing in EOP interpolation stencil" );
        }
    }
    return value;
}

Eigen::Vector2d EarthOrientationDataEvaluator::getCelestialPoleOffsets( const Time& tai ) const
{
    return *evaluate( tai, true ).celestialPoleOffsets;
}

Time EarthOrientationDataEvaluator::convertTaiToScalarUtc( const Time& tai )
{
    // SOFA's UTC quasi-Julian date identifies the civil inserted second. It is
    // used only for this explicit UTC request, never for UT1 or rotation.
    const double taiDay = 2451545.0 + tai.fullDaysSinceEpoch( );
    const double taiFraction = static_cast< double >( tai.secondsSinceNoon( ) / 86400.0L );
    double utc1, utc2;
    if( iauTaiutc( taiDay, taiFraction, &utc1, &utc2 ) < 0 ) throw std::runtime_error( "Cannot convert TAI to UTC" );
    // Identify the civil day with SOFA, then compare the inserted-second
    // interval using Time arithmetic. A rounded civil clock alone can classify
    // a sub-nanosecond instant just before midnight as the following day.
    int year, month, day;
    double fraction, mjdOrigin, utcDay;
    if( iauJd2cal( utc1, utc2, &year, &month, &day, &fraction ) != 0 || iauCal2jd( year, month, day, &mjdOrigin, &utcDay ) != 0 )
        throw std::runtime_error( "Cannot determine UTC civil day" );
    for( int dayOffset : { 0, 1 } )
    {
        int y, m, d;
        double f, datBefore, datAfter;
        iauJd2cal( mjdOrigin, utcDay + dayOffset - 1.0, &y, &m, &d, &f );
        if( iauDat( y, m, d, 0.0, &datBefore ) < 0 ) throw std::runtime_error( "Cannot determine UTC offset" );
        iauJd2cal( mjdOrigin, utcDay + dayOffset, &y, &m, &d, &f );
        if( iauDat( y, m, d, 0.0, &datAfter ) < 0 ) throw std::runtime_error( "Cannot determine UTC offset" );
        if( y >= 1972 && datAfter - datBefore > 0.5 )
        {
            const Time midnightTai = Time( ( utcDay + dayOffset - 51544.5 ) * 86400.0 ) + datAfter;
            if( tai >= midnightTai - ( datAfter - datBefore ) && tai < midnightTai )
                throw std::runtime_error( "Scalar UTC cannot uniquely represent 23:59:60; use TAI/TT for this instant" );
        }
    }
    return sofa_interface::convertTAItoUTC( tai );
}

}  // namespace earth_orientation
}  // namespace tudat
