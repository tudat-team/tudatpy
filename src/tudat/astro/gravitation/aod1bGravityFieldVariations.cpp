/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#include "tudat/astro/gravitation/aod1bGravityFieldVariations.h"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <stdexcept>
#include "tudat/astro/gravitation/tabulatedGravityFieldVariations.h"
#include "tudat/math/basic/legendrePolynomials.h"

namespace tudat
{
namespace gravitation
{
namespace
{

//! Keep the existing linear interpolation and change only its query time and allowed intervals.
class Aod1bGravityFieldVariations : public TabulatedGravityFieldVariations
{
public:
    //! Store scaled potential anomalies and the explicitly selected terrestrial time converter.
    Aod1bGravityFieldVariations( const std::map< double, Eigen::MatrixXd >& cosine,
                                 const std::map< double, Eigen::MatrixXd >& sine,
                                 const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter >& converter ):
        TabulatedGravityFieldVariations(
                cosine,
                sine,
                2,
                0,
                std::make_shared< interpolators::InterpolatorSettings >( interpolators::linear_interpolator,
                                                                         interpolators::huntingAlgorithm,
                                                                         false,
                                                                         interpolators::throw_exception_at_boundary ) ),
        converter_( converter )
    {
        for( const auto& entry : cosine )
        {
            gpsKnots_.push_back( entry.first );
        }
    }

    //! Convert TDB to GPS, accept exact converted knots and reject missing three-hour intervals.
    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > calculateSphericalHarmonicsCorrections( const double time ) override
    {
        if( !std::isfinite( time ) )
        {
            throw std::runtime_error( "AOD1B requires a finite TDB epoch" );
        }

        // GPS is continuous TAI minus 19 seconds; use the rotation's selected conversion to obtain TAI.
        double gps = static_cast< double >( converter_->getTaiTime( basic_astrodynamics::tdb_scale, Time( time ) ) - 19.0 );

        auto next = std::lower_bound( gpsKnots_.begin( ), gpsKnots_.end( ), gps );

        // A rounded double TDB knot can invert just outside coverage. Match only its exact forward image.
        for( int side = 0; side < 2; ++side )
        {
            const auto candidate = side == 0 ? next : ( next == gpsKnots_.begin( ) ? gpsKnots_.end( ) : std::prev( next ) );

            if( candidate != gpsKnots_.end( ) && std::abs( gps - *candidate ) < 1.0e-3 &&
                static_cast< double >( converter_->getCurrentTime< Time >(
                        basic_astrodynamics::tai_scale, basic_astrodynamics::tdb_scale, Time( *candidate ) + 19.0 ) ) == time )
            {
                gps = *candidate;
                next = candidate;
                break;
            }
        }
        if( next == gpsKnots_.end( ) || ( next == gpsKnots_.begin( ) && gps != *next ) )
        {
            throw std::runtime_error( "AOD1B query is outside supplied GPS coverage" );
        }
        if( gps != *next && *next - *std::prev( next ) != 10800.0 )
        {
            throw std::runtime_error( "AOD1B query falls in a missing three-hour interval" );
        }
        return TabulatedGravityFieldVariations::calculateSphericalHarmonicsCorrections( gps );
    }

private:
    //! Continuous source knots identify valid interpolation intervals without another coefficient copy.
    std::vector< double > gpsKnots_;

    //! Model-owned reference to the converter selected by the Earth rotation settings.
    std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > converter_;
};

}  // namespace

//! Convert one selected component to output-field units once, retaining interpolation in GPS seconds.
std::shared_ptr< GravityFieldVariations > createAod1bGravityFieldVariations(
        const std::vector< input_output::Aod1bGravityFieldData >& data,
        const input_output::Aod1bComponent component,
        const int maximumDegree,
        const int maximumOrder,
        const double outputGravitationalParameter,
        const double outputReferenceRadius,
        const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter >& converter )
{
    if( data.empty( ) || !converter || maximumDegree < 2 || maximumOrder < 0 || maximumOrder > maximumDegree ||
        !std::isfinite( outputGravitationalParameter ) || outputGravitationalParameter <= 0.0 || !std::isfinite( outputReferenceRadius ) ||
        outputReferenceRadius <= 0.0 )
    {
        throw std::runtime_error( "AOD1B requires RL07 data, an explicit converter, positive output constants and degree >= 2" );
    }

    // Build the selected output tables once; raw products are not retained by the runtime model.
    std::map< double, Eigen::MatrixXd > cosine, sine;

    const auto& first = data.front( );

    for( const auto& file : data )
    {
        input_output::validateAod1bGravityFieldData( file );
        if( file.release != first.release || file.softwareVersion != first.softwareVersion || file.maximumDegree != first.maximumDegree ||
            file.gravitationalParameter != first.gravitationalParameter || file.referenceRadius != first.referenceRadius ||
            file.fullyNormalized != first.fullyNormalized || maximumDegree > file.maximumDegree )
        {
            throw std::runtime_error( "AOD1B daily files have incompatible release, normalization, constants or degree" );
        }
        for( const auto& sample : file.samples )
        {
            if( sample.component != component )
            {
                continue;
            }

            // Include low degrees temporarily so the existing normalization utility uses physical indices.
            Eigen::MatrixXd c = Eigen::MatrixXd::Zero( maximumDegree + 1, maximumOrder + 1 );

            Eigen::MatrixXd s = c;

            for( int degree = 0; degree <= maximumDegree; ++degree )
            {
                for( int order = 0; order <= std::min( degree, maximumOrder ); ++order )
                {
                    const auto& values = sample.coefficients.at( { degree, order } ).values;

                    c( degree, order ) = values[ 0 ];
                    s( degree, order ) = values[ 1 ];
                }
            }
            if( !file.fullyNormalized )
            {
                const auto normalized = basic_mathematics::convertUnnormalizedToGeodesyNormalizedCoefficients( c, s );

                c = normalized.first;
                s = normalized.second;
            }
            for( int degree = 2; degree <= maximumDegree; ++degree )
            {
                const double scale = file.gravitationalParameter / outputGravitationalParameter *
                        std::pow( file.referenceRadius / outputReferenceRadius, degree );

                c.row( degree ) *= scale;
                s.row( degree ) *= scale;
            }
            if( !c.allFinite( ) || !s.allFinite( ) ||
                !cosine.emplace( sample.gpsSecondsSinceJ2000, c.bottomRows( maximumDegree - 1 ) ).second )
            {
                throw std::runtime_error( "AOD1B scaling is nonfinite or a selected component epoch is duplicated" );
            }
            sine.emplace( sample.gpsSecondsSinceJ2000, s.bottomRows( maximumDegree - 1 ) );
        }
    }
    if( cosine.size( ) < 2 )
    {
        throw std::runtime_error( "AOD1B interpolation requires at least two epochs of the selected component" );
    }
    return std::make_shared< Aod1bGravityFieldVariations >( cosine, sine, converter );
}

}  // namespace gravitation
}  // namespace tudat
