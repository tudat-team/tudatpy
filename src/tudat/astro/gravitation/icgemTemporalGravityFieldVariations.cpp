/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#include "tudat/astro/gravitation/icgemTemporalGravityFieldVariations.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include "tudat/astro/basic_astro/dateTime.h"
#include "tudat/math/basic/legendrePolynomials.h"
#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{
namespace gravitation
{
namespace
{

//! Keep the integer year separate from its fraction to avoid subtracting large decimal-year values.
std::pair< int, long double > calendarYear( const Time& epoch )
{
    // Shift the existing noon-based day count to the civil midnight used by provider dates.
    int day, month, year;

    basic_astrodynamics::convertShiftedJulianDayToCalendarDate(
            ( epoch + 43200.0 ).fullDaysSinceEpoch( ) + basic_astrodynamics::JULIAN_DAY_ON_J2000_INT, day, month, year );

    const Time start = basic_astrodynamics::timeFromDecomposedDateTime< Time >( year, 1, 1, 0, 0, 0.0L );

    const Time end = basic_astrodynamics::timeFromDecomposedDateTime< Time >( year + 1, 1, 1, 0, 0, 0.0L );

    return { year, static_cast< long double >( epoch - start ) / static_cast< long double >( end - start ) };
}

//! Resolve one printed date in the declared source scale; irregular suffixes require an explicit interpretation.
Time sourceDate( const std::string& token, const std::map< std::string, double >& overrides )
{
    const auto override = overrides.find( token );

    if( override != overrides.end( ) )
    {
        return Time( override->second );
    }

    // The neutral reader already validates Gregorian dates and the format-specific token length.
    const int hour = token.size( ) == 13 ? std::stoi( token.substr( 9, 2 ) ) : 0;

    const int minute = token.size( ) == 13 ? std::stoi( token.substr( 11, 2 ) ) : 0;

    if( hour > 23 || minute > 59 )
    {
        throw std::runtime_error( "ICGEM date requires a valid HHMM suffix or an explicit date override: " + token );
    }
    return basic_astrodynamics::timeFromDecomposedDateTime< Time >(
            std::stoi( token.substr( 0, 4 ) ), std::stoi( token.substr( 4, 2 ) ), std::stoi( token.substr( 6, 2 ) ), hour, minute, 0.0L );
}

//! Store only prepared numerical terms, rather than copying the raw provider product into the runtime model.
struct TemporalTerm {
    //! Original term keyword distinguishes reference, trend and sine/cosine amplitudes.
    std::string type;

    //! TAI reference and exclusive validity end; 1.0 terms have unbounded validity.
    Time origin, end;

    bool bounded = false;

    //! Source-calendar reference coordinate and period in the explicitly selected definition of years.
    std::pair< int, long double > calendarOrigin;

    double period = 0.0;

    //! Fully normalized output-field values; reference terms already have the adopted nominal value removed.
    double cosine = 0.0, sine = 0.0;
};

//! Keep each coefficient's base pieces and independent trend/periodic intervals in their original physical order.
struct TemporalCoefficient {
    //! Physical degree/order locate the coefficient inside the degree-two correction block.
    int degree = 0, order = 0;

    //! Base pieces are sorted by TAI origin; independent terms are filtered by their own validity.
    std::vector< TemporalTerm > bases, terms;
};

//! Evaluate prepared coefficient intervals directly, retaining the existing additive gravity update interface.
class IcgemTemporalGravityFieldVariations : public GravityFieldVariations
{
public:
    //! Prepare each selected coefficient once using the adopted nominal field and the selected time converter.
    IcgemTemporalGravityFieldVariations( const input_output::IcgemGravityFieldData& data,
                                         const basic_astrodynamics::TimeScales sourceScale,
                                         const IcgemYearConvention convention,
                                         const int maximumDegree,
                                         const int maximumOrder,
                                         const std::map< std::string, double >& overrides,
                                         const Eigen::MatrixXd& nominalCosine,
                                         const Eigen::MatrixXd& nominalSine,
                                         const double outputMu,
                                         const double outputRadius,
                                         const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter >& converter ):
        GravityFieldVariations( 2, 0, maximumDegree, maximumOrder ), sourceScale_( sourceScale ), convention_( convention ),
        converter_( converter )
    {
        // Reuse a date conversion across all selected coefficients and terms that share a printed token.
        std::map< std::string, std::pair< Time, Time > > dates;

        std::map< std::pair< int, int >, std::vector< const input_output::IcgemTemporalCoefficient* > > records;

        for( const auto& record : data.temporalCoefficients )
        {
            if( record.degree >= 2 && record.degree <= maximumDegree && record.order <= maximumOrder )
            {
                records[ { record.degree, record.order } ].push_back( &record );
            }
        }
        for( int degree = 2; degree <= maximumDegree; ++degree )
        {
            for( int order = 0; order <= std::min( degree, maximumOrder ); ++order )
            {
                const auto key = std::make_pair( degree, order );

                const auto temporal = records.find( key );

                const bool hasStatic = data.coefficients.count( key ) != 0;

                const bool hasTemporal = temporal != records.end( );

                if( ( hasStatic && hasTemporal ) || ( !hasStatic && !hasTemporal ) )
                {
                    throw std::runtime_error( "ICGEM selected coefficient requires exactly one static or temporal baseline" );
                }
                if( hasStatic )
                {
                    continue;
                }

                // Convert source normalization and potential constants once at the physical degree/order.
                const double scale = data.gravitationalParameter / outputMu * std::pow( data.referenceRadius / outputRadius, degree ) /
                        ( data.normalization == input_output::IcgemNormalization::unnormalized
                                  ? basic_mathematics::calculateLegendreGeodesyNormalizationFactor( degree, order )
                                  : 1.0 );

                TemporalCoefficient coefficient;

                coefficient.degree = degree;
                coefficient.order = order;
                for( const auto* record : temporal->second )
                {
                    if( record->type == "gfct" )
                    {
                        coefficient.bases.push_back( prepareTerm( *record, record->referenceDate, data.format, scale, overrides, dates ) );
                        coefficient.bases.back( ).cosine -= nominalCosine( degree, order );
                        coefficient.bases.back( ).sine -= nominalSine( degree, order );
                        if( !std::isfinite( coefficient.bases.back( ).cosine ) || !std::isfinite( coefficient.bases.back( ).sine ) )
                        {
                            throw std::runtime_error( "ICGEM adopted baseline subtraction is nonfinite" );
                        }
                    }
                }
                std::sort( coefficient.bases.begin( ), coefficient.bases.end( ), []( const TemporalTerm& left, const TemporalTerm& right ) {
                    return left.origin < right.origin;
                } );
                for( std::size_t index = 1; index < coefficient.bases.size( ); ++index )
                {
                    if( coefficient.bases[ index ].origin < coefficient.bases[ index - 1 ].end )
                    {
                        throw std::runtime_error( "ICGEM gfct intervals overlap after date interpretation" );
                    }
                }
                // Version 1.0 shares the sole gfct reference even when another record appeared first.
                const auto reference = std::find_if(
                        temporal->second.begin( ), temporal->second.end( ), []( const auto* record ) { return record->type == "gfct"; } );

                for( const auto* record : temporal->second )
                {
                    if( record->type != "gfct" )
                    {
                        coefficient.terms.push_back(
                                prepareTerm( *record,
                                             record->referenceDate.empty( ) ? ( *reference )->referenceDate : record->referenceDate,
                                             data.format,
                                             scale,
                                             overrides,
                                             dates ) );
                    }
                }
                validateTerms( coefficient );
                coefficients_.push_back( std::move( coefficient ) );
            }
        }
    }

    //! Evaluate half-open 2.0 intervals or unbounded 1.0 terms; an exact converted start selects the interval beginning there.
    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > calculateSphericalHarmonicsCorrections( const double time ) override
    {
        if( !std::isfinite( time ) )
        {
            throw std::runtime_error( "ICGEM temporal gravity requires finite TDB time" );
        }

        const auto boundary = boundaries_.find( time );

        const Time tai = boundary == boundaries_.end( ) ? converter_->getTaiTime( basic_astrodynamics::tdb_scale, time ) : boundary->second;

        const auto calendar = convention_ == IcgemYearConvention::fractionalGregorianYear
                ? calendarYear( converter_->getCurrentTime< Time >( basic_astrodynamics::tai_scale, sourceScale_, tai ) )
                : std::pair< int, long double >{ 0, 0.0L };

        Eigen::MatrixXd cosine = Eigen::MatrixXd::Zero( numberOfDegrees_, numberOfOrders_ );

        Eigen::MatrixXd sine = cosine;

        for( const auto& coefficient : coefficients_ )
        {
            // At a shared start and end, upper_bound selects the interval beginning at that epoch.
            auto base = std::upper_bound(
                    coefficient.bases.begin( ), coefficient.bases.end( ), tai, []( const Time& query, const TemporalTerm& term ) {
                        return query < term.origin;
                    } );

            if( !coefficient.bases.front( ).bounded )
            {
                base = coefficient.bases.begin( );
            }
            else
            {
                if( base == coefficient.bases.begin( ) || !( tai < ( --base )->end ) )
                {
                    throw std::runtime_error( "ICGEM temporal query is outside supplied gfct coverage or in a gap" );
                }
            }

            // Accumulate small anomalies directly; the full nominal coefficient was subtracted only at construction.
            long double c = base->cosine, s = base->sine;

            for( const auto& term : coefficient.terms )
            {
                if( term.bounded && ( tai < term.origin || !( tai < term.end ) ) )
                {
                    continue;
                }

                const long double years = convention_ == IcgemYearConvention::julianYear
                        ? static_cast< long double >( tai - term.origin ) / 31557600.0L
                        : static_cast< long double >( calendar.first - term.calendarOrigin.first ) + calendar.second -
                                term.calendarOrigin.second;

                const long double factor = term.type == "trnd" || term.type == "dot"
                        ? years
                        : ( term.type == "asin" ? std::sin( 2.0L * mathematical_constants::LONG_PI * years / term.period )
                                                : std::cos( 2.0L * mathematical_constants::LONG_PI * years / term.period ) );

                c += term.cosine * factor;
                s += term.sine * factor;
            }
            cosine( coefficient.degree - 2, coefficient.order ) = static_cast< double >( c );
            sine( coefficient.degree - 2, coefficient.order ) = static_cast< double >( s );
        }
        if( !cosine.allFinite( ) || !sine.allFinite( ) )
        {
            throw std::runtime_error( "ICGEM temporal correction is nonfinite" );
        }
        return { cosine, sine };
    }

private:
    //! Resolve dates and scale one record, retaining precise TAI epochs for exact converted interval starts and ends.
    TemporalTerm prepareTerm( const input_output::IcgemTemporalCoefficient& record,
                              const std::string& reference,
                              const std::string& format,
                              const double scale,
                              const std::map< std::string, double >& overrides,
                              std::map< std::string, std::pair< Time, Time > >& dates )
    {
        // Resolve each printed date once; the pair retains both source-calendar and TAI coordinates.
        const auto resolve = [ & ]( const std::string& token ) -> const std::pair< Time, Time >& {
            auto found = dates.find( token );

            if( found == dates.end( ) )
            {
                const Time source = sourceDate( token, overrides );

                found = dates.emplace( token, std::make_pair( source, converter_->getTaiTime( sourceScale_, source ) ) ).first;
            }
            return found->second;
        };

        const auto& origin = resolve( reference );

        TemporalTerm term;

        term.type = record.type;
        term.origin = origin.second;
        if( convention_ == IcgemYearConvention::fractionalGregorianYear )
        {
            term.calendarOrigin = calendarYear( origin.first );
        }
        term.bounded = format == "icgem2.0";
        term.period = record.periodYears.value_or( 0.0 );
        term.cosine = record.coefficient.cosine * scale;
        term.sine = record.coefficient.sine * scale;
        if( !std::isfinite( term.cosine ) || !std::isfinite( term.sine ) )
        {
            throw std::runtime_error( "ICGEM normalization or potential scaling is nonfinite" );
        }
        if( term.bounded )
        {
            term.end = resolve( record.validityEndDate ).second;
            if( !( term.origin < term.end ) )
            {
                throw std::runtime_error( "ICGEM interval is empty or reversed after date interpretation" );
            }
            for( const Time& edge : { term.origin, term.end } )
            {
                const double image = static_cast< double >(
                        converter_->getCurrentTime< Time >( basic_astrodynamics::tai_scale, basic_astrodynamics::tdb_scale, edge ) );

                // Distinct starts or ends that round to the same double TDB cannot be distinguished.
                const auto inserted = boundaries_.emplace( image, edge );

                if( !std::isfinite( image ) || ( !inserted.second && inserted.first->second != edge ) )
                {
                    throw std::runtime_error( "Distinct ICGEM interval boundaries have one scalar TDB image" );
                }
            }
        }
        return term;
    }

    //! Validate independently dated amplitudes after overrides, including coverage across adjacent base pieces.
    void validateTerms( const TemporalCoefficient& coefficient ) const
    {
        for( std::size_t index = 0; index < coefficient.terms.size( ); ++index )
        {
            const auto& term = coefficient.terms[ index ];

            for( std::size_t previous = 0; previous < index; ++previous )
            {
                const auto& other = coefficient.terms[ previous ];

                const bool sameKind = ( term.type == "dot" ? "trnd" : term.type ) == ( other.type == "dot" ? "trnd" : other.type );

                if( sameKind && term.period == other.period && ( !term.bounded || ( term.origin < other.end && other.origin < term.end ) ) )
                {
                    throw std::runtime_error( "ICGEM terms of the same kind and period overlap after date interpretation" );
                }
            }
            if( term.bounded )
            {
                Time covered = term.origin;

                for( const auto& base : coefficient.bases )
                {
                    if( covered < base.origin )
                    {
                        break;
                    }
                    if( covered < base.end )
                    {
                        covered = base.end;
                    }
                    if( !( covered < term.end ) )
                    {
                        break;
                    }
                }
                if( covered < term.end )
                {
                    throw std::runtime_error( "ICGEM temporal term extends beyond resolved gfct coverage" );
                }
            }
        }
    }

    //! Prepared selected coefficients and their declared source time/year interpretation.
    std::vector< TemporalCoefficient > coefficients_;

    basic_astrodynamics::TimeScales sourceScale_;

    IcgemYearConvention convention_;

    //! Use the same converter as the selected Earth rotation, with no separate default time model.
    std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > converter_;

    //! Retain TAI starts and ends for their exact converted double TDB values; convert other query times normally.
    std::map< double, Time > boundaries_;
};

}  // namespace

//! Validate the consumption contract and prepare additive temporal corrections against one nominal-field snapshot.
std::shared_ptr< GravityFieldVariations > createIcgemTemporalGravityFieldVariations(
        const input_output::IcgemGravityFieldData& data,
        const basic_astrodynamics::TimeScales sourceTimeScale,
        const IcgemYearConvention yearConvention,
        const int maximumDegree,
        const int maximumOrder,
        const std::map< std::string, double >& dateOverrides,
        const Eigen::MatrixXd& nominalCosine,
        const Eigen::MatrixXd& nominalSine,
        const double outputGravitationalParameter,
        const double outputReferenceRadius,
        const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter >& converter )
{
    input_output::validateIcgemGravityFieldData( data );
    if( !converter || sourceTimeScale < basic_astrodynamics::tai_scale || sourceTimeScale > basic_astrodynamics::ut1_scale ||
        ( yearConvention != IcgemYearConvention::julianYear && yearConvention != IcgemYearConvention::fractionalGregorianYear ) ||
        maximumDegree < 2 || maximumDegree > data.maximumDegree || maximumOrder < 0 || maximumOrder > maximumDegree ||
        nominalCosine.rows( ) <= maximumDegree || nominalSine.rows( ) <= maximumDegree || nominalCosine.cols( ) <= maximumOrder ||
        nominalSine.cols( ) <= maximumOrder || !nominalCosine.allFinite( ) || !nominalSine.allFinite( ) ||
        !std::isfinite( outputGravitationalParameter ) || outputGravitationalParameter <= 0.0 || !std::isfinite( outputReferenceRadius ) ||
        outputReferenceRadius <= 0.0 )
    {
        throw std::runtime_error(
                "ICGEM temporal gravity requires explicit supported time/year conventions, a selected converter and a valid output field "
                "block" );
    }
    for( const auto& entry : dateOverrides )
    {
        if( !std::isfinite( entry.second ) )
        {
            throw std::runtime_error( "ICGEM date overrides require finite source-scale epochs" );
        }
    }
    return std::make_shared< IcgemTemporalGravityFieldVariations >( data,
                                                                    sourceTimeScale,
                                                                    yearConvention,
                                                                    maximumDegree,
                                                                    maximumOrder,
                                                                    dateOverrides,
                                                                    nominalCosine,
                                                                    nominalSine,
                                                                    outputGravitationalParameter,
                                                                    outputReferenceRadius,
                                                                    converter );
}

}  // namespace gravitation
}  // namespace tudat
