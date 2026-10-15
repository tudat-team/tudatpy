/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include "tudat/astro/gravitation/icgemTemporalGravityFieldVariations.h"
#include "tudat/astro/basic_astro/dateTime.h"
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>

using namespace tudat;
using namespace tudat::basic_astrodynamics;
using namespace tudat::gravitation;
namespace
{
//! Locate the same pinned ICGEM extracts used by the neutral reader tests.
std::string fixture( const std::string& name )
{
    const std::string file = __FILE__;

    return file.substr( 0, file.find_last_of( "/\\" ) ) + "/../../../reference_files/icgem_gravity/" + name;
}

//! One independently computed scalar C/S correction and its arithmetic/time-sensitivity bounds.
struct Reference {
    //! Declared source-scale calendar seconds since J2000 noon.
    double sourceSeconds;

    //! Minus one requests the preceding TDB image; zero requests the exact image.
    int side;

    //! Whether this row uses fractional Gregorian years rather than elapsed TAI Julian years.
    bool fractionalYear;

    //! C/S anomalies, absolute arithmetic bounds and absolute coefficient rates per input second.
    Eigen::Vector2d value, rounding, absoluteRate;
};

//! Read the frozen scalar results without evaluating an ICGEM record in the test helper.
std::vector< Reference > references( const std::string& kind )
{
    std::ifstream file( fixture( "temporal_reference.txt" ) );

    BOOST_REQUIRE( file );
    std::vector< Reference > result;

    std::string line;

    std::getline( file, line );
    while( std::getline( file, line ) )
    {
        std::istringstream row( line );

        std::string label;

        Reference value;

        double printedCosine, printedSine;

        BOOST_REQUIRE( row >> label >> value.sourceSeconds >> value.side >> value.fractionalYear >> value.value.x( ) >> value.value.y( ) >>
                       value.rounding.x( ) >> value.rounding.y( ) >> value.absoluteRate.x( ) >> value.absoluteRate.y( ) >> printedCosine >>
                       printedSine );
        if( label == kind )
        {
            result.push_back( value );
        }
    }
    return result;
}

//! Map a declared source epoch to its public double TDB image using the supplied converter.
double tdbImage( const Time& source,
                 TimeScales scale,
                 const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter >& converter )
{
    return static_cast< double >( converter->getContinuousOutputTime( tdb_scale, converter->getTaiTime( scale, source ) ) );
}

//! Check a frozen scalar result; adjacent TDB images use a derivative bound, never an epoch-snapping tolerance.
Eigen::Vector2d checkReference( const std::shared_ptr< GravityFieldVariations >& model,
                                const Reference& expected,
                                TimeScales sourceScale,
                                int order,
                                const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter >& converter )
{
    double tdb = tdbImage( Time( expected.sourceSeconds ), sourceScale, converter );

    if( expected.side < 0 )
    {
        tdb = std::nextafter( tdb, -std::numeric_limits< double >::infinity( ) );
    }
    const auto actual = model->calculateSphericalHarmonicsCorrections( tdb );

    const Eigen::Vector2d value( actual.first( 0, order ), actual.second( 0, order ) );

    // Four input ULPs conservatively cover conversion rounding and an explicitly requested adjacent image.
    const double epochBudget = 4 * ( std::nextafter( tdb, std::numeric_limits< double >::infinity( ) ) - tdb );

    for( int i = 0; i < 2; ++i )
    {
        BOOST_CHECK_LE( std::abs( value[ i ] - expected.value[ i ] ), expected.rounding[ i ] + expected.absoluteRate[ i ] * epochBudget );
    }
    return value;
}
}  // namespace

//! Keep published long-period origins independent of adjacent gfct pieces and classify shared boundaries exactly.
BOOST_AUTO_TEST_CASE( eigenIndependentIntervalsAndOriginalDoubleBoundaryImages )
{
    const auto data = input_output::readIcgemGravityField( fixture( "eigen6s4_v2_c20_1950_2006.gfc" ) );

    const auto converter = std::make_shared< earth_orientation::TerrestrialTimeScaleConverter >( nullptr, nullptr );

    Eigen::MatrixXd nominal = Eigen::MatrixXd::Zero( 3, 1 );

    nominal( 2, 0 ) = -4.84165e-4;
    const Eigen::MatrixXd sine = Eigen::MatrixXd::Zero( 3, 1 );

    // This is an explicit caller choice for an irregular token, not an interpretation of its printed suffix.
    const std::map< std::string, double > overrides{ { "20041226.0060",
                                                       timeFromDecomposedDateTime< double >( 2004, 12, 26, 1, 0, 0.0L ) } };

    BOOST_CHECK_THROW( createIcgemTemporalGravityFieldVariations( data,
                                                                  tai_scale,
                                                                  IcgemYearConvention::julianYear,
                                                                  2,
                                                                  0,
                                                                  {},
                                                                  nominal,
                                                                  sine,
                                                                  data.gravitationalParameter,
                                                                  data.referenceRadius,
                                                                  converter ),
                       std::runtime_error );
    const auto model = createIcgemTemporalGravityFieldVariations( data,
                                                                  tai_scale,
                                                                  IcgemYearConvention::julianYear,
                                                                  2,
                                                                  0,
                                                                  overrides,
                                                                  nominal,
                                                                  sine,
                                                                  data.gravitationalParameter,
                                                                  data.referenceRadius,
                                                                  converter );

    const auto expected = references( "EIGEN" );

    BOOST_REQUIRE_EQUAL( expected.size( ), 7 );
    for( const auto& row : expected )
    {
        checkReference( model, row, tai_scale, 0, converter );
    }
    const double first = tdbImage( timeFromDecomposedDateTime< Time >( 1950, 1, 1, 0, 0, 0.0L ), tai_scale, converter );

    BOOST_CHECK_THROW(
            model->calculateSphericalHarmonicsCorrections( std::nextafter( first, -std::numeric_limits< double >::infinity( ) ) ),
            std::runtime_error );
    BOOST_CHECK_THROW( model->calculateSphericalHarmonicsCorrections(
                               tdbImage( timeFromDecomposedDateTime< Time >( 2006, 1, 1, 0, 0, 0.0L ), tai_scale, converter ) ),
                       std::runtime_error );
    BOOST_CHECK_THROW( createIcgemTemporalGravityFieldVariations( data,
                                                                  tai_scale,
                                                                  IcgemYearConvention::julianYear,
                                                                  2,
                                                                  1,
                                                                  overrides,
                                                                  Eigen::MatrixXd::Zero( 3, 2 ),
                                                                  Eigen::MatrixXd::Zero( 3, 2 ),
                                                                  data.gravitationalParameter,
                                                                  data.referenceRadius,
                                                                  converter ),
                       std::runtime_error );
}

//! Combine trend and annual/semiannual terms with normalization, GM/radius scaling and one nominal subtraction.
BOOST_AUTO_TEST_CASE( normalizedScaledAnomaliesAndDistinctLeapYearConventions )
{
    const auto data = input_output::readIcgemGravityField( fixture( "controlled_temporal_c21.gfc" ) );

    const auto converter = std::make_shared< earth_orientation::TerrestrialTimeScaleConverter >( nullptr, nullptr );

    Eigen::MatrixXd cosine = Eigen::MatrixXd::Zero( 3, 3 ), sine = cosine;

    cosine( 0, 0 ) = 1.0;
    cosine( 2, 0 ) = cosine( 2, 2 ) = 0.001;

    cosine( 2, 1 ) = 1e-6;
    sine( 2, 1 ) = -2e-6;
    const auto expected = references( "CONTROLLED" );

    BOOST_REQUIRE_EQUAL( expected.size( ), 4 );
    std::map< double, Eigen::Vector2d > julian;

    for( bool fractional : { false, true } )
    {
        const auto model = createIcgemTemporalGravityFieldVariations(
                data,
                tt_scale,
                fractional ? IcgemYearConvention::fractionalGregorianYear : IcgemYearConvention::julianYear,
                2,
                2,
                {},
                cosine,
                sine,
                5e14,
                7e6,
                converter );

        for( const auto& row : expected )
        {
            if( row.fractionalYear != fractional )
            {
                continue;
            }
            const auto value = checkReference( model, row, tt_scale, 1, converter );

            if( !fractional )
            {
                julian[ row.sourceSeconds ] = value;
            }
            else
            {
                BOOST_CHECK_GT( ( value - julian.at( row.sourceSeconds ) ).norm( ), 1e-12 );
            }
            const auto actual = model->calculateSphericalHarmonicsCorrections( tdbImage( Time( row.sourceSeconds ), tt_scale, converter ) );

            BOOST_CHECK_EQUAL( actual.first( 0, 0 ), 0.0 );
            BOOST_CHECK_EQUAL( actual.first( 0, 2 ), 0.0 );
        }
    }
}

//! Format 1.0 date-only origins are midnight and are inherited by its undated trend and periodic terms.
BOOST_AUTO_TEST_CASE( dateOnlyOriginsAreInheritedWithoutInventingValidityLimits )
{
    std::istringstream source(
            "product_type gravity_field\nmodelname controlled_date_only\nformat icgem1.0\n"
            "earth_gravity_constant 4e14\nradius 8e6\nmax_degree 2\nnorm fully_normalized\nerrors no\n"
            "tide_system tide_free\nend_of_head\ngfct 2 0 2e-6 0 20000101\ntrnd 2 0 4e-8 0\n"
            "acos 2 0 5e-8 0 1.0\nasin 2 0 -6e-8 0 1.0\n" );

    const auto data = input_output::readIcgemGravityField( source );

    const auto converter = std::make_shared< earth_orientation::TerrestrialTimeScaleConverter >( nullptr, nullptr );

    const auto model = createIcgemTemporalGravityFieldVariations( data,
                                                                  tai_scale,
                                                                  IcgemYearConvention::julianYear,
                                                                  2,
                                                                  0,
                                                                  {},
                                                                  Eigen::MatrixXd::Zero( 3, 1 ),
                                                                  Eigen::MatrixXd::Zero( 3, 1 ),
                                                                  4e14,
                                                                  8e6,
                                                                  converter );

    const Time sourceEpoch = timeFromDecomposedDateTime< Time >( 2000, 1, 1, 0, 0, 0.0L ) + 31557600.0 / 4;

    // At one quarter Julian year, the annual cosine vanishes and the sine is one:2e-6+1e-8-6e-8.
    BOOST_CHECK_SMALL(
            model->calculateSphericalHarmonicsCorrections( tdbImage( sourceEpoch, tai_scale, converter ) ).first( 0, 0 ) - 1.95e-6,
            32 * std::numeric_limits< double >::epsilon( ) * 2.1e-6 );
}
