/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <limits>
#include "tudat/astro/gravitation/aod1bGravityFieldVariations.h"

namespace
{

//! Load the small provider extract used by the raw-reader tests.
tudat::input_output::Aod1bGravityFieldData fixture( )
{
    const std::string source = __FILE__;

    return tudat::input_output::readAod1bGravityField( source.substr( 0, source.find_last_of( "/\\" ) ) +
                                                       "/../../../reference_files/aod1b_gravity/AOD1B_2026-01-01_degree2_extract_07.asc" );
}

//! Use the public double TDB image of a GPS epoch, with the explicitly selected converter.
double tdbFromGps( const double gps, const std::shared_ptr< tudat::earth_orientation::TerrestrialTimeScaleConverter >& converter )
{
    return static_cast< double >( converter->getCurrentTime< tudat::Time >(
            tudat::basic_astrodynamics::tai_scale, tudat::basic_astrodynamics::tdb_scale, tudat::Time( gps ) + 19.0 ) );
}

}  // namespace

//! Test actual converted endpoints, GPS midpoint weights and one GM/radius rescaling of the source potential.
BOOST_AUTO_TEST_CASE( gpsInterpolationAndOutputScaling )
{
    const auto data = fixture( );

    const auto converter = std::make_shared< tudat::earth_orientation::TerrestrialTimeScaleConverter >( nullptr, nullptr );

    const auto model = tudat::gravitation::createAod1bGravityFieldVariations( { data },
                                                                              tudat::input_output::Aod1bComponent::global,
                                                                              2,
                                                                              2,
                                                                              data.gravitationalParameter / 2.0,
                                                                              data.referenceRadius / 2.0,
                                                                              converter );

    const double start = data.samples.at( 2 ).gpsSecondsSinceJ2000;

    const double c0 = data.samples.at( 2 ).coefficients.at( { 2, 0 } ).values[ 0 ];

    const double c1 = data.samples.at( 6 ).coefficients.at( { 2, 0 } ).values[ 0 ];

    BOOST_CHECK_EQUAL( model->getMinimumDegree( ), 2 );
    BOOST_CHECK_EQUAL( model->getMaximumDegree( ), 2 );
    BOOST_CHECK_EQUAL( model->calculateSphericalHarmonicsCorrections( tdbFromGps( start, converter ) ).first( 0, 0 ), 8.0 * c0 );
    BOOST_CHECK_EQUAL( model->calculateSphericalHarmonicsCorrections( tdbFromGps( start + 10800.0, converter ) ).first( 0, 0 ), 8.0 * c1 );
    BOOST_CHECK_NO_THROW( model->calculateSphericalHarmonicsCorrections( tdbFromGps( start + 21600.0, converter ) ) );

    const double midpointTdb = tdbFromGps( start + 5400.0, converter );

    // One public-double epoch ULP bounds time rounding; the remaining term bounds arithmetic rounding.
    const double budget = 8.0 * std::abs( c1 - c0 ) / 10800.0 *
                    ( std::nextafter( midpointTdb, std::numeric_limits< double >::infinity( ) ) - midpointTdb ) +
            16.0 * std::numeric_limits< double >::epsilon( ) * std::abs( c0 + c1 );

    BOOST_CHECK_SMALL( model->calculateSphericalHarmonicsCorrections( midpointTdb ).first( 0, 0 ) - 4.0 * ( c0 + c1 ), budget );
    BOOST_CHECK_THROW( model->calculateSphericalHarmonicsCorrections( tdbFromGps( start - 1.0, converter ) ), std::runtime_error );
    BOOST_CHECK_THROW( model->calculateSphericalHarmonicsCorrections( tdbFromGps( start + 21601.0, converter ) ), std::runtime_error );
}

//! Merge supplied files without interpolating across missing epochs or mixing product metadata.
BOOST_AUTO_TEST_CASE( sourceCoverageAcrossFilesAndGaps )
{
    auto first = fixture( );

    auto second = first;

    first.samples.resize( 8 );
    second.samples.erase( second.samples.begin( ), second.samples.begin( ) + 8 );

    const auto converter = std::make_shared< tudat::earth_orientation::TerrestrialTimeScaleConverter >( nullptr, nullptr );

    const double start = first.samples.front( ).gpsSecondsSinceJ2000;

    const auto merged = tudat::gravitation::createAod1bGravityFieldVariations( { first, second },
                                                                               tudat::input_output::Aod1bComponent::global,
                                                                               2,
                                                                               0,
                                                                               first.gravitationalParameter,
                                                                               first.referenceRadius,
                                                                               converter );

    BOOST_CHECK_NO_THROW( merged->calculateSphericalHarmonicsCorrections( tdbFromGps( start + 16200.0, converter ) ) );
    first.samples.resize( 4 );

    const auto gap = tudat::gravitation::createAod1bGravityFieldVariations( { first, second },
                                                                            tudat::input_output::Aod1bComponent::global,
                                                                            2,
                                                                            0,
                                                                            first.gravitationalParameter,
                                                                            first.referenceRadius,
                                                                            converter );

    BOOST_CHECK_THROW( gap->calculateSphericalHarmonicsCorrections( tdbFromGps( start + 10800.0, converter ) ), std::runtime_error );
    BOOST_CHECK_NO_THROW( gap->calculateSphericalHarmonicsCorrections( tdbFromGps( start + 21600.0, converter ) ) );
    second.softwareVersion = "atm_ocean_dealise.07_v2";
    BOOST_CHECK_THROW( tudat::gravitation::createAod1bGravityFieldVariations( { first, second },
                                                                              tudat::input_output::Aod1bComponent::global,
                                                                              2,
                                                                              0,
                                                                              first.gravitationalParameter,
                                                                              first.referenceRadius,
                                                                              converter ),
                       std::runtime_error );
}
