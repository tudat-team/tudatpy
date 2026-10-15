#define BOOST_TEST_MAIN

#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <limits>

#include "tudat/io/basicInputOutput.h"
#include "tudat/astro/ground_stations/oceanTideEarthDeformation.h"
#include "tudat/astro/basic_astro/dateTime.h"
#include "tudat/astro/basic_astro/oblateSpheroidBodyShapeModel.h"
#include "../../../reference_files/ocean_loading/hardisp_reference_cases.h"
#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/environment_setup/createBodyDeformationModel.h"

namespace tudat
{
namespace unit_tests
{

using namespace basic_astrodynamics;
using namespace earth_orientation;

Time tdbFromTt( double tt )
{
    return getGeocentricTimeConverter( )->getCurrentTime< Time >( tt_scale, tdb_scale, Time( tt ) );
}
Eigen::Vector3d enuDisplacement( OceanTideEarthDeformation& model, const Time& tdb, const std::string& id )
{
    auto station = std::make_shared< ground_stations::GroundStationState >( Eigen::Vector3d( 6378137, 0, 0 ) );
    station->setSiteId( id );
    return station->getRotationFromBodyFixedToTopocentricFrame( static_cast< double >( tdb ) ) *
            model.calculateDisplacement( tdb, station );
}

BOOST_AUTO_TEST_SUITE( test_ocean_tide_site_displacement )

BOOST_AUTO_TEST_CASE( testOceanTideSiteDisplacement )
{
    std::string testFile = paths::getEarthDeformationDataFilesPath( ) + "/oceanTidetest.blq";

    OceanTideEarthDeformation oceanTideDeformationModel = OceanTideEarthDeformation( testFile );

    double timeStep = 3600.0;
    double testTime = timeFromDecomposedDateTime< double >( 2009, 06, 25, 01, 10, 45 );

    Eigen::Matrix< double, 24, 3 > onsalaResults;
    Eigen::Matrix< double, 24, 3 > icelandResults;

    for( int i = 0; i < 24; i++ )
    {
        onsalaResults.block( i, 0, 1, 3 ) =
                enuDisplacement( oceanTideDeformationModel, tdbFromTt( sofa_interface::convertUTCtoTT( testTime ) ), "ONSALA" )
                        .transpose( );
        icelandResults.block( i, 0, 1, 3 ) =
                enuDisplacement( oceanTideDeformationModel, tdbFromTt( sofa_interface::convertUTCtoTT( testTime ) ), "REYKJAVIK" )
                        .transpose( );
        testTime += timeStep;
    }

    Eigen::Matrix< double, 24, 3 > onsalaRawTestMatrix;
    onsalaRawTestMatrix << 0.003094, -0.001538, -0.000895, 0.001812, -0.000950, -0.000193, 0.000218, -0.000248, 0.000421, -0.001104,
            0.000404, 0.000741, -0.001668, 0.000863, 0.000646, -0.001209, 0.001042, 0.000137, 0.000235, 0.000926, -0.000667, 0.002337,
            0.000580, -0.001555, 0.004554, 0.000125, -0.002278, 0.006271, -0.000291, -0.002615, 0.006955, -0.000537, -0.002430, 0.006299,
            -0.000526, -0.001706, 0.004305, -0.000244, -0.000559, 0.001294, 0.000245, 0.000793, -0.002163, 0.000819, 0.002075, -0.005375,
            0.001326, 0.003024, -0.007695, 0.001622, 0.003448, -0.008669, 0.001610, 0.003272, -0.008143, 0.001262, 0.002557, -0.006290,
            0.000633, 0.001477, -0.003566, -0.000155, 0.000282, -0.000593, -0.000941, -0.000766, 0.001992, -0.001561, -0.001457, 0.003689,
            -0.001889, -0.001680;

    Eigen::Matrix< double, 24, 3 > icelandRawTestMatrix;
    icelandRawTestMatrix << -0.005940, -0.001245, -0.000278, 0.013516, -0.001086, 0.003212, 0.029599, -0.000353, 0.005483, 0.038468,
            0.000699, 0.005997, 0.038098, 0.001721, 0.004690, 0.028780, 0.002363, 0.001974, 0.013016, 0.002371, -0.001369, -0.005124,
            0.001653, -0.004390, -0.021047, 0.000310, -0.006225, -0.030799, -0.001383, -0.006313, -0.032056, -0.003048, -0.004549,
            -0.024698, -0.004288, -0.001314, -0.010814, -0.004794, 0.002623, 0.005849, -0.004416, 0.006291, 0.020857, -0.003208, 0.008766,
            0.030226, -0.001413, 0.009402, 0.031437, 0.000594, 0.007996, 0.024079, 0.002389, 0.004844, 0.009945, 0.003606, 0.000663,
            -0.007426, 0.004022, -0.003581, -0.023652, 0.003601, -0.006911, -0.034618, 0.002505, -0.008585, -0.037515, 0.001044, -0.008270,
            -0.031544, -0.000402, -0.006125;

    Eigen::Matrix< double, 24, 3 > onsalaTestMatrix = onsalaRawTestMatrix;
    onsalaTestMatrix.block( 0, 0, 24, 1 ) = -onsalaRawTestMatrix.block( 0, 2, 24, 1 );
    onsalaTestMatrix.block( 0, 1, 24, 1 ) = -onsalaRawTestMatrix.block( 0, 1, 24, 1 );
    onsalaTestMatrix.block( 0, 2, 24, 1 ) = onsalaRawTestMatrix.block( 0, 0, 24, 1 );

    Eigen::Matrix< double, 24, 3 > icelandTestMatrix = icelandRawTestMatrix;
    //
    icelandTestMatrix.block( 0, 0, 24, 1 ) = -icelandRawTestMatrix.block( 0, 2, 24, 1 );
    icelandTestMatrix.block( 0, 1, 24, 1 ) = -icelandRawTestMatrix.block( 0, 1, 24, 1 );
    icelandTestMatrix.block( 0, 2, 24, 1 ) = icelandRawTestMatrix.block( 0, 0, 24, 1 );

    BOOST_CHECK_SMALL( ( onsalaResults - onsalaTestMatrix ).block( 0, 0, 24, 1 ).cwiseAbs( ).maxCoeff( ), 1.0E-6 );
    BOOST_CHECK_SMALL( ( onsalaResults - onsalaTestMatrix ).block( 0, 1, 24, 1 ).cwiseAbs( ).maxCoeff( ), 1.0E-6 );
    BOOST_CHECK_SMALL( ( onsalaResults - onsalaTestMatrix ).block( 0, 2, 24, 1 ).cwiseAbs( ).maxCoeff( ), 1.0E-6 );

    BOOST_CHECK_SMALL( ( icelandResults - icelandTestMatrix ).block( 0, 0, 24, 1 ).cwiseAbs( ).maxCoeff( ), 1.0E-6 );
    BOOST_CHECK_SMALL( ( icelandResults - icelandTestMatrix ).block( 0, 1, 24, 1 ).cwiseAbs( ).maxCoeff( ), 1.0E-6 );
    BOOST_CHECK_SMALL( ( icelandResults - icelandTestMatrix ).block( 0, 2, 24, 1 ).cwiseAbs( ).maxCoeff( ), 1.0E-6 );
}

BOOST_AUTO_TEST_CASE( testOceanLoadingCompatibilityAndPreloadedData )
{
    const std::string file = paths::getEarthDeformationDataFilesPath( ) + "/oceanTidetest.blq";
    const auto data = input_output::readOceanLoadingData( file );
    const auto compatibility = readBlqFile( file );
    const auto& onsala = data.getStation( "ONSALA" );
    BOOST_CHECK_EQUAL( compatibility[ 0 ].at( "ONSALA" )( 0, 0 ), 0.00352 );
    BOOST_CHECK_SMALL( compatibility[ 1 ].at( "ONSALA" )( 0, 0 ) + 1.1292280260403312, 1.0E-14 );
    BOOST_CHECK_EQUAL( onsala.phaseLagsDegrees( 0, 0 ), -64.7 );
    const std::function< Eigen::Vector6d( const double ) > arguments = []( const double time ) {
        return sofa_interface::calculateDoodsonFundamentalArguments( time );
    };
    OceanTideEarthDeformation fromFile( file, arguments );
    OceanTideEarthDeformation preloaded( data, arguments );

    // Both empty construction routes support later append while preserving provider station provenance.
    OceanTideEarthDeformation initiallyEmpty( input_output::OceanLoadingData{}, arguments );

    initiallyEmpty.addOceanLoadingData( data );
    for( const double time : { 0.0, 1.0E8, 3.0E8 } )
    {
        BOOST_CHECK_SMALL( ( enuDisplacement( fromFile, time, "ONSALA" ) - enuDisplacement( preloaded, time, "ONSALA" ) ).norm( ),
                           1.0E-14 );
        BOOST_CHECK_SMALL( ( enuDisplacement( initiallyEmpty, time, "ONSALA" ) - enuDisplacement( preloaded, time, "ONSALA" ) ).norm( ),
                           1.0E-14 );
    }
    BOOST_CHECK_EQUAL( preloaded.getOceanLoadingData( ).getStation( "ONSALA" ).source, file );
    BOOST_CHECK_EQUAL( preloaded.getOceanLoadingData( ).getStation( "ONSALA" ).oceanTideModel, "CSR4.0_f_PP" );
    auto conflicting = data;
    auto additional = onsala;
    additional.stationIdentifier = "ADDED_FIRST";
    conflicting.stations.emplace( additional.stationIdentifier, additional );
    const auto stationCount = preloaded.getOceanLoadingData( ).stations.size( );
    const auto headerCount = preloaded.getOceanLoadingData( ).sourceHeaders.size( );
    const auto before = enuDisplacement( preloaded, 1.0E8, "ONSALA" );
    BOOST_CHECK_THROW( preloaded.addOceanLoadingData( conflicting ), std::runtime_error );
    BOOST_CHECK_EQUAL( preloaded.getOceanLoadingData( ).stations.size( ), stationCount );
    BOOST_CHECK_EQUAL( preloaded.getOceanLoadingData( ).sourceHeaders.size( ), headerCount );
    BOOST_CHECK_EQUAL( preloaded.getOceanLoadingData( ).stations.count( "ADDED_FIRST" ), 0 );
    BOOST_CHECK_SMALL( ( enuDisplacement( preloaded, 1.0E8, "ONSALA" ) - before ).norm( ), 1.0E-14 );
    BOOST_CHECK_THROW( enuDisplacement( preloaded, 0.0, "MISSING" ), std::out_of_range );
}

namespace
{

input_output::OceanLoadingData makeHardispReferenceData( )
{
    input_output::OceanLoadingData data;
    data.source = "Pinned IERS HARDISP/ILRS reference inputs";
    for( const auto& reference : hardispReferenceStations )
    {
        input_output::OceanLoadingStation station;
        station.stationIdentifier = reference.identifier;
        station.source = data.source;
        for( int row = 0; row < 3; ++row )
        {
            for( int column = 0; column < 11; ++column )
            {
                station.amplitudesMetres( row, column ) = reference.coefficients[ row ][ column ];
                station.phaseLagsDegrees( row, column ) = reference.coefficients[ row + 3 ][ column ];
            }
        }
        data.stations.emplace( station.stationIdentifier, station );
    }
    return data;
}

void checkOceanLoadingVector( const Eigen::Vector3d& actual, const Eigen::Vector3d& expected, const double tolerance )
{
    for( int component = 0; component < 3; ++component )
    {
        BOOST_CHECK_SMALL( actual( component ) - expected( component ), tolerance );
    }
}

}  // namespace

BOOST_AUTO_TEST_CASE( testExpandedOceanLoadingAgainstIndependentHardisp )
{
    OceanTideEarthDeformation model( makeHardispReferenceData( ) );
    BOOST_CHECK( model.getSynthesis( ) == OceanTideSynthesis::iersHardisp );
    for( const auto& reference : hardispReferenceCases )
    {
        const auto& station = hardispReferenceStations[ reference.stationIndex ];
        BOOST_TEST_CONTEXT( station.identifier << ", UTC epoch " << reference.year << "-" << reference.month << "-" << reference.day )
        {
            const double utc = timeFromDecomposedDateTime< double >(
                    reference.year, reference.month, reference.day, reference.hour, reference.minute, reference.second );
            const double tt = sofa_interface::convertUTCtoTT( utc );

            // HARDISP writes U/S/W; the public spherical result is E/N/radial.
            const Eigen::Vector3d expected( -reference.upSouthWest[ 2 ], -reference.upSouthWest[ 1 ], reference.upSouthWest[ 0 ] );
            checkOceanLoadingVector( enuDisplacement( model, tdbFromTt( tt ), station.identifier ), expected, 1.0E-6 );
        }
    }
}

BOOST_AUTO_TEST_CASE( testOceanLoadingSphericalAndOblateGeometry )
{
    const auto& reference = hardispReferenceCases[ 35 ];
    const auto& coefficients = hardispReferenceStations[ reference.stationIndex ];
    const double utc = timeFromDecomposedDateTime< double >(
            reference.year, reference.month, reference.day, reference.hour, reference.minute, reference.second );
    const double tt = sofa_interface::convertUTCtoTT( utc );
    OceanTideEarthDeformation model( makeHardispReferenceData( ) );

    // A fixed WGS84 point, independently evaluated from geodetic 48/123 degrees, height 250 m.
    const Eigen::Vector3d point( -2328813.9273189958, 3586058.9739219621, 4717062.1163217751 );
    const auto shape = std::make_shared< OblateSpheroidBodyShapeModel >( 6378137.0, 1.0 / 298.257223563 );
    const auto spherical = std::make_shared< ground_stations::GroundStationState >( point, coordinate_conversions::cartesian_position );
    const auto oblate = std::make_shared< ground_stations::GroundStationState >( point, coordinate_conversions::cartesian_position, shape );
    spherical->setSiteId( coefficients.identifier );
    oblate->setSiteId( coefficients.identifier );
    const double phi = std::atan2( point( 2 ), std::hypot( point( 0 ), point( 1 ) ) );
    const double lambda = 123.0 * mathematical_constants::PI / 180.0;
    const double geodeticLatitude = 48.0 * mathematical_constants::PI / 180.0;
    const Eigen::Vector3d east( -std::sin( lambda ), std::cos( lambda ), 0.0 );
    const Eigen::Vector3d north( -std::sin( phi ) * std::cos( lambda ), -std::sin( phi ) * std::sin( lambda ), std::cos( phi ) );
    const Eigen::Vector3d radial( std::cos( phi ) * std::cos( lambda ), std::cos( phi ) * std::sin( lambda ), std::sin( phi ) );
    const Eigen::Vector3d sphericalEnu( -reference.upSouthWest[ 2 ], -reference.upSouthWest[ 1 ], reference.upSouthWest[ 0 ] );
    const Eigen::Vector3d expectedCartesian = sphericalEnu( 0 ) * east + sphericalEnu( 1 ) * north + sphericalEnu( 2 ) * radial;
    checkOceanLoadingVector( model.calculateDisplacement( tdbFromTt( tt ), spherical ), expectedCartesian, 1.0E-6 );
    checkOceanLoadingVector( model.calculateDisplacement( tdbFromTt( tt ), oblate ), expectedCartesian, 1.0E-6 );
    checkOceanLoadingVector(
            model.calculateDisplacement( tdbFromTt( tt ), spherical ), model.calculateDisplacement( tdbFromTt( tt ), oblate ), 1.0E-14 );
    const Eigen::Vector3d geodeticNorth( -std::sin( geodeticLatitude ) * std::cos( lambda ),
                                         -std::sin( geodeticLatitude ) * std::sin( lambda ),
                                         std::cos( geodeticLatitude ) );
    const Eigen::Vector3d geodeticUp( std::cos( geodeticLatitude ) * std::cos( lambda ),
                                      std::cos( geodeticLatitude ) * std::sin( lambda ),
                                      std::sin( geodeticLatitude ) );
    const Eigen::Vector3d expectedTopocentric(
            east.dot( expectedCartesian ), geodeticNorth.dot( expectedCartesian ), geodeticUp.dot( expectedCartesian ) );
    checkOceanLoadingVector(
            ( oblate->getRotationFromBodyFixedToTopocentricFrame( tt ) * model.calculateDisplacement( tdbFromTt( tt ), oblate ) ),
            expectedTopocentric,
            1.0E-6 );
    const Eigen::Vector3d uws( reference.upSouthWest[ 0 ], reference.upSouthWest[ 2 ], reference.upSouthWest[ 1 ] );
    checkOceanLoadingVector(
            convertLocalOceanTideDisplacementToPlentocentricDisplacement( uws, { east, north, radial } ), expectedCartesian, 1.0E-14 );
}

BOOST_AUTO_TEST_CASE( testOceanLoadingLegacyModeAndZeroAmplitudes )
{
    auto data = makeHardispReferenceData( );
    for( auto& station : data.stations )
    {
        station.second.amplitudesMetres.setZero( );
    }
    OceanTideEarthDeformation zero( data );
    checkOceanLoadingVector( enuDisplacement( zero, 1.0E8, "ONSALA" ), Eigen::Vector3d::Zero( ), 1.0E-15 );
    data.stations.at( "ONSALA" ).amplitudesMetres( 0, 0 ) = 0.01;
    data.stations.at( "ONSALA" ).phaseLagsDegrees( 0, 0 ) = 0.1 * 180.0 / mathematical_constants::PI;
    const std::function< Eigen::Vector6d( const double ) > arguments = []( const double ) {
        Eigen::Vector6d result = Eigen::Vector6d::Zero( );
        result( 0 ) = 0.25;
        return result;
    };
    OceanTideEarthDeformation legacy( data, arguments );
    BOOST_CHECK( legacy.getSynthesis( ) == OceanTideSynthesis::legacyElevenConstituents );

    // Independent single-M2 fixture: 0.01*cos(0.5-0.1), with all other amplitudes zero.
    checkOceanLoadingVector( enuDisplacement( legacy, 0.0, "ONSALA" ), Eigen::Vector3d( 0.0, 0.0, 0.009210609940028852 ), 1.0E-14 );
    OceanTideEarthDeformation explicitLegacy( data, OceanTideSynthesis::legacyElevenConstituents );
    OceanTideEarthDeformation callbackLegacy(
            data, []( const double time ) { return sofa_interface::calculateDoodsonFundamentalArguments( time ); } );
    checkOceanLoadingVector(
            enuDisplacement( explicitLegacy, tdbFromTt( 1.0E8 ), "ONSALA" ), enuDisplacement( callbackLegacy, 1.0E8, "ONSALA" ), 1.0E-14 );
    BOOST_CHECK_THROW( enuDisplacement( zero, 0.0, "MISSING" ), std::out_of_range );
    BOOST_CHECK_THROW( enuDisplacement( zero, std::numeric_limits< double >::quiet_NaN( ), "ONSALA" ), std::invalid_argument );
}

BOOST_AUTO_TEST_CASE( testOceanLoadingFactoriesConsumePreloadedData )
{
    const auto data = makeHardispReferenceData( );
    simulation_setup::SystemOfBodies bodies;
    bodies.createEmptyBody( "Earth" );

    // Empty filename settings normalize to empty provider data, allowing station records to be appended later.
    const auto emptySettings = simulation_setup::oceanTideBodyShapeDeformation( std::vector< std::string >{} );

    const auto empty = std::dynamic_pointer_cast< OceanTideEarthDeformation >(
            simulation_setup::createBodyDeformationModel( emptySettings, "Earth", bodies ) );

    BOOST_REQUIRE( empty != nullptr );
    BOOST_CHECK( empty->getOceanLoadingData( ).stations.empty( ) );
    empty->addOceanLoadingData( data );
    BOOST_CHECK_EQUAL( empty->getOceanLoadingData( ).stations.size( ), data.stations.size( ) );
    const auto fullSettings = simulation_setup::oceanTideBodyShapeDeformationFromData( data );
    const auto full = std::dynamic_pointer_cast< OceanTideEarthDeformation >(
            simulation_setup::createBodyDeformationModel( fullSettings, "Earth", bodies ) );
    BOOST_REQUIRE( full != nullptr );
    BOOST_CHECK( full->getSynthesis( ) == OceanTideSynthesis::iersHardisp );
    const auto legacySettings =
            simulation_setup::oceanTideBodyShapeDeformationFromData( data, OceanTideSynthesis::legacyElevenConstituents );
    const auto legacy = std::dynamic_pointer_cast< OceanTideEarthDeformation >(
            simulation_setup::createBodyDeformationModel( legacySettings, "Earth", bodies ) );
    BOOST_REQUIRE( legacy != nullptr );
    BOOST_CHECK( legacy->getSynthesis( ) == OceanTideSynthesis::legacyElevenConstituents );
    const auto& reference = hardispReferenceCases[ 0 ];
    const double utc = timeFromDecomposedDateTime< double >(
            reference.year, reference.month, reference.day, reference.hour, reference.minute, reference.second );
    const double tt = sofa_interface::convertUTCtoTT( utc );
    const Eigen::Vector3d expected( -reference.upSouthWest[ 2 ], -reference.upSouthWest[ 1 ], reference.upSouthWest[ 0 ] );
    checkOceanLoadingVector( enuDisplacement( *full, tdbFromTt( tt ), "ONSALA" ), expected, 1.0E-6 );
    BOOST_CHECK_EQUAL( full->getOceanLoadingData( ).stations.size( ), data.stations.size( ) );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
