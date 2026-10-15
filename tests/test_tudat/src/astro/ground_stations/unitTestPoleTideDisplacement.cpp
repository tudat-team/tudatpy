#define BOOST_TEST_MAIN

#include <boost/test/included/unit_test.hpp>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "tudat/astro/basic_astro/oblateSpheroidBodyShapeModel.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/ground_stations/poleTideDeformation.h"
#include "tudat/math/interpolators/linearInterpolator.h"
#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"

namespace tudat
{
namespace unit_tests
{

using namespace basic_astrodynamics;
using namespace earth_orientation;
using namespace interpolators;
using namespace coordinate_conversions;
using namespace ground_stations;

namespace
{
Time tdbFromTt( double tt )
{
    return getGeocentricTimeConverter( )->getCurrentTime< Time >( tt_scale, tdb_scale, Time( tt ) );
}
Eigen::Vector3d enuDisplacement( PoleTideDeformation& model, double tt, const std::shared_ptr< GroundStationState >& station )
{
    return station->getRotationFromBodyFixedToTopocentricFrame( tt ) * model.calculateDisplacement( tdbFromTt( tt ), station );
}

// Inject measured pole coordinates without external EOP files or short-period terms.
std::shared_ptr< PolarMotionCalculator > createConstantPolarMotionCalculator( const double tt, const double xpMas, const double ypMas )
{
    const Eigen::Vector2d polarMotion = unit_conversions::convertArcSecondsToRadians( Eigen::Vector2d( xpMas / 1000.0, ypMas / 1000.0 ) );

    // Daily records with a constant pole, on the UTC days around the test epoch.
    const double utcDay = std::floor( sofa_interface::convertTTtoUTC( tt ) / 86400.0 + 51544.5 );
    input_output::EarthOrientationData data;
    for( int i = -1; i <= 2; ++i )
    {
        input_output::EarthOrientationRecord record;
        record.utcModifiedJulianDate = utcDay + i;
        record.polarMotion = { polarMotion.x( ), polarMotion.y( ) };
        record.ut1MinusUtc = 0.0;
        data.records.emplace( record.utcModifiedJulianDate, record );
    }
    const auto zeroShortPeriod = std::make_shared< ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > >(
            1.0, 0.0, std::vector< std::string >( ), std::vector< std::string >( ), []( const double ) -> Eigen::Vector6d {
                return Eigen::Vector6d::Zero( );
            } );
    return std::make_shared< PolarMotionCalculator >( std::make_shared< EarthOrientationDataEvaluator >( data ), zeroShortPeriod );
}

void checkDisplacement( const Eigen::Vector3d& actual, const Eigen::Vector3d& expected )
{
    for( int component = 0; component < 3; ++component )
    {
        BOOST_CHECK_SMALL( actual( component ) - expected( component ), 1.0E-14 );
    }
}

struct PoleTideReferenceCase {
    double yearsSinceJ2000;
    double latitudeDegrees;
    double longitudeDegrees;
    double xpMas;
    double ypMas;
    double enuMetres[ 3 ];
    double cartesianMetres[ 3 ];
};

}  // namespace

BOOST_AUTO_TEST_SUITE( test_pole_tide_site_displacement )

BOOST_AUTO_TEST_CASE( testPoleTideIndependentIersReferenceCases )
{
    // Frozen fixtures evaluated independently from the accepted IERS chapter 7
    // Exact PDF hash and offline scalar reproduction: reference_files/solid_pole_tide/README.md.
    // (1 February 2018), https://iers-conventions.obspm.fr/content/chapter7/icc7.pdf.
    // Eq. (26) was evaluated in its original colatitude/south/east/radial form;
    // Eq. (28) was transposed to produce Cartesian values. Millimetres were
    // converted to metres, and south was negated for spherical ENU.
    // Wobble (m1,m2) in mas by row: (100,-200), (-120,80), (100,200),
    // (100,200), (0,100), (100,0), (-230,-140), (100,200).
    // Input xp/yp include the 2018 secular pole at each specified epoch.
    const PoleTideReferenceCase cases[] = {
        { 0.0,
          45.0,
          30.0,
          155.0,
          520.5,
          { 0.0014204684357863766, 0.0, 0.00044211616751135173 },
          { -0.00043949446353377188, 0.0013864734207246891, 0.00031262334011948438 } },
        { 10.0,
          -35.0,
          -70.0,
          -48.23,
          275.1,
          { -0.00044085860789372865, 0.00035773954020290583, -0.0036038981320963993 },
          { -0.0013537836600300739, 0.0024305056991498031, 0.0023601541232602585 } },
        { 0.0, 0.0, 0.0, 155.0, 120.5, { 0.0, -0.0009, 0.0 }, { 0.0, 0.0, -0.0009 } },
        { 0.0, 0.0, 90.0, 155.0, 120.5, { 0.0, -0.0018, 0.0 }, { 0.0, 0.0, -0.0018 } },
        { 0.0, 45.0, 90.0, 55.0, 220.5, { 0.0, 0.0, -0.0033 }, { 0.0, -0.002333452377915607, -0.0023334523779156074 } },
        { 18.0, -45.0, 180.0, 185.186, 382.78, { 0.0, 0.0, -0.0033 }, { 0.0023334523779156074, 0.0, 0.002333452377915607 } },
        { 26.0,
          67.0,
          145.0,
          -131.398,
          550.46,
          { -0.0020430005044474982, 0.000675859818697258, -0.0025662022104913657 },
          { 0.0025027975769723997, 0.00074156537064727157, -0.0020981221200375476 } },
        { 0.0, 90.0, 40.0, 155.0, 120.5, { -0.00080037114889627511, 0.0018464576962428508, 0.0 }, { -0.0009, -0.0018, 0.0 } }
    };

    for( const auto& reference : cases )
    {
        BOOST_TEST_CONTEXT( "latitude = " << reference.latitudeDegrees << ", longitude = " << reference.longitudeDegrees
                                          << ", years since J2000 = " << reference.yearsSinceJ2000 )
        {
            const double time = reference.yearsSinceJ2000 * physical_constants::JULIAN_YEAR;
            PoleTideDeformation model( createConstantPolarMotionCalculator( time, reference.xpMas, reference.ypMas ) );
            const Eigen::Vector3d sphericalPosition( 6378137.0,
                                                     unit_conversions::convertDegreesToRadians( reference.latitudeDegrees ),
                                                     unit_conversions::convertDegreesToRadians( reference.longitudeDegrees ) );
            const auto station = std::make_shared< GroundStationState >( sphericalPosition, spherical_position );
            const Eigen::Vector3d expectedEnu( reference.enuMetres[ 0 ], reference.enuMetres[ 1 ], reference.enuMetres[ 2 ] );
            const Eigen::Vector3d expectedCartesian(
                    reference.cartesianMetres[ 0 ], reference.cartesianMetres[ 1 ], reference.cartesianMetres[ 2 ] );
            checkDisplacement( enuDisplacement( model, time, station ), expectedEnu );
            checkDisplacement( model.calculateDisplacement( tdbFromTt( time ), station ), expectedCartesian );
        }
    }
}

BOOST_AUTO_TEST_CASE( testPoleTideOblateStationRepresentation )
{
    const double time = 26.0 * physical_constants::JULIAN_YEAR;
    PoleTideDeformation model( createConstantPolarMotionCalculator( time, 198.602, 610.46 ) );
    const auto shape = std::make_shared< OblateSpheroidBodyShapeModel >( 6378137.0, 1.0 / 298.257223563 );

    // WGS84 point: geodetic latitude 48 degrees, east longitude 123 degrees,
    // altitude 250 m. Cartesian coordinates independently evaluated from the
    // ellipsoid equations; IERS Eq. (26)/(28) uses atan2(z,hypot(x,y)).
    const Eigen::Vector3d cartesianPosition( -2328813.9273189958, 3586058.9739219621, 4717062.1163217751 );
    const Eigen::Vector3d geodeticPosition(
            250.0, unit_conversions::convertDegreesToRadians( 48.0 ), unit_conversions::convertDegreesToRadians( 123.0 ) );
    const auto sphericalStation = std::make_shared< GroundStationState >( cartesianPosition, cartesian_position );
    const std::vector< std::shared_ptr< GroundStationState > > oblateStations = {
        std::make_shared< GroundStationState >( cartesianPosition, cartesian_position, shape ),
        std::make_shared< GroundStationState >( geodeticPosition, geodetic_position, shape ),
        std::make_shared< GroundStationState >( sphericalStation->getNominalSphericalPosition( ), spherical_position, shape )
    };

    const Eigen::Vector3d expectedCartesian( -0.0026080881675806561, 0.0043229283367490524, 0.0052751641594906765 );
    const Eigen::Vector3d expectedSphericalEnu( -0.00016710873300951107, -0.00019573999494750939, 0.0072973249023081761 );

    // Independent projection of the same reference Cartesian vector onto the
    // geodetic east/north/normal basis at 48/123 degrees.
    const Eigen::Vector3d expectedGeodeticEnu( -0.00016710873300951118, -0.00022011967763310627, 0.0072966301950275498 );
    checkDisplacement( model.calculateDisplacement( tdbFromTt( time ), sphericalStation ), expectedCartesian );
    checkDisplacement( enuDisplacement( model, time, sphericalStation ), expectedSphericalEnu );
    for( const auto& station : oblateStations )
    {
        checkDisplacement( model.calculateDisplacement( tdbFromTt( time ), station ), expectedCartesian );
        const Eigen::Vector3d enu = enuDisplacement( model, time, station );
        checkDisplacement( enu, expectedGeodeticEnu );
        checkDisplacement( station->getRotationFromBodyFixedToTopocentricFrame( time ).inverse( ) * enu, expectedCartesian );
    }
}

BOOST_AUTO_TEST_CASE( testPoleTideZeroWobbleAnd2018SecularPole )
{
    // Secular pole values from IERS Eq. (21), frozen at J2000, 2018 and 2026.
    const double secularPoleCases[][ 3 ] = { { 0.0, 55.0, 320.5 }, { 18.0, 85.186, 382.78 }, { 26.0, 98.602, 410.46 } };
    const auto shape = std::make_shared< OblateSpheroidBodyShapeModel >( 6378137.0, 1.0 / 298.257223563 );
    const std::vector< std::shared_ptr< GroundStationState > > stations = {
        std::make_shared< GroundStationState >( Eigen::Vector3d( 6378137.0, 0.0, 0.0 ), spherical_position ),
        std::make_shared< GroundStationState >( Eigen::Vector3d( 6378137.0, -0.7, 1.1 ), spherical_position ),
        std::make_shared< GroundStationState >( Eigen::Vector3d( 250.0, 0.8, 2.1 ), geodetic_position, shape )
    };
    for( const auto& reference : secularPoleCases )
    {
        const double time = reference[ 0 ] * physical_constants::JULIAN_YEAR;
        PoleTideDeformation model( createConstantPolarMotionCalculator( time, reference[ 1 ], reference[ 2 ] ) );
        for( const auto& station : stations )
        {
            checkDisplacement( enuDisplacement( model, time, station ), Eigen::Vector3d::Zero( ) );
            checkDisplacement( model.calculateDisplacement( tdbFromTt( time ), station ), Eigen::Vector3d::Zero( ) );
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
