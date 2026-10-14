#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include <complex>
#include <limits>
#include "tudat/astro/ground_stations/atmosphericTideDeformation.h"
#include "tudat/astro/basic_astro/oblateSpheroidBodyShapeModel.h"
#include "tudat/astro/basic_astro/sphericalBodyShapeModel.h"
#include "tudat/simulation/environment_setup/createBodies.h"
#include "tudat/simulation/environment_setup/createBodyDeformationModel.h"
#include "tudat/astro/ephemerides/itrsToGcrsRotationModel.h"

using namespace tudat;
using namespace basic_astrodynamics;
using namespace input_output;
using namespace earth_orientation;
namespace
{
void close( const Eigen::Vector3d& a, const Eigen::Vector3d& b, double tolerance = 1e-13 )
{
    for( int i = 0; i < 3; ++i ) BOOST_CHECK_SMALL( a[ i ] - b[ i ], tolerance );
}
std::shared_ptr< TerrestrialTimeScaleConverter > modelClock( double shift = 0, double maximumGap = 129600 )
{
    EarthOrientationData data;
    for( int i = 0; i < 5; ++i )
    {
        EarthOrientationRecord row;
        row.utcModifiedJulianDate = 57752 + i;
        const Time utc( ( row.utcModifiedJulianDate - 51544.5 ) * 86400.0 );
        const Time tai = sofa_interface::convertUTCtoTAI( utc );
        row.polarMotion = { 0, 0 };
        row.ut1MinusUtc = -36.4 + shift + static_cast< double >( tai - utc );
        row.celestialPoleOffsets = std::array< double, 2 >{ 0, 0 };
        data.records.emplace( row.utcModifiedJulianDate, row );
    }
    EarthOrientationEvaluationSettings settings;
    settings.maximumGapSeconds = maximumGap;
    return std::make_shared< TerrestrialTimeScaleConverter >(
            nullptr, nullptr, nullptr, std::make_shared< EarthOrientationDataEvaluator >( data, settings ) );
}
AtmosphericTideEvaluationSettings explicitGeometry( )
{
    AtmosphericTideEvaluationSettings s;
    s.lookupLatitudeConvention = AtmosphericTideLatitudeConvention::geocentric;
    s.displacementBasis = AtmosphericTideDisplacementBasis::geocentricRadialNorthEast;
    return s;
}
AtmosphericTideGrid grid( bool latitudeDependent = false )
{
    AtmosphericTideGrid g;
    g.longitudesDegrees = { 0, 90 };
    g.latitudesDegrees = { 60, 0 };
    for( int longitude = 0; longitude < 2; ++longitude )
        for( double latitude : g.latitudesDegrees )
        {
            Eigen::Matrix< double, 3, 4 > c = Eigen::Matrix< double, 3, 4 >::Zero( );
            c.col( 0 ) << 1e-3 * ( latitudeDependent ? 1 + latitude : 1 ), 2e-3, 3e-3;
            c.col( 1 ) << 4e-3, 5e-3, 6e-3;
            c.col( 2 ) << 7e-3, 8e-3, 9e-3;
            c.col( 3 ) << 10e-3, 11e-3, 12e-3;
            g.coefficientsMetres.push_back( c );
        }
    return g;
}
Time tdbAtUt1( const std::shared_ptr< TerrestrialTimeScaleConverter >& c, const Time& ut1 )
{
    return c->getCurrentTime< Time >( ut1_scale, tdb_scale, ut1 );
}
Eigen::Vector3d independentLocal( const Eigen::Matrix< double, 3, 4 >& c, double ut1 )
{
    // Complex Fourier synthesis is independent of the production real sine/cosine expression.
    const double phase = 2 * mathematical_constants::PI * ( ( ut1 - 86400.0 * std::floor( ut1 / 86400.0 ) ) / 86400.0 + .5 );
    Eigen::Vector3d result;
    for( int i = 0; i < 3; ++i )
        result[ i ] = std::real( std::complex< double >( c( i, 0 ), -c( i, 1 ) ) * std::exp( std::complex< double >( 0, phase ) ) +
                                 std::complex< double >( c( i, 2 ), -c( i, 3 ) ) * std::exp( std::complex< double >( 0, 2 * phase ) ) );
    return result;
}
Eigen::Vector3d independentLocalRate( const Eigen::Matrix< double, 3, 4 >& c, double ut1 )
{
    const double phase = 2 * mathematical_constants::PI * ( ( ut1 - 86400 * std::floor( ut1 / 86400 ) ) / 86400 + .5 );
    const std::complex< double > iOmega( 0, 2 * mathematical_constants::PI / 86400 );
    Eigen::Vector3d result;
    for( int i = 0; i < 3; ++i )
        result[ i ] =
                std::real( iOmega *
                           ( std::complex< double >( c( i, 0 ), -c( i, 1 ) ) * std::exp( std::complex< double >( 0, phase ) ) +
                             2. * std::complex< double >( c( i, 2 ), -c( i, 3 ) ) * std::exp( std::complex< double >( 0, 2 * phase ) ) ) );
    return result;
}
Eigen::Matrix3d basis( const Eigen::Vector3d& p )
{
    const Eigen::Vector3d r = p.normalized( );
    const Eigen::Vector3d e = Eigen::Vector3d::UnitZ( ).cross( r ).normalized( );
    Eigen::Matrix3d m;
    m.col( 0 ) = r;
    m.col( 1 ) = r.cross( e );
    m.col( 2 ) = e;
    return m;
}
std::shared_ptr< ground_stations::GroundStationState > station( const Eigen::Vector3d& p,
                                                                const std::shared_ptr< BodyShapeModel >& shape = nullptr,
                                                                const std::string& id = "7839" )
{
    auto s = std::make_shared< ground_stations::GroundStationState >( p, coordinate_conversions::cartesian_position, shape );
    s->setSiteId( id );
    return s;
}
AtmosphericTideMetadata rp03( )
{
    AtmosphericTideMetadata m;
    m.modelName = "Ray-Ponte 2003";
    m.greensFunction = "Farrell elastic Green functions";
    m.sourceProfile = AtmosphericTideSourceProfile::iers2010Rp03NoIbFarrell;
    m.referenceOrigin = AtmosphericLoadingReferenceOrigin::solidEarthCentre;
    return m;
}
}  // namespace
BOOST_AUTO_TEST_CASE( nativeContinuousLeapPhaseAndOwnClock )
{
    const auto c = modelClock( );
    const auto g = grid( );
    AtmosphericTideDeformation model( g, c, explicitGeometry( ) );
    const Eigen::Vector3d p( 6378137, 0, 0 );
    const Time midnightTai( ( 57754.0 - 51544.5 ) * 86400.0 + 37 );
    BOOST_CHECK_THROW( EarthOrientationDataEvaluator::convertTaiToScalarUtc( midnightTai - .5 ), std::runtime_error );
    for( double offset : { -2., -1., -.5, -.25, 0., .25, .5, 1., 2. } )
    {
        const Time tai = midnightTai + offset;
        const Time tdb = c->getCurrentTime< Time >( tai_scale, tdb_scale, tai );
        const double ut1 = static_cast< double >( tai - 36.4 );
        const auto local = independentLocal( g.coefficientsMetres.front( ), ut1 );
        const Eigen::Vector3d expected( local.x( ), local.z( ), local.y( ) );
        // The model owns the converter used by this native TDB query.
        close( model.calculateDisplacement( tdb, station( p ) ), expected, 2e-14 );
        close( model.calculateDisplacement( Time( static_cast< double >( tdb ) ), station( p ) ), expected, 3e-13 );
    }
    const auto early = tdbAtUt1( c, Time( ( 57752. - 51544.5 ) * 86400. - .4 ) );
    BOOST_CHECK_THROW( model.calculateDisplacement( early - 1.0, station( p ) ), std::runtime_error );
    AtmosphericTideDeformation gap( grid( ), modelClock( 0, 3600 ), explicitGeometry( ) );
    BOOST_CHECK_THROW( gap.calculateDisplacement( tdbAtUt1( c, Time( ( 57754. - 51544.5 ) * 86400. - .4 ) ), station( p ) ),
                       std::runtime_error );
}
BOOST_AUTO_TEST_CASE( independentLookupAndBasisAndSameCartesianRepresentation )
{
    const auto c = modelClock( );
    const double ut1 = ( 57754. - 51544.5 ) * 86400. + 12345;
    const Time epoch( tdbAtUt1( c, Time( ut1 ) ) );
    const auto oblate = std::make_shared< OblateSpheroidBodyShapeModel >( 6378137, 1. / 298.257223563 );
    const auto sphere = std::make_shared< SphericalBodyShapeModel >( 6378137 );
    // Known ellipsoid surface coordinate 45deg,30deg; normal direction is known independently.
    const double f = 1. / 298.257223563, e2 = f * ( 2 - f ), lat = mathematical_constants::PI / 4, lon = mathematical_constants::PI / 6;
    const double n = 6378137 / std::sqrt( 1 - e2 * std::pow( std::sin( lat ), 2 ) );
    const Eigen::Vector3d p(
            n * std::cos( lat ) * std::cos( lon ), n * std::cos( lat ) * std::sin( lon ), n * ( 1 - e2 ) * std::sin( lat ) );
    const Eigen::Vector3d normal( std::cos( lat ) * std::cos( lon ), std::cos( lat ) * std::sin( lon ), std::sin( lat ) );
    Eigen::Matrix3d normalBasis;
    normalBasis.col( 0 ) = normal;
    normalBasis.col( 2 ) = Eigen::Vector3d::UnitZ( ).cross( normal ).normalized( );
    normalBasis.col( 1 ) = normal.cross( normalBasis.col( 2 ) );
    const auto g = grid( true );
    const double gc = std::asin( p.normalized( ).z( ) ) * 180 / mathematical_constants::PI;
    for( bool geodeticLookup : { false, true } )
        for( bool geodeticBasis : { false, true } )
        {
            auto settings = explicitGeometry( );
            if( geodeticLookup ) settings.lookupLatitudeConvention = AtmosphericTideLatitudeConvention::geodetic;
            if( geodeticBasis ) settings.displacementBasis = AtmosphericTideDisplacementBasis::geodeticUpNorthEast;
            AtmosphericTideDeformation model( g, c, settings, oblate );
            auto coeff = g.coefficientsMetres.front( );
            coeff( 0, 0 ) = 1e-3 * ( 1 + ( geodeticLookup ? 45 : gc ) );
            const Eigen::Vector3d expected = ( geodeticBasis ? normalBasis : basis( p ) ) * independentLocal( coeff, ut1 );
            close( model.calculateDisplacement( epoch, station( p, oblate ) ), expected, 1e-13 );
            close( model.calculateDisplacement( epoch, station( p, sphere ) ), expected, 1e-13 );
            close( model.calculateDisplacement( epoch, station( p ) ), expected, 1e-13 );
            const Eigen::Vector3d expectedRate = ( geodeticBasis ? normalBasis : basis( p ) ) *
                    independentLocalRate( coeff, static_cast< double >( c->getCurrentTime< Time >( tdb_scale, ut1_scale, epoch ) ) ) *
                    c->getGeocentricTtRateWithRespectToTdb( epoch );
            close( model.calculateDisplacementState( epoch, station( p, oblate ) ).tail< 3 >( ), expectedRate, 2e-20 );
            close( model.calculateDisplacementState( epoch, station( p, sphere ) ).tail< 3 >( ), expectedRate, 2e-20 );
            close( model.calculateDisplacementState( epoch, station( p ) ).tail< 3 >( ), expectedRate, 2e-20 );
            const Eigen::Vector3d expectedEnuRate( expectedRate.dot( normalBasis.col( 2 ) ),
                                                   expectedRate.dot( normalBasis.col( 1 ) ),
                                                   expectedRate.dot( normalBasis.col( 0 ) ) );
            const Eigen::Vector3d actualEnuRate =
                    station( p, oblate )->getRotationFromBodyFixedToTopocentricFrame( static_cast< double >( epoch ) ) *
                    model.calculateDisplacementState( epoch, station( p ) ).tail< 3 >( );
            close( actualEnuRate, expectedEnuRate, 2e-20 );
            const auto s = station( p, oblate );
            close( s->getRotationFromBodyFixedToTopocentricFrame( static_cast< double >( epoch ) ) *
                           model.calculateDisplacement( epoch, s ),
                   s->getRotationFromBodyFixedToTopocentricFrame( static_cast< double >( epoch ) ) * expected );
        }
}
BOOST_AUTO_TEST_CASE( requiredGeometryAndSourceConflicts )
{
    const auto c = modelClock( );
    auto g = grid( );
    auto settings = explicitGeometry( );
    BOOST_CHECK_THROW( static_cast< void >( AtmosphericTideDeformation( g, c ) ), std::runtime_error );
    BOOST_CHECK_THROW( static_cast< void >( AtmosphericTideDeformation( g, nullptr, settings ) ), std::runtime_error );
    settings.lookupLatitudeConvention = AtmosphericTideLatitudeConvention::geodetic;
    BOOST_CHECK_THROW( static_cast< void >( AtmosphericTideDeformation( g, c, settings ) ), std::runtime_error );
    settings = explicitGeometry( );
    settings.displacementBasis = AtmosphericTideDisplacementBasis::geodeticUpNorthEast;
    BOOST_CHECK_THROW( static_cast< void >( AtmosphericTideDeformation( g, c, settings ) ), std::runtime_error );
    settings = explicitGeometry( );
    g.metadata.latitudeConvention = AtmosphericTideLatitudeConvention::geodetic;
    BOOST_CHECK_THROW( static_cast< void >( AtmosphericTideDeformation( g, c, settings ) ), std::runtime_error );
    g = grid( );
    g.metadata.displacementBasis = AtmosphericTideDisplacementBasis::geodeticUpNorthEast;
    BOOST_CHECK_THROW( static_cast< void >( AtmosphericTideDeformation( g, c, settings ) ), std::runtime_error );
}
BOOST_AUTO_TEST_CASE( optInRp03OriginExactlyOnceInCartesian )
{
    const auto c = modelClock( );
    auto g = grid( );
    g.metadata = rp03( );
    auto settings = explicitGeometry( );
    AtmosphericTideDeformation ce( g, c, settings );
    settings.originPolicy = AtmosphericLoadingOriginPolicy::convertIdentifiedRp03CeToCm;
    AtmosphericTideDeformation cm( g, c, settings );
    const double ut1 = ( 57754. - 51544.5 ) * 86400. + 12345;
    const Time epoch( tdbAtUt1( c, Time( ut1 ) ) );
    const Eigen::Vector3d p( 4000000, 3000000, 3500000 );
    const auto translation = calculateIers2010Rp03AtmosphericGeocentreTranslation( ut1 );
    close( cm.calculateDisplacement( epoch, station( p ) ) - ce.calculateDisplacement( epoch, station( p ) ), translation );
    BOOST_CHECK( cm.getOutputReferenceOrigin( ) == AtmosphericLoadingReferenceOrigin::wholeEarthCentreOfMass );
    BOOST_CHECK( ce.getOutputReferenceOrigin( ) == AtmosphericLoadingReferenceOrigin::solidEarthCentre );
    g.metadata.referenceOrigin = AtmosphericLoadingReferenceOrigin::wholeEarthCentreOfMass;
    AtmosphericTideDeformation alreadyCm( g, c, settings );
    close( alreadyCm.calculateDisplacement( epoch, station( p ) ), ce.calculateDisplacement( epoch, station( p ) ) );
    for( int wrong = 0; wrong < 4; ++wrong )
    {
        auto bad = g;
        if( wrong == 0 ) bad.metadata.sourceProfile = AtmosphericTideSourceProfile::unspecified;
        if( wrong == 1 ) bad.metadata.modelName = "Vienna APL";
        if( wrong == 2 ) bad.metadata.greensFunction = "other";
        if( wrong == 3 ) bad.metadata.referenceOrigin = AtmosphericLoadingReferenceOrigin::unspecified;
        BOOST_CHECK_THROW( static_cast< void >( AtmosphericTideDeformation( bad, c, settings ) ), std::runtime_error );
    }
}
BOOST_AUTO_TEST_CASE( exactStationCoefficientsAndAnalyticRate )
{
    const auto c = modelClock( );
    AtmosphericTideStationData data;
    AtmosphericTideStation row;
    row.stationIdentifier = "7839";
    row.coefficientsMetres = grid( ).coefficientsMetres.front( );
    data.stations.emplace( row.stationIdentifier, row );
    BOOST_CHECK_THROW( static_cast< void >( AtmosphericTideDeformation( data, c ) ), std::runtime_error );
    AtmosphericTideDeformation model( data, c, explicitGeometry( ) );
    const Eigen::Vector3d p( 6378137, 0, 0 );
    const auto s = station( p );
    const auto epoch = tdbAtUt1( c, Time( ( 57754. - 51544.5 ) * 86400. ) );
    close( model.calculateDisplacement( epoch, s ),
           basis( p ) *
                   independentLocal( row.coefficientsMetres,
                                     static_cast< double >( c->getCurrentTime< Time >( tdb_scale, ut1_scale, epoch ) ) ) );
    BOOST_CHECK_THROW( model.calculateDisplacement( epoch, station( p, nullptr, "783" ) ), std::out_of_range );
    const std::shared_ptr< BodyDeformationModel > base = std::make_shared< AtmosphericTideDeformation >( data, c, explicitGeometry( ) );
    const double ut1 = static_cast< double >( c->getCurrentTime< Time >( tdb_scale, ut1_scale, epoch ) );
    const Eigen::Vector3d expected =
            basis( p ) * independentLocalRate( row.coefficientsMetres, ut1 ) * c->getGeocentricTtRateWithRespectToTdb( epoch );
    close( base->calculateDisplacementState( epoch, s ).tail< 3 >( ), expected, 2e-20 );
    close( base->calculateDisplacementState( Time( static_cast< double >( epoch ) ), s ).tail< 3 >( ), expected, 2e-18 );
    close( base->calculateDisplacementState( epoch, s ).head< 3 >( ), model.calculateDisplacement( epoch, s ) );
}
BOOST_AUTO_TEST_CASE( modelEndpointPoliciesAndRegionalCoverage )
{
    auto g = grid( );
    g.longitudesDegrees = { 0, 180, 360 };
    g.coefficientsMetres = { g.coefficientsMetres[ 0 ], g.coefficientsMetres[ 1 ], g.coefficientsMetres[ 0 ],
                             g.coefficientsMetres[ 1 ], g.coefficientsMetres[ 0 ], g.coefficientsMetres[ 1 ] };
    g.coefficientsMetres.back( )( 0, 0 ) += 8e-7;
    const auto c = modelClock( );
    auto settings = explicitGeometry( );
    BOOST_CHECK_THROW( static_cast< void >( AtmosphericTideDeformation( g, c, settings ) ), std::runtime_error );
    settings.endpointPolicy = AtmosphericTideEndpointPolicy::averageSeamAndPolarCartesianApproximation;
    AtmosphericTideDeformation approximate( g, c, settings );
    BOOST_CHECK_SMALL( approximate.getGridInterpolator( )->getMaximumSeamDiscrepancyMetres( ) - 8e-7, 1e-18 );
    BOOST_CHECK_EQUAL( approximate.getGridInterpolator( )->getProviderGrid( ).coefficientsMetres.back( )( 0, 0 ),
                       g.coefficientsMetres.back( )( 0, 0 ) );
    const Time epoch( tdbAtUt1( c, Time( ( 57754. - 51544.5 ) * 86400. ) ) );
    BOOST_CHECK_THROW( approximate.calculateDisplacement( epoch, station( Eigen::Vector3d( 0, 0, 6378137 ) ) ), std::out_of_range );
    AtmosphericTideDeformation regional( grid( ), c, explicitGeometry( ) );
    BOOST_CHECK_THROW( regional.calculateDisplacement( epoch, station( Eigen::Vector3d( 0, -6378137, 0 ) ) ), std::out_of_range );
}
BOOST_AUTO_TEST_CASE( nativeFactoryPreloadedGridStationFileParityAndRejections )
{
    using namespace simulation_setup;
    const auto c = modelClock( );
    const auto angles = std::make_shared< EarthOrientationAnglesCalculator >( nullptr, nullptr, c );
    SystemOfBodies bodies;
    bodies.createEmptyBody( "Earth" );
    bodies.at( "Earth" )->setRotationalEphemeris( std::make_shared< ephemerides::GcrsToItrsRotationModel >( angles ) );
    bodies.at( "Earth" )->setShapeModel( std::make_shared< OblateSpheroidBodyShapeModel >( 6378137, 1. / 298.257223563 ) );
    const auto base = std::string( __FILE__ ).substr( 0, std::string( __FILE__ ).find_last_of( "/\\" ) ) +
            "/../../../reference_files/atmospheric_tides/";
    const auto path = base + "rp03_ce_grid_extract.txt";
    const auto g = readAtmosphericTideGrid( path );
    const auto settings = explicitGeometry( );
    const auto preloaded = createBodyDeformationModel( atmosphericTideBodyShapeDeformationFromGrid( g, settings ), "Earth", bodies );
    // File settings already own the provider grid before an environment model is constructed.
    const auto fileSettings = atmosphericTideBodyShapeDeformationFromGridFile( path, {}, settings );

    BOOST_CHECK_EQUAL( fileSettings->getGrid( ).metadata.source, g.metadata.source );
    BOOST_CHECK_EQUAL( fileSettings->getGrid( ).coefficientsMetres.size( ), g.coefficientsMetres.size( ) );
    BOOST_CHECK_THROW( atmosphericTideBodyShapeDeformationFromGridFile( path + ".missing" ), std::runtime_error );

    const auto file = createBodyDeformationModel( fileSettings, "Earth", bodies );

    const double lat = -19.5 * mathematical_constants::PI / 180, lon = 47.5 * mathematical_constants::PI / 180;
    const Eigen::Vector3d p =
            6378137 * Eigen::Vector3d( std::cos( lat ) * std::cos( lon ), std::cos( lat ) * std::sin( lon ), std::sin( lat ) );
    const Time midnightTai( ( 57754. - 51544.5 ) * 86400. + 37 );
    for( double offset : { -.5, 0., .5 } )
    {
        const auto tdb = c->getCurrentTime< Time >( tai_scale, tdb_scale, midnightTai + offset );
        close( preloaded->calculateDisplacement( tdb, station( p ) ), file->calculateDisplacement( tdb, station( p ) ) );
        const Eigen::Vector6d a = preloaded->calculateDisplacementState( tdb, station( p ) );
        const Eigen::Vector6d b = file->calculateDisplacementState( tdb, station( p ) );
        close( a.head< 3 >( ), b.head< 3 >( ) );
        close( a.tail< 3 >( ), b.tail< 3 >( ), 2e-20 );
    }
    const auto stations = readAtmosphericTideStationData( base + "bernese_manual_stations.atl" );
    const auto stationPreloaded =
            createBodyDeformationModel( atmosphericTideBodyShapeDeformationFromStations( stations, settings ), "Earth", bodies );
    const auto stationFile = createBodyDeformationModel(
            atmosphericTideBodyShapeDeformationFromStationFile( base + "bernese_manual_stations.atl", {}, settings ), "Earth", bodies );
    const auto namedStation = station( p, nullptr, stations.stations.begin( )->first );
    const Time epoch( c->getCurrentTime< Time >( tai_scale, tdb_scale, midnightTai ) );
    close( stationPreloaded->calculateDisplacement( epoch, namedStation ), stationFile->calculateDisplacement( epoch, namedStation ) );
    close( stationPreloaded->calculateDisplacementState( epoch, namedStation ).tail< 3 >( ),
           stationFile->calculateDisplacementState( epoch, namedStation ).tail< 3 >( ),
           2e-20 );
    BOOST_CHECK_THROW( createBodyDeformationModel( atmosphericTideBodyShapeDeformationFromGrid( g ), "Earth", bodies ),
                       std::runtime_error );
    BOOST_CHECK_THROW( createBodyDeformationModel( atmosphericTideBodyShapeDeformationFromGrid( g, settings ), "Mars", bodies ),
                       std::runtime_error );
    BOOST_CHECK_THROW(
            createBodyDeformationModel( std::make_shared< BodyDeformationSettings >( gravitation::atmospheric_tide ), "Earth", bodies ),
            std::runtime_error );
    bodies.at( "Earth" )->setRotationalEphemeris( nullptr );
    BOOST_CHECK_THROW( createBodyDeformationModel( atmosphericTideBodyShapeDeformationFromGrid( g, settings ), "Earth", bodies ),
                       std::runtime_error );
}
namespace
{
std::shared_ptr< TerrestrialTimeScaleConverter > slopedClock( const std::vector< double >& slopes,
                                                              EarthOrientationEvaluationSettings settings = {} )
{
    EarthOrientationData data;
    Time previous;
    double baseline = -36.4;
    for( std::size_t i = 0; i < slopes.size( ) + 1; ++i )
    {
        EarthOrientationRecord row;
        row.utcModifiedJulianDate = 57752 + i;
        const Time utc( ( row.utcModifiedJulianDate - 51544.5 ) * 86400 );
        const Time tai = sofa_interface::convertUTCtoTAI( utc );
        if( i ) baseline += slopes[ i - 1 ] * static_cast< double >( tai - previous );
        row.polarMotion = { 0, 0 };
        row.ut1MinusUtc = baseline + static_cast< double >( tai - utc );
        row.celestialPoleOffsets = std::array< double, 2 >{ 0, 0 };
        data.records.emplace( row.utcModifiedJulianDate, row );
        previous = tai;
    }
    return std::make_shared< TerrestrialTimeScaleConverter >(
            nullptr, nullptr, nullptr, std::make_shared< EarthOrientationDataEvaluator >( data, settings ) );
}
Time epochAtTai( const std::shared_ptr< TerrestrialTimeScaleConverter >& c, const Time& tai )
{
    return c->getCurrentTime< Time >( tai_scale, tdb_scale, tai );
}
}  // namespace
BOOST_AUTO_TEST_CASE( selectedDailySlopesRightKnotsFinalLeftAndScalarBoundary )
{
    const std::vector< double > slopes{ .01, .02, .03, .04 };
    const auto c = slopedClock( slopes );
    AtmosphericTideDeformation model( grid( ), c, explicitGeometry( ) );
    const Eigen::Vector3d p( 6378137, 0, 0 );
    const auto s = station( p );
    const auto& epochs = c->getEarthOrientationDataEvaluator( )->getEpochsTai( );
    for( std::size_t i = 0; i < epochs.size( ); ++i )
        for( double shift : { 0., 1000. } )
        {
            if( i + 1 == epochs.size( ) && shift > 0 ) continue;
            const Time tai = epochs[ i ] + shift;
            const auto epoch = epochAtTai( c, tai );
            if( shift == 0 )
            {
                BOOST_REQUIRE_MESSAGE( c->getTaiTime( tdb_scale, epoch ) == tai,
                                       "TAI round trip at node " << i << " differs by "
                                                                 << static_cast< long double >( c->getTaiTime( tdb_scale, epoch ) - tai )
                                                                 << " s" );
            }
            else
            {
                // Time stores seconds within an hour as long double. Interior
                // clock round trips need only respect that representation's
                // arithmetic error; exact knot identity above controls side choice.
                BOOST_REQUIRE_SMALL( static_cast< long double >( c->getTaiTime( tdb_scale, epoch ) - tai ),
                                     16.0L * std::numeric_limits< long double >::epsilon( ) * TIME_NORMALIZATION_INTEGER_TERM );
            }
            const auto value = c->getEarthOrientationDataEvaluator( )->evaluate( tai );
            const double slope = slopes[ std::min( i, slopes.size( ) - 1 ) ];
            const double ut1 = static_cast< double >( tai + value.ut1MinusTai );
            const Eigen::Vector3d local = independentLocalRate( grid( ).coefficientsMetres.front( ), ut1 );
            const Eigen::Vector3d expected =
                    Eigen::Vector3d( local.x( ), local.z( ), local.y( ) ) * ( 1 + slope ) * c->getGeocentricTtRateWithRespectToTdb( epoch );
            close( model.calculateDisplacementState( epoch, s ).tail< 3 >( ), expected, 2e-20 );
            // Only the ground-station scalar boundary normalizes rounded endpoint images.
            std::vector< std::shared_ptr< BodyDeformationModel > > models{ std::make_shared< AtmosphericTideDeformation >(
                    grid( ), c, explicitGeometry( ) ) };
            ground_stations::BodyDeformationStationMotionModel adapter( [ & ]( ) -> auto& { return models; }, true );
            const auto scalarState = adapter.getBodyFixedStationMotion( static_cast< double >( epoch ), s );
            if( shift == 0 && ( i == 0 || i + 1 == epochs.size( ) ) )
                close( scalarState.tail< 3 >( ), expected, 2e-20 );
            else
                close( scalarState.tail< 3 >( ),
                       model.calculateDisplacementState( Time( static_cast< double >( epoch ) ), s ).tail< 3 >( ),
                       2e-20 );
        }
}
BOOST_AUTO_TEST_CASE( unavailableSelectedRatePreservesNodePositionAndRejectsState )
{
    EarthOrientationEvaluationSettings policy;
    policy.maximumGapSeconds = 3600;
    const auto c = slopedClock( { .01, .02, .03, .04 }, policy );
    AtmosphericTideDeformation model( grid( ), c, explicitGeometry( ) );
    const Eigen::Vector3d p( 6378137, 0, 0 );
    const auto& nodes = c->getEarthOrientationDataEvaluator( )->getEpochsTai( );
    const auto epoch = epochAtTai( c, nodes[ 2 ] );
    BOOST_CHECK( model.calculateDisplacement( epoch, station( p ) ).allFinite( ) );
    BOOST_CHECK_THROW( model.calculateDisplacementState( epoch, station( p ) ), std::runtime_error );
    BOOST_CHECK_THROW( model.calculateDisplacement( epochAtTai( c, nodes[ 2 ] + 1000. ), station( p ) ), std::runtime_error );
    BOOST_CHECK_THROW( model.calculateDisplacement( epochAtTai( c, nodes.front( ) - 1. ), station( p ) ), std::runtime_error );
    BOOST_CHECK_THROW( model.calculateDisplacementState( epochAtTai( c, nodes.back( ) + 1. ), station( p ) ), std::runtime_error );
}
BOOST_AUTO_TEST_CASE( analyticRatesSmoothConvergenceLeapContinuityAndOriginOnce )
{
    const auto c = slopedClock( { .02, .02, .02, .02 } );
    auto g = grid( );
    g.metadata = rp03( );
    auto options = explicitGeometry( );
    AtmosphericTideDeformation ce( g, c, options );
    options.originPolicy = AtmosphericLoadingOriginPolicy::convertIdentifiedRp03CeToCm;
    AtmosphericTideDeformation cm( g, c, options );
    const Eigen::Vector3d p( 4000000, 3000000, 3500000 );
    const Time midnightTai( ( 57754. - 51544.5 ) * 86400. + 37 );
    for( double offset : { -.75, -.5, -.25, 0., .25, .5, .75, 12345. } )
    {
        const auto epoch = epochAtTai( c, midnightTai + offset );
        const double rate = 1.02 * c->getGeocentricTtRateWithRespectToTdb( epoch );
        const Eigen::Vector3d expected = differentiateIers2010Rp03AtmosphericGeocentreTranslationWrtUt1(
                                                 static_cast< double >( c->getCurrentTime< Time >( tdb_scale, ut1_scale, epoch ) ) ) *
                rate;
        close( cm.calculateDisplacementState( epoch, station( p ) ).tail< 3 >( ) -
                       ce.calculateDisplacementState( epoch, station( p ) ).tail< 3 >( ),
               expected,
               2e-20 );
        const Eigen::Vector3d analytic = cm.calculateDisplacementState( epoch, station( p ) ).tail< 3 >( );
        double prior = 1;
        for( double h : { 32., 16., 8. } )
        {
            const Eigen::Vector3d numeric =
                    ( cm.calculateDisplacement( epoch + h, station( p ) ) - cm.calculateDisplacement( epoch - h, station( p ) ) ) /
                    ( 2 * h );
            const double error = ( numeric - analytic ).norm( );
            BOOST_CHECK_LT( error, prior * .3 );
            prior = error;
        }
        BOOST_CHECK_SMALL( prior, 1e-12 );
    }
}
