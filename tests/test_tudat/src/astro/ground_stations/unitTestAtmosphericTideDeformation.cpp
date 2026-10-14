#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include <complex>
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
            AtmosphericTideCoefficients c = AtmosphericTideCoefficients::Zero( );
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
Eigen::Vector3d independentLocal( const AtmosphericTideCoefficients& c, double ut1 )
{
    // Complex Fourier synthesis is independent of the production real sine/cosine expression.
    const double phase = 2 * mathematical_constants::PI * ( ( ut1 - 86400.0 * std::floor( ut1 / 86400.0 ) ) / 86400.0 + .5 );
    Eigen::Vector3d result;
    for( int i = 0; i < 3; ++i )
        result[ i ] = std::real( std::complex< double >( c( i, 0 ), -c( i, 1 ) ) * std::exp( std::complex< double >( 0, phase ) ) +
                                 std::complex< double >( c( i, 2 ), -c( i, 3 ) ) * std::exp( std::complex< double >( 0, 2 * phase ) ) );
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
    const auto c = modelClock( ), other = modelClock( 3600 );
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
        // Even a context supplied with another clock must use this model's clock.
        close( model.calculateDisplacement( EarthDeformationEpoch( tdb, other ), p ), expected, 2e-14 );
        close( model.calculateDisplacement( static_cast< double >( tdb ), p ), expected, 3e-13 );
    }
    const auto early = tdbAtUt1( c, Time( ( 57752. - 51544.5 ) * 86400. - .4 ) );
    BOOST_CHECK_THROW( model.calculateDisplacement( EarthDeformationEpoch( early - 1.0, c ), p ), std::runtime_error );
    AtmosphericTideDeformation gap( grid( ), modelClock( 0, 3600 ), explicitGeometry( ) );
    BOOST_CHECK_THROW(
            gap.calculateDisplacement( EarthDeformationEpoch( tdbAtUt1( c, Time( ( 57754. - 51544.5 ) * 86400. - .4 ) ), c ), p ),
            std::runtime_error );
}
BOOST_AUTO_TEST_CASE( independentLookupAndBasisAndSameCartesianRepresentation )
{
    const auto c = modelClock( );
    const double ut1 = ( 57754. - 51544.5 ) * 86400. + 12345;
    const EarthDeformationEpoch epoch( tdbAtUt1( c, Time( ut1 ) ), c );
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
            close( model.calculateDisplacement( epoch, p ), expected, 1e-13 );
            const auto s = station( p, oblate );
            close( model.calculateEnuDisplacement( epoch, s ),
                   s->getRotationFromBodyFixedToTopocentricFrame( static_cast< double >( epoch.getTdb( ) ) ) * expected );
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
    const EarthDeformationEpoch epoch( tdbAtUt1( c, Time( ut1 ) ), c );
    const Eigen::Vector3d p( 4000000, 3000000, 3500000 );
    const auto translation = calculateIers2010Rp03AtmosphericGeocentreTranslation( ut1 );
    close( cm.calculateDisplacement( epoch, p ) - ce.calculateDisplacement( epoch, p ), translation );
    close( cm.calculateDisplacement( epoch, p ) - ce.calculateDisplacement( epoch, p ), translation );
    BOOST_CHECK( cm.getOutputReferenceOrigin( ) == AtmosphericLoadingReferenceOrigin::wholeEarthCentreOfMass );
    BOOST_CHECK( ce.getOutputReferenceOrigin( ) == AtmosphericLoadingReferenceOrigin::solidEarthCentre );
    g.metadata.referenceOrigin = AtmosphericLoadingReferenceOrigin::wholeEarthCentreOfMass;
    AtmosphericTideDeformation alreadyCm( g, c, settings );
    close( alreadyCm.calculateDisplacement( epoch, p ), ce.calculateDisplacement( epoch, p ) );
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
BOOST_AUTO_TEST_CASE( exactStationCoefficientsAndUnavailableRates )
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
    const auto epoch = EarthDeformationEpoch( tdbAtUt1( c, Time( ( 57754. - 51544.5 ) * 86400. ) ), c );
    close( model.calculateDisplacement( epoch, s ),
           basis( p ) * independentLocal( row.coefficientsMetres, static_cast< double >( epoch.getUt1( c ) ) ) );
    BOOST_CHECK_THROW( model.calculateDisplacement( epoch, station( p, nullptr, "783" ) ), std::runtime_error );
    BOOST_CHECK_THROW( model.calculateDisplacement( epoch, p ), std::runtime_error );
    const std::shared_ptr< BodyDeformationModel > base = std::make_shared< AtmosphericTideDeformation >( data, c, explicitGeometry( ) );
    BOOST_CHECK_THROW( base->calculateDisplacementVelocity( epoch, s ), std::runtime_error );
    BOOST_CHECK_THROW( base->calculateDisplacementVelocity( epoch, p ), std::runtime_error );
    BOOST_CHECK_THROW( base->calculateDisplacementVelocity( epoch.getTdbSeconds( ), s ), std::runtime_error );
    BOOST_CHECK_THROW( base->calculateDisplacementVelocity( epoch.getTdbSeconds( ), p ), std::runtime_error );
    BOOST_CHECK_THROW( base->calculateDisplacementState( epoch, s ), std::runtime_error );
    BOOST_CHECK_THROW( base->calculateDisplacementState( epoch, p ), std::runtime_error );
    BOOST_CHECK_THROW( base->calculateDisplacementState( epoch.getTdbSeconds( ), s ), std::runtime_error );
    BOOST_CHECK_THROW( base->calculateDisplacementState( epoch.getTdbSeconds( ), p ), std::runtime_error );
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
    const EarthDeformationEpoch epoch( tdbAtUt1( c, Time( ( 57754. - 51544.5 ) * 86400. ) ), c );
    BOOST_CHECK_THROW( approximate.calculateDisplacement( epoch, Eigen::Vector3d( 0, 0, 6378137 ) ), std::runtime_error );
    AtmosphericTideDeformation regional( grid( ), c, explicitGeometry( ) );
    BOOST_CHECK_THROW( regional.calculateDisplacement( epoch, Eigen::Vector3d( 0, -6378137, 0 ) ), std::runtime_error );
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
    const auto file = createBodyDeformationModel( atmosphericTideBodyShapeDeformationFromGridFile( path, {}, settings ), "Earth", bodies );
    const double lat = -19.5 * mathematical_constants::PI / 180, lon = 47.5 * mathematical_constants::PI / 180;
    const Eigen::Vector3d p =
            6378137 * Eigen::Vector3d( std::cos( lat ) * std::cos( lon ), std::cos( lat ) * std::sin( lon ), std::sin( lat ) );
    const Time midnightTai( ( 57754. - 51544.5 ) * 86400. + 37 );
    for( double offset : { -.5, 0., .5 } )
    {
        const auto tdb = c->getCurrentTime< Time >( tai_scale, tdb_scale, midnightTai + offset );
        close( preloaded->calculateDisplacement( EarthDeformationEpoch( tdb, c ), p ),
               file->calculateDisplacement( EarthDeformationEpoch( tdb, c ), p ) );
    }
    const auto stations = readAtmosphericTideStationData( base + "bernese_manual_stations.atl" );
    const auto stationPreloaded =
            createBodyDeformationModel( atmosphericTideBodyShapeDeformationFromStations( stations, settings ), "Earth", bodies );
    const auto stationFile = createBodyDeformationModel(
            atmosphericTideBodyShapeDeformationFromStationFile( base + "bernese_manual_stations.atl", {}, settings ), "Earth", bodies );
    const auto namedStation = station( p, nullptr, stations.stations.begin( )->first );
    const EarthDeformationEpoch epoch( c->getCurrentTime< Time >( tai_scale, tdb_scale, midnightTai ), c );
    close( stationPreloaded->calculateDisplacement( epoch, namedStation ), stationFile->calculateDisplacement( epoch, namedStation ) );
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
