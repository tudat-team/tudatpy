#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include <limits>
#include <Eigen/LU>
#include "tudat/io/basicInputOutput.h"
#include "tudat/astro/reference_frames/terrestrialFrameTransformation.h"
#include "tudat/astro/reference_frames/referenceFrameTransformations.h"
#include "pr1092/reference_files/itrf_transformations/itrf_reference.h"
using namespace tudat;
using namespace reference_frames;
using namespace input_output;
namespace
{

//! Compare positions and velocities with the independent reference values.
void close( const Eigen::Vector6d& actual, const Eigen::Vector6d& expected )
{

    // Compare Cartesian positions in metres and velocities in metres per second.
    for( int i = 0; i < 6; ++i )
    {
        BOOST_CHECK_SMALL( actual[ i ] - expected[ i ], i < 3 ? 2e-9 : 2e-17 );
    }
}

//! Load the frozen Cartesian input state used by both published catalogues.
Eigen::Vector6d state( )
{
    return Eigen::Map< const Eigen::Vector6d >( itrfReferenceInputState );
}
}  // namespace

//! Compare both published ITRF catalogues with independent forward and inverse state references.
BOOST_AUTO_TEST_CASE( allPublishedRowsIndependentHighPrecisionForwardAndInverseStates )
{

    // Read the published tables used to generate the reference states.
    const auto d2014 = readTerrestrialFrameTransformationData( tudat::paths::getTudatTestDataPath( ) +
                                                               "pr1092/reference_files/itrf_transformations/Transfo-ITRF2014_ITRFs.txt" );

    const auto d2020 = readTerrestrialFrameTransformationData( tudat::paths::getTudatTestDataPath( ) +
                                                               "pr1092/reference_files/itrf_transformations/Transfo-ITRF2020_TRFs.txt" );

    for( const auto& fixture : itrfReferenceFixtures )
    {
        const auto& data = fixture.catalogue == 2014 ? d2014 : d2020;

        const auto transformation = evaluateLinearizedHelmert( data.getParameters( fixture.frame ), fixture.year, fixture.yearRate );

        close( transformation.transformState( state( ) ), Eigen::Map< const Eigen::Vector6d >( fixture.forward ) );
        close( transformation.inverse( ).transformState( state( ) ), Eigen::Map< const Eigen::Vector6d >( fixture.inverse ) );
    }
}

//! Check the sign, units and time-rate contribution of all seven Helmert parameters.
BOOST_AUTO_TEST_CASE( individualSevenParametersAndPhysicalRateChain )
{
    TerrestrialFrameTransformationParameters p;
    p.sourceFrame.realization = "A";
    p.targetFrame.realization = "B";
    p.referenceDecimalYear = 2000;
    p.rotationConvention = HelmertRotationConvention::positionVector;
    const Eigen::Vector3d position( 2, 3, 5 );
    const Eigen::Vector3d effects[] = { Eigen::Vector3d( 1, 0, 0 ),  Eigen::Vector3d( 0, 1, 0 ),  Eigen::Vector3d( 0, 0, 1 ), position,
                                        Eigen::Vector3d( 0, -5, 3 ), Eigen::Vector3d( 5, 0, -2 ), Eigen::Vector3d( -3, 2, 0 ) };
    for( int k = 0; k < 7; ++k )
    {
        auto q = p;
        const double amplitude = 1e-3;
        if( k < 3 )
        {
            q.translationMetres[ k ] = amplitude;
        }
        else if( k == 3 )
        {
            q.scale = amplitude;
        }
        else
        {
            q.rotationRadians[ k - 4 ] = amplitude;
        }
        const auto value = evaluateLinearizedHelmert( q, 2000, 1 / 31557600. );

        // The seven independent parameter effects have known Cartesian signs and magnitudes.
        BOOST_CHECK_SMALL( ( value.transformPosition( position ) - position - amplitude * effects[ k ] ).norm( ), 2e-15 );
        auto rate = p;
        if( k < 3 )
        {
            rate.translationRateMetresPerSecond[ k ] = amplitude / 31557600.;
        }
        else if( k == 3 )
        {
            rate.scaleRatePerSecond = amplitude / 31557600.;
        }
        else
        {
            rate.rotationRateRadiansPerSecond[ k - 4 ] = amplitude / 31557600.;
        }

        // Parameter rates contribute to both the displaced position and the Cartesian velocity.
        const auto changing = evaluateLinearizedHelmert( rate, 2002, 2 / 31557600. );

        BOOST_CHECK_SMALL( ( changing.transformPosition( position ) - position - 2 * amplitude * effects[ k ] ).norm( ), 2e-15 );
        const auto stationary = ( Eigen::Vector6d( ) << position, Eigen::Vector3d::Zero( ) ).finished( );
        BOOST_CHECK_SMALL( ( changing.transformState( stationary ).tail< 3 >( ) - ( 2 * amplitude / 31557600. ) * effects[ k ] ).norm( ),
                           1e-23 );

        // Coordinate-frame rotations reverse the position-vector rotation signs.
        if( k >= 4 )
        {
            rate.rotationConvention = HelmertRotationConvention::coordinateFrame;
            const auto reversed = evaluateLinearizedHelmert( rate, 2002, 2 / 31557600. );

            BOOST_CHECK_SMALL( ( reversed.transformPosition( position ) - position + 2 * amplitude * effects[ k ] ).norm( ), 2e-15 );
            BOOST_CHECK_SMALL(
                    ( reversed.transformState( stationary ).tail< 3 >( ) + ( 2 * amplitude / 31557600. ) * effects[ k ] ).norm( ), 1e-23 );
        }
        q.rotationConvention = HelmertRotationConvention::coordinateFrame;
        if( k >= 4 )
        {
            BOOST_CHECK_SMALL( ( evaluateLinearizedHelmert( q, 2000, 1 / 31557600. ).transformPosition( position ) - position +
                                 amplitude * effects[ k ] )
                                       .norm( ),
                               2e-15 );
        }
    }
}

//! Check published provenance, hub composition and the full matrix inverse.
BOOST_AUTO_TEST_CASE( modelProvenanceAndHubComposition )
{

    // Published catalogue and explicit year/rate preserve the original composition case.
    const auto data = readOfficialItrf2020TransformationData( );

    HelmertTransformationModel model( data );

    const double year = 2025.25;

    const double yearRate = 2 / 31557600.;

    const auto source = model.getTransformation( "ITRF93", "ITRF2020", year, yearRate );

    const auto target = model.getTransformation( "ITRF2020", "ITRF88", year, yearRate );

    close( model.transformState( state( ), "ITRF93", "ITRF88", year, yearRate ),
           target.transformState( source.transformState( state( ) ) ) );
    BOOST_CHECK_EQUAL( model.getData( ).referenceFrame.realization, "ITRF2020" );

    // The full matrix inverse, including scale, must recover position and velocity.
    const auto complete = model.getTransformation( "ITRF93", "ITRF88", year, yearRate );

    close( complete.inverse( ).transformState( complete.transformState( state( ) ) ), state( ) );
}

//! Check the corrected legacy table entries and transform moving states at the requested epoch.
BOOST_AUTO_TEST_CASE( legacyTranscriptionDefectsAndMovingStateAtRequestedEpoch )
{
    constexpr double year = 31557600.;
    const double mas = 3.1415926535897932384626433832795 / 648000000.;

    // Check the published translation and angular-rate entries corrected from the legacy transcription.
    BOOST_CHECK_SMALL( getItrf2014ToArbitraryItrfTranslation( "ITRF2005" )[ 3 ] * year - .0003, 1e-18 );
    BOOST_CHECK_SMALL( getItrf2014ToArbitraryItrfTranslation( "ITRF90" )[ 1 ] - .0115, 1e-18 );
    for( const auto frame : { "ITRF90", "ITRF89", "ITRF88" } )
    {
        BOOST_CHECK_SMALL( getItrf2014ToArbitraryItrfTranslation( frame )[ 4 ] * year + .0005, 1e-18 );
    }
    for( const auto frame : { "ITRF97", "ITRF96", "ITRF94" } )
    {
        BOOST_CHECK_SMALL( getItrf2014ToArbitraryItrfRotationMatrixDerivative( frame )( 1, 0 ) * year - .02 * mas, 1e-24 );
    }
    close( convertGroundStationStateBetweenItrfFrames( state( ), 20 * year, "ITRF2014", "ITRF2014" ), state( ) );
    const double scale = -.02e-9 + .03e-9 * 10;
    Eigen::Vector6d expected;
    expected.head< 3 >( ) = Eigen::Vector3d( .0016, .0019, .0014 ) + ( 1 + scale ) * state( ).head< 3 >( );
    expected.tail< 3 >( ) =
            Eigen::Vector3d( 0, 0, -.0001 / year ) + ( .03e-9 / year ) * state( ).head< 3 >( ) + ( 1 + scale ) * state( ).tail< 3 >( );
    close( convertGroundStationStateItrf2014ToArbitraryItrf( state( ), 20 * year, "ITRF2008" ), expected );
    close( convertGroundStationStateArbitraryItrfToItrf2014( expected, 20 * year, "ITRF2008" ), state( ) );
}

//! Recover position and velocity with the complete nonorthogonal transformation inverse.
BOOST_AUTO_TEST_CASE( exactNonorthogonalInverse )
{
    TerrestrialFrameAffineTransformation transformation;
    transformation.matrix = ( Eigen::Matrix3d( ) << 2, .1, .2, .3, 3, .4, .5, .6, 4 ).finished( );
    transformation.matrixRate.setConstant( .001 );
    transformation.translation << 1, 2, 3;
    transformation.translationRate << .01, .02, .03;
    const auto recovered = transformation.inverse( ).transformState( transformation.transformState( state( ) ) );
    const auto normInf = []( const Eigen::Matrix3d& matrix ) { return matrix.cwiseAbs( ).rowwise( ).sum( ).maxCoeff( ); };
    const double inverseNorm = normInf( transformation.matrix.inverse( ) );
    const double condition = normInf( transformation.matrix ) * inverseNorm;
    const double epsilon = std::numeric_limits< double >::epsilon( );

    // gamma_32 covers the small-matrix factorization and scalar products. Bound
    // cancellation using the input magnitudes, not the observed round-trip residual.
    const double gamma = 32 * epsilon / ( 1 - 32 * epsilon );
    const Eigen::Vector6d forward = transformation.transformState( state( ) );
    const double positionMagnitude =
            inverseNorm * ( forward.head< 3 >( ).lpNorm< Eigen::Infinity >( ) + transformation.translation.lpNorm< Eigen::Infinity >( ) ) +
            state( ).head< 3 >( ).lpNorm< Eigen::Infinity >( );
    const double velocityMagnitude = inverseNorm *
                    ( forward.tail< 3 >( ).lpNorm< Eigen::Infinity >( ) + transformation.translationRate.lpNorm< Eigen::Infinity >( ) +
                      normInf( transformation.matrixRate ) * state( ).head< 3 >( ).lpNorm< Eigen::Infinity >( ) ) +
            state( ).tail< 3 >( ).lpNorm< Eigen::Infinity >( );

    // Compare both recovered state blocks with conditioning-based floating-point error bounds.
    BOOST_CHECK_SMALL( ( recovered.head< 3 >( ) - state( ).head< 3 >( ) ).lpNorm< Eigen::Infinity >( ),
                       gamma * condition * positionMagnitude );
    BOOST_CHECK_SMALL( ( recovered.tail< 3 >( ) - state( ).tail< 3 >( ) ).lpNorm< Eigen::Infinity >( ),
                       gamma * condition * velocityMagnitude );
    BOOST_CHECK_GT( ( transformation.inverse( ).matrix - transformation.matrix.transpose( ) ).norm( ), 1 );
}
