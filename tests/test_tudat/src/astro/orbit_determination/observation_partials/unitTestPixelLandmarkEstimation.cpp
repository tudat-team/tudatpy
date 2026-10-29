/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

// Estimation tests for the SUM/LMK pixel-landmark observable, through the OrbitDeterminationManager: recovery of per-image camera
// pointing corrections with the trajectory fixed, and of the spacecraft initial state together with the pointing corrections
// (without and with the SIGMA_PTG a-priori), on synthetic data; and state and joint state/pointing estimation on real Rosetta
// SUM/LMK data of comet 67P.

#define BOOST_TEST_MAIN

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include <boost/test/included/unit_test.hpp>

#include "tudat/astro/basic_astro/keplerPropagator.h"
#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/ephemerides/constantEphemeris.h"
#include "tudat/astro/ephemerides/constantRotationalEphemeris.h"
#include "tudat/astro/ephemerides/simpleRotationalEphemeris.h"
#include "tudat/astro/ephemerides/tabulatedEphemeris.h"
#include "tudat/astro/ephemerides/tabulatedRotationalEphemeris.h"
#include "tudat/astro/gravitation/gravityFieldModel.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/cameraPointingCorrection.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/io/readSumLmkFiles.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/math/integrators/createNumericalIntegrator.h"
#include "tudat/math/interpolators/cubicSplineInterpolator.h"
#include "tudat/simulation/estimation_setup/createEstimatableParametersFactory.h"
#include "tudat/simulation/estimation_setup/createObservationModelFactory.h"
#include "tudat/simulation/estimation_setup/createObservationPartials.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"
#include "tudat/support/sumLmkTestHelpers.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"
#include "tudat/simulation/propagation_setup/createAccelerationModels.h"
#include "tudat/simulation/propagation_setup/propagationSettings.h"
#include "tudat/simulation/propagation_setup/propagationTerminationSettings.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::ephemerides;
using namespace tudat::observation_models;
using namespace tudat::simulation_setup;
using namespace tudat::observation_partials;
using namespace tudat::estimatable_parameters;
using namespace tudat::numerical_integrators;
using namespace tudat::propagators;
using namespace tudat::orbital_element_conversions;
using namespace tudat::basic_astrodynamics;

namespace
{

//! Synthetic target gravitational parameter [m^3 s^-2]: chosen (not 67P) so the orbital period is a
//! few thousand seconds, giving good velocity observability over a short test arc.
constexpr double targetGravitationalParameter = 1.0E6;

//! Landmarks spread across the field of view so that all three rotation degrees of freedom
//! (including rotation about the boresight) are observable.
std::map< std::string, input_output::sum_lmk::LmkLandmarkData > makeLandmarks( )
{
    const std::map< std::string, Eigen::Vector3d > positions = { { "L1", ( Eigen::Vector3d( ) << 600.0, 200.0, 50.0 ).finished( ) },
                                                                 { "L2", ( Eigen::Vector3d( ) << -500.0, 400.0, -30.0 ).finished( ) },
                                                                 { "L3", ( Eigen::Vector3d( ) << 300.0, -550.0, 20.0 ).finished( ) },
                                                                 { "L4", ( Eigen::Vector3d( ) << -400.0, -300.0, -60.0 ).finished( ) } };
    std::map< std::string, input_output::sum_lmk::LmkLandmarkData > landmarks;
    for( const auto& entry : positions )
    {
        input_output::sum_lmk::LmkLandmarkData landmark;
        landmark.landmarkId_ = entry.first;
        landmark.bodyFixedPosition_ = entry.second;
        landmarks[ entry.first ] = landmark;
    }
    return landmarks;
}

//! Build an orthonormal body-fixed -> camera rotation whose +Z (boresight) points from the
//! spacecraft towards the target centre, so landmarks near the centre always project in front of
//! the camera. Rows are CX, CY, CZ as in the SUM camera-axes convention.
Eigen::Matrix3d boresightCameraAxes( const Eigen::Vector3d& spacecraftBodyFixedPosition )
{
    const Eigen::Vector3d cameraZ = ( -spacecraftBodyFixedPosition ).normalized( );
    Eigen::Vector3d reference = Eigen::Vector3d::UnitZ( );
    if( std::fabs( reference.dot( cameraZ ) ) > 0.95 )
    {
        reference = Eigen::Vector3d::UnitX( );
    }
    const Eigen::Vector3d cameraX = reference.cross( cameraZ ).normalized( );
    const Eigen::Vector3d cameraY = cameraZ.cross( cameraX ).normalized( );
    Eigen::Matrix3d cameraAxes;
    cameraAxes.row( 0 ) = cameraX.transpose( );
    cameraAxes.row( 1 ) = cameraY.transpose( );
    cameraAxes.row( 2 ) = cameraZ.transpose( );
    return cameraAxes;
}

//! Format a "YYYY MON DD HH:MM:SS.sss" UTC string for a whole-second offset from 2015 JUN 05
//! 06:00:00, kept within the same day (the arc is short) so no date/leap-second handling is needed.
std::string makeUtcString( const int secondsOffsetFromBase )
{
    const int baseSecondsOfDay = 6 * 3600;
    const int totalSeconds = baseSecondsOfDay + secondsOffsetFromBase;
    const int hours = totalSeconds / 3600;
    const int minutes = ( totalSeconds % 3600 ) / 60;
    const int seconds = totalSeconds % 60;
    char buffer[ 64 ];
    std::snprintf( buffer, sizeof( buffer ), "2015 JUN 05 %02d:%02d:%02d.000", hours, minutes, seconds );
    return std::string( buffer );
}

//! Body-fixed landmark positions on the target [m], spread within a small radius so the full set
//! stays in the field of view across the synthetic orbit arc while still spanning three dimensions.
std::map< std::string, input_output::sum_lmk::LmkLandmarkData > makeOrbitLandmarks( )
{
    const std::map< std::string, Eigen::Vector3d > positions = { { "LMK01", ( Eigen::Vector3d( ) << 700.0, 200.0, 150.0 ).finished( ) },
                                                                 { "LMK02", ( Eigen::Vector3d( ) << -600.0, 300.0, -100.0 ).finished( ) },
                                                                 { "LMK03", ( Eigen::Vector3d( ) << 200.0, -650.0, 250.0 ).finished( ) },
                                                                 { "LMK04", ( Eigen::Vector3d( ) << -300.0, -400.0, -200.0 ).finished( ) },
                                                                 { "LMK05", ( Eigen::Vector3d( ) << 500.0, 500.0, 300.0 ).finished( ) },
                                                                 { "LMK06", ( Eigen::Vector3d( ) << -450.0, 150.0, 400.0 ).finished( ) } };
    std::map< std::string, input_output::sum_lmk::LmkLandmarkData > landmarks;
    for( const auto& entry : positions )
    {
        input_output::sum_lmk::LmkLandmarkData landmark;
        landmark.landmarkId_ = entry.first;
        landmark.bodyFixedPosition_ = entry.second;
        landmarks[ entry.first ] = landmark;
    }
    return landmarks;
}

//! Compute the per-component RMS of a concatenated residual vector.
double residualRms( const Eigen::VectorXd& residuals )
{
    return ( residuals.size( ) > 0 ) ? std::sqrt( residuals.squaredNorm( ) / static_cast< double >( residuals.size( ) ) ) : 0.0;
}

//! Load the committed comet body-fixed attitude (sampled from the SPICE 67P/C-G_CK frame over the
//! data arc) into a tabulated rotational ephemeris. File columns: epoch [TDB s since J2000],
//! quaternion (w,x,y,z) from the comet body-fixed frame to J2000, and body-fixed angular velocity.
std::shared_ptr< RotationalEphemeris > loadCometAttitude( const std::string& attitudeFile )
{
    std::ifstream stream( attitudeFile );
    if( !stream.is_open( ) )
    {
        throw std::runtime_error( "Could not open comet attitude file: " + attitudeFile );
    }
    std::map< double, Eigen::Matrix< double, 7, 1 > > rotationalStateMap;
    double epoch;
    while( stream >> epoch )
    {
        Eigen::Matrix< double, 7, 1 > rotationalState;
        for( int i = 0; i < 7; ++i )
        {
            stream >> rotationalState( i );
        }
        rotationalStateMap[ epoch ] = rotationalState;
    }
    const std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< double, 7, 1 > > > interpolator =
            std::make_shared< interpolators::CubicSplineInterpolator< double, Eigen::Matrix< double, 7, 1 > > >( rotationalStateMap );
    return std::make_shared< TabulatedRotationalEphemeris< double, double > >( interpolator, "J2000", "Comet_Fixed" );
}

//! A synthetic Keplerian orbit arc imaged over roughly one period, with one SUM image per epoch and a
//! boresight pointed at the target centre, and SIGMA_PTG pointingSigmaRadians for each image.
struct SyntheticOrbitScenario {
    SystemOfBodies bodies_;
    std::map< std::string, input_output::sum_lmk::LmkLandmarkData > landmarks_;
    LoadedPixelTestData conversionResult_;
    std::shared_ptr< TranslationalStatePropagatorSettings< double, double > > propagatorSettings_;
    Eigen::Vector6d truthInitialState_;
    std::vector< std::string > imageIds_;
};

SyntheticOrbitScenario buildSyntheticOrbitScenario( const double pointingSigmaRadians )
{
    const int numberOfImages = 24;
    SyntheticOrbitScenario scenario;

    SystemOfBodies bodies( "SSB", "J2000" );
    bodies.createEmptyBody< double, double >( "Target", false );
    bodies.createEmptyBody< double, double >( "Spacecraft", false );

    bodies.at( "Target" )->setEphemeris( std::make_shared< ConstantEphemeris >( Eigen::Vector6d::Zero( ), "SSB", "J2000" ) );
    bodies.at( "Target" )->setGravityFieldModel( std::make_shared< gravitation::GravityFieldModel >( targetGravitationalParameter ) );
    bodies.at( "Target" )
            ->setRotationalEphemeris( std::make_shared< SimpleRotationalEphemeris >(
                    0.3, 1.1, 0.2, 2.0 * mathematical_constants::PI / 12000.0, 0.0, "J2000", "Target_Fixed" ) );
    bodies.at( "Spacecraft" )
            ->setRotationalEphemeris( std::make_shared< ConstantRotationalEphemeris >(
                    Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ), "J2000", "Spacecraft_Fixed" ) );
    bodies.at( "Spacecraft" )
            ->setEphemeris( std::make_shared< TabulatedCartesianEphemeris<> >(
                    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector6d > >( ), "SSB", "J2000" ) );
    bodies.processBodyFrameDefinitions< double, double >( );

    Eigen::Vector6d truthKeplerianElements = Eigen::Vector6d::Zero( );
    truthKeplerianElements( semiMajorAxisIndex ) = 1.0E4;
    truthKeplerianElements( eccentricityIndex ) = 0.05;
    truthKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 30.0 );
    truthKeplerianElements( argumentOfPeriapsisIndex ) = unit_conversions::convertDegreesToRadians( 40.0 );
    truthKeplerianElements( longitudeOfAscendingNodeIndex ) = unit_conversions::convertDegreesToRadians( 25.0 );
    truthKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 10.0 );

    scenario.landmarks_ = makeOrbitLandmarks( );
    std::vector< std::string > landmarkIds;
    for( const auto& entry : scenario.landmarks_ )
    {
        landmarkIds.push_back( entry.first );
    }

    const int epochSpacingSeconds = 250;
    std::vector< input_output::sum_lmk::SumImageData > images;
    std::vector< double > epochs;
    for( int imageIndex = 0; imageIndex < numberOfImages; ++imageIndex )
    {
        input_output::sum_lmk::SumImageData image;
        image.imageId_ = "IMG" + std::to_string( imageIndex );
        image.utcEpochString_ = makeUtcString( imageIndex * epochSpacingSeconds );
        image.imageSize_ = Eigen::Vector2i( 1024, 1024 );
        image.focalLengthMm_ = 100.0;
        image.opticalCenter_ = Eigen::Vector2d( 512.0, 512.0 );
        image.kMatrix_ << 10.0, 0.0, 0.0, 0.0, 10.0, 0.0;
        image.pointingSigma_ = Eigen::Vector3d::Constant( pointingSigmaRadians );

        const double epoch = sumImageTdbTimeForTest( image );
        epochs.push_back( epoch );
        const double timeSinceFirstEpoch = epoch - epochs.front( );
        const Eigen::Vector6d keplerianAtEpoch =
                propagateKeplerOrbit< double >( truthKeplerianElements, timeSinceFirstEpoch, targetGravitationalParameter );
        const Eigen::Vector6d inertialStateAtEpoch = convertKeplerianToCartesianElements( keplerianAtEpoch, targetGravitationalParameter );
        const Eigen::Matrix3d rotationInertialToBodyFixed =
                bodies.at( "Target" )->getRotationalEphemeris( )->getRotationToTargetFrame( epoch ).toRotationMatrix( );
        const Eigen::Vector3d spacecraftBodyFixedPosition = rotationInertialToBodyFixed * inertialStateAtEpoch.head( 3 );

        image.spacecraftObjectVector_ = -spacecraftBodyFixedPosition;
        image.cameraAxes_ = boresightCameraAxes( spacecraftBodyFixedPosition );

        for( const std::string& landmarkId : landmarkIds )
        {
            input_output::sum_lmk::SumLandmarkObservation observation;
            observation.landmarkId_ = landmarkId;
            observation.pixelCoordinates_ = Eigen::Vector2d::Zero( );  // overwritten by simulation
            image.landmarkObservations_.push_back( observation );
        }
        images.push_back( image );
        scenario.imageIds_.push_back( image.imageId_ );
    }

    scenario.truthInitialState_ = convertKeplerianToCartesianElements( truthKeplerianElements, targetGravitationalParameter );

    scenario.conversionResult_ = loadPixelTestData( images, scenario.landmarks_, bodies );

    SelectedAccelerationMap accelerationSettingsMap;
    accelerationSettingsMap[ "Spacecraft" ][ "Target" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    const std::vector< std::string > bodiesToIntegrate = { "Spacecraft" };
    const std::vector< std::string > centralBodies = { "Target" };
    const AccelerationMap accelerationModelMap =
            createAccelerationModelsMap( bodies, accelerationSettingsMap, bodiesToIntegrate, centralBodies );

    scenario.propagatorSettings_ = translationalStatePropagatorSettings< double, double >(
            centralBodies,
            accelerationModelMap,
            bodiesToIntegrate,
            scenario.truthInitialState_,
            epochs.front( ),
            rungeKuttaFixedStepSettings< double >( 10.0, CoefficientSets::rungeKuttaFehlberg78 ),
            propagationTimeTerminationSettings( epochs.back( ) + 100.0 ) );

    scenario.bodies_ = bodies;
    return scenario;
}

//! Get the start index of each image's pointing correction in the estimated parameter vector. The camera has one correction per
//! image, in time order.
std::vector< int > getPointingStartIndices( const std::shared_ptr< EstimatableParameterSet< double > >& parametersToEstimate,
                                            const int numberOfImages )
{
    const EstimatebleParameterIdentifier parameterIdentifier( camera_pointing_correction,
                                                              std::make_pair( "Spacecraft", sumLmkTestCameraName ) );
    const std::vector< std::pair< int, int > > indices = parametersToEstimate->getIndicesForParameterType( parameterIdentifier );
    BOOST_REQUIRE_EQUAL( indices.size( ), 1 );
    BOOST_REQUIRE_EQUAL( indices.at( 0 ).second, 3 * numberOfImages );
    std::vector< int > startIndices;
    for( int imageIndex = 0; imageIndex < numberOfImages; ++imageIndex )
    {
        startIndices.push_back( indices.at( 0 ).first + 3 * imageIndex );
    }
    return startIndices;
}

//! Simulate ideal pixel observations at the epochs of the converted SUM/LMK data, with unit weights.
std::shared_ptr< ObservationCollection< double, double > > simulatePixelObservations(
        const LoadedPixelTestData& conversionResult,
        const OrbitDeterminationManager< double, double >& orbitDeterminationManager,
        const SystemOfBodies& bodies )
{
    const std::shared_ptr< ObservationCollection< double, double > > simulatedObservations = simulateObservations< double, double >(
            getObservationSimulationSettingsFromObservations< double, double >( conversionResult.observationCollection_, bodies ),
            orbitDeterminationManager.getObservationSimulators( ),
            bodies );
    std::map< std::shared_ptr< ObservationCollectionParser >, double > weightsPerObservationParser;
    weightsPerObservationParser[ observationParser( pixel_coordinates ) ] = 1.0;
    simulatedObservations->setConstantWeightPerObservable( weightsPerObservationParser );
    return simulatedObservations;
}

//! A deterministic, distinct pointing offset per image, of the given magnitude [rad].
Eigen::Vector3d syntheticPointingOffset( const int imageIndex, const double magnitude )
{
    return magnitude *
            ( Eigen::Vector3d( ) << std::sin( 0.7 * imageIndex + 0.3 ),
              std::cos( 1.3 * imageIndex + 1.1 ),
              std::sin( 2.1 * imageIndex + 0.5 ) )
                    .finished( );
}

//! The committed real Rosetta/67P scenario: comet with point-mass gravity and the SPC body-fixed
//! attitude, spacecraft on the reconstructed SPICE orbiter arc, and the reduced SUM/LMK dataset
//! ingested into a pixel observation collection. Shared by the real-data estimation tests.
struct RealRosettaScenario {
    SystemOfBodies bodies_;
    LoadedPixelTestData conversionResult_;
    std::shared_ptr< TranslationalStatePropagatorSettings< double, double > > propagatorSettings_;
    Eigen::Vector6d initialStateGuess_;
};

RealRosettaScenario buildRealRosettaScenario( )
{
    RealRosettaScenario scenario;
    const std::string dataPath = paths::getTudatTestDataPath( ) + "/sum_lmk/real_67p";
    spice_interface::loadSpiceKernelInTudat( dataPath + "/rosetta_67p_orbiter_arc.bsp" );

    const double cometGravitationalParameter = 666.2;  // [m^3 s^-2], 67P/Churyumov-Gerasimenko.
    SystemOfBodies bodies( "SSB", "J2000" );
    bodies.createEmptyBody< double, double >( "Comet", false );
    bodies.createEmptyBody< double, double >( "Spacecraft", false );
    bodies.at( "Comet" )->setEphemeris( std::make_shared< ConstantEphemeris >( Eigen::Vector6d::Zero( ), "SSB", "J2000" ) );
    bodies.at( "Comet" )->setGravityFieldModel( std::make_shared< gravitation::GravityFieldModel >( cometGravitationalParameter ) );
    bodies.at( "Comet" )->setRotationalEphemeris( loadCometAttitude( dataPath + "/rosetta_67p_attitude_arc.txt" ) );
    bodies.at( "Spacecraft" )
            ->setRotationalEphemeris( std::make_shared< ConstantRotationalEphemeris >(
                    Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ), "J2000", "Spacecraft_Fixed" ) );
    bodies.at( "Spacecraft" )
            ->setEphemeris( std::make_shared< TabulatedCartesianEphemeris<> >(
                    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector6d > >( ), "SSB", "J2000" ) );
    bodies.processBodyFrameDefinitions< double, double >( );

    std::vector< std::string > sumFiles, lmkFiles;
    for( const std::filesystem::directory_entry& entry : std::filesystem::directory_iterator( dataPath ) )
    {
        const std::string ext = entry.path( ).extension( ).string( );
        if( ext == ".SUM" )
        {
            sumFiles.push_back( entry.path( ).string( ) );
        }
        else if( ext == ".LMK" )
        {
            lmkFiles.push_back( entry.path( ).string( ) );
        }
    }
    if( sumFiles.size( ) < 2 || lmkFiles.empty( ) )
    {
        throw std::runtime_error( "Real Rosetta SUM/LMK test data is incomplete." );
    }

    std::vector< input_output::sum_lmk::SumImageData > sumImages = input_output::sum_lmk::readSumFiles( sumFiles );
    std::sort( sumImages.begin( ), sumImages.end( ), []( const auto& a, const auto& b ) {
        return sumImageTdbTimeForTest( a ) < sumImageTdbTimeForTest( b );
    } );
    const double epoch0 = sumImageTdbTimeForTest( sumImages.front( ) );
    const double epochLast = sumImageTdbTimeForTest( sumImages.back( ) );

    // A-priori spacecraft state = reconstructed Rosetta orbiter state relative to the comet, from the
    // committed SPK (NAIF ids: Rosetta orbiter -226, comet 1000012), in metres.
    scenario.initialStateGuess_ = spice_interface::getBodyCartesianStateAtEpoch( "-226", "1000012", "J2000", "none", epoch0 );

    scenario.conversionResult_ = loadPixelTestData( sumImages, input_output::sum_lmk::readLmkFiles( lmkFiles ), bodies, "Comet" );

    SelectedAccelerationMap accelerationSettingsMap;
    accelerationSettingsMap[ "Spacecraft" ][ "Comet" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    const std::vector< std::string > bodiesToIntegrate = { "Spacecraft" };
    const std::vector< std::string > centralBodies = { "Comet" };
    const AccelerationMap accelerationModelMap =
            createAccelerationModelsMap( bodies, accelerationSettingsMap, bodiesToIntegrate, centralBodies );

    scenario.propagatorSettings_ = translationalStatePropagatorSettings< double, double >(
            centralBodies,
            accelerationModelMap,
            bodiesToIntegrate,
            scenario.initialStateGuess_,
            epoch0,
            rungeKuttaFixedStepSettings< double >( 30.0, CoefficientSets::rungeKuttaFehlberg78 ),
            propagationTimeTerminationSettings( epochLast + 600.0 ) );

    scenario.bodies_ = bodies;
    return scenario;
}

}  // namespace

BOOST_AUTO_TEST_SUITE( test_pixel_landmark_estimation )

//! Recover a distinct injected pointing correction for several images at once (given without SIGMA_PTG), with the trajectory held
//! fixed, through the OrbitDeterminationManager in observation-only mode. Each image's pointing is then determined by its own landmarks alone,
//! so the normal matrix must be block diagonal, with one 3x3 block per image.
BOOST_AUTO_TEST_CASE( testPointingRecovery )
{
    spice_interface::loadStandardSpiceKernels( );

    const std::vector< std::string > landmarkIds = { "L1", "L2", "L3", "L4" };
    const std::vector< std::string > imageIds = { "IMGA", "IMGB", "IMGC" };
    SystemOfBodies bodies = createSyntheticSumLmkBodies( );
    std::vector< input_output::sum_lmk::SumImageData > images;
    for( unsigned int imageIndex = 0; imageIndex < imageIds.size( ); ++imageIndex )
    {
        images.push_back( makeSyntheticSumImage( imageIds.at( imageIndex ), landmarkIds ) );
        images.back( ).utcEpochString_ = "2015 JUN 05 07:2" + std::to_string( imageIndex ) + ":42.053";
        images.back( ).pointingSigma_ = Eigen::Vector3d::Constant( std::numeric_limits< double >::quiet_NaN( ) );
    }
    LoadedPixelTestData conversionResult = loadPixelTestData( images, makeLandmarks( ), bodies );

    // Without SIGMA_PTG, no pointing a-priori is created.
    BOOST_CHECK( conversionResult.inverseAprioriCovarianceDiagonalEntries_.empty( ) );

    const std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double, double >( pointingSettingsForTest( conversionResult ), bodies );
    const std::shared_ptr< PropagatorSettings< double > > nullPropagatorSettings;
    OrbitDeterminationManager< double, double > orbitDeterminationManager(
            bodies, parametersToEstimate, conversionResult.observationModelSettings_, nullPropagatorSettings );
    const std::vector< int > pointingStartIndices = getPointingStartIndices( parametersToEstimate, imageIds.size( ) );

    // Simulate observations with a distinct correction per image, and estimate starting from zero corrections.
    Eigen::VectorXd truthParameters = parametersToEstimate->getFullParameterValues< double >( );
    for( unsigned int imageIndex = 0; imageIndex < imageIds.size( ); ++imageIndex )
    {
        truthParameters.segment( pointingStartIndices.at( imageIndex ), 3 ) = syntheticPointingOffset( imageIndex, 1.5E-3 );
    }
    parametersToEstimate->resetParameterValues( truthParameters );
    const std::shared_ptr< ObservationCollection< double, double > > simulatedObservations =
            simulatePixelObservations( conversionResult, orbitDeterminationManager, bodies );
    parametersToEstimate->resetParameterValues( Eigen::VectorXd( Eigen::VectorXd::Zero( truthParameters.rows( ) ) ) );

    const std::shared_ptr< EstimationInput< double, double > > estimationInput =
            std::make_shared< EstimationInput< double, double > >( simulatedObservations );
    estimationInput->defineEstimationSettings( true, true, true, false, true, false );
    estimationInput->setConvergenceChecker( std::make_shared< EstimationConvergenceChecker >( 8 ) );
    const std::shared_ptr< EstimationOutput< double > > estimationOutput = orbitDeterminationManager.estimateParameters( estimationInput );

    BOOST_CHECK_SMALL( residualRms( estimationOutput->residuals_ ), 1.0E-6 );
    BOOST_CHECK_SMALL( ( estimationOutput->parameterEstimate_ - truthParameters ).norm( ), 1.0E-9 );

    const Eigen::MatrixXd designMatrix = estimationOutput->getUnnormalizedDesignMatrix( );
    const Eigen::MatrixXd normalMatrix = designMatrix.transpose( ) * designMatrix;
    for( unsigned int row = 0; row < imageIds.size( ); ++row )
    {
        const Eigen::Matrix3d diagonalBlock = normalMatrix.block( pointingStartIndices.at( row ), pointingStartIndices.at( row ), 3, 3 );
        BOOST_CHECK_GT( diagonalBlock.norm( ), 0.0 );
        for( unsigned int column = 0; column < imageIds.size( ); ++column )
        {
            if( row != column )
            {
                BOOST_CHECK_SMALL( normalMatrix.block( pointingStartIndices.at( row ), pointingStartIndices.at( column ), 3, 3 ).norm( ) /
                                           diagonalBlock.norm( ),
                                   1.0E-12 );
            }
        }
    }
}


//! Estimate the spacecraft initial state and a per-image camera pointing correction together, through the full
//! OrbitDeterminationManager (propagation, variational equations and pixel partials), without and with the SIGMA_PTG a-priori.
//!
//! A camera rotation and a spacecraft translation perpendicular to the line of sight both shift the image; they are separated
//! only by the parallax across the landmarks in one image. Recovering both from noise-free data without a-priori therefore checks
//! that the pointing partials and the state partials are mutually consistent. The a-priori must then reach the normal equations
//! at the pointing parameters: it pulls the corrections towards zero and reduces their formal errors.
BOOST_AUTO_TEST_CASE( testJointStateAndPointingEstimation )
{
    spice_interface::loadStandardSpiceKernels( );
    const double pointingSigma = 1.0E-3;
    SyntheticOrbitScenario scenario = buildSyntheticOrbitScenario( pointingSigma );

    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames = pointingSettingsForTest( scenario.conversionResult_ );
    parameterNames.insert( parameterNames.begin( ),
                           std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
                                   "Spacecraft", scenario.truthInitialState_, "Target" ) );
    const std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double, double >( parameterNames, scenario.bodies_, scenario.propagatorSettings_ );
    OrbitDeterminationManager< double, double > orbitDeterminationManager(
            scenario.bodies_, parametersToEstimate, scenario.conversionResult_.observationModelSettings_, scenario.propagatorSettings_ );
    const std::vector< int > pointingStartIndices = getPointingStartIndices( parametersToEstimate, scenario.imageIds_.size( ) );

    // Simulate observations with the truth state and a distinct correction per image.
    Eigen::VectorXd truthParameters = parametersToEstimate->getFullParameterValues< double >( );
    for( unsigned int imageIndex = 0; imageIndex < pointingStartIndices.size( ); ++imageIndex )
    {
        truthParameters.segment( pointingStartIndices.at( imageIndex ), 3 ) = syntheticPointingOffset( imageIndex, 1.0E-3 );
    }
    parametersToEstimate->resetParameterValues( truthParameters );
    const std::shared_ptr< ObservationCollection< double, double > > simulatedObservations =
            simulatePixelObservations( scenario.conversionResult_, orbitDeterminationManager, scenario.bodies_ );

    // The a-priori created from the SUM SIGMA_PTG values: 1/sigma^2 on the pointing parameters only.
    const Eigen::MatrixXd inverseAprioriCovariance = createCovarianceFromDiagonalEntries(
            parametersToEstimate, scenario.conversionResult_.inverseAprioriCovarianceDiagonalEntries_, false );
    Eigen::VectorXd expectedDiagonal = Eigen::VectorXd::Zero( truthParameters.rows( ) );
    for( const int startIndex : pointingStartIndices )
    {
        expectedDiagonal.segment( startIndex, 3 ).setConstant( 1.0 / ( pointingSigma * pointingSigma ) );
    }
    BOOST_CHECK_SMALL( ( inverseAprioriCovariance - Eigen::MatrixXd( expectedDiagonal.asDiagonal( ) ) ).norm( ) / expectedDiagonal.norm( ),
                       1.0E-12 );

    // Perturb the state and zero all pointing corrections, then estimate without and with the a-priori.
    Eigen::VectorXd perturbedParameters = Eigen::VectorXd::Zero( truthParameters.rows( ) );
    perturbedParameters.segment( 0, 6 ) = truthParameters.segment( 0, 6 );
    perturbedParameters.segment( 0, 3 ) += Eigen::Vector3d::Constant( 10.0 );
    perturbedParameters.segment( 3, 3 ) += Eigen::Vector3d::Constant( 0.01 );
    std::vector< std::shared_ptr< EstimationOutput< double > > > estimationOutputs;
    for( const Eigen::MatrixXd& aprioriMatrix : { Eigen::MatrixXd( Eigen::MatrixXd::Zero( 0, 0 ) ), inverseAprioriCovariance } )
    {
        parametersToEstimate->resetParameterValues( perturbedParameters );
        const std::shared_ptr< EstimationInput< double, double > > estimationInput =
                std::make_shared< EstimationInput< double, double > >( simulatedObservations, aprioriMatrix );
        estimationInput->defineEstimationSettings( true, true, true, true, true, true );
        estimationInput->setConvergenceChecker( std::make_shared< EstimationConvergenceChecker >( 10 ) );
        estimationOutputs.push_back( orbitDeterminationManager.estimateParameters( estimationInput ) );
    }
    const std::shared_ptr< EstimationOutput< double > > unconstrainedOutput = estimationOutputs.at( 0 );
    const std::shared_ptr< EstimationOutput< double > > constrainedOutput = estimationOutputs.at( 1 );

    // Without a-priori, the truth is recovered.
    const Eigen::VectorXd estimationError = unconstrainedOutput->parameterEstimate_ - truthParameters;
    BOOST_CHECK_SMALL( residualRms( unconstrainedOutput->residuals_ ), 1.0E-6 );
    for( unsigned int i = 0; i < 3; ++i )
    {
        BOOST_CHECK_SMALL( std::fabs( estimationError( i ) ), 1.0E-5 );
        BOOST_CHECK_SMALL( std::fabs( estimationError( i + 3 ) ), 1.0E-8 );
    }
    for( const int startIndex : pointingStartIndices )
    {
        BOOST_CHECK_SMALL( estimationError.segment( startIndex, 3 ).norm( ), 1.0E-10 );
    }

    // With a-priori, the corrections are pulled away from the truth, with smaller formal errors and larger residuals.
    double constrainedDeviation = 0.0;
    for( const int startIndex : pointingStartIndices )
    {
        constrainedDeviation +=
                ( constrainedOutput->parameterEstimate_.segment( startIndex, 3 ) - truthParameters.segment( startIndex, 3 ) ).squaredNorm( );
        for( int component = 0; component < 3; ++component )
        {
            BOOST_CHECK_LT( constrainedOutput->getFormalErrorVector( )( startIndex + component ),
                            unconstrainedOutput->getFormalErrorVector( )( startIndex + component ) );
        }
    }
    BOOST_CHECK_GT( std::sqrt( constrainedDeviation ), 1.0E-6 );
    BOOST_CHECK_GT( residualRms( constrainedOutput->residuals_ ), residualRms( unconstrainedOutput->residuals_ ) );
}

//! Estimation from real Rosetta SUM/LMK data (2015-07-15, comet 67P), propagating the spacecraft under comet point-mass
//! gravity from the reconstructed SPICE orbiter state, with the comet attitude from the SPICE 67P/C-G_CK frame.
//!
//! First, the state alone is estimated: since the SPICE and SPC reconstructions are independent (~1 km apart), the residuals
//! must drop from the a-priori level to the few-pixel level of the SPC solution. Second, per-image pointing corrections are
//! estimated as well, with the SIGMA_PTG a-priori: this must improve the fit further, with bounded corrections. On these data
//! the corrections are nearly degenerate with the spacecraft position (the ~2 km comet seen from ~165 km gives ~1% parallax),
//! so their values are constrained mainly by the a-priori, and are not compared to a reference.
BOOST_AUTO_TEST_CASE( testRealRosettaEstimation )
{
    spice_interface::loadStandardSpiceKernels( );

    // Estimate the state alone, and the state with the pointing corrections, each from its own (unmodified) scenario.
    Eigen::Vector3d aprioriPosition;
    std::vector< std::shared_ptr< EstimationOutput< double > > > estimationOutputs;
    for( const bool estimatePointing : { false, true } )
    {
        RealRosettaScenario scenario = buildRealRosettaScenario( );
        aprioriPosition = scenario.initialStateGuess_.segment( 0, 3 );
        std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
        if( estimatePointing )
        {
            parameterNames = pointingSettingsForTest( scenario.conversionResult_ );
        }
        parameterNames.insert( parameterNames.begin( ),
                               std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
                                       "Spacecraft", scenario.initialStateGuess_, "Comet" ) );
        const std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate =
                createParametersToEstimate< double, double >( parameterNames, scenario.bodies_, scenario.propagatorSettings_ );
        OrbitDeterminationManager< double, double > orbitDeterminationManager( scenario.bodies_,
                                                                               parametersToEstimate,
                                                                               scenario.conversionResult_.observationModelSettings_,
                                                                               scenario.propagatorSettings_ );

        std::map< std::shared_ptr< ObservationCollectionParser >, double > weightsPerObservationParser;
        weightsPerObservationParser[ observationParser( pixel_coordinates ) ] = 1.0;
        scenario.conversionResult_.observationCollection_->setConstantWeightPerObservable( weightsPerObservationParser );

        Eigen::MatrixXd inverseAprioriCovariance = Eigen::MatrixXd::Zero( 0, 0 );
        if( estimatePointing )
        {
            inverseAprioriCovariance = createCovarianceFromDiagonalEntries(
                    parametersToEstimate, scenario.conversionResult_.inverseAprioriCovarianceDiagonalEntries_, false );
        }
        const std::shared_ptr< EstimationInput< double, double > > estimationInput = std::make_shared< EstimationInput< double, double > >(
                scenario.conversionResult_.observationCollection_, inverseAprioriCovariance );
        estimationInput->defineEstimationSettings( true, true, true, true, true, true );
        estimationInput->setConvergenceChecker( std::make_shared< EstimationConvergenceChecker >( 8 ) );
        estimationOutputs.push_back( orbitDeterminationManager.estimateParameters( estimationInput ) );
    }

    const std::shared_ptr< EstimationOutput< double > > stateOnlyOutput = estimationOutputs.at( 0 );
    const double initialRms = residualRms( stateOnlyOutput->getResidualHistoryMatrix( ).col( 0 ) );
    const double stateOnlyRms = residualRms( stateOnlyOutput->residuals_ );

    const std::shared_ptr< EstimationOutput< double > > jointOutput = estimationOutputs.at( 1 );
    const double jointRms = residualRms( jointOutput->residuals_ );
    double largestCorrection = 0.0;
    for( unsigned int i = 6; i < jointOutput->parameterEstimate_.rows( ); i += 3 )
    {
        largestCorrection = std::max( largestCorrection, jointOutput->parameterEstimate_.segment( i, 3 ).norm( ) );
    }

    BOOST_TEST_MESSAGE( "Real-data residual RMS: a-priori " << initialRms << " px, state-only " << stateOnlyRms << " px, joint " << jointRms
                                                            << " px; largest pointing correction " << largestCorrection << " rad" );

    BOOST_CHECK_LT( stateOnlyRms, 3.0 );
    BOOST_CHECK_LT( stateOnlyRms, 0.1 * initialRms );
    BOOST_CHECK_LT( jointRms, stateOnlyRms );
    BOOST_CHECK_LT( largestCorrection, 5.0E-3 );
    BOOST_CHECK_LT( ( stateOnlyOutput->parameterEstimate_.segment( 0, 3 ) - aprioriPosition ).norm( ), 5.0E3 );
    BOOST_CHECK_LT( ( jointOutput->parameterEstimate_.segment( 0, 3 ) - aprioriPosition ).norm( ), 5.0E3 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
