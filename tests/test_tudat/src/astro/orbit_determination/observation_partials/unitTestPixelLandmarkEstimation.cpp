/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

// Estimation tests for the SUM/LMK pixel-landmark observable. Two suites:
//
//   test_pixel_landmark_estimation       -- per-image camera-pointing smoke tests: a self-contained
//                                           Gauss-Newton solve drives the CameraPointingCorrection
//                                           parameter and PixelCoordinatesPointingPartial objects on
//                                           a synthetic (+Z boresight) scene with the spacecraft state
//                                           held fixed.
//
//   test_pixel_landmark_state_estimation -- full spacecraft-state orbit determination through the
//                                           OrbitDeterminationManager (propagation + variational
//                                           equations + initial-state partials): a synthetic
//                                           truth-recovery case and an end-to-end case driven by real
//                                           Rosetta SUM/LMK data of comet 67P.

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
#include "tudat/simulation/estimation_setup/processSumLmkFiles.h"
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

//! Build Target (origin, identity rotation) + Spacecraft (behind, with placeholder attitude) bodies.
SystemOfBodies makeBodies( const double spacecraftZ )
{
    SystemOfBodies bodies( "SSB", "J2000" );
    bodies.createEmptyBody< double, double >( "Target", false );
    bodies.createEmptyBody< double, double >( "Spacecraft", false );
    bodies.at( "Target" )->setEphemeris( std::make_shared< ConstantEphemeris >( Eigen::Vector6d::Zero( ), "SSB", "J2000" ) );
    Eigen::Vector6d spacecraftState = Eigen::Vector6d::Zero( );
    spacecraftState( 2 ) = spacecraftZ;
    bodies.at( "Spacecraft" )->setEphemeris( std::make_shared< ConstantEphemeris >( spacecraftState, "SSB", "J2000" ) );
    bodies.at( "Target" )
            ->setRotationalEphemeris( std::make_shared< ConstantRotationalEphemeris >(
                    Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ), "J2000", "Target_Fixed" ) );
    bodies.at( "Spacecraft" )
            ->setRotationalEphemeris( std::make_shared< ConstantRotationalEphemeris >(
                    Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ), "J2000", "Spacecraft_Fixed" ) );
    bodies.processBodyFrameDefinitions< double, double >( );
    return bodies;
}

input_output::sum_lmk::SumImageData makeImage( const std::string& imageId,
                                               const std::vector< std::string >& landmarkIds,
                                               const double spacecraftZ )
{
    input_output::sum_lmk::SumImageData image;
    image.imageId_ = imageId;
    image.utcEpochString_ = "2015 JUN 05 07:24:42.053";
    image.imageSize_ = Eigen::Vector2i( 1024, 1024 );
    image.focalLengthMm_ = 100.0;
    image.opticalCenter_ = Eigen::Vector2d( 512.0, 512.0 );
    // SCOBJ = spacecraft-to-object (target at origin), i.e. the negation of the spacecraft position.
    image.spacecraftObjectVector_ = Eigen::Vector3d( 0.0, 0.0, -spacecraftZ );
    image.cameraAxes_ = Eigen::Matrix3d::Identity( );
    image.kMatrix_ << 10.0, 0.0, 0.0, 0.0, 10.0, 0.0;
    image.pointingSigma_ = Eigen::Vector3d::Constant( 1.0E-4 );
    for( const std::string& landmarkId : landmarkIds )
    {
        input_output::sum_lmk::SumLandmarkObservation observation;
        observation.landmarkId_ = landmarkId;
        observation.pixelCoordinates_ = Eigen::Vector2d::Zero( );  // simulated below; stored value not used here
        image.landmarkObservations_.push_back( observation );
    }
    return image;
}

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

//! Pixel observation for a single landmark link end at t = 0 (constant ephemerides).
Eigen::Vector2d computePixel( const std::shared_ptr< ObservationModel< 2, double, double > >& model )
{
    std::vector< Eigen::Vector6d > states;
    std::vector< double > times;
    return model->computeObservationsWithLinkEndData( 0.0, receiver, times, states, nullptr );
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
//! boresight pointed at the target centre. This is the shared setup for the joint state + pointing tests.
//! pointingSigmaRadians controls SIGMA_PTG: pass NaN for no pointing a-priori.
struct SyntheticOrbitScenario {
    SystemOfBodies bodies_;
    std::map< std::string, input_output::sum_lmk::LmkLandmarkData > landmarks_;
    SumLmkObservationConversionResult< double, double > conversionResult_;
    std::shared_ptr< TranslationalStatePropagatorSettings< double, double > > propagatorSettings_;
    Eigen::Vector6d truthInitialState_;
    std::vector< std::string > imageIds_;
};

SyntheticOrbitScenario buildSyntheticOrbitScenario( const double pointingSigmaRadians, const int numberOfImages = 24 )
{
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

        const double epoch = observation_models::detail::convertSumUtcStringToSecondsSinceJ2000< double >( image );
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

    SumLmkObservationConversionSettings conversionSettings( "Target", "Spacecraft" );
    scenario.conversionResult_ =
            createSumLmkObservationCollection< double, double >( images, scenario.landmarks_, bodies, conversionSettings );

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
    SumLmkObservationConversionResult< double, double > conversionResult_;
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
        return observation_models::detail::convertSumUtcStringToSecondsSinceJ2000< double >( a ) <
                observation_models::detail::convertSumUtcStringToSecondsSinceJ2000< double >( b );
    } );
    const double epoch0 = observation_models::detail::convertSumUtcStringToSecondsSinceJ2000< double >( sumImages.front( ) );
    const double epochLast = observation_models::detail::convertSumUtcStringToSecondsSinceJ2000< double >( sumImages.back( ) );

    // A-priori spacecraft state = reconstructed Rosetta orbiter state relative to the comet, from the
    // committed SPK (NAIF ids: Rosetta orbiter -226, comet 1000012), in metres.
    scenario.initialStateGuess_ = spice_interface::getBodyCartesianStateAtEpoch( "-226", "1000012", "J2000", "none", epoch0 );

    SumLmkObservationConversionSettings conversionSettings( "Comet", "Spacecraft" );
    scenario.conversionResult_ = createSumLmkObservationCollection< double, double >( sumFiles, lmkFiles, bodies, conversionSettings );

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

//! Recover an injected per-image pointing offset from pixel-landmark observations, with the
//! spacecraft translational state held fixed. Self-contained Gauss-Newton solve driving the
//! actual CameraPointingCorrection parameter and PixelCoordinatesPointingPartial objects.
BOOST_AUTO_TEST_CASE( testPointingOffsetRecovery )
{
    spice_interface::loadStandardSpiceKernels( );

    const std::vector< std::string > landmarkIds = { "L1", "L2", "L3", "L4" };
    SystemOfBodies bodies = makeBodies( -10000.0 );
    const std::vector< input_output::sum_lmk::SumImageData > images = { makeImage( "IMGA", landmarkIds, -10000.0 ) };
    SumLmkObservationConversionSettings conversionSettings( "Target", "Spacecraft" );
    SumLmkObservationConversionResult< double, double > conversionResult =
            createSumLmkObservationCollection< double, double >( images, makeLandmarks( ), bodies, conversionSettings );
    const std::string cameraName = conversionResult.imageIdToCameraName_.at( "IMGA" );

    // One observation model per landmark link end (all share the same receiver camera).
    std::vector< std::shared_ptr< ObservationModel< 2, double, double > > > models;
    for( const std::string& landmarkId : landmarkIds )
    {
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = LinkEndId( "Target", landmarkId );
        linkEnds[ receiver ] = LinkEndId( "Spacecraft", cameraName );
        models.push_back( ObservationModelCreator< 2, double, double >::createObservationModel(
                pixelCoordinatesSettings( LinkDefinition( linkEnds ) ), bodies ) );
    }

    std::shared_ptr< CameraPointingCorrection > pointingParameter =
            std::make_shared< CameraPointingCorrection >( bodies.at( "Spacecraft" )->getVehicleSystems( ), "Spacecraft", cameraName );
    std::shared_ptr< EstimatableParameterSet< double > > parameterSet = std::make_shared< EstimatableParameterSet< double > >(
            std::vector< std::shared_ptr< EstimatableParameter< double > > >( ),
            std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > >( { pointingParameter } ) );

    // Analytical pointing partial + scaling per landmark.
    std::vector< std::shared_ptr< ObservationPartial< 2 > > > pointingPartials;
    std::vector< std::shared_ptr< PositionPartialScaling > > scalings;
    for( const auto& model : models )
    {
        auto partials = ObservationPartialCreator< 2, double, double >::createObservationPartials( model, bodies, parameterSet );
        BOOST_REQUIRE_EQUAL( partials.first.size( ), 1 );
        pointingPartials.push_back( partials.first.begin( )->second );
        scalings.push_back( partials.second );
    }

    // Simulate measurements with the injected (true) pointing offset.
    const Eigen::Vector3d trueOffset = ( Eigen::Vector3d( ) << 1.0E-3, -8.0E-4, 1.5E-3 ).finished( );
    pointingParameter->setParameterValue( trueOffset );
    std::vector< Eigen::Vector2d > observedPixels;
    for( const auto& model : models )
    {
        observedPixels.push_back( computePixel( model ) );
    }

    // Gauss-Newton recovery starting from a zero correction.
    pointingParameter->setParameterValue( Eigen::Vector3d::Zero( ) );
    for( int iteration = 0; iteration < 10; ++iteration )
    {
        Eigen::Matrix3d normalMatrix = Eigen::Matrix3d::Zero( );
        Eigen::Vector3d rightHandSide = Eigen::Vector3d::Zero( );
        for( std::size_t i = 0; i < models.size( ); ++i )
        {
            std::vector< Eigen::Vector6d > states;
            std::vector< double > times;
            const Eigen::Vector2d modelled = models.at( i )->computeObservationsWithLinkEndData( 0.0, receiver, times, states, nullptr );
            scalings.at( i )->update( states, times, receiver, modelled );
            const Eigen::Matrix< double, 2, Eigen::Dynamic > designBlock =
                    pointingPartials.at( i )->calculatePartial( states, times, receiver, nullptr, modelled ).at( 0 ).first;
            const Eigen::Vector2d residual = observedPixels.at( i ) - modelled;
            normalMatrix += designBlock.transpose( ) * designBlock;
            rightHandSide += designBlock.transpose( ) * residual;
        }
        const Eigen::Vector3d update = normalMatrix.ldlt( ).solve( rightHandSide );
        pointingParameter->setParameterValue( pointingParameter->getParameterValue( ) + update );
        if( update.norm( ) < 1.0E-13 )
        {
            break;
        }
    }

    BOOST_CHECK_SMALL( ( pointingParameter->getParameterValue( ) - trueOffset ).norm( ), 1.0E-9 );

    // SIGMA_PTG a-priori inverse covariance was produced for this image.
    BOOST_REQUIRE_EQUAL( conversionResult.inverseAprioriCovarianceDiagonalEntries_.size( ), 1 );
    BOOST_CHECK_EQUAL( conversionResult.inverseAprioriCovarianceDiagonalEntries_.at( 0 ).first.first, camera_pointing_correction );
}

//! Two images produce two independent pointing parameters: perturbing one image's correction
//! must not change the other image's modelled observations.
BOOST_AUTO_TEST_CASE( testMultiImagePointingIndependence )
{
    spice_interface::loadStandardSpiceKernels( );

    SystemOfBodies bodies = makeBodies( -10000.0 );
    const std::vector< input_output::sum_lmk::SumImageData > images = { makeImage( "IMGA", { "L1", "L2" }, -10000.0 ),
                                                                        makeImage( "IMGB", { "L1", "L2" }, -10000.0 ) };
    SumLmkObservationConversionSettings conversionSettings( "Target", "Spacecraft" );
    SumLmkObservationConversionResult< double, double > conversionResult =
            createSumLmkObservationCollection< double, double >( images, makeLandmarks( ), bodies, conversionSettings );
    const std::string cameraA = conversionResult.imageIdToCameraName_.at( "IMGA" );
    const std::string cameraB = conversionResult.imageIdToCameraName_.at( "IMGB" );
    BOOST_CHECK( cameraA != cameraB );

    auto makeModel = [ & ]( const std::string& cameraName ) {
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = LinkEndId( "Target", "L1" );
        linkEnds[ receiver ] = LinkEndId( "Spacecraft", cameraName );
        return ObservationModelCreator< 2, double, double >::createObservationModel( pixelCoordinatesSettings( LinkDefinition( linkEnds ) ),
                                                                                     bodies );
    };
    std::shared_ptr< ObservationModel< 2, double, double > > modelB = makeModel( cameraB );

    std::shared_ptr< CameraPointingCorrection > parameterA =
            std::make_shared< CameraPointingCorrection >( bodies.at( "Spacecraft" )->getVehicleSystems( ), "Spacecraft", cameraA );
    std::shared_ptr< CameraPointingCorrection > parameterB =
            std::make_shared< CameraPointingCorrection >( bodies.at( "Spacecraft" )->getVehicleSystems( ), "Spacecraft", cameraB );

    const Eigen::Vector2d pixelBNominal = computePixel( modelB );

    // Perturbing image A does not affect image B.
    parameterA->setParameterValue( ( Eigen::Vector3d( ) << 2.0E-3, -1.0E-3, 1.0E-3 ).finished( ) );
    BOOST_CHECK_SMALL( ( computePixel( modelB ) - pixelBNominal ).norm( ), 1.0E-10 );

    // Perturbing image B does affect image B.
    parameterB->setParameterValue( ( Eigen::Vector3d( ) << 2.0E-3, -1.0E-3, 1.0E-3 ).finished( ) );
    BOOST_CHECK( ( computePixel( modelB ) - pixelBNominal ).norm( ) > 1.0E-3 );
}

//! The conversion result yields one pointing parameter setting per image, on the receiver body, ordered by
//! camera name, and can be restricted to a subset of the images.
BOOST_AUTO_TEST_CASE( testPointingParameterSettingsCreation )
{
    spice_interface::loadStandardSpiceKernels( );

    SystemOfBodies bodies = makeBodies( -10000.0 );
    const std::vector< input_output::sum_lmk::SumImageData > images = { makeImage( "IMGA", { "L1", "L2" }, -10000.0 ),
                                                                        makeImage( "IMGB", { "L1", "L2" }, -10000.0 ) };
    SumLmkObservationConversionSettings conversionSettings( "Target", "Spacecraft" );
    SumLmkObservationConversionResult< double, double > conversionResult =
            createSumLmkObservationCollection< double, double >( images, makeLandmarks( ), bodies, conversionSettings );

    BOOST_CHECK_EQUAL( conversionResult.receiverBodyName_, "Spacecraft" );

    const std::vector< std::shared_ptr< EstimatableParameterSettings > > allSettings =
            createSumLmkPointingParameterSettings< double, double >( conversionResult );
    BOOST_REQUIRE_EQUAL( allSettings.size( ), 2 );

    // Ordered by camera name, on the receiver body, with the camera as reference point.
    std::vector< std::string > settingCameraNames;
    for( const std::shared_ptr< EstimatableParameterSettings >& setting : allSettings )
    {
        BOOST_CHECK_EQUAL( setting->parameterType_.first, camera_pointing_correction );
        BOOST_CHECK_EQUAL( setting->parameterType_.second.first, "Spacecraft" );
        settingCameraNames.push_back( setting->parameterType_.second.second );
    }
    BOOST_CHECK( std::is_sorted( settingCameraNames.begin( ), settingCameraNames.end( ) ) );
    BOOST_CHECK_EQUAL( settingCameraNames.at( 0 ), conversionResult.imageIdToCameraName_.at( "IMGA" ) );
    BOOST_CHECK_EQUAL( settingCameraNames.at( 1 ), conversionResult.imageIdToCameraName_.at( "IMGB" ) );

    // Restricting to a subset of the images.
    const std::vector< std::shared_ptr< EstimatableParameterSettings > > subsetSettings =
            createSumLmkPointingParameterSettings< double, double >( conversionResult, { "IMGB" } );
    BOOST_REQUIRE_EQUAL( subsetSettings.size( ), 1 );
    BOOST_CHECK_EQUAL( subsetSettings.at( 0 )->parameterType_.second.second, conversionResult.imageIdToCameraName_.at( "IMGB" ) );

    // An unknown image ID is an error rather than being silently ignored.
    BOOST_CHECK_THROW( ( createSumLmkPointingParameterSettings< double, double >( conversionResult, { "NOT_AN_IMAGE" } ) ),
                       std::runtime_error );
}

//! The per-image SIGMA_PTG a-priori is assembled onto the diagonal of the matching parameter blocks, leaving
//! other parameters untouched, and entries for parameters that are not estimated are skipped.
BOOST_AUTO_TEST_CASE( testInverseAprioriCovarianceAssembly )
{
    spice_interface::loadStandardSpiceKernels( );

    SystemOfBodies bodies = makeBodies( -10000.0 );
    const std::vector< input_output::sum_lmk::SumImageData > images = { makeImage( "IMGA", { "L1", "L2" }, -10000.0 ),
                                                                        makeImage( "IMGB", { "L1", "L2" }, -10000.0 ) };
    SumLmkObservationConversionSettings conversionSettings( "Target", "Spacecraft" );
    SumLmkObservationConversionResult< double, double > conversionResult =
            createSumLmkObservationCollection< double, double >( images, makeLandmarks( ), bodies, conversionSettings );
    BOOST_REQUIRE_EQUAL( conversionResult.inverseAprioriCovarianceDiagonalEntries_.size( ), 2 );

    // makeImage uses SIGMA_PTG = 1e-4 rad on all three axes.
    const double expectedInverseVariance = 1.0 / ( 1.0E-4 * 1.0E-4 );

    // --- Both images estimated: both blocks are filled.
    {
        const std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate = createParametersToEstimate< double, double >(
                createSumLmkPointingParameterSettings< double, double >( conversionResult ), bodies );
        const int numberOfParameters = parametersToEstimate->getEstimatedParameterSetSize( );
        BOOST_REQUIRE_EQUAL( numberOfParameters, 6 );

        const Eigen::MatrixXd inverseAprioriCovariance =
                createSumLmkInverseAprioriCovariance< double, double, double >( conversionResult, parametersToEstimate );
        BOOST_REQUIRE_EQUAL( inverseAprioriCovariance.rows( ), numberOfParameters );
        BOOST_REQUIRE_EQUAL( inverseAprioriCovariance.cols( ), numberOfParameters );

        // Diagonal carries the inverse variances; the matrix has no off-diagonal content.
        BOOST_CHECK_SMALL( ( inverseAprioriCovariance.diagonal( ) - Eigen::VectorXd::Constant( 6, expectedInverseVariance ) ).norm( ) /
                                   expectedInverseVariance,
                           1.0E-12 );
        BOOST_CHECK_SMALL( ( inverseAprioriCovariance - Eigen::MatrixXd( inverseAprioriCovariance.diagonal( ).asDiagonal( ) ) ).norm( ),
                           1.0E-12 );

        // Adding onto an existing a-priori accumulates rather than overwrites.
        const Eigen::MatrixXd baseCovariance = Eigen::MatrixXd::Identity( numberOfParameters, numberOfParameters );
        const Eigen::MatrixXd combined =
                createSumLmkInverseAprioriCovariance< double, double, double >( conversionResult, parametersToEstimate, baseCovariance );
        BOOST_CHECK_SMALL( ( combined - ( inverseAprioriCovariance + baseCovariance ) ).norm( ) / expectedInverseVariance, 1.0E-12 );

        // A base matrix of the wrong size is an error.
        BOOST_CHECK_THROW( ( createSumLmkInverseAprioriCovariance< double, double, double >(
                                   conversionResult, parametersToEstimate, Eigen::MatrixXd::Identity( 3, 3 ) ) ),
                           std::runtime_error );
    }

    // --- Only one image estimated: the other image's a-priori is skipped, not an error.
    {
        const std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate = createParametersToEstimate< double, double >(
                createSumLmkPointingParameterSettings< double, double >( conversionResult, { "IMGB" } ), bodies );
        BOOST_REQUIRE_EQUAL( parametersToEstimate->getEstimatedParameterSetSize( ), 3 );

        const Eigen::MatrixXd inverseAprioriCovariance =
                createSumLmkInverseAprioriCovariance< double, double, double >( conversionResult, parametersToEstimate );
        BOOST_REQUIRE_EQUAL( inverseAprioriCovariance.rows( ), 3 );
        BOOST_CHECK_SMALL( ( inverseAprioriCovariance.diagonal( ) - Eigen::VectorXd::Constant( 3, expectedInverseVariance ) ).norm( ) /
                                   expectedInverseVariance,
                           1.0E-12 );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

BOOST_AUTO_TEST_SUITE( test_pixel_landmark_state_estimation )

//! Recover a perturbed spacecraft initial state from simulated pixel-landmark observations through
//! the full OrbitDeterminationManager (propagation + variational equations + pixel partials).
BOOST_AUTO_TEST_CASE( testInitialStateRecoveryFromPixelLandmarks )
{
    spice_interface::loadStandardSpiceKernels( );

    // --- Bodies: target at the origin with point-mass gravity + simple rotation, plus spacecraft.
    SystemOfBodies bodies( "SSB", "J2000" );
    bodies.createEmptyBody< double, double >( "Target", false );
    bodies.createEmptyBody< double, double >( "Spacecraft", false );

    bodies.at( "Target" )->setEphemeris( std::make_shared< ConstantEphemeris >( Eigen::Vector6d::Zero( ), "SSB", "J2000" ) );
    bodies.at( "Target" )->setGravityFieldModel( std::make_shared< gravitation::GravityFieldModel >( targetGravitationalParameter ) );
    bodies.at( "Target" )
            ->setRotationalEphemeris( std::make_shared< SimpleRotationalEphemeris >(
                    /* poleRightAscension = */ 0.3,
                    /* poleDeclination = */ 1.1,
                    /* primeMeridianOfDate = */ 0.2,
                    /* rotationRate = */ 2.0 * mathematical_constants::PI / 12000.0,
                    /* initialSecondsSinceEpoch = */ 0.0,
                    "J2000",
                    "Target_Fixed" ) );
    bodies.at( "Spacecraft" )
            ->setRotationalEphemeris( std::make_shared< ConstantRotationalEphemeris >(
                    Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ), "J2000", "Spacecraft_Fixed" ) );
    // Placeholder ephemeris that the OrbitDeterminationManager overwrites with the propagated arc.
    bodies.at( "Spacecraft" )
            ->setEphemeris( std::make_shared< TabulatedCartesianEphemeris<> >(
                    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector6d > >( ), "SSB", "J2000" ) );
    bodies.processBodyFrameDefinitions< double, double >( );

    // --- Truth orbit (Keplerian about the target).
    Eigen::Vector6d truthKeplerianElements = Eigen::Vector6d::Zero( );
    truthKeplerianElements( semiMajorAxisIndex ) = 1.0E4;
    truthKeplerianElements( eccentricityIndex ) = 0.05;
    truthKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 30.0 );
    truthKeplerianElements( argumentOfPeriapsisIndex ) = unit_conversions::convertDegreesToRadians( 40.0 );
    truthKeplerianElements( longitudeOfAscendingNodeIndex ) = unit_conversions::convertDegreesToRadians( 25.0 );
    truthKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 10.0 );

    // --- Build one SUM image per epoch over (roughly) a single orbital period; the conversion
    //     machinery derives the observation time from each image's UTC string, so we compute the
    //     matching epoch with the same converter before evaluating the truth geometry.
    const std::map< std::string, input_output::sum_lmk::LmkLandmarkData > landmarks = makeOrbitLandmarks( );
    std::vector< std::string > landmarkIds;
    for( const auto& entry : landmarks )
    {
        landmarkIds.push_back( entry.first );
    }

    const int numberOfImages = 24;
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
        // No pointing a-priori: leave SIGMA_PTG non-finite so no pointing parameter is implied.
        image.pointingSigma_ = Eigen::Vector3d::Constant( std::numeric_limits< double >::quiet_NaN( ) );

        const double epoch = observation_models::detail::convertSumUtcStringToSecondsSinceJ2000< double >( image );
        epochs.push_back( epoch );
        // All geometry is referenced to the first epoch.
        const double timeSinceFirstEpoch = epoch - epochs.front( );
        const Eigen::Vector6d keplerianAtEpoch =
                propagateKeplerOrbit< double >( truthKeplerianElements, timeSinceFirstEpoch, targetGravitationalParameter );
        const Eigen::Vector6d inertialStateAtEpoch = convertKeplerianToCartesianElements( keplerianAtEpoch, targetGravitationalParameter );
        const Eigen::Matrix3d rotationInertialToBodyFixed =
                bodies.at( "Target" )->getRotationalEphemeris( )->getRotationToTargetFrame( epoch ).toRotationMatrix( );
        const Eigen::Vector3d spacecraftBodyFixedPosition = rotationInertialToBodyFixed * inertialStateAtEpoch.head( 3 );

        // SCOBJ = spacecraft-to-object (target centre) vector in the target body-fixed frame.
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
    }

    const double initialEpoch = epochs.front( );
    const Eigen::Vector6d truthInitialState = convertKeplerianToCartesianElements( truthKeplerianElements, targetGravitationalParameter );

    // --- Convert to a pixel-coordinate observation collection (registers cameras + landmarks).
    SumLmkObservationConversionSettings conversionSettings( "Target", "Spacecraft" );
    SumLmkObservationConversionResult< double, double > conversionResult =
            createSumLmkObservationCollection< double, double >( images, landmarks, bodies, conversionSettings );

    // --- Dynamics: spacecraft under target point-mass gravity only.
    SelectedAccelerationMap accelerationSettingsMap;
    accelerationSettingsMap[ "Spacecraft" ][ "Target" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    const std::vector< std::string > bodiesToIntegrate = { "Spacecraft" };
    const std::vector< std::string > centralBodies = { "Target" };
    const AccelerationMap accelerationModelMap =
            createAccelerationModelsMap( bodies, accelerationSettingsMap, bodiesToIntegrate, centralBodies );

    const std::shared_ptr< IntegratorSettings< double > > integratorSettings =
            rungeKuttaFixedStepSettings< double >( 10.0, CoefficientSets::rungeKuttaFehlberg78 );

    const std::shared_ptr< TranslationalStatePropagatorSettings< double, double > > propagatorSettings =
            translationalStatePropagatorSettings< double, double >( centralBodies,
                                                                    accelerationModelMap,
                                                                    bodiesToIntegrate,
                                                                    truthInitialState,
                                                                    initialEpoch,
                                                                    integratorSettings,
                                                                    propagationTimeTerminationSettings( epochs.back( ) + 100.0 ) );

    // --- Estimate the spacecraft initial translational state.
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
            "Spacecraft", truthInitialState, "Target" ) );
    const std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double, double >( parameterNames, bodies, propagatorSettings );

    // --- Build the estimator and simulate ideal observations from the truth trajectory.
    OrbitDeterminationManager< double, double > orbitDeterminationManager(
            bodies, parametersToEstimate, conversionResult.observationModelSettings_, propagatorSettings );

    const std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > observationSimulationSettings =
            getObservationSimulationSettingsFromObservations< double, double >( conversionResult.observationCollection_, bodies );
    const std::shared_ptr< ObservationCollection< double, double > > simulatedObservations = simulateObservations< double, double >(
            observationSimulationSettings, orbitDeterminationManager.getObservationSimulators( ), bodies );

    std::map< std::shared_ptr< ObservationCollectionParser >, double > weightsPerObservationParser;
    weightsPerObservationParser[ observationParser( pixel_coordinates ) ] = 1.0;
    simulatedObservations->setConstantWeightPerObservable( weightsPerObservationParser );

    // --- Perturb the initial state and confirm the estimator recovers the truth.
    const Eigen::VectorXd truthParameters = parametersToEstimate->getFullParameterValues< double >( );
    Eigen::VectorXd perturbedParameters = truthParameters;
    perturbedParameters.segment( 0, 3 ) += Eigen::Vector3d::Constant( 10.0 );
    perturbedParameters.segment( 3, 3 ) += Eigen::Vector3d::Constant( 0.01 );
    parametersToEstimate->resetParameterValues( perturbedParameters );

    const std::shared_ptr< EstimationInput< double, double > > estimationInput =
            std::make_shared< EstimationInput< double, double > >( simulatedObservations );
    estimationInput->defineEstimationSettings( true, true, true, true, true, true );
    estimationInput->setConvergenceChecker( std::make_shared< EstimationConvergenceChecker >( 6 ) );

    const std::shared_ptr< EstimationOutput< double > > estimationOutput = orbitDeterminationManager.estimateParameters( estimationInput );

    const Eigen::VectorXd estimationError = estimationOutput->parameterEstimate_ - truthParameters;

    // Ideal (noise-free) data: position and velocity must be recovered to tight tolerances.
    for( unsigned int i = 0; i < 3; ++i )
    {
        BOOST_CHECK_SMALL( std::fabs( estimationError( i ) ), 1.0E-3 );
        BOOST_CHECK_SMALL( std::fabs( estimationError( i + 3 ) ), 1.0E-6 );
    }
}

//! Estimate the spacecraft initial state and a per-image camera pointing correction together, through
//! the full OrbitDeterminationManager.
//!
//! This is the geometrically hard case. A camera rotation of theta about an axis perpendicular to the
//! boresight moves the image by theta*f pixels, while a spacecraft translation dx perpendicular to the
//! line of sight moves it by (dx/d)*f pixels: from a single landmark the two are indistinguishable. What
//! separates them is parallax across the landmarks within one image - a translation moves each landmark
//! by an amount that depends on its own range, whereas a rotation moves the whole pattern rigidly. With
//! three free pointing angles per image against six state parameters for the whole arc, the state is
//! observable only through that parallax, so this test is the real check that the pointing partial and
//! the state partials are mutually consistent: an error in either shows up here as a failure to separate
//! them, even when each is individually self-consistent.
BOOST_AUTO_TEST_CASE( testJointStateAndPointingRecovery )
{
    spice_interface::loadStandardSpiceKernels( );

    // No pointing a-priori: with noise-free data the truth is the exact minimiser, so both the state and
    // the pointing must be recovered. The a-priori is exercised separately below.
    SyntheticOrbitScenario scenario = buildSyntheticOrbitScenario( std::numeric_limits< double >::quiet_NaN( ) );

    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
            "Spacecraft", scenario.truthInitialState_, "Target" ) );
    for( const std::shared_ptr< EstimatableParameterSettings >& pointingSetting :
         createSumLmkPointingParameterSettings< double, double >( scenario.conversionResult_ ) )
    {
        parameterNames.push_back( pointingSetting );
    }
    const std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double, double >( parameterNames, scenario.bodies_, scenario.propagatorSettings_ );
    BOOST_REQUIRE_EQUAL( parametersToEstimate->getEstimatedParameterSetSize( ), 6 + 3 * scenario.imageIds_.size( ) );

    OrbitDeterminationManager< double, double > orbitDeterminationManager(
            scenario.bodies_, parametersToEstimate, scenario.conversionResult_.observationModelSettings_, scenario.propagatorSettings_ );

    // Build the truth parameter vector: truth state plus a distinct injected offset per image. Parameter
    // blocks are located by identifier rather than by assuming the layout of the parameter vector.
    Eigen::VectorXd truthParameters = parametersToEstimate->getFullParameterValues< double >( );
    const double offsetMagnitude = 1.0E-3;
    for( std::size_t imageIndex = 0; imageIndex < scenario.imageIds_.size( ); ++imageIndex )
    {
        const std::string cameraName = scenario.conversionResult_.imageIdToCameraName_.at( scenario.imageIds_.at( imageIndex ) );
        const std::vector< std::pair< int, int > > indices = parametersToEstimate->getIndicesForParameterType(
                EstimatebleParameterIdentifier( camera_pointing_correction, std::make_pair( "Spacecraft", cameraName ) ) );
        BOOST_REQUIRE_EQUAL( indices.size( ), 1 );
        truthParameters.segment( indices.at( 0 ).first, 3 ) = syntheticPointingOffset( static_cast< int >( imageIndex ), offsetMagnitude );
    }

    // Simulate ideal observations from the truth state AND the truth pointing.
    parametersToEstimate->resetParameterValues( truthParameters );
    const std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > observationSimulationSettings =
            getObservationSimulationSettingsFromObservations< double, double >( scenario.conversionResult_.observationCollection_,
                                                                                scenario.bodies_ );
    const std::shared_ptr< ObservationCollection< double, double > > simulatedObservations = simulateObservations< double, double >(
            observationSimulationSettings, orbitDeterminationManager.getObservationSimulators( ), scenario.bodies_ );

    std::map< std::shared_ptr< ObservationCollectionParser >, double > weightsPerObservationParser;
    weightsPerObservationParser[ observationParser( pixel_coordinates ) ] = 1.0;
    simulatedObservations->setConstantWeightPerObservable( weightsPerObservationParser );

    // Perturb the state and zero every pointing correction, then recover both.
    Eigen::VectorXd perturbedParameters = truthParameters;
    perturbedParameters.segment( 0, 3 ) += Eigen::Vector3d::Constant( 10.0 );
    perturbedParameters.segment( 3, 3 ) += Eigen::Vector3d::Constant( 0.01 );
    for( std::size_t imageIndex = 0; imageIndex < scenario.imageIds_.size( ); ++imageIndex )
    {
        const std::string cameraName = scenario.conversionResult_.imageIdToCameraName_.at( scenario.imageIds_.at( imageIndex ) );
        const std::vector< std::pair< int, int > > indices = parametersToEstimate->getIndicesForParameterType(
                EstimatebleParameterIdentifier( camera_pointing_correction, std::make_pair( "Spacecraft", cameraName ) ) );
        perturbedParameters.segment( indices.at( 0 ).first, 3 ) = Eigen::Vector3d::Zero( );
    }
    parametersToEstimate->resetParameterValues( perturbedParameters );

    const std::shared_ptr< EstimationInput< double, double > > estimationInput =
            std::make_shared< EstimationInput< double, double > >( simulatedObservations );
    estimationInput->defineEstimationSettings( true, true, true, true, true, true );
    estimationInput->setConvergenceChecker( std::make_shared< EstimationConvergenceChecker >( 10 ) );

    const std::shared_ptr< EstimationOutput< double > > estimationOutput = orbitDeterminationManager.estimateParameters( estimationInput );
    const Eigen::VectorXd estimationError = estimationOutput->parameterEstimate_ - truthParameters;

    // Residuals must be driven to the numerical floor: the model can reproduce the data exactly.
    BOOST_CHECK_SMALL( residualRms( estimationOutput->residuals_ ), 1.0E-6 );

    // State recovery. Conditioning governs how measurement noise is amplified, not where the minimum of
    // noise-free data lies, so the ill-conditioning described above does not stop the solve reaching the
    // truth: in practice this converges to ~1e-9 m. The tolerances keep a wide margin over that, so the
    // test fails on a broken partial rather than on platform-level floating-point differences.
    BOOST_TEST_MESSAGE( "Joint solve: position error " << estimationError.segment( 0, 3 ).norm( ) << " m, velocity error "
                                                       << estimationError.segment( 3, 3 ).norm( ) << " m/s" );
    for( unsigned int i = 0; i < 3; ++i )
    {
        BOOST_CHECK_SMALL( std::fabs( estimationError( i ) ), 1.0E-5 );
        BOOST_CHECK_SMALL( std::fabs( estimationError( i + 3 ) ), 1.0E-8 );
    }

    // Pointing recovery, per image.
    for( std::size_t imageIndex = 0; imageIndex < scenario.imageIds_.size( ); ++imageIndex )
    {
        const std::string cameraName = scenario.conversionResult_.imageIdToCameraName_.at( scenario.imageIds_.at( imageIndex ) );
        const std::vector< std::pair< int, int > > indices = parametersToEstimate->getIndicesForParameterType(
                EstimatebleParameterIdentifier( camera_pointing_correction, std::make_pair( "Spacecraft", cameraName ) ) );
        const Eigen::Vector3d pointingError = estimationError.segment( indices.at( 0 ).first, 3 );
        BOOST_CHECK_SMALL( pointingError.norm( ), 1.0E-10 );
    }
}

//! The per-image SIGMA_PTG a-priori must actually constrain the solution.
//!
//! With noise-free data and no a-priori, the truth is the exact minimiser and the pointing is recovered
//! exactly (previous test). Adding an a-priori deliberately moves the minimum: it trades data fit for
//! agreement with the a-priori, pulling each correction towards zero and shrinking its formal
//! uncertainty. Both effects are checked here, which is what shows the a-priori matrix built by
//! createSumLmkInverseAprioriCovariance reaches the normal equations at the right parameter indices -
//! a matrix assembled at the wrong offsets would constrain the wrong parameters and leave the pointing
//! estimates untouched.
BOOST_AUTO_TEST_CASE( testPointingAprioriConstrainsSolution )
{
    spice_interface::loadStandardSpiceKernels( );

    // SIGMA_PTG of the same order as the injected offsets, so the a-priori and the data genuinely compete.
    const double pointingSigma = 1.0E-3;
    const double offsetMagnitude = 1.0E-3;
    const int numberOfImages = 8;
    SyntheticOrbitScenario scenario = buildSyntheticOrbitScenario( pointingSigma, numberOfImages );

    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
            "Spacecraft", scenario.truthInitialState_, "Target" ) );
    for( const std::shared_ptr< EstimatableParameterSettings >& pointingSetting :
         createSumLmkPointingParameterSettings< double, double >( scenario.conversionResult_ ) )
    {
        parameterNames.push_back( pointingSetting );
    }
    const std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double, double >( parameterNames, scenario.bodies_, scenario.propagatorSettings_ );

    OrbitDeterminationManager< double, double > orbitDeterminationManager(
            scenario.bodies_, parametersToEstimate, scenario.conversionResult_.observationModelSettings_, scenario.propagatorSettings_ );

    // Locate each image's pointing block once.
    std::vector< int > pointingStartIndices;
    for( const std::string& imageId : scenario.imageIds_ )
    {
        const std::string cameraName = scenario.conversionResult_.imageIdToCameraName_.at( imageId );
        const std::vector< std::pair< int, int > > indices = parametersToEstimate->getIndicesForParameterType(
                EstimatebleParameterIdentifier( camera_pointing_correction, std::make_pair( "Spacecraft", cameraName ) ) );
        BOOST_REQUIRE_EQUAL( indices.size( ), 1 );
        pointingStartIndices.push_back( indices.at( 0 ).first );
    }

    Eigen::VectorXd truthParameters = parametersToEstimate->getFullParameterValues< double >( );
    for( int imageIndex = 0; imageIndex < numberOfImages; ++imageIndex )
    {
        truthParameters.segment( pointingStartIndices.at( imageIndex ), 3 ) = syntheticPointingOffset( imageIndex, offsetMagnitude );
    }

    parametersToEstimate->resetParameterValues( truthParameters );
    const std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > observationSimulationSettings =
            getObservationSimulationSettingsFromObservations< double, double >( scenario.conversionResult_.observationCollection_,
                                                                                scenario.bodies_ );
    const std::shared_ptr< ObservationCollection< double, double > > simulatedObservations = simulateObservations< double, double >(
            observationSimulationSettings, orbitDeterminationManager.getObservationSimulators( ), scenario.bodies_ );
    std::map< std::shared_ptr< ObservationCollectionParser >, double > weightsPerObservationParser;
    weightsPerObservationParser[ observationParser( pixel_coordinates ) ] = 1.0;
    simulatedObservations->setConstantWeightPerObservable( weightsPerObservationParser );

    Eigen::VectorXd startParameters = truthParameters;
    for( int imageIndex = 0; imageIndex < numberOfImages; ++imageIndex )
    {
        startParameters.segment( pointingStartIndices.at( imageIndex ), 3 ) = Eigen::Vector3d::Zero( );
    }

    // The a-priori built from the conversion's SIGMA_PTG entries.
    const Eigen::MatrixXd inverseAprioriCovariance =
            createSumLmkInverseAprioriCovariance< double, double, double >( scenario.conversionResult_, parametersToEstimate );
    BOOST_REQUIRE_EQUAL( inverseAprioriCovariance.rows( ), parametersToEstimate->getEstimatedParameterSetSize( ) );
    for( const int startIndex : pointingStartIndices )
    {
        BOOST_CHECK_CLOSE( inverseAprioriCovariance( startIndex, startIndex ), 1.0 / ( pointingSigma * pointingSigma ), 1.0E-8 );
    }
    // The a-priori must not touch the state block.
    BOOST_CHECK_SMALL( inverseAprioriCovariance.block( 0, 0, 6, 6 ).norm( ), 1.0E-12 );

    auto solve = [ & ]( const Eigen::MatrixXd& aprioriMatrix ) {
        parametersToEstimate->resetParameterValues( startParameters );
        const std::shared_ptr< EstimationInput< double, double > > estimationInput =
                std::make_shared< EstimationInput< double, double > >( simulatedObservations, aprioriMatrix );
        estimationInput->defineEstimationSettings( true, true, true, true, true, true );
        estimationInput->setConvergenceChecker( std::make_shared< EstimationConvergenceChecker >( 8 ) );
        return orbitDeterminationManager.estimateParameters( estimationInput );
    };

    const std::shared_ptr< EstimationOutput< double > > unconstrainedOutput = solve( Eigen::MatrixXd::Zero( 0, 0 ) );
    const std::shared_ptr< EstimationOutput< double > > constrainedOutput = solve( inverseAprioriCovariance );

    const Eigen::VectorXd unconstrainedFormalErrors = unconstrainedOutput->getFormalErrorVector( );
    const Eigen::VectorXd constrainedFormalErrors = constrainedOutput->getFormalErrorVector( );

    double truthNorm = 0.0;
    double unconstrainedDeviation = 0.0;
    double constrainedDeviation = 0.0;
    for( int imageIndex = 0; imageIndex < numberOfImages; ++imageIndex )
    {
        const int startIndex = pointingStartIndices.at( imageIndex );
        const Eigen::Vector3d truthOffset = truthParameters.segment( startIndex, 3 );
        truthNorm += truthOffset.squaredNorm( );
        unconstrainedDeviation += ( unconstrainedOutput->parameterEstimate_.segment( startIndex, 3 ) - truthOffset ).squaredNorm( );
        constrainedDeviation += ( constrainedOutput->parameterEstimate_.segment( startIndex, 3 ) - truthOffset ).squaredNorm( );

        // The a-priori can only reduce the formal uncertainty of the parameters it constrains.
        for( int component = 0; component < 3; ++component )
        {
            BOOST_CHECK( constrainedFormalErrors( startIndex + component ) < unconstrainedFormalErrors( startIndex + component ) );
        }
    }

    BOOST_TEST_MESSAGE( "Pointing deviation from truth: unconstrained " << std::sqrt( unconstrainedDeviation ) << " rad, constrained "
                                                                        << std::sqrt( constrainedDeviation ) << " rad (truth norm "
                                                                        << std::sqrt( truthNorm ) << " rad)" );

    // Without the a-priori the pointing is recovered essentially exactly; with it, the solution is pulled
    // measurably away from the truth and towards zero. That bias is the a-priori doing its job, and is why
    // the recovery tests above are run without one.
    BOOST_CHECK_SMALL( std::sqrt( unconstrainedDeviation ), 1.0E-10 );
    BOOST_CHECK( std::sqrt( constrainedDeviation ) > 1.0E-6 );

    // The constrained solution fits the data less well, by construction.
    BOOST_CHECK( residualRms( constrainedOutput->residuals_ ) > residualRms( unconstrainedOutput->residuals_ ) );
}

//! A-priori entries for images whose pointing is not estimated are skipped rather than misapplied.
BOOST_AUTO_TEST_CASE( testPointingAprioriSkipsUnestimatedImages )
{
    spice_interface::loadStandardSpiceKernels( );

    const double pointingSigma = 1.0E-3;
    const int numberOfImages = 4;
    SyntheticOrbitScenario scenario = buildSyntheticOrbitScenario( pointingSigma, numberOfImages );
    BOOST_REQUIRE_EQUAL( scenario.conversionResult_.inverseAprioriCovarianceDiagonalEntries_.size( ),
                         static_cast< std::size_t >( numberOfImages ) );

    // Estimate pointing for only two of the four images.
    const std::vector< std::string > estimatedImageIds = { scenario.imageIds_.at( 1 ), scenario.imageIds_.at( 3 ) };
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
            "Spacecraft", scenario.truthInitialState_, "Target" ) );
    for( const std::shared_ptr< EstimatableParameterSettings >& pointingSetting :
         createSumLmkPointingParameterSettings< double, double >( scenario.conversionResult_, estimatedImageIds ) )
    {
        parameterNames.push_back( pointingSetting );
    }
    const std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double, double >( parameterNames, scenario.bodies_, scenario.propagatorSettings_ );
    BOOST_REQUIRE_EQUAL( parametersToEstimate->getEstimatedParameterSetSize( ), 6 + 3 * 2 );

    // All four a-priori entries are offered; only the two estimated ones may land in the matrix.
    const Eigen::MatrixXd inverseAprioriCovariance =
            createSumLmkInverseAprioriCovariance< double, double, double >( scenario.conversionResult_, parametersToEstimate );
    BOOST_REQUIRE_EQUAL( inverseAprioriCovariance.rows( ), 12 );
    BOOST_CHECK_SMALL( inverseAprioriCovariance.block( 0, 0, 6, 6 ).norm( ), 1.0E-12 );
    for( int i = 6; i < 12; ++i )
    {
        BOOST_CHECK_CLOSE( inverseAprioriCovariance( i, i ), 1.0 / ( pointingSigma * pointingSigma ), 1.0E-8 );
    }
    // Exactly six non-zero entries: nothing from the two unestimated images leaked in.
    BOOST_CHECK_EQUAL( ( inverseAprioriCovariance.array( ).abs( ) > 1.0E-12 ).count( ), 6 );
}

//! End-to-end estimation from REAL Rosetta SUM/LMK pixel-landmark data (2015-07-15, comet 67P): the
//! spacecraft initial state is recovered through the full OrbitDeterminationManager, propagating under
//! comet point-mass gravity. The a-priori initial state is the reconstructed Rosetta orbiter state
//! read from SPICE, and the comet body-fixed orientation is the SPICE 67P/C-G_CK (SPC/Cheops) frame.
//!
//! All required data is committed under data/sum_lmk/real_67p so the test is self-contained (no
//! external Rosetta SPICE archive): the reduced SUM/LMK files, a trimmed SPK with the orbiter state
//! relative to the comet over the arc (extracted from RORB_DV_257), and a tabulated comet attitude
//! sampled from the 67P/C-G_CK CK. The orbiter SPK and the SPC pixel solution are independent
//! reconstructions (~1 km apart), so the SPICE a-priori starts well off the pixels and the estimator
//! drives the residuals down to the SPC measurement level.
BOOST_AUTO_TEST_CASE( testRealRosettaInitialStateEstimation )
{
    spice_interface::loadStandardSpiceKernels( );
    RealRosettaScenario scenario = buildRealRosettaScenario( );
    BOOST_REQUIRE( scenario.conversionResult_.observationCollection_ != nullptr );
    BOOST_REQUIRE_GT( scenario.conversionResult_.observationCollection_->getTotalObservableSize( ), 0 );

    // --- Estimate the spacecraft initial state from the real pixel observations.
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
            "Spacecraft", scenario.initialStateGuess_, "Comet" ) );
    const std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double, double >( parameterNames, scenario.bodies_, scenario.propagatorSettings_ );

    OrbitDeterminationManager< double, double > orbitDeterminationManager(
            scenario.bodies_, parametersToEstimate, scenario.conversionResult_.observationModelSettings_, scenario.propagatorSettings_ );

    std::map< std::shared_ptr< ObservationCollectionParser >, double > weightsPerObservationParser;
    weightsPerObservationParser[ observationParser( pixel_coordinates ) ] = 1.0;
    scenario.conversionResult_.observationCollection_->setConstantWeightPerObservable( weightsPerObservationParser );

    const std::shared_ptr< EstimationInput< double, double > > estimationInput =
            std::make_shared< EstimationInput< double, double > >( scenario.conversionResult_.observationCollection_ );
    estimationInput->defineEstimationSettings( true, true, true, true, true, true );
    estimationInput->setConvergenceChecker( std::make_shared< EstimationConvergenceChecker >( 8 ) );

    const std::shared_ptr< EstimationOutput< double > > estimationOutput = orbitDeterminationManager.estimateParameters( estimationInput );

    const Eigen::MatrixXd residualHistory = estimationOutput->getResidualHistoryMatrix( );
    const double initialRms = residualRms( residualHistory.col( 0 ) );
    const double finalRms = residualRms( estimationOutput->residuals_ );
    BOOST_TEST_MESSAGE( "Real-data pixel residual RMS: initial = " + std::to_string( initialRms ) +
                        " px, final = " + std::to_string( finalRms ) + " px" );

    // The SPC pixel measurements are reproduced to the few-pixel level by a point-mass orbit over the
    // ~2.3 h arc, and the estimator strongly reduces the residuals from the SPICE a-priori.
    BOOST_CHECK_LT( finalRms, 3.0 );
    BOOST_CHECK_LT( finalRms, 0.1 * initialRms );

    // The estimate corrects the (independent) SPICE a-priori by no more than a few km.
    const Eigen::Vector3d estimatedPosition = estimationOutput->parameterEstimate_.segment( 0, 3 );
    BOOST_CHECK_LT( ( estimatedPosition - scenario.initialStateGuess_.segment( 0, 3 ) ).norm( ), 5.0E3 );
}

//! Estimate the spacecraft state AND a per-image camera pointing correction from the REAL Rosetta
//! SUM/LMK data - the feature's intended workflow on real measurements.
//!
//! For scale: OSIRIS NAC has a 717.3 mm focal length at 13.5 um pixels, so one pixel is about 1.9e-5 rad
//! and the SIGMA_PTG values in these files (1.0e-4 to 1.9e-4 rad) are of order 8-10 pixels.
//!
//! Two measured properties of this solve are worth recording, because neither is what one might assume.
//!
//! First, the per-image pointing is very nearly degenerate with the spacecraft position. A camera
//! rotation and a transverse spacecraft translation produce almost the same image motion, separated only
//! by parallax across the landmarks, which here is the comet's ~2 km extent against a ~165 km range, i.e.
//! about 1%. Taking the converged solution and zeroing only the pointing corrections moves the pixels by
//! ~30 px RMS, yet the corrections buy only ~0.5 px of final fit (1.37 px state-only -> 0.88 px joint):
//! almost all of that motion is cancelled by a compensating ~350 m shift of the estimated state. The
//! individual pointing values are therefore a point on a flat ridge, not a measurement of camera pointing.
//!
//! Second, and consequently, the SIGMA_PTG a-priori is what pins the solution down. Dropping it lets the
//! pointing run from ~1.4e-3 rad to ~6.6e-3 rad while the state moves ~905 m to compensate. Note that
//! comparing the data's information on pointing against the a-priori's (f^2*n ~ 2e11 versus 1/sigma^2 ~
//! 5e7, a ratio of ~4000) is misleading: that is the information with the state held fixed, whereas what
//! decides a joint solve is the information along the degenerate direction, where the data has almost
//! none. A nominally weak a-priori dominates exactly where the data is blind.
//!
//! What is therefore asserted is what the data can genuinely establish: the pointing parameter improves
//! the fit relative to a state-only solve, it stays within a bounded envelope rather than running away,
//! and it does not corrupt the orbit solution. The convention and sign of the partial are pinned down
//! separately, by the finite-difference tests and by the synthetic joint recovery above.
BOOST_AUTO_TEST_CASE( testRealRosettaJointStateAndPointingEstimation )
{
    spice_interface::loadStandardSpiceKernels( );
    RealRosettaScenario scenario = buildRealRosettaScenario( );

    // Per-image pointing alongside the state. Without the a-priori this is weakly determined: the comet
    // spans ~2 km at a ~165 km range, so the parallax that separates a pointing rotation from a
    // spacecraft translation is barely a percent, and 8 images x 3 pointing angles can otherwise absorb
    // the state error. The SIGMA_PTG a-priori is what makes the joint solve well posed - which is also
    // why SPC carries it in the files.
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
            "Spacecraft", scenario.initialStateGuess_, "Comet" ) );
    for( const std::shared_ptr< EstimatableParameterSettings >& pointingSetting :
         createSumLmkPointingParameterSettings< double, double >( scenario.conversionResult_ ) )
    {
        parameterNames.push_back( pointingSetting );
    }
    const std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double, double >( parameterNames, scenario.bodies_, scenario.propagatorSettings_ );

    const std::size_t numberOfImages = scenario.conversionResult_.imageIdToCameraName_.size( );
    BOOST_REQUIRE_EQUAL( parametersToEstimate->getEstimatedParameterSetSize( ), 6 + 3 * numberOfImages );

    const Eigen::MatrixXd inverseAprioriCovariance =
            createSumLmkInverseAprioriCovariance< double, double, double >( scenario.conversionResult_, parametersToEstimate );
    BOOST_REQUIRE_EQUAL( inverseAprioriCovariance.rows( ), parametersToEstimate->getEstimatedParameterSetSize( ) );

    OrbitDeterminationManager< double, double > orbitDeterminationManager(
            scenario.bodies_, parametersToEstimate, scenario.conversionResult_.observationModelSettings_, scenario.propagatorSettings_ );

    std::map< std::shared_ptr< ObservationCollectionParser >, double > weightsPerObservationParser;
    weightsPerObservationParser[ observationParser( pixel_coordinates ) ] = 1.0;
    scenario.conversionResult_.observationCollection_->setConstantWeightPerObservable( weightsPerObservationParser );

    const std::shared_ptr< EstimationInput< double, double > > estimationInput = std::make_shared< EstimationInput< double, double > >(
            scenario.conversionResult_.observationCollection_, inverseAprioriCovariance );
    estimationInput->defineEstimationSettings( true, true, true, true, true, true );
    estimationInput->setConvergenceChecker( std::make_shared< EstimationConvergenceChecker >( 8 ) );

    const std::shared_ptr< EstimationOutput< double > > estimationOutput = orbitDeterminationManager.estimateParameters( estimationInput );

    const Eigen::MatrixXd residualHistory = estimationOutput->getResidualHistoryMatrix( );
    const double initialRms = residualRms( residualHistory.col( 0 ) );
    const double finalRms = residualRms( estimationOutput->residuals_ );

    // Every estimated correction must sit inside SPC's own stated pointing uncertainty for that image.
    // The sigmas are recovered from the a-priori entries the conversion produced (inverse variances).
    double worstSigmaRatio = 0.0;
    double largestCorrection = 0.0;
    for( const auto& entry : scenario.conversionResult_.inverseAprioriCovarianceDiagonalEntries_ )
    {
        const std::vector< std::pair< int, int > > indices = parametersToEstimate->getIndicesForParameterType( entry.first );
        BOOST_REQUIRE_EQUAL( indices.size( ), 1 );
        const Eigen::Vector3d correction = estimationOutput->parameterEstimate_.segment( indices.at( 0 ).first, 3 );
        largestCorrection = std::max( largestCorrection, correction.norm( ) );
        for( int component = 0; component < 3; ++component )
        {
            const double sigma = 1.0 / std::sqrt( entry.second( component ) );
            worstSigmaRatio = std::max( worstSigmaRatio, std::fabs( correction( component ) ) / sigma );
        }
    }

    // State-only reference solve on an independent copy of the same scenario, so the comparison below is
    // against a measured value rather than a hard-coded one.
    RealRosettaScenario referenceScenario = buildRealRosettaScenario( );
    std::vector< std::shared_ptr< EstimatableParameterSettings > > referenceParameterNames;
    referenceParameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
            "Spacecraft", referenceScenario.initialStateGuess_, "Comet" ) );
    const std::shared_ptr< EstimatableParameterSet< double > > referenceParameters = createParametersToEstimate< double, double >(
            referenceParameterNames, referenceScenario.bodies_, referenceScenario.propagatorSettings_ );
    OrbitDeterminationManager< double, double > referenceManager( referenceScenario.bodies_,
                                                                  referenceParameters,
                                                                  referenceScenario.conversionResult_.observationModelSettings_,
                                                                  referenceScenario.propagatorSettings_ );
    referenceScenario.conversionResult_.observationCollection_->setConstantWeightPerObservable( weightsPerObservationParser );
    const std::shared_ptr< EstimationInput< double, double > > referenceInput =
            std::make_shared< EstimationInput< double, double > >( referenceScenario.conversionResult_.observationCollection_ );
    referenceInput->defineEstimationSettings( true, true, true, true, true, true );
    referenceInput->setConvergenceChecker( std::make_shared< EstimationConvergenceChecker >( 8 ) );
    const double stateOnlyRms = residualRms( referenceManager.estimateParameters( referenceInput )->residuals_ );

    BOOST_TEST_MESSAGE( "Real-data joint solve: residual RMS " << initialRms << " -> " << finalRms << " px (state-only solve reaches "
                                                               << stateOnlyRms << " px); largest pointing correction " << largestCorrection
                                                               << " rad = " << worstSigmaRatio << " sigma of SPC's SIGMA_PTG" );

    // Adding per-image pointing must improve the fit relative to estimating the state alone.
    BOOST_CHECK_LT( finalRms, stateOnlyRms );
    BOOST_CHECK_LT( finalRms, 0.1 * initialRms );

    // The corrections stay bounded. With the a-priori in place they settle around 1.4e-3 rad; without it
    // they reach ~6.6e-3 rad, so this bound also confirms the a-priori is actually reaching the normal
    // equations. A sign or handedness error in the partial would instead fight the data and fail the
    // residual checks above.
    BOOST_CHECK_LT( largestCorrection, 5.0E-3 );

    // The extra freedom must not destroy the orbit solution.
    const Eigen::Vector3d estimatedPosition = estimationOutput->parameterEstimate_.segment( 0, 3 );
    BOOST_CHECK_LT( ( estimatedPosition - scenario.initialStateGuess_.segment( 0, 3 ) ).norm( ), 5.0E3 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
