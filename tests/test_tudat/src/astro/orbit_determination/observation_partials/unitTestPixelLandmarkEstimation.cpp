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
//                                           parameter and DirectObservationPartial objects on
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

#include <Eigen/Cholesky>

#include "tudat/astro/basic_astro/keplerPropagator.h"
#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/ephemerides/constantEphemeris.h"
#include "tudat/astro/ephemerides/constantRotationalEphemeris.h"
#include "tudat/astro/ephemerides/simpleRotationalEphemeris.h"
#include "tudat/astro/ephemerides/multiArcEphemeris.h"
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

//! A deterministic, non-identity local frame for a synthetic landmark. UX/UY/UZ are the COLUMNS of the
//! rotation from the landmark local frame to the body-fixed frame, with UZ the maplet plane normal, per
//! the SPC convention (https://web.psi.edu/spc_wiki/MAPFILES). A non-identity frame is what makes
//! local-frame uncertainty reporting distinguishable from the body-fixed covariance.
void setSyntheticLandmarkLocalFrame( input_output::sum_lmk::LmkLandmarkData& landmark, const int index )
{
    const Eigen::Matrix3d localFrame =
            Eigen::AngleAxisd( 0.4 + 0.3 * index, ( Eigen::Vector3d( ) << 0.2, 0.7, -0.5 ).finished( ).normalized( ) ).toRotationMatrix( );
    landmark.localXAxis_ = localFrame.col( 0 );
    landmark.localYAxis_ = localFrame.col( 1 );
    landmark.localZAxis_ = localFrame.col( 2 );
}

//! Body-fixed landmark positions on the target [m], spread within a small radius so the full set
//! stays in the field of view across the synthetic orbit arc while still spanning three dimensions.
//! Pass a finite landmarkSigmaMetres to attach a SIGMA_LMK record, which makes the conversion emit a
//! landmark a-priori; leave it NaN for the tests that must see no landmark a-priori at all.
std::map< std::string, input_output::sum_lmk::LmkLandmarkData > makeOrbitLandmarks(
        const double landmarkSigmaMetres = std::numeric_limits< double >::quiet_NaN( ) )
{
    const std::map< std::string, Eigen::Vector3d > positions = { { "LMK01", ( Eigen::Vector3d( ) << 700.0, 200.0, 150.0 ).finished( ) },
                                                                 { "LMK02", ( Eigen::Vector3d( ) << -600.0, 300.0, -100.0 ).finished( ) },
                                                                 { "LMK03", ( Eigen::Vector3d( ) << 200.0, -650.0, 250.0 ).finished( ) },
                                                                 { "LMK04", ( Eigen::Vector3d( ) << -300.0, -400.0, -200.0 ).finished( ) },
                                                                 { "LMK05", ( Eigen::Vector3d( ) << 500.0, 500.0, 300.0 ).finished( ) },
                                                                 { "LMK06", ( Eigen::Vector3d( ) << -450.0, 150.0, 400.0 ).finished( ) } };
    std::map< std::string, input_output::sum_lmk::LmkLandmarkData > landmarks;
    int landmarkIndex = 0;
    for( const auto& entry : positions )
    {
        input_output::sum_lmk::LmkLandmarkData landmark;
        landmark.landmarkId_ = entry.first;
        landmark.bodyFixedPosition_ = entry.second;
        setSyntheticLandmarkLocalFrame( landmark, landmarkIndex );
        if( std::isfinite( landmarkSigmaMetres ) )
        {
            landmark.landmarkPositionSigma_ = Eigen::Vector3d::Constant( landmarkSigmaMetres );
        }
        landmarks[ entry.first ] = landmark;
        ++landmarkIndex;
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
    LoadedPixelTestData conversionResult_;
    std::shared_ptr< TranslationalStatePropagatorSettings< double, double > > propagatorSettings_;
    Eigen::Vector6d truthInitialState_;
    std::vector< std::string > imageIds_;
};

SyntheticOrbitScenario buildSyntheticOrbitScenario( const double pointingSigmaRadians,
                                                    const int numberOfImages = 24,
                                                    const double landmarkSigmaMetres = std::numeric_limits< double >::quiet_NaN( ) )
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

    scenario.landmarks_ = makeOrbitLandmarks( landmarkSigmaMetres );
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

//! An observation-only scenario for landmark position estimation: the spacecraft follows a known
//! Keplerian arc supplied as a tabulated ephemeris, so nothing is propagated and the trajectory is not
//! estimated, while the target rotates and images are taken along the arc.
//!
//! The varying geometry is the point. From a single fixed viewpoint a landmark only ever yields two
//! numbers no matter how many images are taken, so its depth would stay unobservable; it is the motion of
//! the observer and the rotation of the body that make all three body-fixed components estimable. The
//! semi-major axis is 10 km, the close-orbit regime where landmark positions carry real information.
struct LandmarkObservationScenario {
    SystemOfBodies bodies_;
    std::map< std::string, input_output::sum_lmk::LmkLandmarkData > landmarks_;
    LoadedPixelTestData conversionResult_;
    std::vector< std::string > imageIds_;
    std::vector< double > epochs_;
};

LandmarkObservationScenario buildLandmarkObservationScenario(
        const int numberOfImages = 12,
        const double landmarkSigmaMetres = std::numeric_limits< double >::quiet_NaN( ),
        const int epochSpacingSeconds = 500 )
{
    LandmarkObservationScenario scenario;

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

    Eigen::Vector6d keplerianElements = Eigen::Vector6d::Zero( );
    keplerianElements( semiMajorAxisIndex ) = 1.0E4;
    keplerianElements( eccentricityIndex ) = 0.05;
    keplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 30.0 );
    keplerianElements( argumentOfPeriapsisIndex ) = unit_conversions::convertDegreesToRadians( 40.0 );
    keplerianElements( longitudeOfAscendingNodeIndex ) = unit_conversions::convertDegreesToRadians( 25.0 );
    keplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 10.0 );

    // Image epochs first: the Keplerian arc is defined relative to the first of them.
    for( int imageIndex = 0; imageIndex < numberOfImages; ++imageIndex )
    {
        input_output::sum_lmk::SumImageData epochOnlyImage;
        epochOnlyImage.utcEpochString_ = makeUtcString( imageIndex * epochSpacingSeconds );
        scenario.epochs_.push_back( sumImageTdbTimeForTest( epochOnlyImage ) );
    }

    // Tabulated spacecraft ephemeris sampled from the Keplerian arc. The same object is used to simulate
    // the observations and to evaluate the model during estimation, so interpolation error is common to
    // both and cannot stop the residuals reaching the numerical floor.
    const double sampleStep = 5.0;
    const double firstEpoch = scenario.epochs_.front( );
    const double lastEpoch = scenario.epochs_.back( );
    std::map< double, Eigen::Vector6d > spacecraftStates;
    for( double epoch = firstEpoch - 100.0; epoch <= lastEpoch + 100.0 + 0.5 * sampleStep; epoch += sampleStep )
    {
        const Eigen::Vector6d keplerianAtEpoch =
                propagateKeplerOrbit< double >( keplerianElements, epoch - firstEpoch, targetGravitationalParameter );
        spacecraftStates[ epoch ] = convertKeplerianToCartesianElements( keplerianAtEpoch, targetGravitationalParameter );
    }
    bodies.at( "Spacecraft" )
            ->setEphemeris( std::make_shared< TabulatedCartesianEphemeris<> >(
                    std::make_shared< interpolators::CubicSplineInterpolator< double, Eigen::Vector6d > >( spacecraftStates ),
                    "SSB",
                    "J2000" ) );
    bodies.processBodyFrameDefinitions< double, double >( );

    scenario.landmarks_ = makeOrbitLandmarks( landmarkSigmaMetres );
    std::vector< std::string > landmarkIds;
    for( const auto& entry : scenario.landmarks_ )
    {
        landmarkIds.push_back( entry.first );
    }

    std::vector< input_output::sum_lmk::SumImageData > images;
    for( int imageIndex = 0; imageIndex < numberOfImages; ++imageIndex )
    {
        const double epoch = scenario.epochs_.at( imageIndex );

        input_output::sum_lmk::SumImageData image;
        image.imageId_ = "LIMG" + std::to_string( imageIndex );
        image.utcEpochString_ = makeUtcString( imageIndex * epochSpacingSeconds );
        image.imageSize_ = Eigen::Vector2i( 1024, 1024 );
        image.focalLengthMm_ = 100.0;
        image.opticalCenter_ = Eigen::Vector2d( 512.0, 512.0 );
        image.kMatrix_ << 10.0, 0.0, 0.0, 0.0, 10.0, 0.0;

        const Eigen::Vector3d inertialPosition =
                bodies.at( "Spacecraft" )->getStateInBaseFrameFromEphemeris< double, double >( epoch ).head( 3 );
        const Eigen::Vector3d spacecraftBodyFixedPosition =
                bodies.at( "Target" )->getRotationalEphemeris( )->getRotationToTargetFrame( epoch ).toRotationMatrix( ) * inertialPosition;
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

    scenario.conversionResult_ = loadPixelTestData( images, scenario.landmarks_, bodies );
    scenario.bodies_ = bodies;
    return scenario;
}

//! Thin non-template wrappers around the landmark helpers.
//!
//! Boost's check macros are preprocessor macros, so a comma inside template arguments at the call site is
//! parsed as an argument separator. These wrappers keep the call sites comma-free.
std::vector< std::shared_ptr< EstimatableParameterSettings > > landmarkParameterSettings(
        const LoadedPixelTestData& loaded,
        const std::shared_ptr< ObservationCollection< double, double > >& observationCollection = nullptr,
        const std::vector< std::string >& landmarkIds = std::vector< std::string >( ),
        const std::string& targetBody = "Target" )
{
    return getObservedReferencePointPositionParameterSettings< double, double >(
            observationCollection != nullptr ? observationCollection : loaded.observationCollection_, targetBody, landmarkIds );
}

//! IDs of the landmarks that still have observations, ordered by ID.
std::vector< std::string > observedLandmarkIds( const LoadedPixelTestData& loaded, const std::string& targetBody = "Target" )
{
    std::vector< std::string > landmarkIds;
    for( const std::shared_ptr< EstimatableParameterSettings >& settings : landmarkParameterSettings( loaded, nullptr, {}, targetBody ) )
    {
        landmarkIds.push_back( settings->parameterType_.second.second );
    }
    return landmarkIds;
}

//! Local-frame uncertainties keyed by landmark ID only.
std::map< std::string, Eigen::Vector3d > landmarkLocalFrameUncertainties(
        const LoadedPixelTestData& loaded,
        const std::shared_ptr< EstimatableParameterSet< double > >& parametersToEstimate,
        const Eigen::MatrixXd& covariance )
{
    std::map< std::string, Eigen::Vector3d > uncertainties;
    for( const auto& entry :
         getReferencePointLocalFrameUncertainties< double >( loaded.supplementaryData_, parametersToEstimate, covariance ) )
    {
        uncertainties[ entry.first.second ] = entry.second;
    }
    return uncertainties;
}

//! Pointing (SIGMA_PTG) and landmark (SIGMA_LMK) inverse a-priori entries together.
std::vector< std::pair< EstimatebleParameterIdentifier, Eigen::VectorXd > > landmarkAndPointingAprioriEntries(
        const LoadedPixelTestData& loaded )
{
    std::vector< std::pair< EstimatebleParameterIdentifier, Eigen::VectorXd > > entries = loaded.inverseAprioriCovarianceDiagonalEntries_;
    for( const auto& entry : getReferencePointPositionInverseAprioriCovarianceEntries( loaded.supplementaryData_ ) )
    {
        entries.push_back( entry );
    }
    return entries;
}

Eigen::MatrixXd landmarkInverseAprioriCovariance( const LoadedPixelTestData& loaded,
                                                  const std::shared_ptr< EstimatableParameterSet< double > >& parametersToEstimate )
{
    return createCovarianceFromDiagonalEntries( parametersToEstimate, landmarkAndPointingAprioriEntries( loaded ), false );
}

//! Simulate ideal pixel observations for a scenario and give them unit weight.
std::shared_ptr< ObservationCollection< double, double > > simulateIdealPixelObservations(
        const std::shared_ptr< ObservationCollection< double, double > >& templateCollection,
        OrbitDeterminationManager< double, double >& orbitDeterminationManager,
        const SystemOfBodies& bodies )
{
    const std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > observationSimulationSettings =
            getObservationSimulationSettingsFromObservations< double, double >( templateCollection, bodies );
    const std::shared_ptr< ObservationCollection< double, double > > simulatedObservations = simulateObservations< double, double >(
            observationSimulationSettings, orbitDeterminationManager.getObservationSimulators( ), bodies );
    std::map< std::shared_ptr< ObservationCollectionParser >, double > weightsPerObservationParser;
    weightsPerObservationParser[ observationParser( pixel_coordinates ) ] = 1.0;
    simulatedObservations->setConstantWeightPerObservable( weightsPerObservationParser );
    return simulatedObservations;
}

}  // namespace

BOOST_AUTO_TEST_SUITE( test_pixel_landmark_estimation )

//! Recover an injected per-image pointing offset from pixel-landmark observations, with the
//! spacecraft translational state held fixed. Self-contained Gauss-Newton solve driving the
//! actual CameraPointingCorrection parameter and DirectObservationPartial objects.
BOOST_AUTO_TEST_CASE( testPointingOffsetRecovery )
{
    spice_interface::loadStandardSpiceKernels( );

    const std::vector< std::string > landmarkIds = { "L1", "L2", "L3", "L4" };
    SystemOfBodies bodies = makeBodies( -10000.0 );
    const std::vector< input_output::sum_lmk::SumImageData > images = { makeImage( "IMGA", landmarkIds, -10000.0 ) };
    LoadedPixelTestData conversionResult = loadPixelTestData( images, makeLandmarks( ), bodies );
    const std::string cameraName = input_output::sum_lmk::getSumImageCameraName( "IMGA" );

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

//! Recover an injected per-image pointing offset through the OrbitDeterminationManager while holding the
//! trajectory fixed, using the explicit observation-only/null-propagator mode.
BOOST_AUTO_TEST_CASE( testObservationOnlyPointingRecoveryThroughOrbitDeterminationManager )
{
    spice_interface::loadStandardSpiceKernels( );

    const std::vector< std::string > landmarkIds = { "L1", "L2", "L3", "L4" };
    SystemOfBodies bodies = makeBodies( -10000.0 );
    const std::vector< input_output::sum_lmk::SumImageData > images = { makeImage( "IMGA", landmarkIds, -10000.0 ) };
    LoadedPixelTestData conversionResult = loadPixelTestData( images, makeLandmarks( ), bodies );

    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames = pointingSettingsForTest( conversionResult );
    const std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double, double >( parameterNames, bodies );
    BOOST_REQUIRE_EQUAL( parametersToEstimate->getEstimatedParameterSetSize( ), 3 );

    const std::shared_ptr< PropagatorSettings< double > > nullPropagatorSettings;
    OrbitDeterminationManager< double, double > orbitDeterminationManager(
            bodies, parametersToEstimate, conversionResult.observationModelSettings_, nullPropagatorSettings );

    std::shared_ptr< CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionInterface =
            orbitDeterminationManager.getStateTransitionAndSensitivityMatrixInterface( );
    BOOST_REQUIRE( stateTransitionInterface != nullptr );
    BOOST_CHECK_EQUAL( stateTransitionInterface->getStateTransitionMatrixSize( ), 0 );
    BOOST_CHECK_EQUAL( stateTransitionInterface->getFullParameterVectorSize( ), 3 );
    BOOST_CHECK( orbitDeterminationManager.getVariationalEquationsSolver( ) == nullptr );

    const Eigen::Vector3d trueOffset = ( Eigen::Vector3d( ) << 1.0E-3, -8.0E-4, 1.5E-3 ).finished( );
    Eigen::VectorXd truthParameters = parametersToEstimate->getFullParameterValues< double >( );
    truthParameters.segment( 0, 3 ) = trueOffset;
    parametersToEstimate->resetParameterValues( truthParameters );

    const std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > observationSimulationSettings =
            getObservationSimulationSettingsFromObservations< double, double >( conversionResult.observationCollection_, bodies );
    const std::shared_ptr< ObservationCollection< double, double > > simulatedObservations = simulateObservations< double, double >(
            observationSimulationSettings, orbitDeterminationManager.getObservationSimulators( ), bodies );
    std::map< std::shared_ptr< ObservationCollectionParser >, double > weightsPerObservationParser;
    weightsPerObservationParser[ observationParser( pixel_coordinates ) ] = 1.0;
    simulatedObservations->setConstantWeightPerObservable( weightsPerObservationParser );

    Eigen::VectorXd startParameters = truthParameters;
    startParameters.segment( 0, 3 ) = Eigen::Vector3d::Zero( );
    parametersToEstimate->resetParameterValues( startParameters );

    const std::shared_ptr< EstimationInput< double, double > > estimationInput =
            std::make_shared< EstimationInput< double, double > >( simulatedObservations );
    estimationInput->defineEstimationSettings( true, true, true, false, true, false );
    estimationInput->setConvergenceChecker( std::make_shared< EstimationConvergenceChecker >( 8 ) );

    const std::shared_ptr< EstimationOutput< double > > estimationOutput = orbitDeterminationManager.estimateParameters( estimationInput );
    BOOST_REQUIRE( !estimationOutput->exceptionDuringInversion_ );
    BOOST_REQUIRE( !estimationOutput->exceptionDuringPropagation_ );
    BOOST_CHECK_SMALL( residualRms( estimationOutput->residuals_ ), 1.0E-6 );
    BOOST_CHECK_SMALL( ( estimationOutput->parameterEstimate_ - truthParameters ).norm( ), 1.0E-9 );

    const Eigen::MatrixXd designMatrix = estimationOutput->getUnnormalizedDesignMatrix( );
    BOOST_CHECK_EQUAL( designMatrix.rows( ), simulatedObservations->getTotalObservableSize( ) );
    BOOST_CHECK_EQUAL( designMatrix.cols( ), 3 );
    BOOST_CHECK_GT( designMatrix.norm( ), 0.0 );

    parametersToEstimate->resetParameterValues( startParameters );
    const std::shared_ptr< EstimationInput< double, double > > stateHistoryInput =
            std::make_shared< EstimationInput< double, double > >( simulatedObservations );
    stateHistoryInput->defineEstimationSettings( true, true, true, false, true, true );
    BOOST_CHECK_THROW( orbitDeterminationManager.estimateParameters( stateHistoryInput ), std::runtime_error );
}

//! Recover a DISTINCT injected pointing offset for several images at once, with the trajectory held
//! fixed, through the OrbitDeterminationManager in observation-only mode.
//!
//! This is the case the pointing parameter exists for. With the state fixed there is no
//! pointing/position ridge to slide along, so each image's three angles are determined solely by its own
//! landmarks and the solve should reach the numerical floor without any a-priori. It also exercises the
//! per-image routing through the full estimation stack rather than at the observation-model level: the
//! normal matrix must come out block diagonal, three columns per image, because no observation may
//! depend on another image's pointing.
BOOST_AUTO_TEST_CASE( testObservationOnlyMultiImagePointingRecovery )
{
    spice_interface::loadStandardSpiceKernels( );

    const std::vector< std::string > landmarkIds = { "L1", "L2", "L3", "L4" };
    const std::vector< std::string > imageIds = { "IMGA", "IMGB", "IMGC" };
    SystemOfBodies bodies = makeBodies( -10000.0 );

    std::vector< input_output::sum_lmk::SumImageData > images;
    for( const std::string& imageId : imageIds )
    {
        images.push_back( makeImage( imageId, landmarkIds, -10000.0 ) );
    }

    LoadedPixelTestData conversionResult = loadPixelTestData( images, makeLandmarks( ), bodies );

    const std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames = pointingSettingsForTest( conversionResult );
    const std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double, double >( parameterNames, bodies );
    BOOST_REQUIRE_EQUAL( parametersToEstimate->getEstimatedParameterSetSize( ), 3 * static_cast< int >( imageIds.size( ) ) );

    const std::shared_ptr< PropagatorSettings< double > > nullPropagatorSettings;
    OrbitDeterminationManager< double, double > orbitDeterminationManager(
            bodies, parametersToEstimate, conversionResult.observationModelSettings_, nullPropagatorSettings );

    // Locate each image's parameter block by identifier rather than assuming the vector layout.
    std::vector< int > pointingStartIndices;
    for( const std::string& imageId : imageIds )
    {
        const std::string cameraName = input_output::sum_lmk::getSumImageCameraName( imageId );
        const std::vector< std::pair< int, int > > indices = parametersToEstimate->getIndicesForParameterType(
                EstimatebleParameterIdentifier( camera_pointing_correction, std::make_pair( "Spacecraft", cameraName ) ) );
        BOOST_REQUIRE_EQUAL( indices.size( ), 1 );
        pointingStartIndices.push_back( indices.at( 0 ).first );
    }

    // A distinct offset per image, so a routing error shows up as one image taking another's value.
    Eigen::VectorXd truthParameters = parametersToEstimate->getFullParameterValues< double >( );
    for( std::size_t imageIndex = 0; imageIndex < imageIds.size( ); ++imageIndex )
    {
        truthParameters.segment( pointingStartIndices.at( imageIndex ), 3 ) =
                syntheticPointingOffset( static_cast< int >( imageIndex ), 1.5E-3 );
    }
    parametersToEstimate->resetParameterValues( truthParameters );

    const std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > observationSimulationSettings =
            getObservationSimulationSettingsFromObservations< double, double >( conversionResult.observationCollection_, bodies );
    const std::shared_ptr< ObservationCollection< double, double > > simulatedObservations = simulateObservations< double, double >(
            observationSimulationSettings, orbitDeterminationManager.getObservationSimulators( ), bodies );
    std::map< std::shared_ptr< ObservationCollectionParser >, double > weightsPerObservationParser;
    weightsPerObservationParser[ observationParser( pixel_coordinates ) ] = 1.0;
    simulatedObservations->setConstantWeightPerObservable( weightsPerObservationParser );

    Eigen::VectorXd startParameters = truthParameters;
    for( const int startIndex : pointingStartIndices )
    {
        startParameters.segment( startIndex, 3 ) = Eigen::Vector3d::Zero( );
    }
    parametersToEstimate->resetParameterValues( startParameters );

    const std::shared_ptr< EstimationInput< double, double > > estimationInput =
            std::make_shared< EstimationInput< double, double > >( simulatedObservations );
    estimationInput->defineEstimationSettings( true, true, true, false, true, false );
    estimationInput->setConvergenceChecker( std::make_shared< EstimationConvergenceChecker >( 8 ) );

    const std::shared_ptr< EstimationOutput< double > > estimationOutput = orbitDeterminationManager.estimateParameters( estimationInput );
    BOOST_REQUIRE( !estimationOutput->exceptionDuringInversion_ );

    // No a-priori is used: with the trajectory fixed the data alone determines every angle.
    BOOST_CHECK_SMALL( residualRms( estimationOutput->residuals_ ), 1.0E-6 );
    for( std::size_t imageIndex = 0; imageIndex < imageIds.size( ); ++imageIndex )
    {
        const int startIndex = pointingStartIndices.at( imageIndex );
        const Eigen::Vector3d recovered = estimationOutput->parameterEstimate_.segment( startIndex, 3 );
        const Eigen::Vector3d injected = truthParameters.segment( startIndex, 3 );
        BOOST_CHECK_SMALL( ( recovered - injected ).norm( ), 1.0E-9 );
        // Guard against a routing error that happens to look small: the offsets must differ per image.
        if( imageIndex > 0 )
        {
            const Eigen::Vector3d previous = truthParameters.segment( pointingStartIndices.at( imageIndex - 1 ), 3 );
            BOOST_CHECK_GT( ( injected - previous ).norm( ), 1.0E-4 );
        }
    }

    // Per-image independence through the full stack: no observation may depend on another image's
    // pointing, so the normal matrix must be block diagonal with one 3x3 block per image. Checking the
    // normal matrix rather than the design matrix keeps this independent of observation row ordering.
    const Eigen::MatrixXd designMatrix = estimationOutput->getUnnormalizedDesignMatrix( );
    BOOST_REQUIRE_EQUAL( designMatrix.cols( ), 3 * static_cast< int >( imageIds.size( ) ) );
    const Eigen::MatrixXd normalMatrix = designMatrix.transpose( ) * designMatrix;
    for( std::size_t row = 0; row < imageIds.size( ); ++row )
    {
        const Eigen::Matrix3d diagonalBlock = normalMatrix.block( pointingStartIndices.at( row ), pointingStartIndices.at( row ), 3, 3 );
        BOOST_CHECK_GT( diagonalBlock.norm( ), 0.0 );
        for( std::size_t column = 0; column < imageIds.size( ); ++column )
        {
            if( row == column )
            {
                continue;
            }
            const Eigen::Matrix3d offDiagonalBlock =
                    normalMatrix.block( pointingStartIndices.at( row ), pointingStartIndices.at( column ), 3, 3 );
            BOOST_CHECK_SMALL( offDiagonalBlock.norm( ) / diagonalBlock.norm( ), 1.0E-12 );
        }
    }
}

//! The SIGMA_PTG a-priori must assemble and apply in observation-only mode too, where the parameter
//! vector has no initial-state block and every index belongs to a pointing parameter.
BOOST_AUTO_TEST_CASE( testObservationOnlyPointingApriori )
{
    spice_interface::loadStandardSpiceKernels( );

    const std::vector< std::string > landmarkIds = { "L1", "L2", "L3", "L4" };
    SystemOfBodies bodies = makeBodies( -10000.0 );
    const std::vector< input_output::sum_lmk::SumImageData > images = { makeImage( "IMGA", landmarkIds, -10000.0 ),
                                                                        makeImage( "IMGB", landmarkIds, -10000.0 ) };
    LoadedPixelTestData conversionResult = loadPixelTestData( images, makeLandmarks( ), bodies );

    const std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double, double >( pointingSettingsForTest( conversionResult ), bodies );
    BOOST_REQUIRE_EQUAL( parametersToEstimate->getEstimatedParameterSetSize( ), 6 );

    // makeImage sets SIGMA_PTG = 1e-4 rad on all three axes.
    const Eigen::MatrixXd inverseAprioriCovariance =
            createCovarianceFromDiagonalEntries( parametersToEstimate, conversionResult.inverseAprioriCovarianceDiagonalEntries_, false );
    BOOST_REQUIRE_EQUAL( inverseAprioriCovariance.rows( ), 6 );
    BOOST_REQUIRE_EQUAL( inverseAprioriCovariance.cols( ), 6 );
    const double expectedInverseVariance = 1.0 / ( 1.0E-4 * 1.0E-4 );
    for( int i = 0; i < 6; ++i )
    {
        BOOST_CHECK_CLOSE( inverseAprioriCovariance( i, i ), expectedInverseVariance, 1.0E-8 );
    }
    BOOST_CHECK_SMALL( ( inverseAprioriCovariance - Eigen::MatrixXd( inverseAprioriCovariance.diagonal( ).asDiagonal( ) ) ).norm( ) /
                               expectedInverseVariance,
                       1.0E-12 );

    // And it must be accepted by an observation-only estimation run.
    const std::shared_ptr< PropagatorSettings< double > > nullPropagatorSettings;
    OrbitDeterminationManager< double, double > orbitDeterminationManager(
            bodies, parametersToEstimate, conversionResult.observationModelSettings_, nullPropagatorSettings );

    const std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > observationSimulationSettings =
            getObservationSimulationSettingsFromObservations< double, double >( conversionResult.observationCollection_, bodies );
    const std::shared_ptr< ObservationCollection< double, double > > simulatedObservations = simulateObservations< double, double >(
            observationSimulationSettings, orbitDeterminationManager.getObservationSimulators( ), bodies );
    std::map< std::shared_ptr< ObservationCollectionParser >, double > weightsPerObservationParser;
    weightsPerObservationParser[ observationParser( pixel_coordinates ) ] = 1.0;
    simulatedObservations->setConstantWeightPerObservable( weightsPerObservationParser );

    const std::shared_ptr< EstimationInput< double, double > > estimationInput =
            std::make_shared< EstimationInput< double, double > >( simulatedObservations, inverseAprioriCovariance );
    estimationInput->defineEstimationSettings( true, true, true, false, true, false );
    estimationInput->setConvergenceChecker( std::make_shared< EstimationConvergenceChecker >( 5 ) );

    const std::shared_ptr< EstimationOutput< double > > estimationOutput = orbitDeterminationManager.estimateParameters( estimationInput );
    BOOST_REQUIRE( !estimationOutput->exceptionDuringInversion_ );
    // The data was simulated at zero correction, which is also where the a-priori pulls, so the solution
    // must stay at zero rather than being displaced by a misindexed a-priori.
    BOOST_CHECK_SMALL( estimationOutput->parameterEstimate_.norm( ), 1.0E-9 );
}

//! A covariance analysis (no parameter update) must run in observation-only mode: the path takes a
//! CovarianceAnalysisInput rather than an EstimationInput, so none of the propagation-dependent
//! branches of the estimation loop may be entered.
BOOST_AUTO_TEST_CASE( testObservationOnlyCovarianceAnalysis )
{
    spice_interface::loadStandardSpiceKernels( );

    const std::vector< std::string > landmarkIds = { "L1", "L2", "L3", "L4" };
    SystemOfBodies bodies = makeBodies( -10000.0 );
    const std::vector< input_output::sum_lmk::SumImageData > images = { makeImage( "IMGA", landmarkIds, -10000.0 ) };
    LoadedPixelTestData conversionResult = loadPixelTestData( images, makeLandmarks( ), bodies );

    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames = pointingSettingsForTest( conversionResult );
    const std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double, double >( parameterNames, bodies );
    BOOST_REQUIRE_EQUAL( parametersToEstimate->getEstimatedParameterSetSize( ), 3 );

    const std::shared_ptr< PropagatorSettings< double > > nullPropagatorSettings;
    OrbitDeterminationManager< double, double > orbitDeterminationManager(
            bodies, parametersToEstimate, conversionResult.observationModelSettings_, nullPropagatorSettings );
    BOOST_REQUIRE( orbitDeterminationManager.getVariationalEquationsSolver( ) == nullptr );

    const std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > observationSimulationSettings =
            getObservationSimulationSettingsFromObservations< double, double >( conversionResult.observationCollection_, bodies );
    const std::shared_ptr< ObservationCollection< double, double > > simulatedObservations = simulateObservations< double, double >(
            observationSimulationSettings, orbitDeterminationManager.getObservationSimulators( ), bodies );
    std::map< std::shared_ptr< ObservationCollectionParser >, double > weightsPerObservationParser;
    weightsPerObservationParser[ observationParser( pixel_coordinates ) ] = 1.0;
    simulatedObservations->setConstantWeightPerObservable( weightsPerObservationParser );

    const std::shared_ptr< CovarianceAnalysisInput< double, double > > covarianceInput =
            std::make_shared< CovarianceAnalysisInput< double, double > >( simulatedObservations );
    covarianceInput->defineCovarianceSettings( true, true, true, false );

    const std::shared_ptr< CovarianceAnalysisOutput< double, double > > covarianceOutput =
            orbitDeterminationManager.computeCovariance( covarianceInput );
    BOOST_REQUIRE( covarianceOutput != nullptr );
    BOOST_CHECK( !covarianceOutput->exceptionDuringPropagation_ );

    // The design matrix must be built purely from observation partials: as many rows as observations,
    // one column per pointing parameter, and no zero column (each parameter is observable).
    const Eigen::MatrixXd designMatrix = covarianceOutput->getUnnormalizedDesignMatrix( );
    BOOST_REQUIRE_EQUAL( designMatrix.rows( ), simulatedObservations->getTotalObservableSize( ) );
    BOOST_REQUIRE_EQUAL( designMatrix.cols( ), 3 );
    for( int i = 0; i < designMatrix.cols( ); i++ )
    {
        BOOST_CHECK_GT( designMatrix.col( i ).norm( ), 0.0 );
    }

    // The resulting covariance must be a finite, symmetric, positive-definite 3x3 matrix.
    const Eigen::MatrixXd covariance = covarianceOutput->getUnnormalizedCovarianceMatrix( );
    BOOST_REQUIRE_EQUAL( covariance.rows( ), 3 );
    BOOST_REQUIRE_EQUAL( covariance.cols( ), 3 );
    BOOST_REQUIRE( covariance.allFinite( ) );
    BOOST_CHECK_SMALL( ( covariance - covariance.transpose( ) ).norm( ) / covariance.norm( ), 1.0E-12 );
    const Eigen::LLT< Eigen::MatrixXd > covarianceCholesky( covariance );
    BOOST_CHECK( covarianceCholesky.info( ) == Eigen::Success );
    for( int i = 0; i < covariance.rows( ); i++ )
    {
        BOOST_CHECK_GT( covariance( i, i ), 0.0 );
    }

    // Covariance analysis is the inverse of the weighted normal matrix built from the same design
    // matrix, so recomputing it by hand must reproduce the manager's result.
    const Eigen::MatrixXd expectedInverseCovariance = designMatrix.transpose( ) * designMatrix;
    const Eigen::MatrixXd actualInverseCovariance = covarianceOutput->getUnnormalizedInverseCovarianceMatrix( );
    BOOST_REQUIRE_EQUAL( actualInverseCovariance.rows( ), 3 );
    BOOST_CHECK_SMALL( ( actualInverseCovariance - expectedInverseCovariance ).norm( ) / expectedInverseCovariance.norm( ), 1.0E-10 );
}

//! Two images produce two independent pointing parameters: perturbing one image's correction
//! must not change the other image's modelled observations.
BOOST_AUTO_TEST_CASE( testMultiImagePointingIndependence )
{
    spice_interface::loadStandardSpiceKernels( );

    SystemOfBodies bodies = makeBodies( -10000.0 );
    const std::vector< input_output::sum_lmk::SumImageData > images = { makeImage( "IMGA", { "L1", "L2" }, -10000.0 ),
                                                                        makeImage( "IMGB", { "L1", "L2" }, -10000.0 ) };
    LoadedPixelTestData conversionResult = loadPixelTestData( images, makeLandmarks( ), bodies );
    const std::string cameraA = input_output::sum_lmk::getSumImageCameraName( "IMGA" );
    const std::string cameraB = input_output::sum_lmk::getSumImageCameraName( "IMGB" );
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
    LoadedPixelTestData conversionResult = loadPixelTestData( images, makeLandmarks( ), bodies );

    const std::vector< std::shared_ptr< EstimatableParameterSettings > > allSettings = pointingSettingsForTest( conversionResult );
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
    BOOST_CHECK_EQUAL( settingCameraNames.at( 0 ), input_output::sum_lmk::getSumImageCameraName( "IMGA" ) );
    BOOST_CHECK_EQUAL( settingCameraNames.at( 1 ), input_output::sum_lmk::getSumImageCameraName( "IMGB" ) );

    // Restricting to a subset of the images.
    const std::vector< std::shared_ptr< EstimatableParameterSettings > > subsetSettings = { cameraPointingCorrection( "Spacecraft",
                                                                                                                      "Camera_IMGB" ) };
    BOOST_REQUIRE_EQUAL( subsetSettings.size( ), 1 );
    BOOST_CHECK_EQUAL( subsetSettings.at( 0 )->parameterType_.second.second, input_output::sum_lmk::getSumImageCameraName( "IMGB" ) );

    // A parameter for an unknown camera is rejected by the ordinary parameter factory.
    BOOST_CHECK_THROW(
            ( createParametersToEstimate< double, double >( { cameraPointingCorrection( "Spacecraft", "Camera_NOT_AN_IMAGE" ) }, bodies ) ),
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
    LoadedPixelTestData conversionResult = loadPixelTestData( images, makeLandmarks( ), bodies );
    BOOST_REQUIRE_EQUAL( conversionResult.inverseAprioriCovarianceDiagonalEntries_.size( ), 2 );

    // makeImage uses SIGMA_PTG = 1e-4 rad on all three axes.
    const double expectedInverseVariance = 1.0 / ( 1.0E-4 * 1.0E-4 );

    // --- Both images estimated: both blocks are filled.
    {
        const std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate =
                createParametersToEstimate< double, double >( pointingSettingsForTest( conversionResult ), bodies );
        const int numberOfParameters = parametersToEstimate->getEstimatedParameterSetSize( );
        BOOST_REQUIRE_EQUAL( numberOfParameters, 6 );

        const Eigen::MatrixXd inverseAprioriCovariance = createCovarianceFromDiagonalEntries(
                parametersToEstimate, conversionResult.inverseAprioriCovarianceDiagonalEntries_, false );
        BOOST_REQUIRE_EQUAL( inverseAprioriCovariance.rows( ), numberOfParameters );
        BOOST_REQUIRE_EQUAL( inverseAprioriCovariance.cols( ), numberOfParameters );

        // Diagonal carries the inverse variances; the matrix has no off-diagonal content.
        BOOST_CHECK_SMALL( ( inverseAprioriCovariance.diagonal( ) - Eigen::VectorXd::Constant( 6, expectedInverseVariance ) ).norm( ) /
                                   expectedInverseVariance,
                           1.0E-12 );
        BOOST_CHECK_SMALL( ( inverseAprioriCovariance - Eigen::MatrixXd( inverseAprioriCovariance.diagonal( ).asDiagonal( ) ) ).norm( ),
                           1.0E-12 );

        // The generic helper assigns the requested diagonal and preserves unrelated entries.
        const Eigen::MatrixXd baseCovariance = Eigen::MatrixXd::Identity( numberOfParameters, numberOfParameters );
        const Eigen::MatrixXd combined = addCovarianceDiagonalEntries(
                baseCovariance, parametersToEstimate, conversionResult.inverseAprioriCovarianceDiagonalEntries_ );
        BOOST_CHECK_SMALL( ( combined - inverseAprioriCovariance ).norm( ) / expectedInverseVariance, 1.0E-12 );

        // A base matrix of the wrong size is an error.
        BOOST_CHECK_THROW( ( addCovarianceDiagonalEntries( Eigen::MatrixXd::Identity( 3, 3 ),
                                                           parametersToEstimate,
                                                           conversionResult.inverseAprioriCovarianceDiagonalEntries_ ) ),
                           std::runtime_error );
    }

    // --- Only one image estimated: the other image's a-priori is skipped, not an error.
    {
        const std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate =
                createParametersToEstimate< double, double >( { cameraPointingCorrection( "Spacecraft", "Camera_IMGB" ) }, bodies );
        BOOST_REQUIRE_EQUAL( parametersToEstimate->getEstimatedParameterSetSize( ), 3 );

        const Eigen::MatrixXd inverseAprioriCovariance = createCovarianceFromDiagonalEntries(
                parametersToEstimate, conversionResult.inverseAprioriCovarianceDiagonalEntries_, false );
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

        const double epoch = sumImageTdbTimeForTest( image );
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
    LoadedPixelTestData conversionResult = loadPixelTestData( images, landmarks, bodies );

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
    for( const std::shared_ptr< EstimatableParameterSettings >& pointingSetting : pointingSettingsForTest( scenario.conversionResult_ ) )
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
        const std::string cameraName = input_output::sum_lmk::getSumImageCameraName( scenario.imageIds_.at( imageIndex ) );
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
        const std::string cameraName = input_output::sum_lmk::getSumImageCameraName( scenario.imageIds_.at( imageIndex ) );
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
        const std::string cameraName = input_output::sum_lmk::getSumImageCameraName( scenario.imageIds_.at( imageIndex ) );
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
//! uncertainty. Both effects check that the generic a-priori covariance reaches the normal
//! equations at the right parameter indices -
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
    for( const std::shared_ptr< EstimatableParameterSettings >& pointingSetting : pointingSettingsForTest( scenario.conversionResult_ ) )
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
        const std::string cameraName = input_output::sum_lmk::getSumImageCameraName( imageId );
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
    const Eigen::MatrixXd inverseAprioriCovariance = createCovarianceFromDiagonalEntries(
            parametersToEstimate, scenario.conversionResult_.inverseAprioriCovarianceDiagonalEntries_, false );
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
    for( const std::string& imageId : estimatedImageIds )
    {
        parameterNames.push_back( cameraPointingCorrection( "Spacecraft", input_output::sum_lmk::getSumImageCameraName( imageId ) ) );
    }
    const std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double, double >( parameterNames, scenario.bodies_, scenario.propagatorSettings_ );
    BOOST_REQUIRE_EQUAL( parametersToEstimate->getEstimatedParameterSetSize( ), 6 + 3 * 2 );

    // All four a-priori entries are offered; only the two estimated ones may land in the matrix.
    const Eigen::MatrixXd inverseAprioriCovariance = createCovarianceFromDiagonalEntries(
            parametersToEstimate, scenario.conversionResult_.inverseAprioriCovarianceDiagonalEntries_, false );
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
    for( const std::shared_ptr< EstimatableParameterSettings >& pointingSetting : pointingSettingsForTest( scenario.conversionResult_ ) )
    {
        parameterNames.push_back( pointingSetting );
    }
    const std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double, double >( parameterNames, scenario.bodies_, scenario.propagatorSettings_ );

    const std::size_t numberOfImages = pointingSettingsForTest( scenario.conversionResult_ ).size( );
    BOOST_REQUIRE_EQUAL( parametersToEstimate->getEstimatedParameterSetSize( ), 6 + 3 * numberOfImages );

    const Eigen::MatrixXd inverseAprioriCovariance = createCovarianceFromDiagonalEntries(
            parametersToEstimate, scenario.conversionResult_.inverseAprioriCovarianceDiagonalEntries_, false );
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
        // The conversion also emits landmark SIGMA_LMK entries; only the pointing ones are checked here.
        if( entry.first.first != camera_pointing_correction )
        {
            continue;
        }
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

//! Landmark position estimation: landmarks are registered as body-fixed ground stations by the SUM/LMK
//! conversion, so their positions are estimated through the ordinary ground_station_position parameter.
//! Because that parameter is not arc-wise, each landmark is a single global parameter shared by every arc
//! that observes it, which is what makes a landmark seen from several arcs a tie between them.
BOOST_AUTO_TEST_SUITE( test_pixel_landmark_position_estimation )

//! Recover perturbed landmark positions with the trajectory held fixed, through the
//! OrbitDeterminationManager in observation-only mode.
//!
//! This is the base case the parameter exists for, and the first test to exercise
//! ground_station_position through a pixel observable at all: the partial comes from the generic
//! vector-parameter route and the ground_station_position case of createCartesianStatePartialsWrtParameter,
//! with no pixel-specific code involved. With the trajectory fixed there is no network/orbit degeneracy to
//! slide along, so noise-free data determines all three body-fixed components of every landmark.
BOOST_AUTO_TEST_CASE( testObservationOnlyLandmarkPositionRecovery )
{
    spice_interface::loadStandardSpiceKernels( );

    LandmarkObservationScenario scenario = buildLandmarkObservationScenario( 12 );

    const std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
            landmarkParameterSettings( scenario.conversionResult_ );
    BOOST_REQUIRE_EQUAL( parameterNames.size( ), scenario.landmarks_.size( ) );

    const std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double, double >( parameterNames, scenario.bodies_ );
    BOOST_REQUIRE_EQUAL( parametersToEstimate->getEstimatedParameterSetSize( ), static_cast< int >( 3 * scenario.landmarks_.size( ) ) );

    // Exactly one parameter block per landmark: a landmark is one parameter, however many images see it.
    for( const auto& landmarkEntry : scenario.landmarks_ )
    {
        const std::vector< std::pair< int, int > > indices = parametersToEstimate->getIndicesForParameterType(
                EstimatebleParameterIdentifier( ground_station_position, std::make_pair( "Target", landmarkEntry.first ) ) );
        BOOST_REQUIRE_EQUAL( indices.size( ), 1 );
        BOOST_CHECK_EQUAL( indices.at( 0 ).second, 3 );
    }

    const std::shared_ptr< PropagatorSettings< double > > nullPropagatorSettings;
    OrbitDeterminationManager< double, double > orbitDeterminationManager(
            scenario.bodies_, parametersToEstimate, scenario.conversionResult_.observationModelSettings_, nullPropagatorSettings );

    // Truth = the positions the LMK data declares, which the conversion has already installed.
    const Eigen::VectorXd truthParameters = parametersToEstimate->getFullParameterValues< double >( );
    const std::shared_ptr< ObservationCollection< double, double > > simulatedObservations = simulateIdealPixelObservations(
            scenario.conversionResult_.observationCollection_, orbitDeterminationManager, scenario.bodies_ );

    // Perturb every landmark by a metre-level offset, distinct per landmark and per component.
    Eigen::VectorXd perturbedParameters = truthParameters;
    for( int i = 0; i < perturbedParameters.size( ); ++i )
    {
        perturbedParameters( i ) += 5.0 * std::sin( 0.9 * i + 0.4 );
    }
    parametersToEstimate->resetParameterValues( perturbedParameters );

    const std::shared_ptr< EstimationInput< double, double > > estimationInput =
            std::make_shared< EstimationInput< double, double > >( simulatedObservations );
    estimationInput->defineEstimationSettings( true, true, true, false, true, false );
    estimationInput->setConvergenceChecker( std::make_shared< EstimationConvergenceChecker >( 10 ) );

    const std::shared_ptr< EstimationOutput< double > > estimationOutput = orbitDeterminationManager.estimateParameters( estimationInput );
    BOOST_REQUIRE( !estimationOutput->exceptionDuringInversion_ );
    BOOST_CHECK_SMALL( residualRms( estimationOutput->residuals_ ), 1.0E-6 );

    const Eigen::VectorXd estimationError = estimationOutput->parameterEstimate_ - truthParameters;
    BOOST_TEST_MESSAGE( "Landmark-only solve: worst component error " << estimationError.cwiseAbs( ).maxCoeff( ) << " m" );
    BOOST_CHECK_SMALL( estimationError.cwiseAbs( ).maxCoeff( ), 1.0E-5 );
}

//! Only landmarks that the estimated observations actually observe may get a parameter.
//!
//! Images, and therefore landmarks, are dropped before the estimation starts. This checks the three ways
//! that happens: a landmark no image references, an image dropped by skipObservationsWithMissingLandmarks,
//! and observations removed from the collection afterwards by filterObservations - including the case where
//! filtering empties a set but leaves it in place, which a helper that merely listed link ends would miss.
BOOST_AUTO_TEST_CASE( testLandmarkParameterSettingsExcludeFilteredLandmarks )
{
    spice_interface::loadStandardSpiceKernels( );

    // --- Landmarks that no image references are not estimable, even though LMK data exists for them.
    {
        SystemOfBodies bodies = makeBodies( -10000.0 );
        const std::vector< input_output::sum_lmk::SumImageData > images = { makeImage( "IMGA", { "L1", "L3" }, -10000.0 ) };
        // makeLandmarks( ) defines L1..L4; only L1 and L3 are observed.
        LoadedPixelTestData conversionResult = loadPixelTestData( images, makeLandmarks( ), bodies );

        const std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames = landmarkParameterSettings( conversionResult );
        BOOST_REQUIRE_EQUAL( parameterNames.size( ), 2 );
        // Ordered by landmark ID, so the parameter vector layout is reproducible.
        BOOST_CHECK_EQUAL( parameterNames.at( 0 )->parameterType_.second.second, "L1" );
        BOOST_CHECK_EQUAL( parameterNames.at( 1 )->parameterType_.second.second, "L3" );
        BOOST_CHECK_EQUAL( parameterNames.at( 0 )->parameterType_.first, ground_station_position );
        BOOST_CHECK_EQUAL( parameterNames.at( 0 )->parameterType_.second.first, "Target" );

        // An unobserved landmark, and an entirely unknown one, must both raise rather than be dropped.
        BOOST_CHECK_THROW( landmarkParameterSettings( conversionResult, nullptr, { "L2" } ), std::runtime_error );
        BOOST_CHECK_THROW( landmarkParameterSettings( conversionResult, nullptr, { "NOPE" } ), std::runtime_error );
        // A subset of the observed landmarks is fine.
        BOOST_CHECK_EQUAL( landmarkParameterSettings( conversionResult, nullptr, { "L3" } ).size( ), 1 );
    }

    // --- An image dropped by skipObservationsWithMissingLandmarks takes its landmarks with it.
    {
        SystemOfBodies bodies = makeBodies( -10000.0 );
        // IMGB references only L9, for which no LMK data exists, so the whole image is skipped.
        const std::vector< input_output::sum_lmk::SumImageData > images = { makeImage( "IMGA", { "L1" }, -10000.0 ),
                                                                            makeImage( "IMGB", { "L9" }, -10000.0 ) };
        LoadedPixelTestData conversionResult = loadPixelTestData( images, makeLandmarks( ), bodies, "Target", "Spacecraft", true );

        const std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames = landmarkParameterSettings( conversionResult );
        BOOST_REQUIRE_EQUAL( parameterNames.size( ), 1 );
        BOOST_CHECK_EQUAL( parameterNames.at( 0 )->parameterType_.second.second, "L1" );
    }

    // --- Observations removed from the collection after conversion.
    {
        SystemOfBodies bodies = makeBodies( -10000.0 );
        input_output::sum_lmk::SumImageData earlyImage = makeImage( "IMGEARLY", { "L1" }, -10000.0 );
        input_output::sum_lmk::SumImageData lateImage = makeImage( "IMGLATE", { "L2" }, -10000.0 );
        lateImage.utcEpochString_ = "2015 JUN 05 08:24:42.053";

        LoadedPixelTestData conversionResult = loadPixelTestData( { earlyImage, lateImage }, makeLandmarks( ), bodies );
        BOOST_REQUIRE_EQUAL( landmarkParameterSettings( conversionResult ).size( ), 2 );

        const double earlyTime = sumImageTdbTimeForTest( earlyImage );
        conversionResult.observationCollection_->filterObservations(
                observationFilter( time_bounds_filtering, earlyTime - 1.0, earlyTime + 1.0, true, true ) );

        // Deliberately NOT calling removeEmptySingleObservationSets( ): L2's set is still present but now
        // holds no observations, so a helper that listed link ends rather than counting observations would
        // still hand back a parameter for a landmark the solution cannot see.
        const std::vector< std::shared_ptr< EstimatableParameterSettings > > filteredParameterNames =
                landmarkParameterSettings( conversionResult );
        BOOST_REQUIRE_EQUAL( filteredParameterNames.size( ), 1 );
        BOOST_CHECK_EQUAL( filteredParameterNames.at( 0 )->parameterType_.second.second, "L1" );

        // Requesting the filtered-away landmark explicitly must raise.
        BOOST_CHECK_THROW( landmarkParameterSettings( conversionResult, nullptr, { "L2" } ), std::runtime_error );

        // Passing a separate filtered collection must give the same answer as filtering in place.
        conversionResult.observationCollection_->removeEmptySingleObservationSets( );
        BOOST_CHECK_EQUAL( landmarkParameterSettings( conversionResult, conversionResult.observationCollection_ ).size( ), 1 );
    }
}

//! The SIGMA_LMK a-priori must reach the normal equations at the right parameter indices.
//!
//! SIGMA_LMK is a per-component sigma on the landmark's BODY-FIXED position, which is the frame
//! ground_station_position is expressed in, so the a-priori is a plain diagonal of inverse variances. An
//! a-priori assembled at the wrong offsets would constrain the wrong parameters and leave the landmarks
//! unconstrained, which is what the index checks below rule out.
BOOST_AUTO_TEST_CASE( testLandmarkAprioriFromSigmaLmk )
{
    spice_interface::loadStandardSpiceKernels( );

    const double landmarkSigma = 0.5;
    LandmarkObservationScenario scenario = buildLandmarkObservationScenario( 4, landmarkSigma );
    const double expectedInverseVariance = 1.0 / ( landmarkSigma * landmarkSigma );

    // One landmark a-priori entry per landmark, plus one pointing entry per image only where SIGMA_PTG is
    // given - the scenario leaves SIGMA_PTG unset, so every entry here is a landmark entry.
    std::size_t numberOfLandmarkEntries = 0;
    for( const auto& entry : landmarkAndPointingAprioriEntries( scenario.conversionResult_ ) )
    {
        if( entry.first.first == ground_station_position )
        {
            ++numberOfLandmarkEntries;
            BOOST_CHECK_EQUAL( entry.first.second.first, "Target" );
            BOOST_REQUIRE_EQUAL( entry.second.size( ), 3 );
            for( int component = 0; component < 3; ++component )
            {
                BOOST_CHECK_CLOSE( entry.second( component ), expectedInverseVariance, 1.0E-9 );
            }
        }
    }
    BOOST_CHECK_EQUAL( numberOfLandmarkEntries, scenario.landmarks_.size( ) );

    // --- Every landmark estimated: each gets its own diagonal block.
    {
        const std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate =
                createParametersToEstimate< double, double >( landmarkParameterSettings( scenario.conversionResult_ ), scenario.bodies_ );
        const Eigen::MatrixXd inverseAprioriCovariance =
                landmarkInverseAprioriCovariance( scenario.conversionResult_, parametersToEstimate );

        const int numberOfParameters = parametersToEstimate->getEstimatedParameterSetSize( );
        BOOST_REQUIRE_EQUAL( inverseAprioriCovariance.rows( ), numberOfParameters );
        // Purely diagonal, every entry the same inverse variance.
        BOOST_CHECK_SMALL( ( inverseAprioriCovariance - inverseAprioriCovariance.diagonal( ).asDiagonal( ).toDenseMatrix( ) ).norm( ),
                           1.0E-9 );
        for( int i = 0; i < numberOfParameters; ++i )
        {
            BOOST_CHECK_CLOSE( inverseAprioriCovariance( i, i ), expectedInverseVariance, 1.0E-9 );
        }
    }

    // --- Only a subset estimated: entries for landmarks outside the parameter set are skipped, not
    //     written at whatever index happens to be free.
    {
        std::vector< std::string > landmarkIds;
        for( const auto& entry : scenario.landmarks_ )
        {
            landmarkIds.push_back( entry.first );
        }
        const std::vector< std::string > estimatedLandmarkIds = { landmarkIds.at( 1 ), landmarkIds.at( 3 ) };
        const std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate = createParametersToEstimate< double, double >(
                landmarkParameterSettings( scenario.conversionResult_, nullptr, estimatedLandmarkIds ), scenario.bodies_ );
        BOOST_REQUIRE_EQUAL( parametersToEstimate->getEstimatedParameterSetSize( ), 6 );

        const Eigen::MatrixXd inverseAprioriCovariance =
                landmarkInverseAprioriCovariance( scenario.conversionResult_, parametersToEstimate );
        BOOST_REQUIRE_EQUAL( inverseAprioriCovariance.rows( ), 6 );
        for( int i = 0; i < 6; ++i )
        {
            BOOST_CHECK_CLOSE( inverseAprioriCovariance( i, i ), expectedInverseVariance, 1.0E-9 );
        }
    }

    // --- The a-priori must actually pull the solution: with a tight SIGMA_LMK the estimate stays near the
    //     a-priori value even when the data alone would move it further.
    {
        const std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate =
                createParametersToEstimate< double, double >( landmarkParameterSettings( scenario.conversionResult_ ), scenario.bodies_ );
        const std::shared_ptr< PropagatorSettings< double > > nullPropagatorSettings;
        OrbitDeterminationManager< double, double > orbitDeterminationManager(
                scenario.bodies_, parametersToEstimate, scenario.conversionResult_.observationModelSettings_, nullPropagatorSettings );
        const Eigen::VectorXd truthParameters = parametersToEstimate->getFullParameterValues< double >( );
        const std::shared_ptr< ObservationCollection< double, double > > simulatedObservations = simulateIdealPixelObservations(
                scenario.conversionResult_.observationCollection_, orbitDeterminationManager, scenario.bodies_ );

        const std::shared_ptr< EstimationInput< double, double > > estimationInput = std::make_shared< EstimationInput< double, double > >(
                simulatedObservations, landmarkInverseAprioriCovariance( scenario.conversionResult_, parametersToEstimate ) );
        estimationInput->defineEstimationSettings( true, true, true, false, true, false );
        estimationInput->setConvergenceChecker( std::make_shared< EstimationConvergenceChecker >( 6 ) );

        const std::shared_ptr< EstimationOutput< double > > estimationOutput =
                orbitDeterminationManager.estimateParameters( estimationInput );
        BOOST_REQUIRE( !estimationOutput->exceptionDuringInversion_ );

        // The a-priori shrinks the formal uncertainty below the prior sigma itself: the data adds
        // information on top of it rather than replacing it.
        const Eigen::VectorXd formalErrors = estimationOutput->getFormalErrorVector( );
        BOOST_REQUIRE_EQUAL( formalErrors.size( ), truthParameters.size( ) );
        BOOST_CHECK_LT( formalErrors.maxCoeff( ), landmarkSigma );
        BOOST_CHECK_GT( formalErrors.minCoeff( ), 0.0 );
    }
}

//! The formal uncertainties can be reported along each landmark's own local frame, which is the one
//! genuinely useful property of a local-frame parameterisation - and is available without one, because
//! rotating the covariance block afterwards is equivalent.
BOOST_AUTO_TEST_CASE( testLandmarkLocalFrameUncertainties )
{
    spice_interface::loadStandardSpiceKernels( );

    LandmarkObservationScenario scenario = buildLandmarkObservationScenario( 8, 0.5 );
    const std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double, double >( landmarkParameterSettings( scenario.conversionResult_ ), scenario.bodies_ );
    const int numberOfParameters = parametersToEstimate->getEstimatedParameterSetSize( );

    // A synthetic but non-diagonal covariance, so that the rotation genuinely changes the answer: a
    // diagonal input with equal entries would be invariant and would not test anything.
    Eigen::MatrixXd covariance = Eigen::MatrixXd::Zero( numberOfParameters, numberOfParameters );
    for( int i = 0; i < numberOfParameters; ++i )
    {
        for( int j = 0; j < numberOfParameters; ++j )
        {
            covariance( i, j ) = ( i == j ) ? ( 0.2 + 0.05 * i ) : ( 0.01 * std::cos( 0.5 * ( i + j ) ) );
        }
    }
    covariance = ( covariance * covariance.transpose( ) ).eval( );  // symmetric positive definite

    const std::map< std::string, Eigen::Vector3d > localFrameUncertainties =
            landmarkLocalFrameUncertainties( scenario.conversionResult_, parametersToEstimate, covariance );
    BOOST_REQUIRE_EQUAL( localFrameUncertainties.size( ), scenario.landmarks_.size( ) );

    for( const auto& landmarkEntry : scenario.landmarks_ )
    {
        const std::vector< std::pair< int, int > > indices = parametersToEstimate->getIndicesForParameterType(
                EstimatebleParameterIdentifier( ground_station_position, std::make_pair( "Target", landmarkEntry.first ) ) );
        BOOST_REQUIRE_EQUAL( indices.size( ), 1 );
        const int startIndex = indices.at( 0 ).first;

        Eigen::Matrix3d localFrame;
        localFrame.col( 0 ) = landmarkEntry.second.localXAxis_;
        localFrame.col( 1 ) = landmarkEntry.second.localYAxis_;
        localFrame.col( 2 ) = landmarkEntry.second.localZAxis_;
        const Eigen::Matrix3d expectedLocalCovariance =
                localFrame.transpose( ) * covariance.block( startIndex, startIndex, 3, 3 ) * localFrame;

        const Eigen::Vector3d reported = localFrameUncertainties.at( landmarkEntry.first );
        for( int i = 0; i < 3; ++i )
        {
            BOOST_CHECK_CLOSE( reported( i ), std::sqrt( expectedLocalCovariance( i, i ) ), 1.0E-9 );
        }

        // The local-frame sigmas must differ from the body-fixed ones, otherwise the rotation is a no-op
        // and the test would pass even if the local frame were ignored entirely.
        const Eigen::Vector3d bodyFixedSigmas = covariance.block( startIndex, startIndex, 3, 3 ).diagonal( ).cwiseSqrt( );
        BOOST_CHECK_GT( ( reported - bodyFixedSigmas ).norm( ), 1.0E-6 );

        // A rotation preserves the total variance, whatever axes it is expressed in.
        BOOST_CHECK_CLOSE( reported.squaredNorm( ), bodyFixedSigmas.squaredNorm( ), 1.0E-6 );
    }

    // Landmarks that are not estimated are omitted rather than reported as zero.
    {
        std::vector< std::string > landmarkIds;
        for( const auto& entry : scenario.landmarks_ )
        {
            landmarkIds.push_back( entry.first );
        }
        const std::shared_ptr< EstimatableParameterSet< double > > subsetParameters = createParametersToEstimate< double, double >(
                landmarkParameterSettings( scenario.conversionResult_, nullptr, { landmarkIds.at( 0 ) } ), scenario.bodies_ );
        const std::map< std::string, Eigen::Vector3d > subsetUncertainties =
                landmarkLocalFrameUncertainties( scenario.conversionResult_, subsetParameters, Eigen::MatrixXd::Identity( 3, 3 ) );
        BOOST_REQUIRE_EQUAL( subsetUncertainties.size( ), 1 );
        BOOST_CHECK_EQUAL( subsetUncertainties.begin( )->first, landmarkIds.at( 0 ) );
    }

    // A covariance of the wrong size must be rejected rather than silently indexed out of range.
    BOOST_CHECK_THROW(
            landmarkLocalFrameUncertainties( scenario.conversionResult_, parametersToEstimate, Eigen::MatrixXd::Identity( 2, 2 ) ),
            std::runtime_error );
}

//! The real archived LMK local frames must pass validation at the tolerance the code uses.
//!
//! This pins the tolerance: the SUM camera tolerance of 1e-7 would reject most real landmark files, whose
//! worst orthonormality error is 2.3e-7, so reusing it here would be a silent data-rejection regression.
BOOST_AUTO_TEST_CASE( testRealLandmarkLocalFramesPassValidation )
{
    const std::string dataPath = paths::getTudatTestDataPath( ) + "/sum_lmk/real_67p";
    std::vector< std::string > lmkFiles;
    for( const std::filesystem::directory_entry& entry : std::filesystem::directory_iterator( dataPath ) )
    {
        if( entry.path( ).extension( ).string( ) == ".LMK" )
        {
            lmkFiles.push_back( entry.path( ).string( ) );
        }
    }
    BOOST_REQUIRE_GT( lmkFiles.size( ), 10u );

    const std::map< std::string, input_output::sum_lmk::LmkLandmarkData > landmarks = input_output::sum_lmk::readLmkFiles( lmkFiles );
    double worstOrthogonalityError = 0.0;
    for( const auto& landmarkEntry : landmarks )
    {
        Eigen::Matrix3d localFrame;
        localFrame.col( 0 ) = landmarkEntry.second.localXAxis_;
        localFrame.col( 1 ) = landmarkEntry.second.localYAxis_;
        localFrame.col( 2 ) = landmarkEntry.second.localZAxis_;
        worstOrthogonalityError =
                std::max( worstOrthogonalityError, ( localFrame * localFrame.transpose( ) - Eigen::Matrix3d::Identity( ) ).norm( ) );

        // Must not throw, and must come back as an exact rotation.
        Eigen::Matrix3d validatedFrame;
        BOOST_REQUIRE_NO_THROW( validatedFrame =
                                        input_output::sum_lmk::getLmkLandmarkRotationFromLocalToBodyFixedFrame( landmarkEntry.second ) );
        BOOST_CHECK_SMALL( ( validatedFrame * validatedFrame.transpose( ) - Eigen::Matrix3d::Identity( ) ).norm( ), 1.0E-14 );
        BOOST_CHECK_CLOSE( validatedFrame.determinant( ), 1.0, 1.0E-10 );

        // UZ, the maplet plane normal, points outward for every archived landmark.
        BOOST_CHECK_GT( landmarkEntry.second.localZAxis_.dot( landmarkEntry.second.bodyFixedPosition_.normalized( ) ), 0.0 );
    }

    BOOST_TEST_MESSAGE( "Worst LMK local-frame orthonormality error: " << worstOrthogonalityError );
    BOOST_CHECK_GT( worstOrthogonalityError, input_output::sum_lmk::sumCameraRotationMatrixTolerance );
    BOOST_CHECK_LT( worstOrthogonalityError, input_output::sum_lmk::lmkLandmarkFrameTolerance );

    // A materially non-rotational frame is still rejected.
    input_output::sum_lmk::LmkLandmarkData brokenLandmark = landmarks.begin( )->second;
    brokenLandmark.localZAxis_ = brokenLandmark.localXAxis_;
    BOOST_CHECK_THROW( input_output::sum_lmk::getLmkLandmarkRotationFromLocalToBodyFixedFrame( brokenLandmark ), std::runtime_error );
}

//! A landmark seen from several arcs is ONE global parameter, fed by the observations of every arc.
//!
//! This is the requirement the feature exists for. With the spacecraft orbiting close to the body, the same
//! landmark reappears in images belonging to different arcs. Because ground_station_position is an ordinary
//! (non-arc-wise) parameter, it is not duplicated per arc: the estimated vector holds one arc-wise initial
//! state per arc plus a single position per landmark, and every arc's observations contribute to that one
//! position - which is what ties the arcs together.
//!
//! Three things are checked: the parameter-vector layout (one block per landmark, not one per arc), that the
//! design matrix couples each landmark to observations in BOTH arcs, and that a joint solve recovers both
//! the per-arc states and the landmark positions.
BOOST_AUTO_TEST_CASE( testGlobalLandmarkAcrossMultipleArcs )
{
    spice_interface::loadStandardSpiceKernels( );

    const int numberOfImagesPerArc = 6;
    const int epochSpacingSeconds = 250;
    const int arcGapSeconds = 3000;

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
    // Multi-arc propagated bodies carry a MultiArcEphemeris, with the central body as its origin.
    bodies.at( "Spacecraft" )
            ->setEphemeris(
                    std::make_shared< MultiArcEphemeris >( std::map< double, std::shared_ptr< Ephemeris > >( ), "Target", "J2000" ) );
    bodies.processBodyFrameDefinitions< double, double >( );

    Eigen::Vector6d truthKeplerianElements = Eigen::Vector6d::Zero( );
    truthKeplerianElements( semiMajorAxisIndex ) = 1.0E4;
    truthKeplerianElements( eccentricityIndex ) = 0.05;
    truthKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 30.0 );
    truthKeplerianElements( argumentOfPeriapsisIndex ) = unit_conversions::convertDegreesToRadians( 40.0 );
    truthKeplerianElements( longitudeOfAscendingNodeIndex ) = unit_conversions::convertDegreesToRadians( 25.0 );
    truthKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 10.0 );

    const std::map< std::string, input_output::sum_lmk::LmkLandmarkData > landmarks = makeOrbitLandmarks( );
    std::vector< std::string > landmarkIds;
    for( const auto& entry : landmarks )
    {
        landmarkIds.push_back( entry.first );
    }

    // Two groups of images separated by a gap, becoming two propagation arcs. The underlying trajectory is
    // one continuous Keplerian orbit, so the two arc initial states are mutually consistent.
    std::vector< input_output::sum_lmk::SumImageData > images;
    std::vector< double > epochs;
    std::vector< int > arcIndexPerImage;
    for( int arcIndex = 0; arcIndex < 2; ++arcIndex )
    {
        for( int imageIndex = 0; imageIndex < numberOfImagesPerArc; ++imageIndex )
        {
            const int secondsOffset = arcIndex * arcGapSeconds + imageIndex * epochSpacingSeconds;
            input_output::sum_lmk::SumImageData image;
            image.imageId_ = "MARC" + std::to_string( arcIndex ) + "IMG" + std::to_string( imageIndex );
            image.utcEpochString_ = makeUtcString( secondsOffset );
            image.imageSize_ = Eigen::Vector2i( 1024, 1024 );
            image.focalLengthMm_ = 100.0;
            image.opticalCenter_ = Eigen::Vector2d( 512.0, 512.0 );
            image.kMatrix_ << 10.0, 0.0, 0.0, 0.0, 10.0, 0.0;

            const double epoch = sumImageTdbTimeForTest( image );
            if( epochs.empty( ) )
            {
                epochs.push_back( epoch );
            }
            const double timeSinceStart = epoch - epochs.front( );
            const Eigen::Vector6d inertialState = convertKeplerianToCartesianElements(
                    propagateKeplerOrbit< double >( truthKeplerianElements, timeSinceStart, targetGravitationalParameter ),
                    targetGravitationalParameter );
            const Eigen::Vector3d spacecraftBodyFixedPosition =
                    bodies.at( "Target" )->getRotationalEphemeris( )->getRotationToTargetFrame( epoch ).toRotationMatrix( ) *
                    inertialState.head( 3 );
            image.spacecraftObjectVector_ = -spacecraftBodyFixedPosition;
            image.cameraAxes_ = boresightCameraAxes( spacecraftBodyFixedPosition );

            for( const std::string& landmarkId : landmarkIds )
            {
                input_output::sum_lmk::SumLandmarkObservation observation;
                observation.landmarkId_ = landmarkId;
                observation.pixelCoordinates_ = Eigen::Vector2d::Zero( );
                image.landmarkObservations_.push_back( observation );
            }
            images.push_back( image );
            arcIndexPerImage.push_back( arcIndex );
            if( epochs.size( ) < images.size( ) )
            {
                epochs.push_back( epoch );
            }
        }
    }

    LoadedPixelTestData conversionResult = loadPixelTestData( images, landmarks, bodies );

    // Two propagation arcs, each starting at its own first image and ending after its last.
    SelectedAccelerationMap accelerationSettingsMap;
    accelerationSettingsMap[ "Spacecraft" ][ "Target" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    const std::vector< std::string > bodiesToIntegrate = { "Spacecraft" };
    const std::vector< std::string > centralBodies = { "Target" };
    const AccelerationMap accelerationModelMap =
            createAccelerationModelsMap( bodies, accelerationSettingsMap, bodiesToIntegrate, centralBodies );

    std::vector< double > arcStartTimes;
    Eigen::VectorXd concatenatedArcInitialStates = Eigen::VectorXd::Zero( 12 );
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double, double > > > arcPropagatorSettingsList;
    for( int arcIndex = 0; arcIndex < 2; ++arcIndex )
    {
        const double arcFirstEpoch = epochs.at( arcIndex * numberOfImagesPerArc );
        const double arcLastEpoch = epochs.at( arcIndex * numberOfImagesPerArc + numberOfImagesPerArc - 1 );
        const Eigen::Vector6d arcInitialState = convertKeplerianToCartesianElements(
                propagateKeplerOrbit< double >( truthKeplerianElements, arcFirstEpoch - epochs.front( ), targetGravitationalParameter ),
                targetGravitationalParameter );
        arcStartTimes.push_back( arcFirstEpoch - 10.0 );
        concatenatedArcInitialStates.segment( 6 * arcIndex, 6 ) = arcInitialState;

        arcPropagatorSettingsList.push_back( translationalStatePropagatorSettings< double, double >(
                centralBodies,
                accelerationModelMap,
                bodiesToIntegrate,
                arcInitialState,
                arcFirstEpoch - 10.0,
                rungeKuttaFixedStepSettings< double >( 10.0, CoefficientSets::rungeKuttaFehlberg78 ),
                propagationTimeTerminationSettings( arcLastEpoch + 100.0 ) ) );
    }
    const std::shared_ptr< MultiArcPropagatorSettings< double, double > > multiArcPropagatorSettings =
            std::make_shared< MultiArcPropagatorSettings< double, double > >( arcPropagatorSettingsList );

    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( std::make_shared< ArcWiseInitialTranslationalStateEstimatableParameterSettings< double > >(
            "Spacecraft", concatenatedArcInitialStates, arcStartTimes, "Target" ) );
    const std::vector< std::shared_ptr< EstimatableParameterSettings > > landmarkSettings = landmarkParameterSettings( conversionResult );
    BOOST_REQUIRE_EQUAL( landmarkSettings.size( ), landmarkIds.size( ) );
    for( const std::shared_ptr< EstimatableParameterSettings >& landmarkSetting : landmarkSettings )
    {
        parameterNames.push_back( landmarkSetting );
    }

    const std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double, double >( parameterNames, bodies, multiArcPropagatorSettings );

    // THE LAYOUT CLAIM: 6 per arc for the two arc states, plus 3 per landmark ONCE - not once per arc.
    BOOST_REQUIRE_EQUAL( parametersToEstimate->getEstimatedParameterSetSize( ), static_cast< int >( 12 + 3 * landmarkIds.size( ) ) );
    for( const std::string& landmarkId : landmarkIds )
    {
        const std::vector< std::pair< int, int > > indices = parametersToEstimate->getIndicesForParameterType(
                EstimatebleParameterIdentifier( ground_station_position, std::make_pair( "Target", landmarkId ) ) );
        BOOST_REQUIRE_EQUAL( indices.size( ), 1 );
        BOOST_CHECK_EQUAL( indices.at( 0 ).second, 3 );
    }

    OrbitDeterminationManager< double, double > orbitDeterminationManager(
            bodies, parametersToEstimate, conversionResult.observationModelSettings_, multiArcPropagatorSettings );

    const Eigen::VectorXd truthParameters = parametersToEstimate->getFullParameterValues< double >( );
    const std::shared_ptr< ObservationCollection< double, double > > simulatedObservations =
            simulateIdealPixelObservations( conversionResult.observationCollection_, orbitDeterminationManager, bodies );

    // Perturb both arc states and every landmark, then recover all of them together.
    Eigen::VectorXd perturbedParameters = truthParameters;
    for( int arcIndex = 0; arcIndex < 2; ++arcIndex )
    {
        perturbedParameters.segment( 6 * arcIndex, 3 ) += Eigen::Vector3d::Constant( 5.0 );
        perturbedParameters.segment( 6 * arcIndex + 3, 3 ) += Eigen::Vector3d::Constant( 0.005 );
    }
    for( const std::string& landmarkId : landmarkIds )
    {
        const int startIndex = parametersToEstimate
                                       ->getIndicesForParameterType( EstimatebleParameterIdentifier(
                                               ground_station_position, std::make_pair( "Target", landmarkId ) ) )
                                       .at( 0 )
                                       .first;
        perturbedParameters.segment( startIndex, 3 ) += ( Eigen::Vector3d( ) << 2.0, -3.0, 1.5 ).finished( );
    }
    parametersToEstimate->resetParameterValues( perturbedParameters );

    const std::shared_ptr< EstimationInput< double, double > > estimationInput =
            std::make_shared< EstimationInput< double, double > >( simulatedObservations );
    estimationInput->defineEstimationSettings( true, true, true, true, true, false );
    estimationInput->setConvergenceChecker( std::make_shared< EstimationConvergenceChecker >( 12 ) );

    const std::shared_ptr< EstimationOutput< double > > estimationOutput = orbitDeterminationManager.estimateParameters( estimationInput );
    BOOST_REQUIRE( !estimationOutput->exceptionDuringInversion_ );
    BOOST_REQUIRE( !estimationOutput->exceptionDuringPropagation_ );
    BOOST_CHECK_SMALL( residualRms( estimationOutput->residuals_ ), 1.0E-5 );

    // THE COUPLING CLAIM: each landmark's design-matrix columns must be non-zero in rows belonging to BOTH
    // arcs. If the parameter were somehow arc-local, or only the first arc reached it, one of the two
    // contributions below would vanish.
    const Eigen::MatrixXd designMatrix = estimationOutput->getUnnormalizedDesignMatrix( );
    BOOST_REQUIRE_EQUAL( designMatrix.cols( ), parametersToEstimate->getEstimatedParameterSetSize( ) );
    const std::vector< double > observationTimes = simulatedObservations->getConcatenatedTimeVector( );
    BOOST_REQUIRE_EQUAL( static_cast< int >( observationTimes.size( ) ), designMatrix.rows( ) );
    const double arcBoundaryTime = 0.5 * ( epochs.at( numberOfImagesPerArc - 1 ) + epochs.at( numberOfImagesPerArc ) );

    for( const std::string& landmarkId : landmarkIds )
    {
        const int startIndex = parametersToEstimate
                                       ->getIndicesForParameterType( EstimatebleParameterIdentifier(
                                               ground_station_position, std::make_pair( "Target", landmarkId ) ) )
                                       .at( 0 )
                                       .first;
        double firstArcContribution = 0.0;
        double secondArcContribution = 0.0;
        for( int row = 0; row < designMatrix.rows( ); ++row )
        {
            const double rowNorm = designMatrix.block( row, startIndex, 1, 3 ).norm( );
            if( observationTimes.at( row ) < arcBoundaryTime )
            {
                firstArcContribution += rowNorm;
            }
            else
            {
                secondArcContribution += rowNorm;
            }
        }
        BOOST_CHECK_GT( firstArcContribution, 0.0 );
        BOOST_CHECK_GT( secondArcContribution, 0.0 );
    }

    // THE RECOVERY CLAIM: arc states and landmark positions come back together.
    const Eigen::VectorXd estimationError = estimationOutput->parameterEstimate_ - truthParameters;
    double worstLandmarkError = 0.0;
    for( const std::string& landmarkId : landmarkIds )
    {
        const int startIndex = parametersToEstimate
                                       ->getIndicesForParameterType( EstimatebleParameterIdentifier(
                                               ground_station_position, std::make_pair( "Target", landmarkId ) ) )
                                       .at( 0 )
                                       .first;
        worstLandmarkError = std::max( worstLandmarkError, estimationError.segment( startIndex, 3 ).norm( ) );
    }
    BOOST_TEST_MESSAGE( "Multi-arc joint solve: arc 0 position error " << estimationError.segment( 0, 3 ).norm( ) << " m, arc 1 "
                                                                       << estimationError.segment( 6, 3 ).norm( )
                                                                       << " m, worst landmark error " << worstLandmarkError << " m" );
    BOOST_CHECK_SMALL( estimationError.segment( 0, 3 ).norm( ), 1.0E-3 );
    BOOST_CHECK_SMALL( estimationError.segment( 6, 3 ).norm( ), 1.0E-3 );
    BOOST_CHECK_SMALL( worstLandmarkError, 1.0E-3 );
}

//! Real-data behaviour on the committed Rosetta/67P arc, at 164 km range.
//!
//! What this measures, and it is not what a naive information count predicts. Comparing SIGMA_LMK (about
//! 0.4 m per axis, i.e. about 0.13 px at this range) against the data's information about a single landmark
//! suggests the prior should dominate roughly ten to one and the landmarks should barely move. They move
//! about 1 m in the median and about 5.6 m at worst - several times the prior sigma.
//!
//! The reason is the network/orbit degeneracy: a bulk shift of the whole landmark network is nearly
//! indistinguishable from an error in the spacecraft orbit, and along that direction the data has real
//! leverage while the prior only penalises each landmark separately. Remove the a-priori entirely and the
//! network translates almost rigidly by about 54 m, which is what the second solve below shows. This is the
//! same trap as the "tempting calculation that is wrong" in the SIGMA_PTG discussion: information computed
//! with the orbit held fixed does not decide a joint solve.
//!
//! The practical consequence, which the range table in camera-pointing-observability.md section 9 now
//! records: on a long-range arc, estimating landmark positions jointly with the state is NOT harmless. The
//! landmarks absorb orbit error at the metre level even with SIGMA_LMK applied, and the normal matrix is
//! close to numerically rank deficient (condition number ~1e16). The partial itself is sound - verified by
//! finite differences in unitTestPixelCoordinatesPartials and by the synthetic recovery tests above.
BOOST_AUTO_TEST_CASE( testRealArcLandmarkRegression )
{
    spice_interface::loadStandardSpiceKernels( );

    RealRosettaScenario scenario = buildRealRosettaScenario( );

    const std::vector< std::string > landmarkIds = observedLandmarkIds( scenario.conversionResult_, "Comet" );

    // The real LMK files all carry SIGMA_LMK, so the conversion produces a landmark a-priori for each.
    std::size_t numberOfLandmarkAprioriEntries = 0;
    for( const auto& entry : landmarkAndPointingAprioriEntries( scenario.conversionResult_ ) )
    {
        if( entry.first.first == ground_station_position )
        {
            ++numberOfLandmarkAprioriEntries;
        }
    }
    BOOST_REQUIRE_GT( numberOfLandmarkAprioriEntries, 10u );
    BOOST_CHECK_EQUAL( numberOfLandmarkAprioriEntries, landmarkIds.size( ) );

    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
            "Spacecraft", scenario.initialStateGuess_, "Comet" ) );
    const std::vector< std::shared_ptr< EstimatableParameterSettings > > landmarkSettings =
            landmarkParameterSettings( scenario.conversionResult_, nullptr, {}, "Comet" );
    BOOST_REQUIRE_EQUAL( landmarkSettings.size( ), landmarkIds.size( ) );
    for( const std::shared_ptr< EstimatableParameterSettings >& landmarkSetting : landmarkSettings )
    {
        parameterNames.push_back( landmarkSetting );
    }

    const std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double, double >( parameterNames, scenario.bodies_, scenario.propagatorSettings_ );
    BOOST_REQUIRE_EQUAL( parametersToEstimate->getEstimatedParameterSetSize( ), static_cast< int >( 6 + 3 * landmarkSettings.size( ) ) );

    OrbitDeterminationManager< double, double > orbitDeterminationManager(
            scenario.bodies_, parametersToEstimate, scenario.conversionResult_.observationModelSettings_, scenario.propagatorSettings_ );

    const Eigen::VectorXd aprioriParameters = parametersToEstimate->getFullParameterValues< double >( );

    std::map< std::shared_ptr< ObservationCollectionParser >, double > weightsPerObservationParser;
    weightsPerObservationParser[ observationParser( pixel_coordinates ) ] = 1.0;
    scenario.conversionResult_.observationCollection_->setConstantWeightPerObservable( weightsPerObservationParser );

    // Solve twice: with the SIGMA_LMK a-priori, and without it. Comparing the two is what actually shows
    // the landmark a-priori reaches the normal equations at the right parameter indices - a magnitude bound
    // on its own would pass just as well if the a-priori were silently inert.
    auto runSolve = [ & ]( const bool useLandmarkApriori ) {
        parametersToEstimate->resetParameterValues( aprioriParameters );
        const std::shared_ptr< EstimationInput< double, double > > estimationInput = useLandmarkApriori
                ? std::make_shared< EstimationInput< double, double > >(
                          scenario.conversionResult_.observationCollection_,
                          landmarkInverseAprioriCovariance( scenario.conversionResult_, parametersToEstimate ) )
                : std::make_shared< EstimationInput< double, double > >( scenario.conversionResult_.observationCollection_ );
        estimationInput->defineEstimationSettings( true, true, true, false, true, false );
        estimationInput->setConvergenceChecker( std::make_shared< EstimationConvergenceChecker >( 6 ) );
        return orbitDeterminationManager.estimateParameters( estimationInput );
    };

    // Landmark displacements from their LMK values, as a sorted list so the distribution can be reported
    // rather than only its maximum.
    auto landmarkDisplacements = [ & ]( const std::shared_ptr< EstimationOutput< double > >& output ) {
        std::vector< double > displacements;
        for( const std::string& landmarkId : landmarkIds )
        {
            const std::vector< std::pair< int, int > > indices = parametersToEstimate->getIndicesForParameterType(
                    EstimatebleParameterIdentifier( ground_station_position, std::make_pair( "Comet", landmarkId ) ) );
            BOOST_REQUIRE_EQUAL( indices.size( ), 1 );
            displacements.push_back( ( output->parameterEstimate_.segment( indices.at( 0 ).first, 3 ) -
                                       aprioriParameters.segment( indices.at( 0 ).first, 3 ) )
                                             .norm( ) );
        }
        std::sort( displacements.begin( ), displacements.end( ) );
        return displacements;
    };

    const std::shared_ptr< EstimationOutput< double > > estimationOutput = runSolve( true );
    BOOST_REQUIRE( !estimationOutput->exceptionDuringInversion_ );
    BOOST_REQUIRE( !estimationOutput->exceptionDuringPropagation_ );

    const Eigen::MatrixXd residualHistory = estimationOutput->getResidualHistoryMatrix( );
    const double initialRms = residualRms( residualHistory.col( 0 ) );
    const double finalRms = residualRms( estimationOutput->residuals_ );
    BOOST_CHECK_LT( finalRms, initialRms );

    const std::vector< double > constrainedDisplacements = landmarkDisplacements( estimationOutput );
    const double medianConstrained = constrainedDisplacements.at( constrainedDisplacements.size( ) / 2 );
    const double worstConstrained = constrainedDisplacements.back( );

    const std::shared_ptr< EstimationOutput< double > > unconstrainedOutput = runSolve( false );
    BOOST_REQUIRE( !unconstrainedOutput->exceptionDuringInversion_ );
    const std::vector< double > unconstrainedDisplacements = landmarkDisplacements( unconstrainedOutput );
    const double medianUnconstrained = unconstrainedDisplacements.at( unconstrainedDisplacements.size( ) / 2 );
    const double worstUnconstrained = unconstrainedDisplacements.back( );

    BOOST_TEST_MESSAGE( "Real 164 km arc: pixel RMS " << initialRms << " -> " << finalRms << " px" );
    BOOST_TEST_MESSAGE( "  landmark displacement with SIGMA_LMK a-priori:    median " << medianConstrained << " m, worst "
                                                                                      << worstConstrained << " m" );
    BOOST_TEST_MESSAGE( "  landmark displacement without any a-priori:       median " << medianUnconstrained << " m, worst "
                                                                                      << worstUnconstrained << " m" );

    // THE A-PRIORI CLAIM: removing it must let the landmarks move materially further. If the a-priori were
    // assembled at the wrong indices, or skipped, these two would be identical. Measured: the a-priori cuts
    // the median displacement by roughly a factor of 50 (about 1 m against about 54 m).
    BOOST_CHECK_LT( medianConstrained, medianUnconstrained );
    BOOST_CHECK_LT( worstConstrained, worstUnconstrained );
    BOOST_CHECK_GT( medianUnconstrained, 10.0 * medianConstrained );

    // Without the a-priori the whole network shifts almost rigidly - median and worst displacement are
    // nearly equal - which is the network/orbit degeneracy in its purest form: a bulk translation of the
    // landmarks is nearly indistinguishable from an error in the spacecraft orbit.
    BOOST_CHECK_GT( medianUnconstrained, 0.8 * worstUnconstrained );

    // Regression guard on the constrained solve. NOTE the size of these numbers: SIGMA_LMK is about 0.4 m
    // per axis, yet the landmarks still move about 1 m in the median and about 5.6 m at worst. That is NOT
    // a broken a-priori - it is the joint solve buying a better orbit fit by paying the landmark prior, and
    // it is the same trap as the "tempting calculation that is wrong" in the SIGMA_PTG discussion: the
    // information the data has about a landmark *with the orbit held fixed* is not what decides a joint
    // solve. See the range discussion in camera-pointing-observability.md section 9.
    BOOST_CHECK_LT( medianConstrained, 5.0 );
    BOOST_CHECK_LT( worstConstrained, 15.0 );

    // The local-frame reporting must work on real archived frames, not only synthetic ones.
    const std::map< std::string, Eigen::Vector3d > localFrameUncertainties = landmarkLocalFrameUncertainties(
            scenario.conversionResult_, parametersToEstimate, estimationOutput->getUnnormalizedCovarianceMatrix( ) );
    BOOST_CHECK_EQUAL( localFrameUncertainties.size( ), landmarkIds.size( ) );
    for( const auto& entry : localFrameUncertainties )
    {
        BOOST_CHECK_GT( entry.second.minCoeff( ), 0.0 );
        // Cannot be worse than the prior it started from.
        BOOST_CHECK_LT( entry.second.maxCoeff( ), 1.0 );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
