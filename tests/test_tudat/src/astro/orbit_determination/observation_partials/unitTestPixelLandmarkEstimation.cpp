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

    const std::string dataPath = paths::getTudatTestDataPath( ) + "/sum_lmk/real_67p";
    spice_interface::loadSpiceKernelInTudat( dataPath + "/rosetta_67p_orbiter_arc.bsp" );

    // --- Bodies: comet 67P (point-mass gravity, body-fixed frame from the committed attitude) + s/c.
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

    // --- Gather the committed reduced SUM/LMK dataset.
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
    BOOST_REQUIRE_GE( sumFiles.size( ), 2u );
    BOOST_REQUIRE( !lmkFiles.empty( ) );

    // --- Image epochs (the conversion derives observation times from the SUM UTC strings).
    std::vector< input_output::sum_lmk::SumImageData > sumImages = input_output::sum_lmk::readSumFiles( sumFiles );
    std::sort( sumImages.begin( ), sumImages.end( ), []( const auto& a, const auto& b ) {
        return observation_models::detail::convertSumUtcStringToSecondsSinceJ2000< double >( a ) <
                observation_models::detail::convertSumUtcStringToSecondsSinceJ2000< double >( b );
    } );
    const double epoch0 = observation_models::detail::convertSumUtcStringToSecondsSinceJ2000< double >( sumImages.front( ) );
    const double epochLast = observation_models::detail::convertSumUtcStringToSecondsSinceJ2000< double >( sumImages.back( ) );

    // --- A-priori spacecraft state = reconstructed Rosetta orbiter state relative to the comet, from
    //     the committed SPK (NAIF ids: Rosetta orbiter -226, comet 1000012), in metres.
    const Eigen::Vector6d initialStateGuess = spice_interface::getBodyCartesianStateAtEpoch( "-226", "1000012", "J2000", "none", epoch0 );

    // --- Ingest the real pixel-landmark observations (registers cameras + landmarks on the bodies).
    SumLmkObservationConversionSettings conversionSettings( "Comet", "Spacecraft" );
    SumLmkObservationConversionResult< double, double > conversionResult =
            createSumLmkObservationCollection< double, double >( sumFiles, lmkFiles, bodies, conversionSettings );
    BOOST_REQUIRE( conversionResult.observationCollection_ != nullptr );
    BOOST_REQUIRE_GT( conversionResult.observationCollection_->getTotalObservableSize( ), 0 );

    // --- Dynamics: spacecraft under comet point-mass gravity only.
    SelectedAccelerationMap accelerationSettingsMap;
    accelerationSettingsMap[ "Spacecraft" ][ "Comet" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    const std::vector< std::string > bodiesToIntegrate = { "Spacecraft" };
    const std::vector< std::string > centralBodies = { "Comet" };
    const AccelerationMap accelerationModelMap =
            createAccelerationModelsMap( bodies, accelerationSettingsMap, bodiesToIntegrate, centralBodies );

    const std::shared_ptr< IntegratorSettings< double > > integratorSettings =
            rungeKuttaFixedStepSettings< double >( 30.0, CoefficientSets::rungeKuttaFehlberg78 );
    const std::shared_ptr< TranslationalStatePropagatorSettings< double, double > > propagatorSettings =
            translationalStatePropagatorSettings< double, double >( centralBodies,
                                                                    accelerationModelMap,
                                                                    bodiesToIntegrate,
                                                                    initialStateGuess,
                                                                    epoch0,
                                                                    integratorSettings,
                                                                    propagationTimeTerminationSettings( epochLast + 600.0 ) );

    // --- Estimate the spacecraft initial state from the real pixel observations.
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
            "Spacecraft", initialStateGuess, "Comet" ) );
    const std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double, double >( parameterNames, bodies, propagatorSettings );

    OrbitDeterminationManager< double, double > orbitDeterminationManager(
            bodies, parametersToEstimate, conversionResult.observationModelSettings_, propagatorSettings );

    std::map< std::shared_ptr< ObservationCollectionParser >, double > weightsPerObservationParser;
    weightsPerObservationParser[ observationParser( pixel_coordinates ) ] = 1.0;
    conversionResult.observationCollection_->setConstantWeightPerObservable( weightsPerObservationParser );

    const std::shared_ptr< EstimationInput< double, double > > estimationInput =
            std::make_shared< EstimationInput< double, double > >( conversionResult.observationCollection_ );
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
    BOOST_CHECK_LT( ( estimatedPosition - initialStateGuess.segment( 0, 3 ) ).norm( ), 5.0E3 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
