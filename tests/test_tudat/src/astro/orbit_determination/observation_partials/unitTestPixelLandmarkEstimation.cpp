/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

// Integration smoke tests for per-image camera-pointing estimation from pixel-landmark
// observations. The synthetic scene uses the (+Z boresight) geometry the projection model
// requires; see unitTestPixelCoordinatesPartials / unitTestSumLmkFileReader for the
// characterization tests that document why real SPC SUM/LMK data does not yet flow through.

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <boost/test/unit_test.hpp>

#include "tudat/astro/ephemerides/constantEphemeris.h"
#include "tudat/astro/ephemerides/constantRotationalEphemeris.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/cameraPointingCorrection.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/simulation/estimation_setup/createObservationModelFactory.h"
#include "tudat/simulation/estimation_setup/createObservationPartials.h"
#include "tudat/simulation/estimation_setup/processSumLmkFiles.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::ephemerides;
using namespace tudat::observation_models;
using namespace tudat::simulation_setup;
using namespace tudat::observation_partials;
using namespace tudat::estimatable_parameters;

namespace
{

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

}  // namespace unit_tests
}  // namespace tudat
