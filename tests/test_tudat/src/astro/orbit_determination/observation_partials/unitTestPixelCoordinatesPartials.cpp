/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"

#include "tudat/io/basicInputOutput.h"
#include "tudat/interface/spice/spiceInterface.h"

#include "tudat/simulation/estimation_setup/createObservationModelFactory.h"
#include "tudat/astro/observation_models/oneWayRangeObservationModel.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/constantRotationRate.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/initialTranslationalState.h"
#include "tudat/simulation/estimation_setup/createObservationPartials.h"
#include "tudat/support/numericalObservationPartial.h"
#include "tudat/simulation/environment_setup/createCameras.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"

#include "tudat/support/observationPartialTestFunctions.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::gravitation;
using namespace tudat::ephemerides;
using namespace tudat::observation_models;
using namespace tudat::simulation_setup;
using namespace tudat::spice_interface;
using namespace tudat::observation_partials;
using namespace tudat::estimatable_parameters;

BOOST_AUTO_TEST_SUITE( test_pixel_coordinates_partials )

//! Test partial derivatives of pixel coordinate observable, using general test suite of observation partials.
BOOST_AUTO_TEST_CASE( testPixelCoordinatesPartials )
{
    // Define and create ground stations.
    std::vector< std::pair< std::string, std::string > > groundStations;
    groundStations.resize( 2 );
    groundStations[ 0 ] = std::make_pair( "Earth", "Graz" );
    groundStations[ 1 ] = std::make_pair( "Mars", "MSL" );

    Eigen::VectorXd parameterPerturbationMultipliers = Eigen::VectorXd::Constant( 4, 1.0 );
    parameterPerturbationMultipliers( 2 ) = 10.0;

    // Test partials with constant ephemerides (allows test of position partials).
    {
        const double stateEvaluationTime = 1.1E7;

        // Create environment.
        SystemOfBodies bodies = setupEnvironment( groundStations, 1.0E7, 1.2E7, stateEvaluationTime, true );

        // Create a camera on Earth with non-identity mounting angles, pointed along the body-fixed Earth-to-Mars line of sight.
        const std::string cameraName = "Camera";
        const Eigen::Vector3d cameraBodyFixedPosition = ( Eigen::Vector3d( ) << 1.7E6, -6.2E6, 1.3E5 ).finished( );
        const Eigen::Vector3d marsStationPositionInInertialFrame =
                bodies.at( "Mars" )->getStateInBaseFrameFromEphemeris< double, double >( stateEvaluationTime ).segment( 0, 3 ) +
                bodies.at( "Mars" )->getRotationalEphemeris( )->getRotationToTargetFrame( stateEvaluationTime ).inverse( ) *
                        bodies.at( "Mars" )
                                ->getGroundStation( groundStations[ 1 ].second )
                                ->getNominalStationState( )
                                ->getNominalCartesianPosition( );
        const Eigen::Vector3d relativePositionInInertialFrame = marsStationPositionInInertialFrame -
                bodies.at( "Earth" )->getStateInBaseFrameFromEphemeris< double, double >( stateEvaluationTime ).segment( 0, 3 );
        const Eigen::Vector3d relativePositionInBodyFixedFrame =
                bodies.at( "Earth" )->getRotationalEphemeris( )->getRotationToTargetFrame( stateEvaluationTime ) *
                        relativePositionInInertialFrame -
                cameraBodyFixedPosition;
        const Eigen::Vector3d cameraBoresightDirection = relativePositionInBodyFixedFrame.normalized( );
        const Eigen::Vector3d cameraBoresightEulerAngles =
                ( Eigen::Vector3d( ) << std::atan2( cameraBoresightDirection.y( ), cameraBoresightDirection.x( ) ),
                  std::asin( cameraBoresightDirection.z( ) ),
                  0.3 )
                        .finished( );

        createCamera( bodies.at( "Earth" ),
                      cameraName,
                      cameraBoresightEulerAngles,
                      std::make_pair( 1200.0, 900.0 ),
                      std::make_pair( 320.0, 240.0 ),
                      cameraBodyFixedPosition );

        // Set link ends for observation model.
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = groundStations[ 1 ];
        linkEnds[ receiver ] = std::make_pair( std::string( "Earth" ), cameraName );

        // Generate pixel coordinates model.
        std::vector< std::string > perturbingBodies;
        perturbingBodies.push_back( "Earth" );
        std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionSettings = {
            std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( perturbingBodies )
        };
        std::shared_ptr< ObservationModel< 2 > > pixelCoordinatesModel =
                observation_models::ObservationModelCreator< 2, double, double >::createObservationModel(
                        pixelCoordinatesSettings( linkEnds, lightTimeCorrectionSettings ), bodies );

        // Create parameter objects.
        std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > > initialStateParameters;
        initialStateParameters.push_back( std::make_shared< InitialTranslationalStateParameter< double > >(
                "Earth",
                bodies.at( "Earth" )->getStateInBaseFrameFromEphemeris< double, double >( stateEvaluationTime ),
                "SSB",
                "ECLIPJ2000" ) );
        initialStateParameters.push_back( std::make_shared< InitialTranslationalStateParameter< double > >(
                "Mars",
                bodies.at( "Mars" )->getStateInBaseFrameFromEphemeris< double, double >( stateEvaluationTime ),
                "SSB",
                "ECLIPJ2000" ) );

        std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet =
                std::make_shared< EstimatableParameterSet< double > >(
                        std::vector< std::shared_ptr< EstimatableParameter< double > > >( ),
                        std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > >( ),
                        initialStateParameters );

        testObservationPartials( pixelCoordinatesModel,
                                 bodies,
                                 fullEstimatableParameterSet,
                                 linkEnds,
                                 pixel_coordinates,
                                 1.0E-5,
                                 true,
                                 false,
                                 1.0,
                                 parameterPerturbationMultipliers,
                                 nullptr,
                                 stateEvaluationTime );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
