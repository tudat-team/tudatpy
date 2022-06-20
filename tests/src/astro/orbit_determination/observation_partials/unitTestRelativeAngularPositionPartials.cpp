/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <limits>
#include <string>
#include <vector>

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>
#include <boost/lambda/lambda.hpp>

#include "tudat/basics/testMacros.h"

#include "tudat/io/basicInputOutput.h"
#include "tudat/interface/spice/spiceInterface.h"

#include "tudat/simulation/estimation_setup/createObservationModel.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/constantRotationRate.h"
#include "tudat/simulation/estimation_setup/createObservationPartials.h"
#include "tudat/support/numericalObservationPartial.h"
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

BOOST_AUTO_TEST_SUITE( test_relative_angular_position_partials)

//! Test partial derivatives of relative angular position observable, using general test suite of observation partials.
BOOST_AUTO_TEST_CASE( testRelativeAngularPositionPartials )
{
    std::cout.precision( 16 );

    spice_interface::loadStandardSpiceKernels( { paths::getSpiceKernelPath( ) + "/de430_mar097_small.bsp" } );

    // Define bodies to use.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Mars" );
    bodiesToCreate.push_back( "Phobos" );

    // Specify initial time
    double initialEphemerisTime = 0.0;
    double finalEphemerisTime = initialEphemerisTime + 7.0 * 86400.0;
    double maximumTimeStep = 3600.0;
    double buffer = 10.0 * maximumTimeStep;

    // Create bodies settings needed in simulation
    BodyListSettings bodySettings =
            getDefaultBodySettings(
                    bodiesToCreate, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    bodySettings.addSettings( "Phobos" );
    bodySettings.at( "Phobos" )->ephemerisSettings = getDefaultEphemerisSettings( "Phobos" );

    // Create bodies
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // Define link ends for observations.
    LinkEnds linkEnds;
    linkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Earth" , ""  );
    linkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Mars" , ""  );
    linkEnds[ transmitter2 ] = std::make_pair< std::string, std::string >( "Phobos" , ""  );

    LinkEnds linkEndsAngularPosition;
    linkEndsAngularPosition[ receiver ] = std::make_pair< std::string, std::string >( "Earth" , ""  );
    linkEndsAngularPosition[ transmitter ] = std::make_pair< std::string, std::string >( "Mars" , ""  );

    // Create light-time correction settings
    std::vector< std::string > lightTimePerturbingBodies = { "Sun" };
    std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionSettings;
    lightTimeCorrectionSettings.push_back( std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
            lightTimePerturbingBodies ) );

    // Create observation settings
    std::shared_ptr< ObservationModelSettings > angularPositionSettings = std::make_shared< ObservationModelSettings >
            ( angular_position, linkEndsAngularPosition, lightTimeCorrectionSettings,
              std::make_shared< ConstantObservationBiasSettings >( ( Eigen::Vector2d( ) << 3.2E-9, -1.5E-8 ).finished( ), true ) );

    std::shared_ptr< ObservationModelSettings > relativeAngularPositionSettings = std::make_shared< ObservationModelSettings >
            ( relative_angular_position, linkEnds, lightTimeCorrectionSettings,
              std::make_shared< ConstantObservationBiasSettings >( ( Eigen::Vector2d( ) << 3.2E-9, -1.5E-8 ).finished( ), true ) );

    // Create observation model.
    std::shared_ptr< ObservationModel< 2, double, double > > angularPositionModel =
            ObservationModelCreator< 2, double, double >::createObservationModel( angularPositionSettings, bodies );

    std::shared_ptr< ObservationModel< 2, double, double > > relativeAngularPositionModel =
            ObservationModelCreator< 2, double, double >::createObservationModel( relativeAngularPositionSettings, bodies );


    // Compute observation separately with two functions.
    double receiverObservationTime = ( finalEphemerisTime + initialEphemerisTime ) / 2.0;
    std::vector< double > linkEndTimesRelativeAngularPosition;
    std::vector< Eigen::Vector6d > linkEndStatesRelativeAngularPosition;
    Eigen::Vector2d relativeAngularObservation = relativeAngularPositionModel->computeObservationsWithLinkEndData(
            receiverObservationTime, receiver, linkEndTimesRelativeAngularPosition, linkEndStatesRelativeAngularPosition );

    Eigen::Vector3d positionDifferenceMarsEarth = ( linkEndStatesRelativeAngularPosition[ 0 ] - linkEndStatesRelativeAngularPosition[ 2 ] ).segment( 0, 3 );
    Eigen::Vector3d sphericalCoordinatesMars = coordinate_conversions::convertCartesianToSpherical( positionDifferenceMarsEarth );
    double rightAscensionMars = sphericalCoordinatesMars.z( );
    double declinationMars = mathematical_constants::PI / 2.0 - sphericalCoordinatesMars.y( );

    Eigen::Vector3d modifiedAbsolutePositionMars = linkEndStatesRelativeAngularPosition[ 0 ].segment( 0, 3 ) + ( Eigen::Vector3d( ) << 100.0, 0.0, 0.0 ).finished( );
    Eigen::Vector3d modifiedRelativePositionMarsEarth = modifiedAbsolutePositionMars - linkEndStatesRelativeAngularPosition[ 2 ].segment( 0, 3 );
    Eigen::Vector3d modifiedSphericalCoordinatesMars = coordinate_conversions::convertCartesianToSpherical( modifiedRelativePositionMarsEarth );
    double modifiedRightAscensionMars = modifiedSphericalCoordinatesMars.z( );
    double modifiedDeclinationMars = mathematical_constants::PI / 2.0 - modifiedSphericalCoordinatesMars.y( );

    double partialRightAscensionWrtX = ( modifiedRightAscensionMars - rightAscensionMars ) / 100.0;
    double partialDeclinationWrtX = ( modifiedDeclinationMars - declinationMars ) / 100.0;
    std::cout << "partialRightAscensionWrtX: " << partialRightAscensionWrtX << "\n\n";
    std::cout << "partialDeclinationWrtX: " << partialDeclinationWrtX << "\n\n";

    modifiedAbsolutePositionMars = linkEndStatesRelativeAngularPosition[ 0 ].segment( 0, 3 ) + ( Eigen::Vector3d( ) << 0.0, 100.0, 0.0 ).finished( );
    modifiedRelativePositionMarsEarth = modifiedAbsolutePositionMars - linkEndStatesRelativeAngularPosition[ 2 ].segment( 0, 3 );
    modifiedSphericalCoordinatesMars = coordinate_conversions::convertCartesianToSpherical( modifiedRelativePositionMarsEarth );
    modifiedRightAscensionMars = modifiedSphericalCoordinatesMars.z( );
    modifiedDeclinationMars = mathematical_constants::PI / 2.0 - modifiedSphericalCoordinatesMars.y( );

    double partialRightAscensionWrtY = ( modifiedRightAscensionMars - rightAscensionMars ) / 100.0;
    double partialDeclinationWrtY = ( modifiedDeclinationMars - declinationMars ) / 100.0;
    std::cout << "partialRightAscensionWrtY: " << partialRightAscensionWrtY << "\n\n";
    std::cout << "partialDeclinationWrtY: " << partialDeclinationWrtY << "\n\n";

    modifiedAbsolutePositionMars = linkEndStatesRelativeAngularPosition[ 0 ].segment( 0, 3 ) + ( Eigen::Vector3d( ) << 0.0, 0.0, 100.0 ).finished( );
    modifiedRelativePositionMarsEarth = modifiedAbsolutePositionMars - linkEndStatesRelativeAngularPosition[ 2 ].segment( 0, 3 );
    modifiedSphericalCoordinatesMars = coordinate_conversions::convertCartesianToSpherical( modifiedRelativePositionMarsEarth );
    modifiedRightAscensionMars = modifiedSphericalCoordinatesMars.z( );
    modifiedDeclinationMars = mathematical_constants::PI / 2.0 - modifiedSphericalCoordinatesMars.y( );

    double partialRightAscensionWrtZ = ( modifiedRightAscensionMars - rightAscensionMars ) / 100.0;
    double partialDeclinationWrtZ = ( modifiedDeclinationMars - declinationMars ) / 100.0;
    std::cout << "partialRightAscensionWrtZ: " << partialRightAscensionWrtZ << "\n\n";
    std::cout << "partialDeclinationWrtZ: " << partialDeclinationWrtZ << "\n\n";

    std::vector< double > linkEndTimesAngularPosition;
    std::vector< Eigen::Vector6d > linkEndStatesAngularPosition;
    Eigen::Vector2d angularObservation = angularPositionModel->computeObservationsWithLinkEndData(
            receiverObservationTime, receiver, linkEndTimesAngularPosition, linkEndStatesAngularPosition );

    Eigen::Matrix< double, 2, 3 > partialsAngularPositionWrtPosition = calculatePartialOfAngularPositionWrtLinkEndPosition(
            - positionDifferenceMarsEarth,  false );
    Eigen::Matrix< double, 2, 3 > partialsAngularPositionWrtPosition2 = calculatePartialOfAngularPositionWrtLinkEndPosition2(
            - positionDifferenceMarsEarth,  false );

    std::cout << "partialsAngularPositionWrtPosition: " << partialsAngularPositionWrtPosition << "\n\n";
    std::cout << "partialsAngularPositionWrtPosition - v2: " << partialsAngularPositionWrtPosition2 << "\n\n";

    // Define and create ground stations.
    std::vector< std::pair< std::string, std::string > > groundStations;
    groundStations.resize( 3 );
    groundStations[ 0 ] = std::make_pair( "Earth", "Graz" );
    groundStations[ 1 ] = std::make_pair( "Mars", "MSL" );
    groundStations[ 2 ] = std::make_pair( "Moon", "" );

    Eigen::VectorXd parameterPerturbationMultipliers = Eigen::VectorXd::Constant( 4, 1.0 );
    parameterPerturbationMultipliers( 2 ) = 10.0;
    // Test partials with constant ephemerides (allows test of position partials)
    {
        // Create environment
        SystemOfBodies bodies = setupEnvironment( groundStations, 1.0E7, 1.2E7, 1.1E7, true );

        // Set link ends for observation model
        LinkEnds linkEnds;
        linkEnds[ receiver ] = groundStations[ 0 ];
        linkEnds[ transmitter ] = groundStations[ 1 ];
        linkEnds[ transmitter2 ] = groundStations[ 2 ];


        // Generate one-way range model
        std::vector< std::string > perturbingBodies;
        perturbingBodies.push_back( "Earth" );
        std::shared_ptr< ObservationModel< 2 > > relativeAngularPositionModel =
                observation_models::ObservationModelCreator< 2, double, double >::createObservationModel(
                    std::make_shared< observation_models::ObservationModelSettings >(
                        observation_models::relative_angular_position, linkEnds, std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                            perturbingBodies ) ), bodies  );

        // Create parameter objects.
        std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet =
                createEstimatableParameters( bodies, 1.1E7 );

        testObservationPartials( relativeAngularPositionModel, bodies, fullEstimatableParameterSet, linkEnds, relative_angular_position, 1.0E-4,
                                 true, true, 1.0, parameterPerturbationMultipliers );
    }


//    // Test partials with real ephemerides (without test of position partials)
//    {
//        std::cout << "Test 1" << std::endl;
//        // Create environment
//        SystemOfBodies bodies = setupEnvironment( groundStations, 1.0E7, 1.2E7, 1.1E7, false );

//        // Set link ends for observation model
//        LinkEnds linkEnds;
//        linkEnds[ transmitter ] = groundStations[ 1 ];
//        linkEnds[ receiver ] = groundStations[ 0 ];

//        // Generate one-way range model
//        std::shared_ptr< ObservationModel< 2 > > angularPositionModel =
//                observation_models::ObservationModelCreator< 2, double, double >::createObservationModel(
//                    linkEnds, std::make_shared< observation_models::ObservationModelSettings >(
//                        observation_models::angular_position ), bodies  );

//        // Create parameter objects.
//        std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet =
//                createEstimatableParameters( bodies, 1.1E7 );

//        testObservationPartials( angularPositionModel, bodies, fullEstimatableParameterSet, linkEnds, angular_position, 1.0E-4,
//        false, true, 1.0, parameterPerturbationMultipliers );

//    }
}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat





