/*    Copyright (c) 2010-2023, Delft University of Technology
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

#include <limits>
#include <string>
#include <vector>

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"

#include "tudat/io/basicInputOutput.h"
#include "tudat/interface/spice/spiceInterface.h"

#include "tudat/simulation/estimation_setup/createObservationModelFactory.h"
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

using namespace tudat;
using namespace tudat::unit_tests;
using namespace tudat::gravitation;
using namespace tudat::ephemerides;
using namespace tudat::observation_models;
using namespace tudat::simulation_setup;
using namespace tudat::spice_interface;
using namespace tudat::observation_partials;
using namespace tudat::estimatable_parameters;

BOOST_AUTO_TEST_SUITE( test_one_way_doppler_measured_frequency_partials )

//! Test partial derivatives of one-way Doppler measured frequency observable, using general test suite.
BOOST_AUTO_TEST_CASE( testOneWayDopplerMeasuredFrequencyPartials )
{
    // Define and create ground stations.
    std::vector< std::pair< std::string, std::string > > groundStations;
    groundStations.resize( 2 );
    groundStations[ 0 ] = std::make_pair( "Earth", "Graz" );
    groundStations[ 1 ] = std::make_pair( "Mars", "MSL" );

    Eigen::VectorXd parameterPerturbationMultipliers = Eigen::Vector4d::Constant( 1.0 );

    // Test partials with constant ephemerides (allows test of position partials)
    {
        // Create environment
        SystemOfBodies bodies = setupEnvironment( groundStations, 1.0E7, 1.2E7, 1.1E7, true );

        // Set transmitter frequency on Mars ground station
        bodies.getBody( "Mars" )->getGroundStation( "MSL" )->setTransmittingFrequencyCalculator(
                std::make_shared< ground_stations::ConstantFrequencyInterpolator >( 8.4e9 ) );

        // Set link ends for observation model
        LinkDefinition linkEnds;
        linkEnds[ transmitter ] = groundStations[ 1 ];
        linkEnds[ receiver ] = groundStations[ 0 ];

        std::cout << "Case A: Constant ephemerides " << std::endl;

        // Create one-way Doppler measured frequency model via convenience function
        std::vector< std::string > perturbingBodies;
        perturbingBodies.push_back( "Earth" );
        std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionSettings;
        lightTimeCorrectionSettings.push_back(
                std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( perturbingBodies ) );

        std::shared_ptr< ObservationModel< 1 > > measuredFreqModel =
                observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                        oneWayDopplerMeasuredFrequencySettings(
                                linkEnds,
                                lightTimeCorrectionSettings ),
                        bodies );

        // Create parameter objects.
        std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet =
                createEstimatableParameters( bodies, 1.1E7 );

        printEstimatableParameterEntries( fullEstimatableParameterSet );

        testObservationPartials< 1 >( measuredFreqModel,
                                      bodies,
                                      fullEstimatableParameterSet,
                                      linkEnds,
                                      one_way_doppler_measured_frequency,
                                      5.0E-5,
                                      true,
                                      true,
                                      10.0,
                                      parameterPerturbationMultipliers,
                                      nullptr,
                                      1.1E7,
                                      20.0 );
    }

    // Test partials with real ephemerides (without test of position partials)
    {
        // Create environment
        SystemOfBodies bodies = setupEnvironment( groundStations, 1.0E7, 1.2E7, 1.1E7, false );

        // Set transmitter frequency on Mars ground station
        bodies.getBody( "Mars" )->getGroundStation( "MSL" )->setTransmittingFrequencyCalculator(
                std::make_shared< ground_stations::ConstantFrequencyInterpolator >( 8.4e9 ) );

        // Set link ends for observation model
        LinkDefinition linkEnds;
        linkEnds[ transmitter ] = groundStations[ 1 ];
        linkEnds[ receiver ] = groundStations[ 0 ];

        std::cout << "Case B: Real ephemerides " << std::endl;

        // Create one-way Doppler measured frequency model via convenience function
        std::vector< std::string > perturbingBodies;
        perturbingBodies.push_back( "Earth" );
        std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionSettings;
        lightTimeCorrectionSettings.push_back(
                std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( perturbingBodies ) );

        std::shared_ptr< ObservationModel< 1 > > measuredFreqModel =
                observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                        oneWayDopplerMeasuredFrequencySettings(
                                linkEnds,
                                lightTimeCorrectionSettings ),
                        bodies );

        // Create parameter objects.
        std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet =
                createEstimatableParameters( bodies, 1.1E7 );

        printEstimatableParameterEntries( fullEstimatableParameterSet );

        testObservationPartials< 1 >( measuredFreqModel,
                                      bodies,
                                      fullEstimatableParameterSet,
                                      linkEnds,
                                      one_way_doppler_measured_frequency,
                                      1.0E-4,
                                      false,
                                      true,
                                      1.0,
                                      parameterPerturbationMultipliers,
                                      nullptr,
                                      1.1E7,
                                      20.0 );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
