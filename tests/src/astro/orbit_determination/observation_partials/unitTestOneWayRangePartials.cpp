/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
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
#include "tudat/astro/observation_models/oneWayRangeObservationModel.h"
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

BOOST_AUTO_TEST_SUITE( test_one_way_observation_partials)

//! Test partial derivatives of one-way range observable, using general test suite of observation partials.
BOOST_AUTO_TEST_CASE( testOneWayRangePartials )
{

    // Define and create ground stations.
    std::vector< std::pair< std::string, std::string > > groundStations;
    groundStations.resize( 2 );
    groundStations[ 0 ] = std::make_pair( "Earth", "Graz" );
    groundStations[ 1 ] = std::make_pair( "Mars", "MSL" );

    // Test partials with constant ephemerides (allows test of position partials)
    {
        // Create environment
        SystemOfBodies bodies = setupEnvironment( groundStations, 1.0E7, 1.2E7, 1.1E7, true );

        // Set link ends for observation model
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = groundStations[ 1 ];
        linkEnds[ receiver ] = groundStations[ 0 ];

        // Generate one-way range model
        std::vector< std::string > perturbingBodies;
        perturbingBodies.push_back( "Earth" );
        std::shared_ptr< ObservationModel< 1 > > oneWayRangeModel =
                observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                    std::make_shared< observation_models::ObservationModelSettings >(
                        observation_models::one_way_range, linkEnds,
                        std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
         perturbingBodies ) ), bodies  );

        // Create parameter objects.
        std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet =
                createEstimatableParameters( bodies, 1.1E7 );

        testObservationPartials< 1 >(
                    oneWayRangeModel, bodies, fullEstimatableParameterSet, linkEnds, one_way_range, 1.0E-6, true, true );
    }

    std::cout<<" **************************************************************************************** "<<std::endl;

    // Test partials with real ephemerides (without test of position partials)
    {
        // Create environment
        SystemOfBodies bodies = setupEnvironment( groundStations, 1.0E7, 1.2E7, 1.1E7, false );

        // Set link ends for observation model
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = groundStations[ 1 ];
        linkEnds[ receiver ] = groundStations[ 0 ];

        // Generate one-way range model
        std::vector< std::string > perturbingBodies;
        perturbingBodies.push_back( "Earth" );
        std::shared_ptr< ObservationModel< 1 > > oneWayRangeModel =
                observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                    std::make_shared< observation_models::ObservationModelSettings >(
                        observation_models::one_way_range, linkEnds,
                        std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
         perturbingBodies ) ), bodies  );

        // Create parameter objects.
        std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet =
                createEstimatableParameters( bodies, 1.1E7 );

        testObservationPartials< 1 >(
                    oneWayRangeModel, bodies, fullEstimatableParameterSet, linkEnds, one_way_range, 1.0E-6, false, true );
    }

    std::cout<<" **************************************************************************************** "<<std::endl;

    // Test partials with constant rotational ephemerides (with test of rotation state partial)
    {

        // Create environment
        SystemOfBodies bodies = setupEnvironment( groundStations, 1.0E7, 1.2E7, 1.1E7, false, 1.0, true );

        // Set link ends for observation model
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = groundStations[ 1 ];
        linkEnds[ receiver ] = groundStations[ 0 ];

        // Generate one-way range model
        std::vector< std::string > perturbingBodies;
        perturbingBodies.push_back( "Earth" );
        std::shared_ptr< ObservationModel< 1 > > oneWayRangeModel =
                observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                    std::make_shared< observation_models::ObservationModelSettings >(
                        observation_models::one_way_range, linkEnds,
                        std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
         perturbingBodies ) ), bodies  );

        // Create parameter objects.
        std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet =
                createEstimatableParameters( bodies, 1.1E7, false, true );

        testObservationPartials< 1 >(
                    oneWayRangeModel, bodies, fullEstimatableParameterSet, linkEnds, one_way_range, 1.0E-6, false, true );
    }
}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat




