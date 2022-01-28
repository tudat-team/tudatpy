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

BOOST_AUTO_TEST_SUITE( test_euler_angle_partials)

//! Test partial derivatives of Euler angle observable, using general test suite of observation partials.
BOOST_AUTO_TEST_CASE( testEulerAnglePartials )
{

    std::cout<<" **************************************************************************************** "<<std::endl;

    // Test partials with constant rotational ephemerides (with test of rotation state partial)
    {

        // Create environment
        SystemOfBodies bodies = setupEnvironment( std::vector< std::pair< std::string, std::string > >( ),
                                                 1.0E7, 1.2E7, 1.1E7, false, 1.0, true );

        // Set link ends for observation model
        LinkEnds linkEnds;
        linkEnds[ observed_body ] = std::make_pair( "Mars", "" );

        // Generate one-way range model
        std::shared_ptr< ObservationModel< 3 > > eulerAngleModel =
                observation_models::ObservationModelCreator< 3, double, double >::createObservationModel(
                    std::make_shared< observation_models::ObservationModelSettings >(
                        euler_angle_313_observable, linkEnds ), bodies  );

        // Create parameter objects.
        std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet =
                createEstimatableParameters( bodies, 1.1E7, false, true );

        testObservationPartials< 3 >(
                    eulerAngleModel, bodies, fullEstimatableParameterSet, linkEnds, euler_angle_313_observable, 1.0E-6, false, false );
    }
}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat




