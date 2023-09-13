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

#include <boost/test/unit_test.hpp>

#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/math/integrators/rungeKuttaCoefficients.h"
#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/createBodies.h"
#include "tudat/simulation/estimation_setup/createNumericalSimulator.h"
#include "tudat/math/integrators/createNumericalIntegrator.h"

namespace tudat
{

namespace unit_tests
{

using namespace simulation_setup;

BOOST_AUTO_TEST_SUITE( test_yarkovsky_acceleration )

BOOST_AUTO_TEST_CASE( testYarkovskyAccelerations )
{
    spice_interface::loadStandardSpiceKernels( );
    
    // Set simulation end epoch.
    const double simulationEndEpoch = 2.0 * tudat::physical_constants::JULIAN_DAY;

    // Set numerical integration fixed step size.
    const double fixedStepSize = 15.0;

    // Define body settings for simulation.
    BodyListSettings bodySettings = BodyListSettings(
                std::map< std::string, std::shared_ptr< BodySettings > >( ), "SSB", "J2000" );
    bodySettings.addSettings( "Earth" );
    bodySettings.at( "Earth" )->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >(
                Eigen::Vector6d::Zero( ), "SSB", "J2000" );
    bodySettings.at( "Earth" )->gravityFieldSettings = std::make_shared< GravityFieldSettings >( central_spice );

    // Create Earth object
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // Create spacecraft object.
    bodies.createEmptyBody( "Apophis" );

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // FIXME: This is not implemented yet.
    std::cerr << "The Yarkovski test is not implemented. Do not trust" << std::endl;
}

BOOST_AUTO_TEST_SUITE_END()

}

}
