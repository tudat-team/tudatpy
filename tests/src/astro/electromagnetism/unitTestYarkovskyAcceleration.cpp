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

#include <boost/test/unit_test.hpp>
#include <Eigen/Core>

#include "tudat/basics/testMacros.h"
#include "tudat/astro/electromagnetism/yarkovskyAcceleration.h"
#include "tudat/astro/basic_astro/physicalConstants.h"

#include "tudat/math/basic/linearAlgebra.h"
#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/math/integrators/rungeKuttaCoefficients.h"
#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/simulation/propagation_setup/dynamicsSimulator.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/createBodies.h"
#include "tudat/math/integrators/createNumericalIntegrator.h"


namespace tudat
{

namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_yarkovsky_acceleration )

const double AU = physical_constants::ASTRONOMICAL_UNIT;
const double JD = physical_constants::JULIAN_DAY;
const double yarkovskyParameter = -2.899e-14 * AU / ( JD * JD );

//! Test the implementation of the Yarkovsky Acceleration Model
BOOST_AUTO_TEST_CASE( testYarkovskyAccelerationVerySimple )
{
    const Eigen::Vector6d state = { AU, 0.0, 0.0, 0.0, AU, 0.0 };

    const Eigen::Vector3d expectedYarkovskyAcceleration = Eigen::Vector3d { 0.0, yarkovskyParameter, 0.0
    };
    const Eigen::Vector3d computedYarkovskyAcceleration = electromagnetism::computeYarkovskyAcceleration(
            yarkovskyParameter,
            state );

    TUDAT_CHECK_MATRIX_CLOSE( computedYarkovskyAcceleration, expectedYarkovskyAcceleration, 1.e-10 );
}

//! Test the implementation of the Yarkovsky Acceleration Model
BOOST_AUTO_TEST_CASE( testYarkovskyAccelerationSimple )
{
    const Eigen::Vector6d state = { 2.0 * AU, AU, 1.0e10, 3.0e5, 4.0e5, 12.0e5
    };

    const Eigen::Vector3d expectedYarkovskyAcceleration =
            yarkovskyParameter * 0.19982142476807888 * Eigen::Vector3d { 3.0, 4.0, 12.0 } / 13.0;
    const Eigen::Vector3d computedYarkovskyAcceleration = electromagnetism::computeYarkovskyAcceleration(
            yarkovskyParameter,
            state );

    TUDAT_CHECK_MATRIX_CLOSE( computedYarkovskyAcceleration, expectedYarkovskyAcceleration, 1.0e-10 );
}

Eigen::Vector6d bodyState = Eigen::Vector6d { AU, 0, 0, 12.0e5, 3.0e5, 4.0e5 };

Eigen::Vector6d getBodyState( ) { return bodyState; }

Eigen::Vector6d centralBodyState = Eigen::Vector6d::Zero( );

Eigen::Vector6d getCentralBodyState( ) { return centralBodyState; }

//! Test the class implementation of the Yarkovsky Acceleration model (updateMembers)
BOOST_AUTO_TEST_CASE( testYarkovskyAccelerationModelClassUpdateMembers )
{
    electromagnetism::YarkovskyAccelerationPointer yarkovskyAccelerationModel = std::make_shared< electromagnetism::YarkovskyAcceleration >(
            yarkovskyParameter,
            &getBodyState,
            &getCentralBodyState );

    yarkovskyAccelerationModel->updateMembers( 0.0 );
    const Eigen::Vector3d expectedYarkovskyAcceleration =
            yarkovskyParameter * Eigen::Vector3d { 12.0, 3.0, 4.0 } / 13.0;
    const Eigen::Vector3d computedYarkovskyAcceleration = yarkovskyAccelerationModel->getAcceleration( );

    TUDAT_CHECK_MATRIX_CLOSE( computedYarkovskyAcceleration, expectedYarkovskyAcceleration, 1.0e-10 )
}

//! Test the complete implementation of the yarkovsky Acceleration by comparing with Perez-Hernandez and Benet results
//! Reference: Pérez-Hernández, J. A., & Benet, L. (2022). Non-zero Yarkovsky acceleration for near-Earth asteroid (99942) Apophis. Communications Earth & Environment, 3(1), Article 1. https://doi.org/10.1038/s43247-021-00337-x
//! Reference: Perez-Hernandez et al.: https://static-content.springer.com/esm/art%3A10.1038%2Fs43247-021-00337-x/MediaObjects/43247_2021_337_MOESM1_ESM.pdf
//! Reference: Farnocchia, D., Chesley, S. R., Vokrouhlický, D., Milani, A., Spoto, F., & Bottke, W. F. (2013). Near Earth Asteroids with measurable Yarkovsky effect. Icarus, 224(1), 1–13. https://doi.org/10.1016/j.icarus.2013.02.004
BOOST_AUTO_TEST_CASE( testYarkovskyAccelerationHernandez )
{
    using namespace simulation_setup;

    spice_interface::loadStandardSpiceKernels( );

    // Simulation Setup Constants
    const std::string frameOrigin = "SSB";
    const std::string frameOrientation = "ECLIPJ2000";
    const double simulationStartEpoch = 2459200.5;
    const double simulationEndEpoch = 2459565.5;
    const double fixedStepSize = 1.0;

    // Apophis
    const std::string apophisName = "Apophis";
    const double apophisMass = 0.0;

    std::vector< std::string > bodiesToCreate { "Sun" };
    BodyListSettings bodySettings = getDefaultBodySettings( bodiesToCreate, frameOrigin, frameOrientation );

    // Create Sun object
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // Create Asteroid
    bodies.createEmptyBody( "Apophis" );
    bodies.at( apophisName )->setConstantBodyMass( apophisMass );

    // Define propagator settings variables.
    std::vector< std::string > bodiesToPropagate { apophisName };
    std::vector< std::string > centralBodies { "Sun" };

    // Accelerations
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfApophis = {
            { "Sun", { pointMassGravityAcceleration( ), yarkovskyAcceleration( yarkovskyParameter ) }},
    };

    // Make the acceleration models
    SelectedAccelerationMap accelerationSettings {{ apophisName, accelerationsOfApophis }};
    basic_astrodynamics::AccelerationMap accelerationModels = createAccelerationModelsMap( bodies,
                                                                                           accelerationSettings,
                                                                                           bodiesToPropagate,
                                                                                           centralBodies );

    // Initial State
    Eigen::VectorXd initialState = Eigen::Vector6d { -0.18034829, 0.94069106, 0.34573599, -0.0162659398, 4.39155,
                                                     -0.000395204
    };
    initialState *= physical_constants::ASTRONOMICAL_UNIT;
    initialState.tail( 3 ) /= physical_constants::JULIAN_DAY;

    // Termination Setting
    std::shared_ptr< propagators::PropagationTerminationSettings > terminationSettings = propagators::propagationTimeTerminationSettings(
            simulationEndEpoch,
            true );

    // Integrator Settings
    std::shared_ptr< numerical_integrators::IntegratorSettings< double>> integratorSettings = numerical_integrators::rungeKutta4Settings(
            fixedStepSize );

    // Propagator Settings
    std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double, double>> propagatorSettings = propagators::translationalStatePropagatorSettings(
            centralBodies,
            accelerationModels,
            bodiesToPropagate,
            initialState,
            simulationStartEpoch,
            integratorSettings,
            terminationSettings,
            propagators::cowell );

    // Dynamics Simulator
    propagators::SingleArcDynamicsSimulator< > dynamicsSimulator( bodies, propagatorSettings );

    // State history
    std::map< double, Eigen::VectorXd > stateHist = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

    const double sunGravitationalParameter = bodies.at( "Sun" )->getGravitationalParameter( );
    Eigen::Vector6d keplerInitialState = orbital_element_conversions::convertCartesianToKeplerianElements( Eigen::Vector6d(
            initialState ), sunGravitationalParameter );
    Eigen::Vector6d keplerFinalState = orbital_element_conversions::convertCartesianToKeplerianElements( Eigen::Vector6d(
            stateHist.at( simulationEndEpoch )), sunGravitationalParameter );

    // Check drift in semi-major axis
    const double expectedFinalSemiMajorAxis = keplerInitialState[0] - 199.0;
    BOOST_CHECK_CLOSE( expectedFinalSemiMajorAxis, keplerFinalState[0], 1e-10 );
}

BOOST_AUTO_TEST_SUITE_END( )

}

}