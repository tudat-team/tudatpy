/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Reference: Pérez-Hernández, J. A., & Benet, L. (2022). Non-zero Yarkovsky acceleration for near-Earth asteroid (99942) Apophis. Communications Earth & Environment, 3(1), Article 1. https://doi.org/10.1038/s43247-021-00337-x
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>
#include <Eigen/Core>

#include "tudat/basics/testMacros.h"
#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/astro/electromagnetism/yarkovskyAcceleration.h"
#include "tudat/math/basic/linearAlgebra.h"
#include "tudat/math/basic/leastSquaresEstimation.h"
#include "tudat/math/integrators/rungeKuttaCoefficients.h"
#include "tudat/math/integrators/createNumericalIntegrator.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/simulation/propagation_setup/dynamicsSimulator.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/createBodies.h"


namespace tudat
{

namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_yarkovsky_acceleration )

const double AU = physical_constants::ASTRONOMICAL_UNIT;
const double JD = physical_constants::JULIAN_DAY;
const double yarkovskyParameter = -2.899e-14 * AU / ( JD * JD ); // A2 for Apophis (Pérez-Hernández & Benet, 2022)

//! Compute mean motion of the orbit
double computeMeanMotion( const double semiMajorAxis, const double mu )
{
    return std::sqrt( mu / semiMajorAxis / semiMajorAxis / semiMajorAxis );
}

//! Computer orbital Period
double computeOrbitalPeriod( const double semiMajorAxis, const double mu )
{
    return 2 * M_PI / computeMeanMotion( semiMajorAxis, mu );
}

//! Yarkovsky acceleration is A2 * (AU / rSun)^2 in the direction of the velocity vector (unit length)
Eigen::Vector3d computeExpectedYarkovskyAcceleration( const double A2, const Eigen::Vector6d& state )
{
    const double r0 = AU;
    const double rSun = state.head( 3 ).norm( );
    const Eigen::Vector3d tHat = state.tail( 3 ).normalized( );
    const Eigen::Vector3d acceleration = A2 * ( r0 * r0 ) / ( rSun * rSun ) * tHat;
    return acceleration;
}

//! Expected drift in semi-major axis (m/s) according to a.o. Perez-Hernandez & Benet (2022)
double computeExpectedSemiMajorAxisDrift( const Eigen::Vector6d& keplerElements, const double mu )
{
    const double a = keplerElements[0];
    const double e = keplerElements[1];

    const double p = a * ( 1 - e * e ); // Semi-latus rectum
    const double n = computeMeanMotion( a, mu ); // Mean motion

    const double semiMajorAxisDrift = 2 * yarkovskyParameter * ( 1 - e * e ) * AU * AU / ( n * p * p );
    return semiMajorAxisDrift;
}

//! Test the implementation of the Yarkovsky Acceleration Model
BOOST_AUTO_TEST_CASE( testYarkovskyAccelerationVerySimple )
{
    const Eigen::Vector6d state = { 0.5 * AU, 0.0, 0.0, 0.0, 5e4, 0.0 }; // Arbitrary numbers

    const Eigen::Vector3d computedYarkovskyAcceleration = electromagnetism::computeYarkovskyAcceleration(
            yarkovskyParameter,
            state );
    const Eigen::Vector3d expectedYarkovskyAcceleration = computeExpectedYarkovskyAcceleration( yarkovskyParameter,
                                                                                                state );

    TUDAT_CHECK_MATRIX_CLOSE( computedYarkovskyAcceleration, expectedYarkovskyAcceleration, 1.e-10 );
}

//! Test the implementation of the Yarkovsky Acceleration Model
BOOST_AUTO_TEST_CASE( testYarkovskyAccelerationSimple )
{
    const Eigen::Vector6d state = { 2.0 * AU, AU, 3.0 * AU, 3.0e5, 4.0e5, 12.0e5 };
    const Eigen::Vector3d expectedYarkovskyAcceleration = computeExpectedYarkovskyAcceleration( yarkovskyParameter,
                                                                                                state );
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
    const Eigen::Vector3d expectedYarkovskyAcceleration = computeExpectedYarkovskyAcceleration( yarkovskyParameter,
                                                                                                bodyState );
    const Eigen::Vector3d computedYarkovskyAcceleration = yarkovskyAccelerationModel->getAcceleration( );

    TUDAT_CHECK_MATRIX_CLOSE( computedYarkovskyAcceleration, expectedYarkovskyAcceleration, 1.0e-10 )
}

//! Test the complete implementation of the yarkovsky Acceleration using the drift in semi-major axis for a circular orbit
//! Reference: Pérez-Hernández, J. A., & Benet, L. (2022). Non-zero Yarkovsky acceleration for near-Earth asteroid (99942) Apophis. Communications Earth & Environment, 3(1), Article 1. https://doi.org/10.1038/s43247-021-00337-x
//! For Circular Orbit
BOOST_AUTO_TEST_CASE( testYarkovskyAccelerationCircular )
{
    using namespace simulation_setup;
    spice_interface::loadStandardSpiceKernels( );
    const double MU_SUN = spice_interface::getBodyGravitationalParameter( "Sun" );

    // Initial Conditions
    const double initialSemiMajorAxis = 0.75 * AU;
    const double initialEccentricity = 0.05;
    const double initialInclination = 0.1;
    const double initialRAAN = 0.0;
    const double initialArgOfPeri = 0.0;
    const double initialTrueAnomaly = 0.0;

    // Simulation Setup Constants
    const std::string frameOrigin = "SSB";
    const std::string frameOrientation = "ECLIPJ2000";
    const double simulationStartEpoch = 0.0;
    const double simulationDuration = computeOrbitalPeriod( initialSemiMajorAxis, MU_SUN );
    const double simulationEndEpoch = simulationStartEpoch + simulationDuration;
    const double fixedStepSize = 1000;

    // Bodies
    std::vector< std::string > bodiesToCreate { "Sun" };
    BodyListSettings bodySettings = getDefaultBodySettings( bodiesToCreate, frameOrigin, frameOrientation );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // Asteroid
    bodies.createEmptyBody( "Asteroid" );
    bodies.at( "Asteroid" )->setConstantBodyMass( 0.0 );

    // Define propagator settings variables.
    std::vector< std::string > bodiesToPropagate { "Asteroid" };
    std::vector< std::string > centralBodies { "Sun" };

    // Accelerations
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfAsteroid = {
            { "Sun", { pointMassGravityAcceleration( ), yarkovskyAcceleration( yarkovskyParameter ) }},
    };

    // Make the acceleration models
    SelectedAccelerationMap accelerationSettings {{ "Asteroid", accelerationsOfAsteroid }};
    basic_astrodynamics::AccelerationMap accelerationModels = createAccelerationModelsMap( bodies,
                                                                                           accelerationSettings,
                                                                                           bodiesToPropagate,
                                                                                           centralBodies );

    // Initial State
    Eigen::VectorXd initialState = orbital_element_conversions::convertKeplerianToCartesianElements(
            initialSemiMajorAxis,
            initialEccentricity,
            initialInclination,
            initialRAAN,
            initialArgOfPeri,
            initialTrueAnomaly,
            MU_SUN );

    // Dependent Variables - Storing the semi-major axis (for drift test) and Yarkovsky acceleration
    std::vector< std::shared_ptr< propagators::SingleDependentVariableSaveSettings>> dependentVariablesToSave = {
            propagators::singleAccelerationDependentVariable( basic_astrodynamics::yarkovsky_acceleration,
                                                              "Asteroid",
                                                              "Sun" ),
            propagators::keplerianStateDependentVariable( "Asteroid", "Sun" ),
    };

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
            propagators::cowell,
            dependentVariablesToSave );

    // Dynamics Simulator
    propagators::SingleArcDynamicsSimulator< > dynamicsSimulator( bodies, propagatorSettings );

    // State histories
    std::map< double, Eigen::VectorXd > stateHist = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > depVarHist = dynamicsSimulator.getDependentVariableHistory( );

    // Map of semi-major axes and Yarkokvsky Accelerations;
    std::map< double, double > semiMajorAxes;
    std::map< double, Eigen::Vector3d > yarkovskyAccelerations;

    for ( const auto& pair: depVarHist ) {
        const double& time = pair.first;
        const Eigen::VectorXd& vector = pair.second;
        yarkovskyAccelerations[time] = vector.head( 3 ); // First three elements of the depVar vector
        semiMajorAxes[time] = vector[3]; // Fourth element of the depVar vector
    }


    // Calculate Expected Drift
    const Eigen::Vector6d initialKeplerElements = depVarHist[simulationStartEpoch].segment( 3, 6 );
    const Eigen::Vector6d endKeplerElements = depVarHist[simulationEndEpoch].segment( 3, 6 );

    const double expectedStartSemiMajorAxisDrift = computeExpectedSemiMajorAxisDrift( initialKeplerElements, MU_SUN );
    const double expectedEndSemiMajorAxisDrift = computeExpectedSemiMajorAxisDrift( endKeplerElements, MU_SUN );
    const double expectedSemiMajorAxisDrift = ( expectedStartSemiMajorAxisDrift + expectedEndSemiMajorAxisDrift ) / 2.0;

    // Calculate Simulated Drift
    const double semiMajorAxisDrift = ( endKeplerElements[0] - initialSemiMajorAxis ) / simulationDuration;

    // Check that the change of the drift is small in the simulation time frame
    BOOST_CHECK_SMALL( expectedEndSemiMajorAxisDrift - expectedStartSemiMajorAxisDrift, 1.0e-4 );

    // Check drift in semi-major axis
    BOOST_CHECK_CLOSE( expectedSemiMajorAxisDrift, semiMajorAxisDrift, 1.0e-1 );

    // Check parallel with velocity vector
    const Eigen::Vector3d endYarkovskyAcc = yarkovskyAccelerations[simulationEndEpoch];
    const Eigen::Vector3d endVelocity = stateHist[simulationEndEpoch].tail( 3 );
    BOOST_CHECK_SMALL( endYarkovskyAcc.cross( endVelocity ).norm( ), 1.0e-20 );
}

BOOST_AUTO_TEST_SUITE_END( )

}

}