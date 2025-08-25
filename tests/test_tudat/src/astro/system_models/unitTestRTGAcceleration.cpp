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

#include <memory>

#include "tudat/astro/aerodynamics/testApolloCapsuleCoefficients.h"
#include "tudat/astro/basic_astro/sphericalStateConversions.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/ephemerides/directionBasedRotationalEphemeris.h"
#include "tudat/astro/reference_frames/referenceFrameTransformations.h"
#include "tudat/basics/testMacros.h"
#include "tudat/simulation/propagation_setup/dynamicsSimulator.h"
#include "tudat/interface/spice/spiceEphemeris.h"
#include "tudat/interface/spice/spiceRotationalEphemeris.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/io/multiDimensionalArrayReader.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/estimation_setup/createNumericalSimulator.h"
#include "tudat/simulation/propagation_setup/createMassRateModels.h"
#include "tudat/simulation/estimation_setup/variationalEquationsSolver.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/createSystemModel.h"
#include <limits>
#include <string>

#include <Eigen/Core>

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_rtg_acceleration )


// Test case 1: simple setup of rtgAccelerationModel, Vehicle body with rotation ephemeris and mass function.
// Updating vehicle environment sequentially to isolate effects from environment, then calling updateMembers on Acceleration Model and comparing with expected values.
// First (1a): evaluate acceleration vector at acceleration model referenceTime, then testTime to check correct modelling of decay process.
// Second (1b): update vehicle mass at testTime to follow vehicle mass function, re-evaluate acceleration model current acceleration to check correct interfacing with vehicle mass function.
// Third (1c): update vehicle rotational state at testTime to match rotational ephemeris, re-evaluate acceleration model current acceleration to check correct interfacing with vehicle rotation state.

// Test case 2: Acceleration model and environment setup described above is embedded into a numerical integration and propagated for 24h.
//              First environment dependencies vehicle mass, rotation state) are reset to inital state.
//              Then acceleration model with all dependencies (decay, mass function and rotation ephemeris effects included) are propagated.
//              RTG acceleration is recovered from propagation through dependent variable and compared to expected value.


BOOST_AUTO_TEST_CASE( testRTGAcceleration )
{
    using namespace tudat::simulation_setup;
    using namespace tudat;

    ////////////////////////////////////////////////////////////////
    ///       SETUP                                              ///
    ////////////////////////////////////////////////////////////////

    // Load Spice kernels
    spice_interface::loadStandardSpiceKernels( );

    // Get settings for Earth.
    BodyListSettings bodySettings;
    bodySettings.addSettings( getDefaultSingleBodySettings( "Earth", 0.0, 86400.0 ), "Earth" );
    bodySettings.addSettings( "Vehicle" );

    // Create body objects.
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // Define Relevant Epochs
    double referenceEpoch = 0.0;
    double testTime = 500;

    // Define function describing rotational ephemeris of vehicle
    std::function<Eigen::Matrix3d(double)> timeDependentRotationFunction =
    [](double epoch) {
        double angleRad = 1/7. * epoch * mathematical_constants::PI / 180.0;
        return Eigen::AngleAxisd(angleRad, Eigen::Vector3d::UnitZ()).toRotationMatrix();
    };

    // Assign and update
    bodies.at( "Vehicle" )
            ->setRotationalEphemeris( createRotationModel( std::make_shared< CustomRotationModelSettings >(
                                                                   "ECLIPJ2000",
                                                                   "VehicleFixed",
                                                                   timeDependentRotationFunction,
                                                                   1.0 ),
                                                           "Vehicle",
                                                           bodies ) );
    bodies.at( "Vehicle" )->setCurrentRotationalStateToLocalFrameFromEphemeris( referenceEpoch );

    // Define function describing mass function of vehicle
    double initialVehicleMass = 5000;

    std::function<double(double)> vehicleMassFunction =
        [=](double epoch) {
            //double delta_epoch = epoch - referenceEpoch;
            double delta_test = epoch - testTime;
            if (delta_test > -100. && delta_test <= 200.) {
                return initialVehicleMass - (delta_test+100.);
            } else if (delta_test <= -100.) {
                return initialVehicleMass;
            } else {
                return initialVehicleMass - 100.;
            }
    };

    // Assign and Update
    bodies.at( "Vehicle" )->setBodyMassFunction( vehicleMassFunction );
    bodies.at( "Vehicle" )->updateMass( referenceEpoch );

    // Define settings for accelerations
    SelectedAccelerationMap accelerationSettingsMap;

    Eigen::Vector3d rtgForceVector;
    rtgForceVector << 0.5E-5, 0.5E-5, 0.5E-5;
    double decayScaleFactor = 1.5e-05;       // corresponding to a half-life of roughly half a day

    // Define origin of integration
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    bodiesToPropagate.push_back( "Vehicle" );
    centralBodies.push_back( "Earth" );


    // Alternative: Create Settings object directly via class constructor
    //accelerationSettingsMap[ "Vehicle" ][ "Vehicle" ].push_back(
    //                    std::make_shared< RTGAccelerationSettings >(rtgForceVector, decayScaleFactor, referenceEpoch));

    //  Create Settings object directly via factory function
    accelerationSettingsMap[ "Vehicle" ][ "Vehicle" ].push_back( rtgAcceleration(rtgForceVector, decayScaleFactor, referenceEpoch) );


    // Create accelerations
    basic_astrodynamics::AccelerationMap accelerationsMap = createAccelerationModelsMap( bodies, accelerationSettingsMap, bodiesToPropagate, centralBodies );
    //std::shared_ptr< basic_astrodynamics::AccelerationModel3d > rtgAccelerationModel = accelerationsMap[ "Vehicle"] ["Vehicle"][ 0 ];
    std::shared_ptr< basic_astrodynamics::AccelerationModel3d > accelerationModel = accelerationsMap[ "Vehicle"] ["Vehicle"][ 0 ];


    // Dynamic cast acceleration settings to required type and check consistency.
    std::shared_ptr< system_models::RTGAccelerationModel > rtgAccelerationModel =
            std::dynamic_pointer_cast< system_models::RTGAccelerationModel >( accelerationModel );

    // Test Successful Construction
    BOOST_CHECK_EQUAL( rtgAccelerationModel != nullptr, true );


    ////////////////////////////////////////////////////////////////
    ///       Test 1                                             ///
    ////////////////////////////////////////////////////////////////


    ///////   Test 1a)  ////////////////////////////////////////////

    rtgAccelerationModel->updateMembers( referenceEpoch );
    Eigen::Vector3d expectedAcceleration = rtgForceVector / initialVehicleMass;

    for (int i = 0; i < 3; ++i)
    {
        BOOST_CHECK_CLOSE(expectedAcceleration[i], rtgAccelerationModel->getAcceleration( )[i], 1e-10);  //
    }

    // re-evaluate current acceleration at testTime, without updating mass or rotational state of vehicle --> isolated effect of modelled decay process
    rtgAccelerationModel->updateMembers( testTime );

    // check acceleration equivalence on test epoch values
    double expectedDecay = std::exp(-decayScaleFactor*(testTime-referenceEpoch));
    expectedAcceleration = (rtgForceVector  / initialVehicleMass) * expectedDecay;

    for (int i = 0; i < 3; ++i)
    {
        BOOST_CHECK_CLOSE(expectedAcceleration[i], rtgAccelerationModel->getAcceleration( )[i], 1e-10);  //
    }

    ///////   Test 1b)  ////////////////////////////////////////////

    bodies.at( "Vehicle" )->updateMass( testTime );
    // re-evaluate current acceleration at testTime, after updating mass of vehicle --> effect of modelled decay process + mass decrease
    rtgAccelerationModel->updateMembers( testTime );

    expectedAcceleration *= initialVehicleMass / bodies.at( "Vehicle" )->getBodyMass(  );

    for (int i = 0; i < 3; ++i)
    {
        BOOST_CHECK_CLOSE(expectedAcceleration[i], rtgAccelerationModel->getAcceleration( )[i], 1e-10);  //
    }


    ///////   Test 1c)  ////////////////////////////////////////////

    // re-evaluate current acceleration at testTime, after updating rotational state of vehicle --> effect of modelled decay process + mass decrease + rotation
    bodies.at( "Vehicle" )->setCurrentRotationalStateToLocalFrameFromEphemeris( testTime );
    rtgAccelerationModel->updateMembers( testTime );

    Eigen::Matrix3d rotationMatrixFromFunction = timeDependentRotationFunction( testTime );
    expectedAcceleration = rotationMatrixFromFunction * expectedAcceleration;

    for (int i = 0; i < 3; ++i)
    {
        BOOST_CHECK_CLOSE(expectedAcceleration[i], rtgAccelerationModel->getAcceleration( )[i], 1e-10);  //
    }

    double newTestTime = -500;
    bodies.at( "Vehicle" )->setCurrentRotationalStateToLocalFrameFromEphemeris( newTestTime );
    bodies.at( "Vehicle" )->updateMass( newTestTime );
    rtgAccelerationModel->updateMembers( newTestTime );

    rotationMatrixFromFunction = timeDependentRotationFunction( newTestTime );
    std::exp(-decayScaleFactor*(newTestTime-referenceEpoch));
    expectedAcceleration = rtgForceVector * std::exp(-decayScaleFactor*(newTestTime-referenceEpoch));
    expectedAcceleration = rotationMatrixFromFunction * expectedAcceleration / initialVehicleMass;

    for (int i = 0; i < 3; ++i)
    {
        BOOST_CHECK_CLOSE(expectedAcceleration[i], rtgAccelerationModel->getAcceleration( )[i], 1e-10);  //
    }


    ////////////////////////////////////////////////////////////////
    ///       Test 2                                             ///
    ////////////////////////////////////////////////////////////////

    // Reset environment dependencies
    bodies.at( "Vehicle" )->setCurrentRotationalStateToLocalFrameFromEphemeris( referenceEpoch );
    bodies.at( "Vehicle" )->updateMass( referenceEpoch );


    Eigen::Vector6d systemInitialState = Eigen::Vector6d::Zero( );
    systemInitialState( 0 ) = 8.0E6;
    systemInitialState( 4 ) = 7.5E3;

    std::vector< std::shared_ptr< propagators::SingleDependentVariableSaveSettings > > dependentVariables;
    dependentVariables.push_back(std::make_shared< propagators::SingleAccelerationDependentVariableSaveSettings >(
        basic_astrodynamics::rtg_acceleration, "Vehicle", "Vehicle", 0 ) );


    // Create Propagator, Integrator objects
    std::shared_ptr< propagators::PropagationTimeTerminationSettings > terminationSettings =
        std::make_shared< propagators::PropagationTimeTerminationSettings >( 1000. );
    std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
            std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >(
                    centralBodies, accelerationsMap, bodiesToPropagate, systemInitialState, terminationSettings, propagators::cowell, dependentVariables );
    std::shared_ptr< numerical_integrators::IntegratorSettings<> > integratorSettings = std::make_shared< numerical_integrators::IntegratorSettings<> >( numerical_integrators::rungeKutta4, 0.0, 0.1 );


    // Create simulation object and propagate dynamics.
    propagators::SingleArcDynamicsSimulator<> dynamicsSimulator( bodies, integratorSettings, translationalPropagatorSettings, true, false, false );


    // Retrieve numerical solutions for state and dependent variables
    std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > numericalSolution =
            dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

    std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > dependentVariableSolution =
            dynamicsSimulator.getDependentVariableHistory( );


    Eigen::Vector3d test_acceleration;
    double matchedEpoch;
    for( std::map< double, Eigen::VectorXd >::iterator variableIterator = dependentVariableSolution.begin( );
     variableIterator != dependentVariableSolution.end( );
     variableIterator++ )
    {
        // Retrieve data at test epoch from dependent variables/propagated dynamics

        const double tolerance = 1.0e-5;
        if (std::abs(variableIterator->first - testTime) < tolerance)
        {
            // Retrieve data at test epoch from dependent variables/propagated dynamics

            matchedEpoch = variableIterator->first;
            test_acceleration = variableIterator->second;
        }
    }

    // Match of testEpoch and keys in dependentVariablesHistory is only good to order 1e-7
    // So in order to assert similarity of acceleration at the usual level of 1e-10, we recompute the expected acceleration at the epoch from the dependentVariablesHistory

    double testVehicleMass = vehicleMassFunction( matchedEpoch );
    Eigen::Matrix3d rotationMatrix = timeDependentRotationFunction( matchedEpoch );
    Eigen::Vector3d expectedAccelerationAtMatchedEpoch = (rotationMatrix*rtgForceVector * std::exp(-decayScaleFactor * (matchedEpoch-referenceEpoch)));
    expectedAccelerationAtMatchedEpoch = expectedAccelerationAtMatchedEpoch / testVehicleMass;

    for (int i = 0; i < 3; ++i)
    {
        BOOST_CHECK_CLOSE(expectedAccelerationAtMatchedEpoch[i], test_acceleration[i], 1e-10);  //
    }

}

BOOST_AUTO_TEST_SUITE_END(  )
    }  // namespace unit_tests

// namespace tudat
}


