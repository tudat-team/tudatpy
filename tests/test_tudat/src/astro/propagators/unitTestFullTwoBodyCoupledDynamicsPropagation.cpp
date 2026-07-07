/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/torqueModel.h"
#include "tudat/astro/gravitation/sphericalHarmonicsGravityField.h"
#include "tudat/basics/testMacros.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/math/basic/linearAlgebra.h"
#include "tudat/math/integrators/rungeKutta4Integrator.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/createGravityField.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/propagation_setup/createAccelerationModels.h"
#include "tudat/simulation/propagation_setup/createTorqueModel.h"
#include "tudat/simulation/propagation_setup/dynamicsSimulatorBase.h"
#include "tudat/simulation/propagation_setup/propagationOutputSettings.h"
#include "tudat/simulation/propagation_setup/propagationSettings.h"
#include "tudat/simulation/propagation_setup/singleArcDynamicsSimulator.h"

namespace tudat
{
namespace unit_tests
{

using namespace basic_astrodynamics;
using namespace numerical_integrators;
using namespace propagators;
using namespace simulation_setup;

SystemOfBodies createIoJupiterSystemWithDegreeTwoGravityFields( const double initialEpoch, const double finalEpoch )
{
    const std::vector< std::string > bodiesToCreate = { "Jupiter", "Io" };
    const std::map< std::string, double > gravitationalParameters = { { "Jupiter", 1.266865349218008E17 }, { "Io", 5.959916E12 } };
    const std::map< std::string, double > referenceRadii = { { "Jupiter", 69911.0E3 }, { "Io", 1821.6E3 } };
    BodyListSettings bodySettings =
            getDefaultBodySettings( bodiesToCreate, initialEpoch - 3600.0, finalEpoch + 3600.0, "Jupiter", "ECLIPJ2000" );

    for( const std::string& bodyName : bodiesToCreate )
    {
        bodySettings.at( bodyName )->gravityFieldSettings = nullptr;
        bodySettings.at( bodyName )->rigidBodyPropertiesSettings = nullptr;
    }

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    for( const std::string& bodyName : bodiesToCreate )
    {
        const double gravitationalParameter = gravitationalParameters.at( bodyName );
        const double referenceRadius = referenceRadii.at( bodyName );

        Eigen::MatrixXd cosineCoefficients = Eigen::MatrixXd::Zero( 3, 3 );
        Eigen::MatrixXd sineCoefficients = Eigen::MatrixXd::Zero( 3, 3 );
        cosineCoefficients( 0, 0 ) = 1.0;

        if( bodyName == "Jupiter" )
        {
            cosineCoefficients( 2, 0 ) = 1.20E-2;
            cosineCoefficients( 2, 1 ) = -0.90E-2;
            cosineCoefficients( 2, 2 ) = 1.10E-2;
            sineCoefficients( 2, 1 ) = 0.80E-2;
            sineCoefficients( 2, 2 ) = -1.00E-2;
        }
        else
        {
            cosineCoefficients( 2, 0 ) = -0.65E-2;
            cosineCoefficients( 2, 1 ) = 1.35E-2;
            cosineCoefficients( 2, 2 ) = -1.55E-2;
            sineCoefficients( 2, 1 ) = -1.25E-2;
            sineCoefficients( 2, 2 ) = 0.45E-2;
        }

        addGravityFieldModel( bodies,
                              bodyName,
                              std::make_shared< SphericalHarmonicsGravityFieldSettings >(
                                      gravitationalParameter,
                                      referenceRadius,
                                      cosineCoefficients,
                                      sineCoefficients,
                                      bodies.at( bodyName )->getRotationalEphemeris( )->getTargetFrameOrientation( ),
                                      0.35 ),
                              std::vector< std::shared_ptr< GravityFieldVariationSettings > >( ) );
    }

    return bodies;
}

void setCoupledIoStateInEnvironment( const SystemOfBodies& bodies, const Eigen::VectorXd& propagatedState, const double time )
{
    bodies.at( "Jupiter" )->setState( Eigen::Vector6d::Zero( ) );
    bodies.at( "Jupiter" )->setCurrentRotationalStateToLocalFrameFromEphemeris( time );

    bodies.at( "Io" )->setState( propagatedState.segment( 0, 6 ) );
    bodies.at( "Io" )->setCurrentRotationalStateToLocalFrame( propagatedState.segment( 6, 7 ) );
}

BOOST_AUTO_TEST_SUITE( test_full_two_body_coupled_dynamics_propagation )

BOOST_AUTO_TEST_CASE( testFullTwoBodyAccelerationAndTorqueInCoupledOrbitRotationPropagation )
{
    // Propagate Io's translational and rotational states together with the full two-body spherical-harmonic
    // acceleration and torque. The test verifies that the acceleration/torque values saved during propagation
    // are the same as values obtained by independently evaluating fresh model instances at each propagated state.
    // This checks the coupled propagation environment update, dependent-variable saving, and model consistency.
    spice_interface::loadStandardSpiceKernels( );

    const double initialEpoch = 1.0E8;
    const double integrationStep = 30.0;
    const double finalEpoch = initialEpoch + 5.0 * integrationStep;

    const std::vector< std::string > bodiesToPropagate = { "Io" };
    const std::vector< std::string > centralBodies = { "Jupiter" };
    const std::vector< std::string > bodiesToRotate = { "Io" };
    const std::vector< std::string > baseOrientationFrames = { "ECLIPJ2000" };

    SystemOfBodies bodies = createIoJupiterSystemWithDegreeTwoGravityFields( initialEpoch, finalEpoch );

    SelectedAccelerationMap selectedAccelerationMap;
    selectedAccelerationMap[ "Io" ][ "Jupiter" ].push_back( fullTwoBodySphericalHarmonicAcceleration( 2, 2, 2, 2 ) );
    AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodies, selectedAccelerationMap, bodiesToPropagate, centralBodies );

    SelectedTorqueMap selectedTorqueMap;
    selectedTorqueMap[ "Io" ][ "Jupiter" ].push_back( fullTwoBodySphericalHarmonicGravitationalTorque( 2, 2, 2, 2 ) );
    TorqueModelMap torqueModelMap = createTorqueModelsMap( bodies, selectedTorqueMap, bodiesToRotate );

    Eigen::Vector6d initialTranslationalState = Eigen::Vector6d::Zero( );
    const Eigen::Vector3d initialPosition =
            spice_interface::getBodyCartesianPositionAtEpoch( "Io", "Jupiter", "ECLIPJ2000", "NONE", initialEpoch );
    initialTranslationalState.segment( 0, 3 ) = initialPosition;
    Eigen::Vector3d velocityDirection = Eigen::Vector3d::UnitZ( ).cross( initialPosition );
    if( velocityDirection.norm( ) < std::numeric_limits< double >::epsilon( ) )
    {
        velocityDirection = Eigen::Vector3d::UnitX( ).cross( initialPosition );
    }
    velocityDirection.normalize( );
    initialTranslationalState.segment( 3, 3 ) = velocityDirection * std::sqrt( 1.266865349218008E17 / initialPosition.norm( ) );
    const Eigen::VectorXd initialRotationalState =
            getInitialRotationalStatesOfBodies< double, double >( bodiesToRotate, baseOrientationFrames, bodies, initialEpoch );

    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;
    dependentVariablesList.push_back( singleAccelerationDependentVariable( full_two_body_spherical_harmonic_gravity, "Io", "Jupiter" ) );
    dependentVariablesList.push_back( totalAccelerationDependentVariable( "Io" ) );
    dependentVariablesList.push_back( singleTorqueVariable( full_two_body_spherical_harmonic_gravitational_torque, "Io", "Jupiter" ) );

    std::shared_ptr< IntegratorSettings<> > integratorSettings =
            std::make_shared< IntegratorSettings<> >( rungeKutta4, initialEpoch, integrationStep );
    std::shared_ptr< PropagationTerminationSettings > terminationSettings =
            std::make_shared< PropagationTimeTerminationSettings >( finalEpoch );

    std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >( centralBodies,
                                                                                accelerationModelMap,
                                                                                bodiesToPropagate,
                                                                                initialTranslationalState,
                                                                                initialEpoch,
                                                                                integratorSettings,
                                                                                terminationSettings );

    std::shared_ptr< RotationalStatePropagatorSettings< double > > rotationalPropagatorSettings =
            std::make_shared< RotationalStatePropagatorSettings< double > >( torqueModelMap,
                                                                             bodiesToRotate,
                                                                             initialRotationalState,
                                                                             initialEpoch,
                                                                             integratorSettings,
                                                                             terminationSettings,
                                                                             quaternions );

    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsList;
    propagatorSettingsList.push_back( translationalPropagatorSettings );
    propagatorSettingsList.push_back( rotationalPropagatorSettings );

    std::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings = std::make_shared< MultiTypePropagatorSettings< double > >(
            propagatorSettingsList, integratorSettings, initialEpoch, terminationSettings, dependentVariablesList );

    SingleArcDynamicsSimulator< double > dynamicsSimulator( bodies, propagatorSettings );

    const std::map< double, Eigen::VectorXd >& stateHistory = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    const std::map< double, Eigen::VectorXd >& dependentVariableHistory = dynamicsSimulator.getDependentVariableHistory( );

    // Verify that propagation actually advanced beyond the initial state.
    BOOST_REQUIRE_GT( stateHistory.size( ), 1 );
    // Verify that a dependent-variable vector was saved for every propagated state epoch.
    BOOST_REQUIRE_EQUAL( stateHistory.size( ), dependentVariableHistory.size( ) );

    AccelerationMap independentAccelerationModelMap =
            createAccelerationModelsMap( bodies, selectedAccelerationMap, bodiesToPropagate, centralBodies );
    TorqueModelMap independentTorqueModelMap = createTorqueModelsMap( bodies, selectedTorqueMap, bodiesToRotate );

    std::shared_ptr< AccelerationModel< Eigen::Vector3d > > independentAccelerationModel =
            independentAccelerationModelMap.at( "Io" ).at( "Jupiter" ).at( 0 );
    std::shared_ptr< TorqueModel > independentTorqueModel = independentTorqueModelMap.at( "Io" ).at( "Jupiter" ).at( 0 );

    double maximumRelativeAccelerationDifference = 0.0;
    double maximumRelativeTorqueDifference = 0.0;
    double maximumAccelerationNorm = 0.0;
    double maximumTorqueNorm = 0.0;

    for( const auto& stateEntry : stateHistory )
    {
        const double currentTime = stateEntry.first;
        // Verify that the dependent-variable history contains the same epoch as the propagated state history.
        BOOST_REQUIRE( dependentVariableHistory.count( currentTime ) == 1 );

        setCoupledIoStateInEnvironment( bodies, stateEntry.second, currentTime );

        independentAccelerationModel->resetCurrentTime( );
        independentTorqueModel->resetCurrentTime( );
        independentAccelerationModel->updateMembers( currentTime );
        independentTorqueModel->updateMembers( currentTime );

        const Eigen::Vector3d savedAcceleration = dependentVariableHistory.at( currentTime ).segment( 0, 3 );
        const Eigen::Vector3d savedTorque = dependentVariableHistory.at( currentTime ).segment( 6, 3 );
        const Eigen::Vector3d independentlyEvaluatedAcceleration = independentAccelerationModel->getAcceleration( );
        const Eigen::Vector3d independentlyEvaluatedTorque = independentTorqueModel->getTorque( );

        maximumAccelerationNorm = std::max( maximumAccelerationNorm, savedAcceleration.norm( ) );
        maximumTorqueNorm = std::max( maximumTorqueNorm, savedTorque.norm( ) );

        const double accelerationScale = std::max(
                { savedAcceleration.norm( ), independentlyEvaluatedAcceleration.norm( ), std::numeric_limits< double >::epsilon( ) } );
        const double torqueScale =
                std::max( { savedTorque.norm( ), independentlyEvaluatedTorque.norm( ), std::numeric_limits< double >::epsilon( ) } );

        maximumRelativeAccelerationDifference =
                std::max( maximumRelativeAccelerationDifference,
                          ( savedAcceleration - independentlyEvaluatedAcceleration ).norm( ) / accelerationScale );
        maximumRelativeTorqueDifference =
                std::max( maximumRelativeTorqueDifference, ( savedTorque - independentlyEvaluatedTorque ).norm( ) / torqueScale );
    }

    // Verify that the saved acceleration signal is nonzero, so a trivial zero-output match cannot pass the test.
    BOOST_CHECK_GT( maximumAccelerationNorm, 0.0 );
    // Verify that the saved torque signal is nonzero, so a trivial zero-output match cannot pass the test.
    BOOST_CHECK_GT( maximumTorqueNorm, 0.0 );
    // Verify that the saved full two-body acceleration matches an independent model evaluation at the propagated states.
    BOOST_CHECK_SMALL( maximumRelativeAccelerationDifference, 1.0E-12 );
    // Verify that the saved full two-body torque matches an independent model evaluation at the propagated states.
    BOOST_CHECK_SMALL( maximumRelativeTorqueDifference, 1.0E-11 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
