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

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <boost/test/unit_test.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/ephemerides/keplerEphemeris.h"
#include "tudat/astro/gravitation/sphericalHarmonicsGravityField.h"
#include "tudat/basics/testMacros.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/math/interpolators/linearInterpolator.h"
#include "tudat/simulation/environment_setup/createGravityField.h"
#include "tudat/simulation/environment_setup/createRotationModel.h"
#include "tudat/simulation/estimation_setup/createAccelerationPartials.h"
#include "tudat/simulation/estimation_setup/createEstimatableParametersFactory.h"
#include "tudat/simulation/estimation_setup/singleArcVariationalEquationsSolver.h"
#include "tudat/simulation/propagation_setup/createAccelerationModels.h"
#include "tudat/simulation/propagation_setup/createTorqueModel.h"
#include "tudat/simulation/propagation_setup/propagationSettings.h"

namespace tudat
{
namespace unit_tests
{

using namespace basic_astrodynamics;
using namespace estimatable_parameters;
using namespace interpolators;
using namespace numerical_integrators;
using namespace propagators;
using namespace simulation_setup;

BOOST_AUTO_TEST_SUITE( test_full_two_body_variational_equations )

enum class PhobosGravityModel { mutualSphericalHarmonic, fullTwoBodySinglePointTerms };

struct PhobosRotationSetup {
    SystemOfBodies bodies;
    std::shared_ptr< MultiTypePropagatorSettings< double > > propagatorSettings;
    std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate;
    std::shared_ptr< IntegratorSettings< double > > integratorSettings;
    AccelerationMap accelerationModelMap;
    TorqueModelMap torqueModelMap;
};

PhobosRotationSetup createPhobosRotationSetup( const PhobosGravityModel gravityModel, const double finalEphemerisTime )
{
    const double initialEphemerisTime = 0.0;

    SystemOfBodies bodies( "Mars", "ECLIPJ2000" );
    bodies.createEmptyBody( "Mars", false );
    bodies.at( "Mars" )->setEphemeris(
            std::make_shared< ephemerides::ConstantEphemeris >( [ = ]( ) { return Eigen::Vector6d::Zero( ); } ) );
    bodies.at( "Mars" )->setRotationalEphemeris(
            createRotationModel( getDefaultRotationModelSettings( "Mars", initialEphemerisTime, finalEphemerisTime ), "Mars" ) );
    bodies.at( "Mars" )->setGravityFieldModel(
            createGravityFieldModel( getDefaultGravityFieldSettings( "Mars", initialEphemerisTime, finalEphemerisTime ), "Mars", bodies ) );
    const double marsGravitationalParameter = bodies.at( "Mars" )->getGravityFieldModel( )->getGravitationalParameter( );

    bodies.createEmptyBody( "Phobos" );
    Eigen::Matrix3d phobosInertiaTensor = Eigen::Matrix3d::Zero( );
    phobosInertiaTensor( 0, 0 ) = 0.3615;
    phobosInertiaTensor( 1, 1 ) = 0.4265;
    phobosInertiaTensor( 2, 2 ) = 0.5024;
    phobosInertiaTensor *= ( 11.27E3 * 11.27E3 * 1.0659E16 );

    const double phobosGravitationalParameter = 1.0659E16 * physical_constants::GRAVITATIONAL_CONSTANT;
    const double phobosReferenceRadius = 11.27E3;
    Eigen::MatrixXd phobosCosineGravityFieldCoefficients = Eigen::MatrixXd::Zero( 6, 6 );
    Eigen::MatrixXd phobosSineGravityFieldCoefficients = Eigen::MatrixXd::Zero( 6, 6 );
    double phobosScaledMeanMomentOfInertia = TUDAT_NAN;
    gravitation::getDegreeTwoSphericalHarmonicCoefficients( phobosInertiaTensor,
                                                            phobosGravitationalParameter,
                                                            phobosReferenceRadius,
                                                            true,
                                                            phobosCosineGravityFieldCoefficients,
                                                            phobosSineGravityFieldCoefficients,
                                                            phobosScaledMeanMomentOfInertia );
    bodies.at( "Phobos" )
            ->setGravityFieldModel( std::make_shared< gravitation::SphericalHarmonicsGravityField >( phobosGravitationalParameter,
                                                                                                     phobosReferenceRadius,
                                                                                                     phobosCosineGravityFieldCoefficients,
                                                                                                     phobosSineGravityFieldCoefficients,
                                                                                                     "Phobos_Fixed",
                                                                                                     phobosScaledMeanMomentOfInertia ) );

    Eigen::Vector6d phobosKeplerElements = Eigen::Vector6d::Zero( );
    const double phobosSemiMajorAxis = 9376.0E3;
    phobosKeplerElements( 0 ) = phobosSemiMajorAxis;
    phobosKeplerElements( 2 ) = 0.1;
    bodies.at( "Phobos" )
            ->setEphemeris( ephemerides::getTabulatedEphemeris(
                    std::make_shared< ephemerides::KeplerEphemeris >( phobosKeplerElements, 0.0, marsGravitationalParameter ),
                    initialEphemerisTime - 3600.0,
                    finalEphemerisTime + 3600.0,
                    60.0 ) );

    Eigen::Quaterniond noRotationQuaternion =
            Eigen::Quaterniond( Eigen::AngleAxisd( 1.0, Eigen::Vector3d::UnitZ( ) ) * Eigen::AngleAxisd( 2.0, Eigen::Vector3d::UnitX( ) ) *
                                Eigen::AngleAxisd( -0.5, Eigen::Vector3d::UnitZ( ) ) );
    Eigen::Matrix< double, 7, 1 > unitRotationState = Eigen::Matrix< double, 7, 1 >::Zero( );
    unitRotationState( 0 ) = noRotationQuaternion.w( );
    unitRotationState( 1 ) = noRotationQuaternion.x( );
    unitRotationState( 2 ) = noRotationQuaternion.y( );
    unitRotationState( 3 ) = noRotationQuaternion.z( );
    unitRotationState( 4 ) = 1.0E-4;
    unitRotationState( 5 ) = -1.0E-4;
    unitRotationState( 6 ) = std::sqrt( marsGravitationalParameter / std::pow( phobosSemiMajorAxis, 3.0 ) ) + 1.0E-4;

    std::map< double, Eigen::Matrix< double, 7, 1 > > dummyRotationMap;
    dummyRotationMap[ -1.0E100 ] = unitRotationState;
    dummyRotationMap[ 1.0E100 ] = unitRotationState;
    bodies.at( "Phobos" )
            ->setRotationalEphemeris( std::make_shared< ephemerides::TabulatedRotationalEphemeris< double, double > >(
                    std::make_shared< LinearInterpolator< double, Eigen::Matrix< double, 7, 1 > > >( dummyRotationMap ),
                    "ECLIPJ2000",
                    "Phobos_Fixed" ) );

    SelectedAccelerationMap accelerationMap;
    SelectedTorqueMap torqueMap;
    if( gravityModel == PhobosGravityModel::mutualSphericalHarmonic )
    {
        accelerationMap[ "Phobos" ][ "Mars" ].push_back( std::make_shared< MutualSphericalHarmonicAccelerationSettings >( 2, 2, 2, 2 ) );
        torqueMap[ "Phobos" ][ "Mars" ].push_back( std::make_shared< TorqueSettings >( second_order_gravitational_torque ) );
    }
    else
    {
        accelerationMap[ "Phobos" ][ "Mars" ].push_back( std::make_shared< FullTwoBodySphericalHarmonicAccelerationSettings >(
                getExtendedSinglePointMassInteractions( 2, 2, 2, 2 ) ) );
        torqueMap[ "Phobos" ][ "Mars" ].push_back( fullTwoBodySphericalHarmonicGravitationalTorque( 2, 2, 0, 0 ) );
    }

    const std::vector< std::string > translationalBodiesToIntegrate{ "Phobos" };
    const std::vector< std::string > translationalCentralBodies{ "Mars" };
    AccelerationMap accelerationModelMap =
            createAccelerationModelsMap( bodies, accelerationMap, translationalBodiesToIntegrate, translationalCentralBodies );
    TorqueModelMap torqueModelMap = createTorqueModelsMap( bodies, torqueMap, translationalBodiesToIntegrate );

    std::shared_ptr< PropagationTerminationSettings > terminationSettings =
            std::make_shared< PropagationTimeTerminationSettings >( finalEphemerisTime );
    std::shared_ptr< RotationalStatePropagatorSettings< double > > rotationalPropagatorSettings =
            std::make_shared< RotationalStatePropagatorSettings< double > >(
                    torqueModelMap, translationalBodiesToIntegrate, unitRotationState, terminationSettings );

    Eigen::VectorXd initialTranslationalState =
            getInitialStatesOfBodies( translationalBodiesToIntegrate, translationalCentralBodies, bodies, initialEphemerisTime );
    std::shared_ptr< TranslationalStatePropagatorSettings<> > translationalPropagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings<> >( translationalCentralBodies,
                                                                        accelerationModelMap,
                                                                        translationalBodiesToIntegrate,
                                                                        initialTranslationalState,
                                                                        finalEphemerisTime,
                                                                        cowell );
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsList{ rotationalPropagatorSettings,
                                                                                                    translationalPropagatorSettings };
    std::shared_ptr< MultiTypePropagatorSettings< double > > propagatorSettings =
            std::make_shared< MultiTypePropagatorSettings< double > >( propagatorSettingsList, terminationSettings );

    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back(
            std::make_shared< InitialRotationalStateEstimatableParameterSettings< double > >( "Phobos", unitRotationState, "ECLIPJ2000" ) );
    parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
            "Phobos", initialTranslationalState, "Mars" ) );
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Phobos", mean_moment_of_inertia ) );
    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
            1, 0, 2, 2, "Phobos", spherical_harmonics_cosine_coefficient_block ) );
    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
            2, 1, 2, 2, "Phobos", spherical_harmonics_sine_coefficient_block ) );

    return PhobosRotationSetup{ bodies,
                                propagatorSettings,
                                createParametersToEstimate( parameterNames, bodies ),
                                std::make_shared< IntegratorSettings< double > >( rungeKutta4, initialEphemerisTime, 15.0 ),
                                accelerationModelMap,
                                torqueModelMap };
}

Eigen::Vector3d getSingleAcceleration( const AccelerationMap& accelerationModelMap )
{
    return accelerationModelMap.at( "Phobos" ).at( "Mars" ).at( 0 )->getAcceleration( );
}

Eigen::Vector3d getSingleTorque( const TorqueModelMap& torqueModelMap )
{
    return torqueModelMap.at( "Phobos" ).at( "Mars" ).at( 0 )->getTorque( );
}

struct InitialDerivativeDiagnostics {
    Eigen::Vector3d acceleration;
    Eigen::Vector3d torque;
    Eigen::MatrixXd variationalRhs;
    Eigen::MatrixXd stateAndVariationalRhs;
    Eigen::MatrixXd directAccelerationPartialWrtPhobosRotation;
};

InitialDerivativeDiagnostics evaluateInitialDerivativeDiagnostics( const PhobosGravityModel gravityModel )
{
    PhobosRotationSetup setup = createPhobosRotationSetup( gravityModel, 15.0 );
    SingleArcVariationalEquationsSolver< double, double > solver(
            setup.bodies, setup.integratorSettings, setup.propagatorSettings, setup.parametersToEstimate, true, nullptr, false, false );

    std::shared_ptr< DynamicsStateDerivativeModel< double, double > > derivativeModel =
            solver.getDynamicsSimulator( )->getDynamicsStateDerivative( );
    derivativeModel->setPropagationSettings( std::vector< IntegratedStateType >( ), true, true );

    Eigen::MatrixXd initialState = Eigen::MatrixXd::Zero( 13, setup.parametersToEstimate->getParameterSetSize( ) + 1 );
    initialState.block( 0, 0, 13, 13 ).setIdentity( );
    initialState.col( setup.parametersToEstimate->getParameterSetSize( ) ) =
            derivativeModel->convertFromOutputSolution( setup.propagatorSettings->getInitialStates( ), 0.0 );

    const Eigen::MatrixXd stateAndVariationalRhs = derivativeModel->computeStateDerivative( 0.0, initialState );

    std::shared_ptr< acceleration_partials::AccelerationPartial > accelerationPartial =
            createAnalyticalAccelerationPartial( setup.accelerationModelMap.at( "Phobos" ).at( "Mars" ).at( 0 ),
                                                 std::make_pair( "Phobos", setup.bodies.at( "Phobos" ) ),
                                                 std::make_pair( "Mars", setup.bodies.at( "Mars" ) ),
                                                 setup.bodies,
                                                 setup.parametersToEstimate );
    accelerationPartial->update( 0.0 );
    Eigen::MatrixXd directAccelerationPartialWrtPhobosRotation = Eigen::MatrixXd::Zero( 3, 7 );
    accelerationPartial->wrtNonTranslationalStateOfAdditionalBody(
            directAccelerationPartialWrtPhobosRotation.block( 0, 0, 3, 7 ), std::make_pair( "Phobos", "" ), rotational_state, true );

    return InitialDerivativeDiagnostics{ getSingleAcceleration( setup.accelerationModelMap ),
                                         getSingleTorque( setup.torqueModelMap ),
                                         stateAndVariationalRhs.block( 0, 0, 13, setup.parametersToEstimate->getParameterSetSize( ) ),
                                         stateAndVariationalRhs,
                                         directAccelerationPartialWrtPhobosRotation };
}

BOOST_AUTO_TEST_CASE( testFullTwoBodyEquivalentPhobosRotationVariationalEquationDerivative )
{
    // This test isolates the existing Phobos coupled translational-rotational variational setup at a single time.
    // Both configurations use the same environment construction path; only the acceleration/torque settings are switched.
    spice_interface::loadStandardSpiceKernels( );

    const InitialDerivativeDiagnostics mutualDiagnostics =
            evaluateInitialDerivativeDiagnostics( PhobosGravityModel::mutualSphericalHarmonic );
    const InitialDerivativeDiagnostics fullTwoBodyDiagnostics =
            evaluateInitialDerivativeDiagnostics( PhobosGravityModel::fullTwoBodySinglePointTerms );

    // Verify that the restricted full two-body force model reproduces the legacy mutual spherical-harmonic acceleration value.
    BOOST_CHECK_SMALL( ( fullTwoBodyDiagnostics.acceleration - mutualDiagnostics.acceleration ).norm( ), 1.0E-14 );
    // Verify that the restricted full two-body torque model reproduces the legacy second-degree gravitational torque value.
    BOOST_CHECK_SMALL( ( fullTwoBodyDiagnostics.torque - mutualDiagnostics.torque ).norm( ) / mutualDiagnostics.torque.norm( ), 1.0E-9 );

    const Eigen::MatrixXd assembledAccelerationPartialWrtPhobosRotation = fullTwoBodyDiagnostics.variationalRhs.block( 3, 6, 3, 7 );
    // Verify that the full two-body acceleration partial is inserted into the translational variational-equation rows.
    BOOST_CHECK_SMALL( ( assembledAccelerationPartialWrtPhobosRotation - fullTwoBodyDiagnostics.directAccelerationPartialWrtPhobosRotation )
                               .cwiseAbs( )
                               .maxCoeff( ),
                       1.0E-14 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
