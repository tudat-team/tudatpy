/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <boost/test/included/unit_test.hpp>

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

enum class PhobosGravityModel {
    mutualSphericalHarmonic,
    fullTwoBodySinglePointTerms,
    fullTwoBodyPhobosFigureSinglePointTerms,
    fullTwoBodyMarsFigureSinglePointTerms,
    fullTwoBodyDegreeThree
};

struct PhobosRotationSetup {
    SystemOfBodies bodies;
    std::shared_ptr< MultiTypePropagatorSettings< double > > propagatorSettings;
    std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate;
    std::shared_ptr< IntegratorSettings< double > > integratorSettings;
    AccelerationMap accelerationModelMap;
    TorqueModelMap torqueModelMap;
    Eigen::Matrix< double, 13, 1 > appliedInitialStateDifference;
};

PhobosRotationSetup createPhobosRotationSetup( const PhobosGravityModel gravityModel,
                                               const double finalEphemerisTime,
                                               const Eigen::Matrix< double, 13, 1 >& initialStateDifference,
                                               const bool estimateMarsGravityCoefficients,
                                               const double integratorStep )
{
    const double initialEphemerisTime = 0.0;
    Eigen::Matrix< double, 13, 1 > appliedInitialStateDifference = Eigen::Matrix< double, 13, 1 >::Zero( );

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
    if( gravityModel == PhobosGravityModel::fullTwoBodyDegreeThree )
    {
        phobosCosineGravityFieldCoefficients( 3, 0 ) = 2.0E-3;
        phobosCosineGravityFieldCoefficients( 3, 1 ) = -1.5E-3;
        phobosSineGravityFieldCoefficients( 3, 1 ) = 1.0E-3;
        phobosCosineGravityFieldCoefficients( 3, 2 ) = 0.8E-3;
        phobosSineGravityFieldCoefficients( 3, 2 ) = -0.7E-3;
        phobosCosineGravityFieldCoefficients( 3, 3 ) = 0.5E-3;
        phobosSineGravityFieldCoefficients( 3, 3 ) = 0.4E-3;
    }
    const std::shared_ptr< gravitation::SphericalHarmonicsGravityField > phobosGravityField =
            std::make_shared< gravitation::SphericalHarmonicsGravityField >( phobosGravitationalParameter,
                                                                             phobosReferenceRadius,
                                                                             phobosCosineGravityFieldCoefficients,
                                                                             phobosSineGravityFieldCoefficients,
                                                                             "Phobos_Fixed" );
    bodies.at( "Phobos" )->setGravityFieldModel( phobosGravityField );
    bodies.at( "Phobos" )
            ->setMassProperties( std::make_shared< simulation_setup::FromGravityFieldRigidBodyProperties >(
                    phobosGravityField, phobosScaledMeanMomentOfInertia ) );

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
    if( gravityModel == PhobosGravityModel::fullTwoBodyDegreeThree )
    {
        noRotationQuaternion = Eigen::Quaterniond::Identity( );
    }
    Eigen::Matrix< double, 7, 1 > unitRotationState = Eigen::Matrix< double, 7, 1 >::Zero( );
    unitRotationState( 0 ) = noRotationQuaternion.w( );
    unitRotationState( 1 ) = noRotationQuaternion.x( );
    unitRotationState( 2 ) = noRotationQuaternion.y( );
    unitRotationState( 3 ) = noRotationQuaternion.z( );
    unitRotationState( 4 ) = 1.0E-4;
    unitRotationState( 5 ) = -1.0E-4;
    unitRotationState( 6 ) = std::sqrt( marsGravitationalParameter / std::pow( phobosSemiMajorAxis, 3.0 ) ) + 1.0E-4;

    const Eigen::Matrix< double, 7, 1 > originalRotationState = unitRotationState;
    unitRotationState += initialStateDifference.segment( 6, 7 );
    unitRotationState( 0 ) = originalRotationState( 0 ) / std::fabs( originalRotationState( 0 ) ) *
            std::sqrt( 1.0 - std::pow( unitRotationState.segment( 1, 3 ).norm( ), 2.0 ) );
    appliedInitialStateDifference.segment( 0, 6 ) = initialStateDifference.segment( 0, 6 );
    appliedInitialStateDifference.segment( 6, 7 ) = unitRotationState - originalRotationState;

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
    else if( gravityModel == PhobosGravityModel::fullTwoBodySinglePointTerms ||
             gravityModel == PhobosGravityModel::fullTwoBodyPhobosFigureSinglePointTerms ||
             gravityModel == PhobosGravityModel::fullTwoBodyMarsFigureSinglePointTerms )
    {
        std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > fullTwoBodyInteractions =
                getExtendedSinglePointMassInteractions( 2, 2, 2, 2 );
        if( gravityModel == PhobosGravityModel::fullTwoBodyPhobosFigureSinglePointTerms )
        {
            fullTwoBodyInteractions.erase(
                    std::remove_if( fullTwoBodyInteractions.begin( ),
                                    fullTwoBodyInteractions.end( ),
                                    []( const std::tuple< unsigned int, unsigned int, unsigned int, unsigned int >& interaction ) {
                                        return std::get< 0 >( interaction ) == 0 && std::get< 1 >( interaction ) == 0;
                                    } ),
                    fullTwoBodyInteractions.end( ) );
        }
        else if( gravityModel == PhobosGravityModel::fullTwoBodyMarsFigureSinglePointTerms )
        {
            fullTwoBodyInteractions.erase(
                    std::remove_if( fullTwoBodyInteractions.begin( ),
                                    fullTwoBodyInteractions.end( ),
                                    []( const std::tuple< unsigned int, unsigned int, unsigned int, unsigned int >& interaction ) {
                                        return !( std::get< 0 >( interaction ) == 0 && std::get< 1 >( interaction ) == 0 ) ||
                                                ( std::get< 2 >( interaction ) == 0 && std::get< 3 >( interaction ) == 0 );
                                    } ),
                    fullTwoBodyInteractions.end( ) );
        }
        accelerationMap[ "Phobos" ][ "Mars" ].push_back(
                std::make_shared< FullTwoBodySphericalHarmonicAccelerationSettings >( fullTwoBodyInteractions ) );
        torqueMap[ "Phobos" ][ "Mars" ].push_back( fullTwoBodySphericalHarmonicGravitationalTorque( fullTwoBodyInteractions ) );
    }
    else
    {
        accelerationMap[ "Phobos" ][ "Mars" ].push_back( fullTwoBodySphericalHarmonicAcceleration( 3, 3, 3, 3 ) );
        torqueMap[ "Phobos" ][ "Mars" ].push_back( fullTwoBodySphericalHarmonicGravitationalTorque( 3, 3, 3, 3 ) );
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
    initialTranslationalState += initialStateDifference.segment( 0, 6 );
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
    const int maximumEstimatedGravityDegree = ( gravityModel == PhobosGravityModel::fullTwoBodyDegreeThree ) ? 3 : 2;
    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
            2, 0, maximumEstimatedGravityDegree, maximumEstimatedGravityDegree, "Phobos", spherical_harmonics_cosine_coefficient_block ) );
    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
            2, 1, maximumEstimatedGravityDegree, maximumEstimatedGravityDegree, "Phobos", spherical_harmonics_sine_coefficient_block ) );
    if( estimateMarsGravityCoefficients )
    {
        parameterNames.push_back(
                std::make_shared< SphericalHarmonicEstimatableParameterSettings >( 2,
                                                                                   0,
                                                                                   maximumEstimatedGravityDegree,
                                                                                   maximumEstimatedGravityDegree,
                                                                                   "Mars",
                                                                                   spherical_harmonics_cosine_coefficient_block ) );
        parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                2, 1, maximumEstimatedGravityDegree, maximumEstimatedGravityDegree, "Mars", spherical_harmonics_sine_coefficient_block ) );
    }

    return PhobosRotationSetup{ bodies,
                                propagatorSettings,
                                createParametersToEstimate( parameterNames, bodies ),
                                std::make_shared< IntegratorSettings< double > >( rungeKutta4, initialEphemerisTime, integratorStep ),
                                accelerationModelMap,
                                torqueModelMap,
                                appliedInitialStateDifference };
}

PhobosRotationSetup createPhobosRotationSetup( const PhobosGravityModel gravityModel, const double finalEphemerisTime )
{
    return createPhobosRotationSetup( gravityModel, finalEphemerisTime, Eigen::Matrix< double, 13, 1 >::Zero( ), false, 15.0 );
}

Eigen::Vector3d getSingleAcceleration( const AccelerationMap& accelerationModelMap )
{
    return accelerationModelMap.at( "Phobos" ).at( "Mars" ).at( 0 )->getAcceleration( );
}

Eigen::Vector3d getSingleTorque( const TorqueModelMap& torqueModelMap )
{
    return torqueModelMap.at( "Phobos" ).at( "Mars" ).at( 0 )->getTorque( );
}

struct InitialAccelerationPartialData {
    Eigen::Vector3d acceleration;
    Eigen::Vector3d torque;
    Eigen::MatrixXd stateAndVariationalRhs;
    Eigen::MatrixXd assembledAccelerationPartialWrtPhobosRotation;
    Eigen::MatrixXd directAccelerationPartialWrtPhobosRotation;
    Eigen::MatrixXd directAccelerationPartialWrtPhobosPosition;
    int parameterSetSize;
};

InitialAccelerationPartialData evaluateInitialAccelerationPartialData( const PhobosGravityModel gravityModel )
{
    PhobosRotationSetup setup = createPhobosRotationSetup( gravityModel, 15.0 );
    setup.propagatorSettings->resetInitialTime( setup.integratorSettings->initialTimeDeprecated_ );
    setup.propagatorSettings->setIntegratorSettings( setup.integratorSettings );
    setup.propagatorSettings->getOutputSettings( )->setClearNumericalSolutions( false );
    setup.propagatorSettings->getOutputSettings( )->setIntegratedResult( false );
    setup.propagatorSettings->getOutputSettings( )->setIntegratedVariationalResult( false );
    setup.propagatorSettings->makeOutputSettingsConsistent( );

    SingleArcVariationalEquationsSolver< double, double > solver(
            setup.bodies, setup.propagatorSettings, setup.parametersToEstimate, true, false );

    std::shared_ptr< DynamicsStateDerivativeModel< double, double > > derivativeModel =
            solver.getDynamicsSimulator( )->getDynamicsStateDerivative( );
    derivativeModel->setPropagationSettings( std::vector< IntegratedStateType >( ), true, true );

    const int stateSize = 13;
    const int parameterSetSize = setup.parametersToEstimate->getParameterSetSize( );
    Eigen::MatrixXd initialStateAndVariationalMatrix = Eigen::MatrixXd::Zero( stateSize, parameterSetSize + 1 );
    initialStateAndVariationalMatrix.block( 0, 0, stateSize, stateSize ).setIdentity( );
    initialStateAndVariationalMatrix.col( parameterSetSize ) =
            derivativeModel->convertFromOutputSolution( setup.propagatorSettings->getInitialStates( ), 0.0 );

    simulation_setup::setAreBodiesInPropagation( setup.bodies, true );
    const Eigen::MatrixXd stateAndVariationalRhs = derivativeModel->computeStateDerivative( 0.0, initialStateAndVariationalMatrix );

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

    Eigen::MatrixXd directAccelerationPartialWrtPhobosPosition = Eigen::MatrixXd::Zero( 3, 3 );
    accelerationPartial->wrtPositionOfAcceleratedBody( directAccelerationPartialWrtPhobosPosition.block( 0, 0, 3, 3 ), true );

    const Eigen::Vector3d acceleration = getSingleAcceleration( setup.accelerationModelMap );
    const Eigen::Vector3d torque = getSingleTorque( setup.torqueModelMap );
    simulation_setup::setAreBodiesInPropagation( setup.bodies, false );

    return InitialAccelerationPartialData{ acceleration,
                                           torque,
                                           stateAndVariationalRhs,
                                           stateAndVariationalRhs.block( 3, 6, 3, 7 ),
                                           directAccelerationPartialWrtPhobosRotation,
                                           directAccelerationPartialWrtPhobosPosition,
                                           parameterSetSize };
}

struct FullTwoBodyPropagationOutput {
    Eigen::VectorXd finalState;
    Eigen::MatrixXd combinedStateTransitionAndSensitivityMatrix;
    Eigen::Matrix< double, 13, 1 > appliedInitialStateDifference;
};

struct FullTwoBodyPropagationHistory {
    std::map< double, Eigen::VectorXd > stateHistory;
    std::map< double, Eigen::MatrixXd > combinedStateTransitionAndSensitivityMatrixHistory;
    Eigen::Matrix< double, 13, 1 > appliedInitialStateDifference;
};

FullTwoBodyPropagationHistory executeFullTwoBodyPhobosVariationalHistory(
        const Eigen::Matrix< double, 13, 1 >& initialStateDifference,
        const Eigen::VectorXd& parameterPerturbation = Eigen::VectorXd( ),
        const bool propagateVariationalEquations = true,
        const std::vector< double >& evaluationTimes = std::vector< double >{ 3600.0 },
        const PhobosGravityModel gravityModel = PhobosGravityModel::fullTwoBodyDegreeThree )
{
    const double finalEphemerisTime = evaluationTimes.back( );
    PhobosRotationSetup setup = createPhobosRotationSetup( gravityModel, finalEphemerisTime, initialStateDifference, true, 15.0 );

    Eigen::VectorXd parameterVector = setup.parametersToEstimate->getFullParameterValues< double >( );
    if( parameterPerturbation.rows( ) > 0 )
    {
        parameterVector.segment( 13, parameterPerturbation.rows( ) ) += parameterPerturbation;
        setup.parametersToEstimate->resetParameterValues( parameterVector );
    }

    SingleArcVariationalEquationsSolver< double, double > solver(
            setup.bodies, setup.integratorSettings, setup.propagatorSettings, setup.parametersToEstimate, true, nullptr, false, false );
    if( propagateVariationalEquations )
    {
        solver.integrateVariationalAndDynamicalEquations( setup.propagatorSettings->getInitialStates( ), true );
    }
    else
    {
        solver.integrateDynamicalEquationsOfMotionOnly( setup.propagatorSettings->getInitialStates( ) );
    }

    FullTwoBodyPropagationHistory history;
    for( double evaluationTime : evaluationTimes )
    {
        Eigen::VectorXd state = Eigen::VectorXd::Zero( 13 );
        state.segment( 0, 6 ) = setup.bodies.at( "Phobos" )->getEphemeris( )->getCartesianState( evaluationTime );
        state.segment( 6, 7 ) = setup.bodies.at( "Phobos" )->getRotationalEphemeris( )->getRotationStateVector( evaluationTime );
        history.stateHistory[ evaluationTime ] = state;

        if( propagateVariationalEquations )
        {
            history.combinedStateTransitionAndSensitivityMatrixHistory[ evaluationTime ] =
                    solver.getStateTransitionMatrixInterface( )->getCombinedStateTransitionAndSensitivityMatrix( evaluationTime );
        }
    }
    history.appliedInitialStateDifference = setup.appliedInitialStateDifference;
    return history;
}

FullTwoBodyPropagationOutput executeFullTwoBodyPhobosVariationalSimulation(
        const Eigen::Matrix< double, 13, 1 >& initialStateDifference,
        const Eigen::VectorXd& parameterPerturbation = Eigen::VectorXd( ),
        const bool propagateVariationalEquations = true,
        const PhobosGravityModel gravityModel = PhobosGravityModel::fullTwoBodyDegreeThree )
{
    const std::vector< double > evaluationTimes{ 3600.0 };
    const FullTwoBodyPropagationHistory history = executeFullTwoBodyPhobosVariationalHistory(
            initialStateDifference, parameterPerturbation, propagateVariationalEquations, evaluationTimes, gravityModel );
    return FullTwoBodyPropagationOutput{ history.stateHistory.at( evaluationTimes.back( ) ),
                                         propagateVariationalEquations
                                                 ? history.combinedStateTransitionAndSensitivityMatrixHistory.at( evaluationTimes.back( ) )
                                                 : Eigen::MatrixXd( ),
                                         history.appliedInitialStateDifference };
}

double getRelativeErrorForReferenceNorm( const double absoluteError, const double referenceNorm )
{
    if( referenceNorm > 0.0 )
    {
        return absoluteError / referenceNorm;
    }
    return absoluteError == 0.0 ? 0.0 : std::numeric_limits< double >::infinity( );
}

Eigen::Vector4d getInitialQuaternionVectorForModel( const PhobosGravityModel gravityModel )
{
    if( gravityModel == PhobosGravityModel::fullTwoBodyDegreeThree )
    {
        return Eigen::Vector4d::UnitX( );
    }

    const Eigen::Quaterniond initialQuaternion =
            Eigen::Quaterniond( Eigen::AngleAxisd( 1.0, Eigen::Vector3d::UnitZ( ) ) * Eigen::AngleAxisd( 2.0, Eigen::Vector3d::UnitX( ) ) *
                                Eigen::AngleAxisd( -0.5, Eigen::Vector3d::UnitZ( ) ) );
    return Eigen::Vector4d( initialQuaternion.w( ), initialQuaternion.x( ), initialQuaternion.y( ), initialQuaternion.z( ) );
}

Eigen::Matrix< double, 4, 3 > getPartialOfQuaternionWrtVectorPart( const Eigen::Vector4d& quaternion )
{
    BOOST_REQUIRE_GT( std::fabs( quaternion( 0 ) ), std::numeric_limits< double >::epsilon( ) );
    Eigen::Matrix< double, 4, 3 > partialOfQuaternionWrtVectorPart = Eigen::Matrix< double, 4, 3 >::Zero( );
    partialOfQuaternionWrtVectorPart.row( 0 ) = -quaternion.segment( 1, 3 ).transpose( ) / quaternion( 0 );
    partialOfQuaternionWrtVectorPart.block( 1, 0, 3, 3 ).setIdentity( );
    return partialOfQuaternionWrtVectorPart;
}

Eigen::MatrixXd projectQuaternionColumnsToVectorPart( const Eigen::MatrixXd& matrixWithFullQuaternionColumns,
                                                      const Eigen::Vector4d& quaternion )
{
    BOOST_REQUIRE_EQUAL( matrixWithFullQuaternionColumns.cols( ), 4 );
    return matrixWithFullQuaternionColumns * getPartialOfQuaternionWrtVectorPart( quaternion );
}

void checkColumnsCloseRelative( const Eigen::MatrixXd& computed,
                                const Eigen::MatrixXd& expected,
                                const double relativeTolerance,
                                const std::string& checkIdentifier,
                                const int firstColumnIndex = 0 )
{
    BOOST_REQUIRE_EQUAL( computed.rows( ), expected.rows( ) );
    BOOST_REQUIRE_EQUAL( computed.cols( ), expected.cols( ) );

    for( int column = 0; column < computed.cols( ); ++column )
    {
        const Eigen::VectorXd computedColumn = computed.col( column );
        const Eigen::VectorXd expectedColumn = expected.col( column );
        const double errorNorm = ( computedColumn - expectedColumn ).norm( );
        const double expectedColumnNorm = expectedColumn.norm( );
        const double absoluteTolerance = relativeTolerance * expectedColumnNorm;
        const double relativeError = getRelativeErrorForReferenceNorm( errorNorm, expectedColumnNorm );
        std::ostringstream failureMessage;
        failureMessage << std::setprecision( 17 ) << "Failed " << checkIdentifier << ", column " << firstColumnIndex + column
                       << "\ncomputed value: " << computedColumn.transpose( ) << "\nexpected value: " << expectedColumn.transpose( )
                       << "\nabsolute tolerance: " << absoluteTolerance << "\nrelative tolerance: " << relativeTolerance
                       << "\nabsolute error: " << errorNorm << "\nrelative error: " << relativeError;
        BOOST_TEST_CONTEXT( "column " << firstColumnIndex + column )
        {
            BOOST_CHECK_MESSAGE( errorNorm <= absoluteTolerance, failureMessage.str( ) );
        }
    }
}

struct StateVectorBlock {
    std::string name;
    int startRow;
    int numberOfRows;
};

void checkStateResponseColumnsCloseRelative( const Eigen::MatrixXd& computed,
                                             const Eigen::MatrixXd& expected,
                                             const double relativeTolerance,
                                             const std::string& checkIdentifier,
                                             const int firstColumnIndex = 0 )
{
    const int stateSize = 13;
    BOOST_REQUIRE_EQUAL( computed.rows( ), stateSize );
    BOOST_REQUIRE_EQUAL( expected.rows( ), stateSize );
    BOOST_REQUIRE_EQUAL( computed.cols( ), expected.cols( ) );

    const std::vector< StateVectorBlock > stateVectorBlocks{
        { "position", 0, 3 }, { "velocity", 3, 3 }, { "quaternion", 6, 4 }, { "angular velocity", 10, 3 }
    };

    for( int column = 0; column < computed.cols( ); ++column )
    {
        for( const StateVectorBlock& block : stateVectorBlocks )
        {
            const Eigen::VectorXd computedColumn = computed.block( block.startRow, column, block.numberOfRows, 1 );
            const Eigen::VectorXd expectedColumn = expected.block( block.startRow, column, block.numberOfRows, 1 );
            const double errorNorm = ( computedColumn - expectedColumn ).norm( );
            const double expectedColumnNorm = expectedColumn.norm( );
            const double absoluteTolerance = relativeTolerance * expectedColumnNorm;
            const double relativeError = getRelativeErrorForReferenceNorm( errorNorm, expectedColumnNorm );
            std::ostringstream failureMessage;
            failureMessage << std::setprecision( 17 ) << "Failed " << checkIdentifier << ", " << block.name << ", column "
                           << firstColumnIndex + column << "\ncomputed value: " << computedColumn.transpose( )
                           << "\nexpected value: " << expectedColumn.transpose( ) << "\nabsolute tolerance: " << absoluteTolerance
                           << "\nrelative tolerance: " << relativeTolerance << "\nabsolute error: " << errorNorm
                           << "\nrelative error: " << relativeError;
            BOOST_TEST_CONTEXT( block.name << ", column " << firstColumnIndex + column )
            {
                BOOST_CHECK_MESSAGE( errorNorm <= absoluteTolerance, failureMessage.str( ) );
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( testRestrictedFullTwoBodyInitialAccelerationPartials )
{
    // Compare the restricted full two-body single-point interaction set against the equivalent legacy mutual
    // spherical-harmonic acceleration/torque setup at the initial state only. This is not a propagation test:
    // it directly checks the acceleration, torque, and the assembled acceleration partials that enter the
    // variational-equation RHS.
    spice_interface::loadStandardSpiceKernels( );

    const InitialAccelerationPartialData mutualData = evaluateInitialAccelerationPartialData( PhobosGravityModel::mutualSphericalHarmonic );
    const InitialAccelerationPartialData fullTwoBodyData =
            evaluateInitialAccelerationPartialData( PhobosGravityModel::fullTwoBodySinglePointTerms );
    const InitialAccelerationPartialData fullTwoBodyPhobosFigureData =
            evaluateInitialAccelerationPartialData( PhobosGravityModel::fullTwoBodyPhobosFigureSinglePointTerms );
    const InitialAccelerationPartialData fullTwoBodyMarsFigureData =
            evaluateInitialAccelerationPartialData( PhobosGravityModel::fullTwoBodyMarsFigureSinglePointTerms );

    // Verify that both formulations expose the same estimatable-parameter vector to the variational equations.
    BOOST_REQUIRE_EQUAL( fullTwoBodyData.parameterSetSize, mutualData.parameterSetSize );

    // Verify physical equivalence of the restricted acceleration and torque before inspecting their partials.
    BOOST_CHECK_SMALL( ( fullTwoBodyData.acceleration - mutualData.acceleration ).norm( ), 1.0E-14 );
    BOOST_CHECK_SMALL( ( fullTwoBodyData.torque - mutualData.torque ).norm( ) / mutualData.torque.norm( ), 1.0E-9 );

    // Verify that the acceleration partial object contributes the same Phobos-rotation block as is assembled
    // into the complete state/variational-equation RHS.
    checkColumnsCloseRelative( fullTwoBodyData.assembledAccelerationPartialWrtPhobosRotation,
                               fullTwoBodyData.directAccelerationPartialWrtPhobosRotation,
                               1.0E-12,
                               "restricted full two-body acceleration partial assembly wrt Phobos rotation" );
    // Verify translational-state equivalence against the legacy mutual spherical-harmonic formulation.
    checkColumnsCloseRelative( fullTwoBodyData.directAccelerationPartialWrtPhobosPosition,
                               mutualData.stateAndVariationalRhs.block( 3, 0, 3, 3 ),
                               1.0E-12,
                               "restricted full two-body acceleration partial wrt Phobos position vs legacy mutual model" );
    // Verify additivity of the two restricted single-point interaction groups in the rotation partial.
    checkColumnsCloseRelative( fullTwoBodyPhobosFigureData.directAccelerationPartialWrtPhobosRotation +
                                       fullTwoBodyMarsFigureData.directAccelerationPartialWrtPhobosRotation,
                               fullTwoBodyData.directAccelerationPartialWrtPhobosRotation,
                               1.0E-12,
                               "restricted full two-body split acceleration partials wrt Phobos rotation vs combined model" );
    const Eigen::Vector4d initialQuaternion = getInitialQuaternionVectorForModel( PhobosGravityModel::mutualSphericalHarmonic );
    // Verify that the Mars-figure-only restricted contribution has no physical dependence on the constrained
    // Phobos quaternion-vector perturbation. The projection removes the non-unique raw quaternion column.
    checkColumnsCloseRelative(
            projectQuaternionColumnsToVectorPart( fullTwoBodyMarsFigureData.directAccelerationPartialWrtPhobosRotation.block( 0, 0, 3, 4 ),
                                                  initialQuaternion ),
            Eigen::MatrixXd::Zero( 3, 3 ),
            1.0E-12,
            "restricted full two-body Mars-figure acceleration partial wrt constrained Phobos quaternion vector should be zero",
            7 );
    // Verify that the same Mars-figure-only restricted contribution has no angular-velocity dependence.
    checkColumnsCloseRelative( fullTwoBodyMarsFigureData.directAccelerationPartialWrtPhobosRotation.block( 0, 4, 3, 3 ),
                               Eigen::MatrixXd::Zero( 3, 3 ),
                               1.0E-12,
                               "restricted full two-body Mars-figure acceleration partial wrt Phobos angular velocity should be zero",
                               10 );
    // Verify the physically meaningful constrained quaternion-vector partial against the legacy formulation.
    checkColumnsCloseRelative(
            projectQuaternionColumnsToVectorPart( fullTwoBodyData.directAccelerationPartialWrtPhobosRotation.block( 0, 0, 3, 4 ),
                                                  initialQuaternion ),
            projectQuaternionColumnsToVectorPart( mutualData.assembledAccelerationPartialWrtPhobosRotation.block( 0, 0, 3, 4 ),
                                                  initialQuaternion ),
            1.0E-12,
            "restricted full two-body acceleration partial wrt constrained Phobos quaternion vector vs legacy mutual model",
            7 );
    // Verify angular-velocity partial equivalence against the legacy formulation.
    checkColumnsCloseRelative( fullTwoBodyData.directAccelerationPartialWrtPhobosRotation.block( 0, 4, 3, 3 ),
                               mutualData.assembledAccelerationPartialWrtPhobosRotation.block( 0, 4, 3, 3 ),
                               1.0E-12,
                               "restricted full two-body acceleration partial wrt Phobos angular velocity vs legacy mutual model",
                               10 );
}

BOOST_AUTO_TEST_CASE( testMarsPhobosFullTwoBodyVariationalEquationsFiniteDifference )
{
    // Validate the propagated full two-body translational-rotational STM and selected sensitivities by direct finite differences.
    spice_interface::loadStandardSpiceKernels( );

    const int stateSize = 13;
    const int sensitivityParameterSize = 25;
    const int translationalStateStartColumn = 0;
    const int rotationalStateStartColumn = 6;
    const FullTwoBodyPropagationOutput nominalOutput = executeFullTwoBodyPhobosVariationalSimulation(
            Eigen::Matrix< double, stateSize, 1 >::Zero( ), Eigen::VectorXd::Zero( sensitivityParameterSize ), true );
    const Eigen::MatrixXd& analyticalMatrix = nominalOutput.combinedStateTransitionAndSensitivityMatrix;
    BOOST_CHECK_EQUAL( analyticalMatrix.rows( ), stateSize );
    BOOST_CHECK_EQUAL( analyticalMatrix.cols( ), stateSize + sensitivityParameterSize );

    Eigen::Matrix< double, stateSize, 1 > statePerturbation;
    statePerturbation << 10.0, 10.0, 10.0, 1.0E-2, 1.0E-2, 1.0E-2, 0.0, 3.0E-6, 3.0E-6, 3.0E-6, 1.0E-5, 1.0E-5, 1.0E-5;
    Eigen::MatrixXd numericalMatrix = Eigen::MatrixXd::Zero( stateSize, stateSize + sensitivityParameterSize );

    for( int column = 0; column < 6; ++column )
    {
        Eigen::Matrix< double, stateSize, 1 > perturbation = Eigen::Matrix< double, stateSize, 1 >::Zero( );
        perturbation( column ) = statePerturbation( column );
        const Eigen::VectorXd upState =
                executeFullTwoBodyPhobosVariationalSimulation( perturbation, Eigen::VectorXd::Zero( sensitivityParameterSize ), false )
                        .finalState;
        const Eigen::VectorXd downState =
                executeFullTwoBodyPhobosVariationalSimulation( -perturbation, Eigen::VectorXd::Zero( sensitivityParameterSize ), false )
                        .finalState;
        // Check how the full two-body variational equations map initial translational-state perturbations into the final coupled state.
        numericalMatrix.block( 0, translationalStateStartColumn + column, stateSize, 1 ) =
                ( upState - downState ) / ( 2.0 * statePerturbation( column ) );
    }

    for( int column = 10; column < stateSize; ++column )
    {
        Eigen::Matrix< double, stateSize, 1 > perturbation = Eigen::Matrix< double, stateSize, 1 >::Zero( );
        perturbation( column ) = statePerturbation( column );
        const Eigen::VectorXd upState =
                executeFullTwoBodyPhobosVariationalSimulation( perturbation, Eigen::VectorXd::Zero( sensitivityParameterSize ), false )
                        .finalState;
        const Eigen::VectorXd downState =
                executeFullTwoBodyPhobosVariationalSimulation( -perturbation, Eigen::VectorXd::Zero( sensitivityParameterSize ), false )
                        .finalState;
        // Check how the full two-body variational equations map initial angular-velocity perturbations into the final coupled state.
        numericalMatrix.block( 0, rotationalStateStartColumn + column - 6, stateSize, 1 ) =
                ( upState - downState ) / ( 2.0 * statePerturbation( column ) );
    }

    for( int column = 7; column < 10; ++column )
    {
        Eigen::Matrix< double, stateSize, 1 > perturbation = Eigen::Matrix< double, stateSize, 1 >::Zero( );
        perturbation( column ) = statePerturbation( column );
        const FullTwoBodyPropagationOutput upOutput =
                executeFullTwoBodyPhobosVariationalSimulation( perturbation, Eigen::VectorXd::Zero( sensitivityParameterSize ), false );
        const FullTwoBodyPropagationOutput downOutput =
                executeFullTwoBodyPhobosVariationalSimulation( -perturbation, Eigen::VectorXd::Zero( sensitivityParameterSize ), false );
        const Eigen::VectorXd appliedStateDifference = upOutput.appliedInitialStateDifference - downOutput.appliedInitialStateDifference;
        // Raw quaternion columns are not unique. Check the physically applied constrained quaternion-vector
        // perturbation through the complete STM action instead of comparing individual raw quaternion columns.
        checkStateResponseColumnsCloseRelative( ( analyticalMatrix.block( 0, 0, stateSize, stateSize ) * appliedStateDifference ),
                                                ( upOutput.finalState - downOutput.finalState ),
                                                1.0E-2,
                                                "constrained quaternion-vector STM action",
                                                column );
    }

    Eigen::VectorXd parameterPerturbation = Eigen::VectorXd::Constant( sensitivityParameterSize, 1.0E-4 );
    for( int column = 0; column < sensitivityParameterSize; ++column )
    {
        Eigen::VectorXd perturbation = Eigen::VectorXd::Zero( sensitivityParameterSize );
        perturbation( column ) = parameterPerturbation( column );
        const Eigen::VectorXd upState =
                executeFullTwoBodyPhobosVariationalSimulation( Eigen::Matrix< double, stateSize, 1 >::Zero( ), perturbation, false )
                        .finalState;
        const Eigen::VectorXd downState =
                executeFullTwoBodyPhobosVariationalSimulation( Eigen::Matrix< double, stateSize, 1 >::Zero( ), -perturbation, false )
                        .finalState;
        // Check sensitivities to Phobos mean moment of inertia and degree-2/3 gravity-field coefficients of Phobos and Mars.
        numericalMatrix.block( 0, stateSize + column, stateSize, 1 ) = ( upState - downState ) / ( 2.0 * parameterPerturbation( column ) );
    }

    // Verify the translational-state STM columns against central finite differences.
    checkStateResponseColumnsCloseRelative( numericalMatrix.block( 0, translationalStateStartColumn, stateSize, 6 ),
                                            analyticalMatrix.block( 0, translationalStateStartColumn, stateSize, 6 ),
                                            5.0E-4,
                                            "translational initial-state STM finite difference",
                                            translationalStateStartColumn );
    // Verify the angular-velocity STM columns against central finite differences.
    checkStateResponseColumnsCloseRelative( numericalMatrix.block( 0, rotationalStateStartColumn + 4, stateSize, 3 ),
                                            analyticalMatrix.block( 0, rotationalStateStartColumn + 4, stateSize, 3 ),
                                            1.0E-2,
                                            "angular-velocity initial-state STM finite difference",
                                            rotationalStateStartColumn + 4 );
    // Verify selected gravity and inertia sensitivity columns against central finite differences.
    checkStateResponseColumnsCloseRelative( numericalMatrix.block( 0, stateSize, stateSize, sensitivityParameterSize ),
                                            analyticalMatrix.block( 0, stateSize, stateSize, sensitivityParameterSize ),
                                            2.0E-2,
                                            "sensitivity finite difference",
                                            stateSize );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
