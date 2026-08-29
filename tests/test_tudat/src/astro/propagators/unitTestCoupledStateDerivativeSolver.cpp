/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license.
 */

#define BOOST_TEST_MAIN

#include <boost/test/included/unit_test.hpp>

#include <Eigen/Core>

#include "tudat/astro/propagators/coupledStateDerivativeSolver.h"
#include "tudat/astro/propagators/dynamicsStateDerivativeModel.h"
#include "tudat/basics/testMacros.h"

namespace tudat
{
namespace unit_tests
{

using namespace propagators;

namespace
{

class TestCoupledStateDerivative : public SingleStateTypeDerivative< double, double >
{
public:
    TestCoupledStateDerivative( const IntegratedStateType stateType,
                                const int stateSize,
                                const std::function< Eigen::VectorXd( ) >& derivativeFunction ):
        SingleStateTypeDerivative< double, double >( stateType, { "Body" } ), stateSize_( stateSize ),
        derivativeFunction_( derivativeFunction )
    {}

    void calculateSystemStateDerivative( const double, const Eigen::VectorXd&, Eigen::Block< Eigen::MatrixXd > stateDerivative ) override
    {
        stateDerivative = derivativeFunction_( );
    }

    void clearStateDerivativeModel( ) override {}

    void updateStateDerivativeModel( const double ) override {}

    void convertCurrentStateToGlobalRepresentation( const Eigen::VectorXd& internalSolution,
                                                    const double&,
                                                    Eigen::Block< Eigen::VectorXd > currentState ) override
    {
        currentState = internalSolution;
    }

    Eigen::MatrixXd convertFromOutputSolution( const Eigen::MatrixXd& outputSolution, const double& ) override
    {
        return outputSolution;
    }

    void convertToOutputSolution( const Eigen::MatrixXd& internalSolution,
                                  const double&,
                                  Eigen::Block< Eigen::VectorXd > currentState ) override
    {
        currentState = internalSolution;
    }

    int getConventionalStateSize( ) override
    {
        return stateSize_;
    }

private:
    int stateSize_;
    std::function< Eigen::VectorXd( ) > derivativeFunction_;
};

}  // namespace

BOOST_AUTO_TEST_CASE( testDirectScaledAffineCoupledDerivativeSolution )
{
    const Eigen::Vector2d componentScales( 1.0e-12, 1.0e-3 );
    Eigen::Matrix2d scaledCouplingMatrix;
    scaledCouplingMatrix << 0.1, 0.2, -0.3, 0.05;
    const Eigen::Matrix2d couplingMatrix =
            componentScales.asDiagonal( ) * scaledCouplingMatrix * componentScales.cwiseInverse( ).asDiagonal( );
    const Eigen::Vector2d constantTerm( 2.0e-15, -4.0e-4 );
    const auto derivativeMapping = [ = ]( const Eigen::VectorXd& derivative ) { return constantTerm + couplingMatrix * derivative; };

    const CoupledStateDerivativeSolution< double > solution =
            solveCoupledStateDerivative< double >( derivativeMapping,
                                                   Eigen::Vector2d::Zero( ),
                                                   componentScales,
                                                   CoupledStateDerivativeSolverSettings( true, 1.0e-12, 1.0e-12, 10 ) );

    const Eigen::Vector2d expectedDerivative = ( Eigen::Matrix2d::Identity( ) - couplingMatrix ).inverse( ) * constantTerm;
    const Eigen::Matrix2d expectedMultiplier = ( Eigen::Matrix2d::Identity( ) - couplingMatrix ).inverse( );
    BOOST_CHECK( solution.usedDirectSolution_ );
    BOOST_CHECK( solution.converged_ );
    BOOST_CHECK_EQUAL( solution.numberOfIterations_, 0 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( solution.stateDerivative_, expectedDerivative, 2.0e-14 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( solution.implicitDerivativeMultiplier_, expectedMultiplier, 2.0e-14 );
}

BOOST_AUTO_TEST_CASE( testConfiguredIterativeCoupledDerivativeSolution )
{
    Eigen::Matrix2d couplingMatrix;
    couplingMatrix << 0.2, 0.0, 0.0, -0.1;
    const Eigen::Vector2d constantTerm( 1.0e-9, -2.0e-6 );
    const auto derivativeMapping = [ = ]( const Eigen::VectorXd& derivative ) { return constantTerm + couplingMatrix * derivative; };
    const Eigen::Vector2d componentScales( 1.0e-9, 1.0e-6 );

    const CoupledStateDerivativeSolution< double > solution =
            solveCoupledStateDerivative< double >( derivativeMapping,
                                                   Eigen::Vector2d::Zero( ),
                                                   componentScales,
                                                   CoupledStateDerivativeSolverSettings( false, 1.0e-13, 1.0e-13, 30 ) );

    const Eigen::Vector2d expectedDerivative = ( Eigen::Matrix2d::Identity( ) - couplingMatrix ).inverse( ) * constantTerm;
    BOOST_CHECK( !solution.usedDirectSolution_ );
    BOOST_CHECK( solution.converged_ );
    BOOST_CHECK_GT( solution.numberOfIterations_, 1 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( solution.stateDerivative_, expectedDerivative, 2.0e-13 );
}

BOOST_AUTO_TEST_CASE( testCoupledDerivativeNonConvergencePolicies )
{
    const auto divergentMapping = []( const Eigen::VectorXd& derivative ) { return ( derivative.array( ) + 1.0 ).matrix( ).eval( ); };
    const Eigen::VectorXd initialDerivative = Eigen::VectorXd::Zero( 1 );
    const Eigen::VectorXd componentScales = Eigen::VectorXd::Ones( 1 );

    BOOST_CHECK_THROW( solveCoupledStateDerivative< double >(
                               divergentMapping,
                               initialDerivative,
                               componentScales,
                               CoupledStateDerivativeSolverSettings( false, 0.0, 0.0, 3, throw_exception_on_coupled_derivative_failure ) ),
                       std::runtime_error );

    const CoupledStateDerivativeSolution< double > acceptedSolution = solveCoupledStateDerivative< double >(
            divergentMapping,
            initialDerivative,
            componentScales,
            CoupledStateDerivativeSolverSettings( false, 0.0, 0.0, 3, accept_last_coupled_derivative_iteration ) );
    BOOST_CHECK( !acceptedSolution.converged_ );
    BOOST_CHECK_EQUAL( acceptedSolution.numberOfIterations_, 3 );
    BOOST_CHECK_EQUAL( acceptedSolution.stateDerivative_( 0 ), 3.0 );
}

BOOST_AUTO_TEST_CASE( testCoupledDerivativeDynamicsModelIntegration )
{
    Eigen::Vector5d gravityRateInput = Eigen::Vector5d::Zero( );
    Eigen::Vector3d angularAccelerationInput = Eigen::Vector3d::Zero( );
    const Eigen::Vector5d uncoupledGravityRate = ( Eigen::Vector5d( ) << 2.0e-12, -3.0e-12, 1.0e-12, 4.0e-12, -2.0e-12 ).finished( );
    const Eigen::Vector3d uncoupledAngularAcceleration( 4.0e-7, -2.0e-7, 1.0e-7 );
    Eigen::Matrix< double, 5, 3 > gravityFromAcceleration;
    gravityFromAcceleration << 2.0e-6, 0.0, 0.0, 0.0, -1.0e-6, 0.0, 0.0, 0.0, 3.0e-6, 1.0e-6, 0.0, 0.0, 0.0, 2.0e-6, -1.0e-6;
    Eigen::Matrix< double, 3, 5 > accelerationFromGravity;
    accelerationFromGravity << 2.0e4, 0.0, -1.0e4, 0.0, 0.0, 0.0, -3.0e4, 0.0, 1.0e4, 0.0, 0.0, 0.0, 2.0e4, 0.0, -1.0e4;

    const std::shared_ptr< TestCoupledStateDerivative > gravityModel =
            std::make_shared< TestCoupledStateDerivative >( gravity_deformation_state, 5, [ & ]( ) {
                return ( uncoupledGravityRate + gravityFromAcceleration * angularAccelerationInput ).eval( );
            } );
    const Eigen::Vector4d quaternionDerivative( 0.1, 0.2, 0.3, 0.4 );
    const std::shared_ptr< TestCoupledStateDerivative > rotationModel =
            std::make_shared< TestCoupledStateDerivative >( rotational_state, 7, [ & ]( ) {
                Eigen::VectorXd derivative( 7 );
                derivative << quaternionDerivative, uncoupledAngularAcceleration + accelerationFromGravity * gravityRateInput;
                return derivative;
            } );

    std::unordered_map< IntegratedStateType, std::vector< std::shared_ptr< SingleStateTypeDerivative< double, double > > > > modelsToUpdate;
    modelsToUpdate[ gravity_deformation_state ] = { gravityModel };
    modelsToUpdate[ rotational_state ] = { rotationModel };
    std::map< StateDerivativeDependency, std::vector< std::string > > updateSettings = {
        { inertia_tensor_derivative_dependency, { "Body" } }, { rotation_rate_derivative_dependency, { "Body" } }
    };
    using EnvironmentFunction = std::function< void( const Eigen::VectorXd& ) >;
    std::map< StateDerivativeDependency, std::vector< std::pair< std::string, EnvironmentFunction > > > environmentFunctions;
    environmentFunctions[ inertia_tensor_derivative_dependency ].push_back(
            { "Body", [ & ]( const Eigen::VectorXd& derivative ) { gravityRateInput = derivative; } } );
    environmentFunctions[ rotation_rate_derivative_dependency ].push_back(
            { "Body", [ & ]( const Eigen::VectorXd& derivative ) { angularAccelerationInput = derivative; } } );

    const std::shared_ptr< StateDerivativeUpdater< double, double > > updater =
            std::make_shared< StateDerivativeUpdater< double, double > >( simulation_setup::SystemOfBodies( ),
                                                                          updateSettings,
                                                                          std::vector< std::function< void( double ) > >( ),
                                                                          modelsToUpdate,
                                                                          environmentFunctions );
    DynamicsStateDerivativeModel< double, double > dynamicsModel(
            { gravityModel, rotationModel },
            []( const double,
                const std::unordered_map< IntegratedStateType, Eigen::VectorXd >&,
                const std::vector< IntegratedStateType > ) {},
            nullptr,
            updater );
    dynamicsModel.setPropagationSettings( {}, true, false );

    const Eigen::Vector3d expectedAngularAcceleration =
            ( Eigen::Matrix3d::Identity( ) - accelerationFromGravity * gravityFromAcceleration ).inverse( ) *
            ( uncoupledAngularAcceleration + accelerationFromGravity * uncoupledGravityRate );
    const Eigen::Vector5d expectedGravityRate = uncoupledGravityRate + gravityFromAcceleration * expectedAngularAcceleration;
    const Eigen::MatrixXd derivative = dynamicsModel.computeStateDerivative( 0.0, Eigen::VectorXd::Zero( 12 ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( derivative.block( 0, 0, 5, 1 ), expectedGravityRate, 5.0e-13 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( derivative.block( 5, 0, 4, 1 ), quaternionDerivative, 5.0e-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( derivative.block( 9, 0, 3, 1 ), expectedAngularAcceleration, 5.0e-13 );

    Eigen::MatrixXd derivativeBeforeModifier;
    dynamicsModel.setStateDerivativeModifierFunction( [ & ]( const Eigen::MatrixXd&, const Eigen::MatrixXd& physicalDerivative ) {
        derivativeBeforeModifier = physicalDerivative;
        return 2.0 * physicalDerivative;
    } );
    const Eigen::MatrixXd modifiedDerivative = dynamicsModel.computeStateDerivative( 0.0, Eigen::VectorXd::Zero( 12 ) );
    const Eigen::MatrixXd expectedModifiedDerivative = 2.0 * derivative;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( derivativeBeforeModifier, derivative, 5.0e-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( modifiedDerivative, expectedModifiedDerivative, 5.0e-15 );
}

BOOST_AUTO_TEST_CASE( testOneWayDerivativeDependencyIsNotTreatedAsFeedback )
{
    Eigen::Vector3d angularAccelerationInput = Eigen::Vector3d::Zero( );
    const Eigen::Vector5d uncoupledGravityRate = Eigen::Vector5d::Constant( 1.0e-12 );
    const Eigen::Vector3d angularAcceleration( 2.0e-7, -3.0e-7, 4.0e-7 );
    const Eigen::Matrix< double, 5, 3 > gravityFromAcceleration = ( Eigen::Matrix< double, 5, 3 >( ) << 1.0e-6,
                                                                    0.0,
                                                                    0.0,
                                                                    0.0,
                                                                    2.0e-6,
                                                                    0.0,
                                                                    0.0,
                                                                    0.0,
                                                                    -1.0e-6,
                                                                    1.0e-6,
                                                                    1.0e-6,
                                                                    0.0,
                                                                    0.0,
                                                                    -1.0e-6,
                                                                    1.0e-6 )
                                                                          .finished( );

    const std::shared_ptr< TestCoupledStateDerivative > gravityModel =
            std::make_shared< TestCoupledStateDerivative >( gravity_deformation_state, 5, [ & ]( ) {
                return ( uncoupledGravityRate + gravityFromAcceleration * angularAccelerationInput ).eval( );
            } );
    const std::shared_ptr< TestCoupledStateDerivative > rotationModel =
            std::make_shared< TestCoupledStateDerivative >( rotational_state, 7, [ & ]( ) {
                Eigen::VectorXd derivative = Eigen::VectorXd::Zero( 7 );
                derivative.tail( 3 ) = angularAcceleration;
                return derivative;
            } );

    std::unordered_map< IntegratedStateType, std::vector< std::shared_ptr< SingleStateTypeDerivative< double, double > > > > modelsToUpdate;
    modelsToUpdate[ gravity_deformation_state ] = { gravityModel };
    std::map< StateDerivativeDependency, std::vector< std::string > > updateSettings = { { rotation_rate_derivative_dependency,
                                                                                           { "Body" } } };
    using EnvironmentFunction = std::function< void( const Eigen::VectorXd& ) >;
    std::map< StateDerivativeDependency, std::vector< std::pair< std::string, EnvironmentFunction > > > environmentFunctions;
    environmentFunctions[ rotation_rate_derivative_dependency ].push_back(
            { "Body", [ & ]( const Eigen::VectorXd& derivative ) { angularAccelerationInput = derivative; } } );
    const std::shared_ptr< StateDerivativeUpdater< double, double > > updater =
            std::make_shared< StateDerivativeUpdater< double, double > >( simulation_setup::SystemOfBodies( ),
                                                                          updateSettings,
                                                                          std::vector< std::function< void( double ) > >( ),
                                                                          modelsToUpdate,
                                                                          environmentFunctions );
    DynamicsStateDerivativeModel< double, double > dynamicsModel(
            { gravityModel, rotationModel },
            []( const double,
                const std::unordered_map< IntegratedStateType, Eigen::VectorXd >&,
                const std::vector< IntegratedStateType > ) {},
            nullptr,
            updater );
    dynamicsModel.setPropagationSettings( {}, true, false );

    const Eigen::MatrixXd derivative = dynamicsModel.computeStateDerivative( 0.0, Eigen::VectorXd::Zero( 12 ) );
    const Eigen::Vector5d expectedGravityRate = uncoupledGravityRate + gravityFromAcceleration * angularAcceleration;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( derivative.block( 0, 0, 5, 1 ), expectedGravityRate, 5.0e-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( derivative.block( 9, 0, 3, 1 ), angularAcceleration, 5.0e-15 );
}

}  // namespace unit_tests
}  // namespace tudat
