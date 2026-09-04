/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <limits>
#include <sstream>
#include <string>

#include <boost/test/included/unit_test.hpp>

#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/ephemerides/simpleRotationalEphemeris.h"
#include "tudat/astro/gravitation/periodicGravityFieldVariations.h"
#include "tudat/astro/gravitation/polynomialGravityFieldVariations.h"
#include "tudat/astro/orbit_determination/acceleration_partials/numericalAccelerationPartial.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/constantRotationalOrientation.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/constantRotationRate.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameterSet.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/gravityFieldVariationParameters.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/gravitationalParameter.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/initialMassState.h"
#include "tudat/astro/orbit_determination/observation_partials/rotationMatrixPartial.h"
#include "tudat/basics/testMacros.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/rigidBodyProperties.h"
#include "tudat/simulation/estimation_setup/createEstimatableParametersFactory.h"
#include "tudat/simulation/estimation_setup/createTorquePartials.h"
#include "tudat/simulation/propagation_setup/createTorqueModel.h"

namespace tudat
{

namespace unit_tests
{

using namespace tudat::simulation_setup;
using namespace tudat::estimatable_parameters;
using namespace tudat::acceleration_partials;
using namespace tudat::basic_astrodynamics;

namespace
{

observation_partials::RotationMatrixPartialNamedList createSimpleRotationPartialMap( const std::shared_ptr< Body >& body )
{
    std::shared_ptr< ephemerides::SimpleRotationalEphemeris > rotationModel =
            std::dynamic_pointer_cast< ephemerides::SimpleRotationalEphemeris >( body->getRotationalEphemeris( ) );
    observation_partials::RotationMatrixPartialNamedList partials;
    partials[ std::make_pair( constant_rotation_rate, "" ) ] =
            std::make_shared< observation_partials::RotationMatrixPartialWrtConstantRotationRate >( rotationModel );
    partials[ std::make_pair( rotation_pole_position, "" ) ] =
            std::make_shared< observation_partials::RotationMatrixPartialWrtPoleOrientation >( rotationModel );
    return partials;
}

SystemOfBodies createTwoBodyTorquePartialTestSystem( const double testTime,
                                                     const std::string& bodyUndergoingTorqueName,
                                                     const std::string& bodyExertingTorqueName,
                                                     const bool useArbitraryRotationStates,
                                                     const bool useZeroScaledMeanMomentOfInertia = false )
{
    SystemOfBodies bodies = SystemOfBodies( "SSB", "ECLIPJ2000" );
    bodies.createEmptyBody( bodyUndergoingTorqueName );
    bodies.createEmptyBody( bodyExertingTorqueName );

    std::shared_ptr< Body > bodyUndergoingTorque = bodies.at( bodyUndergoingTorqueName );
    std::shared_ptr< Body > bodyExertingTorque = bodies.at( bodyExertingTorqueName );

    double scaledMeanMomentOfInertiaBodyUndergoingTorque = useZeroScaledMeanMomentOfInertia ? 0.0 : 0.3;
    double scaledMeanMomentOfInertiaBodyExertingTorque = useZeroScaledMeanMomentOfInertia ? 0.0 : 0.26;
    const double massBodyUndergoingTorque = 6.0E11;
    const double massBodyExertingTorque = 4.0E11;
    const double gravitationalParameterBodyUndergoingTorque = massBodyUndergoingTorque * physical_constants::GRAVITATIONAL_CONSTANT;
    const double gravitationalParameterBodyExertingTorque = massBodyExertingTorque * physical_constants::GRAVITATIONAL_CONSTANT;
    const double referenceRadiusBodyUndergoingTorque = 520.0;
    const double referenceRadiusBodyExertingTorque = 430.0;

    Eigen::Matrix3d normalizedInertiaTensorOfBodyUndergoingTorque = Eigen::Matrix3d::Zero( );
    normalizedInertiaTensorOfBodyUndergoingTorque( 0, 0 ) = 0.30;
    normalizedInertiaTensorOfBodyUndergoingTorque( 1, 1 ) = 0.33;
    normalizedInertiaTensorOfBodyUndergoingTorque( 2, 2 ) = 0.37;
    normalizedInertiaTensorOfBodyUndergoingTorque( 0, 1 ) = 0.01;
    normalizedInertiaTensorOfBodyUndergoingTorque( 1, 0 ) = normalizedInertiaTensorOfBodyUndergoingTorque( 0, 1 );
    normalizedInertiaTensorOfBodyUndergoingTorque( 0, 2 ) = -0.008;
    normalizedInertiaTensorOfBodyUndergoingTorque( 2, 0 ) = normalizedInertiaTensorOfBodyUndergoingTorque( 0, 2 );

    Eigen::Matrix3d normalizedInertiaTensorOfBodyExertingTorque = Eigen::Matrix3d::Zero( );
    normalizedInertiaTensorOfBodyExertingTorque( 0, 0 ) = 0.28;
    normalizedInertiaTensorOfBodyExertingTorque( 1, 1 ) = 0.35;
    normalizedInertiaTensorOfBodyExertingTorque( 2, 2 ) = 0.37;
    normalizedInertiaTensorOfBodyExertingTorque( 0, 1 ) = -0.012;
    normalizedInertiaTensorOfBodyExertingTorque( 1, 0 ) = normalizedInertiaTensorOfBodyExertingTorque( 0, 1 );
    normalizedInertiaTensorOfBodyExertingTorque( 1, 2 ) = 0.009;
    normalizedInertiaTensorOfBodyExertingTorque( 2, 1 ) = normalizedInertiaTensorOfBodyExertingTorque( 1, 2 );

    const Eigen::Matrix3d inertiaTensorOfBodyUndergoingTorque = normalizedInertiaTensorOfBodyUndergoingTorque * massBodyUndergoingTorque *
            referenceRadiusBodyUndergoingTorque * referenceRadiusBodyUndergoingTorque;
    const Eigen::Matrix3d inertiaTensorOfBodyExertingTorque = normalizedInertiaTensorOfBodyExertingTorque * massBodyExertingTorque *
            referenceRadiusBodyExertingTorque * referenceRadiusBodyExertingTorque;

    Eigen::MatrixXd cosineCoefficientsOfBodyUndergoingTorque = Eigen::MatrixXd::Zero( 6, 6 );
    Eigen::MatrixXd sineCoefficientsOfBodyUndergoingTorque = Eigen::MatrixXd::Zero( 6, 6 );
    Eigen::MatrixXd cosineCoefficientsOfBodyExertingTorque = Eigen::MatrixXd::Zero( 6, 6 );
    Eigen::MatrixXd sineCoefficientsOfBodyExertingTorque = Eigen::MatrixXd::Zero( 6, 6 );

    gravitation::getDegreeTwoSphericalHarmonicCoefficients( inertiaTensorOfBodyUndergoingTorque,
                                                            gravitationalParameterBodyUndergoingTorque,
                                                            referenceRadiusBodyUndergoingTorque,
                                                            true,
                                                            cosineCoefficientsOfBodyUndergoingTorque,
                                                            sineCoefficientsOfBodyUndergoingTorque,
                                                            scaledMeanMomentOfInertiaBodyUndergoingTorque );

    gravitation::getDegreeTwoSphericalHarmonicCoefficients( inertiaTensorOfBodyExertingTorque,
                                                            gravitationalParameterBodyExertingTorque,
                                                            referenceRadiusBodyExertingTorque,
                                                            true,
                                                            cosineCoefficientsOfBodyExertingTorque,
                                                            sineCoefficientsOfBodyExertingTorque,
                                                            scaledMeanMomentOfInertiaBodyExertingTorque );

    bodyUndergoingTorque->setGravityFieldModel(
            std::make_shared< gravitation::SphericalHarmonicsGravityField >( gravitationalParameterBodyUndergoingTorque,
                                                                             referenceRadiusBodyUndergoingTorque,
                                                                             cosineCoefficientsOfBodyUndergoingTorque,
                                                                             sineCoefficientsOfBodyUndergoingTorque,
                                                                             bodyUndergoingTorqueName + "_Fixed",
                                                                             scaledMeanMomentOfInertiaBodyUndergoingTorque ) );
    bodyExertingTorque->setGravityFieldModel(
            std::make_shared< gravitation::SphericalHarmonicsGravityField >( gravitationalParameterBodyExertingTorque,
                                                                             referenceRadiusBodyExertingTorque,
                                                                             cosineCoefficientsOfBodyExertingTorque,
                                                                             sineCoefficientsOfBodyExertingTorque,
                                                                             bodyExertingTorqueName + "_Fixed",
                                                                             scaledMeanMomentOfInertiaBodyExertingTorque ) );

    bodyUndergoingTorque->getMassProperties( )->update( testTime );
    bodyExertingTorque->getMassProperties( )->update( testTime );

    Eigen::Matrix< double, 7, 1 > unitRotationState = Eigen::Matrix< double, 7, 1 >::Zero( );
    unitRotationState( 0 ) = 1.0;
    std::map< double, Eigen::Matrix< double, 7, 1 > > dummyRotationMap;
    dummyRotationMap[ -1.0E100 ] = unitRotationState;
    dummyRotationMap[ 1.0E100 ] = unitRotationState;
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< double, 7, 1 > > > dummyInterpolator =
            std::make_shared< interpolators::LinearInterpolator< double, Eigen::Matrix< double, 7, 1 > > >( dummyRotationMap );
    bodyUndergoingTorque->setRotationalEphemeris( std::make_shared< ephemerides::TabulatedRotationalEphemeris< double, double > >(
            dummyInterpolator, "ECLIPJ2000", bodyUndergoingTorqueName + "_Fixed" ) );
    bodyExertingTorque->setRotationalEphemeris( std::make_shared< ephemerides::TabulatedRotationalEphemeris< double, double > >(
            dummyInterpolator, "ECLIPJ2000", bodyExertingTorqueName + "_Fixed" ) );

    Eigen::Vector6d translationalStateOfBodyUndergoingTorque = Eigen::Vector6d::Zero( );
    translationalStateOfBodyUndergoingTorque << 120.0, -45.0, 300.0, 0.02, -0.03, 0.01;
    Eigen::Vector6d translationalStateOfBodyExertingTorque = Eigen::Vector6d::Zero( );
    translationalStateOfBodyExertingTorque << 3780.0, 2510.0, -1640.0, -0.01, 0.04, -0.02;
    bodyUndergoingTorque->setState( translationalStateOfBodyUndergoingTorque );
    bodyExertingTorque->setState( translationalStateOfBodyExertingTorque );

    if( useArbitraryRotationStates )
    {
        Eigen::Vector7d arbitraryRotationStateOfBodyUndergoingTorque = Eigen::Vector7d::Zero( );
        arbitraryRotationStateOfBodyUndergoingTorque.segment( 0, 4 ) =
                tudat::linear_algebra::convertQuaternionToVectorFormat( Eigen::Quaterniond(
                        Eigen::AngleAxisd( 0.6, Eigen::Vector3d::UnitX( ) ) * Eigen::AngleAxisd( -0.25, Eigen::Vector3d::UnitY( ) ) *
                        Eigen::AngleAxisd( 0.4, Eigen::Vector3d::UnitZ( ) ) ) );
        Eigen::Vector7d arbitraryRotationStateOfBodyExertingTorque = Eigen::Vector7d::Zero( );
        arbitraryRotationStateOfBodyExertingTorque.segment( 0, 4 ) =
                tudat::linear_algebra::convertQuaternionToVectorFormat( Eigen::Quaterniond(
                        Eigen::AngleAxisd( -0.2, Eigen::Vector3d::UnitX( ) ) * Eigen::AngleAxisd( 0.45, Eigen::Vector3d::UnitY( ) ) *
                        Eigen::AngleAxisd( -0.5, Eigen::Vector3d::UnitZ( ) ) ) );

        bodyUndergoingTorque->setCurrentRotationalStateToLocalFrame( arbitraryRotationStateOfBodyUndergoingTorque );
        bodyExertingTorque->setCurrentRotationalStateToLocalFrame( arbitraryRotationStateOfBodyExertingTorque );
    }
    else
    {
        bodyUndergoingTorque->setCurrentRotationalStateToLocalFrame( unitRotationState );
        bodyExertingTorque->setCurrentRotationalStateToLocalFrame( unitRotationState );
    }

    return bodies;
}

void checkMatrixClosePerElement( const Eigen::MatrixXd& analyticalMatrix,
                                 const Eigen::MatrixXd& numericalMatrix,
                                 const double tolerance,
                                 const std::string& label = "" )
{
    // Require equal dimensions so every analytical partial entry is compared to the corresponding numerical entry.
    BOOST_REQUIRE_EQUAL( analyticalMatrix.rows( ), numericalMatrix.rows( ) );
    // Require equal dimensions so no coefficient/state partial column is silently skipped.
    BOOST_REQUIRE_EQUAL( analyticalMatrix.cols( ), numericalMatrix.cols( ) );

    for( int i = 0; i < analyticalMatrix.rows( ); i++ )
    {
        for( int j = 0; j < analyticalMatrix.cols( ); j++ )
        {
            const double relativeDifference =
                    std::fabs( analyticalMatrix( i, j ) - numericalMatrix( i, j ) ) / std::max( 1.0, std::fabs( numericalMatrix( i, j ) ) );
            // Verify the analytical partial entry agrees with its central finite-difference reference.
            BOOST_CHECK_MESSAGE( relativeDifference < tolerance,
                                 "Torque partial mismatch"
                                         << ( label.empty( ) ? "" : " in " + label ) << " at (" << i << "," << j
                                         << "): analytical=" << analyticalMatrix( i, j ) << ", numerical=" << numericalMatrix( i, j )
                                         << ", relativeDifference=" << relativeDifference << ", tolerance=" << tolerance );
        }
    }
}

Eigen::Matrix< double, 6, 1 > getAuxiliaryFunctionVector(
        const acceleration_partials::detail::FourthDegreeTorqueAuxiliaryQuantities& auxiliaryQuantities )
{
    Eigen::Matrix< double, 6, 1 > auxiliaryFunctionVector;
    auxiliaryFunctionVector << auxiliaryQuantities.fyzFunction, auxiliaryQuantities.fxzFunction, auxiliaryQuantities.fxyFunction,
            auxiliaryQuantities.gyzFunction, auxiliaryQuantities.gxzFunction, auxiliaryQuantities.gxyFunction;
    return auxiliaryFunctionVector;
}

void updateBodyMassDistributions( const std::shared_ptr< Body >& bodyUndergoingTorque,
                                  const std::shared_ptr< Body >& bodyExertingTorque,
                                  const double testTime )
{
    bodyUndergoingTorque->getMassProperties( )->resetCurrentTime( );
    bodyExertingTorque->getMassProperties( )->resetCurrentTime( );
    bodyUndergoingTorque->updateMassDistribution( testTime );
    bodyExertingTorque->updateMassDistribution( testTime );
}

void applyGravityFieldVariation( const std::shared_ptr< gravitation::SphericalHarmonicsGravityField >& gravityField,
                                 const Eigen::MatrixXd& baseCosineCoefficients,
                                 const Eigen::MatrixXd& baseSineCoefficients,
                                 const std::shared_ptr< gravitation::GravityFieldVariations >& variationModel,
                                 const double testTime )
{
    Eigen::MatrixXd cosineCoefficients = baseCosineCoefficients;
    Eigen::MatrixXd sineCoefficients = baseSineCoefficients;
    variationModel->addSphericalHarmonicsCorrections( testTime, sineCoefficients, cosineCoefficients );
    gravityField->setCosineCoefficients( cosineCoefficients );
    gravityField->setSineCoefficients( sineCoefficients );
}

void makeBodyMassIndependentOfGravityField( const std::shared_ptr< Body >& body )
{
    std::ostringstream discardedWarning;
    std::streambuf* originalWarningBuffer = std::cerr.rdbuf( discardedWarning.rdbuf( ) );
    body->setMassProperties( std::make_shared< TimeDependentRigidBodyProperties >(
            body->getBodyMass( ), body->getBodyFixedCenterOfMass( ), body->getBodyInertiaTensor( ) ) );
    std::cerr.rdbuf( originalWarningBuffer );
}

std::shared_ptr< InitialMassStateParameter< double > > makeInitialMassParameter( const std::shared_ptr< Body >& body,
                                                                                 const std::string& bodyName )
{
    Eigen::VectorXd massState( 1 );
    massState( 0 ) = body->getBodyMass( );
    std::shared_ptr< InitialMassStateParameter< double > > massParameter =
            std::make_shared< InitialMassStateParameter< double > >( bodyName, massState );
    massParameter->addStateClosureFunctions(
            [ body ]( ) {
                Eigen::VectorXd currentMass( 1 );
                currentMass( 0 ) = body->getBodyMass( );
                return currentMass;
            },
            [ body ]( const Eigen::VectorXd& newMass ) { body->setCurrentPropagatedBodyMass( newMass( 0 ) ); } );
    return massParameter;
}

}  // namespace

BOOST_AUTO_TEST_SUITE( test_full_two_body_torque_partials )

BOOST_AUTO_TEST_CASE( testFourthDegreeTorqueAuxiliaryFunctionPartials )
{
    // Test rationale:
    // Verify partial derivatives of the fourth-degree torque auxiliary functions (f/g terms and contracted
    // inertia tensor contributions used in Schutz-style closed-form torque expressions) against central
    // finite differences. This is a low-level unit check that isolates the algebraic building blocks
    // before they are used inside the full torque-partial model.
    const Eigen::Vector3d relativePositionInBodyFixedFrame( 3275.0, -1840.0, 2510.0 );
    const double massOfBodyExertingTorque = 4.0E11;
    const Eigen::Matrix3d inertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque =
            ( Eigen::Matrix3d( ) << 2.6E16, 1.2E15, -9.0E14, 1.2E15, 3.1E16, -1.6E15, -9.0E14, -1.6E15, 3.5E16 ).finished( );
    const Eigen::Matrix< double, 6, 1 > independentInertiaTensorComponentsOfBodyExertingTorque =
            acceleration_partials::detail::getIndependentInertiaTensorComponentsFromMatrix(
                    inertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque );
    const acceleration_partials::detail::FourthDegreeTorqueAuxiliaryQuantities auxiliaryQuantities =
            acceleration_partials::detail::computeFourthDegreeTorqueAuxiliaryQuantities(
                    relativePositionInBodyFixedFrame, massOfBodyExertingTorque, independentInertiaTensorComponentsOfBodyExertingTorque );

    const double positionPerturbation = 1.0E-3;
    for( int coordinateIndex = 0; coordinateIndex < 3; coordinateIndex++ )
    {
        Eigen::Vector3d upPerturbedPosition = relativePositionInBodyFixedFrame;
        upPerturbedPosition( coordinateIndex ) += positionPerturbation;
        const auto upPerturbedAuxiliaryQuantities = acceleration_partials::detail::computeFourthDegreeTorqueAuxiliaryQuantities(
                upPerturbedPosition, massOfBodyExertingTorque, independentInertiaTensorComponentsOfBodyExertingTorque );

        Eigen::Vector3d downPerturbedPosition = relativePositionInBodyFixedFrame;
        downPerturbedPosition( coordinateIndex ) -= positionPerturbation;
        const auto downPerturbedAuxiliaryQuantities = acceleration_partials::detail::computeFourthDegreeTorqueAuxiliaryQuantities(
                downPerturbedPosition, massOfBodyExertingTorque, independentInertiaTensorComponentsOfBodyExertingTorque );

        const double numericalPartialOfContractedInertiaTensorWrtCoordinate =
                ( upPerturbedAuxiliaryQuantities.contractedInertiaTensorOfBodyExertingTorque -
                  downPerturbedAuxiliaryQuantities.contractedInertiaTensorOfBodyExertingTorque ) /
                ( 2.0 * positionPerturbation );
        const double analyticalPartialOfContractedInertiaTensorWrtCoordinate =
                acceleration_partials::detail::computePartialOfContractedInertiaTensorOfBodyExertingTorqueWrtCoordinate(
                        auxiliaryQuantities, coordinateIndex );

        // Verify the analytical derivative of the contracted inertia scalar matches a central difference.
        BOOST_CHECK_SMALL( std::fabs( analyticalPartialOfContractedInertiaTensorWrtCoordinate -
                                      numericalPartialOfContractedInertiaTensorWrtCoordinate ) /
                                   std::max( 1.0, std::fabs( numericalPartialOfContractedInertiaTensorWrtCoordinate ) ),
                           5.0E-9 );

        const Eigen::Matrix< double, 6, 1 > numericalPartialOfAuxiliaryFunctionsWrtCoordinate =
                ( getAuxiliaryFunctionVector( upPerturbedAuxiliaryQuantities ) -
                  getAuxiliaryFunctionVector( downPerturbedAuxiliaryQuantities ) ) /
                ( 2.0 * positionPerturbation );
        const Eigen::Matrix< double, 6, 1 > analyticalPartialOfAuxiliaryFunctionsWrtCoordinate =
                acceleration_partials::detail::computePartialOfAuxiliaryFunctionsWrtPositionCoordinate( auxiliaryQuantities,
                                                                                                        coordinateIndex );

        // Verify the analytical position derivative of the Schutz auxiliary f/g vector.
        checkMatrixClosePerElement( analyticalPartialOfAuxiliaryFunctionsWrtCoordinate,
                                    numericalPartialOfAuxiliaryFunctionsWrtCoordinate,
                                    5.0E-7,
                                    "Schutz auxiliary functions wrt position coordinate " + std::to_string( coordinateIndex ) );
    }

    const double inertiaTensorComponentPerturbation = 1.0E7;
    Eigen::Matrix< double, 6, 6 > numericalPartialOfAuxiliaryFunctionsWrtInertiaTensorComponents = Eigen::Matrix< double, 6, 6 >::Zero( );
    for( int componentIndex = 0; componentIndex < 6; componentIndex++ )
    {
        Eigen::Matrix< double, 6, 1 > upPerturbedIndependentInertiaTensorComponents =
                independentInertiaTensorComponentsOfBodyExertingTorque;
        upPerturbedIndependentInertiaTensorComponents( componentIndex ) += inertiaTensorComponentPerturbation;
        const auto upPerturbedAuxiliaryQuantities = acceleration_partials::detail::computeFourthDegreeTorqueAuxiliaryQuantities(
                relativePositionInBodyFixedFrame, massOfBodyExertingTorque, upPerturbedIndependentInertiaTensorComponents );

        Eigen::Matrix< double, 6, 1 > downPerturbedIndependentInertiaTensorComponents =
                independentInertiaTensorComponentsOfBodyExertingTorque;
        downPerturbedIndependentInertiaTensorComponents( componentIndex ) -= inertiaTensorComponentPerturbation;
        const auto downPerturbedAuxiliaryQuantities = acceleration_partials::detail::computeFourthDegreeTorqueAuxiliaryQuantities(
                relativePositionInBodyFixedFrame, massOfBodyExertingTorque, downPerturbedIndependentInertiaTensorComponents );

        numericalPartialOfAuxiliaryFunctionsWrtInertiaTensorComponents.col( componentIndex ) =
                ( getAuxiliaryFunctionVector( upPerturbedAuxiliaryQuantities ) -
                  getAuxiliaryFunctionVector( downPerturbedAuxiliaryQuantities ) ) /
                ( 2.0 * inertiaTensorComponentPerturbation );
    }

    const Eigen::Matrix< double, 6, 6 > analyticalPartialOfAuxiliaryFunctionsWrtInertiaTensorComponents =
            acceleration_partials::detail::computePartialOfAuxiliaryFunctionsWrtIndependentInertiaTensorComponentsOfBodyExertingTorque(
                    auxiliaryQuantities );

    // Verify the analytical inertia-component derivative of the Schutz auxiliary f/g vector.
    checkMatrixClosePerElement( analyticalPartialOfAuxiliaryFunctionsWrtInertiaTensorComponents,
                                numericalPartialOfAuxiliaryFunctionsWrtInertiaTensorComponents,
                                5.0E-5,
                                "Schutz auxiliary functions wrt body-2 inertia components" );
}

BOOST_AUTO_TEST_CASE( testFourthDegreeFullTwoBodyGravitationalTorquePartials )
{
    // Test rationale:
    // Validate analytical partials of the fourth-degree full-two-body torque model against numerical
    // finite differences for orientation, translational state, and degree-2 coefficient blocks of both bodies.
    //
    // This test covers both arbitrary and aligned attitude configurations, and both non-zero and zero
    // scaled mean moments of inertia, ensuring robustness of the closed-form partial implementation.
    const std::string bodyUndergoingTorqueName = "BodyA";
    const std::string bodyExertingTorqueName = "BodyB";
    const double testTime = 1250.0;
    const double quaternionPerturbationValue = 1.0E-9;
    const double positionPerturbationValue = 5.0E-3;
    const double coefficientPerturbationValue = 1.0E-2;

    const Eigen::Vector4d quaternionPerturbation = Eigen::Vector4d::Constant( quaternionPerturbationValue );

    for( const bool useZeroScaledMeanMomentOfInertia : { false, true } )
    {
        for( const bool useArbitraryRotationStates : { false, true } )
        {
            SystemOfBodies bodies = createTwoBodyTorquePartialTestSystem( testTime,
                                                                          bodyUndergoingTorqueName,
                                                                          bodyExertingTorqueName,
                                                                          useArbitraryRotationStates,
                                                                          useZeroScaledMeanMomentOfInertia );
            std::shared_ptr< Body > bodyUndergoingTorque = bodies.at( bodyUndergoingTorqueName );
            std::shared_ptr< Body > bodyExertingTorque = bodies.at( bodyExertingTorqueName );

            std::shared_ptr< TorqueModel > torqueModel =
                    createFourthDegreeFullTwoBodyGravitationalTorqueModel( bodyUndergoingTorque,
                                                                           bodyExertingTorque,
                                                                           fourthDegreeFullTwoBodyGravitationalTorque( ),
                                                                           bodyUndergoingTorqueName,
                                                                           bodyExertingTorqueName );

            std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterSettings;
            parameterSettings.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                    2, 0, 2, 2, bodyUndergoingTorqueName, spherical_harmonics_cosine_coefficient_block ) );
            parameterSettings.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                    2, 1, 2, 2, bodyUndergoingTorqueName, spherical_harmonics_sine_coefficient_block ) );
            parameterSettings.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                    2, 0, 2, 2, bodyExertingTorqueName, spherical_harmonics_cosine_coefficient_block ) );
            parameterSettings.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                    2, 1, 2, 2, bodyExertingTorqueName, spherical_harmonics_sine_coefficient_block ) );
            std::shared_ptr< EstimatableParameterSet< double > > parameterSet = createParametersToEstimate( parameterSettings, bodies );

            std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > cosineCoefficientsOfBodyUndergoingTorqueParameter =
                    parameterSet->getEstimatedVectorParameters( ).at( 0 );
            std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > sineCoefficientsOfBodyUndergoingTorqueParameter =
                    parameterSet->getEstimatedVectorParameters( ).at( 1 );
            std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > cosineCoefficientsOfBodyExertingTorqueParameter =
                    parameterSet->getEstimatedVectorParameters( ).at( 2 );
            std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > sineCoefficientsOfBodyExertingTorqueParameter =
                    parameterSet->getEstimatedVectorParameters( ).at( 3 );
            std::shared_ptr< GravitationalParameter > gravitationalParameterOfBodyUndergoingTorque =
                    std::make_shared< GravitationalParameter >( bodyUndergoingTorque->getGravityFieldModel( ), bodyUndergoingTorqueName );
            std::shared_ptr< GravitationalParameter > gravitationalParameterOfBodyExertingTorque =
                    std::make_shared< GravitationalParameter >( bodyExertingTorque->getGravityFieldModel( ), bodyExertingTorqueName );

            std::shared_ptr< TorquePartial > torquePartial =
                    createAnalyticalTorquePartial( torqueModel,
                                                   std::make_pair( bodyUndergoingTorqueName, bodyUndergoingTorque ),
                                                   std::make_pair( bodyExertingTorqueName, bodyExertingTorque ),
                                                   basic_astrodynamics::SingleBodyTorqueModelMap( ),
                                                   bodies,
                                                   parameterSet );
            torquePartial->update( testTime );

            Eigen::MatrixXd analyticalPartialWrtOrientationOfBodyUndergoingTorque = Eigen::MatrixXd::Zero( 3, 4 );
            torquePartial->wrtOrientationOfAcceleratedBody( analyticalPartialWrtOrientationOfBodyUndergoingTorque.block( 0, 0, 3, 4 ) );
            Eigen::MatrixXd analyticalPartialWrtOrientationOfBodyExertingTorque = Eigen::MatrixXd::Zero( 3, 4 );
            torquePartial->wrtOrientationOfAcceleratingBody( analyticalPartialWrtOrientationOfBodyExertingTorque.block( 0, 0, 3, 4 ) );

            Eigen::MatrixXd analyticalPartialWrtTranslationalStateOfBodyUndergoingTorque = Eigen::MatrixXd::Zero( 3, 6 );
            torquePartial->wrtNonRotationalStateOfAdditionalBody(
                    analyticalPartialWrtTranslationalStateOfBodyUndergoingTorque.block( 0, 0, 3, 6 ),
                    std::make_pair( bodyUndergoingTorqueName, "" ),
                    propagators::translational_state );
            Eigen::MatrixXd analyticalPartialWrtTranslationalStateOfBodyExertingTorque = Eigen::MatrixXd::Zero( 3, 6 );
            torquePartial->wrtNonRotationalStateOfAdditionalBody(
                    analyticalPartialWrtTranslationalStateOfBodyExertingTorque.block( 0, 0, 3, 6 ),
                    std::make_pair( bodyExertingTorqueName, "" ),
                    propagators::translational_state );

            const Eigen::MatrixXd analyticalPartialWrtCosineCoefficientsOfBodyUndergoingTorque =
                    torquePartial->wrtParameter( cosineCoefficientsOfBodyUndergoingTorqueParameter );
            const Eigen::MatrixXd analyticalPartialWrtSineCoefficientsOfBodyUndergoingTorque =
                    torquePartial->wrtParameter( sineCoefficientsOfBodyUndergoingTorqueParameter );
            const Eigen::MatrixXd analyticalPartialWrtCosineCoefficientsOfBodyExertingTorque =
                    torquePartial->wrtParameter( cosineCoefficientsOfBodyExertingTorqueParameter );
            const Eigen::MatrixXd analyticalPartialWrtSineCoefficientsOfBodyExertingTorque =
                    torquePartial->wrtParameter( sineCoefficientsOfBodyExertingTorqueParameter );

            const Eigen::Vector3d positionPerturbation = Eigen::Vector3d::Constant( positionPerturbationValue );
            const Eigen::VectorXd cosineCoefficientPerturbationBodyUndergoingTorque = Eigen::VectorXd::Constant(
                    cosineCoefficientsOfBodyUndergoingTorqueParameter->getParameterSize( ), coefficientPerturbationValue );
            const Eigen::VectorXd sineCoefficientPerturbationBodyUndergoingTorque = Eigen::VectorXd::Constant(
                    sineCoefficientsOfBodyUndergoingTorqueParameter->getParameterSize( ), coefficientPerturbationValue );
            const Eigen::VectorXd cosineCoefficientPerturbationBodyExertingTorque = Eigen::VectorXd::Constant(
                    cosineCoefficientsOfBodyExertingTorqueParameter->getParameterSize( ), coefficientPerturbationValue );
            const Eigen::VectorXd sineCoefficientPerturbationBodyExertingTorque = Eigen::VectorXd::Constant(
                    sineCoefficientsOfBodyExertingTorqueParameter->getParameterSize( ), coefficientPerturbationValue );

            const std::function< void( Eigen::Vector7d ) > setRotationalStateOfBodyUndergoingTorque =
                    std::bind( &Body::setCurrentRotationalStateToLocalFrame, bodyUndergoingTorque, std::placeholders::_1 );
            const std::function< void( Eigen::Vector7d ) > setRotationalStateOfBodyExertingTorque =
                    std::bind( &Body::setCurrentRotationalStateToLocalFrame, bodyExertingTorque, std::placeholders::_1 );
            const std::function< void( Eigen::Vector6d ) > setTranslationalStateOfBodyUndergoingTorque =
                    std::bind( &Body::setState, bodyUndergoingTorque, std::placeholders::_1 );
            const std::function< void( Eigen::Vector6d ) > setTranslationalStateOfBodyExertingTorque =
                    std::bind( &Body::setState, bodyExertingTorque, std::placeholders::_1 );

            std::vector< Eigen::Vector4d > appliedQuaternionPerturbationOfBodyUndergoingTorque;
            const Eigen::MatrixXd torqueDeviationDueToOrientationChangeOfBodyUndergoingTorque =
                    calculateTorqueDeviationDueToOrientationChange( setRotationalStateOfBodyUndergoingTorque,
                                                                    torqueModel,
                                                                    bodyUndergoingTorque->getRotationalStateVector( ),
                                                                    quaternionPerturbation,
                                                                    appliedQuaternionPerturbationOfBodyUndergoingTorque,
                                                                    emptyFunction,
                                                                    testTime );

            std::vector< Eigen::Vector4d > appliedQuaternionPerturbationOfBodyExertingTorque;
            const Eigen::MatrixXd torqueDeviationDueToOrientationChangeOfBodyExertingTorque =
                    calculateTorqueDeviationDueToOrientationChange( setRotationalStateOfBodyExertingTorque,
                                                                    torqueModel,
                                                                    bodyExertingTorque->getRotationalStateVector( ),
                                                                    quaternionPerturbation,
                                                                    appliedQuaternionPerturbationOfBodyExertingTorque,
                                                                    emptyFunction,
                                                                    testTime );

            const Eigen::MatrixXd numericalPartialWrtPositionOfBodyUndergoingTorque =
                    calculateTorqueWrtTranslationalStatePartials( setTranslationalStateOfBodyUndergoingTorque,
                                                                  torqueModel,
                                                                  bodyUndergoingTorque->getState( ),
                                                                  positionPerturbation,
                                                                  0,
                                                                  emptyFunction,
                                                                  testTime );
            const Eigen::MatrixXd numericalPartialWrtPositionOfBodyExertingTorque =
                    calculateTorqueWrtTranslationalStatePartials( setTranslationalStateOfBodyExertingTorque,
                                                                  torqueModel,
                                                                  bodyExertingTorque->getState( ),
                                                                  positionPerturbation,
                                                                  0,
                                                                  emptyFunction,
                                                                  testTime );
            const std::function< void( ) > updateFunction =
                    std::bind( &updateBodyMassDistributions, bodyUndergoingTorque, bodyExertingTorque, testTime );

            const Eigen::MatrixXd numericalPartialWrtCosineCoefficientsOfBodyUndergoingTorque =
                    calculateTorqueWrtParameterPartials( cosineCoefficientsOfBodyUndergoingTorqueParameter,
                                                         torqueModel,
                                                         cosineCoefficientPerturbationBodyUndergoingTorque,
                                                         updateFunction,
                                                         testTime );
            const Eigen::MatrixXd numericalPartialWrtSineCoefficientsOfBodyUndergoingTorque =
                    calculateTorqueWrtParameterPartials( sineCoefficientsOfBodyUndergoingTorqueParameter,
                                                         torqueModel,
                                                         sineCoefficientPerturbationBodyUndergoingTorque,
                                                         updateFunction,
                                                         testTime );
            const Eigen::MatrixXd numericalPartialWrtCosineCoefficientsOfBodyExertingTorque =
                    calculateTorqueWrtParameterPartials( cosineCoefficientsOfBodyExertingTorqueParameter,
                                                         torqueModel,
                                                         cosineCoefficientPerturbationBodyExertingTorque,
                                                         updateFunction,
                                                         testTime );
            const Eigen::MatrixXd numericalPartialWrtSineCoefficientsOfBodyExertingTorque =
                    calculateTorqueWrtParameterPartials( sineCoefficientsOfBodyExertingTorqueParameter,
                                                         torqueModel,
                                                         sineCoefficientPerturbationBodyExertingTorque,
                                                         updateFunction,
                                                         testTime );

            if( useArbitraryRotationStates && !useZeroScaledMeanMomentOfInertia )
            {
                std::shared_ptr< gravitation::SphericalHarmonicsGravityField > gravityFieldOfBodyUndergoingTorque =
                        std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >(
                                bodyUndergoingTorque->getGravityFieldModel( ) );
                std::shared_ptr< gravitation::SphericalHarmonicsGravityField > gravityFieldOfBodyExertingTorque =
                        std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >(
                                bodyExertingTorque->getGravityFieldModel( ) );
                const Eigen::MatrixXd baseCosineCoefficientsOfBodyUndergoingTorque =
                        gravityFieldOfBodyUndergoingTorque->getCosineCoefficients( );
                const Eigen::MatrixXd baseSineCoefficientsOfBodyUndergoingTorque =
                        gravityFieldOfBodyUndergoingTorque->getSineCoefficients( );
                const Eigen::MatrixXd baseCosineCoefficientsOfBodyExertingTorque =
                        gravityFieldOfBodyExertingTorque->getCosineCoefficients( );
                const Eigen::MatrixXd baseSineCoefficientsOfBodyExertingTorque = gravityFieldOfBodyExertingTorque->getSineCoefficients( );
                const Eigen::MatrixXd zeroVariationBlock = Eigen::MatrixXd::Zero( 3, 5 );
                const std::map< int, std::vector< std::pair< int, int > > > variationCosineIndices = {
                    { 1, { { 2, 0 }, { 3, 1 }, { 4, 2 } } }
                };
                const std::map< int, std::vector< std::pair< int, int > > > variationSineIndices = { { 1,
                                                                                                       { { 2, 1 }, { 3, 2 }, { 4, 4 } } } };

                auto resetGravityFields = [ & ]( ) {
                    gravityFieldOfBodyUndergoingTorque->setCosineCoefficients( baseCosineCoefficientsOfBodyUndergoingTorque );
                    gravityFieldOfBodyUndergoingTorque->setSineCoefficients( baseSineCoefficientsOfBodyUndergoingTorque );
                    gravityFieldOfBodyExertingTorque->setCosineCoefficients( baseCosineCoefficientsOfBodyExertingTorque );
                    gravityFieldOfBodyExertingTorque->setSineCoefficients( baseSineCoefficientsOfBodyExertingTorque );
                };
                auto checkPolynomialVariationPartial =
                        [ & ]( const std::string& bodyName,
                               const std::shared_ptr< gravitation::SphericalHarmonicsGravityField >& gravityField,
                               const Eigen::MatrixXd& baseCosineCoefficients,
                               const Eigen::MatrixXd& baseSineCoefficients,
                               const std::string& label ) {
                            std::shared_ptr< gravitation::PolynomialGravityFieldVariations > variationModel =
                                    std::make_shared< gravitation::PolynomialGravityFieldVariations >(
                                            std::map< int, Eigen::MatrixXd >{ { 1, zeroVariationBlock } },
                                            std::map< int, Eigen::MatrixXd >{ { 1, zeroVariationBlock } },
                                            testTime - 0.3,
                                            2,
                                            0 );
                            std::shared_ptr< PolynomialGravityFieldVariationsParameters > parameter =
                                    std::make_shared< PolynomialGravityFieldVariationsParameters >(
                                            variationModel, variationCosineIndices, variationSineIndices, bodyName );

                            resetGravityFields( );
                            applyGravityFieldVariation(
                                    gravityField, baseCosineCoefficients, baseSineCoefficients, variationModel, testTime );
                            updateFunction( );
                            torqueModel->resetCurrentTime( );
                            torquePartial->resetCurrentTime( );
                            torqueModel->updateMembers( testTime );
                            torquePartial->update( testTime );
                            const Eigen::MatrixXd analyticalPartial = torquePartial->wrtParameter( parameter );
                            const std::function< void( ) > updateVariationEnvironment = [ & ]( ) {
                                resetGravityFields( );
                                applyGravityFieldVariation(
                                        gravityField, baseCosineCoefficients, baseSineCoefficients, variationModel, testTime );
                                updateFunction( );
                            };
                            const Eigen::MatrixXd numericalPartial = calculateTorqueWrtParameterPartials(
                                    parameter,
                                    torqueModel,
                                    Eigen::VectorXd::Constant( parameter->getParameterSize( ), 2.0E-5 ),
                                    updateVariationEnvironment,
                                    testTime );
                            checkMatrixClosePerElement( analyticalPartial, numericalPartial, 5.0E-5, label );
                            resetGravityFields( );
                            updateFunction( );
                        };

                checkPolynomialVariationPartial( bodyUndergoingTorqueName,
                                                 gravityFieldOfBodyUndergoingTorque,
                                                 baseCosineCoefficientsOfBodyUndergoingTorque,
                                                 baseSineCoefficientsOfBodyUndergoingTorque,
                                                 "schutz polynomial variations body undergoing" );
                checkPolynomialVariationPartial( bodyExertingTorqueName,
                                                 gravityFieldOfBodyExertingTorque,
                                                 baseCosineCoefficientsOfBodyExertingTorque,
                                                 baseSineCoefficientsOfBodyExertingTorque,
                                                 "schutz polynomial variations body exerting" );
            }

            for( int index = 1; index < 4; index++ )
            {
                const Eigen::Vector3d analyticalTorqueDeviationOfBodyUndergoingTorque =
                        analyticalPartialWrtOrientationOfBodyUndergoingTorque *
                        appliedQuaternionPerturbationOfBodyUndergoingTorque.at( index );
                const Eigen::Vector3d numericalTorqueDeviationOfBodyUndergoingTorque =
                        torqueDeviationDueToOrientationChangeOfBodyUndergoingTorque.col( index - 1 );
                // Verify the quaternion partial for the body undergoing torque through direct attitude perturbation.
                checkMatrixClosePerElement( analyticalTorqueDeviationOfBodyUndergoingTorque,
                                            numericalTorqueDeviationOfBodyUndergoingTorque,
                                            5.0E-9,
                                            "Schutz orientation body undergoing, perturbation " + std::to_string( index ) );

                const Eigen::Vector3d analyticalTorqueDeviationOfBodyExertingTorque =
                        analyticalPartialWrtOrientationOfBodyExertingTorque * appliedQuaternionPerturbationOfBodyExertingTorque.at( index );
                const Eigen::Vector3d numericalTorqueDeviationOfBodyExertingTorque =
                        torqueDeviationDueToOrientationChangeOfBodyExertingTorque.col( index - 1 );
                // Verify the quaternion partial for the body exerting torque through direct attitude perturbation.
                checkMatrixClosePerElement( analyticalTorqueDeviationOfBodyExertingTorque,
                                            numericalTorqueDeviationOfBodyExertingTorque,
                                            5.0E-9,
                                            "Schutz orientation body exerting, perturbation " + std::to_string( index ) );
            }

            const double analyticalPositionPartialNorm =
                    analyticalPartialWrtTranslationalStateOfBodyUndergoingTorque.block( 0, 0, 3, 3 ).norm( ) +
                    analyticalPartialWrtTranslationalStateOfBodyExertingTorque.block( 0, 0, 3, 3 ).norm( );
            const double numericalPositionPartialNorm =
                    numericalPartialWrtPositionOfBodyUndergoingTorque.norm( ) + numericalPartialWrtPositionOfBodyExertingTorque.norm( );

            // Verify the analytical position partials are active for this non-degenerate torque geometry.
            BOOST_CHECK_GT( analyticalPositionPartialNorm, 1.0E-20 );
            // Verify the numerical position partials are active so the comparison is not vacuous.
            BOOST_CHECK_GT( numericalPositionPartialNorm, 1.0E-20 );

            // Verify the body-undergoing position Jacobian against central finite differences.
            checkMatrixClosePerElement( analyticalPartialWrtTranslationalStateOfBodyUndergoingTorque.block( 0, 0, 3, 3 ),
                                        numericalPartialWrtPositionOfBodyUndergoingTorque,
                                        5.0E-9,
                                        "Schutz position body undergoing" );
            // Verify the body-exerting position Jacobian against central finite differences.
            checkMatrixClosePerElement( analyticalPartialWrtTranslationalStateOfBodyExertingTorque.block( 0, 0, 3, 3 ),
                                        numericalPartialWrtPositionOfBodyExertingTorque,
                                        5.0E-9,
                                        "Schutz position body exerting" );
            // Verify there is no explicit velocity dependence in the fourth-degree torque model.
            BOOST_CHECK_SMALL( analyticalPartialWrtTranslationalStateOfBodyUndergoingTorque.block( 0, 3, 3, 3 ).norm( ), 1.0E-30 );
            // Verify there is no explicit velocity dependence in the fourth-degree torque model.
            BOOST_CHECK_SMALL( analyticalPartialWrtTranslationalStateOfBodyExertingTorque.block( 0, 3, 3, 3 ).norm( ), 1.0E-30 );

            const Eigen::MatrixXd numericalPartialWrtVelocityOfBodyUndergoingTorque =
                    calculateTorqueWrtTranslationalStatePartials( setTranslationalStateOfBodyUndergoingTorque,
                                                                  torqueModel,
                                                                  bodyUndergoingTorque->getState( ),
                                                                  positionPerturbation,
                                                                  3,
                                                                  emptyFunction,
                                                                  testTime );
            const Eigen::MatrixXd numericalPartialWrtVelocityOfBodyExertingTorque =
                    calculateTorqueWrtTranslationalStatePartials( setTranslationalStateOfBodyExertingTorque,
                                                                  torqueModel,
                                                                  bodyExertingTorque->getState( ),
                                                                  positionPerturbation,
                                                                  3,
                                                                  emptyFunction,
                                                                  testTime );
            // Verify finite differencing also sees no velocity dependence for the body undergoing torque.
            BOOST_CHECK_SMALL( numericalPartialWrtVelocityOfBodyUndergoingTorque.norm( ), 1.0E-18 );
            // Verify finite differencing also sees no velocity dependence for the body exerting torque.
            BOOST_CHECK_SMALL( numericalPartialWrtVelocityOfBodyExertingTorque.norm( ), 1.0E-18 );

            const double analyticalCoefficientPartialNorm = analyticalPartialWrtCosineCoefficientsOfBodyUndergoingTorque.norm( ) +
                    analyticalPartialWrtSineCoefficientsOfBodyUndergoingTorque.norm( ) +
                    analyticalPartialWrtCosineCoefficientsOfBodyExertingTorque.norm( ) +
                    analyticalPartialWrtSineCoefficientsOfBodyExertingTorque.norm( );
            const double numericalCoefficientPartialNorm = numericalPartialWrtCosineCoefficientsOfBodyUndergoingTorque.norm( ) +
                    numericalPartialWrtSineCoefficientsOfBodyUndergoingTorque.norm( ) +
                    numericalPartialWrtCosineCoefficientsOfBodyExertingTorque.norm( ) +
                    numericalPartialWrtSineCoefficientsOfBodyExertingTorque.norm( );

            // Verify coefficient partials are active for the non-spherical degree-2 bodies.
            BOOST_CHECK_GT( analyticalCoefficientPartialNorm, 1.0E-20 );
            // Verify numerical coefficient partials are active so coefficient comparisons are meaningful.
            BOOST_CHECK_GT( numericalCoefficientPartialNorm, 1.0E-20 );

            // Verify cosine-coefficient partials of the body undergoing torque against finite differences.
            checkMatrixClosePerElement( analyticalPartialWrtCosineCoefficientsOfBodyUndergoingTorque,
                                        numericalPartialWrtCosineCoefficientsOfBodyUndergoingTorque,
                                        1.0E-10,
                                        "Schutz cosine coefficients body undergoing" );
            // Verify sine-coefficient partials of the body undergoing torque against finite differences.
            checkMatrixClosePerElement( analyticalPartialWrtSineCoefficientsOfBodyUndergoingTorque,
                                        numericalPartialWrtSineCoefficientsOfBodyUndergoingTorque,
                                        1.0E-10,
                                        "Schutz sine coefficients body undergoing" );
            // Verify cosine-coefficient partials of the body exerting torque against finite differences.
            checkMatrixClosePerElement( analyticalPartialWrtCosineCoefficientsOfBodyExertingTorque,
                                        numericalPartialWrtCosineCoefficientsOfBodyExertingTorque,
                                        1.0E-10,
                                        "Schutz cosine coefficients body exerting" );
            // Verify sine-coefficient partials of the body exerting torque against finite differences.
            checkMatrixClosePerElement( analyticalPartialWrtSineCoefficientsOfBodyExertingTorque,
                                        numericalPartialWrtSineCoefficientsOfBodyExertingTorque,
                                        1.0E-10,
                                        "Schutz sine coefficients body exerting" );

            makeBodyMassIndependentOfGravityField( bodyUndergoingTorque );
            makeBodyMassIndependentOfGravityField( bodyExertingTorque );
            torqueModel->resetCurrentTime( );
            torquePartial->resetCurrentTime( );
            torqueModel->updateMembers( testTime );
            torquePartial->update( testTime );

            std::shared_ptr< InitialMassStateParameter< double > > massOfBodyUndergoingTorqueParameter =
                    makeInitialMassParameter( bodyUndergoingTorque, bodyUndergoingTorqueName );
            std::shared_ptr< InitialMassStateParameter< double > > massOfBodyExertingTorqueParameter =
                    makeInitialMassParameter( bodyExertingTorque, bodyExertingTorqueName );
            const Eigen::MatrixXd analyticalPartialWrtMassOfBodyUndergoingTorque =
                    torquePartial->wrtParameter( massOfBodyUndergoingTorqueParameter );
            const Eigen::MatrixXd analyticalPartialWrtMassOfBodyExertingTorque =
                    torquePartial->wrtParameter( massOfBodyExertingTorqueParameter );
            const Eigen::MatrixXd numericalPartialWrtMassOfBodyUndergoingTorque = calculateTorqueWrtParameterPartials(
                    massOfBodyUndergoingTorqueParameter,
                    torqueModel,
                    Eigen::VectorXd::Constant( massOfBodyUndergoingTorqueParameter->getParameterSize( ), 1.0E6 ),
                    emptyFunction,
                    testTime );
            const Eigen::MatrixXd numericalPartialWrtMassOfBodyExertingTorque = calculateTorqueWrtParameterPartials(
                    massOfBodyExertingTorqueParameter,
                    torqueModel,
                    Eigen::VectorXd::Constant( massOfBodyExertingTorqueParameter->getParameterSize( ), 1.0E6 ),
                    emptyFunction,
                    testTime );
            const Eigen::Vector3d analyticalPartialWrtGravitationalParameterOfBodyUndergoingTorque =
                    torquePartial->wrtParameter( gravitationalParameterOfBodyUndergoingTorque );
            const Eigen::Vector3d analyticalPartialWrtGravitationalParameterOfBodyExertingTorque =
                    torquePartial->wrtParameter( gravitationalParameterOfBodyExertingTorque );
            const Eigen::Vector3d numericalPartialWrtGravitationalParameterOfBodyUndergoingTorque = calculateTorqueWrtParameterPartials(
                    gravitationalParameterOfBodyUndergoingTorque, torqueModel, 1.0E-2, emptyFunction, testTime );
            const Eigen::Vector3d numericalPartialWrtGravitationalParameterOfBodyExertingTorque = calculateTorqueWrtParameterPartials(
                    gravitationalParameterOfBodyExertingTorque, torqueModel, 1.0E-2, emptyFunction, testTime );

            // Verify the Schutz torque is independent of body-1 mass when inertia is held fixed.
            BOOST_CHECK_SMALL( analyticalPartialWrtMassOfBodyUndergoingTorque.norm( ), 1.0E-20 );
            // Verify the numerical mass perturbation of body 1 also gives zero torque response.
            BOOST_CHECK_SMALL( numericalPartialWrtMassOfBodyUndergoingTorque.norm( ), 1.0E-20 );
            // Verify the analytical body-2 mass partial against a central finite-difference reference.
            checkMatrixClosePerElement( analyticalPartialWrtMassOfBodyExertingTorque,
                                        numericalPartialWrtMassOfBodyExertingTorque,
                                        1.0E-10,
                                        "Schutz mass body exerting" );
            // Verify the torque is independent of body-1 gravitational parameter when inertia is held fixed.
            BOOST_CHECK_SMALL( analyticalPartialWrtGravitationalParameterOfBodyUndergoingTorque.norm( ), 1.0E-20 );
            // Verify finite differencing also gives zero response to body-1 gravitational parameter.
            BOOST_CHECK_SMALL( numericalPartialWrtGravitationalParameterOfBodyUndergoingTorque.norm( ), 1.0E-20 );
            // Verify the Schutz formulation has no direct body-2 gravitational-parameter dependence when mass is fixed separately.
            BOOST_CHECK_SMALL( analyticalPartialWrtGravitationalParameterOfBodyExertingTorque.norm( ), 1.0E-20 );
            // Verify finite differencing also gives no direct body-2 gravitational-parameter response.
            BOOST_CHECK_SMALL( numericalPartialWrtGravitationalParameterOfBodyExertingTorque.norm( ), 1.0E-20 );
        }
    }
}

BOOST_AUTO_TEST_CASE( testFullTwoBodyTorqueRotationModelParameterPartials )
{
    // Test rationale:
    // Verify that rotation-model parameter partials are correctly chained through quaternion/orientation partials
    // for both the Schutz fourth-degree model and the DMR full spherical-harmonic torque model.
    const std::string bodyUndergoingTorqueName = "BodyA";
    const std::string bodyExertingTorqueName = "BodyB";
    const double testTime = 1250.0;

    SystemOfBodies bodies =
            createTwoBodyTorquePartialTestSystem( testTime, bodyUndergoingTorqueName, bodyExertingTorqueName, false, false );
    std::shared_ptr< Body > bodyUndergoingTorque = bodies.at( bodyUndergoingTorqueName );
    std::shared_ptr< Body > bodyExertingTorque = bodies.at( bodyExertingTorqueName );

    std::shared_ptr< ephemerides::SimpleRotationalEphemeris > rotationModelOfBodyUndergoingTorque =
            std::make_shared< ephemerides::SimpleRotationalEphemeris >(
                    0.42, 1.02, -0.25, 2.0E-5, 1000.0, "ECLIPJ2000", bodyUndergoingTorqueName + "_Fixed" );
    std::shared_ptr< ephemerides::SimpleRotationalEphemeris > rotationModelOfBodyExertingTorque =
            std::make_shared< ephemerides::SimpleRotationalEphemeris >(
                    -0.31, 0.91, 0.37, -1.5E-5, 1000.0, "ECLIPJ2000", bodyExertingTorqueName + "_Fixed" );
    bodyUndergoingTorque->setRotationalEphemeris( rotationModelOfBodyUndergoingTorque );
    bodyExertingTorque->setRotationalEphemeris( rotationModelOfBodyExertingTorque );
    bodyUndergoingTorque->setCurrentRotationToLocalFrameFromEphemeris( testTime );
    bodyExertingTorque->setCurrentRotationToLocalFrameFromEphemeris( testTime );

    const observation_partials::RotationMatrixPartialNamedList rotationPartialsOfBodyUndergoingTorque =
            createSimpleRotationPartialMap( bodyUndergoingTorque );
    const observation_partials::RotationMatrixPartialNamedList rotationPartialsOfBodyExertingTorque =
            createSimpleRotationPartialMap( bodyExertingTorque );

    std::shared_ptr< TorqueModel > fourthDegreeTorqueModel =
            createFourthDegreeFullTwoBodyGravitationalTorqueModel( bodyUndergoingTorque,
                                                                   bodyExertingTorque,
                                                                   fourthDegreeFullTwoBodyGravitationalTorque( ),
                                                                   bodyUndergoingTorqueName,
                                                                   bodyExertingTorqueName );
    std::shared_ptr< FourthDegreeFullTwoBodyGravitationalTorquePartial > fourthDegreeTorquePartial =
            std::make_shared< FourthDegreeFullTwoBodyGravitationalTorquePartial >(
                    std::dynamic_pointer_cast< gravitation::FourthDegreeFullTwoBodyGravitationalTorqueModel >( fourthDegreeTorqueModel ),
                    std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >(
                            bodyUndergoingTorque->getGravityFieldModel( ) ),
                    std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >( bodyExertingTorque->getGravityFieldModel( ) ),
                    bodyUndergoingTorqueName,
                    bodyExertingTorqueName,
                    rotationPartialsOfBodyUndergoingTorque,
                    rotationPartialsOfBodyExertingTorque );

    std::shared_ptr< TorqueModel > dmrTorqueModel =
            createFullTwoBodySphericalHarmonicGravitationalTorqueModel( bodyUndergoingTorque,
                                                                        bodyExertingTorque,
                                                                        fullTwoBodySphericalHarmonicGravitationalTorque( 4, 4, 4, 4 ),
                                                                        bodyUndergoingTorqueName,
                                                                        bodyExertingTorqueName );
    std::shared_ptr< gravitation::FullTwoBodySphericalHarmonicTorque > dmrTorqueModelTyped =
            std::dynamic_pointer_cast< gravitation::FullTwoBodySphericalHarmonicTorque >( dmrTorqueModel );
    std::shared_ptr< FullTwoBodySphericalHarmonicsGravityPartial > dmrAccelerationPartial =
            std::make_shared< FullTwoBodySphericalHarmonicsGravityPartial >( bodyUndergoingTorqueName,
                                                                             bodyExertingTorqueName,
                                                                             dmrTorqueModelTyped->getAccelerationBetweenBodies( ),
                                                                             rotationPartialsOfBodyUndergoingTorque,
                                                                             rotationPartialsOfBodyExertingTorque );
    std::shared_ptr< FullTwoBodySphericalHarmonicGravitationalTorquePartial > dmrTorquePartial =
            std::make_shared< FullTwoBodySphericalHarmonicGravitationalTorquePartial >( dmrTorqueModelTyped,
                                                                                        dmrAccelerationPartial,
                                                                                        bodyUndergoingTorqueName,
                                                                                        bodyExertingTorqueName,
                                                                                        rotationPartialsOfBodyUndergoingTorque,
                                                                                        rotationPartialsOfBodyExertingTorque );

    auto updateRotationsAndTorque = [ & ]( const std::shared_ptr< TorqueModel >& torqueModel,
                                           const std::shared_ptr< TorquePartial >& torquePartial ) {
        bodyUndergoingTorque->setCurrentRotationToLocalFrameFromEphemeris( testTime );
        bodyExertingTorque->setCurrentRotationToLocalFrameFromEphemeris( testTime );
        updateBodyMassDistributions( bodyUndergoingTorque, bodyExertingTorque, testTime );
        torqueModel->resetCurrentTime( );
        torquePartial->resetCurrentTime( );
        torqueModel->updateMembers( testTime );
        torquePartial->update( testTime );
        return torqueModel->getTorque( );
    };

    auto checkRotationRatePartial = [ & ]( const std::shared_ptr< TorqueModel >& torqueModel,
                                           const std::shared_ptr< TorquePartial >& torquePartial,
                                           const std::shared_ptr< RotationRate >& rotationRateParameter,
                                           const std::string& label ) {
        updateRotationsAndTorque( torqueModel, torquePartial );
        const Eigen::Vector3d analyticalPartial = torquePartial->wrtParameter( rotationRateParameter );
        const Eigen::Vector3d numericalPartial = calculateTorqueWrtParameterPartials(
                rotationRateParameter, torqueModel, 1.0E-6, emptyFunction, testTime, [ & ]( const double currentTime ) {
                    bodyUndergoingTorque->setCurrentRotationToLocalFrameFromEphemeris( currentTime );
                    bodyExertingTorque->setCurrentRotationToLocalFrameFromEphemeris( currentTime );
                    updateBodyMassDistributions( bodyUndergoingTorque, bodyExertingTorque, currentTime );
                } );
        updateRotationsAndTorque( torqueModel, torquePartial );
        checkMatrixClosePerElement( analyticalPartial, numericalPartial, 1.0E-5, label );
    };

    auto checkPolePositionPartial = [ & ]( const std::shared_ptr< TorqueModel >& torqueModel,
                                           const std::shared_ptr< TorquePartial >& torquePartial,
                                           const std::shared_ptr< ConstantRotationalOrientation >& polePositionParameter,
                                           const std::string& label ) {
        updateRotationsAndTorque( torqueModel, torquePartial );
        const Eigen::MatrixXd analyticalPartial = torquePartial->wrtParameter( polePositionParameter );
        const Eigen::MatrixXd numericalPartial = calculateTorqueWrtParameterPartials(
                polePositionParameter,
                torqueModel,
                Eigen::VectorXd::Constant( polePositionParameter->getParameterSize( ), 3.0E-5 ),
                emptyFunction,
                testTime,
                [ & ]( const double currentTime ) {
                    bodyUndergoingTorque->setCurrentRotationToLocalFrameFromEphemeris( currentTime );
                    bodyExertingTorque->setCurrentRotationToLocalFrameFromEphemeris( currentTime );
                    updateBodyMassDistributions( bodyUndergoingTorque, bodyExertingTorque, currentTime );
                } );
        updateRotationsAndTorque( torqueModel, torquePartial );
        checkMatrixClosePerElement( analyticalPartial, numericalPartial, 1.0E-5, label );
    };

    std::shared_ptr< RotationRate > rotationRateBodyUndergoingTorque =
            std::make_shared< RotationRate >( rotationModelOfBodyUndergoingTorque, bodyUndergoingTorqueName );
    std::shared_ptr< RotationRate > rotationRateBodyExertingTorque =
            std::make_shared< RotationRate >( rotationModelOfBodyExertingTorque, bodyExertingTorqueName );
    std::shared_ptr< ConstantRotationalOrientation > polePositionBodyUndergoingTorque =
            std::make_shared< ConstantRotationalOrientation >( rotationModelOfBodyUndergoingTorque, bodyUndergoingTorqueName );
    std::shared_ptr< ConstantRotationalOrientation > polePositionBodyExertingTorque =
            std::make_shared< ConstantRotationalOrientation >( rotationModelOfBodyExertingTorque, bodyExertingTorqueName );

    // Verify Schutz torque sensitivity to body-1 constant rotation rate.
    checkRotationRatePartial( fourthDegreeTorqueModel, fourthDegreeTorquePartial, rotationRateBodyUndergoingTorque, "Schutz rate body 1" );
    // Verify Schutz torque sensitivity to body-2 constant rotation rate.
    checkRotationRatePartial( fourthDegreeTorqueModel, fourthDegreeTorquePartial, rotationRateBodyExertingTorque, "Schutz rate body 2" );
    // Verify Schutz torque sensitivity to body-1 pole orientation.
    checkPolePositionPartial( fourthDegreeTorqueModel, fourthDegreeTorquePartial, polePositionBodyUndergoingTorque, "Schutz pole body 1" );
    // Verify Schutz torque sensitivity to body-2 pole orientation.
    checkPolePositionPartial( fourthDegreeTorqueModel, fourthDegreeTorquePartial, polePositionBodyExertingTorque, "Schutz pole body 2" );

    // Verify DMR torque sensitivity to body-1 constant rotation rate.
    checkRotationRatePartial( dmrTorqueModel, dmrTorquePartial, rotationRateBodyUndergoingTorque, "DMR rate body 1" );
    // Verify DMR torque sensitivity to body-2 constant rotation rate.
    checkRotationRatePartial( dmrTorqueModel, dmrTorquePartial, rotationRateBodyExertingTorque, "DMR rate body 2" );
    // Verify DMR torque sensitivity to body-1 pole orientation.
    checkPolePositionPartial( dmrTorqueModel, dmrTorquePartial, polePositionBodyUndergoingTorque, "DMR pole body 1" );
    // Verify DMR torque sensitivity to body-2 pole orientation.
    checkPolePositionPartial( dmrTorqueModel, dmrTorquePartial, polePositionBodyExertingTorque, "DMR pole body 2" );
}

BOOST_AUTO_TEST_CASE( testFullTwoBodySphericalHarmonicGravitationalTorquePartials )
{
    // Test rationale:
    // Validate analytical partial derivatives of the full two-body spherical-harmonic torque model
    // against independent finite-difference references for:
    // 1) both body orientations (quaternion perturbations),
    // 2) both body translational states (position blocks),
    // 3) degree-2 through degree-4 spherical-harmonic coefficients of both bodies (cosine/sine).
    //
    // The torque model is based on the mutual-potential coupling expansion (Dirkx et al., 2019). This test
    // verifies that the implemented Jacobians are consistent with the implemented torque itself for both aligned
    // and arbitrary body-frame orientations.
    const std::string bodyUndergoingTorqueName = "BodyA";
    const std::string bodyExertingTorqueName = "BodyB";
    const double testTime = 1250.0;

    for( const bool useArbitraryRotationStates : { false, true } )
    {
        SystemOfBodies bodies = createTwoBodyTorquePartialTestSystem(
                testTime, bodyUndergoingTorqueName, bodyExertingTorqueName, useArbitraryRotationStates );
        std::shared_ptr< Body > bodyUndergoingTorque = bodies.at( bodyUndergoingTorqueName );
        std::shared_ptr< Body > bodyExertingTorque = bodies.at( bodyExertingTorqueName );

        std::shared_ptr< TorqueModel > torqueModel =
                createFullTwoBodySphericalHarmonicGravitationalTorqueModel( bodyUndergoingTorque,
                                                                            bodyExertingTorque,
                                                                            fullTwoBodySphericalHarmonicGravitationalTorque( 4, 4, 4, 4 ),
                                                                            bodyUndergoingTorqueName,
                                                                            bodyExertingTorqueName );

        std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterSettings;
        parameterSettings.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                2, 0, 4, 4, bodyUndergoingTorqueName, spherical_harmonics_cosine_coefficient_block ) );
        parameterSettings.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                2, 1, 4, 4, bodyUndergoingTorqueName, spherical_harmonics_sine_coefficient_block ) );
        parameterSettings.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                2, 0, 4, 4, bodyExertingTorqueName, spherical_harmonics_cosine_coefficient_block ) );
        parameterSettings.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                2, 1, 4, 4, bodyExertingTorqueName, spherical_harmonics_sine_coefficient_block ) );
        std::shared_ptr< EstimatableParameterSet< double > > parameterSet = createParametersToEstimate( parameterSettings, bodies );

        std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > cosineCoefficientsOfBodyUndergoingTorqueParameter =
                parameterSet->getEstimatedVectorParameters( ).at( 0 );
        std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > sineCoefficientsOfBodyUndergoingTorqueParameter =
                parameterSet->getEstimatedVectorParameters( ).at( 1 );
        std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > cosineCoefficientsOfBodyExertingTorqueParameter =
                parameterSet->getEstimatedVectorParameters( ).at( 2 );
        std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > sineCoefficientsOfBodyExertingTorqueParameter =
                parameterSet->getEstimatedVectorParameters( ).at( 3 );
        std::shared_ptr< GravitationalParameter > gravitationalParameterOfBodyUndergoingTorque =
                std::make_shared< GravitationalParameter >( bodyUndergoingTorque->getGravityFieldModel( ), bodyUndergoingTorqueName );
        std::shared_ptr< GravitationalParameter > gravitationalParameterOfBodyExertingTorque =
                std::make_shared< GravitationalParameter >( bodyExertingTorque->getGravityFieldModel( ), bodyExertingTorqueName );

        std::shared_ptr< TorquePartial > torquePartial =
                createAnalyticalTorquePartial( torqueModel,
                                               std::make_pair( bodyUndergoingTorqueName, bodyUndergoingTorque ),
                                               std::make_pair( bodyExertingTorqueName, bodyExertingTorque ),
                                               basic_astrodynamics::SingleBodyTorqueModelMap( ),
                                               bodies,
                                               parameterSet );
        torquePartial->update( testTime );

        Eigen::MatrixXd analyticalPartialWrtOrientationOfBodyUndergoingTorque = Eigen::MatrixXd::Zero( 3, 4 );
        torquePartial->wrtOrientationOfAcceleratedBody( analyticalPartialWrtOrientationOfBodyUndergoingTorque.block( 0, 0, 3, 4 ) );
        Eigen::MatrixXd analyticalPartialWrtOrientationOfBodyExertingTorque = Eigen::MatrixXd::Zero( 3, 4 );
        torquePartial->wrtOrientationOfAcceleratingBody( analyticalPartialWrtOrientationOfBodyExertingTorque.block( 0, 0, 3, 4 ) );

        Eigen::MatrixXd analyticalPartialWrtTranslationalStateOfBodyUndergoingTorque = Eigen::MatrixXd::Zero( 3, 6 );
        torquePartial->wrtNonRotationalStateOfAdditionalBody(
                analyticalPartialWrtTranslationalStateOfBodyUndergoingTorque.block( 0, 0, 3, 6 ),
                std::make_pair( bodyUndergoingTorqueName, "" ),
                propagators::translational_state );
        Eigen::MatrixXd analyticalPartialWrtTranslationalStateOfBodyExertingTorque = Eigen::MatrixXd::Zero( 3, 6 );
        torquePartial->wrtNonRotationalStateOfAdditionalBody(
                analyticalPartialWrtTranslationalStateOfBodyExertingTorque.block( 0, 0, 3, 6 ),
                std::make_pair( bodyExertingTorqueName, "" ),
                propagators::translational_state );

        const Eigen::MatrixXd analyticalPartialWrtCosineCoefficientsOfBodyUndergoingTorque =
                torquePartial->wrtParameter( cosineCoefficientsOfBodyUndergoingTorqueParameter );
        const Eigen::MatrixXd analyticalPartialWrtSineCoefficientsOfBodyUndergoingTorque =
                torquePartial->wrtParameter( sineCoefficientsOfBodyUndergoingTorqueParameter );
        const Eigen::MatrixXd analyticalPartialWrtCosineCoefficientsOfBodyExertingTorque =
                torquePartial->wrtParameter( cosineCoefficientsOfBodyExertingTorqueParameter );
        const Eigen::MatrixXd analyticalPartialWrtSineCoefficientsOfBodyExertingTorque =
                torquePartial->wrtParameter( sineCoefficientsOfBodyExertingTorqueParameter );

        // These finite-difference steps were selected from a perturbation sweep: smaller translational/coefficient
        // steps are roundoff dominated for dimensional torques, while larger translational steps show truncation error.
        const Eigen::Vector4d orientationPerturbation = Eigen::Vector4d::Constant( 1.0E-10 );
        const Eigen::Vector3d positionPerturbation = Eigen::Vector3d::Constant( 1.0E-1 );
        const Eigen::VectorXd cosineCoefficientPerturbationBodyUndergoingTorque =
                Eigen::VectorXd::Constant( cosineCoefficientsOfBodyUndergoingTorqueParameter->getParameterSize( ), 100.0 );
        const Eigen::VectorXd sineCoefficientPerturbationBodyUndergoingTorque =
                Eigen::VectorXd::Constant( sineCoefficientsOfBodyUndergoingTorqueParameter->getParameterSize( ), 100.0 );
        const Eigen::VectorXd cosineCoefficientPerturbationBodyExertingTorque =
                Eigen::VectorXd::Constant( cosineCoefficientsOfBodyExertingTorqueParameter->getParameterSize( ), 100.0 );
        const Eigen::VectorXd sineCoefficientPerturbationBodyExertingTorque =
                Eigen::VectorXd::Constant( sineCoefficientsOfBodyExertingTorqueParameter->getParameterSize( ), 100.0 );

        const std::function< void( Eigen::Vector7d ) > setRotationalStateOfBodyUndergoingTorque =
                std::bind( &Body::setCurrentRotationalStateToLocalFrame, bodyUndergoingTorque, std::placeholders::_1 );
        const std::function< void( Eigen::Vector7d ) > setRotationalStateOfBodyExertingTorque =
                std::bind( &Body::setCurrentRotationalStateToLocalFrame, bodyExertingTorque, std::placeholders::_1 );
        const std::function< void( Eigen::Vector6d ) > setTranslationalStateOfBodyUndergoingTorque =
                std::bind( &Body::setState, bodyUndergoingTorque, std::placeholders::_1 );
        const std::function< void( Eigen::Vector6d ) > setTranslationalStateOfBodyExertingTorque =
                std::bind( &Body::setState, bodyExertingTorque, std::placeholders::_1 );

        // Numerical references are computed with central finite differences and compared against analytical
        // partials from FullTwoBodySphericalHarmonicGravitationalTorquePartial. Agreement here confirms that
        // coefficient/attitude/state sensitivities are implemented correctly and consistently with the model output.
        std::vector< Eigen::Vector4d > appliedQuaternionPerturbationOfBodyUndergoingTorque;
        const Eigen::MatrixXd torqueDeviationDueToOrientationChangeOfBodyUndergoingTorque =
                calculateTorqueDeviationDueToOrientationChange( setRotationalStateOfBodyUndergoingTorque,
                                                                torqueModel,
                                                                bodyUndergoingTorque->getRotationalStateVector( ),
                                                                orientationPerturbation,
                                                                appliedQuaternionPerturbationOfBodyUndergoingTorque,
                                                                emptyFunction,
                                                                testTime );
        std::vector< Eigen::Vector4d > appliedQuaternionPerturbationOfBodyExertingTorque;
        const Eigen::MatrixXd torqueDeviationDueToOrientationChangeOfBodyExertingTorque =
                calculateTorqueDeviationDueToOrientationChange( setRotationalStateOfBodyExertingTorque,
                                                                torqueModel,
                                                                bodyExertingTorque->getRotationalStateVector( ),
                                                                orientationPerturbation,
                                                                appliedQuaternionPerturbationOfBodyExertingTorque,
                                                                emptyFunction,
                                                                testTime );
        const Eigen::MatrixXd numericalPartialWrtPositionOfBodyUndergoingTorque =
                calculateTorqueWrtTranslationalStatePartials( setTranslationalStateOfBodyUndergoingTorque,
                                                              torqueModel,
                                                              bodyUndergoingTorque->getState( ),
                                                              positionPerturbation,
                                                              0,
                                                              emptyFunction,
                                                              testTime );
        const Eigen::MatrixXd numericalPartialWrtPositionOfBodyExertingTorque =
                calculateTorqueWrtTranslationalStatePartials( setTranslationalStateOfBodyExertingTorque,
                                                              torqueModel,
                                                              bodyExertingTorque->getState( ),
                                                              positionPerturbation,
                                                              0,
                                                              emptyFunction,
                                                              testTime );

        const Eigen::MatrixXd numericalPartialWrtCosineCoefficientsOfBodyUndergoingTorque =
                calculateTorqueWrtParameterPartials( cosineCoefficientsOfBodyUndergoingTorqueParameter,
                                                     torqueModel,
                                                     cosineCoefficientPerturbationBodyUndergoingTorque,
                                                     emptyFunction,
                                                     testTime );
        const Eigen::MatrixXd numericalPartialWrtSineCoefficientsOfBodyUndergoingTorque =
                calculateTorqueWrtParameterPartials( sineCoefficientsOfBodyUndergoingTorqueParameter,
                                                     torqueModel,
                                                     sineCoefficientPerturbationBodyUndergoingTorque,
                                                     emptyFunction,
                                                     testTime );
        const Eigen::MatrixXd numericalPartialWrtCosineCoefficientsOfBodyExertingTorque =
                calculateTorqueWrtParameterPartials( cosineCoefficientsOfBodyExertingTorqueParameter,
                                                     torqueModel,
                                                     cosineCoefficientPerturbationBodyExertingTorque,
                                                     emptyFunction,
                                                     testTime );
        const Eigen::MatrixXd numericalPartialWrtSineCoefficientsOfBodyExertingTorque =
                calculateTorqueWrtParameterPartials( sineCoefficientsOfBodyExertingTorqueParameter,
                                                     torqueModel,
                                                     sineCoefficientPerturbationBodyExertingTorque,
                                                     emptyFunction,
                                                     testTime );

        if( useArbitraryRotationStates )
        {
            std::shared_ptr< gravitation::SphericalHarmonicsGravityField > gravityFieldOfBodyUndergoingTorque =
                    std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >(
                            bodyUndergoingTorque->getGravityFieldModel( ) );
            std::shared_ptr< gravitation::SphericalHarmonicsGravityField > gravityFieldOfBodyExertingTorque =
                    std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >( bodyExertingTorque->getGravityFieldModel( ) );
            const Eigen::MatrixXd baseCosineCoefficientsOfBodyUndergoingTorque =
                    gravityFieldOfBodyUndergoingTorque->getCosineCoefficients( );
            const Eigen::MatrixXd baseSineCoefficientsOfBodyUndergoingTorque = gravityFieldOfBodyUndergoingTorque->getSineCoefficients( );
            const Eigen::MatrixXd baseCosineCoefficientsOfBodyExertingTorque = gravityFieldOfBodyExertingTorque->getCosineCoefficients( );
            const Eigen::MatrixXd baseSineCoefficientsOfBodyExertingTorque = gravityFieldOfBodyExertingTorque->getSineCoefficients( );
            const Eigen::MatrixXd zeroVariationBlock = Eigen::MatrixXd::Zero( 3, 5 );
            const std::map< int, std::vector< std::pair< int, int > > > polynomialCosineIndices = { { 1, { { 3, 1 }, { 4, 2 } } },
                                                                                                    { 2, { { 4, 4 } } } };
            const std::map< int, std::vector< std::pair< int, int > > > polynomialSineIndices = { { 1, { { 3, 1 }, { 4, 4 } } } };
            const std::map< int, std::vector< std::pair< int, int > > > periodicCosineIndices = { { 0, { { 3, 0 }, { 4, 4 } } },
                                                                                                  { 1, { { 4, 2 } } } };
            const std::map< int, std::vector< std::pair< int, int > > > periodicSineIndices = { { 0, { { 3, 1 } } }, { 1, { { 4, 4 } } } };
            const std::vector< double > variationFrequencies = { 1.0E-3, 1.7E-3 };
            auto resetGravityFields = [ & ]( ) {
                gravityFieldOfBodyUndergoingTorque->setCosineCoefficients( baseCosineCoefficientsOfBodyUndergoingTorque );
                gravityFieldOfBodyUndergoingTorque->setSineCoefficients( baseSineCoefficientsOfBodyUndergoingTorque );
                gravityFieldOfBodyExertingTorque->setCosineCoefficients( baseCosineCoefficientsOfBodyExertingTorque );
                gravityFieldOfBodyExertingTorque->setSineCoefficients( baseSineCoefficientsOfBodyExertingTorque );
            };
            auto checkGravityFieldVariationPartial =
                    [ & ]( const std::shared_ptr< EstimatableParameter< Eigen::VectorXd > >& parameter,
                           const std::shared_ptr< gravitation::GravityFieldVariations >& variationModel,
                           const std::shared_ptr< gravitation::SphericalHarmonicsGravityField >& gravityField,
                           const Eigen::MatrixXd& baseCosineCoefficients,
                           const Eigen::MatrixXd& baseSineCoefficients,
                           const std::string& label ) {
                        resetGravityFields( );
                        applyGravityFieldVariation( gravityField, baseCosineCoefficients, baseSineCoefficients, variationModel, testTime );
                        torqueModel->resetCurrentTime( );
                        torquePartial->resetCurrentTime( );
                        torqueModel->updateMembers( testTime );
                        torquePartial->update( testTime );
                        const Eigen::MatrixXd analyticalPartial = torquePartial->wrtParameter( parameter );
                        const std::function< void( ) > updateVariationEnvironment = [ & ]( ) {
                            resetGravityFields( );
                            applyGravityFieldVariation(
                                    gravityField, baseCosineCoefficients, baseSineCoefficients, variationModel, testTime );
                        };
                        const Eigen::MatrixXd numericalPartial =
                                calculateTorqueWrtParameterPartials( parameter,
                                                                     torqueModel,
                                                                     Eigen::VectorXd::Constant( parameter->getParameterSize( ), 5.0E-2 ),
                                                                     updateVariationEnvironment,
                                                                     testTime );
                        checkMatrixClosePerElement( analyticalPartial, numericalPartial, 1.0E-5, label );
                        resetGravityFields( );
                    };
            auto createPolynomialVariationParameter =
                    [ & ]( const std::string& bodyName ) -> std::shared_ptr< PolynomialGravityFieldVariationsParameters > {
                std::shared_ptr< gravitation::PolynomialGravityFieldVariations > variationModel =
                        std::make_shared< gravitation::PolynomialGravityFieldVariations >(
                                std::map< int, Eigen::MatrixXd >{ { 1, zeroVariationBlock }, { 2, zeroVariationBlock } },
                                std::map< int, Eigen::MatrixXd >{ { 1, zeroVariationBlock } },
                                testTime - 0.25,
                                2,
                                0 );
                return std::make_shared< PolynomialGravityFieldVariationsParameters >(
                        variationModel, polynomialCosineIndices, polynomialSineIndices, bodyName );
            };
            auto createPeriodicVariationParameter =
                    [ & ]( const std::string& bodyName ) -> std::shared_ptr< PeriodicGravityFieldVariationsParameters > {
                const std::vector< Eigen::MatrixXd > zeroPeriodicBlocks( 2, zeroVariationBlock );
                std::shared_ptr< gravitation::PeriodicGravityFieldVariations > variationModel =
                        std::make_shared< gravitation::PeriodicGravityFieldVariations >( zeroPeriodicBlocks,
                                                                                         zeroPeriodicBlocks,
                                                                                         zeroPeriodicBlocks,
                                                                                         zeroPeriodicBlocks,
                                                                                         variationFrequencies,
                                                                                         testTime - 0.4,
                                                                                         2,
                                                                                         0 );
                return std::make_shared< PeriodicGravityFieldVariationsParameters >(
                        variationModel, periodicCosineIndices, periodicSineIndices, bodyName );
            };

            std::shared_ptr< PolynomialGravityFieldVariationsParameters > bodyUndergoingPolynomialParameter =
                    createPolynomialVariationParameter( bodyUndergoingTorqueName );
            checkGravityFieldVariationPartial( bodyUndergoingPolynomialParameter,
                                               bodyUndergoingPolynomialParameter->getPolynomialVariationModel( ),
                                               gravityFieldOfBodyUndergoingTorque,
                                               baseCosineCoefficientsOfBodyUndergoingTorque,
                                               baseSineCoefficientsOfBodyUndergoingTorque,
                                               "dmr polynomial variations body undergoing" );
            std::shared_ptr< PeriodicGravityFieldVariationsParameters > bodyUndergoingPeriodicParameter =
                    createPeriodicVariationParameter( bodyUndergoingTorqueName );
            checkGravityFieldVariationPartial( bodyUndergoingPeriodicParameter,
                                               bodyUndergoingPeriodicParameter->getPeriodicVariationModel( ),
                                               gravityFieldOfBodyUndergoingTorque,
                                               baseCosineCoefficientsOfBodyUndergoingTorque,
                                               baseSineCoefficientsOfBodyUndergoingTorque,
                                               "dmr periodic variations body undergoing" );
            std::shared_ptr< PolynomialGravityFieldVariationsParameters > bodyExertingPolynomialParameter =
                    createPolynomialVariationParameter( bodyExertingTorqueName );
            checkGravityFieldVariationPartial( bodyExertingPolynomialParameter,
                                               bodyExertingPolynomialParameter->getPolynomialVariationModel( ),
                                               gravityFieldOfBodyExertingTorque,
                                               baseCosineCoefficientsOfBodyExertingTorque,
                                               baseSineCoefficientsOfBodyExertingTorque,
                                               "dmr polynomial variations body exerting" );
            std::shared_ptr< PeriodicGravityFieldVariationsParameters > bodyExertingPeriodicParameter =
                    createPeriodicVariationParameter( bodyExertingTorqueName );
            checkGravityFieldVariationPartial( bodyExertingPeriodicParameter,
                                               bodyExertingPeriodicParameter->getPeriodicVariationModel( ),
                                               gravityFieldOfBodyExertingTorque,
                                               baseCosineCoefficientsOfBodyExertingTorque,
                                               baseSineCoefficientsOfBodyExertingTorque,
                                               "dmr periodic variations body exerting" );
        }

        // This check verifies that the analytical quaternion partial of the body undergoing torque is not identically zero.
        BOOST_CHECK_GT( analyticalPartialWrtOrientationOfBodyUndergoingTorque.norm( ), 1.0E-20 );
        // This check verifies that the analytical quaternion partial of the body exerting torque is not identically zero.
        BOOST_CHECK_GT( analyticalPartialWrtOrientationOfBodyExertingTorque.norm( ), 1.0E-20 );
        // This check verifies that the analytical position partial of the body undergoing torque is not identically zero.
        BOOST_CHECK_GT( analyticalPartialWrtTranslationalStateOfBodyUndergoingTorque.block( 0, 0, 3, 3 ).norm( ), 1.0E-20 );
        // This check verifies that the analytical position partial of the body exerting torque is not identically zero.
        BOOST_CHECK_GT( analyticalPartialWrtTranslationalStateOfBodyExertingTorque.block( 0, 0, 3, 3 ).norm( ), 1.0E-20 );
        // This check verifies that the numerical quaternion reference for the body undergoing torque is not identically zero.
        BOOST_CHECK_GT( torqueDeviationDueToOrientationChangeOfBodyUndergoingTorque.norm( ), 1.0E-20 );
        // This check verifies that the numerical quaternion reference for the body exerting torque is not identically zero.
        BOOST_CHECK_GT( torqueDeviationDueToOrientationChangeOfBodyExertingTorque.norm( ), 1.0E-20 );
        // This check verifies that the numerical position partial of the body undergoing torque is not identically zero.
        BOOST_CHECK_GT( numericalPartialWrtPositionOfBodyUndergoingTorque.norm( ), 1.0E-20 );
        // This check verifies that the numerical position partial of the body exerting torque is not identically zero.
        BOOST_CHECK_GT( numericalPartialWrtPositionOfBodyExertingTorque.norm( ), 1.0E-20 );
        // This check verifies that the analytical cosine-coefficient partial of the body undergoing torque is not identically zero.
        BOOST_CHECK_GT( analyticalPartialWrtCosineCoefficientsOfBodyUndergoingTorque.norm( ), 1.0E-20 );
        // This check verifies that the analytical sine-coefficient partial of the body undergoing torque is not identically zero.
        BOOST_CHECK_GT( analyticalPartialWrtSineCoefficientsOfBodyUndergoingTorque.norm( ), 1.0E-20 );
        // This check verifies that the analytical cosine-coefficient partial of the body exerting torque is not identically zero.
        BOOST_CHECK_GT( analyticalPartialWrtCosineCoefficientsOfBodyExertingTorque.norm( ), 1.0E-20 );
        // This check verifies that the analytical sine-coefficient partial of the body exerting torque is not identically zero.
        BOOST_CHECK_GT( analyticalPartialWrtSineCoefficientsOfBodyExertingTorque.norm( ), 1.0E-20 );
        // This check verifies that the numerical cosine-coefficient partial of the body undergoing torque is not identically zero.
        BOOST_CHECK_GT( numericalPartialWrtCosineCoefficientsOfBodyUndergoingTorque.norm( ), 1.0E-20 );
        // This check verifies that the numerical sine-coefficient partial of the body undergoing torque is not identically zero.
        BOOST_CHECK_GT( numericalPartialWrtSineCoefficientsOfBodyUndergoingTorque.norm( ), 1.0E-20 );
        // This check verifies that the numerical cosine-coefficient partial of the body exerting torque is not identically zero.
        BOOST_CHECK_GT( numericalPartialWrtCosineCoefficientsOfBodyExertingTorque.norm( ), 1.0E-20 );
        // This check verifies that the numerical sine-coefficient partial of the body exerting torque is not identically zero.
        BOOST_CHECK_GT( numericalPartialWrtSineCoefficientsOfBodyExertingTorque.norm( ), 1.0E-20 );

        for( int index = 1; index < 4; index++ )
        {
            // Quaternion-element perturbation mapping:
            // analytical delta_tau = (d tau / d q) * delta_q must reproduce the directly perturbed torque delta.
            const Eigen::Vector3d analyticalTorqueDeviationOfBodyUndergoingTorque =
                    analyticalPartialWrtOrientationOfBodyUndergoingTorque * appliedQuaternionPerturbationOfBodyUndergoingTorque.at( index );
            const Eigen::Vector3d analyticalTorqueDeviationOfBodyExertingTorque =
                    analyticalPartialWrtOrientationOfBodyExertingTorque * appliedQuaternionPerturbationOfBodyExertingTorque.at( index );

            checkMatrixClosePerElement( analyticalTorqueDeviationOfBodyUndergoingTorque,
                                        torqueDeviationDueToOrientationChangeOfBodyUndergoingTorque.col( index - 1 ),
                                        1.0E-5,
                                        "dmr orientation body undergoing, perturbation " + std::to_string( index ) );
            checkMatrixClosePerElement( analyticalTorqueDeviationOfBodyExertingTorque,
                                        torqueDeviationDueToOrientationChangeOfBodyExertingTorque.col( index - 1 ),
                                        1.0E-5,
                                        "dmr orientation body exerting, perturbation " + std::to_string( index ) );
        }
        // Verify body-undergoing position partials against finite differences.
        checkMatrixClosePerElement( analyticalPartialWrtTranslationalStateOfBodyUndergoingTorque.block( 0, 0, 3, 3 ),
                                    numericalPartialWrtPositionOfBodyUndergoingTorque,
                                    3.0E-7,
                                    "dmr position body undergoing" );
        // Verify body-exerting position partials against finite differences.
        checkMatrixClosePerElement( analyticalPartialWrtTranslationalStateOfBodyExertingTorque.block( 0, 0, 3, 3 ),
                                    numericalPartialWrtPositionOfBodyExertingTorque,
                                    3.0E-7,
                                    "dmr position body exerting" );
        // Verify body-undergoing cosine coefficient partials against finite differences.
        checkMatrixClosePerElement( analyticalPartialWrtCosineCoefficientsOfBodyUndergoingTorque,
                                    numericalPartialWrtCosineCoefficientsOfBodyUndergoingTorque,
                                    1.0E-8,
                                    "dmr cosine coefficients body undergoing" );
        // Verify body-undergoing sine coefficient partials against finite differences.
        checkMatrixClosePerElement( analyticalPartialWrtSineCoefficientsOfBodyUndergoingTorque,
                                    numericalPartialWrtSineCoefficientsOfBodyUndergoingTorque,
                                    1.0E-8,
                                    "dmr sine coefficients body undergoing" );
        // Verify body-exerting cosine coefficient partials against finite differences.
        checkMatrixClosePerElement( analyticalPartialWrtCosineCoefficientsOfBodyExertingTorque,
                                    numericalPartialWrtCosineCoefficientsOfBodyExertingTorque,
                                    1.0E-8,
                                    "dmr cosine coefficients body exerting" );
        // Verify body-exerting sine coefficient partials against finite differences.
        checkMatrixClosePerElement( analyticalPartialWrtSineCoefficientsOfBodyExertingTorque,
                                    numericalPartialWrtSineCoefficientsOfBodyExertingTorque,
                                    1.0E-8,
                                    "dmr sine coefficients body exerting" );

        makeBodyMassIndependentOfGravityField( bodyUndergoingTorque );
        makeBodyMassIndependentOfGravityField( bodyExertingTorque );
        torqueModel->resetCurrentTime( );
        torquePartial->resetCurrentTime( );
        torqueModel->updateMembers( testTime );
        torquePartial->update( testTime );

        const Eigen::Vector3d analyticalPartialWrtGravitationalParameterOfBodyUndergoingTorque =
                torquePartial->wrtParameter( gravitationalParameterOfBodyUndergoingTorque );
        const Eigen::Vector3d analyticalPartialWrtGravitationalParameterOfBodyExertingTorque =
                torquePartial->wrtParameter( gravitationalParameterOfBodyExertingTorque );
        const Eigen::Vector3d numericalPartialWrtGravitationalParameterOfBodyUndergoingTorque = calculateTorqueWrtParameterPartials(
                gravitationalParameterOfBodyUndergoingTorque, torqueModel, 1.0E-2, emptyFunction, testTime );
        const Eigen::Vector3d numericalPartialWrtGravitationalParameterOfBodyExertingTorque = calculateTorqueWrtParameterPartials(
                gravitationalParameterOfBodyExertingTorque, torqueModel, 1.0E-2, emptyFunction, testTime );

        // Verify body-1 gravitational parameter has no direct effect once mass/inertia are decoupled.
        BOOST_CHECK_SMALL( analyticalPartialWrtGravitationalParameterOfBodyUndergoingTorque.norm( ), 1.0E-20 );
        // Verify the same zero body-1 gravitational-parameter response by finite differencing.
        BOOST_CHECK_SMALL( numericalPartialWrtGravitationalParameterOfBodyUndergoingTorque.norm( ), 1.0E-20 );
        // Verify body-2 gravitational parameter partial against a finite-difference reference.
        checkMatrixClosePerElement( analyticalPartialWrtGravitationalParameterOfBodyExertingTorque,
                                    numericalPartialWrtGravitationalParameterOfBodyExertingTorque,
                                    1.0E-9,
                                    "dmr GM body exerting" );

        std::shared_ptr< InitialMassStateParameter< double > > massOfBodyUndergoingTorqueParameter =
                makeInitialMassParameter( bodyUndergoingTorque, bodyUndergoingTorqueName );
        std::shared_ptr< InitialMassStateParameter< double > > massOfBodyExertingTorqueParameter =
                makeInitialMassParameter( bodyExertingTorque, bodyExertingTorqueName );
        const Eigen::MatrixXd analyticalPartialWrtMassOfBodyUndergoingTorque =
                torquePartial->wrtParameter( massOfBodyUndergoingTorqueParameter );
        const Eigen::MatrixXd analyticalPartialWrtMassOfBodyExertingTorque =
                torquePartial->wrtParameter( massOfBodyExertingTorqueParameter );
        const Eigen::MatrixXd numericalPartialWrtMassOfBodyUndergoingTorque = calculateTorqueWrtParameterPartials(
                massOfBodyUndergoingTorqueParameter,
                torqueModel,
                Eigen::VectorXd::Constant( massOfBodyUndergoingTorqueParameter->getParameterSize( ), 1.0E6 ),
                emptyFunction,
                testTime );
        const Eigen::MatrixXd numericalPartialWrtMassOfBodyExertingTorque = calculateTorqueWrtParameterPartials(
                massOfBodyExertingTorqueParameter,
                torqueModel,
                Eigen::VectorXd::Constant( massOfBodyExertingTorqueParameter->getParameterSize( ), 1.0E6 ),
                emptyFunction,
                testTime );

        // Verify body-1 mass partial against a finite-difference reference.
        checkMatrixClosePerElement( analyticalPartialWrtMassOfBodyUndergoingTorque,
                                    numericalPartialWrtMassOfBodyUndergoingTorque,
                                    1.0E-10,
                                    "dmr mass body undergoing" );
        // Verify body-2 mass has no direct effect in this parameterization when gravitational parameter is separate.
        BOOST_CHECK_SMALL( analyticalPartialWrtMassOfBodyExertingTorque.norm( ), 1.0E-20 );
        // Verify the same zero body-2 mass response by finite differencing.
        BOOST_CHECK_SMALL( numericalPartialWrtMassOfBodyExertingTorque.norm( ), 1.0E-20 );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
