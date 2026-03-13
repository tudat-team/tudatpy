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

#include <limits>
#include <string>

#include <boost/test/unit_test.hpp>

#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/orbit_determination/acceleration_partials/numericalAccelerationPartial.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameterSet.h"
#include "tudat/basics/testMacros.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
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

SystemOfBodies createTwoBodyTorquePartialTestSystem(
        const double testTime,
        const std::string& bodyUndergoingTorqueName,
        const std::string& bodyExertingTorqueName,
        const bool useArbitraryRotationStates )
{
    SystemOfBodies bodies = SystemOfBodies( "SSB", "ECLIPJ2000" );
    bodies.createEmptyBody( bodyUndergoingTorqueName );
    bodies.createEmptyBody( bodyExertingTorqueName );

    std::shared_ptr< Body > bodyUndergoingTorque = bodies.at( bodyUndergoingTorqueName );
    std::shared_ptr< Body > bodyExertingTorque = bodies.at( bodyExertingTorqueName );

    double scaledMeanMomentOfInertiaBodyUndergoingTorque = 0.3;
    double scaledMeanMomentOfInertiaBodyExertingTorque = 0.26;
    const double massBodyUndergoingTorque = 6.0E11;
    const double massBodyExertingTorque = 4.0E11;
    const double gravitationalParameterBodyUndergoingTorque =
            massBodyUndergoingTorque * physical_constants::GRAVITATIONAL_CONSTANT;
    const double gravitationalParameterBodyExertingTorque =
            massBodyExertingTorque * physical_constants::GRAVITATIONAL_CONSTANT;
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

    const Eigen::Matrix3d inertiaTensorOfBodyUndergoingTorque =
            normalizedInertiaTensorOfBodyUndergoingTorque * massBodyUndergoingTorque *
            referenceRadiusBodyUndergoingTorque * referenceRadiusBodyUndergoingTorque;
    const Eigen::Matrix3d inertiaTensorOfBodyExertingTorque =
            normalizedInertiaTensorOfBodyExertingTorque * massBodyExertingTorque *
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

    bodyUndergoingTorque->setGravityFieldModel( std::make_shared< gravitation::SphericalHarmonicsGravityField >(
            gravitationalParameterBodyUndergoingTorque,
            referenceRadiusBodyUndergoingTorque,
            cosineCoefficientsOfBodyUndergoingTorque,
            sineCoefficientsOfBodyUndergoingTorque,
            bodyUndergoingTorqueName + "_Fixed",
            scaledMeanMomentOfInertiaBodyUndergoingTorque ) );
    bodyExertingTorque->setGravityFieldModel( std::make_shared< gravitation::SphericalHarmonicsGravityField >(
            gravitationalParameterBodyExertingTorque,
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
                        Eigen::AngleAxisd( 0.6, Eigen::Vector3d::UnitX( ) ) *
                        Eigen::AngleAxisd( -0.25, Eigen::Vector3d::UnitY( ) ) *
                        Eigen::AngleAxisd( 0.4, Eigen::Vector3d::UnitZ( ) ) ) );
        Eigen::Vector7d arbitraryRotationStateOfBodyExertingTorque = Eigen::Vector7d::Zero( );
        arbitraryRotationStateOfBodyExertingTorque.segment( 0, 4 ) =
                tudat::linear_algebra::convertQuaternionToVectorFormat( Eigen::Quaterniond(
                        Eigen::AngleAxisd( -0.2, Eigen::Vector3d::UnitX( ) ) *
                        Eigen::AngleAxisd( 0.45, Eigen::Vector3d::UnitY( ) ) *
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

void checkMatrixClosePerElement( const Eigen::MatrixXd& analyticalMatrix, const Eigen::MatrixXd& numericalMatrix, const double tolerance )
{
    BOOST_REQUIRE_EQUAL( analyticalMatrix.rows( ), numericalMatrix.rows( ) );
    BOOST_REQUIRE_EQUAL( analyticalMatrix.cols( ), numericalMatrix.cols( ) );

    for( int i = 0; i < analyticalMatrix.rows( ); i++ )
    {
        for( int j = 0; j < analyticalMatrix.cols( ); j++ )
        {
            BOOST_CHECK_SMALL( std::fabs( analyticalMatrix( i, j ) - numericalMatrix( i, j ) ) /
                                       std::max( 1.0, std::fabs( numericalMatrix( i, j ) ) ),
                               tolerance );
        }
    }
}

Eigen::Matrix< double, 6, 1 > getAuxiliaryFunctionVector(
        const acceleration_partials::detail::FourthDegreeTorqueAuxiliaryQuantities& auxiliaryQuantities )
{
    Eigen::Matrix< double, 6, 1 > auxiliaryFunctionVector;
    auxiliaryFunctionVector << auxiliaryQuantities.fyzFunction,
            auxiliaryQuantities.fxzFunction,
            auxiliaryQuantities.fxyFunction,
            auxiliaryQuantities.gyzFunction,
            auxiliaryQuantities.gxzFunction,
            auxiliaryQuantities.gxyFunction;
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

} // namespace

BOOST_AUTO_TEST_SUITE( test_full_two_body_torque_partials )

BOOST_AUTO_TEST_CASE( testFourthDegreeTorqueAuxiliaryFunctionPartials )
{
    const Eigen::Vector3d relativePositionInBodyFixedFrame( 3275.0, -1840.0, 2510.0 );
    const double massOfBodyExertingTorque = 4.0E11;
    const Eigen::Matrix3d inertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque =
            ( Eigen::Matrix3d( ) << 2.6E16, 1.2E15, -9.0E14,
              1.2E15, 3.1E16, -1.6E15,
             -9.0E14, -1.6E15, 3.5E16 ).finished( );
    const Eigen::Matrix< double, 6, 1 > independentInertiaTensorComponentsOfBodyExertingTorque =
            acceleration_partials::detail::getIndependentInertiaTensorComponentsFromMatrix(
                    inertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque );
    const acceleration_partials::detail::FourthDegreeTorqueAuxiliaryQuantities auxiliaryQuantities =
            acceleration_partials::detail::computeFourthDegreeTorqueAuxiliaryQuantities(
                    relativePositionInBodyFixedFrame,
                    massOfBodyExertingTorque,
                    independentInertiaTensorComponentsOfBodyExertingTorque );

    const double positionPerturbation = 1.0E-3;
    for( int coordinateIndex = 0; coordinateIndex < 3; coordinateIndex++ )
    {
        Eigen::Vector3d upperturbedPosition = relativePositionInBodyFixedFrame;
        upperturbedPosition( coordinateIndex ) += positionPerturbation;
        const auto upPerturbedAuxiliaryQuantities =
                acceleration_partials::detail::computeFourthDegreeTorqueAuxiliaryQuantities(
                        upperturbedPosition,
                        massOfBodyExertingTorque,
                        independentInertiaTensorComponentsOfBodyExertingTorque );

        Eigen::Vector3d downperturbedPosition = relativePositionInBodyFixedFrame;
        downperturbedPosition( coordinateIndex ) -= positionPerturbation;
        const auto downPerturbedAuxiliaryQuantities =
                acceleration_partials::detail::computeFourthDegreeTorqueAuxiliaryQuantities(
                        downperturbedPosition,
                        massOfBodyExertingTorque,
                        independentInertiaTensorComponentsOfBodyExertingTorque );

        const double numericalPartialOfContractedInertiaTensorWrtCoordinate =
                ( upPerturbedAuxiliaryQuantities.contractedInertiaTensorOfBodyExertingTorque -
                  downPerturbedAuxiliaryQuantities.contractedInertiaTensorOfBodyExertingTorque ) /
                ( 2.0 * positionPerturbation );
        const double analyticalPartialOfContractedInertiaTensorWrtCoordinate =
                acceleration_partials::detail::computePartialOfContractedInertiaTensorOfBodyExertingTorqueWrtCoordinate(
                        auxiliaryQuantities, coordinateIndex );

        BOOST_CHECK_SMALL(
                std::fabs( analyticalPartialOfContractedInertiaTensorWrtCoordinate -
                           numericalPartialOfContractedInertiaTensorWrtCoordinate ) /
                        std::max( 1.0, std::fabs( numericalPartialOfContractedInertiaTensorWrtCoordinate ) ),
                5.0E-9 );

        const Eigen::Matrix< double, 6, 1 > numericalPartialOfAuxiliaryFunctionsWrtCoordinate =
                ( getAuxiliaryFunctionVector( upPerturbedAuxiliaryQuantities ) -
                  getAuxiliaryFunctionVector( downPerturbedAuxiliaryQuantities ) ) /
                ( 2.0 * positionPerturbation );
        const Eigen::Matrix< double, 6, 1 > analyticalPartialOfAuxiliaryFunctionsWrtCoordinate =
                acceleration_partials::detail::computePartialOfAuxiliaryFunctionsWrtPositionCoordinate(
                        auxiliaryQuantities, coordinateIndex );

        checkMatrixClosePerElement(
                analyticalPartialOfAuxiliaryFunctionsWrtCoordinate,
                numericalPartialOfAuxiliaryFunctionsWrtCoordinate,
                5.0E-7 );
    }

    const double inertiaTensorComponentPerturbation = 1.0E7;
    Eigen::Matrix< double, 6, 6 > numericalPartialOfAuxiliaryFunctionsWrtInertiaTensorComponents =
            Eigen::Matrix< double, 6, 6 >::Zero( );
    for( int componentIndex = 0; componentIndex < 6; componentIndex++ )
    {
        Eigen::Matrix< double, 6, 1 > upPerturbedIndependentInertiaTensorComponents =
                independentInertiaTensorComponentsOfBodyExertingTorque;
        upPerturbedIndependentInertiaTensorComponents( componentIndex ) += inertiaTensorComponentPerturbation;
        const auto upPerturbedAuxiliaryQuantities =
                acceleration_partials::detail::computeFourthDegreeTorqueAuxiliaryQuantities(
                        relativePositionInBodyFixedFrame,
                        massOfBodyExertingTorque,
                        upPerturbedIndependentInertiaTensorComponents );

        Eigen::Matrix< double, 6, 1 > downPerturbedIndependentInertiaTensorComponents =
                independentInertiaTensorComponentsOfBodyExertingTorque;
        downPerturbedIndependentInertiaTensorComponents( componentIndex ) -= inertiaTensorComponentPerturbation;
        const auto downPerturbedAuxiliaryQuantities =
                acceleration_partials::detail::computeFourthDegreeTorqueAuxiliaryQuantities(
                        relativePositionInBodyFixedFrame,
                        massOfBodyExertingTorque,
                        downPerturbedIndependentInertiaTensorComponents );

        numericalPartialOfAuxiliaryFunctionsWrtInertiaTensorComponents.col( componentIndex ) =
                ( getAuxiliaryFunctionVector( upPerturbedAuxiliaryQuantities ) -
                  getAuxiliaryFunctionVector( downPerturbedAuxiliaryQuantities ) ) /
                ( 2.0 * inertiaTensorComponentPerturbation );
    }

    const Eigen::Matrix< double, 6, 6 > analyticalPartialOfAuxiliaryFunctionsWrtInertiaTensorComponents =
            acceleration_partials::detail::computePartialOfAuxiliaryFunctionsWrtIndependentInertiaTensorComponentsOfBodyExertingTorque(
                    auxiliaryQuantities );

    checkMatrixClosePerElement(
            analyticalPartialOfAuxiliaryFunctionsWrtInertiaTensorComponents,
            numericalPartialOfAuxiliaryFunctionsWrtInertiaTensorComponents,
            5.0E-5 );
}

BOOST_AUTO_TEST_CASE( testFullTwoBodySphericalHarmonicGravitationalTorquePartials )
{
    const std::string bodyUndergoingTorqueName = "BodyA";
    const std::string bodyExertingTorqueName = "BodyB";
    const double testTime = 1250.0;

    for( const bool useArbitraryRotationStates : { false, true } )
    {
        SystemOfBodies bodies = createTwoBodyTorquePartialTestSystem(
                testTime, bodyUndergoingTorqueName, bodyExertingTorqueName, useArbitraryRotationStates );
        std::shared_ptr< Body > bodyUndergoingTorque = bodies.at( bodyUndergoingTorqueName );
        std::shared_ptr< Body > bodyExertingTorque = bodies.at( bodyExertingTorqueName );

        std::shared_ptr< TorqueModel > torqueModel = createFullTwoBodySphericalHarmonicGravitationalTorqueModel(
                bodyUndergoingTorque,
                bodyExertingTorque,
                fullTwoBodySphericalHarmonicGravitationalTorque( 2, 2, 2, 2 ),
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

        std::shared_ptr< TorquePartial > torquePartial = createAnalyticalTorquePartial(
                torqueModel,
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

        const Eigen::Vector4d orientationPerturbation = Eigen::Vector4d::Constant( 5.0E-8 );
        const Eigen::Vector3d positionPerturbation = Eigen::Vector3d::Constant( 2.0 );
        const Eigen::VectorXd cosineCoefficientPerturbationBodyUndergoingTorque =
                Eigen::VectorXd::Constant( cosineCoefficientsOfBodyUndergoingTorqueParameter->getParameterSize( ), 2.0E-8 );
        const Eigen::VectorXd sineCoefficientPerturbationBodyUndergoingTorque =
                Eigen::VectorXd::Constant( sineCoefficientsOfBodyUndergoingTorqueParameter->getParameterSize( ), 2.0E-8 );
        const Eigen::VectorXd cosineCoefficientPerturbationBodyExertingTorque =
                Eigen::VectorXd::Constant( cosineCoefficientsOfBodyExertingTorqueParameter->getParameterSize( ), 2.0E-8 );
        const Eigen::VectorXd sineCoefficientPerturbationBodyExertingTorque =
                Eigen::VectorXd::Constant( sineCoefficientsOfBodyExertingTorqueParameter->getParameterSize( ), 2.0E-8 );

        const std::function< void( Eigen::Vector7d ) > setRotationalStateOfBodyUndergoingTorque =
                std::bind( &Body::setCurrentRotationalStateToLocalFrame, bodyUndergoingTorque, std::placeholders::_1 );
        const std::function< void( Eigen::Vector7d ) > setRotationalStateOfBodyExertingTorque =
                std::bind( &Body::setCurrentRotationalStateToLocalFrame, bodyExertingTorque, std::placeholders::_1 );
        const std::function< void( Eigen::Vector6d ) > setTranslationalStateOfBodyUndergoingTorque =
                std::bind( &Body::setState, bodyUndergoingTorque, std::placeholders::_1 );
        const std::function< void( Eigen::Vector6d ) > setTranslationalStateOfBodyExertingTorque =
                std::bind( &Body::setState, bodyExertingTorque, std::placeholders::_1 );

        const Eigen::MatrixXd numericalPartialWrtOrientationOfBodyUndergoingTorque = calculateTorqueWrtRotationalStatePartials(
                setRotationalStateOfBodyUndergoingTorque,
                torqueModel,
                bodyUndergoingTorque->getRotationalStateVector( ),
                orientationPerturbation,
                0,
                4,
                emptyFunction,
                testTime );
        const Eigen::MatrixXd numericalPartialWrtOrientationOfBodyExertingTorque = calculateTorqueWrtRotationalStatePartials(
                setRotationalStateOfBodyExertingTorque,
                torqueModel,
                bodyExertingTorque->getRotationalStateVector( ),
                orientationPerturbation,
                0,
                4,
                emptyFunction,
                testTime );
        const Eigen::MatrixXd numericalPartialWrtPositionOfBodyUndergoingTorque = calculateTorqueWrtTranslationalStatePartials(
                setTranslationalStateOfBodyUndergoingTorque,
                torqueModel,
                bodyUndergoingTorque->getState( ),
                positionPerturbation,
                0,
                emptyFunction,
                testTime );
        const Eigen::MatrixXd numericalPartialWrtPositionOfBodyExertingTorque = calculateTorqueWrtTranslationalStatePartials(
                setTranslationalStateOfBodyExertingTorque,
                torqueModel,
                bodyExertingTorque->getState( ),
                positionPerturbation,
                0,
                emptyFunction,
                testTime );

        const Eigen::MatrixXd numericalPartialWrtCosineCoefficientsOfBodyUndergoingTorque = calculateTorqueWrtParameterPartials(
                cosineCoefficientsOfBodyUndergoingTorqueParameter,
                torqueModel,
                cosineCoefficientPerturbationBodyUndergoingTorque,
                emptyFunction,
                testTime );
        const Eigen::MatrixXd numericalPartialWrtSineCoefficientsOfBodyUndergoingTorque = calculateTorqueWrtParameterPartials(
                sineCoefficientsOfBodyUndergoingTorqueParameter,
                torqueModel,
                sineCoefficientPerturbationBodyUndergoingTorque,
                emptyFunction,
                testTime );
        const Eigen::MatrixXd numericalPartialWrtCosineCoefficientsOfBodyExertingTorque = calculateTorqueWrtParameterPartials(
                cosineCoefficientsOfBodyExertingTorqueParameter,
                torqueModel,
                cosineCoefficientPerturbationBodyExertingTorque,
                emptyFunction,
                testTime );
        const Eigen::MatrixXd numericalPartialWrtSineCoefficientsOfBodyExertingTorque = calculateTorqueWrtParameterPartials(
                sineCoefficientsOfBodyExertingTorqueParameter,
                torqueModel,
                sineCoefficientPerturbationBodyExertingTorque,
                emptyFunction,
                testTime );

        const double analyticalStatePartialNorm =
                analyticalPartialWrtOrientationOfBodyUndergoingTorque.norm( ) +
                analyticalPartialWrtOrientationOfBodyExertingTorque.norm( ) +
                analyticalPartialWrtTranslationalStateOfBodyUndergoingTorque.block( 0, 0, 3, 3 ).norm( ) +
                analyticalPartialWrtTranslationalStateOfBodyExertingTorque.block( 0, 0, 3, 3 ).norm( );
        const double numericalStatePartialNorm =
                numericalPartialWrtOrientationOfBodyUndergoingTorque.norm( ) +
                numericalPartialWrtOrientationOfBodyExertingTorque.norm( ) +
                numericalPartialWrtPositionOfBodyUndergoingTorque.norm( ) +
                numericalPartialWrtPositionOfBodyExertingTorque.norm( );
        const double analyticalCoefficientPartialNorm =
                analyticalPartialWrtCosineCoefficientsOfBodyUndergoingTorque.norm( ) +
                analyticalPartialWrtSineCoefficientsOfBodyUndergoingTorque.norm( ) +
                analyticalPartialWrtCosineCoefficientsOfBodyExertingTorque.norm( ) +
                analyticalPartialWrtSineCoefficientsOfBodyExertingTorque.norm( );
        const double numericalCoefficientPartialNorm =
                numericalPartialWrtCosineCoefficientsOfBodyUndergoingTorque.norm( ) +
                numericalPartialWrtSineCoefficientsOfBodyUndergoingTorque.norm( ) +
                numericalPartialWrtCosineCoefficientsOfBodyExertingTorque.norm( ) +
                numericalPartialWrtSineCoefficientsOfBodyExertingTorque.norm( );

        BOOST_CHECK_GT( analyticalStatePartialNorm, 1.0E-20 );
        BOOST_CHECK_GT( numericalStatePartialNorm, 1.0E-20 );
        BOOST_CHECK_GT( analyticalCoefficientPartialNorm, 1.0E-20 );
        BOOST_CHECK_GT( numericalCoefficientPartialNorm, 1.0E-20 );

        checkMatrixClosePerElement( analyticalPartialWrtOrientationOfBodyUndergoingTorque,
                                    numericalPartialWrtOrientationOfBodyUndergoingTorque,
                                    5.0E-10 );
        checkMatrixClosePerElement( analyticalPartialWrtOrientationOfBodyExertingTorque,
                                    numericalPartialWrtOrientationOfBodyExertingTorque,
                                    5.0E-10 );
        checkMatrixClosePerElement( analyticalPartialWrtTranslationalStateOfBodyUndergoingTorque.block( 0, 0, 3, 3 ),
                                    numericalPartialWrtPositionOfBodyUndergoingTorque,
                                    5.0E-10 );
        checkMatrixClosePerElement( analyticalPartialWrtTranslationalStateOfBodyExertingTorque.block( 0, 0, 3, 3 ),
                                    numericalPartialWrtPositionOfBodyExertingTorque,
                                    5.0E-10 );
        checkMatrixClosePerElement( analyticalPartialWrtCosineCoefficientsOfBodyUndergoingTorque,
                                    numericalPartialWrtCosineCoefficientsOfBodyUndergoingTorque,
                                    5.0E-10 );
        checkMatrixClosePerElement( analyticalPartialWrtSineCoefficientsOfBodyUndergoingTorque,
                                    numericalPartialWrtSineCoefficientsOfBodyUndergoingTorque,
                                    5.0E-10 );
        checkMatrixClosePerElement( analyticalPartialWrtCosineCoefficientsOfBodyExertingTorque,
                                    numericalPartialWrtCosineCoefficientsOfBodyExertingTorque,
                                    5.0E-10 );
        checkMatrixClosePerElement( analyticalPartialWrtSineCoefficientsOfBodyExertingTorque,
                                    numericalPartialWrtSineCoefficientsOfBodyExertingTorque,
                                    5.0E-10 );
    }
}

BOOST_AUTO_TEST_CASE( testFourthDegreeFullTwoBodyGravitationalTorquePositionPartials )
{
    const std::string bodyUndergoingTorqueName = "BodyA";
    const std::string bodyExertingTorqueName = "BodyB";
    const double testTime = 1250.0;

    for( const bool useArbitraryRotationStates : { false, true } )
    {
        SystemOfBodies bodies = createTwoBodyTorquePartialTestSystem(
                testTime, bodyUndergoingTorqueName, bodyExertingTorqueName, useArbitraryRotationStates );
        std::shared_ptr< Body > bodyUndergoingTorque = bodies.at( bodyUndergoingTorqueName );
        std::shared_ptr< Body > bodyExertingTorque = bodies.at( bodyExertingTorqueName );

        std::shared_ptr< TorqueModel > torqueModel = createFourthDegreeFullTwoBodyGravitationalTorqueModel(
                bodyUndergoingTorque,
                bodyExertingTorque,
                fourthDegreeFullTwoBodyGravitationalTorque( ),
                bodyUndergoingTorqueName,
                bodyExertingTorqueName );

        std::shared_ptr< TorquePartial > torquePartial = createAnalyticalTorquePartial(
                torqueModel,
                std::make_pair( bodyUndergoingTorqueName, bodyUndergoingTorque ),
                std::make_pair( bodyExertingTorqueName, bodyExertingTorque ),
                basic_astrodynamics::SingleBodyTorqueModelMap( ),
                bodies );
        torquePartial->update( testTime );

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

        const Eigen::Vector3d positionPerturbation = Eigen::Vector3d::Constant( 5.0E-1 );
        const std::function< void( Eigen::Vector6d ) > setTranslationalStateOfBodyUndergoingTorque =
                std::bind( &Body::setState, bodyUndergoingTorque, std::placeholders::_1 );
        const std::function< void( Eigen::Vector6d ) > setTranslationalStateOfBodyExertingTorque =
                std::bind( &Body::setState, bodyExertingTorque, std::placeholders::_1 );

        const Eigen::MatrixXd numericalPartialWrtPositionOfBodyUndergoingTorque = calculateTorqueWrtTranslationalStatePartials(
                setTranslationalStateOfBodyUndergoingTorque,
                torqueModel,
                bodyUndergoingTorque->getState( ),
                positionPerturbation,
                0,
                emptyFunction,
                testTime );
        const Eigen::MatrixXd numericalPartialWrtPositionOfBodyExertingTorque = calculateTorqueWrtTranslationalStatePartials(
                setTranslationalStateOfBodyExertingTorque,
                torqueModel,
                bodyExertingTorque->getState( ),
                positionPerturbation,
                0,
                emptyFunction,
                testTime );

        const double analyticalPositionPartialNorm =
                analyticalPartialWrtTranslationalStateOfBodyUndergoingTorque.block( 0, 0, 3, 3 ).norm( ) +
                analyticalPartialWrtTranslationalStateOfBodyExertingTorque.block( 0, 0, 3, 3 ).norm( );
        const double numericalPositionPartialNorm =
                numericalPartialWrtPositionOfBodyUndergoingTorque.norm( ) +
                numericalPartialWrtPositionOfBodyExertingTorque.norm( );

        BOOST_CHECK_GT( analyticalPositionPartialNorm, 1.0E-20 );
        BOOST_CHECK_GT( numericalPositionPartialNorm, 1.0E-20 );

        checkMatrixClosePerElement( analyticalPartialWrtTranslationalStateOfBodyUndergoingTorque.block( 0, 0, 3, 3 ),
                                    numericalPartialWrtPositionOfBodyUndergoingTorque,
                                    5.0E-6 );
        checkMatrixClosePerElement( analyticalPartialWrtTranslationalStateOfBodyExertingTorque.block( 0, 0, 3, 3 ),
                                    numericalPartialWrtPositionOfBodyExertingTorque,
                                    5.0E-6 );
        BOOST_CHECK_SMALL( analyticalPartialWrtTranslationalStateOfBodyUndergoingTorque.block( 0, 3, 3, 3 ).norm( ), 1.0E-30 );
        BOOST_CHECK_SMALL( analyticalPartialWrtTranslationalStateOfBodyExertingTorque.block( 0, 3, 3, 3 ).norm( ), 1.0E-30 );

        const Eigen::MatrixXd numericalPartialWrtVelocityOfBodyUndergoingTorque = calculateTorqueWrtTranslationalStatePartials(
                setTranslationalStateOfBodyUndergoingTorque,
                torqueModel,
                bodyUndergoingTorque->getState( ),
                positionPerturbation,
                3,
                emptyFunction,
                testTime );
        const Eigen::MatrixXd numericalPartialWrtVelocityOfBodyExertingTorque = calculateTorqueWrtTranslationalStatePartials(
                setTranslationalStateOfBodyExertingTorque,
                torqueModel,
                bodyExertingTorque->getState( ),
                positionPerturbation,
                3,
                emptyFunction,
                testTime );
        BOOST_CHECK_SMALL( numericalPartialWrtVelocityOfBodyUndergoingTorque.norm( ), 1.0E-18 );
        BOOST_CHECK_SMALL( numericalPartialWrtVelocityOfBodyExertingTorque.norm( ), 1.0E-18 );
    }
}

BOOST_AUTO_TEST_CASE( testFourthDegreeFullTwoBodyGravitationalTorqueQuaternionPartials )
{
    const std::string bodyUndergoingTorqueName = "BodyA";
    const std::string bodyExertingTorqueName = "BodyB";
    const double testTime = 1250.0;
    const Eigen::Vector4d quaternionPerturbation = Eigen::Vector4d::Constant( 1.0E-9 );

    for( const bool useArbitraryRotationStates : { false, true } )
    {
        SystemOfBodies bodies = createTwoBodyTorquePartialTestSystem(
                testTime, bodyUndergoingTorqueName, bodyExertingTorqueName, useArbitraryRotationStates );
        std::shared_ptr< Body > bodyUndergoingTorque = bodies.at( bodyUndergoingTorqueName );
        std::shared_ptr< Body > bodyExertingTorque = bodies.at( bodyExertingTorqueName );

        std::shared_ptr< TorqueModel > torqueModel = createFourthDegreeFullTwoBodyGravitationalTorqueModel(
                bodyUndergoingTorque,
                bodyExertingTorque,
                fourthDegreeFullTwoBodyGravitationalTorque( ),
                bodyUndergoingTorqueName,
                bodyExertingTorqueName );

        std::shared_ptr< TorquePartial > torquePartial = createAnalyticalTorquePartial(
                torqueModel,
                std::make_pair( bodyUndergoingTorqueName, bodyUndergoingTorque ),
                std::make_pair( bodyExertingTorqueName, bodyExertingTorque ),
                basic_astrodynamics::SingleBodyTorqueModelMap( ),
                bodies );
        torquePartial->update( testTime );

        Eigen::MatrixXd analyticalPartialWrtOrientationOfBodyUndergoingTorque = Eigen::MatrixXd::Zero( 3, 4 );
        torquePartial->wrtOrientationOfAcceleratedBody( analyticalPartialWrtOrientationOfBodyUndergoingTorque.block( 0, 0, 3, 4 ) );
        Eigen::MatrixXd analyticalPartialWrtOrientationOfBodyExertingTorque = Eigen::MatrixXd::Zero( 3, 4 );
        torquePartial->wrtOrientationOfAcceleratingBody( analyticalPartialWrtOrientationOfBodyExertingTorque.block( 0, 0, 3, 4 ) );

        const std::function< void( Eigen::Vector7d ) > setRotationalStateOfBodyUndergoingTorque =
                std::bind( &Body::setCurrentRotationalStateToLocalFrame, bodyUndergoingTorque, std::placeholders::_1 );
        const std::function< void( Eigen::Vector7d ) > setRotationalStateOfBodyExertingTorque =
                std::bind( &Body::setCurrentRotationalStateToLocalFrame, bodyExertingTorque, std::placeholders::_1 );

        std::vector< Eigen::Vector4d > appliedQuaternionPerturbationOfBodyUndergoingTorque;
        const Eigen::MatrixXd torqueDeviationDueToOrientationChangeOfBodyUndergoingTorque =
                calculateTorqueDeviationDueToOrientationChange(
                        setRotationalStateOfBodyUndergoingTorque,
                        torqueModel,
                        bodyUndergoingTorque->getRotationalStateVector( ),
                        quaternionPerturbation,
                        appliedQuaternionPerturbationOfBodyUndergoingTorque,
                        emptyFunction,
                        testTime );

        std::vector< Eigen::Vector4d > appliedQuaternionPerturbationOfBodyExertingTorque;
        const Eigen::MatrixXd torqueDeviationDueToOrientationChangeOfBodyExertingTorque =
                calculateTorqueDeviationDueToOrientationChange(
                        setRotationalStateOfBodyExertingTorque,
                        torqueModel,
                        bodyExertingTorque->getRotationalStateVector( ),
                        quaternionPerturbation,
                        appliedQuaternionPerturbationOfBodyExertingTorque,
                        emptyFunction,
                        testTime );

        for( int index = 1; index < 4; index++ )
        {
            const Eigen::Vector3d analyticalTorqueDeviationOfBodyUndergoingTorque =
                    analyticalPartialWrtOrientationOfBodyUndergoingTorque *
                    appliedQuaternionPerturbationOfBodyUndergoingTorque.at( index );
            const Eigen::Vector3d numericalTorqueDeviationOfBodyUndergoingTorque =
                    torqueDeviationDueToOrientationChangeOfBodyUndergoingTorque.col( index - 1 );
            checkMatrixClosePerElement(
                    analyticalTorqueDeviationOfBodyUndergoingTorque,
                    numericalTorqueDeviationOfBodyUndergoingTorque,
                    5.0E-6 );

            const Eigen::Vector3d analyticalTorqueDeviationOfBodyExertingTorque =
                    analyticalPartialWrtOrientationOfBodyExertingTorque *
                    appliedQuaternionPerturbationOfBodyExertingTorque.at( index );
            const Eigen::Vector3d numericalTorqueDeviationOfBodyExertingTorque =
                    torqueDeviationDueToOrientationChangeOfBodyExertingTorque.col( index - 1 );
            checkMatrixClosePerElement(
                    analyticalTorqueDeviationOfBodyExertingTorque,
                    numericalTorqueDeviationOfBodyExertingTorque,
                    5.0E-6 );
        }
    }
}

BOOST_AUTO_TEST_CASE( testFourthDegreeFullTwoBodyGravitationalTorqueCoefficientPartials )
{
    const std::string bodyUndergoingTorqueName = "BodyA";
    const std::string bodyExertingTorqueName = "BodyB";
    const double testTime = 1250.0;

    for( const bool useArbitraryRotationStates : { false, true } )
    {
        SystemOfBodies bodies = createTwoBodyTorquePartialTestSystem(
                testTime, bodyUndergoingTorqueName, bodyExertingTorqueName, useArbitraryRotationStates );
        std::shared_ptr< Body > bodyUndergoingTorque = bodies.at( bodyUndergoingTorqueName );
        std::shared_ptr< Body > bodyExertingTorque = bodies.at( bodyExertingTorqueName );

        std::shared_ptr< TorqueModel > torqueModel = createFourthDegreeFullTwoBodyGravitationalTorqueModel(
                bodyUndergoingTorque,
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

        std::shared_ptr< TorquePartial > torquePartial = createAnalyticalTorquePartial(
                torqueModel,
                std::make_pair( bodyUndergoingTorqueName, bodyUndergoingTorque ),
                std::make_pair( bodyExertingTorqueName, bodyExertingTorque ),
                basic_astrodynamics::SingleBodyTorqueModelMap( ),
                bodies,
                parameterSet );
        torquePartial->update( testTime );

        const Eigen::MatrixXd analyticalPartialWrtCosineCoefficientsOfBodyUndergoingTorque =
                torquePartial->wrtParameter( cosineCoefficientsOfBodyUndergoingTorqueParameter );
        const Eigen::MatrixXd analyticalPartialWrtSineCoefficientsOfBodyUndergoingTorque =
                torquePartial->wrtParameter( sineCoefficientsOfBodyUndergoingTorqueParameter );
        const Eigen::MatrixXd analyticalPartialWrtCosineCoefficientsOfBodyExertingTorque =
                torquePartial->wrtParameter( cosineCoefficientsOfBodyExertingTorqueParameter );
        const Eigen::MatrixXd analyticalPartialWrtSineCoefficientsOfBodyExertingTorque =
                torquePartial->wrtParameter( sineCoefficientsOfBodyExertingTorqueParameter );

        const Eigen::VectorXd cosineCoefficientPerturbationBodyUndergoingTorque =
                Eigen::VectorXd::Constant( cosineCoefficientsOfBodyUndergoingTorqueParameter->getParameterSize( ), 2.0E-8 );
        const Eigen::VectorXd sineCoefficientPerturbationBodyUndergoingTorque =
                Eigen::VectorXd::Constant( sineCoefficientsOfBodyUndergoingTorqueParameter->getParameterSize( ), 2.0E-8 );
        const Eigen::VectorXd cosineCoefficientPerturbationBodyExertingTorque =
                Eigen::VectorXd::Constant( cosineCoefficientsOfBodyExertingTorqueParameter->getParameterSize( ), 2.0E-8 );
        const Eigen::VectorXd sineCoefficientPerturbationBodyExertingTorque =
                Eigen::VectorXd::Constant( sineCoefficientsOfBodyExertingTorqueParameter->getParameterSize( ), 2.0E-8 );

        const std::function< void( ) > updateFunction =
                std::bind( &updateBodyMassDistributions, bodyUndergoingTorque, bodyExertingTorque, testTime );

        const Eigen::MatrixXd numericalPartialWrtCosineCoefficientsOfBodyUndergoingTorque = calculateTorqueWrtParameterPartials(
                cosineCoefficientsOfBodyUndergoingTorqueParameter,
                torqueModel,
                cosineCoefficientPerturbationBodyUndergoingTorque,
                updateFunction,
                testTime );
        const Eigen::MatrixXd numericalPartialWrtSineCoefficientsOfBodyUndergoingTorque = calculateTorqueWrtParameterPartials(
                sineCoefficientsOfBodyUndergoingTorqueParameter,
                torqueModel,
                sineCoefficientPerturbationBodyUndergoingTorque,
                updateFunction,
                testTime );
        const Eigen::MatrixXd numericalPartialWrtCosineCoefficientsOfBodyExertingTorque = calculateTorqueWrtParameterPartials(
                cosineCoefficientsOfBodyExertingTorqueParameter,
                torqueModel,
                cosineCoefficientPerturbationBodyExertingTorque,
                updateFunction,
                testTime );
        const Eigen::MatrixXd numericalPartialWrtSineCoefficientsOfBodyExertingTorque = calculateTorqueWrtParameterPartials(
                sineCoefficientsOfBodyExertingTorqueParameter,
                torqueModel,
                sineCoefficientPerturbationBodyExertingTorque,
                updateFunction,
                testTime );

        const double analyticalCoefficientPartialNorm =
                analyticalPartialWrtCosineCoefficientsOfBodyUndergoingTorque.norm( ) +
                analyticalPartialWrtSineCoefficientsOfBodyUndergoingTorque.norm( ) +
                analyticalPartialWrtCosineCoefficientsOfBodyExertingTorque.norm( ) +
                analyticalPartialWrtSineCoefficientsOfBodyExertingTorque.norm( );
        const double numericalCoefficientPartialNorm =
                numericalPartialWrtCosineCoefficientsOfBodyUndergoingTorque.norm( ) +
                numericalPartialWrtSineCoefficientsOfBodyUndergoingTorque.norm( ) +
                numericalPartialWrtCosineCoefficientsOfBodyExertingTorque.norm( ) +
                numericalPartialWrtSineCoefficientsOfBodyExertingTorque.norm( );

        BOOST_CHECK_GT( analyticalCoefficientPartialNorm, 1.0E-20 );
        BOOST_CHECK_GT( numericalCoefficientPartialNorm, 1.0E-20 );

        checkMatrixClosePerElement( analyticalPartialWrtCosineCoefficientsOfBodyUndergoingTorque,
                                    numericalPartialWrtCosineCoefficientsOfBodyUndergoingTorque,
                                    2.0E-5 );
        checkMatrixClosePerElement( analyticalPartialWrtSineCoefficientsOfBodyUndergoingTorque,
                                    numericalPartialWrtSineCoefficientsOfBodyUndergoingTorque,
                                    2.0E-5 );
        checkMatrixClosePerElement( analyticalPartialWrtCosineCoefficientsOfBodyExertingTorque,
                                    numericalPartialWrtCosineCoefficientsOfBodyExertingTorque,
                                    2.0E-5 );
        checkMatrixClosePerElement( analyticalPartialWrtSineCoefficientsOfBodyExertingTorque,
                                    numericalPartialWrtSineCoefficientsOfBodyExertingTorque,
                                    2.0E-5 );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
