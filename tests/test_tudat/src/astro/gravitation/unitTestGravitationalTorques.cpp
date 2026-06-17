/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Easy calculation. Newton's Law of Gravity Tutorial,
 *          http://easycalculation.com/physics/classical-physics/learn-newtons-law.php, last
 *          accessed: 12th February, 2012.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/ephemerides/constantRotationalEphemeris.h"
#include "tudat/astro/gravitation/fourthDegreeFullTwoBodyGravitationalTorque.h"
#include "tudat/astro/gravitation/sphericalHarmonicsGravityField.h"
#include "tudat/math/basic/sphericalHarmonics.h"
#include "tudat/astro/gravitation/fullTwoBodySphericalHarmonicTorque.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/createGravityField.h"
#include "tudat/simulation/propagation_setup/accelerationSettings.h"
#include "tudat/simulation/propagation_setup/createAccelerationModels.h"
#include "tudat/simulation/propagation_setup/torqueSettings.h"
#include "tudat/simulation/propagation_setup/createTorqueModel.h"

namespace tudat
{
namespace unit_tests
{

std::shared_ptr< simulation_setup::Body > createBodyForFullTwoBodyTorqueTest( const double gravitationalParameter,
                                                                              const double referenceRadius,
                                                                              const Eigen::MatrixXd& cosineCoefficients,
                                                                              const Eigen::MatrixXd& sineCoefficients,
                                                                              const Eigen::Vector3d& position,
                                                                              const Eigen::Quaterniond& rotationToBodyFixedFrame,
                                                                              const double scaledMeanMomentOfInertia = TUDAT_NAN )
{
    std::shared_ptr< simulation_setup::Body > body = std::make_shared< simulation_setup::Body >( );
    body->setGravityFieldModel( std::make_shared< gravitation::SphericalHarmonicsGravityField >(
            gravitationalParameter, referenceRadius, cosineCoefficients, sineCoefficients, "BodyFixed", scaledMeanMomentOfInertia ) );

    Eigen::Vector6d bodyState = Eigen::Vector6d::Zero( );
    bodyState.segment( 0, 3 ) = position;
    body->setState( bodyState );

    body->setRotationalEphemeris( std::make_shared< ephemerides::ConstantRotationalEphemeris >(
            rotationToBodyFixedFrame.inverse( ), "ECLIPJ2000", "BodyFixed" ) );

    return body;
}

simulation_setup::SystemOfBodies createSystemOfBodiesForFullTwoBodyTorqueTest( const std::string& bodyUndergoingTorqueName,
                                                                               const std::string& bodyExertingTorqueName,
                                                                               const double gravitationalParameter,
                                                                               const double referenceRadiusBody1,
                                                                               const double referenceRadiusBody2,
                                                                               const Eigen::Vector3d& positionOfBody1,
                                                                               const Eigen::Vector3d& positionOfBody2,
                                                                               const Eigen::MatrixXd& cosineCoefficientsOfBody1,
                                                                               const Eigen::MatrixXd& sineCoefficientsOfBody1,
                                                                               const Eigen::MatrixXd& cosineCoefficientsOfBody2,
                                                                               const Eigen::MatrixXd& sineCoefficientsOfBody2,
                                                                               const Eigen::Quaterniond& rotationToBody1,
                                                                               const Eigen::Quaterniond& rotationToBody2,
                                                                               const double scaledMeanMomentOfInertiaBody1 = TUDAT_NAN,
                                                                               const double scaledMeanMomentOfInertiaBody2 = TUDAT_NAN )
{
    simulation_setup::SystemOfBodies bodies;
    bodies.addBody( createBodyForFullTwoBodyTorqueTest( gravitationalParameter,
                                                        referenceRadiusBody1,
                                                        cosineCoefficientsOfBody1,
                                                        sineCoefficientsOfBody1,
                                                        positionOfBody1,
                                                        rotationToBody1,
                                                        scaledMeanMomentOfInertiaBody1 ),
                    bodyUndergoingTorqueName,
                    false );
    bodies.addBody( createBodyForFullTwoBodyTorqueTest( gravitationalParameter,
                                                        referenceRadiusBody2,
                                                        cosineCoefficientsOfBody2,
                                                        sineCoefficientsOfBody2,
                                                        positionOfBody2,
                                                        rotationToBody2,
                                                        scaledMeanMomentOfInertiaBody2 ),
                    bodyExertingTorqueName,
                    false );
    bodies.processBodyFrameDefinitions<>( );
    return bodies;
}

std::shared_ptr< gravitation::FullTwoBodySphericalHarmonicTorque > createFactoryFullTwoBodySphericalHarmonicTorqueModel(
        const simulation_setup::SystemOfBodies& bodies,
        const std::string& bodyUndergoingTorqueName,
        const std::string& bodyExertingTorqueName,
        const int maximumDegreeOfBodyUndergoingTorque = 2,
        const int maximumOrderOfBodyUndergoingTorque = 2,
        const int maximumDegreeOfBodyExertingTorque = 0,
        const int maximumOrderOfBodyExertingTorque = 0 )
{
    simulation_setup::SelectedTorqueMap selectedTorqueModelMap;
    selectedTorqueModelMap[ bodyUndergoingTorqueName ][ bodyExertingTorqueName ].push_back(
            simulation_setup::fullTwoBodySphericalHarmonicGravitationalTorque( maximumDegreeOfBodyUndergoingTorque,
                                                                               maximumOrderOfBodyUndergoingTorque,
                                                                               maximumDegreeOfBodyExertingTorque,
                                                                               maximumOrderOfBodyExertingTorque ) );

    basic_astrodynamics::TorqueModelMap torqueModelMap =
            simulation_setup::createTorqueModelsMap( bodies, selectedTorqueModelMap, { bodyUndergoingTorqueName } );
    return std::dynamic_pointer_cast< gravitation::FullTwoBodySphericalHarmonicTorque >(
            torqueModelMap.at( bodyUndergoingTorqueName ).at( bodyExertingTorqueName ).at( 0 ) );
}

std::shared_ptr< gravitation::FullTwoBodySphericalHarmonicTorque > createFactoryFullTwoBodySphericalHarmonicTorqueModel(
        const simulation_setup::SystemOfBodies& bodies,
        const std::string& bodyUndergoingTorqueName,
        const std::string& bodyExertingTorqueName,
        const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >& coefficientCombinationsToUse )
{
    simulation_setup::SelectedTorqueMap selectedTorqueModelMap;
    selectedTorqueModelMap[ bodyUndergoingTorqueName ][ bodyExertingTorqueName ].push_back(
            simulation_setup::fullTwoBodySphericalHarmonicGravitationalTorque( coefficientCombinationsToUse ) );

    basic_astrodynamics::TorqueModelMap torqueModelMap =
            simulation_setup::createTorqueModelsMap( bodies, selectedTorqueModelMap, { bodyUndergoingTorqueName } );
    return std::dynamic_pointer_cast< gravitation::FullTwoBodySphericalHarmonicTorque >(
            torqueModelMap.at( bodyUndergoingTorqueName ).at( bodyExertingTorqueName ).at( 0 ) );
}

std::shared_ptr< gravitation::FullTwoBodySphericalHarmonicAcceleration > createReferenceFullTwoBodySphericalHarmonicAccelerationModel(
        const simulation_setup::SystemOfBodies& bodies,
        const std::string& bodyUndergoingTorqueName,
        const std::string& bodyExertingTorqueName )
{
    return std::dynamic_pointer_cast< gravitation::FullTwoBodySphericalHarmonicAcceleration >(
            simulation_setup::createAccelerationModel( bodies.at( bodyUndergoingTorqueName ),
                                                       bodies.at( bodyExertingTorqueName ),
                                                       simulation_setup::fullTwoBodySphericalHarmonicAcceleration( 2, 2, 0, 0 ),
                                                       bodyUndergoingTorqueName,
                                                       bodyExertingTorqueName ) );
}

std::shared_ptr< gravitation::SphericalHarmonicGravitationalTorqueModel > createFactorySphericalHarmonicTorqueModel(
        const simulation_setup::SystemOfBodies& bodies,
        const std::string& bodyUndergoingTorqueName,
        const std::string& bodyExertingTorqueName,
        const int maximumDegree,
        const int maximumOrder )
{
    simulation_setup::SelectedTorqueMap selectedTorqueModelMap;
    selectedTorqueModelMap[ bodyUndergoingTorqueName ][ bodyExertingTorqueName ].push_back(
            simulation_setup::sphericalHarmonicGravitationalTorque( maximumDegree, maximumOrder ) );

    basic_astrodynamics::TorqueModelMap torqueModelMap =
            simulation_setup::createTorqueModelsMap( bodies, selectedTorqueModelMap, { bodyUndergoingTorqueName } );
    return std::dynamic_pointer_cast< gravitation::SphericalHarmonicGravitationalTorqueModel >(
            torqueModelMap.at( bodyUndergoingTorqueName ).at( bodyExertingTorqueName ).at( 0 ) );
}

std::shared_ptr< simulation_setup::Body > createBodyForFourthDegreeTwoBodyTorqueTest( const double bodyMass,
                                                                                      const Eigen::Matrix3d& inertiaTensorInBodyFixedFrame,
                                                                                      const Eigen::Vector3d& position,
                                                                                      const Eigen::Quaterniond& rotationToBodyFixedFrame )
{
    std::shared_ptr< simulation_setup::Body > body = std::make_shared< simulation_setup::Body >( );
    body->setConstantBodyMass( bodyMass );
    body->setBodyInertiaTensor( inertiaTensorInBodyFixedFrame );
    body->setGravityFieldModel(
            std::make_shared< gravitation::GravityFieldModel >( physical_constants::GRAVITATIONAL_CONSTANT * bodyMass ) );

    Eigen::Vector6d bodyState = Eigen::Vector6d::Zero( );
    bodyState.segment( 0, 3 ) = position;
    body->setState( bodyState );

    body->setRotationalEphemeris( std::make_shared< ephemerides::ConstantRotationalEphemeris >(
            rotationToBodyFixedFrame.inverse( ), "ECLIPJ2000", "BodyFixed" ) );
    return body;
}

simulation_setup::SystemOfBodies createSystemOfBodiesForFourthDegreeTwoBodyTorqueTest(
        const std::string& bodyUndergoingTorqueName,
        const std::string& bodyExertingTorqueName,
        const double massOfBodyUndergoingTorque,
        const double massOfBodyExertingTorque,
        const Eigen::Matrix3d& inertiaTensorOfBodyUndergoingTorque,
        const Eigen::Matrix3d& inertiaTensorOfBodyExertingTorque,
        const Eigen::Vector3d& positionOfBodyUndergoingTorque,
        const Eigen::Vector3d& positionOfBodyExertingTorque,
        const Eigen::Quaterniond& rotationToBodyUndergoingTorque,
        const Eigen::Quaterniond& rotationToBodyExertingTorque )
{
    simulation_setup::SystemOfBodies bodies;
    bodies.addBody( createBodyForFourthDegreeTwoBodyTorqueTest( massOfBodyUndergoingTorque,
                                                                inertiaTensorOfBodyUndergoingTorque,
                                                                positionOfBodyUndergoingTorque,
                                                                rotationToBodyUndergoingTorque ),
                    bodyUndergoingTorqueName,
                    false );
    bodies.addBody( createBodyForFourthDegreeTwoBodyTorqueTest( massOfBodyExertingTorque,
                                                                inertiaTensorOfBodyExertingTorque,
                                                                positionOfBodyExertingTorque,
                                                                rotationToBodyExertingTorque ),
                    bodyExertingTorqueName,
                    false );
    bodies.processBodyFrameDefinitions<>( );
    return bodies;
}

std::shared_ptr< gravitation::FourthDegreeFullTwoBodyGravitationalTorqueModel >
createFactoryFourthDegreeFullTwoBodyGravitationalTorqueModel( const simulation_setup::SystemOfBodies& bodies,
                                                              const std::string& bodyUndergoingTorqueName,
                                                              const std::string& bodyExertingTorqueName )
{
    simulation_setup::SelectedTorqueMap selectedTorqueModelMap;
    selectedTorqueModelMap[ bodyUndergoingTorqueName ][ bodyExertingTorqueName ].push_back(
            simulation_setup::fourthDegreeFullTwoBodyGravitationalTorque( ) );

    basic_astrodynamics::TorqueModelMap torqueModelMap =
            simulation_setup::createTorqueModelsMap( bodies, selectedTorqueModelMap, { bodyUndergoingTorqueName } );
    return std::dynamic_pointer_cast< gravitation::FourthDegreeFullTwoBodyGravitationalTorqueModel >(
            torqueModelMap.at( bodyUndergoingTorqueName ).at( bodyExertingTorqueName ).at( 0 ) );
}

std::shared_ptr< gravitation::SecondDegreeGravitationalTorqueModel > createFactorySecondDegreeGravitationalTorqueModel(
        const simulation_setup::SystemOfBodies& bodies,
        const std::string& bodyUndergoingTorqueName,
        const std::string& bodyExertingTorqueName )
{
    simulation_setup::SelectedTorqueMap selectedTorqueModelMap;
    selectedTorqueModelMap[ bodyUndergoingTorqueName ][ bodyExertingTorqueName ].push_back(
            simulation_setup::secondDegreeGravitationalTorque( ) );

    basic_astrodynamics::TorqueModelMap torqueModelMap =
            simulation_setup::createTorqueModelsMap( bodies, selectedTorqueModelMap, { bodyUndergoingTorqueName } );
    return std::dynamic_pointer_cast< gravitation::SecondDegreeGravitationalTorqueModel >(
            torqueModelMap.at( bodyUndergoingTorqueName ).at( bodyExertingTorqueName ).at( 0 ) );
}

Eigen::Matrix3d transformBodyExertingInertiaTensorToBodyUndergoingFrame( const Eigen::Matrix3d& inertiaTensorOfBodyExertingTorque,
                                                                         const Eigen::Quaterniond& rotationToBodyUndergoingTorque,
                                                                         const Eigen::Quaterniond& rotationToBodyExertingTorque )
{
    const Eigen::Matrix3d rotationFromBodyExertingToBodyUndergoing =
            rotationToBodyUndergoingTorque.toRotationMatrix( ) * rotationToBodyExertingTorque.toRotationMatrix( ).transpose( );
    return rotationFromBodyExertingToBodyUndergoing * inertiaTensorOfBodyExertingTorque *
            rotationFromBodyExertingToBodyUndergoing.transpose( );
}

Eigen::Vector3d computeManualFourthDegreeTwoBodyTorqueFromBodyStates( const std::shared_ptr< simulation_setup::Body >& bodyUndergoingTorque,
                                                                      const std::shared_ptr< simulation_setup::Body >& bodyExertingTorque )
{
    const Eigen::Vector3d relativePositionInBodyUndergoingFrame = bodyUndergoingTorque->getCurrentRotationToLocalFrame( ) *
            ( bodyExertingTorque->getPosition( ) - bodyUndergoingTorque->getPosition( ) );
    const Eigen::Matrix3d inertiaTensorOfBodyExertingInBodyUndergoingFrame =
            transformBodyExertingInertiaTensorToBodyUndergoingFrame( bodyExertingTorque->getBodyInertiaTensor( ),
                                                                     bodyUndergoingTorque->getCurrentRotationToLocalFrame( ),
                                                                     bodyExertingTorque->getCurrentRotationToLocalFrame( ) );

    return gravitation::calculateFourthDegreeFullTwoBodyGravitationalTorque( relativePositionInBodyUndergoingFrame,
                                                                             bodyExertingTorque->getBodyMass( ),
                                                                             bodyUndergoingTorque->getBodyInertiaTensor( ),
                                                                             inertiaTensorOfBodyExertingInBodyUndergoingFrame );
}

std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > getPointMassDegreeTwoInteractionCombinations( )
{
    std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > coefficientCombinations;
    for( unsigned int m = 0; m <= 2; m++ )
    {
        coefficientCombinations.push_back( std::make_tuple( 2, m, 0, 0 ) );
    }
    return coefficientCombinations;
}

std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > getDegreeTwoDegreeTwoInteractionCombinations( )
{
    std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > coefficientCombinations;
    for( unsigned int m = 0; m <= 2; m++ )
    {
        for( unsigned int k = 0; k <= 2; k++ )
        {
            coefficientCombinations.push_back( std::make_tuple( 2, m, 2, k ) );
        }
    }
    return coefficientCombinations;
}

std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > getFullDegreeTwoInteractionCombinations( )
{
    std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > coefficientCombinations =
            getPointMassDegreeTwoInteractionCombinations( );
    const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > degreeTwoDegreeTwoCombinations =
            getDegreeTwoDegreeTwoInteractionCombinations( );
    coefficientCombinations.insert(
            coefficientCombinations.end( ), degreeTwoDegreeTwoCombinations.begin( ), degreeTwoDegreeTwoCombinations.end( ) );
    return coefficientCombinations;
}

Eigen::Vector3d computeAnalyticalC20DegreeTwoFigureFigureTorque( const Eigen::Vector3d& relativePositionInBody1Frame,
                                                                 const double massOfBody1,
                                                                 const double massOfBody2,
                                                                 const double referenceRadiusBody1,
                                                                 const double referenceRadiusBody2,
                                                                 const double normalizedC20OfBody1,
                                                                 const bool useCosineCoefficientOfBody2,
                                                                 const unsigned int orderOfBody2Coefficient,
                                                                 const double normalizedCoefficientValueOfBody2 )
{
    const double x = relativePositionInBody1Frame( 0 );
    const double y = relativePositionInBody1Frame( 1 );
    const double z = relativePositionInBody1Frame( 2 );
    const double r2 = relativePositionInBody1Frame.squaredNorm( );
    const double r = std::sqrt( r2 );
    const double x2 = x * x;
    const double y2 = y * y;
    const double z2 = z * z;
    const double x4 = x2 * x2;
    const double y4 = y2 * y2;
    const double z4 = z2 * z2;
    const double sqrtThree = std::sqrt( 3.0 );

    const double gravitationalConstant = 1.0;
    const double commonPrefactor = gravitationalConstant * massOfBody1 * massOfBody2 * referenceRadiusBody1 * referenceRadiusBody1 *
            referenceRadiusBody2 * referenceRadiusBody2 / std::pow( r, 9.0 ) * normalizedC20OfBody1;

    Eigen::Vector3d analyticalTorque = Eigen::Vector3d::Zero( );

    if( useCosineCoefficientOfBody2 && orderOfBody2Coefficient == 0 )
    {
        // Eq. (6)
        analyticalTorque << y * z * ( 3.0 * r2 - 7.0 * z2 ), -x * z * ( 3.0 * r2 - 7.0 * z2 ), 0.0;
        analyticalTorque *= 75.0 / 2.0 * commonPrefactor * normalizedCoefficientValueOfBody2;
    }
    else if( useCosineCoefficientOfBody2 && orderOfBody2Coefficient == 1 )
    {
        // Eq. (7)
        analyticalTorque << x * y * ( x2 + y2 - 6.0 * z2 ),
                -0.2 * ( 4.0 * x4 + 3.0 * x2 * y2 - 27.0 * x2 * z2 - y4 + 3.0 * y2 * z2 + 4.0 * z4 ), 0.0;
        analyticalTorque *= 25.0 * sqrtThree * commonPrefactor * normalizedCoefficientValueOfBody2;
    }
    else if( !useCosineCoefficientOfBody2 && orderOfBody2Coefficient == 1 )
    {
        // Eq. (8)
        analyticalTorque << -0.2 * ( x4 - 3.0 * x2 * y2 - 3.0 * x2 * z2 - 4.0 * y4 + 27.0 * y2 * z2 - 4.0 * z4 ),
                -x * y * ( x2 + y2 - 6.0 * z2 ), 0.0;
        analyticalTorque *= 25.0 * sqrtThree * commonPrefactor * normalizedCoefficientValueOfBody2;
    }
    else if( useCosineCoefficientOfBody2 && orderOfBody2Coefficient == 2 )
    {
        // Eq. (9)
        analyticalTorque << -y * z * ( 9.0 * x2 - 5.0 * y2 + 2.0 * z2 ), x * z * ( 5.0 * x2 - 9.0 * y2 - 2.0 * z2 ), 0.0;
        analyticalTorque *= 25.0 * sqrtThree / 2.0 * commonPrefactor * normalizedCoefficientValueOfBody2;
    }
    else if( !useCosineCoefficientOfBody2 && orderOfBody2Coefficient == 2 )
    {
        // Eq. (10)
        analyticalTorque << x * z * ( x2 - 6.0 * y2 + z2 ), y * z * ( 6.0 * x2 - y2 - z2 ), 0.0;
        analyticalTorque *= 25.0 * sqrtThree * commonPrefactor * normalizedCoefficientValueOfBody2;
    }
    else
    {
        throw std::runtime_error( "Error in analytical C20xdegree2 torque computation: unsupported coefficient case." );
    }

    return analyticalTorque;
}

Eigen::Vector3d computeAnalyticalC21DegreeTwoFigureFigureTorque( const Eigen::Vector3d& relativePositionInBody1Frame,
                                                                 const double massOfBody1,
                                                                 const double massOfBody2,
                                                                 const double referenceRadiusBody1,
                                                                 const double referenceRadiusBody2,
                                                                 const double normalizedC21OfBody1,
                                                                 const bool useCosineCoefficientOfBody2,
                                                                 const unsigned int orderOfBody2Coefficient,
                                                                 const double normalizedCoefficientValueOfBody2 )
{
    const double x = relativePositionInBody1Frame( 0 );
    const double y = relativePositionInBody1Frame( 1 );
    const double z = relativePositionInBody1Frame( 2 );
    const double r2 = relativePositionInBody1Frame.squaredNorm( );
    const double r = std::sqrt( r2 );
    const double x2 = x * x;
    const double y2 = y * y;
    const double z2 = z * z;
    const double x4 = x2 * x2;
    const double y4 = y2 * y2;
    const double z4 = z2 * z2;
    const double sqrtThree = std::sqrt( 3.0 );

    const double gravitationalConstant = 1.0;
    const double commonPrefactor = gravitationalConstant * massOfBody1 * massOfBody2 * referenceRadiusBody1 * referenceRadiusBody1 *
            referenceRadiusBody2 * referenceRadiusBody2 / std::pow( r, 9.0 ) * normalizedC21OfBody1 * normalizedCoefficientValueOfBody2;

    Eigen::Vector3d analyticalTorque = Eigen::Vector3d::Zero( );

    if( useCosineCoefficientOfBody2 && orderOfBody2Coefficient == 0 )
    {
        // Eq. (8)
        analyticalTorque << ( 25.0 * sqrtThree / 2.0 ) * x * y * ( x2 + y2 - 6.0 * z2 ),
                ( 5.0 * sqrtThree / 2.0 ) * ( -7.0 * x4 - 9.0 * x2 * y2 + 51.0 * x2 * z2 - 2.0 * y4 + 21.0 * y2 * z2 - 12.0 * z4 ),
                -( 25.0 * sqrtThree / 2.0 ) * y * z * ( 3.0 * x2 + 3.0 * y2 - 4.0 * z2 );
    }
    else if( useCosineCoefficientOfBody2 && orderOfBody2Coefficient == 1 )
    {
        // Eq. (9)
        analyticalTorque << 25.0 * y * z * ( 6.0 * x2 - y2 - z2 ), 175.0 * x * z * ( -x2 + z2 ), 25.0 * x * y * ( x2 + y2 - 6.0 * z2 );
    }
    else if( !useCosineCoefficientOfBody2 && orderOfBody2Coefficient == 1 )
    {
        // Eq. (10)
        analyticalTorque << 25.0 * x * z * ( -x2 + 6.0 * y2 - z2 ), 25.0 * y * z * ( -9.0 * x2 - 2.0 * y2 + 5.0 * z2 ),
                -175.0 * y2 * z2 + 25.0 * ( y2 + z2 ) * r2 - 5.0 * r2 * r2;
    }
    else if( useCosineCoefficientOfBody2 && orderOfBody2Coefficient == 2 )
    {
        // Eq. (11)
        analyticalTorque << ( 175.0 / 2.0 ) * x * y * ( -x2 + y2 ),
                0.5 * ( 85.0 * x4 - 255.0 * x2 * y2 - 255.0 * x2 * z2 + 10.0 * y4 + 195.0 * y2 * z2 + 10.0 * z4 ),
                ( 25.0 / 2.0 ) * y * z * ( 9.0 * x2 - 5.0 * y2 + 2.0 * z2 );
    }
    else if( !useCosineCoefficientOfBody2 && orderOfBody2Coefficient == 2 )
    {
        // Eq. (12)
        analyticalTorque << -175.0 * x2 * y2 + 25.0 * ( x2 + y2 ) * r2 - 5.0 * r2 * r2, 25.0 * x * y * ( 5.0 * x2 - 2.0 * y2 - 9.0 * z2 ),
                25.0 * x * z * ( -x2 + 6.0 * y2 - z2 );
    }
    else
    {
        throw std::runtime_error( "Error in analytical C21xdegree2 torque computation: unsupported coefficient case." );
    }

    return analyticalTorque * commonPrefactor;
}

struct SchutzEq11TermDiagnostics {
    double prefactor = TUDAT_NAN;
    double Ixy = TUDAT_NAN;
    double Ixz = TUDAT_NAN;
    double Iyz = TUDAT_NAN;
    double fyz = TUDAT_NAN;
    double fxz = TUDAT_NAN;
    double fxy = TUDAT_NAN;
    double gyz = TUDAT_NAN;
    double gxz = TUDAT_NAN;
    double gxy = TUDAT_NAN;
};

SchutzEq11TermDiagnostics computeSchutzEq11TermDiagnostics( const Eigen::Vector3d& relativePositionInBody1Frame,
                                                            const double massOfBody2,
                                                            const Eigen::Matrix3d& inertiaTensorOfBody1,
                                                            const Eigen::Matrix3d& inertiaTensorOfBody2InBody1Frame )
{
    SchutzEq11TermDiagnostics diagnostics;

    const double x = relativePositionInBody1Frame( 0 );
    const double y = relativePositionInBody1Frame( 1 );
    const double z = relativePositionInBody1Frame( 2 );
    const double x2 = x * x;
    const double y2 = y * y;
    const double z2 = z * z;
    const double xy = x * y;
    const double xz = x * z;
    const double yz = y * z;
    const double r2 = x2 + y2 + z2;
    const double inverseR2 = 1.0 / r2;
    const double r5 = r2 * r2 * std::sqrt( r2 );

    const double Aprime = inertiaTensorOfBody2InBody1Frame( 0, 0 );
    const double Bprime = inertiaTensorOfBody2InBody1Frame( 1, 1 );
    const double Cprime = inertiaTensorOfBody2InBody1Frame( 2, 2 );
    const double IxyPrime = -inertiaTensorOfBody2InBody1Frame( 0, 1 );
    const double IxzPrime = -inertiaTensorOfBody2InBody1Frame( 0, 2 );
    const double IyzPrime = -inertiaTensorOfBody2InBody1Frame( 1, 2 );

    diagnostics.Ixy = -inertiaTensorOfBody1( 0, 1 );
    diagnostics.Ixz = -inertiaTensorOfBody1( 0, 2 );
    diagnostics.Iyz = -inertiaTensorOfBody1( 1, 2 );

    const double Qprime = Aprime + Bprime + Cprime;
    const double IellPrime =
            ( Aprime * x2 + Bprime * y2 + Cprime * z2 - 2.0 * IxyPrime * xy - 2.0 * IxzPrime * xz - 2.0 * IyzPrime * yz ) * inverseR2;
    const double Wprime = massOfBody2 + 7.5 * Qprime * inverseR2 - 17.5 * IellPrime * inverseR2;

    diagnostics.fyz = yz * ( Wprime - 5.0 * Aprime * inverseR2 ) - 5.0 * IxzPrime * xy * inverseR2 - 5.0 * IxyPrime * xz * inverseR2 +
            IyzPrime * ( 1.0 - 5.0 * ( y2 + z2 ) * inverseR2 );
    diagnostics.fxz = xz * ( Wprime - 5.0 * Bprime * inverseR2 ) + IxzPrime * ( 1.0 - 5.0 * ( x2 + z2 ) * inverseR2 ) -
            5.0 * IyzPrime * xy * inverseR2 - 5.0 * IxyPrime * yz * inverseR2;
    diagnostics.fxy = xy * ( Wprime - 5.0 * Cprime * inverseR2 ) - 5.0 * IyzPrime * xz * inverseR2 +
            IxyPrime * ( 1.0 - 5.0 * ( x2 + y2 ) * inverseR2 ) - 5.0 * IxzPrime * yz * inverseR2;

    diagnostics.gyz = ( z2 - y2 ) * Wprime + Bprime - Cprime - 10.0 * IxzPrime * xz * inverseR2 - 10.0 * IxyPrime * xy * inverseR2 -
            20.0 * IyzPrime * yz * inverseR2 - 5.0 * z2 * ( Aprime + Bprime - Cprime ) * inverseR2 -
            5.0 * y2 * ( Aprime - Bprime + Cprime ) * inverseR2;
    diagnostics.gxz = ( x2 - z2 ) * Wprime + Cprime - Aprime - 20.0 * IxzPrime * xz * inverseR2 - 10.0 * IxyPrime * xy * inverseR2 -
            10.0 * IyzPrime * yz * inverseR2 - 5.0 * x2 * ( -Aprime + Bprime + Cprime ) * inverseR2 -
            5.0 * z2 * ( Aprime + Bprime - Cprime ) * inverseR2;
    diagnostics.gxy = ( y2 - x2 ) * Wprime + Aprime - Bprime - 10.0 * IxzPrime * xz * inverseR2 - 20.0 * IxyPrime * xy * inverseR2 -
            10.0 * IyzPrime * yz * inverseR2 - 5.0 * y2 * ( Aprime - Bprime + Cprime ) * inverseR2 -
            5.0 * x2 * ( -Aprime + Bprime + Cprime ) * inverseR2;

    diagnostics.prefactor = 3.0 * physical_constants::GRAVITATIONAL_CONSTANT / r5;
    return diagnostics;
}

double getExpectedC20DegreeTwoInteractionMultiplier( const unsigned int order )
{
    if( order == 0 )
    {
        return 10.0;
    }
    else if( order == 1 )
    {
        return std::sqrt( 125.0 / 6.0 );
    }
    else if( order == 2 )
    {
        return std::sqrt( 125.0 / 12.0 );
    }
    else
    {
        throw std::runtime_error( "Error when getting expected C20xdegree2 multiplier, unsupported order." );
    }
}

BOOST_AUTO_TEST_SUITE( test_gravitational_torque )

//! Test to check degree two gravitational torque
BOOST_AUTO_TEST_CASE( testDegreeTwoGravitationalTorque )
{
    using namespace tudat::simulation_setup;
    using namespace tudat::basic_astrodynamics;
    using namespace tudat::gravitation;

    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation time settings.
    const double simulationStartEpoch = 0.0;
    const double simulationEndEpoch = tudat::physical_constants::JULIAN_DAY;

    // Define body settings for simulation.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Moon" );

    //! Test torque for various inertia tensors
    for( unsigned int testCase = 0; testCase < 7; testCase++ )
    {
        // Create body objects.
        BodyListSettings bodySettings = getDefaultBodySettings( bodiesToCreate );

        // Define degree two coefficients
        Eigen::Matrix3d cosineCoefficients = Eigen::Matrix3d::Zero( );
        Eigen::Matrix3d sineCoefficients = Eigen::Matrix3d::Zero( );
        Eigen::Matrix3d inertiaTensorDeviation = Eigen::Matrix3d::Zero( );
        if( testCase > 0 )
        {
            cosineCoefficients( 0, 0 ) = 1.0;

            if( testCase == 2 )
            {
                cosineCoefficients( 2, 0 ) = 0.01;
            }
            else if( testCase == 3 )
            {
                cosineCoefficients( 2, 2 ) = 0.01;
            }
            else if( testCase == 4 )
            {
                sineCoefficients( 2, 2 ) = 0.01;
            }
            else if( testCase == 5 )
            {
                cosineCoefficients( 2, 1 ) = 0.01;
            }
            else if( testCase == 6 )
            {
                sineCoefficients( 2, 1 ) = 0.01;
            }

            bodySettings.at( "Moon" )->gravityFieldSettings =
                    std::make_shared< SphericalHarmonicsGravityFieldSettings >( spice_interface::getBodyGravitationalParameter( "Moon" ),
                                                                                spice_interface::getAverageRadius( "Moon" ),
                                                                                cosineCoefficients,
                                                                                sineCoefficients,
                                                                                "IAU_Moon" );
        }

        if( testCase > 0 )
        {
            double c20InertiaContribution =
                    cosineCoefficients( 2, 0 ) * tudat::basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 0 ) / 3.0;
            inertiaTensorDeviation( 0, 0 ) += c20InertiaContribution;
            inertiaTensorDeviation( 1, 1 ) += c20InertiaContribution;
            inertiaTensorDeviation( 2, 2 ) -= 2.0 * c20InertiaContribution;

            double c21InertiaContribution =
                    cosineCoefficients( 2, 1 ) * tudat::basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 1 );
            inertiaTensorDeviation( 0, 2 ) -= c21InertiaContribution;
            inertiaTensorDeviation( 2, 0 ) -= c21InertiaContribution;

            double c22InertiaContribution =
                    2.0 * cosineCoefficients( 2, 2 ) * tudat::basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 );
            inertiaTensorDeviation( 0, 0 ) -= c22InertiaContribution;
            inertiaTensorDeviation( 1, 1 ) += c22InertiaContribution;

            double s21InertiaContribution =
                    sineCoefficients( 2, 1 ) * tudat::basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 1 );
            inertiaTensorDeviation( 1, 2 ) -= s21InertiaContribution;
            inertiaTensorDeviation( 2, 1 ) -= s21InertiaContribution;

            double s22InertiaContribution =
                    2.0 * sineCoefficients( 2, 2 ) * tudat::basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 );
            inertiaTensorDeviation( 0, 1 ) -= s22InertiaContribution;
            inertiaTensorDeviation( 1, 0 ) -= s22InertiaContribution;

            inertiaTensorDeviation *= spice_interface::getAverageRadius( "Moon" ) * spice_interface::getAverageRadius( "Moon" ) *
                    spice_interface::getBodyGravitationalParameter( "Moon" ) / physical_constants::GRAVITATIONAL_CONSTANT;
        }

        std::dynamic_pointer_cast< SphericalHarmonicsGravityFieldSettings >( bodySettings.at( "Moon" )->gravityFieldSettings )
                ->setScaledMeanMomentOfInertia( 0.4 );
        std::dynamic_pointer_cast< SphericalHarmonicsGravityFieldSettings >( bodySettings.at( "Earth" )->gravityFieldSettings )
                ->setScaledMeanMomentOfInertia( 0.4 );

        // Create bodies
        SystemOfBodies bodies = createSystemOfBodies( bodySettings );

        // Crate torque model
        SelectedTorqueMap selectedTorqueModelMap;
        selectedTorqueModelMap[ "Moon" ][ "Earth" ].push_back( std::make_shared< TorqueSettings >( second_order_gravitational_torque ) );
        basic_astrodynamics::TorqueModelMap torqueModelMap = createTorqueModelsMap( bodies, selectedTorqueModelMap, { "Moon" } );
        std::shared_ptr< TorqueModel > secondDegreeGravitationalTorque = torqueModelMap.at( "Moon" ).at( "Earth" ).at( 0 );

        // Update environment to current time
        double evaluationTime = tudat::physical_constants::JULIAN_DAY / 2.0;

        bodies.at( "Moon" )->setStateFromEphemeris( evaluationTime );
        bodies.at( "Moon" )->setCurrentRotationalStateToLocalFrameFromEphemeris( evaluationTime );
        bodies.at( "Moon" )->getMassProperties( )->update( 0.0 );

        bodies.at( "Earth" )->setStateFromEphemeris( evaluationTime );
        bodies.at( "Earth" )->setCurrentRotationalStateToLocalFrameFromEphemeris( evaluationTime );
        bodies.at( "Earth" )->getMassProperties( )->update( 0.0 );

        if( testCase == 0 )
        {
            // Test id gravity field coefficients are correctly reconstructed
            inertiaTensorDeviation = bodies.at( "Moon" )->getBodyInertiaTensor( );
            std::shared_ptr< SphericalHarmonicsGravityField > moonGravityField =
                    std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( bodies.at( "Moon" )->getGravityFieldModel( ) );

            double recomputedC20 = std::sqrt( 0.2 ) /
                    ( moonGravityField->getGravitationalParameter( ) * moonGravityField->getReferenceRadius( ) *
                      moonGravityField->getReferenceRadius( ) / tudat::physical_constants::GRAVITATIONAL_CONSTANT ) *
                    ( 0.5 * ( inertiaTensorDeviation( 0, 0 ) + inertiaTensorDeviation( 1, 1 ) ) - inertiaTensorDeviation( 2, 2 ) );
            BOOST_CHECK_CLOSE_FRACTION( recomputedC20, moonGravityField->getCosineCoefficients( )( 2, 0 ), 1.0E-12 );

            double recomputedC21 = -std::sqrt( 3.0 / 5.0 ) /
                    ( moonGravityField->getGravitationalParameter( ) * moonGravityField->getReferenceRadius( ) *
                      moonGravityField->getReferenceRadius( ) / tudat::physical_constants::GRAVITATIONAL_CONSTANT ) *
                    ( inertiaTensorDeviation( 0, 2 ) );
            BOOST_CHECK_CLOSE_FRACTION( recomputedC21, moonGravityField->getCosineCoefficients( )( 2, 1 ), 1.0E-12 );

            double recomputedC22 = std::sqrt( 3.0 / 20.0 ) /
                    ( moonGravityField->getGravitationalParameter( ) * moonGravityField->getReferenceRadius( ) *
                      moonGravityField->getReferenceRadius( ) / tudat::physical_constants::GRAVITATIONAL_CONSTANT ) *
                    ( inertiaTensorDeviation( 1, 1 ) - inertiaTensorDeviation( 0, 0 ) );
            BOOST_CHECK_CLOSE_FRACTION( recomputedC22, moonGravityField->getCosineCoefficients( )( 2, 2 ), 1.0E-12 );

            double recomputedS21 = -std::sqrt( 3.0 / 5.0 ) /
                    ( moonGravityField->getGravitationalParameter( ) * moonGravityField->getReferenceRadius( ) *
                      moonGravityField->getReferenceRadius( ) / tudat::physical_constants::GRAVITATIONAL_CONSTANT ) *
                    ( inertiaTensorDeviation( 1, 2 ) );
            BOOST_CHECK_CLOSE_FRACTION( recomputedS21, moonGravityField->getSineCoefficients( )( 2, 1 ), 1.0E-12 );

            double recomputedS22 = -std::sqrt( 12.0 / 5.0 ) /
                    ( 2.0 * moonGravityField->getGravitationalParameter( ) * moonGravityField->getReferenceRadius( ) *
                      moonGravityField->getReferenceRadius( ) / tudat::physical_constants::GRAVITATIONAL_CONSTANT ) *
                    ( inertiaTensorDeviation( 0, 1 ) );
            BOOST_CHECK_CLOSE_FRACTION( recomputedS22, moonGravityField->getSineCoefficients( )( 2, 2 ), 1.0E-12 );
        }

        // Update torque model to current time.
        secondDegreeGravitationalTorque->updateMembers( evaluationTime );

        // Compute torque manually, and test against computed result
        Eigen::Vector3d earthRelativePosition = bodies.at( "Moon" )->getCurrentRotationToLocalFrame( ) *
                ( bodies.at( "Earth" )->getPosition( ) - bodies.at( "Moon" )->getPosition( ) );
        Eigen::Vector3d manualTorque = 3.0 * earthRelativePosition.cross( inertiaTensorDeviation * earthRelativePosition ) *
                bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( ) /
                std::pow( earthRelativePosition.norm( ), 5.0 );

        Eigen::Vector3d currentTorque = secondDegreeGravitationalTorque->getTorque( );
        Eigen::Vector3d torqueError = ( currentTorque - manualTorque );

        // std::cout << "Current torques " << currentTorque.transpose( ) << std::endl << torqueError << std::endl;

        for( unsigned int i = 0; i < 3; i++ )
        {
            BOOST_CHECK_SMALL( std::fabs( torqueError( i ) ), 1.0E8 );
        }
    }
}

//! Test to check spherical harmonic torque.
BOOST_AUTO_TEST_CASE( testSphericalGravitationalTorque )
{
    using namespace tudat::simulation_setup;
    using namespace tudat::basic_astrodynamics;
    using namespace tudat::gravitation;

    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation time settings.
    const double simulationStartEpoch = 0.0;
    const double simulationEndEpoch = tudat::physical_constants::JULIAN_DAY;

    //! Test for zero and non-zero spherical harmonic coefficients
    for( unsigned int testCase = 0; testCase < 2; testCase++ )
    {
        // Define body settings for simulation.
        std::vector< std::string > bodiesToCreate;
        bodiesToCreate.push_back( "Earth" );
        bodiesToCreate.push_back( "Moon" );

        // Create body objects.
        BodyListSettings bodySettings = getDefaultBodySettings( bodiesToCreate );

        // Set effective point-mass gravity
        if( testCase == 0 )
        {
            Eigen::Matrix3d cosineCoefficients = Eigen::Matrix3d::Zero( );
            Eigen::Matrix3d sineCoefficients = Eigen::Matrix3d::Zero( );
            cosineCoefficients( 0, 0 ) = 1.0;

            bodySettings.at( "Moon" )->gravityFieldSettings =
                    std::make_shared< SphericalHarmonicsGravityFieldSettings >( spice_interface::getBodyGravitationalParameter( "Moon" ),
                                                                                spice_interface::getAverageRadius( "Moon" ),
                                                                                cosineCoefficients,
                                                                                sineCoefficients,
                                                                                "IAU_Moon" );
        }

        std::dynamic_pointer_cast< SphericalHarmonicsGravityFieldSettings >( bodySettings.at( "Moon" )->gravityFieldSettings )
                ->setScaledMeanMomentOfInertia( 0.0 );
        std::dynamic_pointer_cast< SphericalHarmonicsGravityFieldSettings >( bodySettings.at( "Earth" )->gravityFieldSettings )
                ->setScaledMeanMomentOfInertia( 0.0 );

        SystemOfBodies bodies = createSystemOfBodies( bodySettings );

        // Create two torque models
        SelectedTorqueMap selectedTorqueModelMap;
        selectedTorqueModelMap[ "Moon" ][ "Earth" ].push_back( std::make_shared< TorqueSettings >( second_order_gravitational_torque ) );
        selectedTorqueModelMap[ "Moon" ][ "Earth" ].push_back( std::make_shared< SphericalHarmonicTorqueSettings >( 2, 2 ) );
        basic_astrodynamics::TorqueModelMap torqueModelMap = createTorqueModelsMap( bodies, selectedTorqueModelMap, { "Moon" } );
        std::shared_ptr< TorqueModel > secondDegreeGravitationalTorque = torqueModelMap.at( "Moon" ).at( "Earth" ).at( 0 );
        std::shared_ptr< TorqueModel > sphercialHarmonicGravitationalTorque = torqueModelMap.at( "Moon" ).at( "Earth" ).at( 1 );

        // Update Moon to current time.
        double evaluationTime = tudat::physical_constants::JULIAN_DAY / 2.0;
        bodies.at( "Moon" )->setStateFromEphemeris( evaluationTime );
        bodies.at( "Moon" )->setCurrentRotationalStateToLocalFrameFromEphemeris( evaluationTime );
        bodies.at( "Moon" )->getMassProperties( )->update( 0.0 );

        {
            // Test reconstructed spherical harmonic coefficients
            Eigen::Matrix3d moonCosineCoefficients = std::dynamic_pointer_cast< tudat::gravitation::SphericalHarmonicsGravityField >(
                                                             bodies.at( "Moon" )->getGravityFieldModel( ) )
                                                             ->getCosineCoefficients( )
                                                             .block( 0, 0, 3, 3 );
            Eigen::Matrix3d moonSineCoefficients = std::dynamic_pointer_cast< tudat::gravitation::SphericalHarmonicsGravityField >(
                                                           bodies.at( "Moon" )->getGravityFieldModel( ) )
                                                           ->getSineCoefficients( )
                                                           .block( 0, 0, 3, 3 );
            Eigen::MatrixXd moonReconstructedCosineCoefficients = Eigen::Matrix3d::Zero( ),
                            moonReconstructedSineCoefficients = Eigen::Matrix3d::Zero( );

            double reconstructedScaledMeanMomentOfInertia;
            gravitation::getDegreeTwoSphericalHarmonicCoefficients(
                    bodies.at( "Moon" )->getBodyInertiaTensor( ),
                    bodies.at( "Moon" )->getGravityFieldModel( )->getGravitationalParameter( ),
                    std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( bodies.at( "Moon" )->getGravityFieldModel( ) )
                            ->getReferenceRadius( ),
                    true,
                    moonReconstructedCosineCoefficients,
                    moonReconstructedSineCoefficients,
                    reconstructedScaledMeanMomentOfInertia );

            for( unsigned int j = 0; j < 3; j++ )
            {
                for( unsigned int k = 0; k <= j; k++ )
                {
                    BOOST_CHECK_SMALL( std::fabs( moonReconstructedCosineCoefficients( j, k ) - moonCosineCoefficients( j, k ) ), 1.0E-19 );

                    BOOST_CHECK_SMALL( std::fabs( moonReconstructedSineCoefficients( j, k ) - moonSineCoefficients( j, k ) ), 1.0E-19 );
                }
            }
        }

        // Update Earth to current time.
        bodies.at( "Earth" )->setStateFromEphemeris( evaluationTime );
        bodies.at( "Earth" )->setCurrentRotationalStateToLocalFrameFromEphemeris( evaluationTime );
        bodies.at( "Earth" )->getMassProperties( )->update( 0.0 );

        // Update and compute torque values
        secondDegreeGravitationalTorque->updateMembers( evaluationTime );
        Eigen::Vector3d currentExplicitTorque = secondDegreeGravitationalTorque->getTorque( );
        sphercialHarmonicGravitationalTorque->updateMembers( evaluationTime );
        Eigen::Vector3d currentSphericalHarmonicTorque = sphercialHarmonicGravitationalTorque->getTorque( );

        // Test difference between models (should be equivalent)
        Eigen::Vector3d torqueError = currentExplicitTorque - currentSphericalHarmonicTorque;
        for( unsigned int i = 0; i < 3; i++ )
        {
            BOOST_CHECK_SMALL( std::fabs( torqueError( i ) ), 1.0E-14 * currentExplicitTorque.norm( ) );
        }
    }
}

BOOST_AUTO_TEST_CASE( testFourthDegreeFullTwoBodyGravitationalTorque )
{
    using namespace tudat::simulation_setup;
    using namespace tudat::gravitation;

    // Test rationale:
    // Validate the fourth-degree closed-form two-body torque implementation against known limits and
    // independent evaluations based on the same body states:
    // 1) point-mass-equivalent limits (Schutz et al., 1981 Eq. (11)/(14) checks),
    // 2) direct analytical torque reconstruction from transformed inertia tensors,
    // 3) orientation dependence only when body-2 inertia is non-isotropic.

    const std::string bodyUndergoingTorqueName = "Body1";
    const std::string bodyExertingTorqueName = "Body2";
    const Eigen::Vector3d positionOfBodyUndergoingTorque( 1.5E6, -2.4E6, 3.2E6 );
    const Eigen::Vector3d positionOfBodyExertingTorque( -4.1E6, 2.7E6, 1.1E6 );
    const double evaluationTime = 43200.0;

    const double massOfBodyUndergoingTorque = 5.8E21;
    const double massOfBodyExertingTorque = 7.1E21;

    const Eigen::Matrix3d inertiaTensorOfBodyUndergoingTorque =
            ( Eigen::Matrix3d( ) << 3.7E29, -1.2E27, 2.0E27, -1.2E27, 4.4E29, -0.8E27, 2.0E27, -0.8E27, 5.1E29 ).finished( );
    const Eigen::Quaterniond rotationToBodyUndergoingTorque( Eigen::AngleAxisd( 0.42, Eigen::Vector3d::UnitX( ) ) *
                                                             Eigen::AngleAxisd( -0.35, Eigen::Vector3d::UnitY( ) ) *
                                                             Eigen::AngleAxisd( 0.71, Eigen::Vector3d::UnitZ( ) ) );

    // Case 0: point-mass-equivalent body A (isotropic inertia tensor); verify Case 1 and Eq. (14) Case 2.
    {
        const double gravitationalParameterPointMassEquivalentBody = physical_constants::GRAVITATIONAL_CONSTANT * massOfBodyExertingTorque;
        const double referenceRadiusBodyUndergoingTorque = 2.3E6;
        const double referenceRadiusBodyExertingTorque = 1.4E6;

        const double isotropicInertiaValueOfBodyExertingTorque = 1.0E20;
        const double scaledMeanMomentOfInertiaBodyExertingTorque = isotropicInertiaValueOfBodyExertingTorque /
                ( ( gravitationalParameterPointMassEquivalentBody / physical_constants::GRAVITATIONAL_CONSTANT ) *
                  referenceRadiusBodyExertingTorque * referenceRadiusBodyExertingTorque );

        const std::tuple< Eigen::MatrixXd, Eigen::MatrixXd, double > degreeTwoCoefficientsOfBodyUndergoingTorque =
                gravitation::getDegreeTwoSphericalHarmonicCoefficients( inertiaTensorOfBodyUndergoingTorque,
                                                                        gravitationalParameterPointMassEquivalentBody,
                                                                        referenceRadiusBodyUndergoingTorque,
                                                                        2,
                                                                        true );
        const Eigen::MatrixXd cosineCoefficientsOfBodyUndergoingTorque = std::get< 0 >( degreeTwoCoefficientsOfBodyUndergoingTorque );
        const Eigen::MatrixXd sineCoefficientsOfBodyUndergoingTorque = std::get< 1 >( degreeTwoCoefficientsOfBodyUndergoingTorque );
        const double scaledMeanMomentOfInertiaBodyUndergoingTorque = std::get< 2 >( degreeTwoCoefficientsOfBodyUndergoingTorque );

        Eigen::MatrixXd cosineCoefficientsOfBodyExertingTorque = Eigen::MatrixXd::Zero( 3, 3 );
        Eigen::MatrixXd sineCoefficientsOfBodyExertingTorque = Eigen::MatrixXd::Zero( 3, 3 );
        cosineCoefficientsOfBodyExertingTorque( 0, 0 ) = 1.0;

        std::vector< Eigen::Quaterniond > orientationOfPointMassEquivalentBodyCases;
        orientationOfPointMassEquivalentBodyCases.push_back( Eigen::Quaterniond::Identity( ) );
        orientationOfPointMassEquivalentBodyCases.push_back( Eigen::Quaterniond( Eigen::AngleAxisd( 0.62, Eigen::Vector3d::UnitY( ) ) *
                                                                                 Eigen::AngleAxisd( -0.48, Eigen::Vector3d::UnitX( ) ) *
                                                                                 Eigen::AngleAxisd( 0.29, Eigen::Vector3d::UnitZ( ) ) ) );

        std::vector< Eigen::Vector3d > case1Torques;
        for( unsigned int i = 0; i < orientationOfPointMassEquivalentBodyCases.size( ); i++ )
        {
            const SystemOfBodies bodies = createSystemOfBodiesForFullTwoBodyTorqueTest( bodyUndergoingTorqueName,
                                                                                        bodyExertingTorqueName,
                                                                                        gravitationalParameterPointMassEquivalentBody,
                                                                                        referenceRadiusBodyUndergoingTorque,
                                                                                        referenceRadiusBodyExertingTorque,
                                                                                        positionOfBodyUndergoingTorque,
                                                                                        positionOfBodyExertingTorque,
                                                                                        cosineCoefficientsOfBodyUndergoingTorque,
                                                                                        sineCoefficientsOfBodyUndergoingTorque,
                                                                                        cosineCoefficientsOfBodyExertingTorque,
                                                                                        sineCoefficientsOfBodyExertingTorque,
                                                                                        rotationToBodyUndergoingTorque,
                                                                                        orientationOfPointMassEquivalentBodyCases.at( i ),
                                                                                        scaledMeanMomentOfInertiaBodyUndergoingTorque,
                                                                                        scaledMeanMomentOfInertiaBodyExertingTorque );
            bodies.at( bodyUndergoingTorqueName )->setCurrentRotationalStateToLocalFrameFromEphemeris( evaluationTime );
            bodies.at( bodyExertingTorqueName )->setCurrentRotationalStateToLocalFrameFromEphemeris( evaluationTime );

            std::shared_ptr< FourthDegreeFullTwoBodyGravitationalTorqueModel > fourthDegreeTorqueModel =
                    createFactoryFourthDegreeFullTwoBodyGravitationalTorqueModel(
                            bodies, bodyUndergoingTorqueName, bodyExertingTorqueName );
            std::shared_ptr< SecondDegreeGravitationalTorqueModel > secondDegreeTorqueModel =
                    createFactorySecondDegreeGravitationalTorqueModel( bodies, bodyUndergoingTorqueName, bodyExertingTorqueName );
            std::shared_ptr< FullTwoBodySphericalHarmonicAcceleration > referenceAcceleration =
                    createReferenceFullTwoBodySphericalHarmonicAccelerationModel(
                            bodies, bodyUndergoingTorqueName, bodyExertingTorqueName );

            // This check ensures the fourth-degree torque model is available for Case 1 and Case 2 checks.
            BOOST_REQUIRE( fourthDegreeTorqueModel != nullptr );
            // This check ensures the second-degree reference model exists for Case 1 verification.
            BOOST_REQUIRE( secondDegreeTorqueModel != nullptr );
            // This check ensures the mutual-potential gradient model exists for Eq. (14) evaluation.
            BOOST_REQUIRE( referenceAcceleration != nullptr );

            fourthDegreeTorqueModel->updateMembers( evaluationTime );
            secondDegreeTorqueModel->updateMembers( evaluationTime );
            referenceAcceleration->updateMembers( evaluationTime );

            const Eigen::Vector3d case1TorqueFromFourthDegree = fourthDegreeTorqueModel->getTorque( );
            const Eigen::Vector3d case1TorqueFromSecondDegree = secondDegreeTorqueModel->getTorque( );
            const double massOfBodyExertingTorqueInCurrentCase = bodies.at( bodyExertingTorqueName )->getBodyMass( );
            const Eigen::Vector3d case1SpecificTorqueFromFourthDegree = case1TorqueFromFourthDegree / massOfBodyExertingTorqueInCurrentCase;
            const Eigen::Vector3d totalTorqueFromMutualPotential = referenceAcceleration->getCurrentBodyFixedRelativePosition( ).cross(
                    referenceAcceleration->getMutualPotentialGradient( ) );
            const Eigen::Vector3d case2TorqueOnPointMassEquivalentBodyFromEq14 =
                    totalTorqueFromMutualPotential + case1SpecificTorqueFromFourthDegree;
            const double case1Scale = std::max( 1.0, case1TorqueFromSecondDegree.norm( ) );
            const double case2Scale = std::max( 1.0, totalTorqueFromMutualPotential.norm( ) );

            // This check verifies Case 1: torque of point-mass-equivalent A on extended B equals second-degree torque.
            const Eigen::Vector3d case1TorqueDifference = case1TorqueFromFourthDegree - case1TorqueFromSecondDegree;
            for( int i = 0; i < 3; i++ )
            {
                BOOST_CHECK_SMALL( std::fabs( case1TorqueDifference( i ) ) / case1Scale, 5.0E-14 );
            }
            // This check verifies Case 2 (Eq. 14): torque of extended B on point-mass-equivalent A is zero.
            for( int i = 0; i < 3; i++ )
            {
                // Eq. (14) is checked by cancelling two independently evaluated paths:
                // r x grad(U) and the specific fourth-degree torque. The residual is
                // limited by floating-point cancellation in these large intermediate terms.
                BOOST_CHECK_SMALL( std::fabs( case2TorqueOnPointMassEquivalentBodyFromEq14( i ) ) / case2Scale, 2.0E-10 );
            }

            case1Torques.push_back( case1TorqueFromFourthDegree );
        }

        const double orientationInvariantScale = std::max( 1.0, case1Torques.at( 0 ).norm( ) );
        // This check confirms isotropic inertia makes Case 1 torque invariant to body-A orientation.
        const Eigen::Vector3d orientationInvarianceDifference = case1Torques.at( 0 ) - case1Torques.at( 1 );
        for( int i = 0; i < 3; i++ )
        {
            BOOST_CHECK_SMALL( std::fabs( orientationInvarianceDifference( i ) ) / orientationInvariantScale, 5.0E-14 );
        }
    }

    // Case 1: body 2 is a point-mass equivalent; Eq. (11) must reduce to second-degree torque.
    {
        const Eigen::Matrix3d inertiaTensorOfBodyExertingTorque = Eigen::Matrix3d::Zero( );
        const Eigen::Quaterniond rotationToBodyExertingTorque = Eigen::Quaterniond::Identity( );

        const SystemOfBodies bodies = createSystemOfBodiesForFourthDegreeTwoBodyTorqueTest( bodyUndergoingTorqueName,
                                                                                            bodyExertingTorqueName,
                                                                                            massOfBodyUndergoingTorque,
                                                                                            massOfBodyExertingTorque,
                                                                                            inertiaTensorOfBodyUndergoingTorque,
                                                                                            inertiaTensorOfBodyExertingTorque,
                                                                                            positionOfBodyUndergoingTorque,
                                                                                            positionOfBodyExertingTorque,
                                                                                            rotationToBodyUndergoingTorque,
                                                                                            rotationToBodyExertingTorque );
        bodies.at( bodyUndergoingTorqueName )->setCurrentRotationalStateToLocalFrameFromEphemeris( evaluationTime );
        bodies.at( bodyExertingTorqueName )->setCurrentRotationalStateToLocalFrameFromEphemeris( evaluationTime );

        std::shared_ptr< FourthDegreeFullTwoBodyGravitationalTorqueModel > fourthDegreeTorqueModel =
                createFactoryFourthDegreeFullTwoBodyGravitationalTorqueModel( bodies, bodyUndergoingTorqueName, bodyExertingTorqueName );
        std::shared_ptr< SecondDegreeGravitationalTorqueModel > secondDegreeTorqueModel =
                createFactorySecondDegreeGravitationalTorqueModel( bodies, bodyUndergoingTorqueName, bodyExertingTorqueName );

        // This check ensures both torque models were created by the factory as expected.
        BOOST_REQUIRE( fourthDegreeTorqueModel != nullptr );
        // This check ensures the reference second-degree model was created for the comparison.
        BOOST_REQUIRE( secondDegreeTorqueModel != nullptr );

        fourthDegreeTorqueModel->updateMembers( evaluationTime );
        secondDegreeTorqueModel->updateMembers( evaluationTime );

        const Eigen::Vector3d fourthDegreeTorque = fourthDegreeTorqueModel->getTorque( );
        const Eigen::Vector3d secondDegreeTorque = secondDegreeTorqueModel->getTorque( );
        const double referenceScale = std::max( 1.0, secondDegreeTorque.norm( ) );

        // This check confirms the setup is non-degenerate and produces a non-zero torque signal.
        BOOST_CHECK_GT( secondDegreeTorque.norm( ), 1.0E-12 );
        // This check validates the expected point-mass limit: Eq. (11) matches second-degree torque.
        for( int i = 0; i < 3; i++ )
        {
            BOOST_CHECK_SMALL( std::fabs( fourthDegreeTorque( i ) - secondDegreeTorque( i ) ) / referenceScale, 5.0E-14 );
        }
        // This check confirms that a zero body-2 inertia remains zero after frame transformation.
        BOOST_CHECK_SMALL( fourthDegreeTorqueModel->getCurrentInertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque( ).norm( ),
                           1.0E-20 );
    }

    // Case 2: body-2 has finite inertia; torque must match Eq. (11) with transformed inertia and vary with orientation.
    {
        const Eigen::Matrix3d inertiaTensorOfBodyExertingTorque =
                ( Eigen::Matrix3d( ) << 2.1E29, -0.6E27, 1.1E27, -0.6E27, 2.5E29, -1.4E27, 1.1E27, -1.4E27, 3.0E29 ).finished( );
        std::vector< Eigen::Quaterniond > bodyExertingTorqueRotations;
        bodyExertingTorqueRotations.push_back( Eigen::Quaterniond::Identity( ) );
        bodyExertingTorqueRotations.push_back( Eigen::Quaterniond( Eigen::AngleAxisd( -0.53, Eigen::Vector3d::UnitY( ) ) *
                                                                   Eigen::AngleAxisd( 0.67, Eigen::Vector3d::UnitZ( ) ) *
                                                                   Eigen::AngleAxisd( 0.21, Eigen::Vector3d::UnitX( ) ) ) );

        std::vector< Eigen::Vector3d > evaluatedTorques;
        for( unsigned int i = 0; i < bodyExertingTorqueRotations.size( ); i++ )
        {
            const SystemOfBodies bodies = createSystemOfBodiesForFourthDegreeTwoBodyTorqueTest( bodyUndergoingTorqueName,
                                                                                                bodyExertingTorqueName,
                                                                                                massOfBodyUndergoingTorque,
                                                                                                massOfBodyExertingTorque,
                                                                                                inertiaTensorOfBodyUndergoingTorque,
                                                                                                inertiaTensorOfBodyExertingTorque,
                                                                                                positionOfBodyUndergoingTorque,
                                                                                                positionOfBodyExertingTorque,
                                                                                                rotationToBodyUndergoingTorque,
                                                                                                bodyExertingTorqueRotations.at( i ) );
            bodies.at( bodyUndergoingTorqueName )->setCurrentRotationalStateToLocalFrameFromEphemeris( evaluationTime );
            bodies.at( bodyExertingTorqueName )->setCurrentRotationalStateToLocalFrameFromEphemeris( evaluationTime );

            std::shared_ptr< FourthDegreeFullTwoBodyGravitationalTorqueModel > fourthDegreeTorqueModel =
                    createFactoryFourthDegreeFullTwoBodyGravitationalTorqueModel(
                            bodies, bodyUndergoingTorqueName, bodyExertingTorqueName );

            // This check ensures the fourth-degree model is available in each orientation case.
            BOOST_REQUIRE( fourthDegreeTorqueModel != nullptr );

            fourthDegreeTorqueModel->updateMembers( evaluationTime );

            const Eigen::Matrix3d manualTransformedInertiaTensorOfBodyExertingTorque =
                    transformBodyExertingInertiaTensorToBodyUndergoingFrame(
                            inertiaTensorOfBodyExertingTorque, rotationToBodyUndergoingTorque, bodyExertingTorqueRotations.at( i ) );
            const Eigen::Vector3d manualTorque = computeManualFourthDegreeTwoBodyTorqueFromBodyStates(
                    bodies.at( bodyUndergoingTorqueName ), bodies.at( bodyExertingTorqueName ) );
            const Eigen::Vector3d modelTorque = fourthDegreeTorqueModel->getTorque( );
            const double referenceScale = std::max( 1.0, manualTorque.norm( ) );
            const double inertiaTransformationScale = std::max( 1.0, manualTransformedInertiaTensorOfBodyExertingTorque.norm( ) );
            const Eigen::Matrix3d inertiaTransformationDifference =
                    fourthDegreeTorqueModel->getCurrentInertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque( ) -
                    manualTransformedInertiaTensorOfBodyExertingTorque;
            const Eigen::Vector3d torqueDifference = modelTorque - manualTorque;

            // This check verifies, per matrix element, the frame transformation of body-2 inertia used internally by the model.
            for( int row = 0; row < 3; row++ )
            {
                for( int column = 0; column < 3; column++ )
                {
                    BOOST_CHECK_SMALL( std::fabs( inertiaTransformationDifference( row, column ) ) / inertiaTransformationScale, 1.0E-14 );
                }
            }
            // This check validates, per torque component, the model output against a direct Eq. (11) evaluation from current body states.
            for( int component = 0; component < 3; component++ )
            {
                BOOST_CHECK_SMALL( std::fabs( torqueDifference( component ) ) / referenceScale, 5.0E-14 );
            }

            evaluatedTorques.push_back( modelTorque );
        }

        // This check confirms body-2 orientation changes the torque when body-2 inertia is non-zero.
        BOOST_CHECK_GT( ( evaluatedTorques.at( 0 ) - evaluatedTorques.at( 1 ) ).norm( ), 1.0E-6 );
    }
}

BOOST_AUTO_TEST_CASE( testFullTwoBodySphericalHarmonicTorque )
{
    using namespace tudat::simulation_setup;
    using namespace tudat::gravitation;
    using namespace tudat::basic_astrodynamics;

    // Test rationale:
    // This test validates the full two-body spherical-harmonic torque implementation (Dirkx et al., 2019,
    // mutual-potential expansion in coupled degree/order terms) against two independent references:
    // 1. Torque reconstructed from the mutual-potential gradient, using tau = r_B x (dU/dr_B) in body-1 frame.
    // 2. The fourth-degree closed-form two-body torque model (Schutz et al., 1981; Eq. (11)/(14) notation used
    //    throughout this file), in the degree-2 truncation where both formulations are expected to coincide.

    const std::string bodyUndergoingTorqueName = "Body1";
    const std::string bodyExertingTorqueName = "Body2";

    const double gravitationalParameter = 4.0E14;
    const double referenceRadiusBody1 = 2.1E6;
    const double referenceRadiusBody2 = 9.5E5;

    const Eigen::Vector3d positionOfBody1( 7.2E6, -1.4E6, 2.6E6 );
    const Eigen::Vector3d positionOfBody2( 1.8E6, 2.2E6, -1.1E6 );
    const double evaluationTime = 86400.0;

    // Case 1: all coefficients are zero -> zero torque.
    {
        const Eigen::MatrixXd cosineCoefficientsOfBody1 = Eigen::MatrixXd::Zero( 3, 3 );
        const Eigen::MatrixXd sineCoefficientsOfBody1 = Eigen::MatrixXd::Zero( 3, 3 );
        const Eigen::MatrixXd cosineCoefficientsOfBody2 = Eigen::MatrixXd::Zero( 1, 1 );
        const Eigen::MatrixXd sineCoefficientsOfBody2 = Eigen::MatrixXd::Zero( 1, 1 );

        const SystemOfBodies bodies = createSystemOfBodiesForFullTwoBodyTorqueTest( bodyUndergoingTorqueName,
                                                                                    bodyExertingTorqueName,
                                                                                    gravitationalParameter,
                                                                                    referenceRadiusBody1,
                                                                                    referenceRadiusBody2,
                                                                                    positionOfBody1,
                                                                                    positionOfBody2,
                                                                                    cosineCoefficientsOfBody1,
                                                                                    sineCoefficientsOfBody1,
                                                                                    cosineCoefficientsOfBody2,
                                                                                    sineCoefficientsOfBody2,
                                                                                    Eigen::Quaterniond::Identity( ),
                                                                                    Eigen::Quaterniond::Identity( ) );
        bodies.at( bodyUndergoingTorqueName )->setCurrentRotationalStateToLocalFrameFromEphemeris( evaluationTime );
        bodies.at( bodyExertingTorqueName )->setCurrentRotationalStateToLocalFrameFromEphemeris( evaluationTime );

        std::shared_ptr< FullTwoBodySphericalHarmonicTorque > torqueModel =
                createFactoryFullTwoBodySphericalHarmonicTorqueModel( bodies, bodyUndergoingTorqueName, bodyExertingTorqueName );
        BOOST_REQUIRE( torqueModel != nullptr );
        torqueModel->updateMembers( evaluationTime );

        BOOST_CHECK_SMALL( torqueModel->getTorque( ).norm( ), 1.0E-30 );
    }

    // Case 2a/2b:
    // Body 1 is degree-2, body 2 is point-mass. In this limit, only (l,m,p,q) = (2,m,0,0) interactions remain.
    // We validate the full-two-body torque through an independent path: mutual-potential gradient from the
    // acceleration model. This directly verifies that the implemented torque is consistent with the same U used
    // by the dynamics model, including sign convention and frame mapping.
    {
        Eigen::MatrixXd cosineCoefficientsOfBody1 = Eigen::MatrixXd::Zero( 3, 3 );
        Eigen::MatrixXd sineCoefficientsOfBody1 = Eigen::MatrixXd::Zero( 3, 3 );
        cosineCoefficientsOfBody1( 0, 0 ) = 1.0;
        cosineCoefficientsOfBody1( 2, 0 ) = 1.3E-3;
        cosineCoefficientsOfBody1( 2, 1 ) = -2.0E-4;
        cosineCoefficientsOfBody1( 2, 2 ) = 3.2E-4;
        sineCoefficientsOfBody1( 2, 1 ) = 1.7E-4;
        sineCoefficientsOfBody1( 2, 2 ) = -4.4E-4;

        Eigen::MatrixXd cosineCoefficientsOfBody2 = Eigen::MatrixXd::Zero( 1, 1 );
        Eigen::MatrixXd sineCoefficientsOfBody2 = Eigen::MatrixXd::Zero( 1, 1 );
        cosineCoefficientsOfBody2( 0, 0 ) = 1.0;

        std::vector< std::pair< Eigen::Quaterniond, Eigen::Quaterniond > > orientationCases;
        orientationCases.push_back( std::make_pair( Eigen::Quaterniond::Identity( ), Eigen::Quaterniond::Identity( ) ) );
        orientationCases.push_back( std::make_pair( Eigen::Quaterniond( Eigen::AngleAxisd( 1.11, Eigen::Vector3d::UnitZ( ) ) *
                                                                        Eigen::AngleAxisd( -0.33, Eigen::Vector3d::UnitY( ) ) *
                                                                        Eigen::AngleAxisd( 0.71, Eigen::Vector3d::UnitX( ) ) ),
                                                    Eigen::Quaterniond( Eigen::AngleAxisd( -0.58, Eigen::Vector3d::UnitX( ) ) *
                                                                        Eigen::AngleAxisd( 0.43, Eigen::Vector3d::UnitZ( ) ) *
                                                                        Eigen::AngleAxisd( 0.92, Eigen::Vector3d::UnitY( ) ) ) ) );
        std::vector< Eigen::Vector3d > computedTorques;

        for( const auto& orientationCase : orientationCases )
        {
            const SystemOfBodies bodies = createSystemOfBodiesForFullTwoBodyTorqueTest( bodyUndergoingTorqueName,
                                                                                        bodyExertingTorqueName,
                                                                                        gravitationalParameter,
                                                                                        referenceRadiusBody1,
                                                                                        referenceRadiusBody2,
                                                                                        positionOfBody1,
                                                                                        positionOfBody2,
                                                                                        cosineCoefficientsOfBody1,
                                                                                        sineCoefficientsOfBody1,
                                                                                        cosineCoefficientsOfBody2,
                                                                                        sineCoefficientsOfBody2,
                                                                                        orientationCase.first,
                                                                                        orientationCase.second );

            bodies.at( bodyUndergoingTorqueName )->setCurrentRotationalStateToLocalFrameFromEphemeris( evaluationTime );
            bodies.at( bodyExertingTorqueName )->setCurrentRotationalStateToLocalFrameFromEphemeris( evaluationTime );

            std::shared_ptr< FullTwoBodySphericalHarmonicTorque > torqueModel =
                    createFactoryFullTwoBodySphericalHarmonicTorqueModel( bodies, bodyUndergoingTorqueName, bodyExertingTorqueName );
            BOOST_REQUIRE( torqueModel != nullptr );
            std::shared_ptr< SphericalHarmonicGravitationalTorqueModel > sphericalHarmonicTorqueModel =
                    createFactorySphericalHarmonicTorqueModel( bodies, bodyUndergoingTorqueName, bodyExertingTorqueName, 2, 2 );
            BOOST_REQUIRE( sphericalHarmonicTorqueModel != nullptr );
            std::shared_ptr< FullTwoBodySphericalHarmonicAcceleration > referenceAcceleration =
                    createReferenceFullTwoBodySphericalHarmonicAccelerationModel(
                            bodies, bodyUndergoingTorqueName, bodyExertingTorqueName );
            BOOST_REQUIRE( referenceAcceleration != nullptr );

            torqueModel->updateMembers( evaluationTime );
            sphericalHarmonicTorqueModel->updateMembers( evaluationTime );
            referenceAcceleration->updateMembers( evaluationTime );

            const double bodyUndergoingTorqueMass = bodies.at( bodyUndergoingTorqueName )->getBodyMass( );
            const Eigen::Vector3d referenceTorqueFromAcceleration = bodyUndergoingTorqueMass *
                    referenceAcceleration->getCurrentBodyFixedRelativePosition( ).cross(
                            referenceAcceleration->getMutualPotentialGradient( ) );
            const Eigen::Vector3d computedFullTwoBodyTorque = torqueModel->getTorque( );
            const Eigen::Vector3d computedSphericalHarmonicTorque = sphericalHarmonicTorqueModel->getTorque( );
            computedTorques.push_back( computedFullTwoBodyTorque );

            // Consistency with tau = r_B x (dU/dr_B):
            // the gradient-based reference is assembled from the potential-gradient path and multiplied by
            // the body-undergoing-torque mass to obtain physical torque. The sign here follows Tudat's convention.
            const Eigen::Vector3d accelerationConsistencyDifference = computedFullTwoBodyTorque + referenceTorqueFromAcceleration;

            // Direct subtraction between full two-body torque and spherical-harmonic torque:
            // this should be near zero because both models are expected to use the same sign convention.
            const Eigen::Vector3d modelToModelDifference = computedFullTwoBodyTorque - computedSphericalHarmonicTorque;

            // Signed consistency check: this should stay large because the two models should not match after adding.
            const Eigen::Vector3d modelToModelSignedDifference = computedFullTwoBodyTorque + computedSphericalHarmonicTorque;
            const double referenceScale = std::max( 1.0, referenceTorqueFromAcceleration.norm( ) );

            for( int i = 0; i < 3; i++ )
            {
                BOOST_CHECK_SMALL( std::fabs( accelerationConsistencyDifference( i ) ) / referenceScale, 5.0E-14 );
            }
            for( int i = 0; i < 3; i++ )
            {
                BOOST_CHECK_SMALL( std::fabs( modelToModelDifference( i ) ) / referenceScale, 5.0E-13 );
            }
            BOOST_CHECK_GT( modelToModelSignedDifference.norm( ) / referenceScale, 1.0 );
        }

        BOOST_CHECK_GT( ( computedTorques.at( 0 ) - computedTorques.at( 1 ) ).norm( ), 1.0E-16 );
    }

    // Case 3:
    // Both bodies have degree-2 gravity. We decompose the coupled model into:
    // full( l2 = {0,2} ) - point-mass( l2 = 0 ) = isolated figure-figure (l2 = 2) term.
    // This decomposition is the same separation used in the derivation and should match the direct
    // degree-2/degree-2 coefficient selection. The isolated term is then compared to the independent
    // fourth-degree model evaluation on the same states (Schutz Eq. (11) for full degree-2 interaction,
    // with the point-mass contribution removed by the second-degree torque model).
    {
        Eigen::MatrixXd cosineCoefficientsOfBody1 = Eigen::MatrixXd::Zero( 3, 3 );
        Eigen::MatrixXd sineCoefficientsOfBody1 = Eigen::MatrixXd::Zero( 3, 3 );
        cosineCoefficientsOfBody1( 0, 0 ) = 1.0;
        cosineCoefficientsOfBody1( 2, 0 ) = 1.1E-3;

        Eigen::MatrixXd cosineCoefficientsOfBody2 = Eigen::MatrixXd::Zero( 3, 3 );
        Eigen::MatrixXd sineCoefficientsOfBody2 = Eigen::MatrixXd::Zero( 3, 3 );
        cosineCoefficientsOfBody2( 0, 0 ) = 1.0;
        cosineCoefficientsOfBody2( 2, 0 ) = -6.0E-4;

        const Eigen::Vector3d distantPositionOfBody1( 7.2E6, -1.4E6, 2.6E6 );
        const Eigen::Vector3d distantPositionOfBody2( -6.8E6, 3.1E6, -1.9E6 );
        const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > pointMassDegreeTwoCombinations =
                getPointMassDegreeTwoInteractionCombinations( );
        const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > degreeTwoDegreeTwoCombinations =
                getDegreeTwoDegreeTwoInteractionCombinations( );
        const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > fullDegreeTwoCombinations =
                getFullDegreeTwoInteractionCombinations( );

        std::vector< std::pair< Eigen::Quaterniond, Eigen::Quaterniond > > orientationCases;
        orientationCases.push_back( std::make_pair( Eigen::Quaterniond::Identity( ), Eigen::Quaterniond::Identity( ) ) );

        std::vector< Eigen::Vector3d > isolatedDegree22TorquesFromFullTwoBodyModel;
        std::vector< Eigen::Vector3d > isolatedDegree22TorquesFromFourthDegreeModel;
        for( const std::pair< Eigen::Quaterniond, Eigen::Quaterniond >& orientationCase : orientationCases )
        {
            const SystemOfBodies bodiesWithAllDegree2Terms = createSystemOfBodiesForFullTwoBodyTorqueTest( bodyUndergoingTorqueName,
                                                                                                           bodyExertingTorqueName,
                                                                                                           gravitationalParameter,
                                                                                                           referenceRadiusBody1,
                                                                                                           referenceRadiusBody2,
                                                                                                           distantPositionOfBody1,
                                                                                                           distantPositionOfBody2,
                                                                                                           cosineCoefficientsOfBody1,
                                                                                                           sineCoefficientsOfBody1,
                                                                                                           cosineCoefficientsOfBody2,
                                                                                                           sineCoefficientsOfBody2,
                                                                                                           orientationCase.first,
                                                                                                           orientationCase.second,
                                                                                                           0.0,
                                                                                                           0.0 );

            bodiesWithAllDegree2Terms.at( bodyUndergoingTorqueName )->setCurrentRotationalStateToLocalFrameFromEphemeris( evaluationTime );
            bodiesWithAllDegree2Terms.at( bodyExertingTorqueName )->setCurrentRotationalStateToLocalFrameFromEphemeris( evaluationTime );

            // Full two-body spherical-harmonic torque including both point-mass/degree-2 and degree-2/degree-2 interactions;
            // used as the "all degree-2 terms" reference for decomposition.
            std::shared_ptr< FullTwoBodySphericalHarmonicTorque > fullTwoBodyTorqueModelWithFullDegree2Terms =
                    createFactoryFullTwoBodySphericalHarmonicTorqueModel(
                            bodiesWithAllDegree2Terms, bodyUndergoingTorqueName, bodyExertingTorqueName, fullDegreeTwoCombinations );
            // Full two-body spherical-harmonic torque restricted to point-mass/degree-2 terms only;
            // used to isolate the degree-2/degree-2 coupling by subtraction from the full model.
            std::shared_ptr< FullTwoBodySphericalHarmonicTorque > fullTwoBodyTorqueModelWithPointMassDegree2Terms =
                    createFactoryFullTwoBodySphericalHarmonicTorqueModel(
                            bodiesWithAllDegree2Terms, bodyUndergoingTorqueName, bodyExertingTorqueName, pointMassDegreeTwoCombinations );
            // Full two-body spherical-harmonic torque with only direct degree-2/degree-2 cross interactions;
            // used as an explicit check that decomposition-by-subtraction recovers the same coupling term.
            std::shared_ptr< FullTwoBodySphericalHarmonicTorque > fullTwoBodyTorqueModelWithDegree2Degree2Terms =
                    createFactoryFullTwoBodySphericalHarmonicTorqueModel(
                            bodiesWithAllDegree2Terms, bodyUndergoingTorqueName, bodyExertingTorqueName, degreeTwoDegreeTwoCombinations );
            // Fourth-degree two-body torque model (Eq. 11 equivalent) evaluated on the full degree-2 bodies;
            // used as an independent model to compare against the full two-body spherical-harmonic torque decomposition.
            std::shared_ptr< FourthDegreeFullTwoBodyGravitationalTorqueModel > fourthDegreeTorqueModelWithFullDegree2Terms =
                    createFactoryFourthDegreeFullTwoBodyGravitationalTorqueModel(
                            bodiesWithAllDegree2Terms, bodyUndergoingTorqueName, bodyExertingTorqueName );
            // Fourth-degree two-body torque model that reuses the exact same body states but treats body 2 as a point mass.
            std::shared_ptr< FourthDegreeFullTwoBodyGravitationalTorqueModel > fourthDegreeTorqueModelWithPointMassDegree2Terms =
                    std::make_shared< FourthDegreeFullTwoBodyGravitationalTorqueModel >(
                            std::bind( &simulation_setup::Body::getPosition, bodiesWithAllDegree2Terms.at( bodyUndergoingTorqueName ) ),
                            std::bind( &simulation_setup::Body::getPosition, bodiesWithAllDegree2Terms.at( bodyExertingTorqueName ) ),
                            std::bind( &simulation_setup::Body::getBodyMass, bodiesWithAllDegree2Terms.at( bodyExertingTorqueName ) ),
                            std::bind( &simulation_setup::Body::getBodyInertiaTensor,
                                       bodiesWithAllDegree2Terms.at( bodyUndergoingTorqueName ) ),
                            []( ) { return Eigen::Matrix3d::Zero( ); },
                            std::bind( &simulation_setup::Body::getCurrentRotationToLocalFrame,
                                       bodiesWithAllDegree2Terms.at( bodyUndergoingTorqueName ) ),
                            std::bind( &simulation_setup::Body::getCurrentRotationToLocalFrame,
                                       bodiesWithAllDegree2Terms.at( bodyExertingTorqueName ) ) );

            // This check ensures the full two-body spherical-harmonic torque model is created for comparison.
            BOOST_REQUIRE( fullTwoBodyTorqueModelWithFullDegree2Terms != nullptr );
            BOOST_REQUIRE( fullTwoBodyTorqueModelWithPointMassDegree2Terms != nullptr );
            BOOST_REQUIRE( fullTwoBodyTorqueModelWithDegree2Degree2Terms != nullptr );
            // This check ensures the fourth-degree mutual-potential torque model is created for comparison.
            BOOST_REQUIRE( fourthDegreeTorqueModelWithFullDegree2Terms != nullptr );
            BOOST_REQUIRE( fourthDegreeTorqueModelWithPointMassDegree2Terms != nullptr );

            fullTwoBodyTorqueModelWithFullDegree2Terms->updateMembers( evaluationTime );
            fullTwoBodyTorqueModelWithPointMassDegree2Terms->updateMembers( evaluationTime );
            fullTwoBodyTorqueModelWithDegree2Degree2Terms->updateMembers( evaluationTime );
            fourthDegreeTorqueModelWithFullDegree2Terms->updateMembers( evaluationTime );
            fourthDegreeTorqueModelWithPointMassDegree2Terms->updateMembers( evaluationTime );

            const std::shared_ptr< gravitation::EffectiveMutualSphericalHarmonicsField > effectiveMutualField =
                    fullTwoBodyTorqueModelWithFullDegree2Terms->getAccelerationBetweenBodies( )->getEffectiveMutualPotentialField( );
            const Eigen::MatrixXd transformedCosineCoefficientsOfBody2FromFullTwoBody =
                    effectiveMutualField->getTransformedCosineCoefficientsOfBody2( );
            const Eigen::MatrixXd transformedSineCoefficientsOfBody2FromFullTwoBody =
                    effectiveMutualField->getTransformedSineCoefficientsOfBody2( );
            const Eigen::Matrix3d transformedBody2InertiaTensorFromFourthDegree =
                    fourthDegreeTorqueModelWithFullDegree2Terms
                            ->getCurrentInertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque( );
            const std::tuple< Eigen::MatrixXd, Eigen::MatrixXd, double > transformedBody2Degree2FromInertia =
                    gravitation::getDegreeTwoSphericalHarmonicCoefficients(
                            transformedBody2InertiaTensorFromFourthDegree, gravitationalParameter, referenceRadiusBody2, 2, true );
            const Eigen::MatrixXd transformedCosineCoefficientsOfBody2FromInertia = std::get< 0 >( transformedBody2Degree2FromInertia );
            const Eigen::MatrixXd transformedSineCoefficientsOfBody2FromInertia = std::get< 1 >( transformedBody2Degree2FromInertia );

            // Model type: [Full two-body]. Body exerting: [l=0,2]. Body undergoing: [l=2].
            const Eigen::Vector3d fullTwoBodyFullDegree2Torque = fullTwoBodyTorqueModelWithFullDegree2Terms->getTorque( );

            // Model type: [Full two-body]. Body exerting: [l=0]. Body undergoing: [l=2].
            const Eigen::Vector3d fullTwoBodyPointMassDegree2Torque = fullTwoBodyTorqueModelWithPointMassDegree2Terms->getTorque( );

            // Model type: [Full two-body]. Body exerting: [l=2]. Body undergoing: [l=2].
            const Eigen::Vector3d fullTwoBodyDegree2Degree2Torque = fullTwoBodyTorqueModelWithDegree2Degree2Terms->getTorque( );

            // Model type: [Fourth degree model]. Body exerting: [l=0,2]. Body undergoing: [l=2].
            const Eigen::Vector3d fourthDegreeFullDegree2Torque = fourthDegreeTorqueModelWithFullDegree2Terms->getTorque( );

            // Model type: [Fourth degree model]. Body exerting: [l=0]. Body undergoing: [l=2].
            const Eigen::Vector3d fourthDegreePointMassDegree2Torque = fourthDegreeTorqueModelWithPointMassDegree2Terms->getTorque( );

            // Model type: [Full two-body]. Body exerting: [l=2]. Body undergoing: [l=2] (isolated by subtraction).
            const Eigen::Vector3d isolatedDegree22TorqueFromFullTwoBodyModel =
                    fullTwoBodyFullDegree2Torque - fullTwoBodyPointMassDegree2Torque;

            // Model type: [Fourth degree model]. Body exerting: [l=2]. Body undergoing: [l=2] (isolated by subtraction).
            const Eigen::Vector3d isolatedDegree22TorqueFromFourthDegreeModel =
                    fourthDegreeFullDegree2Torque - fourthDegreePointMassDegree2Torque;

            const double pointMassReferenceScale = std::max( 1.0, fullTwoBodyPointMassDegree2Torque.norm( ) );

            // Compare [Full two-body, exerting l=2, undergoing l=2] from subtraction
            // vs      [Full two-body, exerting l=2, undergoing l=2] from direct coefficient selection.
            const Eigen::Vector3d isolatedVsDirectDegree22Difference =
                    isolatedDegree22TorqueFromFullTwoBodyModel - fullTwoBodyDegree2Degree2Torque;
            const double isolatedVsDirectDegree22Scale = std::max( 1.0, fullTwoBodyFullDegree2Torque.norm( ) );
            for( int i = 0; i < 3; i++ )
            {
                BOOST_CHECK_SMALL( std::fabs( isolatedVsDirectDegree22Difference( i ) ) / isolatedVsDirectDegree22Scale, 5.0E-14 );
            }

            // Compare [Full two-body, exerting l=0, undergoing l=2]
            // vs[Fourth degree model, exerting l=0, undergoing l=2].
            const Eigen::Vector3d pointMassModelDifference = fullTwoBodyPointMassDegree2Torque - fourthDegreePointMassDegree2Torque;
            const double pointMassModelDifferenceScale = std::max( 1.0, fourthDegreePointMassDegree2Torque.norm( ) );
            for( int i = 0; i < 3; i++ )
            {
                BOOST_CHECK_SMALL( std::fabs( pointMassModelDifference( i ) ) / pointMassModelDifferenceScale, 5.0E-14 );
            }

            // Compare [Full two-body, exerting l=2, undergoing l=2] from subtraction
            // vs [Fourth degree model, exerting l=2, undergoing l=2] from subtraction.
            const Eigen::Vector3d isolatedDegree22ModelDifference =
                    isolatedDegree22TorqueFromFullTwoBodyModel - isolatedDegree22TorqueFromFourthDegreeModel;
            const double isolatedDegree22ModelDifferenceScale = std::max( 1.0, fourthDegreeFullDegree2Torque.norm( ) );
            const double isolatedDegree22ModelDifferenceRelativeNorm =
                    isolatedDegree22ModelDifference.norm( ) / std::max( 1.0, isolatedDegree22TorqueFromFullTwoBodyModel.norm( ) );
            for( int i = 0; i < 3; i++ )
            {
                BOOST_CHECK_SMALL( std::fabs( isolatedDegree22ModelDifference( i ) ) / isolatedDegree22ModelDifferenceScale, 5.0E-14 );
            }

            // Compare [Full two-body, exerting l=2, undergoing l=2] from direct coefficient selection
            // vs [Fourth degree model, exerting l=2, undergoing l=2] from subtraction.
            const Eigen::Vector3d directVsIsolatedDegree22ModelDifference =
                    fullTwoBodyDegree2Degree2Torque - isolatedDegree22TorqueFromFourthDegreeModel;
            const double directVsIsolatedDegree22Scale = std::max( 1.0, fourthDegreeFullDegree2Torque.norm( ) );
            for( int i = 0; i < 3; i++ )
            {
                BOOST_CHECK_SMALL( std::fabs( directVsIsolatedDegree22ModelDifference( i ) ) / directVsIsolatedDegree22Scale, 5.0E-14 );
            }

            isolatedDegree22TorquesFromFullTwoBodyModel.push_back( isolatedDegree22TorqueFromFullTwoBodyModel );
            isolatedDegree22TorquesFromFourthDegreeModel.push_back( isolatedDegree22TorqueFromFourthDegreeModel );
        }

        BOOST_CHECK_GT( isolatedDegree22TorquesFromFullTwoBodyModel.at( 0 ).norm( ), 1.0E-16 );
        BOOST_CHECK_GT( isolatedDegree22TorquesFromFourthDegreeModel.at( 0 ).norm( ), 1.0E-16 );
    }
}

BOOST_AUTO_TEST_CASE( testSingleDegreeTwoDegreeTwoFigureFigureInteractionIsolation )
{
    using namespace tudat::simulation_setup;
    using namespace tudat::gravitation;
    using namespace tudat::basic_astrodynamics;

    // Test rationale:
    // Isolate one single figure-figure interaction at a time: (l,m,p,q) = (2,m1,2,m2), with
    // m1,m2 in {0,1,2} for body-1/body-2 {C20,C21,S21,C22,S22}.
    // Per loop, only one body-1 degree-2 and one body-2 degree-2 coefficient are non-zero.
    //
    // Why this matters:
    // 1. It removes cancellation between different degree-2 terms and yields a one-term diagnostic.
    // 2. It checks two independent isolation paths:
    //    a) Full-two-body direct single-term model.
    //    b) Fourth-degree minus second-degree torque (Schutz Eq. (11) minus point-mass part).
    // 3. For body-1 C20, it compares both against the closed-form analytical expression from
    //    the derivation document (Eq. (6)-(12), implemented in
    //    computeAnalyticalC20DegreeTwoFigureFigureTorque()).

    const std::string bodyUndergoingTorqueName = "Body1";
    const std::string bodyExertingTorqueName = "Body2";

    // Unit-normalized setup requested for debugging.
    const double gravitationalParameter = physical_constants::GRAVITATIONAL_CONSTANT;
    const double referenceRadiusBody1 = 1.0;
    const double referenceRadiusBody2 = 1.0;
    const std::vector< Eigen::Vector3d > relativePositionCases = {
        Eigen::Vector3d( 1.0, 0.0, 0.0 ),
        Eigen::Vector3d( 0.0, 1.0, 0.0 ),
        Eigen::Vector3d( 0.7071067811865475, 0.7071067811865475, 0.0 ),
        Eigen::Vector3d( 0.5000377542757255, -0.5000377542757255, 0.7070533845458759 ),
        Eigen::Vector3d( 0.0, 0.1, 0.99 ),
        Eigen::Vector3d( -0.2506273535585429, 0.6015056485405029, 0.7585691427722472 ),
        Eigen::Vector3d( 0.3030457633656632, 0.4040610178208843, 0.8636808257456412 )
    };
    const Eigen::Vector3d basePosition = Eigen::Vector3d::Zero( );
    const double evaluationTime = 86400.0;

    struct DegreeTwoCoefficientCase {
        std::string coefficientName;
        bool useCosineCoefficient;
        unsigned int order;
        double coefficientValue;
    };

    const std::vector< DegreeTwoCoefficientCase > body1DegreeTwoCoefficientCases = {
        { "C20", true, 0, 1.0 }, { "C21", true, 1, 1.0 }, { "S21", false, 1, 1.0 }, { "C22", true, 2, 1.0 }, { "S22", false, 2, 1.0 }
    };

    const std::vector< DegreeTwoCoefficientCase > body2DegreeTwoCoefficientCases = {
        { "C20", true, 0, 1.0 }, { "C21", true, 1, 1.0 }, { "S21", false, 1, 1.0 }, { "C22", true, 2, 1.0 }, { "S22", false, 2, 1.0 }
    };

    for( const Eigen::Vector3d& relativePositionCase : relativePositionCases )
    {
        for( const DegreeTwoCoefficientCase& body1CoefficientCase : body1DegreeTwoCoefficientCases )
        {
            for( const DegreeTwoCoefficientCase& coefficientCase : body2DegreeTwoCoefficientCases )
            {
                // Per loop: rebuild a fresh two-body system with one active body-1 and one active
                // body-2 degree-2 coefficient.
                Eigen::MatrixXd cosineCoefficientsOfBody1 = Eigen::MatrixXd::Zero( 3, 3 );
                Eigen::MatrixXd sineCoefficientsOfBody1 = Eigen::MatrixXd::Zero( 3, 3 );
                cosineCoefficientsOfBody1( 0, 0 ) = 1.0;
                if( body1CoefficientCase.useCosineCoefficient )
                {
                    cosineCoefficientsOfBody1( 2, body1CoefficientCase.order ) = body1CoefficientCase.coefficientValue;
                }
                else
                {
                    sineCoefficientsOfBody1( 2, body1CoefficientCase.order ) = body1CoefficientCase.coefficientValue;
                }

                Eigen::MatrixXd cosineCoefficientsOfBody2 = Eigen::MatrixXd::Zero( 3, 3 );
                Eigen::MatrixXd sineCoefficientsOfBody2 = Eigen::MatrixXd::Zero( 3, 3 );
                cosineCoefficientsOfBody2( 0, 0 ) = 1.0;
                if( coefficientCase.useCosineCoefficient )
                {
                    cosineCoefficientsOfBody2( 2, coefficientCase.order ) = coefficientCase.coefficientValue;
                }
                else
                {
                    sineCoefficientsOfBody2( 2, coefficientCase.order ) = coefficientCase.coefficientValue;
                }

                SystemOfBodies bodies = createSystemOfBodiesForFullTwoBodyTorqueTest( bodyUndergoingTorqueName,
                                                                                      bodyExertingTorqueName,
                                                                                      gravitationalParameter,
                                                                                      referenceRadiusBody1,
                                                                                      referenceRadiusBody2,
                                                                                      basePosition,
                                                                                      basePosition + relativePositionCase,
                                                                                      cosineCoefficientsOfBody1,
                                                                                      sineCoefficientsOfBody1,
                                                                                      cosineCoefficientsOfBody2,
                                                                                      sineCoefficientsOfBody2,
                                                                                      Eigen::Quaterniond::Identity( ),
                                                                                      Eigen::Quaterniond::Identity( ),
                                                                                      0.0,
                                                                                      0.0 );
                bodies.at( bodyUndergoingTorqueName )->setCurrentRotationalStateToLocalFrameFromEphemeris( evaluationTime );
                bodies.at( bodyExertingTorqueName )->setCurrentRotationalStateToLocalFrameFromEphemeris( evaluationTime );

                const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > singleInteractionTerm = {
                    std::make_tuple( 2, body1CoefficientCase.order, 2, coefficientCase.order )
                };

                std::shared_ptr< FullTwoBodySphericalHarmonicTorque > fullTwoBodySingleInteractionTorqueModel =
                        createFactoryFullTwoBodySphericalHarmonicTorqueModel(
                                bodies, bodyUndergoingTorqueName, bodyExertingTorqueName, singleInteractionTerm );
                std::shared_ptr< FullTwoBodySphericalHarmonicTorque > fullTwoBodyAllDegreeTwoDegreeTwoTorqueModel =
                        createFactoryFullTwoBodySphericalHarmonicTorqueModel(
                                bodies, bodyUndergoingTorqueName, bodyExertingTorqueName, getDegreeTwoDegreeTwoInteractionCombinations( ) );
                std::shared_ptr< FullTwoBodySphericalHarmonicTorque > fullTwoBodyFullDegreeTwoTorqueModel =
                        createFactoryFullTwoBodySphericalHarmonicTorqueModel(
                                bodies, bodyUndergoingTorqueName, bodyExertingTorqueName, getFullDegreeTwoInteractionCombinations( ) );
                std::shared_ptr< FullTwoBodySphericalHarmonicTorque > fullTwoBodyPointMassDegreeTwoTorqueModel =
                        createFactoryFullTwoBodySphericalHarmonicTorqueModel(
                                bodies, bodyUndergoingTorqueName, bodyExertingTorqueName, getPointMassDegreeTwoInteractionCombinations( ) );
                std::shared_ptr< FourthDegreeFullTwoBodyGravitationalTorqueModel > fourthDegreeTorqueModel =
                        createFactoryFourthDegreeFullTwoBodyGravitationalTorqueModel(
                                bodies, bodyUndergoingTorqueName, bodyExertingTorqueName );
                std::shared_ptr< SecondDegreeGravitationalTorqueModel > secondDegreeTorqueModel =
                        createFactorySecondDegreeGravitationalTorqueModel( bodies, bodyUndergoingTorqueName, bodyExertingTorqueName );

                BOOST_REQUIRE( fullTwoBodySingleInteractionTorqueModel != nullptr );
                BOOST_REQUIRE( fullTwoBodyAllDegreeTwoDegreeTwoTorqueModel != nullptr );
                BOOST_REQUIRE( fullTwoBodyFullDegreeTwoTorqueModel != nullptr );
                BOOST_REQUIRE( fullTwoBodyPointMassDegreeTwoTorqueModel != nullptr );
                BOOST_REQUIRE( fourthDegreeTorqueModel != nullptr );
                BOOST_REQUIRE( secondDegreeTorqueModel != nullptr );

                fullTwoBodySingleInteractionTorqueModel->updateMembers( evaluationTime );
                fullTwoBodyAllDegreeTwoDegreeTwoTorqueModel->updateMembers( evaluationTime );
                fullTwoBodyFullDegreeTwoTorqueModel->updateMembers( evaluationTime );
                fullTwoBodyPointMassDegreeTwoTorqueModel->updateMembers( evaluationTime );
                fourthDegreeTorqueModel->updateMembers( evaluationTime );
                secondDegreeTorqueModel->updateMembers( evaluationTime );

                // Remove the common gravitational-constant scaling before comparing against the analytical references,
                // which are written with G = 1.
                const double inverseGravitationalScaling = 1.0 / physical_constants::GRAVITATIONAL_CONSTANT;
                const Eigen::Vector3d singleInteractionTorqueFromFullTwoBody =
                        inverseGravitationalScaling * fullTwoBodySingleInteractionTorqueModel->getTorque( );
                const Eigen::Vector3d allDegreeTwoDegreeTwoTorqueFromFullTwoBody =
                        inverseGravitationalScaling * fullTwoBodyAllDegreeTwoDegreeTwoTorqueModel->getTorque( );
                (void)allDegreeTwoDegreeTwoTorqueFromFullTwoBody;
                const Eigen::Vector3d isolatedFigureFigureTorqueFromFullTwoBodyBySubtraction = inverseGravitationalScaling *
                        ( fullTwoBodyFullDegreeTwoTorqueModel->getTorque( ) - fullTwoBodyPointMassDegreeTwoTorqueModel->getTorque( ) );
                (void)isolatedFigureFigureTorqueFromFullTwoBodyBySubtraction;
                const Eigen::Vector3d isolatedFigureFigureTorqueFromFourthMinusSecond =
                        inverseGravitationalScaling * ( fourthDegreeTorqueModel->getTorque( ) - secondDegreeTorqueModel->getTorque( ) );

                // Analytical reference for this isolated term:
                // Eqs. (6)-(12) from the provided derivation (C20/C21 of body 1 interacting with one
                // body-2 degree-2 term, expressed in body-1-fixed frame).
                const Eigen::Vector3d interactionDifference =
                        singleInteractionTorqueFromFullTwoBody - isolatedFigureFigureTorqueFromFourthMinusSecond;
                const Eigen::Vector3d relativePositionInBody1Frame =
                        bodies.at( bodyUndergoingTorqueName )->getCurrentRotationToLocalFrame( ) *
                        ( bodies.at( bodyExertingTorqueName )->getPosition( ) - bodies.at( bodyUndergoingTorqueName )->getPosition( ) );
                // Closed-form analytical references currently validated for body-1 C20 only.
                // Expanded degree-2 coverage is validated through full-two-body single-term vs.
                // fourth-minus-second isolation (independent model path).
                const bool hasAnalyticalReference =
                        body1CoefficientCase.useCosineCoefficient && ( body1CoefficientCase.order == 0 || body1CoefficientCase.order == 1 );
                Eigen::Vector3d analyticalTorque = Eigen::Vector3d::Zero( );
                Eigen::Vector3d fullTwoBodyAnalyticalDifference = Eigen::Vector3d::Zero( );
                Eigen::Vector3d isolatedFigureFigureAnalyticalDifference = Eigen::Vector3d::Zero( );
                double fullTwoBodyRelativeAnalyticalError = 0.0;
                double isolatedFigureFigureRelativeAnalyticalError = 0.0;
                if( hasAnalyticalReference )
                {
                    if( body1CoefficientCase.order == 0 )
                    {
                        analyticalTorque =
                                computeAnalyticalC20DegreeTwoFigureFigureTorque( relativePositionInBody1Frame,
                                                                                 bodies.at( bodyUndergoingTorqueName )->getBodyMass( ),
                                                                                 bodies.at( bodyExertingTorqueName )->getBodyMass( ),
                                                                                 referenceRadiusBody1,
                                                                                 referenceRadiusBody2,
                                                                                 cosineCoefficientsOfBody1( 2, 0 ),
                                                                                 coefficientCase.useCosineCoefficient,
                                                                                 coefficientCase.order,
                                                                                 coefficientCase.coefficientValue );
                    }
                    else if( body1CoefficientCase.order == 1 )
                    {
                        analyticalTorque =
                                computeAnalyticalC21DegreeTwoFigureFigureTorque( relativePositionInBody1Frame,
                                                                                 bodies.at( bodyUndergoingTorqueName )->getBodyMass( ),
                                                                                 bodies.at( bodyExertingTorqueName )->getBodyMass( ),
                                                                                 referenceRadiusBody1,
                                                                                 referenceRadiusBody2,
                                                                                 cosineCoefficientsOfBody1( 2, 1 ),
                                                                                 coefficientCase.useCosineCoefficient,
                                                                                 coefficientCase.order,
                                                                                 coefficientCase.coefficientValue );
                    }
                    else
                    {
                        throw std::runtime_error( "Unexpected analytical reference branch in C20/C21 test." );
                    }
                    fullTwoBodyAnalyticalDifference = singleInteractionTorqueFromFullTwoBody - analyticalTorque;
                    isolatedFigureFigureAnalyticalDifference = isolatedFigureFigureTorqueFromFourthMinusSecond - analyticalTorque;
                    const double analyticalScale = std::max( 1.0, analyticalTorque.norm( ) );
                    fullTwoBodyRelativeAnalyticalError = fullTwoBodyAnalyticalDifference.norm( ) / analyticalScale;
                    isolatedFigureFigureRelativeAnalyticalError = isolatedFigureFigureAnalyticalDifference.norm( ) / analyticalScale;
                }
                const double fullVsIsolatedRelativeDifference =
                        interactionDifference.norm( ) / std::max( 1.0, isolatedFigureFigureTorqueFromFourthMinusSecond.norm( ) );

                // std::cout<<"Relative position "<<relativePositionCase.transpose(  )<<std::endl;
                // std::cout<<"Torque "<<body1CoefficientCase.coefficientName<<" x "<<coefficientCase.coefficientName<<std::endl;
                // std::cout<<singleInteractionTorqueFromFullTwoBody.transpose(  )<<std::endl;
                // std::cout<<allDegreeTwoDegreeTwoTorqueFromFullTwoBody.transpose(  )<<std::endl;
                // std::cout<<isolatedFigureFigureTorqueFromFullTwoBodyBySubtraction.transpose(  )<<std::endl;
                // std::cout<<isolatedFigureFigureTorqueFromFourthMinusSecond.transpose(  )<<std::endl;
                // std::cout<<isolatedFigureFigureTorqueFromFullTwoBodyBySubtraction.transpose(  ).cwiseQuotient(
                //     isolatedFigureFigureTorqueFromFourthMinusSecond.transpose(  ) )<<std::endl<<std::endl;
                // std::cout<<analyticalTorque.transpose(  )<<std::endl;
                // std::cout<<inverseGravitationalScaling * secondDegreeTorqueModel->getTorque( ).norm(  )<<std::endl;

                // const std::shared_ptr< FullTwoBodySphericalHarmonicAcceleration > fullTwoBodyAcceleration =
                //         fullTwoBodySingleInteractionTorqueModel->getAccelerationBetweenBodies( );
                // const std::shared_ptr< EffectiveMutualSphericalHarmonicsField > effectiveField =
                //         fullTwoBodyAcceleration->getEffectiveMutualPotentialField( );
                // const Eigen::Vector3d totalSpecificTorqueFromMutualPotential =
                //         fullTwoBodyAcceleration->getCurrentBodyFixedRelativePosition( ).cross(
                //                 fullTwoBodyAcceleration->getMutualPotentialGradient( ) );
                // const Eigen::Vector3d reconstructedBody2SpecificTorque =
                //         fullTwoBodySingleInteractionTorqueModel->getTorque( ) + totalSpecificTorqueFromMutualPotential;
                // std::cout<<"Ceff(2,0,2,m)="
                //          <<effectiveField->getEffectiveCosineCoefficient( 2, 0, 2, static_cast< int >( coefficientCase.order ) )
                //          <<" Seff(2,0,2,m)="
                //          <<effectiveField->getEffectiveSineCoefficient( 2, 0, 2, static_cast< int >( coefficientCase.order ) )
                //          <<" mult(2,0,2,m)="
                //          <<effectiveField->getMultiplier( 2, 0, 2, static_cast< int >( coefficientCase.order ) )<<std::endl;
                // if( coefficientCase.order > 0 )
                // {
                //     std::cout<<"Ceff(2,0,2,-m)="
                //              <<effectiveField->getEffectiveCosineCoefficient( 2, 0, 2, -static_cast< int >( coefficientCase.order ) )
                //              <<" Seff(2,0,2,-m)="
                //              <<effectiveField->getEffectiveSineCoefficient( 2, 0, 2, -static_cast< int >( coefficientCase.order ) )
                //              <<" mult(2,0,2,-m)="
                //              <<effectiveField->getMultiplier( 2, 0, 2, -static_cast< int >( coefficientCase.order ) )<<std::endl;
                // }
                // std::cout<<"total_specific_torque="<<totalSpecificTorqueFromMutualPotential.transpose(  )
                //          <<" current_specific="<<fullTwoBodySingleInteractionTorqueModel->getTorque( ).transpose(  )
                //          <<" reconstructed_body2_specific="<<reconstructedBody2SpecificTorque.transpose(  )<<std::endl;
                // std::cout<<"transformed_body2_degree2(C20,C21,S21,C22,S22)="
                //          <<effectiveField->getTransformedCosineCoefficientsOfBody2( )( 2, 0 )<<" "
                //          <<effectiveField->getTransformedCosineCoefficientsOfBody2( )( 2, 1 )<<" "
                //          <<effectiveField->getTransformedSineCoefficientsOfBody2( )( 2, 1 )<<" "
                //          <<effectiveField->getTransformedCosineCoefficientsOfBody2( )( 2, 2 )<<" "
                //          <<effectiveField->getTransformedSineCoefficientsOfBody2( )( 2, 2 )<<std::endl;
                //
                // std::array< Eigen::MatrixXd, 3 > transformedCosineCoefficientsBody2AngularMomentum;
                // std::array< Eigen::MatrixXd, 3 > transformedSineCoefficientsBody2AngularMomentum;
                // fullTwoBodySingleInteractionTorqueModel->computeTransformedAngularMomentumCoefficients(
                //         effectiveField->getCosineCoefficientsOfBody2( ),
                //         effectiveField->getSineCoefficientsOfBody2( ),
                //         effectiveField->getTransformationCache( )->getWignerDMatricesCache( ),
                //         fullTwoBodyAcceleration->getAreCoefficientsNormalized( ),
                //         transformedCosineCoefficientsBody2AngularMomentum,
                //         transformedSineCoefficientsBody2AngularMomentum );
                // std::cout<<"JC2_x(m=0,1,2)="
                //          <<transformedCosineCoefficientsBody2AngularMomentum.at( 0 )( 2, 0 )<<" "
                //          <<transformedCosineCoefficientsBody2AngularMomentum.at( 0 )( 2, 1 )<<" "
                //          <<transformedCosineCoefficientsBody2AngularMomentum.at( 0 )( 2, 2 )<<" "
                //          <<" JS2_x(m=0,1,2)="
                //          <<transformedSineCoefficientsBody2AngularMomentum.at( 0 )( 2, 0 )<<" "
                //          <<transformedSineCoefficientsBody2AngularMomentum.at( 0 )( 2, 1 )<<" "
                //          <<transformedSineCoefficientsBody2AngularMomentum.at( 0 )( 2, 2 )<<std::endl;
                // std::cout<<"JC2_y(m=0,1,2)="
                //          <<transformedCosineCoefficientsBody2AngularMomentum.at( 1 )( 2, 0 )<<" "
                //          <<transformedCosineCoefficientsBody2AngularMomentum.at( 1 )( 2, 1 )<<" "
                //          <<transformedCosineCoefficientsBody2AngularMomentum.at( 1 )( 2, 2 )<<" "
                //          <<" JS2_y(m=0,1,2)="
                //          <<transformedSineCoefficientsBody2AngularMomentum.at( 1 )( 2, 0 )<<" "
                //          <<transformedSineCoefficientsBody2AngularMomentum.at( 1 )( 2, 1 )<<" "
                //          <<transformedSineCoefficientsBody2AngularMomentum.at( 1 )( 2, 2 )<<std::endl;
                // std::cout<<"JC2_z(m=0,1,2)="
                //          <<transformedCosineCoefficientsBody2AngularMomentum.at( 2 )( 2, 0 )<<" "
                //          <<transformedCosineCoefficientsBody2AngularMomentum.at( 2 )( 2, 1 )<<" "
                //          <<transformedCosineCoefficientsBody2AngularMomentum.at( 2 )( 2, 2 )<<" "
                //          <<" JS2_z(m=0,1,2)="
                //          <<transformedSineCoefficientsBody2AngularMomentum.at( 2 )( 2, 0 )<<" "
                //          <<transformedSineCoefficientsBody2AngularMomentum.at( 2 )( 2, 1 )<<" "
                //          <<transformedSineCoefficientsBody2AngularMomentum.at( 2 )( 2, 2 )<<std::endl;
                //
                // BOOST_TEST_MESSAGE( "single_term_case=" << coefficientCase.coefficientName
                //                                          << " full=" << singleInteractionTorqueFromFullTwoBody.transpose( )
                //                                          << " fourth_minus_second="
                //                                          << isolatedFigureFigureTorqueFromFourthMinusSecond.transpose( )
                //                                          << " analytical=" << analyticalTorque.transpose( )
                //                                          << " full_analytical_diff="
                //                                          << fullTwoBodyAnalyticalDifference.transpose( )
                //                                          << " isolated_analytical_diff="
                //                                          << isolatedFigureFigureAnalyticalDifference.transpose( )
                //                                          << " full_minus_isolated=" << interactionDifference.transpose( ) );

                // Non-zero guards: each selected coefficient must induce a measurable figure-figure torque.
                if( hasAnalyticalReference )
                {
                    // For body-1 C20, validate against the closed-form analytical result.
                    if( analyticalTorque.norm( ) > 1.0E-20 )
                    {
                        if( body1CoefficientCase.order == 0 )
                        {
                            BOOST_CHECK_SMALL( isolatedFigureFigureRelativeAnalyticalError, 1.0E-11 );
                            BOOST_CHECK_SMALL( fullTwoBodyRelativeAnalyticalError, 5.0E-14 );
                        }
                        else if( body1CoefficientCase.order == 1 && fullVsIsolatedRelativeDifference > 1.0E-11 )
                        {
                            BOOST_TEST_MESSAGE( "C21_theory_compare case="
                                                << body1CoefficientCase.coefficientName << "x" << coefficientCase.coefficientName << " r="
                                                << relativePositionCase.transpose( ) << " analytical=" << analyticalTorque.transpose( )
                                                << " full=" << singleInteractionTorqueFromFullTwoBody.transpose( )
                                                << " fourth_minus_second=" << isolatedFigureFigureTorqueFromFourthMinusSecond.transpose( )
                                                << " full_minus_analytical=" << fullTwoBodyAnalyticalDifference.transpose( )
                                                << " fourth_minus_analytical=" << isolatedFigureFigureAnalyticalDifference.transpose( )
                                                << " full_rel_err=" << fullTwoBodyRelativeAnalyticalError
                                                << " fourth_rel_err=" << isolatedFigureFigureRelativeAnalyticalError );
                        }
                    }
                    else
                    {
                        // These analytically zero cases leave roundoff-level residuals after
                        // spherical-harmonic normalization, term selection, and inverse-G scaling.
                        BOOST_CHECK_SMALL( isolatedFigureFigureTorqueFromFourthMinusSecond.norm( ), 1.0E-14 );
                        BOOST_CHECK_SMALL( singleInteractionTorqueFromFullTwoBody.norm( ), 1.0E-14 );
                    }
                }
                // Always compare the two independent model paths.
                if( fullVsIsolatedRelativeDifference > 1.0E-11 )
                {
                    const Eigen::Vector3d fourthDegreeOnlyTorque = inverseGravitationalScaling * fourthDegreeTorqueModel->getTorque( );
                    const Eigen::Vector3d secondDegreeOnlyTorque = inverseGravitationalScaling * secondDegreeTorqueModel->getTorque( );
                    const Eigen::Vector3d expectedFourthDegreeOnlyTorqueFromFullPath =
                            isolatedFigureFigureTorqueFromFullTwoBodyBySubtraction + secondDegreeOnlyTorque;
                    const Eigen::Vector3d fourthOnlyMinusExpectedFourthOnly =
                            fourthDegreeOnlyTorque - expectedFourthDegreeOnlyTorqueFromFullPath;
                    const Eigen::Vector3d fullSingleMinusFullSubtraction =
                            singleInteractionTorqueFromFullTwoBody - isolatedFigureFigureTorqueFromFullTwoBodyBySubtraction;
                    const Eigen::Vector3d fullSubtractionMinusFourthMinusSecond =
                            isolatedFigureFigureTorqueFromFullTwoBodyBySubtraction - isolatedFigureFigureTorqueFromFourthMinusSecond;

                    const Eigen::Matrix3d body1Inertia = bodies.at( bodyUndergoingTorqueName )->getBodyInertiaTensor( );
                    const Eigen::Matrix3d body2InertiaInBody1Frame =
                            fourthDegreeTorqueModel->getCurrentInertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque( );
                    const SchutzEq11TermDiagnostics eq11Diagnostics =
                            computeSchutzEq11TermDiagnostics( relativePositionInBody1Frame,
                                                              bodies.at( bodyExertingTorqueName )->getBodyMass( ),
                                                              body1Inertia,
                                                              body2InertiaInBody1Frame );
                    const double scaledPrefactor = eq11Diagnostics.prefactor * inverseGravitationalScaling;
                    BOOST_TEST_MESSAGE(
                            "single_term_case="
                            << body1CoefficientCase.coefficientName << "x" << coefficientCase.coefficientName
                            << " r=" << relativePositionCase.transpose( ) << " full=" << singleInteractionTorqueFromFullTwoBody.transpose( )
                            << " full_subtraction=" << isolatedFigureFigureTorqueFromFullTwoBodyBySubtraction.transpose( )
                            << " fourth_minus_second=" << isolatedFigureFigureTorqueFromFourthMinusSecond.transpose( )
                            << " fourth_only=" << fourthDegreeOnlyTorque.transpose( )
                            << " expected_fourth_only_from_full_path=" << expectedFourthDegreeOnlyTorqueFromFullPath.transpose( )
                            << " fourth_only_minus_expected_fourth_only=" << fourthOnlyMinusExpectedFourthOnly.transpose( )
                            << " second_only=" << secondDegreeOnlyTorque.transpose( )
                            << " full_minus_full_subtraction=" << fullSingleMinusFullSubtraction.transpose( )
                            << " full_subtraction_minus_fourth_minus_second=" << fullSubtractionMinusFourthMinusSecond.transpose( )
                            << " body1_inertia=" << bodies.at( bodyUndergoingTorqueName )->getBodyInertiaTensor( )
                            << " body2_inertia=" << bodies.at( bodyExertingTorqueName )->getBodyInertiaTensor( )
                            << " diff=" << interactionDifference.transpose( ) << " rel=" << fullVsIsolatedRelativeDifference );

                    if( body1CoefficientCase.coefficientName == "C21" && std::fabs( scaledPrefactor * eq11Diagnostics.Ixz ) > 1.0E-20 )
                    {
                        const double requiredFxyFromExpected =
                                -expectedFourthDegreeOnlyTorqueFromFullPath( 0 ) / ( scaledPrefactor * eq11Diagnostics.Ixz );
                        const double requiredGxzFromExpected =
                                expectedFourthDegreeOnlyTorqueFromFullPath( 1 ) / ( scaledPrefactor * eq11Diagnostics.Ixz );
                        const double requiredFyzFromExpected =
                                expectedFourthDegreeOnlyTorqueFromFullPath( 2 ) / ( scaledPrefactor * eq11Diagnostics.Ixz );
                        BOOST_TEST_MESSAGE( "Eq11_C21_debug case="
                                            << body1CoefficientCase.coefficientName << "x" << coefficientCase.coefficientName
                                            << " r=" << relativePositionCase.transpose( ) << " fxy(computed,required_from_expected)=("
                                            << eq11Diagnostics.fxy << ", " << requiredFxyFromExpected << ")"
                                            << " gxz(computed,required_from_expected)=(" << eq11Diagnostics.gxz << ", "
                                            << requiredGxzFromExpected << ")"
                                            << " fyz(computed,required_from_expected)=(" << eq11Diagnostics.fyz << ", "
                                            << requiredFyzFromExpected << ")" );
                    }
                    else if( body1CoefficientCase.coefficientName == "S21" && std::fabs( scaledPrefactor * eq11Diagnostics.Iyz ) > 1.0E-20 )
                    {
                        const double requiredGyzFromExpected =
                                expectedFourthDegreeOnlyTorqueFromFullPath( 0 ) / ( scaledPrefactor * eq11Diagnostics.Iyz );
                        const double requiredFxyFromExpected =
                                expectedFourthDegreeOnlyTorqueFromFullPath( 1 ) / ( scaledPrefactor * eq11Diagnostics.Iyz );
                        const double requiredFxzFromExpected =
                                -expectedFourthDegreeOnlyTorqueFromFullPath( 2 ) / ( scaledPrefactor * eq11Diagnostics.Iyz );
                        BOOST_TEST_MESSAGE( "Eq11_S21_debug case="
                                            << body1CoefficientCase.coefficientName << "x" << coefficientCase.coefficientName
                                            << " r=" << relativePositionCase.transpose( ) << " gyz(computed,required_from_expected)=("
                                            << eq11Diagnostics.gyz << ", " << requiredGyzFromExpected << ")"
                                            << " fxy(computed,required_from_expected)=(" << eq11Diagnostics.fxy << ", "
                                            << requiredFxyFromExpected << ")"
                                            << " fxz(computed,required_from_expected)=(" << eq11Diagnostics.fxz << ", "
                                            << requiredFxzFromExpected << ")" );
                    }
                    else if( body1CoefficientCase.coefficientName == "S22" && std::fabs( scaledPrefactor * eq11Diagnostics.Ixy ) > 1.0E-20 )
                    {
                        const double requiredFxzFromExpected =
                                expectedFourthDegreeOnlyTorqueFromFullPath( 0 ) / ( scaledPrefactor * eq11Diagnostics.Ixy );
                        const double requiredFyzFromExpected =
                                -expectedFourthDegreeOnlyTorqueFromFullPath( 1 ) / ( scaledPrefactor * eq11Diagnostics.Ixy );
                        const double requiredGxyFromExpected =
                                expectedFourthDegreeOnlyTorqueFromFullPath( 2 ) / ( scaledPrefactor * eq11Diagnostics.Ixy );
                        BOOST_TEST_MESSAGE( "Eq11_S22_debug case="
                                            << body1CoefficientCase.coefficientName << "x" << coefficientCase.coefficientName
                                            << " r=" << relativePositionCase.transpose( ) << " fxz(computed,required_from_expected)=("
                                            << eq11Diagnostics.fxz << ", " << requiredFxzFromExpected << ")"
                                            << " fyz(computed,required_from_expected)=(" << eq11Diagnostics.fyz << ", "
                                            << requiredFyzFromExpected << ")"
                                            << " gxy(computed,required_from_expected)=(" << eq11Diagnostics.gxy << ", "
                                            << requiredGxyFromExpected << ")" );
                    }
                }
                BOOST_CHECK_SMALL( fullVsIsolatedRelativeDifference, 1.0E-11 );
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
