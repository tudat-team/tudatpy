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

std::shared_ptr< simulation_setup::Body > createBodyForFullTwoBodyTorqueTest(
        const double gravitationalParameter,
        const double referenceRadius,
        const Eigen::MatrixXd& cosineCoefficients,
        const Eigen::MatrixXd& sineCoefficients,
        const Eigen::Vector3d& position,
        const Eigen::Quaterniond& rotationToBodyFixedFrame,
        const double scaledMeanMomentOfInertia = TUDAT_NAN )
{
    std::shared_ptr< simulation_setup::Body > body = std::make_shared< simulation_setup::Body >( );
    body->setGravityFieldModel( std::make_shared< gravitation::SphericalHarmonicsGravityField >(
            gravitationalParameter,
            referenceRadius,
            cosineCoefficients,
            sineCoefficients,
            "BodyFixed",
            scaledMeanMomentOfInertia ) );

    Eigen::Vector6d bodyState = Eigen::Vector6d::Zero( );
    bodyState.segment( 0, 3 ) = position;
    body->setState( bodyState );

    body->setRotationalEphemeris( std::make_shared< ephemerides::ConstantRotationalEphemeris >(
            rotationToBodyFixedFrame.inverse( ), "ECLIPJ2000", "BodyFixed" ) );

    return body;
}

simulation_setup::SystemOfBodies createSystemOfBodiesForFullTwoBodyTorqueTest(
        const std::string& bodyUndergoingTorqueName,
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
    bodies.addBody(
            createBodyForFullTwoBodyTorqueTest( gravitationalParameter,
                                                referenceRadiusBody1,
                                                cosineCoefficientsOfBody1,
                                                sineCoefficientsOfBody1,
                                                positionOfBody1,
                                                rotationToBody1,
                                                scaledMeanMomentOfInertiaBody1 ),
            bodyUndergoingTorqueName,
            false );
    bodies.addBody(
            createBodyForFullTwoBodyTorqueTest( gravitationalParameter,
                                                referenceRadiusBody2,
                                                cosineCoefficientsOfBody2,
                                                sineCoefficientsOfBody2,
                                                positionOfBody2,
                                                rotationToBody2,
                                                scaledMeanMomentOfInertiaBody2 ),
            bodyExertingTorqueName,
            false );
    bodies.processBodyFrameDefinitions< >( );
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
            simulation_setup::fullTwoBodySphericalHarmonicGravitationalTorque(
                    maximumDegreeOfBodyUndergoingTorque,
                    maximumOrderOfBodyUndergoingTorque,
                    maximumDegreeOfBodyExertingTorque,
                    maximumOrderOfBodyExertingTorque ) );

    basic_astrodynamics::TorqueModelMap torqueModelMap = simulation_setup::createTorqueModelsMap(
            bodies, selectedTorqueModelMap, { bodyUndergoingTorqueName } );
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

    basic_astrodynamics::TorqueModelMap torqueModelMap = simulation_setup::createTorqueModelsMap(
            bodies, selectedTorqueModelMap, { bodyUndergoingTorqueName } );
    return std::dynamic_pointer_cast< gravitation::FullTwoBodySphericalHarmonicTorque >(
            torqueModelMap.at( bodyUndergoingTorqueName ).at( bodyExertingTorqueName ).at( 0 ) );
}

std::shared_ptr< gravitation::FullTwoBodySphericalHarmonicAcceleration > createReferenceFullTwoBodySphericalHarmonicAccelerationModel(
        const simulation_setup::SystemOfBodies& bodies,
        const std::string& bodyUndergoingTorqueName,
        const std::string& bodyExertingTorqueName )
{
    return std::dynamic_pointer_cast< gravitation::FullTwoBodySphericalHarmonicAcceleration >(
            simulation_setup::createAccelerationModel(
                    bodies.at( bodyUndergoingTorqueName ),
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

    basic_astrodynamics::TorqueModelMap torqueModelMap = simulation_setup::createTorqueModelsMap(
            bodies, selectedTorqueModelMap, { bodyUndergoingTorqueName } );
    return std::dynamic_pointer_cast< gravitation::SphericalHarmonicGravitationalTorqueModel >(
            torqueModelMap.at( bodyUndergoingTorqueName ).at( bodyExertingTorqueName ).at( 0 ) );
}

std::shared_ptr< simulation_setup::Body > createBodyForFourthDegreeTwoBodyTorqueTest(
        const double bodyMass,
        const Eigen::Matrix3d& inertiaTensorInBodyFixedFrame,
        const Eigen::Vector3d& position,
        const Eigen::Quaterniond& rotationToBodyFixedFrame )
{
    std::shared_ptr< simulation_setup::Body > body = std::make_shared< simulation_setup::Body >( );
    body->setConstantBodyMass( bodyMass );
    body->setBodyInertiaTensor( inertiaTensorInBodyFixedFrame );
    body->setGravityFieldModel( std::make_shared< gravitation::GravityFieldModel >(
            physical_constants::GRAVITATIONAL_CONSTANT * bodyMass ) );

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
    bodies.addBody(
            createBodyForFourthDegreeTwoBodyTorqueTest( massOfBodyUndergoingTorque,
                                                        inertiaTensorOfBodyUndergoingTorque,
                                                        positionOfBodyUndergoingTorque,
                                                        rotationToBodyUndergoingTorque ),
            bodyUndergoingTorqueName,
            false );
    bodies.addBody(
            createBodyForFourthDegreeTwoBodyTorqueTest( massOfBodyExertingTorque,
                                                        inertiaTensorOfBodyExertingTorque,
                                                        positionOfBodyExertingTorque,
                                                        rotationToBodyExertingTorque ),
            bodyExertingTorqueName,
            false );
    bodies.processBodyFrameDefinitions< >( );
    return bodies;
}

std::shared_ptr< gravitation::FourthDegreeFullTwoBodyGravitationalTorqueModel >
createFactoryFourthDegreeFullTwoBodyGravitationalTorqueModel(
        const simulation_setup::SystemOfBodies& bodies,
        const std::string& bodyUndergoingTorqueName,
        const std::string& bodyExertingTorqueName )
{
    simulation_setup::SelectedTorqueMap selectedTorqueModelMap;
    selectedTorqueModelMap[ bodyUndergoingTorqueName ][ bodyExertingTorqueName ].push_back(
            simulation_setup::fourthDegreeFullTwoBodyGravitationalTorque( ) );

    basic_astrodynamics::TorqueModelMap torqueModelMap = simulation_setup::createTorqueModelsMap(
            bodies, selectedTorqueModelMap, { bodyUndergoingTorqueName } );
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

    basic_astrodynamics::TorqueModelMap torqueModelMap = simulation_setup::createTorqueModelsMap(
            bodies, selectedTorqueModelMap, { bodyUndergoingTorqueName } );
    return std::dynamic_pointer_cast< gravitation::SecondDegreeGravitationalTorqueModel >(
            torqueModelMap.at( bodyUndergoingTorqueName ).at( bodyExertingTorqueName ).at( 0 ) );
}

Eigen::Matrix3d transformBodyExertingInertiaTensorToBodyUndergoingFrame(
        const Eigen::Matrix3d& inertiaTensorOfBodyExertingTorque,
        const Eigen::Quaterniond& rotationToBodyUndergoingTorque,
        const Eigen::Quaterniond& rotationToBodyExertingTorque )
{
    const Eigen::Matrix3d rotationFromBodyExertingToBodyUndergoing =
            rotationToBodyUndergoingTorque.toRotationMatrix( ) *
            rotationToBodyExertingTorque.toRotationMatrix( ).transpose( );
    return rotationFromBodyExertingToBodyUndergoing * inertiaTensorOfBodyExertingTorque *
            rotationFromBodyExertingToBodyUndergoing.transpose( );
}

Eigen::Vector3d computeManualFourthDegreeTwoBodyTorqueFromBodyStates(
        const std::shared_ptr< simulation_setup::Body >& bodyUndergoingTorque,
        const std::shared_ptr< simulation_setup::Body >& bodyExertingTorque )
{
    const Eigen::Vector3d relativePositionInBodyUndergoingFrame = bodyUndergoingTorque->getCurrentRotationToLocalFrame( ) *
            ( bodyExertingTorque->getPosition( ) - bodyUndergoingTorque->getPosition( ) );
    const Eigen::Matrix3d inertiaTensorOfBodyExertingInBodyUndergoingFrame =
            transformBodyExertingInertiaTensorToBodyUndergoingFrame(
                    bodyExertingTorque->getBodyInertiaTensor( ),
                    bodyUndergoingTorque->getCurrentRotationToLocalFrame( ),
                    bodyExertingTorque->getCurrentRotationToLocalFrame( ) );

    return gravitation::calculateFourthDegreeFullTwoBodyGravitationalTorque(
            relativePositionInBodyUndergoingFrame,
            bodyExertingTorque->getBodyMass( ),
            bodyUndergoingTorque->getBodyInertiaTensor( ),
            inertiaTensorOfBodyExertingInBodyUndergoingFrame );
}

std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >
getPointMassDegreeTwoInteractionCombinations( )
{
    std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > coefficientCombinations;
    for( unsigned int m = 0; m <= 2; m++ )
    {
        coefficientCombinations.push_back( std::make_tuple( 2, m, 0, 0 ) );
    }
    return coefficientCombinations;
}

std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >
getDegreeTwoDegreeTwoInteractionCombinations( )
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

std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >
getFullDegreeTwoInteractionCombinations( )
{
    std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > coefficientCombinations =
            getPointMassDegreeTwoInteractionCombinations( );
    const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >
            degreeTwoDegreeTwoCombinations = getDegreeTwoDegreeTwoInteractionCombinations( );
    coefficientCombinations.insert(
            coefficientCombinations.end( ), degreeTwoDegreeTwoCombinations.begin( ), degreeTwoDegreeTwoCombinations.end( ) );
    return coefficientCombinations;
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

    const std::string bodyUndergoingTorqueName = "Body1";
    const std::string bodyExertingTorqueName = "Body2";
    const Eigen::Vector3d positionOfBodyUndergoingTorque( 1.5E6, -2.4E6, 3.2E6 );
    const Eigen::Vector3d positionOfBodyExertingTorque( -4.1E6, 2.7E6, 1.1E6 );
    const double evaluationTime = 43200.0;

    const double massOfBodyUndergoingTorque = 5.8E21;
    const double massOfBodyExertingTorque = 7.1E21;

    const Eigen::Matrix3d inertiaTensorOfBodyUndergoingTorque = ( Eigen::Matrix3d( ) << 3.7E29, -1.2E27, 2.0E27,
                                                                   -1.2E27, 4.4E29, -0.8E27,
                                                                   2.0E27, -0.8E27, 5.1E29 )
                                                                          .finished( );
    const Eigen::Quaterniond rotationToBodyUndergoingTorque(
            Eigen::AngleAxisd( 0.42, Eigen::Vector3d::UnitX( ) ) *
            Eigen::AngleAxisd( -0.35, Eigen::Vector3d::UnitY( ) ) *
            Eigen::AngleAxisd( 0.71, Eigen::Vector3d::UnitZ( ) ) );

    // Case 0: point-mass-equivalent body A (isotropic inertia tensor); verify Case 1 and Eq. (14) Case 2.
    {
        const double gravitationalParameterPointMassEquivalentBody =
                physical_constants::GRAVITATIONAL_CONSTANT * massOfBodyExertingTorque;
        const double referenceRadiusBodyUndergoingTorque = 2.3E6;
        const double referenceRadiusBodyExertingTorque = 1.4E6;

        const double isotropicInertiaValueOfBodyExertingTorque = 1.0E20;
        const double scaledMeanMomentOfInertiaBodyExertingTorque =
                isotropicInertiaValueOfBodyExertingTorque /
                ( ( gravitationalParameterPointMassEquivalentBody / physical_constants::GRAVITATIONAL_CONSTANT ) *
                  referenceRadiusBodyExertingTorque * referenceRadiusBodyExertingTorque );

        const std::tuple< Eigen::MatrixXd, Eigen::MatrixXd, double > degreeTwoCoefficientsOfBodyUndergoingTorque =
                gravitation::getDegreeTwoSphericalHarmonicCoefficients(
                        inertiaTensorOfBodyUndergoingTorque,
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
        orientationOfPointMassEquivalentBodyCases.push_back(
                Eigen::Quaterniond( Eigen::AngleAxisd( 0.62, Eigen::Vector3d::UnitY( ) ) *
                                    Eigen::AngleAxisd( -0.48, Eigen::Vector3d::UnitX( ) ) *
                                    Eigen::AngleAxisd( 0.29, Eigen::Vector3d::UnitZ( ) ) ) );

        std::vector< Eigen::Vector3d > case1Torques;
        for( unsigned int i = 0; i < orientationOfPointMassEquivalentBodyCases.size( ); i++ )
        {
            const SystemOfBodies bodies = createSystemOfBodiesForFullTwoBodyTorqueTest(
                    bodyUndergoingTorqueName,
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
            const Eigen::Vector3d case1SpecificTorqueFromFourthDegree =
                    case1TorqueFromFourthDegree / massOfBodyExertingTorqueInCurrentCase;
            const Eigen::Vector3d totalTorqueFromMutualPotential =
                    referenceAcceleration->getCurrentBodyFixedRelativePosition( ).cross(
                            referenceAcceleration->getMutualPotentialGradient( ) );
            const Eigen::Vector3d case2TorqueOnPointMassEquivalentBodyFromEq14 =
                    totalTorqueFromMutualPotential + case1SpecificTorqueFromFourthDegree;
            const double case1Scale = std::max( 1.0, case1TorqueFromSecondDegree.norm( ) );
            const double case2Scale = std::max( 1.0, totalTorqueFromMutualPotential.norm( ) );

            // This check verifies Case 1: torque of point-mass-equivalent A on extended B equals second-degree torque.
            BOOST_CHECK_SMALL(
                    ( case1TorqueFromFourthDegree - case1TorqueFromSecondDegree ).norm( ) / case1Scale,
                    5.0E-13 );
            // This check verifies Case 2 (Eq. 14): torque of extended B on point-mass-equivalent A is zero.
            BOOST_CHECK_SMALL( case2TorqueOnPointMassEquivalentBodyFromEq14.norm( ) / case2Scale, 5.0E-12 );

            case1Torques.push_back( case1TorqueFromFourthDegree );
        }

        const double orientationInvariantScale = std::max( 1.0, case1Torques.at( 0 ).norm( ) );
        // This check confirms isotropic inertia makes Case 1 torque invariant to body-A orientation.
        BOOST_CHECK_SMALL( ( case1Torques.at( 0 ) - case1Torques.at( 1 ) ).norm( ) / orientationInvariantScale, 5.0E-13 );
    }

    // Case 1: body 2 is a point-mass equivalent; Eq. (11) must reduce to second-degree torque.
    {
        const Eigen::Matrix3d inertiaTensorOfBodyExertingTorque = Eigen::Matrix3d::Zero( );
        const Eigen::Quaterniond rotationToBodyExertingTorque = Eigen::Quaterniond::Identity( );

        const SystemOfBodies bodies = createSystemOfBodiesForFourthDegreeTwoBodyTorqueTest(
                bodyUndergoingTorqueName,
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
                createFactoryFourthDegreeFullTwoBodyGravitationalTorqueModel(
                        bodies, bodyUndergoingTorqueName, bodyExertingTorqueName );
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
        BOOST_CHECK_SMALL( ( fourthDegreeTorque - secondDegreeTorque ).norm( ) / referenceScale, 5.0E-14 );
        // This check confirms that a zero body-2 inertia remains zero after frame transformation.
        BOOST_CHECK_SMALL(
                fourthDegreeTorqueModel->getCurrentInertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque( ).norm( ),
                1.0E-20 );
    }

    // Case 2: body-2 has finite inertia; torque must match Eq. (11) with transformed inertia and vary with orientation.
    {
        const Eigen::Matrix3d inertiaTensorOfBodyExertingTorque = ( Eigen::Matrix3d( ) << 2.1E29, -0.6E27, 1.1E27,
                                                                     -0.6E27, 2.5E29, -1.4E27,
                                                                     1.1E27, -1.4E27, 3.0E29 )
                                                                            .finished( );
        std::vector< Eigen::Quaterniond > bodyExertingTorqueRotations;
        bodyExertingTorqueRotations.push_back( Eigen::Quaterniond::Identity( ) );
        bodyExertingTorqueRotations.push_back(
                Eigen::Quaterniond( Eigen::AngleAxisd( -0.53, Eigen::Vector3d::UnitY( ) ) *
                                    Eigen::AngleAxisd( 0.67, Eigen::Vector3d::UnitZ( ) ) *
                                    Eigen::AngleAxisd( 0.21, Eigen::Vector3d::UnitX( ) ) ) );

        std::vector< Eigen::Vector3d > evaluatedTorques;
        for( unsigned int i = 0; i < bodyExertingTorqueRotations.size( ); i++ )
        {
            const SystemOfBodies bodies = createSystemOfBodiesForFourthDegreeTwoBodyTorqueTest(
                    bodyUndergoingTorqueName,
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
                            inertiaTensorOfBodyExertingTorque,
                            rotationToBodyUndergoingTorque,
                            bodyExertingTorqueRotations.at( i ) );
            const Eigen::Vector3d manualTorque = computeManualFourthDegreeTwoBodyTorqueFromBodyStates(
                    bodies.at( bodyUndergoingTorqueName ), bodies.at( bodyExertingTorqueName ) );
            const Eigen::Vector3d modelTorque = fourthDegreeTorqueModel->getTorque( );
            const double referenceScale = std::max( 1.0, manualTorque.norm( ) );

            // This check verifies the frame transformation of body-2 inertia used internally by the model.
            BOOST_CHECK_SMALL(
                    ( fourthDegreeTorqueModel->getCurrentInertiaTensorOfBodyExertingTorqueInFrameOfBodyUndergoingTorque( ) -
                      manualTransformedInertiaTensorOfBodyExertingTorque )
                            .norm( ) /
                            std::max( 1.0, manualTransformedInertiaTensorOfBodyExertingTorque.norm( ) ),
                    1.0E-14 );
            // This check validates the model output against a direct Eq. (11) evaluation from current body states.
            BOOST_CHECK_SMALL( ( modelTorque - manualTorque ).norm( ) / referenceScale, 5.0E-14 );

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

    // Case 2a/2b: degree-2 undergoing body with point-mass perturber, compare against acceleration model reference.
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
        orientationCases.push_back( std::make_pair(
                Eigen::Quaterniond( Eigen::AngleAxisd( 1.11, Eigen::Vector3d::UnitZ( ) ) *
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

            const Eigen::Vector3d referenceTorqueFromAcceleration =
                    referenceAcceleration->getCurrentBodyFixedRelativePosition( ).cross(
                            referenceAcceleration->getMutualPotentialGradient( ) );
            const Eigen::Vector3d computedFullTwoBodyTorque = torqueModel->getTorque( );
            const Eigen::Vector3d computedSphericalHarmonicTorque = sphericalHarmonicTorqueModel->getTorque( );
            const double bodyExertingTorqueMass = bodies.at( bodyExertingTorqueName )->getBodyMass( );
            const Eigen::Vector3d specificSphericalHarmonicTorque =
                    computedSphericalHarmonicTorque / bodyExertingTorqueMass;
            computedTorques.push_back( computedFullTwoBodyTorque );

            // Consistency of the full two-body torque with the independent acceleration-based torque reference
            // r x (dU/dr), all expressed in body-1-fixed coordinates.
            const Eigen::Vector3d accelerationConsistencyDifference =
                    computedFullTwoBodyTorque - referenceTorqueFromAcceleration;

            // Direct subtraction between full two-body torque and spherical-harmonic torque (specific form):
            // this should stay large because the two models use opposite sign conventions.
            const Eigen::Vector3d modelToModelDifference =
                    computedFullTwoBodyTorque - specificSphericalHarmonicTorque;

            // Signed consistency check for the sign convention: the two models should match after adding
            // (i.e. full-two-body torque ≈ - specific spherical-harmonic torque).
            const Eigen::Vector3d modelToModelSignedDifference =
                    computedFullTwoBodyTorque + specificSphericalHarmonicTorque;
            const double referenceScale = std::max( 1.0, referenceTorqueFromAcceleration.norm( ) );

            for( int i = 0; i < 3; i++ )
            {
                BOOST_CHECK_SMALL( std::fabs( accelerationConsistencyDifference( i ) ) / referenceScale, 5.0E-14 );
            }
            for( int i = 0; i < 3; i++ )
            {
                BOOST_CHECK_SMALL( std::fabs( modelToModelSignedDifference( i ) ) / referenceScale, 5.0E-13 );
            }
            BOOST_CHECK_GT( modelToModelDifference.norm( ) / referenceScale, 1.0 );
        }

        BOOST_CHECK_GT( ( computedTorques.at( 0 ) - computedTorques.at( 1 ) ).norm( ), 1.0E-16 );
    }

    // Case 3: both bodies degree-2; isolate the degree-2/degree-2 coupling and compare both independent models.
    {
        Eigen::MatrixXd cosineCoefficientsOfBody1 = Eigen::MatrixXd::Zero( 3, 3 );
        Eigen::MatrixXd sineCoefficientsOfBody1 = Eigen::MatrixXd::Zero( 3, 3 );
        cosineCoefficientsOfBody1( 0, 0 ) = 1.0;
        cosineCoefficientsOfBody1( 2, 0 ) = 1.1E-3;
        cosineCoefficientsOfBody1( 2, 1 ) = -2.1E-4;
        cosineCoefficientsOfBody1( 2, 2 ) = 3.4E-4;
        sineCoefficientsOfBody1( 2, 1 ) = 1.3E-4;
        sineCoefficientsOfBody1( 2, 2 ) = -2.8E-4;

        Eigen::MatrixXd cosineCoefficientsOfBody2 = Eigen::MatrixXd::Zero( 3, 3 );
        Eigen::MatrixXd sineCoefficientsOfBody2 = Eigen::MatrixXd::Zero( 3, 3 );
        cosineCoefficientsOfBody2( 0, 0 ) = 1.0;
        cosineCoefficientsOfBody2( 2, 0 ) = -6.0E-4;
        cosineCoefficientsOfBody2( 2, 1 ) = 2.3E-4;
        cosineCoefficientsOfBody2( 2, 2 ) = -4.1E-4;
        sineCoefficientsOfBody2( 2, 1 ) = -1.7E-4;
        sineCoefficientsOfBody2( 2, 2 ) = 3.4E-4;

        const Eigen::Vector3d distantPositionOfBody1( 7.2E6, -1.4E6, 2.6E6 );
        const Eigen::Vector3d distantPositionOfBody2( -6.8E6, 3.1E6, -1.9E6 );
        const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >
                pointMassDegreeTwoCombinations = getPointMassDegreeTwoInteractionCombinations( );
        const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >
                degreeTwoDegreeTwoCombinations = getDegreeTwoDegreeTwoInteractionCombinations( );
        const std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >
                fullDegreeTwoCombinations = getFullDegreeTwoInteractionCombinations( );

        std::vector< std::pair< Eigen::Quaterniond, Eigen::Quaterniond > > orientationCases;
        orientationCases.push_back( std::make_pair( Eigen::Quaterniond::Identity( ), Eigen::Quaterniond::Identity( ) ) );
        orientationCases.push_back( std::make_pair(
                Eigen::Quaterniond( Eigen::AngleAxisd( 0.62, Eigen::Vector3d::UnitX( ) ) *
                                    Eigen::AngleAxisd( -0.28, Eigen::Vector3d::UnitZ( ) ) *
                                    Eigen::AngleAxisd( 0.45, Eigen::Vector3d::UnitY( ) ) ),
                Eigen::Quaterniond( Eigen::AngleAxisd( -0.39, Eigen::Vector3d::UnitY( ) ) *
                                    Eigen::AngleAxisd( 0.77, Eigen::Vector3d::UnitX( ) ) *
                                    Eigen::AngleAxisd( 0.21, Eigen::Vector3d::UnitZ( ) ) ) ) );

        std::vector< Eigen::Vector3d > isolatedDegree22TorquesFromFullTwoBodyModel;
        std::vector< Eigen::Vector3d > isolatedDegree22TorquesFromFourthDegreeModel;
        for( const std::pair< Eigen::Quaterniond, Eigen::Quaterniond >& orientationCase : orientationCases )
        {
            const SystemOfBodies bodiesWithAllDegree2Terms = createSystemOfBodiesForFullTwoBodyTorqueTest(
                    bodyUndergoingTorqueName,
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
                    0.3,
                    0.3 );

            bodiesWithAllDegree2Terms.at( bodyUndergoingTorqueName )->setCurrentRotationalStateToLocalFrameFromEphemeris( evaluationTime );
            bodiesWithAllDegree2Terms.at( bodyExertingTorqueName )->setCurrentRotationalStateToLocalFrameFromEphemeris( evaluationTime );

            // Full two-body spherical-harmonic torque including both point-mass/degree-2 and degree-2/degree-2 interactions;
            // used as the "all degree-2 terms" reference for decomposition.
            std::shared_ptr< FullTwoBodySphericalHarmonicTorque > fullTwoBodyTorqueModelWithFullDegree2Terms =
                    createFactoryFullTwoBodySphericalHarmonicTorqueModel(
                            bodiesWithAllDegree2Terms,
                            bodyUndergoingTorqueName,
                            bodyExertingTorqueName,
                            fullDegreeTwoCombinations );
            // Full two-body spherical-harmonic torque restricted to point-mass/degree-2 terms only;
            // used to isolate the degree-2/degree-2 coupling by subtraction from the full model.
            std::shared_ptr< FullTwoBodySphericalHarmonicTorque > fullTwoBodyTorqueModelWithPointMassDegree2Terms =
                    createFactoryFullTwoBodySphericalHarmonicTorqueModel(
                            bodiesWithAllDegree2Terms,
                            bodyUndergoingTorqueName,
                            bodyExertingTorqueName,
                            pointMassDegreeTwoCombinations );
            // Full two-body spherical-harmonic torque with only direct degree-2/degree-2 cross interactions;
            // used as an explicit check that decomposition-by-subtraction recovers the same coupling term.
            std::shared_ptr< FullTwoBodySphericalHarmonicTorque > fullTwoBodyTorqueModelWithDegree2Degree2Terms =
                    createFactoryFullTwoBodySphericalHarmonicTorqueModel(
                            bodiesWithAllDegree2Terms,
                            bodyUndergoingTorqueName,
                            bodyExertingTorqueName,
                            degreeTwoDegreeTwoCombinations );
            // Fourth-degree two-body torque model (Eq. 11 equivalent) evaluated on the full degree-2 bodies;
            // used as an independent model to compare against the full two-body spherical-harmonic torque decomposition.
            std::shared_ptr< FourthDegreeFullTwoBodyGravitationalTorqueModel > fourthDegreeTorqueModelWithFullDegree2Terms =
                    createFactoryFourthDegreeFullTwoBodyGravitationalTorqueModel(
                            bodiesWithAllDegree2Terms, bodyUndergoingTorqueName, bodyExertingTorqueName );

            const double massOfBodyUndergoingTorque =
                    bodiesWithAllDegree2Terms.at( bodyUndergoingTorqueName )->getBodyMass( );
            const double massOfBodyExertingTorque =
                    bodiesWithAllDegree2Terms.at( bodyExertingTorqueName )->getBodyMass( );
            const Eigen::Matrix3d inertiaTensorOfBodyUndergoingTorque =
                    bodiesWithAllDegree2Terms.at( bodyUndergoingTorqueName )->getBodyInertiaTensor( );
            const SystemOfBodies fourthDegreeBodiesWithPointMassBody2 =
                    createSystemOfBodiesForFourthDegreeTwoBodyTorqueTest(
                            bodyUndergoingTorqueName,
                            bodyExertingTorqueName,
                            massOfBodyUndergoingTorque,
                            massOfBodyExertingTorque,
                            inertiaTensorOfBodyUndergoingTorque,
                            Eigen::Matrix3d::Zero( ),
                            distantPositionOfBody1,
                            distantPositionOfBody2,
                            orientationCase.first,
                            orientationCase.second );
            fourthDegreeBodiesWithPointMassBody2.at( bodyUndergoingTorqueName )->setCurrentRotationalStateToLocalFrameFromEphemeris(
                    evaluationTime );
            fourthDegreeBodiesWithPointMassBody2.at( bodyExertingTorqueName )->setCurrentRotationalStateToLocalFrameFromEphemeris(
                    evaluationTime );
            // Fourth-degree two-body torque model with a point-mass-equivalent exerting body (zero inertia tensor);
            // used to extract the point-mass/degree-2 contribution within the fourth-degree model for subtraction.
            std::shared_ptr< FourthDegreeFullTwoBodyGravitationalTorqueModel > fourthDegreeTorqueModelWithPointMassDegree2Terms =
                    createFactoryFourthDegreeFullTwoBodyGravitationalTorqueModel(
                            fourthDegreeBodiesWithPointMassBody2, bodyUndergoingTorqueName, bodyExertingTorqueName );

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

            // Convert full two-body model output from specific torque to torque by multiplying with body-2 mass.
            const double bodyExertingTorqueMass = bodiesWithAllDegree2Terms.at( bodyExertingTorqueName )->getBodyMass( );

            // Model type: [Full two-body]. Body exerting: [l=0,2]. Body undergoing: [l=2].
            const Eigen::Vector3d fullTwoBodyFullDegree2Torque =
                    -bodyExertingTorqueMass * fullTwoBodyTorqueModelWithFullDegree2Terms->getTorque( );

            // Model type: [Full two-body]. Body exerting: [l=0]. Body undergoing: [l=2].
            const Eigen::Vector3d fullTwoBodyPointMassDegree2Torque =
                    -bodyExertingTorqueMass * fullTwoBodyTorqueModelWithPointMassDegree2Terms->getTorque( );

            // Model type: [Full two-body]. Body exerting: [l=2]. Body undergoing: [l=2].
            const Eigen::Vector3d fullTwoBodyDegree2Degree2Torque =
                    -bodyExertingTorqueMass * fullTwoBodyTorqueModelWithDegree2Degree2Terms->getTorque( );

            // Model type: [Fourth degree model]. Body exerting: [l=0,2]. Body undergoing: [l=2].
            const Eigen::Vector3d fourthDegreeFullDegree2Torque = fourthDegreeTorqueModelWithFullDegree2Terms->getTorque( );

            // Model type: [Fourth degree model]. Body exerting: [l=0]. Body undergoing: [l=2].
            const Eigen::Vector3d fourthDegreePointMassDegree2Torque =
                    fourthDegreeTorqueModelWithPointMassDegree2Terms->getTorque( );

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
            const Eigen::Vector3d pointMassModelDifference =
                    fullTwoBodyPointMassDegree2Torque - fourthDegreePointMassDegree2Torque;
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

        // This check confirms the isolated degree-2/degree-2 coupling from the full two-body model varies with orientation.
        BOOST_CHECK_GT(
                ( isolatedDegree22TorquesFromFullTwoBodyModel.at( 0 ) - isolatedDegree22TorquesFromFullTwoBodyModel.at( 1 ) ).norm( ),
                1.0E-1 );

        // This check confirms the isolated degree-2/degree-2 coupling from the fourth-degree model varies with orientation.
        BOOST_CHECK_GT(
                ( isolatedDegree22TorquesFromFourthDegreeModel.at( 0 ) -
                  isolatedDegree22TorquesFromFourthDegreeModel.at( 1 ) )
                        .norm( ),
                1.0E-1 );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
