/*    Copyright (c) 2010-2026, Delft University of Technology
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
#include <map>

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"

#include "tudat/astro/ephemerides/simpleRotationalEphemeris.h"
#include "tudat/astro/gravitation/centralGravityModel.h"
#include "tudat/astro/gravitation/mutualSphericalHarmonicGravityModel.h"
#include "tudat/astro/gravitation/periodicGravityFieldVariations.h"
#include "tudat/astro/gravitation/polynomialGravityFieldVariations.h"
#include "tudat/astro/gravitation/sphericalHarmonicsGravityField.h"
#include "tudat/astro/gravitation/sphericalHarmonicsGravityModel.h"
#include "tudat/astro/orbit_determination/acceleration_partials/centralGravityAccelerationPartial.h"
#include "tudat/astro/orbit_determination/acceleration_partials/fullTwoBodySphericalHarmonicGravityPartial.h"
#include "tudat/astro/orbit_determination/acceleration_partials/mutualSphericalHarmonicGravityPartial.h"
#include "tudat/astro/orbit_determination/acceleration_partials/numericalAccelerationPartial.h"
#include "tudat/astro/orbit_determination/acceleration_partials/sphericalHarmonicAccelerationPartial.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/constantRotationalOrientation.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/constantRotationRate.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameterSet.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/gravityFieldVariationParameters.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/gravitationalParameter.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/sphericalHarmonicCosineCoefficients.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/sphericalHarmonicSineCoefficients.h"
#include "tudat/astro/orbit_determination/observation_partials/rotationMatrixPartial.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/estimation_setup/createAccelerationPartials.h"
#include "tudat/simulation/propagation_setup/accelerationSettings.h"
#include "tudat/simulation/propagation_setup/createAccelerationModels.h"

namespace tudat
{

namespace unit_tests
{

using namespace tudat::acceleration_partials;
using namespace tudat::basic_astrodynamics;
using namespace tudat::estimatable_parameters;
using namespace tudat::gravitation;
using namespace tudat::simulation_setup;

BOOST_AUTO_TEST_SUITE( test_mutual_extended_sh_acceleration_partials )

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

Eigen::Vector4d getSignConsistentQuaternionVector( const Eigen::Quaterniond& quaternion, const Eigen::Vector4d& referenceQuaternion )
{
    Eigen::Vector4d quaternionVector = linear_algebra::convertQuaternionToVectorFormat( quaternion );
    if( quaternionVector.dot( referenceQuaternion ) < 0.0 )
    {
        quaternionVector *= -1.0;
    }
    return quaternionVector;
}

BOOST_AUTO_TEST_CASE( testRotationMatrixPartialToQuaternionPartialChain )
{
    // Test rationale:
    // Verify the low-level chain rule that converts rotation-matrix parameter partials into quaternion
    // partials before it is used by full two-body orientation partials.
    const double rotationRate = 1.9E-5;
    const double initialTime = 300.0;
    const double testTime = 2500.0;
    const double rightAscension = 0.42;
    const double declination = 1.02;
    const double primeMeridian = -0.25;
    std::shared_ptr< ephemerides::SimpleRotationalEphemeris > rotationalEphemeris =
            std::make_shared< ephemerides::SimpleRotationalEphemeris >(
                    rightAscension, declination, primeMeridian, rotationRate, initialTime, "ECLIPJ2000", "Body_Fixed" );

    const Eigen::Vector4d nominalQuaternion =
            linear_algebra::convertQuaternionToVectorFormat( rotationalEphemeris->getRotationToBaseFrame( testTime ) );

    {
        std::shared_ptr< observation_partials::RotationMatrixPartialWrtConstantRotationRate > rotationMatrixPartialObject =
                std::make_shared< observation_partials::RotationMatrixPartialWrtConstantRotationRate >( rotationalEphemeris );
        const Eigen::MatrixXd analyticalQuaternionPartial =
                acceleration_partials::detail::computePartialOfQuaternionWrtRotationMatrixParameter(
                        rotationalEphemeris->getRotationToBaseFrame( testTime ),
                        rotationMatrixPartialObject->calculatePartialOfRotationMatrixToBaseFrameWrParameter( testTime ) );

        const double perturbation = 1.0E-8;
        rotationalEphemeris->resetRotationRate( rotationRate + perturbation );
        const Eigen::Vector4d upQuaternion =
                getSignConsistentQuaternionVector( rotationalEphemeris->getRotationToBaseFrame( testTime ), nominalQuaternion );
        rotationalEphemeris->resetRotationRate( rotationRate - perturbation );
        const Eigen::Vector4d downQuaternion =
                getSignConsistentQuaternionVector( rotationalEphemeris->getRotationToBaseFrame( testTime ), nominalQuaternion );
        rotationalEphemeris->resetRotationRate( rotationRate );

        const Eigen::MatrixXd numericalQuaternionPartial = ( upQuaternion - downQuaternion ) / ( 2.0 * perturbation );
        // Verify the rotation-rate quaternion chain rule against a central finite-difference reference.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( analyticalQuaternionPartial, numericalQuaternionPartial, 1.0E-7 );
    }

    {
        std::shared_ptr< observation_partials::RotationMatrixPartialWrtPoleOrientation > rotationMatrixPartialObject =
                std::make_shared< observation_partials::RotationMatrixPartialWrtPoleOrientation >( rotationalEphemeris );
        const Eigen::MatrixXd analyticalQuaternionPartial =
                acceleration_partials::detail::computePartialOfQuaternionWrtRotationMatrixParameter(
                        rotationalEphemeris->getRotationToBaseFrame( testTime ),
                        rotationMatrixPartialObject->calculatePartialOfRotationMatrixToBaseFrameWrParameter( testTime ) );

        const double perturbation = 1.0E-7;
        Eigen::MatrixXd numericalQuaternionPartial = Eigen::MatrixXd::Zero( 4, 2 );
        for( int i = 0; i < 2; i++ )
        {
            rotationalEphemeris->resetInitialPoleRightAscensionAndDeclination( rightAscension + ( i == 0 ? perturbation : 0.0 ),
                                                                               declination + ( i == 1 ? perturbation : 0.0 ) );
            const Eigen::Vector4d upQuaternion =
                    getSignConsistentQuaternionVector( rotationalEphemeris->getRotationToBaseFrame( testTime ), nominalQuaternion );
            rotationalEphemeris->resetInitialPoleRightAscensionAndDeclination( rightAscension - ( i == 0 ? perturbation : 0.0 ),
                                                                               declination - ( i == 1 ? perturbation : 0.0 ) );
            const Eigen::Vector4d downQuaternion =
                    getSignConsistentQuaternionVector( rotationalEphemeris->getRotationToBaseFrame( testTime ), nominalQuaternion );
            numericalQuaternionPartial.col( i ) = ( upQuaternion - downQuaternion ) / ( 2.0 * perturbation );
        }
        rotationalEphemeris->resetInitialPoleRightAscensionAndDeclination( rightAscension, declination );

        // Verify the pole-orientation quaternion chain rule against central finite differences.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( analyticalQuaternionPartial, numericalQuaternionPartial, 1.0E-7 );
    }
}

BOOST_AUTO_TEST_CASE( testFullTwoBodySphericalHarmonicGravityPartials )
{
    // Test rationale:
    // This test validates analytical partials of the full two-body spherical-harmonic acceleration model
    // (Dirkx et al., 2019 mutual-potential interaction summation) against central finite differences.
    //
    // We run three interaction subsets:
    // 1) "regular": body-1 harmonics acting on point-mass body 2,
    // 2) "mutualSinglePoint": mutual interactions where one side contributes only degree 0,
    // 3) "mutual": full degree-2 mutual interactions.
    // This structure checks both correctness of each subset and that enabling figure-figure terms changes
    // the partials in physically expected directions.
    struct PartialSet {
        Eigen::MatrixXd wrtBody1Position;
        Eigen::MatrixXd wrtBody2Position;
        Eigen::MatrixXd wrtBody1Velocity;
        Eigen::MatrixXd wrtBody2Velocity;
        Eigen::MatrixXd wrtBody1Cosine;
        Eigen::MatrixXd wrtBody1Sine;
        Eigen::MatrixXd wrtBody2Cosine;
        Eigen::MatrixXd wrtBody2Sine;
    };

    struct CaseDefinition {
        std::string name;
        std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > coefficientCombinations;
        bool hasBody2ShapeTerms;
        bool hasFigureFigureTerms;
    };

    const double gravitationalParameter = 5.0E5;
    const double equatorialRadiusOfBody1 = 1300.0;
    const double equatorialRadiusOfBody2 = 1100.0;

    Eigen::MatrixXd cosineCoefficientsOfBody1Base = Eigen::MatrixXd::Zero( 5, 5 );
    Eigen::MatrixXd sineCoefficientsOfBody1Base = Eigen::MatrixXd::Zero( 5, 5 );
    Eigen::MatrixXd cosineCoefficientsOfBody2Base = Eigen::MatrixXd::Zero( 5, 5 );
    Eigen::MatrixXd sineCoefficientsOfBody2Base = Eigen::MatrixXd::Zero( 5, 5 );

    cosineCoefficientsOfBody1Base( 0, 0 ) = 1.0;
    cosineCoefficientsOfBody2Base( 0, 0 ) = 1.0;

    cosineCoefficientsOfBody1Base( 2, 0 ) = 0.23;
    cosineCoefficientsOfBody1Base( 2, 1 ) = -0.11;
    sineCoefficientsOfBody1Base( 2, 1 ) = 0.14;
    cosineCoefficientsOfBody1Base( 2, 2 ) = 0.09;
    sineCoefficientsOfBody1Base( 2, 2 ) = -0.06;
    cosineCoefficientsOfBody1Base( 3, 0 ) = -0.034;
    cosineCoefficientsOfBody1Base( 3, 1 ) = 0.028;
    sineCoefficientsOfBody1Base( 3, 1 ) = -0.021;
    cosineCoefficientsOfBody1Base( 3, 3 ) = 0.017;
    sineCoefficientsOfBody1Base( 3, 3 ) = 0.013;
    cosineCoefficientsOfBody1Base( 4, 0 ) = 0.019;
    cosineCoefficientsOfBody1Base( 4, 2 ) = -0.014;
    sineCoefficientsOfBody1Base( 4, 2 ) = 0.011;
    cosineCoefficientsOfBody1Base( 4, 4 ) = 0.008;
    sineCoefficientsOfBody1Base( 4, 4 ) = -0.006;

    cosineCoefficientsOfBody2Base( 2, 0 ) = -0.19;
    cosineCoefficientsOfBody2Base( 2, 1 ) = 0.16;
    sineCoefficientsOfBody2Base( 2, 1 ) = -0.08;
    cosineCoefficientsOfBody2Base( 2, 2 ) = 0.13;
    sineCoefficientsOfBody2Base( 2, 2 ) = 0.12;
    cosineCoefficientsOfBody2Base( 3, 0 ) = 0.027;
    cosineCoefficientsOfBody2Base( 3, 1 ) = -0.023;
    sineCoefficientsOfBody2Base( 3, 1 ) = 0.018;
    cosineCoefficientsOfBody2Base( 3, 3 ) = -0.016;
    sineCoefficientsOfBody2Base( 3, 3 ) = 0.012;
    cosineCoefficientsOfBody2Base( 4, 0 ) = -0.017;
    cosineCoefficientsOfBody2Base( 4, 2 ) = 0.015;
    sineCoefficientsOfBody2Base( 4, 2 ) = -0.010;
    cosineCoefficientsOfBody2Base( 4, 4 ) = -0.007;
    sineCoefficientsOfBody2Base( 4, 4 ) = 0.005;

    std::shared_ptr< Body > body1 = std::make_shared< Body >( );
    std::shared_ptr< Body > body2 = std::make_shared< Body >( );

    std::shared_ptr< SphericalHarmonicsGravityField > body1GravityField = std::make_shared< SphericalHarmonicsGravityField >(
            gravitationalParameter, equatorialRadiusOfBody1, cosineCoefficientsOfBody1Base, sineCoefficientsOfBody1Base, "IAU_Body1" );
    std::shared_ptr< SphericalHarmonicsGravityField > body2GravityField = std::make_shared< SphericalHarmonicsGravityField >(
            gravitationalParameter, equatorialRadiusOfBody2, cosineCoefficientsOfBody2Base, sineCoefficientsOfBody2Base, "IAU_Body2" );

    body1->setGravityFieldModel( body1GravityField );
    body2->setGravityFieldModel( body2GravityField );

    const std::vector< std::pair< int, int > > cosineIndices = { { 2, 0 }, { 2, 1 }, { 2, 2 }, { 3, 0 }, { 3, 1 },
                                                                 { 3, 3 }, { 4, 0 }, { 4, 2 }, { 4, 4 } };
    const std::vector< std::pair< int, int > > sineIndices = { { 2, 1 }, { 2, 2 }, { 3, 1 }, { 3, 3 }, { 4, 2 }, { 4, 4 } };

    const std::vector< CaseDefinition > testCases = {
        { "regular", FullTwoBodySphericalHarmonicAccelerationSettings( 4, 4, 0, 0 ).coefficientCombinationsToUse_, false, false },
        { "mutualSinglePoint", getExtendedSinglePointMassInteractions( 4, 4, 4, 4 ), true, false },
        { "mutual", FullTwoBodySphericalHarmonicAccelerationSettings( 4, 4, 4, 4 ).coefficientCombinationsToUse_, true, true }
    };

    const Eigen::Vector3d positionPerturbation = Eigen::Vector3d::Constant( 10.0 );
    const Eigen::Vector3d velocityPerturbation = Eigen::Vector3d::Constant( 1.0 );

    for( unsigned int rotationCase = 0; rotationCase < 2; rotationCase++ )
    {
        Eigen::Vector6d stateOfBody1 = Eigen::Vector6d::Zero( );
        Eigen::Vector6d stateOfBody2 = Eigen::Vector6d::Zero( );
        Eigen::Quaterniond rotationToBody1 = Eigen::Quaterniond::Identity( );
        Eigen::Quaterniond rotationToBody2 = Eigen::Quaterniond::Identity( );

        if( rotationCase == 0 )
        {
            stateOfBody1.segment( 0, 3 ) = ( Eigen::Vector3d( ) << 5100.0, -2200.0, 3600.0 ).finished( );
            stateOfBody2.segment( 0, 3 ) = Eigen::Vector3d::Zero( );
        }
        else
        {
            rotationToBody1 = Eigen::Quaterniond( Eigen::AngleAxisd( 0.7, Eigen::Vector3d::UnitZ( ) ) *
                                                  Eigen::AngleAxisd( -0.4, Eigen::Vector3d::UnitX( ) ) *
                                                  Eigen::AngleAxisd( 0.2, Eigen::Vector3d::UnitY( ) ) );
            rotationToBody2 = Eigen::Quaterniond( Eigen::AngleAxisd( -0.3, Eigen::Vector3d::UnitZ( ) ) *
                                                  Eigen::AngleAxisd( 0.5, Eigen::Vector3d::UnitY( ) ) *
                                                  Eigen::AngleAxisd( 0.25, Eigen::Vector3d::UnitX( ) ) );
            stateOfBody1.segment( 0, 3 ) = ( Eigen::Vector3d( ) << -4300.0, 3100.0, 2800.0 ).finished( );
            stateOfBody2.segment( 0, 3 ) = Eigen::Vector3d::Zero( );
        }

        const double evaluationTime = 1000.0 + static_cast< double >( rotationCase ) * 100.0;

        body1->setRotationalEphemeris( std::make_shared< ephemerides::SimpleRotationalEphemeris >(
                rotationToBody1, 0.0, evaluationTime, "ECLIPJ2000", "IAU_Body1" ) );
        body2->setRotationalEphemeris( std::make_shared< ephemerides::SimpleRotationalEphemeris >(
                rotationToBody2, 0.0, evaluationTime, "ECLIPJ2000", "IAU_Body2" ) );

        body1->setCurrentRotationToLocalFrameFromEphemeris( evaluationTime );
        body2->setCurrentRotationToLocalFrameFromEphemeris( evaluationTime );
        body1->setState( stateOfBody1 );
        body2->setState( stateOfBody2 );

        std::map< std::string, PartialSet > analyticalPartialsByCase;

        for( const CaseDefinition& testCase : testCases )
        {
            body1GravityField->setCosineCoefficients( cosineCoefficientsOfBody1Base );
            body1GravityField->setSineCoefficients( sineCoefficientsOfBody1Base );
            body2GravityField->setCosineCoefficients( cosineCoefficientsOfBody2Base );
            body2GravityField->setSineCoefficients( sineCoefficientsOfBody2Base );

            std::function< void( Eigen::Vector6d ) > body1StateSetFunction = std::bind( &Body::setState, body1, std::placeholders::_1 );
            std::function< void( Eigen::Vector6d ) > body2StateSetFunction = std::bind( &Body::setState, body2, std::placeholders::_1 );

            std::shared_ptr< AccelerationSettings > accelerationSettings =
                    std::make_shared< FullTwoBodySphericalHarmonicAccelerationSettings >( testCase.coefficientCombinations );
            std::shared_ptr< FullTwoBodySphericalHarmonicAcceleration > accelerationModel =
                    std::dynamic_pointer_cast< FullTwoBodySphericalHarmonicAcceleration >(
                            createAccelerationModel( body1, body2, accelerationSettings, "Body1", "Body2" ) );

            // Require successful creation so the test exercises the intended full two-body acceleration path.
            BOOST_REQUIRE( accelerationModel != nullptr );

            std::shared_ptr< FullTwoBodySphericalHarmonicsGravityPartial > accelerationPartial =
                    std::make_shared< FullTwoBodySphericalHarmonicsGravityPartial >( "Body1",
                                                                                     "Body2",
                                                                                     accelerationModel,
                                                                                     createSimpleRotationPartialMap( body1 ),
                                                                                     createSimpleRotationPartialMap( body2 ) );

            std::shared_ptr< SphericalHarmonicsCosineCoefficients > body1CosineCoefficientsParameter =
                    std::make_shared< SphericalHarmonicsCosineCoefficients >(
                            std::bind( &SphericalHarmonicsGravityField::getCosineCoefficients, body1GravityField ),
                            std::bind( &SphericalHarmonicsGravityField::setCosineCoefficients, body1GravityField, std::placeholders::_1 ),
                            cosineIndices,
                            "Body1" );
            std::shared_ptr< SphericalHarmonicsSineCoefficients > body1SineCoefficientsParameter =
                    std::make_shared< SphericalHarmonicsSineCoefficients >(
                            std::bind( &SphericalHarmonicsGravityField::getSineCoefficients, body1GravityField ),
                            std::bind( &SphericalHarmonicsGravityField::setSineCoefficients, body1GravityField, std::placeholders::_1 ),
                            sineIndices,
                            "Body1" );
            std::shared_ptr< SphericalHarmonicsCosineCoefficients > body2CosineCoefficientsParameter =
                    std::make_shared< SphericalHarmonicsCosineCoefficients >(
                            std::bind( &SphericalHarmonicsGravityField::getCosineCoefficients, body2GravityField ),
                            std::bind( &SphericalHarmonicsGravityField::setCosineCoefficients, body2GravityField, std::placeholders::_1 ),
                            cosineIndices,
                            "Body2" );
            std::shared_ptr< SphericalHarmonicsSineCoefficients > body2SineCoefficientsParameter =
                    std::make_shared< SphericalHarmonicsSineCoefficients >(
                            std::bind( &SphericalHarmonicsGravityField::getSineCoefficients, body2GravityField ),
                            std::bind( &SphericalHarmonicsGravityField::setSineCoefficients, body2GravityField, std::placeholders::_1 ),
                            sineIndices,
                            "Body2" );
            std::shared_ptr< GravitationalParameter > body1GravitationalParameter =
                    std::make_shared< GravitationalParameter >( body1GravityField, "Body1" );
            std::shared_ptr< GravitationalParameter > body2GravitationalParameter =
                    std::make_shared< GravitationalParameter >( body2GravityField, "Body2" );

            body1->setCurrentRotationToLocalFrameFromEphemeris( evaluationTime );
            body2->setCurrentRotationToLocalFrameFromEphemeris( evaluationTime );
            accelerationModel->updateMembers( evaluationTime );
            accelerationPartial->update( evaluationTime );

            PartialSet analyticalPartials;
            analyticalPartials.wrtBody1Position = Eigen::MatrixXd::Zero( 3, 3 );
            analyticalPartials.wrtBody2Position = Eigen::MatrixXd::Zero( 3, 3 );
            analyticalPartials.wrtBody1Velocity = Eigen::MatrixXd::Zero( 3, 3 );
            analyticalPartials.wrtBody2Velocity = Eigen::MatrixXd::Zero( 3, 3 );

            accelerationPartial->wrtPositionOfAcceleratedBody( analyticalPartials.wrtBody1Position.block( 0, 0, 3, 3 ) );
            accelerationPartial->wrtPositionOfAcceleratingBody( analyticalPartials.wrtBody2Position.block( 0, 0, 3, 3 ) );
            accelerationPartial->wrtVelocityOfAcceleratedBody( analyticalPartials.wrtBody1Velocity.block( 0, 0, 3, 3 ) );
            accelerationPartial->wrtVelocityOfAcceleratingBody( analyticalPartials.wrtBody2Velocity.block( 0, 0, 3, 3 ) );

            analyticalPartials.wrtBody1Cosine = accelerationPartial->wrtParameter( body1CosineCoefficientsParameter );
            analyticalPartials.wrtBody1Sine = accelerationPartial->wrtParameter( body1SineCoefficientsParameter );
            analyticalPartials.wrtBody2Cosine = accelerationPartial->wrtParameter( body2CosineCoefficientsParameter );
            analyticalPartials.wrtBody2Sine = accelerationPartial->wrtParameter( body2SineCoefficientsParameter );
            const Eigen::Vector3d analyticalPartialWrtBody1GravitationalParameter =
                    accelerationPartial->wrtParameter( body1GravitationalParameter );
            const Eigen::Vector3d analyticalPartialWrtBody2GravitationalParameter =
                    accelerationPartial->wrtParameter( body2GravitationalParameter );

            if( testCase.hasFigureFigureTerms )
            {
                const double rotationParameterInitialTime = evaluationTime - 400.0;
                const double rotationParameterEvaluationTime = evaluationTime + 250.0;
                const double body1RotationRate = 2.0E-5;
                const double body2RotationRate = -1.5E-5;
                const double body1RightAscension = 0.35 + 0.08 * static_cast< double >( rotationCase );
                const double body1Declination = 0.95 - 0.04 * static_cast< double >( rotationCase );
                const double body1PrimeMeridian = -0.22 + 0.11 * static_cast< double >( rotationCase );
                const double body2RightAscension = -0.18 + 0.06 * static_cast< double >( rotationCase );
                const double body2Declination = 1.08 - 0.03 * static_cast< double >( rotationCase );
                const double body2PrimeMeridian = 0.31 - 0.09 * static_cast< double >( rotationCase );
                body1->setRotationalEphemeris( std::make_shared< ephemerides::SimpleRotationalEphemeris >( body1RightAscension,
                                                                                                           body1Declination,
                                                                                                           body1PrimeMeridian,
                                                                                                           body1RotationRate,
                                                                                                           rotationParameterInitialTime,
                                                                                                           "ECLIPJ2000",
                                                                                                           "IAU_Body1" ) );
                body2->setRotationalEphemeris( std::make_shared< ephemerides::SimpleRotationalEphemeris >( body2RightAscension,
                                                                                                           body2Declination,
                                                                                                           body2PrimeMeridian,
                                                                                                           body2RotationRate,
                                                                                                           rotationParameterInitialTime,
                                                                                                           "ECLIPJ2000",
                                                                                                           "IAU_Body2" ) );
                std::shared_ptr< FullTwoBodySphericalHarmonicsGravityPartial > rotationAccelerationPartial =
                        std::make_shared< FullTwoBodySphericalHarmonicsGravityPartial >( "Body1",
                                                                                         "Body2",
                                                                                         accelerationModel,
                                                                                         createSimpleRotationPartialMap( body1 ),
                                                                                         createSimpleRotationPartialMap( body2 ) );

                auto evaluateAccelerationForRotationParameter = [ & ]( ) {
                    body1->setCurrentRotationToLocalFrameFromEphemeris( rotationParameterEvaluationTime );
                    body2->setCurrentRotationToLocalFrameFromEphemeris( rotationParameterEvaluationTime );
                    accelerationModel->resetCurrentTime( );
                    rotationAccelerationPartial->resetCurrentTime( );
                    accelerationModel->updateMembers( rotationParameterEvaluationTime );
                    rotationAccelerationPartial->update( rotationParameterEvaluationTime );
                    return accelerationModel->getAcceleration( );
                };
                auto checkCentralQuaternionStatePartial = [ & ]( const std::shared_ptr< Body >& body, const std::string& bodyName ) {
                    evaluateAccelerationForRotationParameter( );
                    Eigen::MatrixXd analyticalPartialWrtRotationalState = Eigen::MatrixXd::Zero( 3, 7 );
                    rotationAccelerationPartial->wrtNonTranslationalStateOfAdditionalBody(
                            analyticalPartialWrtRotationalState.block( 0, 0, 3, 7 ),
                            std::make_pair( bodyName, "" ),
                            propagators::rotational_state );

                    std::vector< Eigen::Vector4d > appliedQuaternionPerturbations;
                    const Eigen::MatrixXd accelerationDeviation = calculateAccelerationDeviationDueToOrientationChange(
                            std::bind( &Body::setCurrentRotationalStateToLocalFrame, body, std::placeholders::_1 ),
                            accelerationModel,
                            body->getRotationalStateVector( ),
                            Eigen::Vector4d::Constant( 1.0E-7 ),
                            appliedQuaternionPerturbations,
                            [ & ]( ) { accelerationModel->resetCurrentTime( ); },
                            rotationParameterEvaluationTime );
                    accelerationModel->resetCurrentTime( );
                    evaluateAccelerationForRotationParameter( );

                    for( int index = 1; index < 4; index++ )
                    {
                        const Eigen::Vector3d analyticalAccelerationDeviation =
                                analyticalPartialWrtRotationalState.block( 0, 0, 3, 4 ) * appliedQuaternionPerturbations.at( index );
                        // Verify the quaternion-state partial against the shared orientation-deviation finite-difference helper.
                        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                                analyticalAccelerationDeviation, accelerationDeviation.col( index - 1 ), 1.0E-5 );
                    }
                };
                checkCentralQuaternionStatePartial( body1, "Body1" );
                checkCentralQuaternionStatePartial( body2, "Body2" );

                auto checkRotationRatePartial = [ & ]( const std::shared_ptr< RotationRate >& rotationRateParameter ) {
                    evaluateAccelerationForRotationParameter( );
                    const Eigen::Vector3d analyticalPartial = rotationAccelerationPartial->wrtParameter( rotationRateParameter );
                    const Eigen::Vector3d numericalPartial = calculateAccelerationWrtParameterPartials(
                            rotationRateParameter,
                            accelerationModel,
                            1.0E-8,
                            emptyFunction,
                            rotationParameterEvaluationTime,
                            [ & ]( const double currentTime ) {
                                body1->setCurrentRotationToLocalFrameFromEphemeris( currentTime );
                                body2->setCurrentRotationToLocalFrameFromEphemeris( currentTime );
                            } );
                    evaluateAccelerationForRotationParameter( );
                    // Verify the constant-rotation-rate partial against direct parameter finite differences.
                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( analyticalPartial, numericalPartial, 1.0E-5 );
                };
                auto checkPolePositionPartial = [ & ]( const std::shared_ptr< ConstantRotationalOrientation >& polePositionParameter ) {
                    evaluateAccelerationForRotationParameter( );
                    const Eigen::MatrixXd analyticalPartial = rotationAccelerationPartial->wrtParameter( polePositionParameter );
                    const Eigen::MatrixXd numericalPartial = calculateAccelerationWrtParameterPartials(
                            polePositionParameter,
                            accelerationModel,
                            Eigen::VectorXd::Constant( polePositionParameter->getParameterSize( ), 1.0E-7 ),
                            emptyFunction,
                            rotationParameterEvaluationTime,
                            [ & ]( const double currentTime ) {
                                body1->setCurrentRotationToLocalFrameFromEphemeris( currentTime );
                                body2->setCurrentRotationToLocalFrameFromEphemeris( currentTime );
                            } );
                    evaluateAccelerationForRotationParameter( );
                    // Verify the pole-position partial against direct parameter finite differences.
                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( analyticalPartial, numericalPartial, 1.0E-5 );
                };

                checkRotationRatePartial( std::make_shared< RotationRate >(
                        std::dynamic_pointer_cast< ephemerides::SimpleRotationalEphemeris >( body1->getRotationalEphemeris( ) ),
                        "Body1" ) );
                checkRotationRatePartial( std::make_shared< RotationRate >(
                        std::dynamic_pointer_cast< ephemerides::SimpleRotationalEphemeris >( body2->getRotationalEphemeris( ) ),
                        "Body2" ) );
                checkPolePositionPartial( std::make_shared< ConstantRotationalOrientation >(
                        std::dynamic_pointer_cast< ephemerides::SimpleRotationalEphemeris >( body1->getRotationalEphemeris( ) ),
                        "Body1" ) );
                checkPolePositionPartial( std::make_shared< ConstantRotationalOrientation >(
                        std::dynamic_pointer_cast< ephemerides::SimpleRotationalEphemeris >( body2->getRotationalEphemeris( ) ),
                        "Body2" ) );

                body1->setRotationalEphemeris( std::make_shared< ephemerides::SimpleRotationalEphemeris >(
                        rotationToBody1, 0.0, evaluationTime, "ECLIPJ2000", "IAU_Body1" ) );
                body2->setRotationalEphemeris( std::make_shared< ephemerides::SimpleRotationalEphemeris >(
                        rotationToBody2, 0.0, evaluationTime, "ECLIPJ2000", "IAU_Body2" ) );
                body1->setCurrentRotationToLocalFrameFromEphemeris( evaluationTime );
                body2->setCurrentRotationToLocalFrameFromEphemeris( evaluationTime );
                accelerationModel->resetCurrentTime( );
                accelerationPartial->resetCurrentTime( );
                accelerationModel->updateMembers( evaluationTime );
                accelerationPartial->update( evaluationTime );
            }

            const std::function< void( ) > updateAccelerationFiniteDifferenceEnvironment = [ & ]( ) {
                body1->setCurrentRotationToLocalFrameFromEphemeris( evaluationTime );
                body2->setCurrentRotationToLocalFrameFromEphemeris( evaluationTime );
            };

            const Eigen::Matrix3d numericalPartialWrtBody1Position =
                    calculateAccelerationWrtStatePartials( body1StateSetFunction,
                                                           accelerationModel,
                                                           body1->getState( ),
                                                           positionPerturbation,
                                                           0,
                                                           updateAccelerationFiniteDifferenceEnvironment,
                                                           evaluationTime );
            const Eigen::Matrix3d numericalPartialWrtBody1Velocity =
                    calculateAccelerationWrtStatePartials( body1StateSetFunction,
                                                           accelerationModel,
                                                           body1->getState( ),
                                                           velocityPerturbation,
                                                           3,
                                                           updateAccelerationFiniteDifferenceEnvironment,
                                                           evaluationTime );
            const Eigen::Matrix3d numericalPartialWrtBody2Position =
                    calculateAccelerationWrtStatePartials( body2StateSetFunction,
                                                           accelerationModel,
                                                           body2->getState( ),
                                                           positionPerturbation,
                                                           0,
                                                           updateAccelerationFiniteDifferenceEnvironment,
                                                           evaluationTime );
            const Eigen::Matrix3d numericalPartialWrtBody2Velocity =
                    calculateAccelerationWrtStatePartials( body2StateSetFunction,
                                                           accelerationModel,
                                                           body2->getState( ),
                                                           velocityPerturbation,
                                                           3,
                                                           updateAccelerationFiniteDifferenceEnvironment,
                                                           evaluationTime );

            const Eigen::MatrixXd numericalPartialWrtBody1Cosine = calculateAccelerationWrtParameterPartials(
                    body1CosineCoefficientsParameter,
                    accelerationModel,
                    Eigen::VectorXd::Constant( body1CosineCoefficientsParameter->getParameterSize( ), 1.0E-2 ),
                    updateAccelerationFiniteDifferenceEnvironment,
                    evaluationTime );
            const Eigen::MatrixXd numericalPartialWrtBody1Sine = calculateAccelerationWrtParameterPartials(
                    body1SineCoefficientsParameter,
                    accelerationModel,
                    Eigen::VectorXd::Constant( body1SineCoefficientsParameter->getParameterSize( ), 1.0E-2 ),
                    updateAccelerationFiniteDifferenceEnvironment,
                    evaluationTime );
            const Eigen::MatrixXd numericalPartialWrtBody2Cosine = calculateAccelerationWrtParameterPartials(
                    body2CosineCoefficientsParameter,
                    accelerationModel,
                    Eigen::VectorXd::Constant( body2CosineCoefficientsParameter->getParameterSize( ), 1.0E-2 ),
                    updateAccelerationFiniteDifferenceEnvironment,
                    evaluationTime );
            const Eigen::MatrixXd numericalPartialWrtBody2Sine = calculateAccelerationWrtParameterPartials(
                    body2SineCoefficientsParameter,
                    accelerationModel,
                    Eigen::VectorXd::Constant( body2SineCoefficientsParameter->getParameterSize( ), 1.0E-2 ),
                    updateAccelerationFiniteDifferenceEnvironment,
                    evaluationTime );
            const Eigen::Vector3d numericalPartialWrtBody1GravitationalParameter = calculateAccelerationWrtParameterPartials(
                    body1GravitationalParameter, accelerationModel, 10.0, updateAccelerationFiniteDifferenceEnvironment, evaluationTime );
            const Eigen::Vector3d numericalPartialWrtBody2GravitationalParameter = calculateAccelerationWrtParameterPartials(
                    body2GravitationalParameter, accelerationModel, 10.0, updateAccelerationFiniteDifferenceEnvironment, evaluationTime );

            // Verify body-1 position partial against central finite differences.
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( numericalPartialWrtBody1Position, analyticalPartials.wrtBody1Position, 5.0E-5 );
            // Verify body-2 position partial against central finite differences.
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( numericalPartialWrtBody2Position, analyticalPartials.wrtBody2Position, 5.0E-5 );
            // Verify body-1 gravitational parameter has no direct response in this acceleration convention.
            BOOST_CHECK_SMALL( analyticalPartialWrtBody1GravitationalParameter.norm( ), 1.0E-20 );
            // Verify finite differencing also gives zero body-1 gravitational-parameter response.
            BOOST_CHECK_SMALL( numericalPartialWrtBody1GravitationalParameter.norm( ), 1.0E-20 );
            // Verify body-2 gravitational-parameter partial against central finite differences.
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    numericalPartialWrtBody2GravitationalParameter, analyticalPartialWrtBody2GravitationalParameter, 1.0E-10 );

            // Verify there is no analytical dependence on body-1 velocity.
            BOOST_CHECK_SMALL( analyticalPartials.wrtBody1Velocity.norm( ), std::numeric_limits< double >::epsilon( ) );
            // Verify there is no analytical dependence on body-2 velocity.
            BOOST_CHECK_SMALL( analyticalPartials.wrtBody2Velocity.norm( ), std::numeric_limits< double >::epsilon( ) );
            // Verify finite differencing also sees no body-1 velocity dependence.
            BOOST_CHECK_SMALL( numericalPartialWrtBody1Velocity.norm( ), std::numeric_limits< double >::epsilon( ) );
            // Verify finite differencing also sees no body-2 velocity dependence.
            BOOST_CHECK_SMALL( numericalPartialWrtBody2Velocity.norm( ), std::numeric_limits< double >::epsilon( ) );

            // Verify body-1 cosine coefficient partials against central finite differences.
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( numericalPartialWrtBody1Cosine, analyticalPartials.wrtBody1Cosine, 1.0E-8 );
            // Verify body-1 sine coefficient partials against central finite differences.
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( numericalPartialWrtBody1Sine, analyticalPartials.wrtBody1Sine, 1.0E-8 );

            if( testCase.hasBody2ShapeTerms )
            {
                // Verify active body-2 cosine coefficient partials against central finite differences.
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( numericalPartialWrtBody2Cosine, analyticalPartials.wrtBody2Cosine, 1.0E-8 );
                // Verify active body-2 sine coefficient partials against central finite differences.
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( numericalPartialWrtBody2Sine, analyticalPartials.wrtBody2Sine, 1.0E-8 );
            }
            else
            {
                // Verify inactive body-2 cosine coefficients produce no analytical partial.
                BOOST_CHECK_SMALL( analyticalPartials.wrtBody2Cosine.norm( ), 1.0E-20 );
                // Verify inactive body-2 sine coefficients produce no analytical partial.
                BOOST_CHECK_SMALL( analyticalPartials.wrtBody2Sine.norm( ), 1.0E-20 );
                // Verify finite differencing confirms zero response to inactive body-2 cosine coefficients.
                BOOST_CHECK_SMALL( numericalPartialWrtBody2Cosine.norm( ), 1.0E-20 );
                // Verify finite differencing confirms zero response to inactive body-2 sine coefficients.
                BOOST_CHECK_SMALL( numericalPartialWrtBody2Sine.norm( ), 1.0E-20 );
            }

            if( testCase.hasFigureFigureTerms && rotationCase == 1 )
            {
                const double variationEvaluationTime = evaluationTime + 37.0;
                const Eigen::MatrixXd zeroVariationBlock = Eigen::MatrixXd::Zero( 3, 5 );
                const std::map< int, std::vector< std::pair< int, int > > > polynomialCosineIndices = { { 1, { { 3, 1 }, { 4, 2 } } },
                                                                                                        { 2, { { 2, 0 } } } };
                const std::map< int, std::vector< std::pair< int, int > > > polynomialSineIndices = { { 1, { { 3, 1 }, { 4, 4 } } } };
                const std::map< int, std::vector< std::pair< int, int > > > periodicCosineIndices = { { 0, { { 3, 0 }, { 4, 4 } } },
                                                                                                      { 1, { { 4, 2 } } } };
                const std::map< int, std::vector< std::pair< int, int > > > periodicSineIndices = { { 0, { { 3, 1 } } },
                                                                                                    { 1, { { 4, 4 } } } };
                const std::vector< double > variationFrequencies = { 1.0E-3, 1.7E-3 };

                auto resetGravityFields = [ & ]( ) {
                    body1GravityField->setCosineCoefficients( cosineCoefficientsOfBody1Base );
                    body1GravityField->setSineCoefficients( sineCoefficientsOfBody1Base );
                    body2GravityField->setCosineCoefficients( cosineCoefficientsOfBody2Base );
                    body2GravityField->setSineCoefficients( sineCoefficientsOfBody2Base );
                };

                auto evaluateAccelerationForVariationParameter = [ & ]( ) {
                    body1->setCurrentRotationToLocalFrameFromEphemeris( variationEvaluationTime );
                    body2->setCurrentRotationToLocalFrameFromEphemeris( variationEvaluationTime );
                    accelerationModel->resetCurrentTime( );
                    accelerationPartial->resetCurrentTime( );
                    accelerationModel->updateMembers( variationEvaluationTime );
                    accelerationPartial->update( variationEvaluationTime );
                    return accelerationModel->getAcceleration( );
                };

                auto checkGravityFieldVariationPartial = [ & ]( const std::shared_ptr< EstimatableParameter< Eigen::VectorXd > >& parameter,
                                                                const std::shared_ptr< GravityFieldVariations >& variationModel,
                                                                const std::shared_ptr< SphericalHarmonicsGravityField >& gravityField,
                                                                const Eigen::MatrixXd& baseCosineCoefficients,
                                                                const Eigen::MatrixXd& baseSineCoefficients ) {
                    auto applyVariation = [ & ]( ) {
                        Eigen::MatrixXd cosineCoefficients = baseCosineCoefficients;
                        Eigen::MatrixXd sineCoefficients = baseSineCoefficients;
                        variationModel->addSphericalHarmonicsCorrections( variationEvaluationTime, sineCoefficients, cosineCoefficients );
                        gravityField->setCosineCoefficients( cosineCoefficients );
                        gravityField->setSineCoefficients( sineCoefficients );
                    };
                    const std::function< void( ) > updateVariationEnvironment = [ & ]( ) {
                        resetGravityFields( );
                        applyVariation( );
                    };

                    updateVariationEnvironment( );
                    evaluateAccelerationForVariationParameter( );
                    const Eigen::MatrixXd analyticalPartial = accelerationPartial->wrtParameter( parameter );

                    const Eigen::VectorXd nominalParameter = parameter->getParameterValue( );
                    const Eigen::MatrixXd numericalPartial = calculateAccelerationWrtParameterPartials(
                            parameter,
                            accelerationModel,
                            Eigen::VectorXd::Constant( nominalParameter.size( ), 1.0E-5 ),
                            updateVariationEnvironment,
                            variationEvaluationTime,
                            [ & ]( const double currentTime ) {
                                body1->setCurrentRotationToLocalFrameFromEphemeris( currentTime );
                                body2->setCurrentRotationToLocalFrameFromEphemeris( currentTime );
                            } );
                    updateVariationEnvironment( );
                    evaluateAccelerationForVariationParameter( );

                    // Verify the gravity-field variation chain rule against direct variation-amplitude finite differences.
                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( analyticalPartial, numericalPartial, 3.0E-3 );
                    resetGravityFields( );
                };

                auto createPolynomialVariationParameter =
                        [ & ]( const std::string& bodyName ) -> std::shared_ptr< PolynomialGravityFieldVariationsParameters > {
                    std::shared_ptr< PolynomialGravityFieldVariations > variationModel =
                            std::make_shared< PolynomialGravityFieldVariations >(
                                    std::map< int, Eigen::MatrixXd >{ { 1, zeroVariationBlock }, { 2, zeroVariationBlock } },
                                    std::map< int, Eigen::MatrixXd >{ { 1, zeroVariationBlock } },
                                    variationEvaluationTime - 0.25,
                                    2,
                                    0 );
                    return std::make_shared< PolynomialGravityFieldVariationsParameters >(
                            variationModel, polynomialCosineIndices, polynomialSineIndices, bodyName );
                };
                auto createPeriodicVariationParameter =
                        [ & ]( const std::string& bodyName ) -> std::shared_ptr< PeriodicGravityFieldVariationsParameters > {
                    const std::vector< Eigen::MatrixXd > zeroPeriodicBlocks( 2, zeroVariationBlock );
                    std::shared_ptr< PeriodicGravityFieldVariations > variationModel =
                            std::make_shared< PeriodicGravityFieldVariations >( zeroPeriodicBlocks,
                                                                                zeroPeriodicBlocks,
                                                                                zeroPeriodicBlocks,
                                                                                zeroPeriodicBlocks,
                                                                                variationFrequencies,
                                                                                variationEvaluationTime - 0.4,
                                                                                2,
                                                                                0 );
                    return std::make_shared< PeriodicGravityFieldVariationsParameters >(
                            variationModel, periodicCosineIndices, periodicSineIndices, bodyName );
                };

                std::shared_ptr< PolynomialGravityFieldVariationsParameters > body1PolynomialParameter =
                        createPolynomialVariationParameter( "Body1" );
                checkGravityFieldVariationPartial( body1PolynomialParameter,
                                                   body1PolynomialParameter->getPolynomialVariationModel( ),
                                                   body1GravityField,
                                                   cosineCoefficientsOfBody1Base,
                                                   sineCoefficientsOfBody1Base );
                std::shared_ptr< PeriodicGravityFieldVariationsParameters > body1PeriodicParameter =
                        createPeriodicVariationParameter( "Body1" );
                checkGravityFieldVariationPartial( body1PeriodicParameter,
                                                   body1PeriodicParameter->getPeriodicVariationModel( ),
                                                   body1GravityField,
                                                   cosineCoefficientsOfBody1Base,
                                                   sineCoefficientsOfBody1Base );
                std::shared_ptr< PolynomialGravityFieldVariationsParameters > body2PolynomialParameter =
                        createPolynomialVariationParameter( "Body2" );
                checkGravityFieldVariationPartial( body2PolynomialParameter,
                                                   body2PolynomialParameter->getPolynomialVariationModel( ),
                                                   body2GravityField,
                                                   cosineCoefficientsOfBody2Base,
                                                   sineCoefficientsOfBody2Base );
                std::shared_ptr< PeriodicGravityFieldVariationsParameters > body2PeriodicParameter =
                        createPeriodicVariationParameter( "Body2" );
                checkGravityFieldVariationPartial( body2PeriodicParameter,
                                                   body2PeriodicParameter->getPeriodicVariationModel( ),
                                                   body2GravityField,
                                                   cosineCoefficientsOfBody2Base,
                                                   sineCoefficientsOfBody2Base );
            }

            analyticalPartialsByCase[ testCase.name ] = analyticalPartials;
        }

        // Verify that adding mutual figure-figure terms changes the position Jacobian relative to the regular model.
        BOOST_CHECK_GT(
                ( analyticalPartialsByCase.at( "mutual" ).wrtBody1Position - analyticalPartialsByCase.at( "regular" ).wrtBody1Position )
                        .norm( ),
                1.0E-7 );
        // Verify that even the single-point mutual model is distinct from the one-way regular model.
        BOOST_CHECK_GT( ( analyticalPartialsByCase.at( "mutualSinglePoint" ).wrtBody1Position -
                          analyticalPartialsByCase.at( "regular" ).wrtBody1Position )
                                .norm( ),
                        1.0E-7 );
        // Verify figure-figure coupling changes body-1 coefficient partials beyond the single-point mutual subset.
        BOOST_CHECK_GT( ( analyticalPartialsByCase.at( "mutual" ).wrtBody1Cosine -
                          analyticalPartialsByCase.at( "mutualSinglePoint" ).wrtBody1Cosine )
                                .norm( ),
                        1.0E-14 );
        // Verify body-2 cosine coefficients are active only in the full mutual model.
        BOOST_CHECK_GT( analyticalPartialsByCase.at( "mutual" ).wrtBody2Cosine.norm( ), 1.0E-8 );
        // Verify body-2 sine coefficients are active only in the full mutual model.
        BOOST_CHECK_GT( analyticalPartialsByCase.at( "mutual" ).wrtBody2Sine.norm( ), 1.0E-8 );
    }
}

BOOST_AUTO_TEST_CASE( testFullTwoBodySphericalHarmonicPartialsAgainstEquivalentSimplerModels )
{
    // Test rationale:
    // Verify that full-two-body partials reduce exactly to simpler known models when interaction sets are
    // constrained. This is a model-equivalence test (not only analytical-vs-numerical):
    // 1) point-mass <-> central gravity,
    // 2) extended body-1 on point-mass body-2 <-> one-way spherical harmonics,
    // 3) point-mass body-1 on extended body-2 <-> opposite one-way spherical harmonics,
    // 4) mutual single-point interactions <-> legacy mutual spherical-harmonic model.
    //
    // The equivalences follow directly from selecting subsets of the Dirkx interaction sum over coefficient pairs.
    enum EquivalentModelType {
        pointMassEquivalent = 0,
        sphericalBody1OnBody2Equivalent = 1,
        sphericalBody2OnBody1Equivalent = 2,
        mutualSphericalEquivalent = 3
    };

    struct EquivalentCaseDefinition {
        std::string name;
        std::shared_ptr< AccelerationSettings > mutualExtendedSettings;
        EquivalentModelType equivalentModelType;
        bool compareBody1Coefficients;
        bool compareBody2Coefficients;
    };

    const double gravitationalParameter = 5.0E5;
    const double equatorialRadiusOfBody1 = 1300.0;
    const double equatorialRadiusOfBody2 = 1100.0;

    Eigen::MatrixXd cosineCoefficientsOfBody1Base = Eigen::MatrixXd::Zero( 3, 3 );
    Eigen::MatrixXd sineCoefficientsOfBody1Base = Eigen::MatrixXd::Zero( 3, 3 );
    Eigen::MatrixXd cosineCoefficientsOfBody2Base = Eigen::MatrixXd::Zero( 3, 3 );
    Eigen::MatrixXd sineCoefficientsOfBody2Base = Eigen::MatrixXd::Zero( 3, 3 );

    cosineCoefficientsOfBody1Base( 0, 0 ) = 1.0;
    cosineCoefficientsOfBody2Base( 0, 0 ) = 1.0;

    cosineCoefficientsOfBody1Base( 2, 0 ) = 0.23;
    cosineCoefficientsOfBody1Base( 2, 1 ) = -0.11;
    sineCoefficientsOfBody1Base( 2, 1 ) = 0.14;
    cosineCoefficientsOfBody1Base( 2, 2 ) = 0.09;
    sineCoefficientsOfBody1Base( 2, 2 ) = -0.06;

    cosineCoefficientsOfBody2Base( 2, 0 ) = -0.19;
    cosineCoefficientsOfBody2Base( 2, 1 ) = 0.16;
    sineCoefficientsOfBody2Base( 2, 1 ) = -0.08;
    cosineCoefficientsOfBody2Base( 2, 2 ) = 0.13;
    sineCoefficientsOfBody2Base( 2, 2 ) = 0.12;

    std::shared_ptr< Body > body1 = std::make_shared< Body >( );
    std::shared_ptr< Body > body2 = std::make_shared< Body >( );

    std::shared_ptr< SphericalHarmonicsGravityField > body1GravityField = std::make_shared< SphericalHarmonicsGravityField >(
            gravitationalParameter, equatorialRadiusOfBody1, cosineCoefficientsOfBody1Base, sineCoefficientsOfBody1Base, "IAU_Body1" );
    std::shared_ptr< SphericalHarmonicsGravityField > body2GravityField = std::make_shared< SphericalHarmonicsGravityField >(
            gravitationalParameter, equatorialRadiusOfBody2, cosineCoefficientsOfBody2Base, sineCoefficientsOfBody2Base, "IAU_Body2" );

    body1->setGravityFieldModel( body1GravityField );
    body2->setGravityFieldModel( body2GravityField );

    const std::vector< std::pair< int, int > > cosineIndices = { { 2, 0 }, { 2, 1 }, { 2, 2 } };
    const std::vector< std::pair< int, int > > sineIndices = { { 2, 1 }, { 2, 2 } };

    const std::vector< EquivalentCaseDefinition > testCases = {
        { "pointMassOnPointMass",
          std::make_shared< FullTwoBodySphericalHarmonicAccelerationSettings >( 0, 0, 0, 0 ),
          pointMassEquivalent,
          false,
          false },
        { "extendedBody1OnPointMass",
          std::make_shared< FullTwoBodySphericalHarmonicAccelerationSettings >( 2, 2, 0, 0 ),
          sphericalBody1OnBody2Equivalent,
          true,
          false },
        { "pointMassOnExtendedBody2",
          std::make_shared< FullTwoBodySphericalHarmonicAccelerationSettings >( 0, 0, 2, 2 ),
          sphericalBody2OnBody1Equivalent,
          false,
          true },
        { "mutualSinglePointInteractions",
          std::make_shared< FullTwoBodySphericalHarmonicAccelerationSettings >( getExtendedSinglePointMassInteractions( 2, 2, 2, 2 ) ),
          mutualSphericalEquivalent,
          true,
          true }
    };

    for( unsigned int rotationCase = 0; rotationCase < 2; rotationCase++ )
    {
        Eigen::Vector6d stateOfBody1 = Eigen::Vector6d::Zero( );
        Eigen::Vector6d stateOfBody2 = Eigen::Vector6d::Zero( );
        Eigen::Quaterniond rotationToBody1 = Eigen::Quaterniond::Identity( );
        Eigen::Quaterniond rotationToBody2 = Eigen::Quaterniond::Identity( );

        if( rotationCase == 0 )
        {
            stateOfBody1.segment( 0, 3 ) = ( Eigen::Vector3d( ) << 5100.0, -2200.0, 3600.0 ).finished( );
            stateOfBody2.segment( 0, 3 ) = Eigen::Vector3d::Zero( );
        }
        else
        {
            rotationToBody1 = Eigen::Quaterniond( Eigen::AngleAxisd( 0.7, Eigen::Vector3d::UnitZ( ) ) *
                                                  Eigen::AngleAxisd( -0.4, Eigen::Vector3d::UnitX( ) ) *
                                                  Eigen::AngleAxisd( 0.2, Eigen::Vector3d::UnitY( ) ) );
            rotationToBody2 = Eigen::Quaterniond( Eigen::AngleAxisd( -0.3, Eigen::Vector3d::UnitZ( ) ) *
                                                  Eigen::AngleAxisd( 0.5, Eigen::Vector3d::UnitY( ) ) *
                                                  Eigen::AngleAxisd( 0.25, Eigen::Vector3d::UnitX( ) ) );
            stateOfBody1.segment( 0, 3 ) = ( Eigen::Vector3d( ) << -4300.0, 3100.0, 2800.0 ).finished( );
            stateOfBody2.segment( 0, 3 ) = Eigen::Vector3d::Zero( );
        }

        const double evaluationTime = 1000.0 + static_cast< double >( rotationCase ) * 100.0;

        body1->setRotationalEphemeris( std::make_shared< ephemerides::SimpleRotationalEphemeris >(
                rotationToBody1, 0.0, evaluationTime, "ECLIPJ2000", "IAU_Body1" ) );
        body2->setRotationalEphemeris( std::make_shared< ephemerides::SimpleRotationalEphemeris >(
                rotationToBody2, 0.0, evaluationTime, "ECLIPJ2000", "IAU_Body2" ) );

        body1->setCurrentRotationToLocalFrameFromEphemeris( evaluationTime );
        body2->setCurrentRotationToLocalFrameFromEphemeris( evaluationTime );
        body1->setState( stateOfBody1 );
        body2->setState( stateOfBody2 );

        for( const EquivalentCaseDefinition& testCase : testCases )
        {
            body1GravityField->setCosineCoefficients( cosineCoefficientsOfBody1Base );
            body1GravityField->setSineCoefficients( sineCoefficientsOfBody1Base );
            body2GravityField->setCosineCoefficients( cosineCoefficientsOfBody2Base );
            body2GravityField->setSineCoefficients( sineCoefficientsOfBody2Base );

            std::shared_ptr< SphericalHarmonicsCosineCoefficients > body1CosineCoefficientsParameter =
                    std::make_shared< SphericalHarmonicsCosineCoefficients >(
                            std::bind( &SphericalHarmonicsGravityField::getCosineCoefficients, body1GravityField ),
                            std::bind( &SphericalHarmonicsGravityField::setCosineCoefficients, body1GravityField, std::placeholders::_1 ),
                            cosineIndices,
                            "Body1" );
            std::shared_ptr< SphericalHarmonicsSineCoefficients > body1SineCoefficientsParameter =
                    std::make_shared< SphericalHarmonicsSineCoefficients >(
                            std::bind( &SphericalHarmonicsGravityField::getSineCoefficients, body1GravityField ),
                            std::bind( &SphericalHarmonicsGravityField::setSineCoefficients, body1GravityField, std::placeholders::_1 ),
                            sineIndices,
                            "Body1" );
            std::shared_ptr< SphericalHarmonicsCosineCoefficients > body2CosineCoefficientsParameter =
                    std::make_shared< SphericalHarmonicsCosineCoefficients >(
                            std::bind( &SphericalHarmonicsGravityField::getCosineCoefficients, body2GravityField ),
                            std::bind( &SphericalHarmonicsGravityField::setCosineCoefficients, body2GravityField, std::placeholders::_1 ),
                            cosineIndices,
                            "Body2" );
            std::shared_ptr< SphericalHarmonicsSineCoefficients > body2SineCoefficientsParameter =
                    std::make_shared< SphericalHarmonicsSineCoefficients >(
                            std::bind( &SphericalHarmonicsGravityField::getSineCoefficients, body2GravityField ),
                            std::bind( &SphericalHarmonicsGravityField::setSineCoefficients, body2GravityField, std::placeholders::_1 ),
                            sineIndices,
                            "Body2" );

            std::shared_ptr< FullTwoBodySphericalHarmonicAcceleration > mutualExtendedModel =
                    std::dynamic_pointer_cast< FullTwoBodySphericalHarmonicAcceleration >(
                            createAccelerationModel( body1, body2, testCase.mutualExtendedSettings, "Body1", "Body2" ) );
            // Require successful creation of the full-two-body model for this equivalence case.
            BOOST_REQUIRE( mutualExtendedModel != nullptr );

            std::shared_ptr< FullTwoBodySphericalHarmonicsGravityPartial > mutualExtendedPartial =
                    std::make_shared< FullTwoBodySphericalHarmonicsGravityPartial >( "Body1", "Body2", mutualExtendedModel );

            mutualExtendedModel->updateMembers( evaluationTime );
            mutualExtendedPartial->update( evaluationTime );

            Eigen::MatrixXd mutualPartialWrtBody1Position = Eigen::MatrixXd::Zero( 3, 3 );
            Eigen::MatrixXd mutualPartialWrtBody2Position = Eigen::MatrixXd::Zero( 3, 3 );
            mutualExtendedPartial->wrtPositionOfAcceleratedBody( mutualPartialWrtBody1Position.block( 0, 0, 3, 3 ) );
            mutualExtendedPartial->wrtPositionOfAcceleratingBody( mutualPartialWrtBody2Position.block( 0, 0, 3, 3 ) );
            Eigen::MatrixXd mutualPartialWrtBody1Cosine = mutualExtendedPartial->wrtParameter( body1CosineCoefficientsParameter );
            Eigen::MatrixXd mutualPartialWrtBody1Sine = mutualExtendedPartial->wrtParameter( body1SineCoefficientsParameter );
            Eigen::MatrixXd mutualPartialWrtBody2Cosine = mutualExtendedPartial->wrtParameter( body2CosineCoefficientsParameter );
            Eigen::MatrixXd mutualPartialWrtBody2Sine = mutualExtendedPartial->wrtParameter( body2SineCoefficientsParameter );

            std::shared_ptr< AccelerationModel< Eigen::Vector3d > > equivalentModel;
            std::shared_ptr< AccelerationPartial > equivalentPartial;
            double scaleToMutual = 1.0;

            if( testCase.equivalentModelType == pointMassEquivalent )
            {
                std::shared_ptr< CentralGravitationalAccelerationModel3d > pointMassModel =
                        std::dynamic_pointer_cast< CentralGravitationalAccelerationModel3d >( createAccelerationModel(
                                body1, body2, std::make_shared< AccelerationSettings >( point_mass_gravity ), "Body1", "Body2" ) );
                // Require the equivalent point-mass model needed for the reduction check.
                BOOST_REQUIRE( pointMassModel != nullptr );
                equivalentModel = pointMassModel;
                equivalentPartial = std::make_shared< CentralGravitationPartial >( pointMassModel, "Body1", "Body2" );
            }
            else if( testCase.equivalentModelType == sphericalBody1OnBody2Equivalent )
            {
                std::shared_ptr< SphericalHarmonicsGravitationalAccelerationModel > sphericalHarmonicModel =
                        std::dynamic_pointer_cast< SphericalHarmonicsGravitationalAccelerationModel >( createAccelerationModel(
                                body2, body1, std::make_shared< SphericalHarmonicAccelerationSettings >( 2, 2 ), "Body2", "Body1" ) );
                // Require the equivalent one-way spherical-harmonic model for body-1 shape terms.
                BOOST_REQUIRE( sphericalHarmonicModel != nullptr );
                equivalentModel = sphericalHarmonicModel;
                equivalentPartial = std::make_shared< SphericalHarmonicsGravityPartial >( "Body2", "Body1", sphericalHarmonicModel );
                // Sign/scale conversion: this equivalent model computes acceleration of body 2 due to body 1.
                // The mutual model here is defined for body 1 due to body 2, hence the factor below.
                scaleToMutual = -body2GravityField->getGravitationalParameter( ) / body1GravityField->getGravitationalParameter( );
            }
            else if( testCase.equivalentModelType == sphericalBody2OnBody1Equivalent )
            {
                std::shared_ptr< SphericalHarmonicsGravitationalAccelerationModel > sphericalHarmonicModel =
                        std::dynamic_pointer_cast< SphericalHarmonicsGravitationalAccelerationModel >( createAccelerationModel(
                                body1, body2, std::make_shared< SphericalHarmonicAccelerationSettings >( 2, 2 ), "Body1", "Body2" ) );
                // Require the equivalent one-way spherical-harmonic model for body-2 shape terms.
                BOOST_REQUIRE( sphericalHarmonicModel != nullptr );
                equivalentModel = sphericalHarmonicModel;
                equivalentPartial = std::make_shared< SphericalHarmonicsGravityPartial >( "Body1", "Body2", sphericalHarmonicModel );
            }
            else
            {
                std::shared_ptr< MutualSphericalHarmonicsGravitationalAccelerationModel > mutualSphericalModel =
                        std::dynamic_pointer_cast< MutualSphericalHarmonicsGravitationalAccelerationModel >(
                                createAccelerationModel( body1,
                                                         body2,
                                                         std::make_shared< MutualSphericalHarmonicAccelerationSettings >( 2, 2, 2, 2 ),
                                                         "Body1",
                                                         "Body2" ) );
                // Require the legacy mutual spherical-harmonic model for the single-point mutual reduction check.
                BOOST_REQUIRE( mutualSphericalModel != nullptr );
                equivalentModel = mutualSphericalModel;

                std::shared_ptr< SphericalHarmonicsGravityPartial > body2OnBody1Partial =
                        std::make_shared< SphericalHarmonicsGravityPartial >(
                                "Body1", "Body2", mutualSphericalModel->getAccelerationModelFromShExpansionOfBodyExertingAcceleration( ) );
                std::shared_ptr< SphericalHarmonicsGravityPartial > body1OnBody2Partial =
                        std::make_shared< SphericalHarmonicsGravityPartial >(
                                "Body2",
                                "Body1",
                                mutualSphericalModel->getAccelerationModelFromShExpansionOfBodyUndergoingAcceleration( ) );
                equivalentPartial = std::make_shared< MutualSphericalHarmonicsGravityPartial >(
                        body2OnBody1Partial, body1OnBody2Partial, "Body1", "Body2", mutualSphericalModel );
            }
            // Require a fully constructed equivalent model before comparing model-reduction partials.
            BOOST_REQUIRE( equivalentModel != nullptr );
            // Require the analytical partial object for the equivalent model.
            BOOST_REQUIRE( equivalentPartial != nullptr );

            equivalentModel->updateMembers( evaluationTime );
            equivalentPartial->update( evaluationTime );

            Eigen::MatrixXd equivalentPartialWrtBody1Position = Eigen::MatrixXd::Zero( 3, 3 );
            Eigen::MatrixXd equivalentPartialWrtBody2Position = Eigen::MatrixXd::Zero( 3, 3 );

            if( testCase.equivalentModelType == sphericalBody1OnBody2Equivalent )
            {
                equivalentPartial->wrtPositionOfAcceleratingBody( equivalentPartialWrtBody1Position.block( 0, 0, 3, 3 ) );
                equivalentPartial->wrtPositionOfAcceleratedBody( equivalentPartialWrtBody2Position.block( 0, 0, 3, 3 ) );
            }
            else
            {
                equivalentPartial->wrtPositionOfAcceleratedBody( equivalentPartialWrtBody1Position.block( 0, 0, 3, 3 ) );
                equivalentPartial->wrtPositionOfAcceleratingBody( equivalentPartialWrtBody2Position.block( 0, 0, 3, 3 ) );
            }

            equivalentPartialWrtBody1Position *= scaleToMutual;
            equivalentPartialWrtBody2Position *= scaleToMutual;

            // Verify full-two-body body-1 position partial is active in this equivalence case.
            BOOST_CHECK_GT( mutualPartialWrtBody1Position.norm( ), 1.0E-12 );
            // Verify full-two-body body-2 position partial is active in this equivalence case.
            BOOST_CHECK_GT( mutualPartialWrtBody2Position.norm( ), 1.0E-12 );
            // Verify equivalent-model body-1 position partial is active.
            BOOST_CHECK_GT( equivalentPartialWrtBody1Position.norm( ), 1.0E-12 );
            // Verify equivalent-model body-2 position partial is active.
            BOOST_CHECK_GT( equivalentPartialWrtBody2Position.norm( ), 1.0E-12 );

            // Position Jacobians from the full-two-body model must collapse to the equivalent model Jacobians.
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( mutualPartialWrtBody1Position, equivalentPartialWrtBody1Position, 2.0E-14 );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( mutualPartialWrtBody2Position, equivalentPartialWrtBody2Position, 2.0E-14 );

            if( testCase.compareBody1Coefficients )
            {
                Eigen::MatrixXd equivalentPartialWrtBody1Cosine =
                        equivalentPartial->wrtParameter( body1CosineCoefficientsParameter ) * scaleToMutual;
                Eigen::MatrixXd equivalentPartialWrtBody1Sine =
                        equivalentPartial->wrtParameter( body1SineCoefficientsParameter ) * scaleToMutual;

                // Verify full-two-body body-1 cosine coefficient partial is active before comparing values.
                BOOST_CHECK_GT( mutualPartialWrtBody1Cosine.norm( ), 1.0E-16 );
                // Verify full-two-body body-1 sine coefficient partial is active before comparing values.
                BOOST_CHECK_GT( mutualPartialWrtBody1Sine.norm( ), 1.0E-16 );
                // Verify equivalent-model body-1 cosine coefficient partial is active.
                BOOST_CHECK_GT( equivalentPartialWrtBody1Cosine.norm( ), 1.0E-16 );
                // Verify equivalent-model body-1 sine coefficient partial is active.
                BOOST_CHECK_GT( equivalentPartialWrtBody1Sine.norm( ), 1.0E-16 );

                // Verify body-1 cosine coefficient partial reduction to the equivalent model.
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( mutualPartialWrtBody1Cosine, equivalentPartialWrtBody1Cosine, 2.0E-14 );
                // Verify body-1 sine coefficient partial reduction to the equivalent model.
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( mutualPartialWrtBody1Sine, equivalentPartialWrtBody1Sine, 2.0E-14 );
            }
            else
            {
                // Verify body-1 cosine coefficients are inactive when absent from the interaction subset.
                BOOST_CHECK_SMALL( mutualPartialWrtBody1Cosine.norm( ), 1.0E-20 );
                // Verify body-1 sine coefficients are inactive when absent from the interaction subset.
                BOOST_CHECK_SMALL( mutualPartialWrtBody1Sine.norm( ), 1.0E-20 );
            }

            if( testCase.compareBody2Coefficients )
            {
                Eigen::MatrixXd equivalentPartialWrtBody2Cosine =
                        equivalentPartial->wrtParameter( body2CosineCoefficientsParameter ) * scaleToMutual;
                Eigen::MatrixXd equivalentPartialWrtBody2Sine =
                        equivalentPartial->wrtParameter( body2SineCoefficientsParameter ) * scaleToMutual;

                // Verify full-two-body body-2 cosine coefficient partial is active before comparing values.
                BOOST_CHECK_GT( mutualPartialWrtBody2Cosine.norm( ), 1.0E-16 );
                // Verify full-two-body body-2 sine coefficient partial is active before comparing values.
                BOOST_CHECK_GT( mutualPartialWrtBody2Sine.norm( ), 1.0E-16 );
                // Verify equivalent-model body-2 cosine coefficient partial is active.
                BOOST_CHECK_GT( equivalentPartialWrtBody2Cosine.norm( ), 1.0E-16 );
                // Verify equivalent-model body-2 sine coefficient partial is active.
                BOOST_CHECK_GT( equivalentPartialWrtBody2Sine.norm( ), 1.0E-16 );

                // Verify body-2 cosine coefficient partial reduction to the equivalent model.
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( mutualPartialWrtBody2Cosine, equivalentPartialWrtBody2Cosine, 2.0E-14 );
                // Verify body-2 sine coefficient partial reduction to the equivalent model.
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( mutualPartialWrtBody2Sine, equivalentPartialWrtBody2Sine, 2.0E-14 );
            }
            else
            {
                // Verify body-2 cosine coefficients are inactive when absent from the interaction subset.
                BOOST_CHECK_SMALL( mutualPartialWrtBody2Cosine.norm( ), 1.0E-20 );
                // Verify body-2 sine coefficients are inactive when absent from the interaction subset.
                BOOST_CHECK_SMALL( mutualPartialWrtBody2Sine.norm( ), 1.0E-20 );
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( testFullTwoBodySphericalHarmonicGravityPartialsThirdBody )
{
    // Test rationale:
    // Validate third-body partials for the full-two-body acceleration model in a non-trivial geometry where all
    // three bodies have degree-2 terms. This checks:
    // 1) partials w.r.t. accelerated, accelerating, and central-body positions;
    // 2) partials w.r.t. degree-2 coefficients of all three bodies.
    //
    // The formulation corresponds to third-body application of the same mutual-potential interaction framework
    // used in Dirkx et al. (2019), now with central-body state dependence included.
    const double gravitationalParameterBody1 = 5.0E5;
    const double gravitationalParameterBody2 = 6.0E5;
    const double gravitationalParameterCentralBody = 7.0E5;
    const double equatorialRadiusBody1 = 1300.0;
    const double equatorialRadiusBody2 = 1100.0;
    const double equatorialRadiusCentralBody = 1700.0;

    Eigen::MatrixXd cosineCoefficientsOfBody1 = Eigen::MatrixXd::Zero( 3, 3 );
    Eigen::MatrixXd sineCoefficientsOfBody1 = Eigen::MatrixXd::Zero( 3, 3 );
    Eigen::MatrixXd cosineCoefficientsOfBody2 = Eigen::MatrixXd::Zero( 3, 3 );
    Eigen::MatrixXd sineCoefficientsOfBody2 = Eigen::MatrixXd::Zero( 3, 3 );
    Eigen::MatrixXd cosineCoefficientsOfCentralBody = Eigen::MatrixXd::Zero( 3, 3 );
    Eigen::MatrixXd sineCoefficientsOfCentralBody = Eigen::MatrixXd::Zero( 3, 3 );

    cosineCoefficientsOfBody1( 0, 0 ) = 1.0;
    cosineCoefficientsOfBody2( 0, 0 ) = 1.0;
    cosineCoefficientsOfCentralBody( 0, 0 ) = 1.0;

    cosineCoefficientsOfBody1( 2, 0 ) = 0.23;
    cosineCoefficientsOfBody1( 2, 1 ) = -0.11;
    sineCoefficientsOfBody1( 2, 1 ) = 0.14;
    cosineCoefficientsOfBody1( 2, 2 ) = 0.09;
    sineCoefficientsOfBody1( 2, 2 ) = -0.06;

    cosineCoefficientsOfBody2( 2, 0 ) = -0.19;
    cosineCoefficientsOfBody2( 2, 1 ) = 0.16;
    sineCoefficientsOfBody2( 2, 1 ) = -0.08;
    cosineCoefficientsOfBody2( 2, 2 ) = 0.13;
    sineCoefficientsOfBody2( 2, 2 ) = 0.12;

    cosineCoefficientsOfCentralBody( 2, 0 ) = 0.21;
    cosineCoefficientsOfCentralBody( 2, 1 ) = -0.05;
    sineCoefficientsOfCentralBody( 2, 1 ) = 0.07;
    cosineCoefficientsOfCentralBody( 2, 2 ) = -0.09;
    sineCoefficientsOfCentralBody( 2, 2 ) = 0.11;

    std::shared_ptr< Body > body1 = std::make_shared< Body >( );
    std::shared_ptr< Body > body2 = std::make_shared< Body >( );
    std::shared_ptr< Body > centralBody = std::make_shared< Body >( );

    body1->setGravityFieldModel( std::make_shared< SphericalHarmonicsGravityField >(
            gravitationalParameterBody1, equatorialRadiusBody1, cosineCoefficientsOfBody1, sineCoefficientsOfBody1, "IAU_Body1" ) );
    body2->setGravityFieldModel( std::make_shared< SphericalHarmonicsGravityField >(
            gravitationalParameterBody2, equatorialRadiusBody2, cosineCoefficientsOfBody2, sineCoefficientsOfBody2, "IAU_Body2" ) );
    centralBody->setGravityFieldModel( std::make_shared< SphericalHarmonicsGravityField >( gravitationalParameterCentralBody,
                                                                                           equatorialRadiusCentralBody,
                                                                                           cosineCoefficientsOfCentralBody,
                                                                                           sineCoefficientsOfCentralBody,
                                                                                           "IAU_CentralBody" ) );

    const double evaluationTime = 1000.0;
    body1->setRotationalEphemeris( std::make_shared< ephemerides::SimpleRotationalEphemeris >(
            Eigen::Quaterniond::Identity( ), 0.0, evaluationTime, "ECLIPJ2000", "IAU_Body1" ) );
    body2->setRotationalEphemeris( std::make_shared< ephemerides::SimpleRotationalEphemeris >(
            Eigen::Quaterniond( Eigen::AngleAxisd( 0.3, Eigen::Vector3d::UnitZ( ) ) *
                                Eigen::AngleAxisd( -0.2, Eigen::Vector3d::UnitX( ) ) ),
            0.0,
            evaluationTime,
            "ECLIPJ2000",
            "IAU_Body2" ) );
    centralBody->setRotationalEphemeris( std::make_shared< ephemerides::SimpleRotationalEphemeris >(
            Eigen::Quaterniond( Eigen::AngleAxisd( -0.4, Eigen::Vector3d::UnitY( ) ) *
                                Eigen::AngleAxisd( 0.15, Eigen::Vector3d::UnitZ( ) ) ),
            0.0,
            evaluationTime,
            "ECLIPJ2000",
            "IAU_CentralBody" ) );

    body1->setCurrentRotationToLocalFrameFromEphemeris( evaluationTime );
    body2->setCurrentRotationToLocalFrameFromEphemeris( evaluationTime );
    centralBody->setCurrentRotationToLocalFrameFromEphemeris( evaluationTime );

    Eigen::Vector6d body1State = Eigen::Vector6d::Zero( );
    body1State.segment( 0, 3 ) = ( Eigen::Vector3d( ) << 5100.0, -2200.0, 3600.0 ).finished( );
    body1->setState( body1State );

    Eigen::Vector6d body2State = Eigen::Vector6d::Zero( );
    body2State.segment( 0, 3 ) = ( Eigen::Vector3d( ) << -700.0, 1300.0, 900.0 ).finished( );
    body2->setState( body2State );

    Eigen::Vector6d centralBodyState = Eigen::Vector6d::Zero( );
    centralBodyState.segment( 0, 3 ) = ( Eigen::Vector3d( ) << 4200.0, -3600.0, 2500.0 ).finished( );
    centralBody->setState( centralBodyState );

    SystemOfBodies bodies;
    bodies.addBody( body1, "Body1" );
    bodies.addBody( body2, "Body2" );
    bodies.addBody( centralBody, "CentralBody" );

    std::shared_ptr< AccelerationSettings > mutualExtendedSettings =
            std::make_shared< FullTwoBodySphericalHarmonicAccelerationSettings >( 2, 2, 2, 2, 2, 2 );
    std::shared_ptr< ThirdBodyFullTwoBodySphericalHarmonicsGravitationalAccelerationModel > mutualExtendedModel =
            std::dynamic_pointer_cast< ThirdBodyFullTwoBodySphericalHarmonicsGravitationalAccelerationModel >(
                    createAccelerationModel( body1, body2, mutualExtendedSettings, "Body1", "Body2", centralBody, "CentralBody", bodies ) );
    // Require the third-body full-two-body acceleration model before testing its partials.
    BOOST_REQUIRE( mutualExtendedModel != nullptr );

    mutualExtendedModel->updateMembers( evaluationTime );
    // Verify the selected three-body geometry produces a nonzero acceleration response.
    BOOST_CHECK_GT( mutualExtendedModel->getAcceleration( ).norm( ), 1.0E-16 );

    std::vector< std::shared_ptr< EstimatableParameter< double > > > doubleParameters;
    std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > > vectorParameters;
    std::shared_ptr< EstimatableParameterSet< double > > parameterSet =
            std::make_shared< EstimatableParameterSet< double > >( doubleParameters, vectorParameters );

    std::shared_ptr< AccelerationPartial > mutualExtendedPartial = createAnalyticalAccelerationPartial(
            mutualExtendedModel, std::make_pair( "Body1", body1 ), std::make_pair( "Body2", body2 ), bodies, parameterSet );

    // Require creation of the analytical third-body partial object.
    BOOST_REQUIRE( mutualExtendedPartial != nullptr );

    mutualExtendedPartial->update( evaluationTime );

    Eigen::MatrixXd partialWrtBody1PositionExtended = Eigen::MatrixXd::Zero( 3, 3 );
    Eigen::MatrixXd partialWrtBody2PositionExtended = Eigen::MatrixXd::Zero( 3, 3 );
    Eigen::MatrixXd partialWrtCentralBodyPositionExtended = Eigen::MatrixXd::Zero( 3, 3 );

    mutualExtendedPartial->wrtPositionOfAcceleratedBody( partialWrtBody1PositionExtended.block( 0, 0, 3, 3 ) );
    mutualExtendedPartial->wrtPositionOfAcceleratingBody( partialWrtBody2PositionExtended.block( 0, 0, 3, 3 ) );
    mutualExtendedPartial->wrtPositionOfAdditionalBody( "CentralBody", partialWrtCentralBodyPositionExtended.block( 0, 0, 3, 3 ) );

    // Verify body-1 position partial is active in the third-body geometry.
    BOOST_CHECK_GT( partialWrtBody1PositionExtended.norm( ), 1.0E-16 );
    // Verify body-2 position partial is active in the third-body geometry.
    BOOST_CHECK_GT( partialWrtBody2PositionExtended.norm( ), 1.0E-16 );
    // Verify central-body position partial is active in the third-body geometry.
    BOOST_CHECK_GT( partialWrtCentralBodyPositionExtended.norm( ), 1.0E-16 );

    std::shared_ptr< SphericalHarmonicsCosineCoefficients > body1CosineCoefficientsParameter =
            std::make_shared< SphericalHarmonicsCosineCoefficients >(
                    std::bind( &SphericalHarmonicsGravityField::getCosineCoefficients,
                               std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( body1->getGravityFieldModel( ) ) ),
                    std::bind( &SphericalHarmonicsGravityField::setCosineCoefficients,
                               std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( body1->getGravityFieldModel( ) ),
                               std::placeholders::_1 ),
                    std::vector< std::pair< int, int > >( { { 2, 0 }, { 2, 1 }, { 2, 2 } } ),
                    "Body1" );
    std::shared_ptr< SphericalHarmonicsSineCoefficients > body1SineCoefficientsParameter =
            std::make_shared< SphericalHarmonicsSineCoefficients >(
                    std::bind( &SphericalHarmonicsGravityField::getSineCoefficients,
                               std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( body1->getGravityFieldModel( ) ) ),
                    std::bind( &SphericalHarmonicsGravityField::setSineCoefficients,
                               std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( body1->getGravityFieldModel( ) ),
                               std::placeholders::_1 ),
                    std::vector< std::pair< int, int > >( { { 2, 1 }, { 2, 2 } } ),
                    "Body1" );

    std::shared_ptr< SphericalHarmonicsCosineCoefficients > body2CosineCoefficientsParameter =
            std::make_shared< SphericalHarmonicsCosineCoefficients >(
                    std::bind( &SphericalHarmonicsGravityField::getCosineCoefficients,
                               std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( body2->getGravityFieldModel( ) ) ),
                    std::bind( &SphericalHarmonicsGravityField::setCosineCoefficients,
                               std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( body2->getGravityFieldModel( ) ),
                               std::placeholders::_1 ),
                    std::vector< std::pair< int, int > >( { { 2, 0 }, { 2, 1 }, { 2, 2 } } ),
                    "Body2" );
    std::shared_ptr< SphericalHarmonicsSineCoefficients > body2SineCoefficientsParameter =
            std::make_shared< SphericalHarmonicsSineCoefficients >(
                    std::bind( &SphericalHarmonicsGravityField::getSineCoefficients,
                               std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( body2->getGravityFieldModel( ) ) ),
                    std::bind( &SphericalHarmonicsGravityField::setSineCoefficients,
                               std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( body2->getGravityFieldModel( ) ),
                               std::placeholders::_1 ),
                    std::vector< std::pair< int, int > >( { { 2, 1 }, { 2, 2 } } ),
                    "Body2" );

    std::shared_ptr< SphericalHarmonicsCosineCoefficients > centralCosineCoefficientsParameter =
            std::make_shared< SphericalHarmonicsCosineCoefficients >(
                    std::bind( &SphericalHarmonicsGravityField::getCosineCoefficients,
                               std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( centralBody->getGravityFieldModel( ) ) ),
                    std::bind( &SphericalHarmonicsGravityField::setCosineCoefficients,
                               std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( centralBody->getGravityFieldModel( ) ),
                               std::placeholders::_1 ),
                    std::vector< std::pair< int, int > >( { { 2, 0 }, { 2, 1 }, { 2, 2 } } ),
                    "CentralBody" );
    std::shared_ptr< SphericalHarmonicsSineCoefficients > centralSineCoefficientsParameter =
            std::make_shared< SphericalHarmonicsSineCoefficients >(
                    std::bind( &SphericalHarmonicsGravityField::getSineCoefficients,
                               std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( centralBody->getGravityFieldModel( ) ) ),
                    std::bind( &SphericalHarmonicsGravityField::setSineCoefficients,
                               std::dynamic_pointer_cast< SphericalHarmonicsGravityField >( centralBody->getGravityFieldModel( ) ),
                               std::placeholders::_1 ),
                    std::vector< std::pair< int, int > >( { { 2, 1 }, { 2, 2 } } ),
                    "CentralBody" );

    const Eigen::MatrixXd partialWrtBody1CosineExtended = mutualExtendedPartial->wrtParameter( body1CosineCoefficientsParameter );
    const Eigen::MatrixXd partialWrtBody1SineExtended = mutualExtendedPartial->wrtParameter( body1SineCoefficientsParameter );
    const Eigen::MatrixXd partialWrtBody2CosineExtended = mutualExtendedPartial->wrtParameter( body2CosineCoefficientsParameter );
    const Eigen::MatrixXd partialWrtBody2SineExtended = mutualExtendedPartial->wrtParameter( body2SineCoefficientsParameter );
    const Eigen::MatrixXd partialWrtCentralCosineExtended = mutualExtendedPartial->wrtParameter( centralCosineCoefficientsParameter );
    const Eigen::MatrixXd partialWrtCentralSineExtended = mutualExtendedPartial->wrtParameter( centralSineCoefficientsParameter );

    // Verify body-1 cosine coefficient partial is active.
    BOOST_CHECK_GT( partialWrtBody1CosineExtended.norm( ), 1.0E-16 );
    // Verify body-1 sine coefficient partial is active.
    BOOST_CHECK_GT( partialWrtBody1SineExtended.norm( ), 1.0E-16 );
    // Verify body-2 cosine coefficient partial is active.
    BOOST_CHECK_GT( partialWrtBody2CosineExtended.norm( ), 1.0E-16 );
    // Verify body-2 sine coefficient partial is active.
    BOOST_CHECK_GT( partialWrtBody2SineExtended.norm( ), 1.0E-16 );
    // Verify central-body cosine coefficient partial is active.
    BOOST_CHECK_GT( partialWrtCentralCosineExtended.norm( ), 1.0E-16 );
    // Verify central-body sine coefficient partial is active.
    BOOST_CHECK_GT( partialWrtCentralSineExtended.norm( ), 1.0E-16 );

    const Eigen::Vector3d positionPerturbation = Eigen::Vector3d::Constant( 10.0 );
    std::function< void( Eigen::Vector6d ) > body1StateSetFunction = std::bind( &Body::setState, body1, std::placeholders::_1 );
    std::function< void( Eigen::Vector6d ) > body2StateSetFunction = std::bind( &Body::setState, body2, std::placeholders::_1 );
    std::function< void( Eigen::Vector6d ) > centralStateSetFunction = std::bind( &Body::setState, centralBody, std::placeholders::_1 );

    Eigen::Matrix3d numericalPartialWrtBody1Position = calculateAccelerationWrtStatePartials(
            body1StateSetFunction, mutualExtendedModel, body1->getState( ), positionPerturbation, 0 );
    Eigen::Matrix3d numericalPartialWrtBody2Position = calculateAccelerationWrtStatePartials(
            body2StateSetFunction, mutualExtendedModel, body2->getState( ), positionPerturbation, 0 );
    Eigen::Matrix3d numericalPartialWrtCentralPosition = calculateAccelerationWrtStatePartials(
            centralStateSetFunction, mutualExtendedModel, centralBody->getState( ), positionPerturbation, 0 );

    // Verify body-1 third-body position partial against central finite differences.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( numericalPartialWrtBody1Position, partialWrtBody1PositionExtended, 5.0E-5 );
    // Verify body-2 third-body position partial against central finite differences.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( numericalPartialWrtBody2Position, partialWrtBody2PositionExtended, 5.0E-5 );
    // Verify central-body third-body position partial against central finite differences.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( numericalPartialWrtCentralPosition, partialWrtCentralBodyPositionExtended, 5.0E-5 );

    Eigen::MatrixXd numericalPartialWrtBody1Cosine = calculateAccelerationWrtParameterPartials(
            body1CosineCoefficientsParameter,
            mutualExtendedModel,
            Eigen::VectorXd::Constant( body1CosineCoefficientsParameter->getParameterSize( ), 1.0E-2 ) );
    Eigen::MatrixXd numericalPartialWrtBody1Sine = calculateAccelerationWrtParameterPartials(
            body1SineCoefficientsParameter,
            mutualExtendedModel,
            Eigen::VectorXd::Constant( body1SineCoefficientsParameter->getParameterSize( ), 1.0E-2 ) );
    Eigen::MatrixXd numericalPartialWrtBody2Cosine = calculateAccelerationWrtParameterPartials(
            body2CosineCoefficientsParameter,
            mutualExtendedModel,
            Eigen::VectorXd::Constant( body2CosineCoefficientsParameter->getParameterSize( ), 1.0E-2 ) );
    Eigen::MatrixXd numericalPartialWrtBody2Sine = calculateAccelerationWrtParameterPartials(
            body2SineCoefficientsParameter,
            mutualExtendedModel,
            Eigen::VectorXd::Constant( body2SineCoefficientsParameter->getParameterSize( ), 1.0E-2 ) );
    Eigen::MatrixXd numericalPartialWrtCentralCosine = calculateAccelerationWrtParameterPartials(
            centralCosineCoefficientsParameter,
            mutualExtendedModel,
            Eigen::VectorXd::Constant( centralCosineCoefficientsParameter->getParameterSize( ), 1.0E-2 ) );
    Eigen::MatrixXd numericalPartialWrtCentralSine = calculateAccelerationWrtParameterPartials(
            centralSineCoefficientsParameter,
            mutualExtendedModel,
            Eigen::VectorXd::Constant( centralSineCoefficientsParameter->getParameterSize( ), 1.0E-2 ) );

    // Verify body-1 cosine coefficient third-body partial against finite differences.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( numericalPartialWrtBody1Cosine, partialWrtBody1CosineExtended, 1.0E-8 );
    // Verify body-1 sine coefficient third-body partial against finite differences.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( numericalPartialWrtBody1Sine, partialWrtBody1SineExtended, 1.0E-8 );
    // Verify body-2 cosine coefficient third-body partial against finite differences.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( numericalPartialWrtBody2Cosine, partialWrtBody2CosineExtended, 1.0E-8 );
    // Verify body-2 sine coefficient third-body partial against finite differences.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( numericalPartialWrtBody2Sine, partialWrtBody2SineExtended, 1.0E-8 );
    // Verify central-body cosine coefficient third-body partial against finite differences.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( numericalPartialWrtCentralCosine, partialWrtCentralCosineExtended, 1.0E-8 );
    // Verify central-body sine coefficient third-body partial against finite differences.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( numericalPartialWrtCentralSine, partialWrtCentralSineExtended, 1.0E-8 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
