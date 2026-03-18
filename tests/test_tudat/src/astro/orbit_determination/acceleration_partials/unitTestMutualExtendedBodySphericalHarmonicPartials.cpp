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
#include "tudat/astro/gravitation/sphericalHarmonicsGravityField.h"
#include "tudat/astro/gravitation/sphericalHarmonicsGravityModel.h"
#include "tudat/astro/orbit_determination/acceleration_partials/centralGravityAccelerationPartial.h"
#include "tudat/astro/orbit_determination/acceleration_partials/fullTwoBodySphericalHarmonicGravityPartial.h"
#include "tudat/astro/orbit_determination/acceleration_partials/mutualSphericalHarmonicGravityPartial.h"
#include "tudat/astro/orbit_determination/acceleration_partials/numericalAccelerationPartial.h"
#include "tudat/astro/orbit_determination/acceleration_partials/sphericalHarmonicAccelerationPartial.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameterSet.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/sphericalHarmonicCosineCoefficients.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/sphericalHarmonicSineCoefficients.h"
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
    struct PartialSet
    {
        Eigen::MatrixXd wrtBody1Position;
        Eigen::MatrixXd wrtBody2Position;
        Eigen::MatrixXd wrtBody1Velocity;
        Eigen::MatrixXd wrtBody2Velocity;
        Eigen::MatrixXd wrtBody1Cosine;
        Eigen::MatrixXd wrtBody1Sine;
        Eigen::MatrixXd wrtBody2Cosine;
        Eigen::MatrixXd wrtBody2Sine;
    };

    struct CaseDefinition
    {
        std::string name;
        std::vector< std::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > coefficientCombinations;
        bool hasBody2ShapeTerms;
        bool hasFigureFigureTerms;
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

    const std::vector< CaseDefinition > testCases = {
        { "regular",
          FullTwoBodySphericalHarmonicAccelerationSettings( 2, 2, 0, 0 ).coefficientCombinationsToUse_,
          false,
          false },
        { "mutualSinglePoint",
          getExtendedSinglePointMassInteractions( 2, 2, 2, 2 ),
          true,
          false },
        { "mutual",
          FullTwoBodySphericalHarmonicAccelerationSettings( 2, 2, 2, 2 ).coefficientCombinationsToUse_,
          true,
          true } };

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
            rotationToBody1 = Eigen::Quaterniond(
                    Eigen::AngleAxisd( 0.7, Eigen::Vector3d::UnitZ( ) ) *
                    Eigen::AngleAxisd( -0.4, Eigen::Vector3d::UnitX( ) ) *
                    Eigen::AngleAxisd( 0.2, Eigen::Vector3d::UnitY( ) ) );
            rotationToBody2 = Eigen::Quaterniond(
                    Eigen::AngleAxisd( -0.3, Eigen::Vector3d::UnitZ( ) ) *
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

            std::function< void( Eigen::Vector6d ) > body1StateSetFunction =
                    std::bind( &Body::setState, body1, std::placeholders::_1 );
            std::function< void( Eigen::Vector6d ) > body2StateSetFunction =
                    std::bind( &Body::setState, body2, std::placeholders::_1 );

            std::shared_ptr< AccelerationSettings > accelerationSettings =
                    std::make_shared< FullTwoBodySphericalHarmonicAccelerationSettings >(
                            testCase.coefficientCombinations );
            std::shared_ptr< FullTwoBodySphericalHarmonicAcceleration > accelerationModel =
                    std::dynamic_pointer_cast< FullTwoBodySphericalHarmonicAcceleration >(
                            createAccelerationModel( body1, body2, accelerationSettings, "Body1", "Body2" ) );

            BOOST_REQUIRE( accelerationModel != nullptr );

            std::shared_ptr< FullTwoBodySphericalHarmonicsGravityPartial > accelerationPartial =
                    std::make_shared< FullTwoBodySphericalHarmonicsGravityPartial >(
                            "Body1", "Body2", accelerationModel );

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

            // Finite-difference reference:
            // each acceleration evaluation advances the time tag by 1 second so the bodies refresh their current
            // rotational states through the ephemeris interface before model evaluation.
            double finiteDifferenceTimeTag = evaluationTime;
            auto evaluateAcceleration = [ & ]( )
            {
                finiteDifferenceTimeTag += 1.0;
                body1->setCurrentRotationToLocalFrameFromEphemeris( finiteDifferenceTimeTag );
                body2->setCurrentRotationToLocalFrameFromEphemeris( finiteDifferenceTimeTag );
                accelerationModel->updateMembers( finiteDifferenceTimeTag );
                return accelerationModel->getAcceleration( );
            };

            auto calculateStatePartial = [ & ]( const std::function< void( Eigen::Vector6d ) >& setState,
                                                Eigen::Vector6d nominalState,
                                                const Eigen::Vector3d& perturbation,
                                                const int startIndex )
            {
                Eigen::Matrix3d partial = Eigen::Matrix3d::Zero( );
                Eigen::Vector6d perturbedState = nominalState;

                for( int i = 0; i < 3; i++ )
                {
                    perturbedState = nominalState;
                    perturbedState( startIndex + i ) += perturbation( i );
                    setState( perturbedState );
                    const Eigen::Vector3d upAcceleration = evaluateAcceleration( );

                    perturbedState = nominalState;
                    perturbedState( startIndex + i ) -= perturbation( i );
                    setState( perturbedState );
                    const Eigen::Vector3d downAcceleration = evaluateAcceleration( );

                    partial.col( i ) = ( upAcceleration - downAcceleration ) / ( 2.0 * perturbation( i ) );
                }

                setState( nominalState );
                evaluateAcceleration( );
                return partial;
            };

            auto calculateVectorParameterPartial = [ & ](
                                                          const std::shared_ptr< EstimatableParameter< Eigen::VectorXd > >& parameter,
                                                          const Eigen::VectorXd& perturbation )
            {
                const Eigen::VectorXd nominalValue = parameter->getParameterValue( );
                Eigen::MatrixXd partial = Eigen::MatrixXd::Zero( 3, nominalValue.size( ) );
                Eigen::VectorXd perturbedValue = nominalValue;

                for( int i = 0; i < nominalValue.size( ); i++ )
                {
                    perturbedValue = nominalValue;
                    perturbedValue( i ) += perturbation( i );
                    parameter->setParameterValue( perturbedValue );
                    const Eigen::Vector3d upAcceleration = evaluateAcceleration( );

                    perturbedValue = nominalValue;
                    perturbedValue( i ) -= perturbation( i );
                    parameter->setParameterValue( perturbedValue );
                    const Eigen::Vector3d downAcceleration = evaluateAcceleration( );

                    partial.col( i ) = ( upAcceleration - downAcceleration ) / ( 2.0 * perturbation( i ) );
                }

                parameter->setParameterValue( nominalValue );
                evaluateAcceleration( );
                return partial;
            };

            const Eigen::Matrix3d numericalPartialWrtBody1Position =
                    calculateStatePartial( body1StateSetFunction, body1->getState( ), positionPerturbation, 0 );
            const Eigen::Matrix3d numericalPartialWrtBody1Velocity =
                    calculateStatePartial( body1StateSetFunction, body1->getState( ), velocityPerturbation, 3 );
            const Eigen::Matrix3d numericalPartialWrtBody2Position =
                    calculateStatePartial( body2StateSetFunction, body2->getState( ), positionPerturbation, 0 );
            const Eigen::Matrix3d numericalPartialWrtBody2Velocity =
                    calculateStatePartial( body2StateSetFunction, body2->getState( ), velocityPerturbation, 3 );

            const Eigen::MatrixXd numericalPartialWrtBody1Cosine = calculateVectorParameterPartial(
                    body1CosineCoefficientsParameter,
                    Eigen::VectorXd::Constant( body1CosineCoefficientsParameter->getParameterSize( ), 1.0E-8 ) );
            const Eigen::MatrixXd numericalPartialWrtBody1Sine = calculateVectorParameterPartial(
                    body1SineCoefficientsParameter,
                    Eigen::VectorXd::Constant( body1SineCoefficientsParameter->getParameterSize( ), 1.0E-8 ) );
            const Eigen::MatrixXd numericalPartialWrtBody2Cosine = calculateVectorParameterPartial(
                    body2CosineCoefficientsParameter,
                    Eigen::VectorXd::Constant( body2CosineCoefficientsParameter->getParameterSize( ), 1.0E-8 ) );
            const Eigen::MatrixXd numericalPartialWrtBody2Sine = calculateVectorParameterPartial(
                    body2SineCoefficientsParameter,
                    Eigen::VectorXd::Constant( body2SineCoefficientsParameter->getParameterSize( ), 1.0E-8 ) );

            // State Jacobians: analytical partial blocks should match numerical references.
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    numericalPartialWrtBody1Position, analyticalPartials.wrtBody1Position, 5.0E-5 );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    numericalPartialWrtBody2Position, analyticalPartials.wrtBody2Position, 5.0E-5 );

            // Gravitational acceleration has no explicit velocity dependence in this model.
            BOOST_CHECK_SMALL( analyticalPartials.wrtBody1Velocity.norm( ), std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_SMALL( analyticalPartials.wrtBody2Velocity.norm( ), std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_SMALL( numericalPartialWrtBody1Velocity.norm( ), std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_SMALL( numericalPartialWrtBody2Velocity.norm( ), std::numeric_limits< double >::epsilon( ) );

            // Coefficient Jacobians: active coefficient sets must match finite-difference references.
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    numericalPartialWrtBody1Cosine, analyticalPartials.wrtBody1Cosine, 1.0E-3 );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    numericalPartialWrtBody1Sine, analyticalPartials.wrtBody1Sine, 1.0E-3 );

            if( testCase.hasBody2ShapeTerms )
            {
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        numericalPartialWrtBody2Cosine, analyticalPartials.wrtBody2Cosine, 1.0E-3 );
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        numericalPartialWrtBody2Sine, analyticalPartials.wrtBody2Sine, 1.0E-3 );
            }
            else
            {
                BOOST_CHECK_SMALL( analyticalPartials.wrtBody2Cosine.norm( ), 1.0E-20 );
                BOOST_CHECK_SMALL( analyticalPartials.wrtBody2Sine.norm( ), 1.0E-20 );
                BOOST_CHECK_SMALL( numericalPartialWrtBody2Cosine.norm( ), 1.0E-20 );
                BOOST_CHECK_SMALL( numericalPartialWrtBody2Sine.norm( ), 1.0E-20 );
            }

            analyticalPartialsByCase[ testCase.name ] = analyticalPartials;
        }

        BOOST_CHECK_GT( ( analyticalPartialsByCase.at( "mutual" ).wrtBody1Position -
                          analyticalPartialsByCase.at( "regular" ).wrtBody1Position )
                                .norm( ),
                        1.0E-7 );
        BOOST_CHECK_GT( ( analyticalPartialsByCase.at( "mutualSinglePoint" ).wrtBody1Position -
                          analyticalPartialsByCase.at( "regular" ).wrtBody1Position )
                                .norm( ),
                        1.0E-7 );
        BOOST_CHECK_GT( ( analyticalPartialsByCase.at( "mutual" ).wrtBody1Cosine -
                          analyticalPartialsByCase.at( "mutualSinglePoint" ).wrtBody1Cosine )
                                .norm( ),
                        1.0E-14 );
        BOOST_CHECK_GT( analyticalPartialsByCase.at( "mutual" ).wrtBody2Cosine.norm( ), 1.0E-8 );
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
    enum EquivalentModelType
    {
        pointMassEquivalent = 0,
        sphericalBody1OnBody2Equivalent = 1,
        sphericalBody2OnBody1Equivalent = 2,
        mutualSphericalEquivalent = 3
    };

    struct EquivalentCaseDefinition
    {
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
          std::make_shared< FullTwoBodySphericalHarmonicAccelerationSettings >(
                  getExtendedSinglePointMassInteractions( 2, 2, 2, 2 ) ),
          mutualSphericalEquivalent,
          true,
          true } };

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
            rotationToBody1 = Eigen::Quaterniond(
                    Eigen::AngleAxisd( 0.7, Eigen::Vector3d::UnitZ( ) ) *
                    Eigen::AngleAxisd( -0.4, Eigen::Vector3d::UnitX( ) ) *
                    Eigen::AngleAxisd( 0.2, Eigen::Vector3d::UnitY( ) ) );
            rotationToBody2 = Eigen::Quaterniond(
                    Eigen::AngleAxisd( -0.3, Eigen::Vector3d::UnitZ( ) ) *
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
            BOOST_REQUIRE( mutualExtendedModel != nullptr );

            std::shared_ptr< FullTwoBodySphericalHarmonicsGravityPartial > mutualExtendedPartial =
                    std::make_shared< FullTwoBodySphericalHarmonicsGravityPartial >(
                            "Body1", "Body2", mutualExtendedModel );

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
                        std::dynamic_pointer_cast< CentralGravitationalAccelerationModel3d >(
                                createAccelerationModel(
                                        body1,
                                        body2,
                                        std::make_shared< AccelerationSettings >( point_mass_gravity ),
                                        "Body1",
                                        "Body2" ) );
                BOOST_REQUIRE( pointMassModel != nullptr );
                equivalentModel = pointMassModel;
                equivalentPartial = std::make_shared< CentralGravitationPartial >( pointMassModel, "Body1", "Body2" );
            }
            else if( testCase.equivalentModelType == sphericalBody1OnBody2Equivalent )
            {
                std::shared_ptr< SphericalHarmonicsGravitationalAccelerationModel > sphericalHarmonicModel =
                        std::dynamic_pointer_cast< SphericalHarmonicsGravitationalAccelerationModel >(
                                createAccelerationModel(
                                        body2,
                                        body1,
                                        std::make_shared< SphericalHarmonicAccelerationSettings >( 2, 2 ),
                                        "Body2",
                                        "Body1" ) );
                BOOST_REQUIRE( sphericalHarmonicModel != nullptr );
                equivalentModel = sphericalHarmonicModel;
                equivalentPartial =
                        std::make_shared< SphericalHarmonicsGravityPartial >( "Body2", "Body1", sphericalHarmonicModel );
                // Sign/scale conversion: this equivalent model computes acceleration of body 2 due to body 1.
                // The mutual model here is defined for body 1 due to body 2, hence the factor below.
                scaleToMutual = -body2GravityField->getGravitationalParameter( ) /
                        body1GravityField->getGravitationalParameter( );
            }
            else if( testCase.equivalentModelType == sphericalBody2OnBody1Equivalent )
            {
                std::shared_ptr< SphericalHarmonicsGravitationalAccelerationModel > sphericalHarmonicModel =
                        std::dynamic_pointer_cast< SphericalHarmonicsGravitationalAccelerationModel >(
                                createAccelerationModel(
                                        body1,
                                        body2,
                                        std::make_shared< SphericalHarmonicAccelerationSettings >( 2, 2 ),
                                        "Body1",
                                        "Body2" ) );
                BOOST_REQUIRE( sphericalHarmonicModel != nullptr );
                equivalentModel = sphericalHarmonicModel;
                equivalentPartial =
                        std::make_shared< SphericalHarmonicsGravityPartial >( "Body1", "Body2", sphericalHarmonicModel );
            }
            else
            {
                std::shared_ptr< MutualSphericalHarmonicsGravitationalAccelerationModel > mutualSphericalModel =
                        std::dynamic_pointer_cast< MutualSphericalHarmonicsGravitationalAccelerationModel >(
                                createAccelerationModel(
                                        body1,
                                        body2,
                                        std::make_shared< MutualSphericalHarmonicAccelerationSettings >( 2, 2, 2, 2 ),
                                        "Body1",
                                        "Body2" ) );
                BOOST_REQUIRE( mutualSphericalModel != nullptr );
                equivalentModel = mutualSphericalModel;

                std::shared_ptr< SphericalHarmonicsGravityPartial > body2OnBody1Partial =
                        std::make_shared< SphericalHarmonicsGravityPartial >(
                                "Body1",
                                "Body2",
                                mutualSphericalModel->getAccelerationModelFromShExpansionOfBodyExertingAcceleration( ) );
                std::shared_ptr< SphericalHarmonicsGravityPartial > body1OnBody2Partial =
                        std::make_shared< SphericalHarmonicsGravityPartial >(
                                "Body2",
                                "Body1",
                                mutualSphericalModel->getAccelerationModelFromShExpansionOfBodyUndergoingAcceleration( ) );
                equivalentPartial = std::make_shared< MutualSphericalHarmonicsGravityPartial >(
                        body2OnBody1Partial, body1OnBody2Partial, "Body1", "Body2", mutualSphericalModel );
            }
            BOOST_REQUIRE( equivalentModel != nullptr );
            BOOST_REQUIRE( equivalentPartial != nullptr );

            equivalentModel->updateMembers( evaluationTime );
            equivalentPartial->update( evaluationTime );

            Eigen::MatrixXd equivalentPartialWrtBody1Position = Eigen::MatrixXd::Zero( 3, 3 );
            Eigen::MatrixXd equivalentPartialWrtBody2Position = Eigen::MatrixXd::Zero( 3, 3 );

            if( testCase.equivalentModelType == sphericalBody1OnBody2Equivalent )
            {
                equivalentPartial->wrtPositionOfAcceleratingBody(
                        equivalentPartialWrtBody1Position.block( 0, 0, 3, 3 ) );
                equivalentPartial->wrtPositionOfAcceleratedBody(
                        equivalentPartialWrtBody2Position.block( 0, 0, 3, 3 ) );
            }
            else
            {
                equivalentPartial->wrtPositionOfAcceleratedBody(
                        equivalentPartialWrtBody1Position.block( 0, 0, 3, 3 ) );
                equivalentPartial->wrtPositionOfAcceleratingBody(
                        equivalentPartialWrtBody2Position.block( 0, 0, 3, 3 ) );
            }

            equivalentPartialWrtBody1Position *= scaleToMutual;
            equivalentPartialWrtBody2Position *= scaleToMutual;

            BOOST_CHECK_GT( mutualPartialWrtBody1Position.norm( ), 1.0E-12 );
            BOOST_CHECK_GT( mutualPartialWrtBody2Position.norm( ), 1.0E-12 );
            BOOST_CHECK_GT( equivalentPartialWrtBody1Position.norm( ), 1.0E-12 );
            BOOST_CHECK_GT( equivalentPartialWrtBody2Position.norm( ), 1.0E-12 );

            // Position Jacobians from the full-two-body model must collapse to the equivalent model Jacobians.
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    mutualPartialWrtBody1Position, equivalentPartialWrtBody1Position, 2.0E-14 );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    mutualPartialWrtBody2Position, equivalentPartialWrtBody2Position, 2.0E-14 );

            if( testCase.compareBody1Coefficients )
            {
                Eigen::MatrixXd equivalentPartialWrtBody1Cosine =
                        equivalentPartial->wrtParameter( body1CosineCoefficientsParameter ) * scaleToMutual;
                Eigen::MatrixXd equivalentPartialWrtBody1Sine =
                        equivalentPartial->wrtParameter( body1SineCoefficientsParameter ) * scaleToMutual;

                BOOST_CHECK_GT( mutualPartialWrtBody1Cosine.norm( ), 1.0E-16 );
                BOOST_CHECK_GT( mutualPartialWrtBody1Sine.norm( ), 1.0E-16 );
                BOOST_CHECK_GT( equivalentPartialWrtBody1Cosine.norm( ), 1.0E-16 );
                BOOST_CHECK_GT( equivalentPartialWrtBody1Sine.norm( ), 1.0E-16 );

                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        mutualPartialWrtBody1Cosine, equivalentPartialWrtBody1Cosine, 2.0E-14 );
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        mutualPartialWrtBody1Sine, equivalentPartialWrtBody1Sine, 2.0E-14 );
            }
            else
            {
                BOOST_CHECK_SMALL( mutualPartialWrtBody1Cosine.norm( ), 1.0E-20 );
                BOOST_CHECK_SMALL( mutualPartialWrtBody1Sine.norm( ), 1.0E-20 );
            }

            if( testCase.compareBody2Coefficients )
            {
                Eigen::MatrixXd equivalentPartialWrtBody2Cosine =
                        equivalentPartial->wrtParameter( body2CosineCoefficientsParameter ) * scaleToMutual;
                Eigen::MatrixXd equivalentPartialWrtBody2Sine =
                        equivalentPartial->wrtParameter( body2SineCoefficientsParameter ) * scaleToMutual;

                BOOST_CHECK_GT( mutualPartialWrtBody2Cosine.norm( ), 1.0E-16 );
                BOOST_CHECK_GT( mutualPartialWrtBody2Sine.norm( ), 1.0E-16 );
                BOOST_CHECK_GT( equivalentPartialWrtBody2Cosine.norm( ), 1.0E-16 );
                BOOST_CHECK_GT( equivalentPartialWrtBody2Sine.norm( ), 1.0E-16 );

                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        mutualPartialWrtBody2Cosine, equivalentPartialWrtBody2Cosine, 2.0E-14 );
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        mutualPartialWrtBody2Sine, equivalentPartialWrtBody2Sine, 2.0E-14 );
            }
            else
            {
                BOOST_CHECK_SMALL( mutualPartialWrtBody2Cosine.norm( ), 1.0E-20 );
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
            gravitationalParameterBody1,
            equatorialRadiusBody1,
            cosineCoefficientsOfBody1,
            sineCoefficientsOfBody1,
            "IAU_Body1" ) );
    body2->setGravityFieldModel( std::make_shared< SphericalHarmonicsGravityField >(
            gravitationalParameterBody2,
            equatorialRadiusBody2,
            cosineCoefficientsOfBody2,
            sineCoefficientsOfBody2,
            "IAU_Body2" ) );
    centralBody->setGravityFieldModel( std::make_shared< SphericalHarmonicsGravityField >(
            gravitationalParameterCentralBody,
            equatorialRadiusCentralBody,
            cosineCoefficientsOfCentralBody,
            sineCoefficientsOfCentralBody,
            "IAU_CentralBody" ) );

    const double evaluationTime = 1000.0;
    body1->setRotationalEphemeris( std::make_shared< ephemerides::SimpleRotationalEphemeris >(
            Eigen::Quaterniond::Identity( ), 0.0, evaluationTime, "ECLIPJ2000", "IAU_Body1" ) );
    body2->setRotationalEphemeris( std::make_shared< ephemerides::SimpleRotationalEphemeris >(
            Eigen::Quaterniond(
                    Eigen::AngleAxisd( 0.3, Eigen::Vector3d::UnitZ( ) ) *
                    Eigen::AngleAxisd( -0.2, Eigen::Vector3d::UnitX( ) ) ),
            0.0,
            evaluationTime,
            "ECLIPJ2000",
            "IAU_Body2" ) );
    centralBody->setRotationalEphemeris( std::make_shared< ephemerides::SimpleRotationalEphemeris >(
            Eigen::Quaterniond(
                    Eigen::AngleAxisd( -0.4, Eigen::Vector3d::UnitY( ) ) *
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
            std::make_shared< FullTwoBodySphericalHarmonicAccelerationSettings >(
                    getExtendedSinglePointMassInteractions( 2, 2, 2, 2 ) );
    std::shared_ptr< ThirdBodyFullTwoBodySphericalHarmonicsGravitationalAccelerationModel > mutualExtendedModel =
            std::dynamic_pointer_cast< ThirdBodyFullTwoBodySphericalHarmonicsGravitationalAccelerationModel >(
                    createAccelerationModel(
                            body1, body2, mutualExtendedSettings, "Body1", "Body2", centralBody, "CentralBody", bodies ) );
    BOOST_REQUIRE( mutualExtendedModel != nullptr );

    mutualExtendedModel->updateMembers( evaluationTime );
    BOOST_CHECK_GT( mutualExtendedModel->getAcceleration( ).norm( ), 1.0E-16 );

    std::vector< std::shared_ptr< EstimatableParameter< double > > > doubleParameters;
    std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > > vectorParameters;
    std::shared_ptr< EstimatableParameterSet< double > > parameterSet =
            std::make_shared< EstimatableParameterSet< double > >( doubleParameters, vectorParameters );

    std::shared_ptr< AccelerationPartial > mutualExtendedPartial = createAnalyticalAccelerationPartial(
            mutualExtendedModel, std::make_pair( "Body1", body1 ), std::make_pair( "Body2", body2 ), bodies, parameterSet );

    BOOST_REQUIRE( mutualExtendedPartial != nullptr );

    mutualExtendedPartial->update( evaluationTime );

    Eigen::MatrixXd partialWrtBody1PositionExtended = Eigen::MatrixXd::Zero( 3, 3 );
    Eigen::MatrixXd partialWrtBody2PositionExtended = Eigen::MatrixXd::Zero( 3, 3 );
    Eigen::MatrixXd partialWrtCentralBodyPositionExtended = Eigen::MatrixXd::Zero( 3, 3 );

    mutualExtendedPartial->wrtPositionOfAcceleratedBody( partialWrtBody1PositionExtended.block( 0, 0, 3, 3 ) );
    mutualExtendedPartial->wrtPositionOfAcceleratingBody( partialWrtBody2PositionExtended.block( 0, 0, 3, 3 ) );
    mutualExtendedPartial->wrtPositionOfAdditionalBody(
            "CentralBody", partialWrtCentralBodyPositionExtended.block( 0, 0, 3, 3 ) );

    BOOST_CHECK_GT( partialWrtBody1PositionExtended.norm( ), 1.0E-16 );
    BOOST_CHECK_GT( partialWrtBody2PositionExtended.norm( ), 1.0E-16 );
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

    BOOST_CHECK_GT( partialWrtBody1CosineExtended.norm( ), 1.0E-16 );
    BOOST_CHECK_GT( partialWrtBody1SineExtended.norm( ), 1.0E-16 );
    BOOST_CHECK_GT( partialWrtBody2CosineExtended.norm( ), 1.0E-16 );
    BOOST_CHECK_GT( partialWrtBody2SineExtended.norm( ), 1.0E-16 );
    BOOST_CHECK_GT( partialWrtCentralCosineExtended.norm( ), 1.0E-16 );
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

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( numericalPartialWrtBody1Position, partialWrtBody1PositionExtended, 5.0E-5 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( numericalPartialWrtBody2Position, partialWrtBody2PositionExtended, 5.0E-5 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( numericalPartialWrtCentralPosition, partialWrtCentralBodyPositionExtended, 5.0E-5 );

    Eigen::MatrixXd numericalPartialWrtBody1Cosine = calculateAccelerationWrtParameterPartials(
            body1CosineCoefficientsParameter,
            mutualExtendedModel,
            Eigen::VectorXd::Constant( body1CosineCoefficientsParameter->getParameterSize( ), 1.0E-8 ) );
    Eigen::MatrixXd numericalPartialWrtBody1Sine = calculateAccelerationWrtParameterPartials(
            body1SineCoefficientsParameter,
            mutualExtendedModel,
            Eigen::VectorXd::Constant( body1SineCoefficientsParameter->getParameterSize( ), 1.0E-8 ) );
    Eigen::MatrixXd numericalPartialWrtBody2Cosine = calculateAccelerationWrtParameterPartials(
            body2CosineCoefficientsParameter,
            mutualExtendedModel,
            Eigen::VectorXd::Constant( body2CosineCoefficientsParameter->getParameterSize( ), 1.0E-8 ) );
    Eigen::MatrixXd numericalPartialWrtBody2Sine = calculateAccelerationWrtParameterPartials(
            body2SineCoefficientsParameter,
            mutualExtendedModel,
            Eigen::VectorXd::Constant( body2SineCoefficientsParameter->getParameterSize( ), 1.0E-8 ) );
    Eigen::MatrixXd numericalPartialWrtCentralCosine = calculateAccelerationWrtParameterPartials(
            centralCosineCoefficientsParameter,
            mutualExtendedModel,
            Eigen::VectorXd::Constant( centralCosineCoefficientsParameter->getParameterSize( ), 1.0E-8 ) );
    Eigen::MatrixXd numericalPartialWrtCentralSine = calculateAccelerationWrtParameterPartials(
            centralSineCoefficientsParameter,
            mutualExtendedModel,
            Eigen::VectorXd::Constant( centralSineCoefficientsParameter->getParameterSize( ), 1.0E-8 ) );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( numericalPartialWrtBody1Cosine, partialWrtBody1CosineExtended, 1.0E-3 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( numericalPartialWrtBody1Sine, partialWrtBody1SineExtended, 1.0E-3 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( numericalPartialWrtBody2Cosine, partialWrtBody2CosineExtended, 1.0E-3 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( numericalPartialWrtBody2Sine, partialWrtBody2SineExtended, 1.0E-3 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( numericalPartialWrtCentralCosine, partialWrtCentralCosineExtended, 1.0E-3 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( numericalPartialWrtCentralSine, partialWrtCentralSineExtended, 1.0E-3 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
