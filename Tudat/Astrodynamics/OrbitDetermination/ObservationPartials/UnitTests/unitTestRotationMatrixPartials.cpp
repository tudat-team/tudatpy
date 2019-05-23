/*    Copyright (c) 2010-2018, Delft University of Technology
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
#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>
#include <boost/lambda/lambda.hpp>

#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/SimulationSetup/EstimationSetup/createCartesianStatePartials.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/rotationMatrixPartial.h"
#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"
#include "Tudat/InputOutput/basicInputOutput.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::ephemerides;
using namespace tudat::estimatable_parameters;
using namespace tudat::observation_partials;

BOOST_AUTO_TEST_SUITE( test_rotation_matrix_partaisl )

//! Test whether partial derivatives of rotation matrix computed by SimpleRotationalEphemeris works correctly
BOOST_AUTO_TEST_CASE( testSimpleRotationalEphemerisPartials )
{
    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Create rotation model
    double nominalRotationRate = 2.0 * mathematical_constants::PI / 86400.0;
    std::shared_ptr< SimpleRotationalEphemeris > rotationalEphemeris =
            std::make_shared< SimpleRotationalEphemeris >(
                spice_interface::computeRotationQuaternionBetweenFrames( "ECLIPJ2000", "IAU_Earth", 1.0E7 ),
                nominalRotationRate, 1.0E7, "ECLIPJ2000", "IAU_Earth" );

    {
        // Create partial object.
        std::shared_ptr< RotationMatrixPartialWrtConstantRotationRate > rotationMatrixPartialObject =
                std::make_shared< RotationMatrixPartialWrtConstantRotationRate >( rotationalEphemeris );

        // Compute partial analytically
        double testTime = 1.0E6;
        Eigen::Matrix3d rotationMatrixPartial =
                rotationMatrixPartialObject->calculatePartialOfRotationMatrixToBaseFrameWrParameter(
                    testTime ).at( 0 );

        Eigen::Matrix3d rotationMatrixDerivativePartial =
                rotationMatrixPartialObject->calculatePartialOfRotationMatrixDerivativeToBaseFrameWrParameter(
                    testTime ).at( 0 );

        // Compute partial numerically.
        double perturbation = 1.0E-12;
        rotationalEphemeris->resetRotationRate( nominalRotationRate + perturbation );
        Eigen::Matrix3d upperturbedRotationMatrix = rotationalEphemeris->getRotationToBaseFrame(
                    testTime).toRotationMatrix( );
        Eigen::Matrix3d upperturbedRotationMatrixDerivative = rotationalEphemeris->getDerivativeOfRotationToBaseFrame(
                    testTime );

        rotationalEphemeris->resetRotationRate( nominalRotationRate - perturbation );
        Eigen::Matrix3d downperturbedRotationMatrix = rotationalEphemeris->getRotationToBaseFrame(
                    testTime).toRotationMatrix( );
        Eigen::Matrix3d downperturbedRotationMatrixDerivative = rotationalEphemeris->getDerivativeOfRotationToBaseFrame(
                    testTime );

        Eigen::Matrix3d numericalRotationMatrixPartial =
                ( upperturbedRotationMatrix - downperturbedRotationMatrix ) / ( 2.0 * perturbation );
        Eigen::Matrix3d numericalRotationMatrixDerivativePartial =
                ( upperturbedRotationMatrixDerivative - downperturbedRotationMatrixDerivative ) / ( 2.0 * perturbation );

        Eigen::Matrix3d matrixDifference = rotationMatrixPartial - numericalRotationMatrixPartial;

        // Compare analytical and numerical result.
        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL( std::fabs( matrixDifference( i, j ) ), 0.1 );
            }
        }

        matrixDifference = rotationMatrixDerivativePartial - numericalRotationMatrixDerivativePartial;

        // Compare analytical and numerical result.
        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL( std::fabs( matrixDifference( i, j ) ), 1.0E-5 );
            }
        }
    }

    {
        // Create partial object.
        std::shared_ptr< RotationMatrixPartialWrtPoleOrientation > rotationMatrixPartialObject =
                std::make_shared< RotationMatrixPartialWrtPoleOrientation >( rotationalEphemeris );

        // Compute partial analytically
        double testTime = 1.0E6;
        std::vector< Eigen::Matrix3d > rotationMatrixPartials =
                rotationMatrixPartialObject->calculatePartialOfRotationMatrixToBaseFrameWrParameter(
                    testTime );

        std::vector< Eigen::Matrix3d > rotationMatrixDerivativePartials =
                rotationMatrixPartialObject->calculatePartialOfRotationMatrixDerivativeToBaseFrameWrParameter(
                    testTime );

        Eigen::Vector3d nominalEulerAngles = rotationalEphemeris->getInitialEulerAngles( );
        double perturbedAngle;

        // Compute partial numerically.
        double perturbation = 1.0E-6;
        {


            // Compute partial for right ascension numerically.
            {
                perturbedAngle = nominalEulerAngles( 0 ) + perturbation;
                rotationalEphemeris->resetInitialPoleRightAscensionAndDeclination( perturbedAngle, nominalEulerAngles( 1 ) );
                Eigen::Matrix3d upperturbedRotationMatrix = rotationalEphemeris->getRotationToBaseFrame(
                            testTime).toRotationMatrix( );
                Eigen::Matrix3d upperturbedRotationMatrixDerivative = rotationalEphemeris->getDerivativeOfRotationToBaseFrame(
                            testTime );

                perturbedAngle = nominalEulerAngles( 0 ) - perturbation;
                rotationalEphemeris->resetInitialPoleRightAscensionAndDeclination( perturbedAngle, nominalEulerAngles( 1 ) );
                Eigen::Matrix3d downperturbedRotationMatrix = rotationalEphemeris->getRotationToBaseFrame(
                            testTime).toRotationMatrix( );
                Eigen::Matrix3d downperturbedRotationMatrixDerivative =
                        rotationalEphemeris->getDerivativeOfRotationToBaseFrame(
                            testTime );

                Eigen::Matrix3d numericalRotationMatrixPartial =
                        ( upperturbedRotationMatrix - downperturbedRotationMatrix ) / ( 2.0 * perturbation );
                Eigen::Matrix3d numericalRotationMatrixDerivativePartial =
                        ( upperturbedRotationMatrixDerivative - downperturbedRotationMatrixDerivative ) /
                        ( 2.0 * perturbation );

                Eigen::Matrix3d matrixDifference = rotationMatrixPartials.at( 0 ) - numericalRotationMatrixPartial;

                // Compare analytical and numerical result.
                for( unsigned int i = 0; i < 3; i++ )
                {
                    for( unsigned int j = 0; j < 3; j++ )
                    {
                        BOOST_CHECK_SMALL( std::fabs( matrixDifference( i, j ) ), 1.0E-8 );
                    }
                }

                matrixDifference = rotationMatrixDerivativePartials.at( 0 ) - numericalRotationMatrixDerivativePartial;

                // Compare analytical and numerical result.
                for( unsigned int i = 0; i < 3; i++ )
                {
                    for( unsigned int j = 0; j < 3; j++ )
                    {
                        BOOST_CHECK_SMALL( std::fabs( matrixDifference( i, j ) ), 1.0E-13 );
                    }
                }
            }

            // Compute partial for declination numerically.
            {
                perturbedAngle = nominalEulerAngles( 1 ) + perturbation;
                rotationalEphemeris->resetInitialPoleRightAscensionAndDeclination( nominalEulerAngles( 0 ), perturbedAngle );
                Eigen::Matrix3d upperturbedRotationMatrix = rotationalEphemeris->getRotationToBaseFrame(
                            testTime).toRotationMatrix( );
                Eigen::Matrix3d upperturbedRotationMatrixDerivative = rotationalEphemeris->getDerivativeOfRotationToBaseFrame(
                            testTime );

                perturbedAngle = nominalEulerAngles( 1 ) - perturbation;
                rotationalEphemeris->resetInitialPoleRightAscensionAndDeclination( nominalEulerAngles( 0 ), perturbedAngle );
                Eigen::Matrix3d downperturbedRotationMatrix = rotationalEphemeris->getRotationToBaseFrame(
                            testTime).toRotationMatrix( );
                Eigen::Matrix3d downperturbedRotationMatrixDerivative =
                        rotationalEphemeris->getDerivativeOfRotationToBaseFrame( testTime );

                Eigen::Matrix3d numericalRotationMatrixPartial =
                        ( upperturbedRotationMatrix - downperturbedRotationMatrix ) / ( 2.0 * perturbation );
                Eigen::Matrix3d numericalRotationMatrixDerivativePartial =
                        ( upperturbedRotationMatrixDerivative -
                          downperturbedRotationMatrixDerivative ) / ( 2.0 * perturbation );

                Eigen::Matrix3d matrixDifference = rotationMatrixPartials.at( 1 ) - numericalRotationMatrixPartial;

                // Compare analytical and numerical result.
                for( unsigned int i = 0; i < 3; i++ )
                {
                    for( unsigned int j = 0; j < 3; j++ )
                    {
                        BOOST_CHECK_SMALL( std::fabs( matrixDifference( i, j ) ), 1.0E-8 );
                    }
                }

                matrixDifference = rotationMatrixDerivativePartials.at( 1 ) - numericalRotationMatrixDerivativePartial;

                // Compare analytical and numerical result.
                for( unsigned int i = 0; i < 3; i++ )
                {
                    for( unsigned int j = 0; j < 3; j++ )
                    {
                        BOOST_CHECK_SMALL( std::fabs( matrixDifference( i, j ) ), 1.0E-13 );
                    }
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( testTidallyLockedRotationPartials )
{
    Eigen::Vector6d currentState =
            tudat::spice_interface::getBodyCartesianStateAtEpoch(
                "Mercury", "SSB", "J2000", "None", 0.0 );
    Eigen::Vector3d positionVector = currentState.segment( 0, 3 );
    Eigen::Vector3d velocityVector = currentState.segment( 3, 3 );

    double positionNorm = positionVector.norm( );

    Eigen::Matrix3d currentRotationMatrix = tudat::reference_frames::getInertialToRswSatelliteCenteredFrameRotationMatrix(
                currentState ).transpose( );

    Eigen::Vector3d rVector = currentRotationMatrix.block( 0, 0, 3, 1 );
    Eigen::Vector3d sVector = currentRotationMatrix.block( 0, 1, 3, 1 );
    Eigen::Vector3d wVector = currentRotationMatrix.block( 0, 2, 3, 1 );
    Eigen::Vector3d unnormalizedWVector = currentState.segment< 3 >( 0 ).cross(
                currentState.segment< 3 >( 3 ) );
    double unnormalizedWVectorNorm = unnormalizedWVector.norm( );

    std::cout<<"pos norm "<<positionNorm<<std::endl;

    Eigen::Matrix3d rVectorDerivativeWrtPosition =
            Eigen::Matrix3d::Identity( ) / positionNorm - positionVector * positionVector.transpose( ) / (
                positionNorm * positionNorm * positionNorm );
    Eigen::Matrix3d unnormalizedWVectorDerivativeWrtPosition =
            -linear_algebra::getCrossProductMatrix( currentState.segment( 3, 3 ) );
    Eigen::Matrix3d unnormalizedWVectorDerivativeWrtVelocity =
            linear_algebra::getCrossProductMatrix( positionVector );
    Eigen::Matrix3d wPartialScalingTerm =
            ( Eigen::Matrix3d::Identity( ) / unnormalizedWVectorNorm -
              unnormalizedWVector * unnormalizedWVector.transpose( ) /
              ( unnormalizedWVectorNorm * unnormalizedWVectorNorm * unnormalizedWVectorNorm ) );

    Eigen::Matrix3d wVectorDerivativeWrtPosition =
            wPartialScalingTerm * unnormalizedWVectorDerivativeWrtPosition;
    Eigen::Matrix3d wVectorDerivativeWrtVelocity =
            wPartialScalingTerm * unnormalizedWVectorDerivativeWrtVelocity;

    Eigen::Matrix3d sVectorDerivativeWrtPosition =
            linear_algebra::getCrossProductMatrix( wVector ) * rVectorDerivativeWrtPosition -
            linear_algebra::getCrossProductMatrix( rVector ) * wVectorDerivativeWrtPosition;
    Eigen::Matrix3d sVectorDerivativeWrtVelocity =
            -linear_algebra::getCrossProductMatrix( positionVector ) * wVectorDerivativeWrtVelocity;

    double positionPerturbation = 10000.0;
    double velocityPerturbation = 1.0;
    Eigen::Vector6d nominalState = currentState;
    Eigen::Matrix3d nominalRotationMatrix = currentRotationMatrix;

    std::cout<<"R partials "<<std::endl<<
               rVectorDerivativeWrtPosition<<std::endl<<std::endl;
    std::cout<<"S partials "<<std::endl<<
               sVectorDerivativeWrtPosition<<std::endl<<std::endl;
    std::cout<<"W partials "<<std::endl<<
               wVectorDerivativeWrtPosition<<std::endl<<std::endl;
    std::cout<<"Unnormalized W partials "<<std::endl<<
               unnormalizedWVectorDerivativeWrtPosition<<std::endl<<std::endl;
    for( int i = 0; i < 3; i++ )
    {
        currentState = nominalState;
        currentState( i ) += positionPerturbation;
        Eigen::Matrix3d upPerturbedRotationMatrix =
                tudat::reference_frames::getInertialToRswSatelliteCenteredFrameRotationMatrix(
                                currentState ).transpose( );
        Eigen::Vector3d upperturbedUnnormalizedWVector = currentState.segment< 3 >( 0 ).cross(
                    currentState.segment< 3 >( 3 ) );

        currentState = nominalState;
        currentState( i ) -= positionPerturbation;
        Eigen::Matrix3d downPerturbedRotationMatrix =
                tudat::reference_frames::getInertialToRswSatelliteCenteredFrameRotationMatrix(
                                currentState ).transpose( );
        Eigen::Vector3d downperturbedUnnormalizedWVector = currentState.segment< 3 >( 0 ).cross(
                    currentState.segment< 3 >( 3 ) );


        std::cout<<"Numerical partial"<<std::endl<<
                   ( upPerturbedRotationMatrix - downPerturbedRotationMatrix ) / ( 2.0 * positionPerturbation )
                   <<std::endl<<std::endl;
//        std::cout<<"Unnormalized w partial"<<std::endl<<
//                   ( upperturbedUnnormalizedWVector - downperturbedUnnormalizedWVector ) / ( 2.0 * positionPerturbation )
//                   <<std::endl<<std::endl;
    }

}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat





