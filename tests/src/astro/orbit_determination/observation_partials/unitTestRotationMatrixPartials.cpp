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
#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>
#include <boost/lambda/lambda.hpp>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/simulation/estimation_setup/createCartesianStatePartials.h"
#include "tudat/astro/orbit_determination/observation_partials/rotationMatrixPartial.h"
#include "tudat/astro/reference_frames/referenceFrameTransformations.h"
#include "tudat/io/basicInputOutput.h"

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

//! Test whether partial derivatives of rotation matrix computed by SynchronousRotationalEphemeris works correctly
BOOST_AUTO_TEST_CASE( testSynchronousRotationPartials )
{
    // Define nominal state
    Eigen::Vector6d nominalState =
            tudat::spice_interface::getBodyCartesianStateAtEpoch(
                                       "Mercury", "SSB", "ECLIPJ2000", "None", 1.0E7 );

    // Define nominal state function
    Eigen::Vector6d currentState = nominalState;
    std::function< Eigen::Vector6d( const double, bool ) > relativeStateFunction =
            [ & ]( const double, bool ){ return currentState; };

    // Create rotation model
    std::shared_ptr< tudat::ephemerides::SynchronousRotationalEphemeris > synchronousRotationModel =
            std::make_shared< ephemerides::SynchronousRotationalEphemeris >(
                relativeStateFunction, "SSB", "Mercury_Fixed", "ECLIPJ2000" );

    // Create rotation partial model
    std::shared_ptr< RotationMatrixPartial > rotationMatrixPartialObject =
            std::make_shared< SynchronousRotationMatrixPartialWrtTranslationalState >( synchronousRotationModel );

    // Define test settings
    double testTime = 1.0E7;
    double positionPerturbation = 10000.0;
    double velocityPerturbation = 0.1;

    // Test partials w.r.t. position and velocity components
    std::vector< Eigen::Matrix3d > rotationMatrixPartials =
            rotationMatrixPartialObject->calculatePartialOfRotationMatrixToBaseFrameWrParameter( testTime );
    for( int i = 0; i < 3; i++ )
    {
        currentState = nominalState;
        currentState( i ) += positionPerturbation;
        Eigen::Matrix3d upPerturbedRotationMatrix = synchronousRotationModel->getRotationToBaseFrame(
                    1.0E7 ).toRotationMatrix( );

        currentState = nominalState;
        currentState( i ) -= positionPerturbation;
        Eigen::Matrix3d downPerturbedRotationMatrix = synchronousRotationModel->getRotationToBaseFrame(
                    1.0E7 ).toRotationMatrix( );

        Eigen::Matrix3d relativePartialError =
                ( ( upPerturbedRotationMatrix - downPerturbedRotationMatrix ) / ( 2.0 * positionPerturbation ) -
                rotationMatrixPartials.at( i ) ) / rotationMatrixPartials.at( i ).norm( );

        for( int j = 0; j < 3; j++ )
        {
            for( int k = 0; k < 3; k++ )
            {
                BOOST_CHECK_SMALL( std::fabs( relativePartialError( j, k ) ), 1.0E-8 );
            }
        }

        currentState = nominalState;
        currentState( i + 3 ) += velocityPerturbation;
        upPerturbedRotationMatrix = synchronousRotationModel->getRotationToBaseFrame(
                    1.0E7 ).toRotationMatrix( );

        currentState = nominalState;
        currentState( i + 3 ) -= velocityPerturbation;
        downPerturbedRotationMatrix = synchronousRotationModel->getRotationToBaseFrame(
                    1.0E7 ).toRotationMatrix( );

        relativePartialError =
                        ( ( upPerturbedRotationMatrix - downPerturbedRotationMatrix ) / ( 2.0 * velocityPerturbation ) -
                        rotationMatrixPartials.at( i + 3 ) ) / rotationMatrixPartials.at( i + 3 ).norm( );

        for( int j = 0; j < 3; j++ )
        {
            for( int k = 0; k < 3; k++ )
            {
                BOOST_CHECK_SMALL( std::fabs( relativePartialError( j, k ) ), 1.0E-8 );
            }
        }
    }
}


//! Test whether partial derivatives of rotation matrix computed by SynchronousRotationalEphemeris works correctly
BOOST_AUTO_TEST_CASE( testIauRotationPartials )
{
    spice_interface::loadSpiceKernelInTudat( paths::getSpiceKernelPath( ) + "/pck00010.tpc" );

    std::string baseFrameOrientation = "J2000";
    std::string targetFrameOrientation = "IAU_Jupiter";

    double degreeToRadian = unit_conversions::convertDegreesToRadians( 1.0 );
    double nominalMeridian = 284.95 * degreeToRadian;
    Eigen::Vector2d nominalPole = ( Eigen::Vector2d( ) << 268.056595, 64.495303 ).finished( ) * degreeToRadian;
    double rotationRate = 870.5360000 * degreeToRadian / physical_constants::JULIAN_DAY;
    Eigen::Vector2d polePrecession = ( Eigen::Vector2d( ) << -0.006499, 0.002413 ).finished( ) * degreeToRadian / physical_constants::JULIAN_CENTURY;


    std::map< double, std::pair< double, double > > meridianPeriodicTerms;
    std::map< double, std::pair< Eigen::Vector2d, double > > polePeriodicTerms;
    polePeriodicTerms[ 4850.4046 * degreeToRadian / physical_constants::JULIAN_CENTURY ] =
            std::make_pair( degreeToRadian * ( Eigen::Vector2d( ) << 0.000117, 0.000050 ).finished( ), degreeToRadian * 99.360714 );
    polePeriodicTerms[ 1191.9605 * degreeToRadian / physical_constants::JULIAN_CENTURY ] =
            std::make_pair( degreeToRadian * ( Eigen::Vector2d( ) << 0.000938, 0.000404 ).finished( ), degreeToRadian * 175.895369 );
    polePeriodicTerms[ 262.5475 * degreeToRadian / physical_constants::JULIAN_CENTURY ] =
            std::make_pair( degreeToRadian * ( Eigen::Vector2d( ) << 0.001432, 0.000617 ).finished( ), degreeToRadian * 300.323162 );
    polePeriodicTerms[ 6070.2476 * degreeToRadian / physical_constants::JULIAN_CENTURY ] =
            std::make_pair( degreeToRadian * ( Eigen::Vector2d( ) << 0.000030, -0.000013 ).finished( ), degreeToRadian * 114.012305 );
    polePeriodicTerms[ 64.3000 * degreeToRadian / physical_constants::JULIAN_CENTURY  ] =
            std::make_pair( degreeToRadian * ( Eigen::Vector2d( ) << 0.002150, 0.000926 ).finished( ), degreeToRadian * 49.511251 );


    // Create rotation model
    std::shared_ptr< tudat::ephemerides::IauRotationModel > iauRotationModel =
            std::make_shared< ephemerides::IauRotationModel >(
                baseFrameOrientation, targetFrameOrientation,
                nominalMeridian, nominalPole, rotationRate, polePrecession,
                meridianPeriodicTerms, polePeriodicTerms );

    {
        // Create partial object.
        std::shared_ptr< RotationMatrixPartialWrtNominalPolePosition > rotationMatrixPartialObject =
                std::make_shared< RotationMatrixPartialWrtNominalPolePosition >( iauRotationModel );

        // Compute partial analytically
        double testTime = 1.0E9;
        std::vector< Eigen::Matrix3d > rotationMatrixPartials =
                rotationMatrixPartialObject->calculatePartialOfRotationMatrixToBaseFrameWrParameter(
                    testTime );

//        std::vector< Eigen::Matrix3d > rotationMatrixDerivativePartials =
//                rotationMatrixPartialObject->calculatePartialOfRotationMatrixDerivativeToBaseFrameWrParameter(
//                    testTime );

        Eigen::Vector2d unperturbedPole = iauRotationModel->getNominalPole( );
        Eigen::Vector2d perturbedPole = unperturbedPole;

        // Compute partial numerically.
        double perturbation = 1.0E-6;
        {
            // Compute partial for right ascension numerically.
            {
                perturbedPole( 0 ) = unperturbedPole( 0 ) + perturbation;
                iauRotationModel->setNominalPole( perturbedPole );
                Eigen::Matrix3d upperturbedRotationMatrix = iauRotationModel->getRotationToBaseFrame(
                            testTime).toRotationMatrix( );
                Eigen::Matrix3d upperturbedRotationMatrixDerivative = iauRotationModel->getDerivativeOfRotationToBaseFrame(
                            testTime );

                perturbedPole( 0 ) = unperturbedPole( 0 ) - perturbation;
                iauRotationModel->setNominalPole( perturbedPole );
                Eigen::Matrix3d downperturbedRotationMatrix = iauRotationModel->getRotationToBaseFrame(
                            testTime).toRotationMatrix( );
                Eigen::Matrix3d downperturbedRotationMatrixDerivative =
                        iauRotationModel->getDerivativeOfRotationToBaseFrame(
                            testTime );
                perturbedPole = unperturbedPole;

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

//                matrixDifference = rotationMatrixDerivativePartials.at( 0 ) - numericalRotationMatrixDerivativePartial;

//                // Compare analytical and numerical result.
//                for( unsigned int i = 0; i < 3; i++ )
//                {
//                    for( unsigned int j = 0; j < 3; j++ )
//                    {
//                        BOOST_CHECK_SMALL( std::fabs( matrixDifference( i, j ) ), 1.0E-13 );
//                    }
//                }
            }

            // Compute partial for declination numerically.
            {
                perturbedPole( 1 ) = unperturbedPole( 1 ) + perturbation;
                iauRotationModel->setNominalPole( perturbedPole );
                Eigen::Matrix3d upperturbedRotationMatrix = iauRotationModel->getRotationToBaseFrame(
                            testTime).toRotationMatrix( );
                Eigen::Matrix3d upperturbedRotationMatrixDerivative = iauRotationModel->getDerivativeOfRotationToBaseFrame(
                            testTime );

                perturbedPole( 1 ) = unperturbedPole( 1 ) - perturbation;
                iauRotationModel->setNominalPole( perturbedPole );
                Eigen::Matrix3d downperturbedRotationMatrix = iauRotationModel->getRotationToBaseFrame(
                            testTime).toRotationMatrix( );
                Eigen::Matrix3d downperturbedRotationMatrixDerivative =
                        iauRotationModel->getDerivativeOfRotationToBaseFrame( testTime );

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

//                matrixDifference = rotationMatrixDerivativePartials.at( 1 ) - numericalRotationMatrixDerivativePartial;

//                // Compare analytical and numerical result.
//                for( unsigned int i = 0; i < 3; i++ )
//                {
//                    for( unsigned int j = 0; j < 3; j++ )
//                    {
//                        BOOST_CHECK_SMALL( std::fabs( matrixDifference( i, j ) ), 1.0E-13 );
//                    }
//                }
            }
        }
    }

    {
        // Create partial object.
        std::shared_ptr< RotationMatrixPartialWrtPolePositionRate > rotationMatrixPartialObject =
                std::make_shared< RotationMatrixPartialWrtPolePositionRate >( iauRotationModel );

        // Compute partial analytically
        double testTime = 1.0E9;
        std::vector< Eigen::Matrix3d > rotationMatrixPartials =
                rotationMatrixPartialObject->calculatePartialOfRotationMatrixToBaseFrameWrParameter(
                    testTime );

//        std::vector< Eigen::Matrix3d > rotationMatrixDerivativePartials =
//                rotationMatrixPartialObject->calculatePartialOfRotationMatrixDerivativeToBaseFrameWrParameter(
//                    testTime );

        Eigen::Vector2d unperturbedPoleRate = iauRotationModel->getPolePrecession( );
        Eigen::Vector2d perturbedPoleRate = unperturbedPoleRate;

        std::cout<<unperturbedPoleRate<<std::endl;
        // Compute partial numerically.
        {
            // Compute partial for right ascension numerically.
            {
                double perturbation = 1.0E-14;
                perturbedPoleRate( 0 ) = unperturbedPoleRate( 0 ) + perturbation;
                iauRotationModel->setPolePrecession( perturbedPoleRate );
                Eigen::Matrix3d upperturbedRotationMatrix = iauRotationModel->getRotationToBaseFrame(
                            testTime).toRotationMatrix( );
                Eigen::Matrix3d upperturbedRotationMatrixDerivative = iauRotationModel->getDerivativeOfRotationToBaseFrame(
                            testTime );

                perturbedPoleRate( 0 ) = unperturbedPoleRate( 0 ) - perturbation;
                iauRotationModel->setPolePrecession( perturbedPoleRate );
                Eigen::Matrix3d downperturbedRotationMatrix = iauRotationModel->getRotationToBaseFrame(
                            testTime).toRotationMatrix( );
                Eigen::Matrix3d downperturbedRotationMatrixDerivative =
                        iauRotationModel->getDerivativeOfRotationToBaseFrame(
                            testTime );
                iauRotationModel->setPolePrecession( unperturbedPoleRate );

                Eigen::Matrix3d numericalRotationMatrixPartial =
                        ( upperturbedRotationMatrix - downperturbedRotationMatrix ) / ( 2.0 * perturbation );
                Eigen::Matrix3d numericalRotationMatrixDerivativePartial =
                        ( upperturbedRotationMatrixDerivative - downperturbedRotationMatrixDerivative ) /
                        ( 2.0 * perturbation );

                Eigen::Matrix3d matrixDifference = rotationMatrixPartials.at( 0 ) - numericalRotationMatrixPartial;

                std::cout<<"Matrices: "<<std::endl<<
                           numericalRotationMatrixPartial<<std::endl<<std::endl<<
                           rotationMatrixPartials.at( 0 )<<std::endl<<std::endl<<
                           matrixDifference<<std::endl<<std::endl;

                   // Compare analytical and numerical result.
                for( unsigned int i = 0; i < 3; i++ )
                {
                    for( unsigned int j = 0; j < 3; j++ )
                    {
                        BOOST_CHECK_SMALL( std::fabs( matrixDifference( i, j ) ), 0.1 );
                    }
                }

//                matrixDifference = rotationMatrixDerivativePartials.at( 0 ) - numericalRotationMatrixDerivativePartial;

//                // Compare analytical and numerical result.
//                for( unsigned int i = 0; i < 3; i++ )
//                {
//                    for( unsigned int j = 0; j < 3; j++ )
//                    {
//                        BOOST_CHECK_SMALL( std::fabs( matrixDifference( i, j ) ), 1.0E-13 );
//                    }
//                }
            }

            // Compute partial for declination numerically.
            {
                double perturbation = 1.0E-14;

                perturbedPoleRate = unperturbedPoleRate;
                perturbedPoleRate( 1 ) = unperturbedPoleRate( 1 ) + perturbation;
                iauRotationModel->setPolePrecession( perturbedPoleRate );
                Eigen::Matrix3d upperturbedRotationMatrix = iauRotationModel->getRotationToBaseFrame(
                            testTime).toRotationMatrix( );
                Eigen::Matrix3d upperturbedRotationMatrixDerivative = iauRotationModel->getDerivativeOfRotationToBaseFrame(
                            testTime );

                perturbedPoleRate( 1 ) = unperturbedPoleRate( 1 ) - perturbation;
                iauRotationModel->setPolePrecession( perturbedPoleRate );
                Eigen::Matrix3d downperturbedRotationMatrix = iauRotationModel->getRotationToBaseFrame(
                            testTime).toRotationMatrix( );
                Eigen::Matrix3d downperturbedRotationMatrixDerivative =
                        iauRotationModel->getDerivativeOfRotationToBaseFrame( testTime );
                perturbedPoleRate = unperturbedPoleRate;

                Eigen::Matrix3d numericalRotationMatrixPartial =
                        ( upperturbedRotationMatrix - downperturbedRotationMatrix ) / ( 2.0 * perturbation );
                Eigen::Matrix3d numericalRotationMatrixDerivativePartial =
                        ( upperturbedRotationMatrixDerivative -
                          downperturbedRotationMatrixDerivative ) / ( 2.0 * perturbation );

                Eigen::Matrix3d matrixDifference = rotationMatrixPartials.at( 1 ) - numericalRotationMatrixPartial;

                std::cout<<"Matrices: "<<std::endl<<
                           numericalRotationMatrixPartial<<std::endl<<std::endl<<
                           rotationMatrixPartials.at( 1 )<<std::endl<<std::endl<<
                           matrixDifference<<std::endl<<std::endl;

                // Compare analytical and numerical result.
                for( unsigned int i = 0; i < 3; i++ )
                {
                    for( unsigned int j = 0; j < 3; j++ )
                    {
                        BOOST_CHECK_SMALL( std::fabs( matrixDifference( i, j ) ), 0.1 );
                    }
                }

//                matrixDifference = rotationMatrixDerivativePartials.at( 1 ) - numericalRotationMatrixDerivativePartial;

//                // Compare analytical and numerical result.
//                for( unsigned int i = 0; i < 3; i++ )
//                {
//                    for( unsigned int j = 0; j < 3; j++ )
//                    {
//                        BOOST_CHECK_SMALL( std::fabs( matrixDifference( i, j ) ), 1.0E-13 );
//                    }
//                }
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat





