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
#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>
#include <boost/lambda/lambda.hpp>

#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/simulation/estimation_setup/createCartesianStatePartials.h"
#include "tudat/astro/orbit_determination/observation_partials/rotationMatrixPartial.h"
#include "tudat/astro/ephemerides/fullPlanetaryRotationModel.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/createRotationModel.h"
#include "tudat/io/basicInputOutput.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::ephemerides;
using namespace tudat::estimatable_parameters;
using namespace tudat::observation_partials;
using namespace tudat::simulation_setup;

BOOST_AUTO_TEST_SUITE( test_planetary_rotation_model_rotation_matrix_partials )

//! Test whether partial derivatives of rotation matrix computed by FullPlanetaryRotationMode works correctly
BOOST_AUTO_TEST_CASE( testPlanetaryRotationModelEphemerisPartials )
{
    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Create rotation model
    double initialTime = 0.0;
    double finalTime = 1.0E8;

    spice_interface::loadStandardSpiceKernels( );

    NamedBodyMap bodyMap;
    bodyMap.addNewBody( "Mars" );

    std::shared_ptr< RotationModelSettings > defaultMarsRotationSettings =
            getHighAccuracyMarsRotationModel( initialTime, finalTime );

    std::shared_ptr< RotationalEphemeris > marsRotationModel =
            createRotationModel( defaultMarsRotationSettings, "Mars" );

    std::shared_ptr< PlanetaryRotationModelSettings > marsRotationModelSettings =
            std::dynamic_pointer_cast< PlanetaryRotationModelSettings >( defaultMarsRotationSettings );

    double coreFactor = marsRotationModelSettings->getCoreFactor( );

    double freeCoreNutationRate = marsRotationModelSettings->getFreeCoreNutationRate( );

    std::map< double, std::pair< double, double > > rotationRateCorrections =
            marsRotationModelSettings -> getRotationRateCorrections( );

    std::map< double, std::pair< double, double > > xPolarMotionCoefficients =
            marsRotationModelSettings -> getxPolarMotionCoefficients( );

    std::map< double, std::pair< double, double > > yPolarMotionCoefficients =
            marsRotationModelSettings -> getyPolarMotionCoefficients( );

    //Test for Rotation Matrix Partial wrt Periodic Spin Variations
    {

        // Create partial object.

        std::shared_ptr< RotationMatrixPartialWrtPeriodicSpinVariations > rotationMatrixPartialObject =
                std::make_shared< RotationMatrixPartialWrtPeriodicSpinVariations >(
                    std::dynamic_pointer_cast< PlanetaryRotationModel >( marsRotationModel ) );

        // Compute partial analytically

        double testTime = 9.0E7;

        std::vector< Eigen::Matrix3d > rotationMatrixPartials =
                rotationMatrixPartialObject->calculatePartialOfRotationMatrixToBaseFrameWrParameter( testTime );

        std::vector< Eigen::Matrix3d > rotationMatrixDerivativePartials =
                rotationMatrixPartialObject->calculatePartialOfRotationMatrixDerivativeToBaseFrameWrParameter( testTime );

        // Compute partial numerically
        {
            // Compute partial for Phi_cj numerically
            {
//                std::cout << std::endl<< std::endl << "Phi_cj:" << std::endl;

                double perturbation = 1.0E-4;

                std::map< double, std::pair< double, double > > upperturbateRotationRateCorrections;

                std::map< double, std::pair< double, double > > downperturbateRotationRateCorrections;

                for( std::map< double, std::pair< double, double > >::iterator correctionIterator = rotationRateCorrections.begin( );
                     correctionIterator != rotationRateCorrections.end( ); correctionIterator++ )
                {
                    upperturbateRotationRateCorrections[ correctionIterator->first ] = std::make_pair(
                                correctionIterator->second.first + perturbation, correctionIterator->second.second );

                    downperturbateRotationRateCorrections[ correctionIterator->first ] = std::make_pair(
                                correctionIterator->second.first - perturbation, correctionIterator->second.second );

                    if ( correctionIterator->first > 1.0 )
                    {
                        upperturbateRotationRateCorrections[ correctionIterator->first ] = std::make_pair(
                                    correctionIterator->second.first, correctionIterator->second.second  );

                        downperturbateRotationRateCorrections[ correctionIterator->first ] = std::make_pair(
                                    correctionIterator->second.first, correctionIterator->second.second);
                    }
                }

                marsRotationModelSettings->setPeriodTerms(
                            marsRotationModelSettings->getMeanMotionDirectNutationCorrections( ),
                            marsRotationModelSettings->getMeanMotionTimeDependentPhaseNutationCorrections( ),
                            upperturbateRotationRateCorrections,
                            xPolarMotionCoefficients,
                            yPolarMotionCoefficients );

                std::shared_ptr< RotationalEphemeris > upperturbatedMarsRotationModel =
                        createRotationModel( marsRotationModelSettings, "Mars" );

                Eigen::Matrix3d upperturbedRotationMatrix =
                        upperturbatedMarsRotationModel->getRotationToBaseFrame( testTime ).toRotationMatrix( );

                Eigen::Matrix3d upperturbedRotationMatrixDerivative =
                        upperturbatedMarsRotationModel->getDerivativeOfRotationToBaseFrame( testTime );

                marsRotationModelSettings->setPeriodTerms(
                            marsRotationModelSettings->getMeanMotionDirectNutationCorrections( ),
                            marsRotationModelSettings->getMeanMotionTimeDependentPhaseNutationCorrections( ),
                            downperturbateRotationRateCorrections,
                            xPolarMotionCoefficients,
                            yPolarMotionCoefficients );

                std::shared_ptr< RotationalEphemeris > downperturbatedMarsRotationModel =
                        createRotationModel( marsRotationModelSettings, "Mars" );

                Eigen::Matrix3d downperturbedRotationMatrix =
                        downperturbatedMarsRotationModel->getRotationToBaseFrame( testTime ).toRotationMatrix( );

                Eigen::Matrix3d downperturbedRotationMatrixDerivative =
                        downperturbatedMarsRotationModel->getDerivativeOfRotationToBaseFrame( testTime );

                Eigen::Matrix3d numericalRotationMatrixPartial =
                        ( upperturbedRotationMatrix - downperturbedRotationMatrix ) / ( 2.0 * perturbation );

                Eigen::Matrix3d numericalRotationMatrixDerivativePartial =
                        ( upperturbedRotationMatrixDerivative - downperturbedRotationMatrixDerivative ) /
                        ( 2.0 * perturbation );

                Eigen::Matrix3d matrixDifference = rotationMatrixPartials.at( 0 ) - numericalRotationMatrixPartial;

//                std::cout << "analitical Rotation Matrix Partial"<< std::endl << rotationMatrixPartials.at(0) << std::endl;

//                std::cout << "numerical Rotation Matrix Partial"<< std::endl << numericalRotationMatrixPartial << std::endl;

//                std::cout << "rotation Matrix Partial difference"<< std::endl << matrixDifference << std::endl;


                // Compare analytical and numerical result.

                for( unsigned int i = 0; i < 3; i++ )
                {
                    for( unsigned int j = 0; j < 3; j++ )
                    {
                        BOOST_CHECK_SMALL( std::fabs( matrixDifference( i, j ) ), 1.0E-9 );
                    }
                }

                matrixDifference = rotationMatrixDerivativePartials.at( 0 ) - numericalRotationMatrixDerivativePartial;

//                std::cout << "analitical Rotation Matrix derivative Partial"<< std::endl <<
//                             rotationMatrixDerivativePartials.at(0) << std::endl;

//                std::cout << "numerical Rotation Matrix derivative Partial"<< std::endl <<
//                             numericalRotationMatrixDerivativePartial << std::endl;

//                std::cout << "rotation Matrix derivative Partial difference"<< std::endl <<
//                             matrixDifference << std::endl << std::endl<< std::endl;

                //Compare analytical and numerical result.

                for( unsigned int i = 0; i < 3; i++ )
                {
                    for( unsigned int j = 0; j < 3; j++ )
                    {
                        BOOST_CHECK_SMALL( std::fabs( matrixDifference( i, j ) ), 1.0E-13 );
                    }
                }

                marsRotationModelSettings->setPeriodTerms(
                            marsRotationModelSettings->getMeanMotionDirectNutationCorrections( ),
                            marsRotationModelSettings->getMeanMotionTimeDependentPhaseNutationCorrections( ),
                            rotationRateCorrections,
                            xPolarMotionCoefficients,
                            yPolarMotionCoefficients );
            }

            //Compute partial for Phi_sj numerically
            {
//                std::cout << std::endl<< std::endl<< "Phi_sj:" << std::endl;

                double perturbation = 9.0E-5;

                std::map< double, std::pair< double, double > > upperturbateRotationRateCorrections;

                std::map< double, std::pair< double, double > > downperturbateRotationRateCorrections;

                for( std::map< double, std::pair< double, double > >::iterator correctionIterator = rotationRateCorrections.begin( );
                     correctionIterator != rotationRateCorrections.end( ); correctionIterator++ )
                {
                    upperturbateRotationRateCorrections[ correctionIterator->first ] = std::make_pair(
                                correctionIterator->second.first, correctionIterator->second.second + perturbation );

                    downperturbateRotationRateCorrections[ correctionIterator->first ] = std::make_pair(
                                correctionIterator->second.first, correctionIterator->second.second - perturbation );

                    if ( correctionIterator->first > 1.0 )
                    {
                        upperturbateRotationRateCorrections[ correctionIterator->first ] = std::make_pair(
                                    correctionIterator->second.first, correctionIterator->second.second  );

                        downperturbateRotationRateCorrections[ correctionIterator->first ] = std::make_pair(
                                    correctionIterator->second.first, correctionIterator->second.second);
                    }
                }

                marsRotationModelSettings->setPeriodTerms(
                            marsRotationModelSettings->getMeanMotionDirectNutationCorrections( ),
                            marsRotationModelSettings->getMeanMotionTimeDependentPhaseNutationCorrections( ),
                            upperturbateRotationRateCorrections,
                            xPolarMotionCoefficients,
                            yPolarMotionCoefficients );

                std::shared_ptr< RotationalEphemeris > upperturbatedMarsRotationModel =
                        createRotationModel( marsRotationModelSettings, "Mars" );

                Eigen::Matrix3d upperturbedRotationMatrix =
                        upperturbatedMarsRotationModel->getRotationToBaseFrame( testTime ).toRotationMatrix( );


                Eigen::Matrix3d upperturbedRotationMatrixDerivative =
                        upperturbatedMarsRotationModel->getDerivativeOfRotationToBaseFrame( testTime );

                marsRotationModelSettings->setPeriodTerms(
                            marsRotationModelSettings->getMeanMotionDirectNutationCorrections( ),
                            marsRotationModelSettings->getMeanMotionTimeDependentPhaseNutationCorrections( ),
                            downperturbateRotationRateCorrections,
                            xPolarMotionCoefficients,
                            yPolarMotionCoefficients );

                std::shared_ptr< RotationalEphemeris > downperturbatedMarsRotationModel =
                        createRotationModel( marsRotationModelSettings, "Mars" );

                Eigen::Matrix3d downperturbedRotationMatrix =
                        downperturbatedMarsRotationModel->getRotationToBaseFrame( testTime ).toRotationMatrix( );

                Eigen::Matrix3d downperturbedRotationMatrixDerivative =
                        downperturbatedMarsRotationModel->getDerivativeOfRotationToBaseFrame( testTime );

                Eigen::Matrix3d numericalRotationMatrixPartial =
                        ( upperturbedRotationMatrix - downperturbedRotationMatrix ) / ( 2.0 * perturbation );

                Eigen::Matrix3d numericalRotationMatrixDerivativePartial =
                        ( upperturbedRotationMatrixDerivative - downperturbedRotationMatrixDerivative ) /
                        ( 2.0 * perturbation );

                Eigen::Matrix3d matrixDifference = rotationMatrixPartials.at( 1 ) - numericalRotationMatrixPartial;

//                std::cout << "analitical Rotation Matrix Partial"<< std::endl << rotationMatrixPartials.at( 1 ) << std::endl;

//                std::cout << "numerical Rotation Matrix Partial"<< std::endl << numericalRotationMatrixPartial << std::endl;

//                std::cout << "rotation Matrix Partial difference"<< std::endl << matrixDifference << std::endl;

                // Compare analytical and numerical result.

                for( unsigned int i = 0; i < 3; i++ )
                {
                    for( unsigned int j = 0; j < 3; j++ )
                    {
                        BOOST_CHECK_SMALL( std::fabs( matrixDifference( i, j ) ), 1.0E-9 );
                    }
                }

                matrixDifference = rotationMatrixDerivativePartials.at( 1 ) - numericalRotationMatrixDerivativePartial;

//                std::cout << "analitical Rotation Matrix derivative Partial"<< std::endl <<
//                             rotationMatrixDerivativePartials.at( 1 ) << std::endl;

//                std::cout << "numerical Rotation Matrix derivative Partial"<< std::endl <<
//                             numericalRotationMatrixDerivativePartial << std::endl;

//                std::cout << "rotation Matrix derivative Partial difference"<< std::endl <<
//                             matrixDifference << std::endl << std::endl<< std::endl;

                //Compare analytical and numerical result.

                for( unsigned int i = 0; i < 3; i++ )
                {
                    for( unsigned int j = 0; j < 3; j++ )
                    {
                        BOOST_CHECK_SMALL( std::fabs( matrixDifference( i, j ) ), 1.0E-13 );
                    }
                }

                marsRotationModelSettings->setPeriodTerms(
                            marsRotationModelSettings->getMeanMotionDirectNutationCorrections( ),
                            marsRotationModelSettings->getMeanMotionTimeDependentPhaseNutationCorrections( ),
                            rotationRateCorrections,
                            xPolarMotionCoefficients,
                            yPolarMotionCoefficients );
            }
        }
    }

    //Test for Rotation Matrix Partial Wrt Polar motion amplitude
    {
        // Create partial object.

        std::shared_ptr< RotationMatrixPartialWrtPolarMotionAmplitude > rotationMatrixPartialObject =
                std::make_shared< RotationMatrixPartialWrtPolarMotionAmplitude >(
                    std::dynamic_pointer_cast< PlanetaryRotationModel >( marsRotationModel ) );

        // Compute partial analytically

        double testTime = 0.9E8;

        std::vector< Eigen::Matrix3d > rotationMatrixPartials =
                rotationMatrixPartialObject->calculatePartialOfRotationMatrixToBaseFrameWrParameter( testTime );

        std::vector< Eigen::Matrix3d > rotationMatrixDerivativePartials =
                rotationMatrixPartialObject->calculatePartialOfRotationMatrixDerivativeToBaseFrameWrParameter( testTime );

        // Compute partial numerically
        {
            //Compute partial wrt X_cj numerically
            {
//                std::cout << std::endl<< std::endl << "X_cj:" << std::endl;

                double perturbation = 1.0E-5;

                std::map< double, std::pair< double, double > > upperturbatePolarMotionCoefficients;

                std::map< double, std::pair< double, double > > downperturbatePolarMotionCoefficients;

                for( std::map< double, std::pair< double, double > >::iterator correctionIterator = xPolarMotionCoefficients.begin( );
                     correctionIterator != xPolarMotionCoefficients.end( ); correctionIterator++ )
                {
                    upperturbatePolarMotionCoefficients[ correctionIterator->first ] = std::make_pair(
                                correctionIterator->second.first + perturbation , correctionIterator->second.second );

                    downperturbatePolarMotionCoefficients[ correctionIterator->first ] = std::make_pair(
                                correctionIterator->second.first - perturbation , correctionIterator->second.second );

                    if ( correctionIterator->first > 1.0 )
                    {
                        upperturbatePolarMotionCoefficients[ correctionIterator->first ] = std::make_pair(
                                    correctionIterator->second.first, correctionIterator->second.second  );

                        downperturbatePolarMotionCoefficients[ correctionIterator->first ] = std::make_pair(
                                    correctionIterator->second.first, correctionIterator->second.second);
                    }
                }

                marsRotationModelSettings->setPeriodTerms(
                            marsRotationModelSettings->getMeanMotionDirectNutationCorrections( ),
                            marsRotationModelSettings->getMeanMotionTimeDependentPhaseNutationCorrections( ),
                            rotationRateCorrections,
                            upperturbatePolarMotionCoefficients,
                            yPolarMotionCoefficients );

                std::shared_ptr< RotationalEphemeris > upperturbatedMarsRotationModel =
                        createRotationModel( marsRotationModelSettings, "Mars" );

                Eigen::Matrix3d upperturbedRotationMatrix = upperturbatedMarsRotationModel->getRotationToBaseFrame(
                            testTime).toRotationMatrix( );

                Eigen::Matrix3d upperturbedRotationMatrixDerivative = upperturbatedMarsRotationModel->getDerivativeOfRotationToBaseFrame(
                            testTime );

                marsRotationModelSettings->setPeriodTerms(
                            marsRotationModelSettings->getMeanMotionDirectNutationCorrections( ),
                            marsRotationModelSettings->getMeanMotionTimeDependentPhaseNutationCorrections( ),
                            rotationRateCorrections,
                            downperturbatePolarMotionCoefficients,
                            yPolarMotionCoefficients );

                std::shared_ptr< RotationalEphemeris > downperturbatedMarsRotationModel =
                        createRotationModel( marsRotationModelSettings, "Mars" );

                Eigen::Matrix3d downperturbedRotationMatrix =
                        downperturbatedMarsRotationModel->getRotationToBaseFrame( testTime ).toRotationMatrix( );

                Eigen::Matrix3d downperturbedRotationMatrixDerivative =
                        downperturbatedMarsRotationModel->getDerivativeOfRotationToBaseFrame( testTime );

                Eigen::Matrix3d numericalRotationMatrixPartial =
                        ( upperturbedRotationMatrix - downperturbedRotationMatrix ) / ( 2.0 * perturbation );

                Eigen::Matrix3d numericalRotationMatrixDerivativePartial =
                        ( upperturbedRotationMatrixDerivative - downperturbedRotationMatrixDerivative ) /
                        ( 2.0 * perturbation );

                Eigen::Matrix3d matrixDifference = rotationMatrixPartials.at( 0 ) - numericalRotationMatrixPartial;

//                std::cout << "analitical Rotation Matrix Partial"<< std::endl << rotationMatrixPartials.at( 0 ) << std::endl;

//                std::cout << "numerical Rotation Matrix Partial"<< std::endl << numericalRotationMatrixPartial << std::endl;

//                std::cout << "rotation Matrix Partial difference"<< std::endl << matrixDifference << std::endl;

                // Compare analytical and numerical result.

                for( unsigned int i = 0; i < 3; i++ )
                {
                    for( unsigned int j = 0; j < 3; j++ )
                    {
                        BOOST_CHECK_SMALL( std::fabs( matrixDifference( i, j ) ), 1.0E-10 );
                    }
                }

                matrixDifference = rotationMatrixDerivativePartials.at( 0 ) - numericalRotationMatrixDerivativePartial;

//                std::cout << "analitical Rotation Matrix derivative Partial"<< std::endl <<
//                             rotationMatrixDerivativePartials.at( 0 ) << std::endl;

//                std::cout << "numerical Rotation Matrix derivative Partial"<< std::endl <<
//                             numericalRotationMatrixDerivativePartial << std::endl;

//                std::cout << "rotation Matrix derivative Partial difference"<< std::endl << matrixDifference << std::endl;

                //Compare analytical and numerical result.

                for( unsigned int i = 0; i < 3; i++ )
                {
                    for( unsigned int j = 0; j < 3; j++ )
                    {
                        BOOST_CHECK_SMALL( std::fabs( matrixDifference( i, j ) ), 1.0E-15 );
                    }
                }

                marsRotationModelSettings->setPeriodTerms(
                            marsRotationModelSettings->getMeanMotionDirectNutationCorrections( ),
                            marsRotationModelSettings->getMeanMotionTimeDependentPhaseNutationCorrections( ),
                            rotationRateCorrections,
                            xPolarMotionCoefficients,
                            yPolarMotionCoefficients );
            }

            // Compute partial wrt X_sj numerically
            {
//                std::cout << std::endl<< std::endl << "X_sj:" << std::endl;

                double perturbation = 9.0E-5;

                std::map< double, std::pair< double, double > > upperturbatePolarMotionCoefficients;

                std::map< double, std::pair< double, double > > downperturbatePolarMotionCoefficients;

                for( std::map< double, std::pair< double, double > >::iterator correctionIterator = xPolarMotionCoefficients.begin( );
                     correctionIterator != xPolarMotionCoefficients.end( ); correctionIterator++ )
                {
                    upperturbatePolarMotionCoefficients[ correctionIterator->first ] = std::make_pair(
                                correctionIterator->second.first , correctionIterator->second.second + perturbation );

                    downperturbatePolarMotionCoefficients[ correctionIterator->first ] = std::make_pair(
                                correctionIterator->second.first , correctionIterator->second.second - perturbation );

                    if ( correctionIterator->first > 1.0 )
                    {
                        upperturbatePolarMotionCoefficients[ correctionIterator->first ] = std::make_pair(
                                    correctionIterator->second.first, correctionIterator->second.second  );

                        downperturbatePolarMotionCoefficients[ correctionIterator->first ] = std::make_pair(
                                    correctionIterator->second.first, correctionIterator->second.second);
                    }
                }

                marsRotationModelSettings->setPeriodTerms(
                            marsRotationModelSettings->getMeanMotionDirectNutationCorrections( ),
                            marsRotationModelSettings->getMeanMotionTimeDependentPhaseNutationCorrections( ),
                            rotationRateCorrections,
                            upperturbatePolarMotionCoefficients,
                            yPolarMotionCoefficients );

                std::shared_ptr< RotationalEphemeris > upperturbatedMarsRotationModel =
                        createRotationModel( marsRotationModelSettings, "Mars" );

                Eigen::Matrix3d upperturbedRotationMatrix = upperturbatedMarsRotationModel->getRotationToBaseFrame(
                            testTime).toRotationMatrix( );

                Eigen::Matrix3d upperturbedRotationMatrixDerivative = upperturbatedMarsRotationModel->getDerivativeOfRotationToBaseFrame(
                            testTime );

                marsRotationModelSettings->setPeriodTerms(
                            marsRotationModelSettings->getMeanMotionDirectNutationCorrections( ),
                            marsRotationModelSettings->getMeanMotionTimeDependentPhaseNutationCorrections( ),
                            rotationRateCorrections,
                            downperturbatePolarMotionCoefficients,
                            yPolarMotionCoefficients );

                std::shared_ptr< RotationalEphemeris > downperturbatedMarsRotationModel =
                        createRotationModel( marsRotationModelSettings, "Mars" );

                Eigen::Matrix3d downperturbedRotationMatrix =
                        downperturbatedMarsRotationModel->getRotationToBaseFrame( testTime ).toRotationMatrix( );

                Eigen::Matrix3d downperturbedRotationMatrixDerivative =
                        downperturbatedMarsRotationModel->getDerivativeOfRotationToBaseFrame( testTime );

                Eigen::Matrix3d numericalRotationMatrixPartial =
                        ( upperturbedRotationMatrix - downperturbedRotationMatrix ) / ( 2.0 * perturbation );

                Eigen::Matrix3d numericalRotationMatrixDerivativePartial =
                        ( upperturbedRotationMatrixDerivative - downperturbedRotationMatrixDerivative ) /
                        ( 2.0 * perturbation );

                Eigen::Matrix3d matrixDifference = rotationMatrixPartials.at( 1 ) - numericalRotationMatrixPartial;

//                std::cout << "analitical Rotation Matrix Partial"<< std::endl << rotationMatrixPartials.at( 1 ) << std::endl;

//                std::cout << "numerical Rotation Matrix Partial"<< std::endl << numericalRotationMatrixPartial << std::endl;

//                std::cout << "rotation Matrix Partial difference"<< std::endl << matrixDifference << std::endl;

                // Compare analytical and numerical result.

                for( unsigned int i = 0; i < 3; i++ )
                {
                    for( unsigned int j = 0; j < 3; j++ )
                    {
                        BOOST_CHECK_SMALL( std::fabs( matrixDifference( i, j ) ), 1.0E-10 );
                    }
                }

                matrixDifference = rotationMatrixDerivativePartials.at( 1 ) - numericalRotationMatrixDerivativePartial;

//                std::cout << "analitical Rotation Matrix derivative Partial"<< std::endl <<
//                             rotationMatrixDerivativePartials.at( 1 ) << std::endl;

//                std::cout << "numerical Rotation Matrix derivative Partial"<< std::endl <<
//                             numericalRotationMatrixDerivativePartial << std::endl;

//                std::cout << "rotation Matrix derivative Partial difference"<< std::endl <<
//                             matrixDifference << std::endl;

                //Compare analytical and numerical result.

                for( unsigned int i = 0; i < 3; i++ )
                {
                    for( unsigned int j = 0; j < 3; j++ )
                    {
                        BOOST_CHECK_SMALL( std::fabs( matrixDifference( i, j ) ), 1.0E-14 );
                    }
                }

                marsRotationModelSettings->setPeriodTerms(
                            marsRotationModelSettings->getMeanMotionDirectNutationCorrections( ),
                            marsRotationModelSettings->getMeanMotionTimeDependentPhaseNutationCorrections( ),
                            rotationRateCorrections,
                            xPolarMotionCoefficients,
                            yPolarMotionCoefficients );
            }

            // Compute partial wrt Y_cj numerically
            {
//                std::cout << std::endl<< std::endl << "Y_cj:" << std::endl;

                double perturbation = 9.0E-7;

                std::map< double, std::pair< double, double > > upperturbatePolarMotionCoefficients;

                std::map< double, std::pair< double, double > > downperturbatePolarMotionCoefficients;

                for( std::map< double, std::pair< double, double > >::iterator correctionIterator = yPolarMotionCoefficients.begin( );
                     correctionIterator != yPolarMotionCoefficients.end( ); correctionIterator++ )
                {
                    upperturbatePolarMotionCoefficients[ correctionIterator->first ] = std::make_pair(
                                correctionIterator->second.first + perturbation , correctionIterator->second.second );

                    downperturbatePolarMotionCoefficients[ correctionIterator->first ] = std::make_pair(
                                correctionIterator->second.first - perturbation , correctionIterator->second.second );

                    if ( correctionIterator->first > 1.0 )
                    {
                        upperturbatePolarMotionCoefficients[ correctionIterator->first ] = std::make_pair(
                                    correctionIterator->second.first, correctionIterator->second.second  );

                        downperturbatePolarMotionCoefficients[ correctionIterator->first ] = std::make_pair(
                                    correctionIterator->second.first, correctionIterator->second.second);
                    }
                }

                marsRotationModelSettings->setPeriodTerms(
                            marsRotationModelSettings->getMeanMotionDirectNutationCorrections( ),
                            marsRotationModelSettings->getMeanMotionTimeDependentPhaseNutationCorrections( ),
                            rotationRateCorrections,
                            xPolarMotionCoefficients,
                            upperturbatePolarMotionCoefficients );

                std::shared_ptr< RotationalEphemeris > upperturbatedMarsRotationModel =
                        createRotationModel( marsRotationModelSettings, "Mars" );

                Eigen::Matrix3d upperturbedRotationMatrix =
                        upperturbatedMarsRotationModel->getRotationToBaseFrame( testTime ).toRotationMatrix( );

                Eigen::Matrix3d upperturbedRotationMatrixDerivative =
                        upperturbatedMarsRotationModel->getDerivativeOfRotationToBaseFrame( testTime );

                marsRotationModelSettings->setPeriodTerms(
                            marsRotationModelSettings->getMeanMotionDirectNutationCorrections( ),
                            marsRotationModelSettings->getMeanMotionTimeDependentPhaseNutationCorrections( ),
                            rotationRateCorrections,
                            xPolarMotionCoefficients,
                            downperturbatePolarMotionCoefficients );

                std::shared_ptr< RotationalEphemeris > downperturbatedMarsRotationModel =
                        createRotationModel( marsRotationModelSettings, "Mars" );

                Eigen::Matrix3d downperturbedRotationMatrix =
                        downperturbatedMarsRotationModel->getRotationToBaseFrame( testTime ).toRotationMatrix( );

                Eigen::Matrix3d downperturbedRotationMatrixDerivative =
                        downperturbatedMarsRotationModel->getDerivativeOfRotationToBaseFrame( testTime );

                Eigen::Matrix3d numericalRotationMatrixPartial =
                        ( upperturbedRotationMatrix - downperturbedRotationMatrix ) / ( 2.0 * perturbation );

                Eigen::Matrix3d numericalRotationMatrixDerivativePartial =
                        ( upperturbedRotationMatrixDerivative - downperturbedRotationMatrixDerivative ) /
                        ( 2.0 * perturbation );

                Eigen::Matrix3d matrixDifference = rotationMatrixPartials.at( 2 ) - numericalRotationMatrixPartial;

////                std::cout << "analitical Rotation Matrix Partial"<< std::endl << rotationMatrixPartials.at( 2 ) << std::endl;

////                std::cout << "numerical Rotation Matrix Partial"<< std::endl << numericalRotationMatrixPartial << std::endl;

////                std::cout << "rotation Matrix Partial difference"<< std::endl << matrixDifference << std::endl;

                // Compare analytical and numerical result.

                for( unsigned int i = 0; i < 3; i++ )
                {
                    for( unsigned int j = 0; j < 3; j++ )
                    {
                        BOOST_CHECK_SMALL( std::fabs( matrixDifference( i, j ) ), 1.0E-10 );
                    }
                }

                matrixDifference = rotationMatrixDerivativePartials.at( 2 ) - numericalRotationMatrixDerivativePartial;

////                std::cout << "analitical Rotation Matrix derivative Partial"<< std::endl <<
////                             rotationMatrixDerivativePartials.at( 2 ) << std::endl;

////                std::cout << "numerical Rotation Matrix derivative Partial"<< std::endl <<
////                             numericalRotationMatrixDerivativePartial << std::endl;

////                std::cout << "rotation Matrix derivative Partial difference"<< std::endl << matrixDifference << std::endl;

                //Compare analytical and numerical result.

                for( unsigned int i = 0; i < 3; i++ )
                {
                    for( unsigned int j = 0; j < 3; j++ )
                    {
                        BOOST_CHECK_SMALL( std::fabs( matrixDifference( i, j ) ), 1.0E-14 );
                    }
                }

                marsRotationModelSettings->setPeriodTerms(
                            marsRotationModelSettings->getMeanMotionDirectNutationCorrections( ),
                            marsRotationModelSettings->getMeanMotionTimeDependentPhaseNutationCorrections( ),
                            rotationRateCorrections,
                            xPolarMotionCoefficients,
                            yPolarMotionCoefficients );
            }

            // Compute partial wrt Y_sj numerically
            {
//                std::cout << std::endl<< std::endl << "Y_sj:" << std::endl;

                double perturbation = 5.0E-5;

                std::map< double, std::pair< double, double > > upperturbatePolarMotionCoefficients;

                std::map< double, std::pair< double, double > > downperturbatePolarMotionCoefficients;

                for( std::map< double, std::pair< double, double > >::iterator correctionIterator = yPolarMotionCoefficients.begin( );
                     correctionIterator != yPolarMotionCoefficients.end( ); correctionIterator++ )
                {
                    upperturbatePolarMotionCoefficients[ correctionIterator->first ] = std::make_pair(
                                correctionIterator->second.first , correctionIterator->second.second + perturbation );

                    downperturbatePolarMotionCoefficients[ correctionIterator->first ] = std::make_pair(
                                correctionIterator->second.first , correctionIterator->second.second - perturbation );

                    if ( correctionIterator->first > 1.0 )
                    {
                        upperturbatePolarMotionCoefficients[ correctionIterator->first ] = std::make_pair(
                                    correctionIterator->second.first, correctionIterator->second.second  );

                        downperturbatePolarMotionCoefficients[ correctionIterator->first ] = std::make_pair(
                                    correctionIterator->second.first, correctionIterator->second.second);
                    }
                }

                marsRotationModelSettings->setPeriodTerms(
                            marsRotationModelSettings->getMeanMotionDirectNutationCorrections( ),
                            marsRotationModelSettings->getMeanMotionTimeDependentPhaseNutationCorrections( ),
                            rotationRateCorrections,
                            xPolarMotionCoefficients,
                            upperturbatePolarMotionCoefficients );

                std::shared_ptr< RotationalEphemeris > upperturbatedMarsRotationModel =
                        createRotationModel( marsRotationModelSettings, "Mars" );

                Eigen::Matrix3d upperturbedRotationMatrix =
                        upperturbatedMarsRotationModel->getRotationToBaseFrame( testTime ).toRotationMatrix( );

                Eigen::Matrix3d upperturbedRotationMatrixDerivative =
                        upperturbatedMarsRotationModel->getDerivativeOfRotationToBaseFrame( testTime );

                marsRotationModelSettings->setPeriodTerms(
                            marsRotationModelSettings->getMeanMotionDirectNutationCorrections( ),
                            marsRotationModelSettings->getMeanMotionTimeDependentPhaseNutationCorrections( ),
                            rotationRateCorrections,
                            xPolarMotionCoefficients,
                            downperturbatePolarMotionCoefficients );

                std::shared_ptr< RotationalEphemeris > downperturbatedMarsRotationModel =
                        createRotationModel( marsRotationModelSettings, "Mars" );

                Eigen::Matrix3d downperturbedRotationMatrix =
                        downperturbatedMarsRotationModel->getRotationToBaseFrame( testTime ).toRotationMatrix( );

                Eigen::Matrix3d downperturbedRotationMatrixDerivative =
                        downperturbatedMarsRotationModel->getDerivativeOfRotationToBaseFrame( testTime );

                Eigen::Matrix3d numericalRotationMatrixPartial =
                        ( upperturbedRotationMatrix - downperturbedRotationMatrix ) / ( 2.0 * perturbation );

                Eigen::Matrix3d numericalRotationMatrixDerivativePartial =
                        ( upperturbedRotationMatrixDerivative - downperturbedRotationMatrixDerivative ) /
                        ( 2.0 * perturbation );

                Eigen::Matrix3d matrixDifference = rotationMatrixPartials.at( 3 ) - numericalRotationMatrixPartial;

//                std::cout << "analitical Rotation Matrix Partial"<< std::endl << rotationMatrixPartials.at( 3 ) << std::endl;

//                std::cout << "numerical Rotation Matrix Partial"<< std::endl << numericalRotationMatrixPartial << std::endl;

//                std::cout << "rotation Matrix Partial difference"<< std::endl << matrixDifference << std::endl;

                // Compare analytical and numerical result.

                for( unsigned int i = 0; i < 3; i++ )
                {
                    for( unsigned int j = 0; j < 3; j++ )
                    {
                        BOOST_CHECK_SMALL( std::fabs( matrixDifference( i, j ) ), 1.0E-10 );
                    }
                }

                matrixDifference = rotationMatrixDerivativePartials.at( 3 ) - numericalRotationMatrixDerivativePartial;

//                std::cout << "analitical Rotation Matrix derivative Partial"<< std::endl <<
//                             rotationMatrixDerivativePartials.at( 3 ) << std::endl;

//                std::cout << "numerical Rotation Matrix derivative Partial"<< std::endl <<
//                             numericalRotationMatrixDerivativePartial << std::endl;

//                std::cout << "rotation Matrix derivative Partial difference"<< std::endl << matrixDifference << std::endl;

                //Compare analytical and numerical result.

                for( unsigned int i = 0; i < 3; i++ )
                {
                    for( unsigned int j = 0; j < 3; j++ )
                    {
                        BOOST_CHECK_SMALL( std::fabs( matrixDifference( i, j ) ), 1.0E-14 );
                    }
                }

                marsRotationModelSettings->setPeriodTerms(
                            marsRotationModelSettings->getMeanMotionDirectNutationCorrections( ),
                            marsRotationModelSettings->getMeanMotionTimeDependentPhaseNutationCorrections( ),
                            rotationRateCorrections,
                            xPolarMotionCoefficients,
                            yPolarMotionCoefficients );
            }
        }
    }

    //Test for Rotation Matrix Partial wrt Core Factor
    {
//        std::cout << std::endl << std::endl << "Core Factor" << std::endl;

        // Create partial object.

        std::shared_ptr< RotationMatrixPartialWrtCoreFactor > rotationMatrixPartialObject =
                std::make_shared< RotationMatrixPartialWrtCoreFactor >(
                    std::dynamic_pointer_cast< PlanetaryRotationModel >( marsRotationModel ) );

        // Compute partial analytically

        double testTime = 9.0E7;

        std::vector< Eigen::Matrix3d > rotationMatrixPartials =
                rotationMatrixPartialObject->calculatePartialOfRotationMatrixToBaseFrameWrParameter( testTime );

        std::vector< Eigen::Matrix3d > rotationMatrixDerivativePartials =
                rotationMatrixPartialObject->calculatePartialOfRotationMatrixDerivativeToBaseFrameWrParameter( testTime );

        // Compute partial numerically
        {
            double perturbation = 9.0E-2;

            double upperturbateCoreFactor = coreFactor + perturbation ;

            double downperturbateCoreFactor = coreFactor - perturbation ;

            marsRotationModelSettings->setCoreFactorAndFreeCoreNutation(
                        upperturbateCoreFactor, freeCoreNutationRate );

            std::shared_ptr< RotationalEphemeris > upperturbatedMarsRotationModel =
                    createRotationModel( marsRotationModelSettings, "Mars" );

            Eigen::Matrix3d upperturbedRotationMatrix =
                    upperturbatedMarsRotationModel->getRotationToBaseFrame( testTime ).toRotationMatrix( );

            Eigen::Matrix3d upperturbedRotationMatrixDerivative =
                    upperturbatedMarsRotationModel->getDerivativeOfRotationToBaseFrame(  testTime );

            marsRotationModelSettings->setCoreFactorAndFreeCoreNutation(
                        downperturbateCoreFactor, freeCoreNutationRate );

            std::shared_ptr< RotationalEphemeris > downperturbatedMarsRotationModel =
                    createRotationModel( marsRotationModelSettings, "Mars" );

            Eigen::Matrix3d downperturbedRotationMatrix =
                    downperturbatedMarsRotationModel->getRotationToBaseFrame( testTime ).toRotationMatrix( );

            Eigen::Matrix3d downperturbedRotationMatrixDerivative =
                    downperturbatedMarsRotationModel->getDerivativeOfRotationToBaseFrame( testTime );

            Eigen::Matrix3d numericalRotationMatrixPartial =
                    ( upperturbedRotationMatrix - downperturbedRotationMatrix ) / ( 2.0 * perturbation );

            Eigen::Matrix3d numericalRotationMatrixDerivativePartial =
                    ( upperturbedRotationMatrixDerivative - downperturbedRotationMatrixDerivative ) /
                    ( 2.0 * perturbation );

            Eigen::Matrix3d matrixDifference = rotationMatrixPartials.at( 0 ) - numericalRotationMatrixPartial;

//            std::cout << "analitical Rotation Matrix Partial"<< std::endl << rotationMatrixPartials.at( 0 ) << std::endl;

//            std::cout << "numerical Rotation Matrix Partial"<< std::endl << numericalRotationMatrixPartial << std::endl;

//            std::cout << "rotation Matrix Partial difference"<< std::endl << matrixDifference << std::endl;

            // Compare analytical and numerical result.

            for( unsigned int i = 0; i < 3; i++ )
            {
                for( unsigned int j = 0; j < 3; j++ )
                {
                    BOOST_CHECK_SMALL( std::fabs( matrixDifference( i, j ) ), 1.0E-12 );
                }
            }

            matrixDifference = rotationMatrixDerivativePartials.at( 0 ) - numericalRotationMatrixDerivativePartial;

//            std::cout << "analitical Rotation Matrix derivative Partial"<< std::endl <<
//                         rotationMatrixDerivativePartials.at( 0 ) << std::endl;

//            std::cout << "numerical Rotation Matrix derivative Partial"<< std::endl <<
//                         numericalRotationMatrixDerivativePartial << std::endl;

//            std::cout << "rotation Matrix derivative Partial difference"<< std::endl << matrixDifference << std::endl;

            //Compare analytical and numerical result.

            for( unsigned int i = 0; i < 3; i++ )
            {
                for( unsigned int j = 0; j < 3; j++ )
                {
                    BOOST_CHECK_SMALL( std::fabs( matrixDifference( i, j ) ), 1.0E-16 );
                }
            }
        }
    }

    //Test for Rotation Matrix Partial wrt Free Core Nutation
    {
//        std::cout << std::endl << std::endl << "Free Core Nutation" << std::endl;

        // Create partial object.

        std::shared_ptr< RotationMatrixPartialWrtFreeCoreNutationRate > rotationMatrixPartialObject =
                std::make_shared< RotationMatrixPartialWrtFreeCoreNutationRate >(
                    std::dynamic_pointer_cast< PlanetaryRotationModel >( marsRotationModel ) );

        // Compute partial analytically

        double testTime = 9.0E7;

        std::vector< Eigen::Matrix3d > rotationMatrixPartials =
                rotationMatrixPartialObject->calculatePartialOfRotationMatrixToBaseFrameWrParameter( testTime );

        std::vector< Eigen::Matrix3d > rotationMatrixDerivativePartials =
                rotationMatrixPartialObject->calculatePartialOfRotationMatrixDerivativeToBaseFrameWrParameter( testTime );

        // Compute partial numerically
        {
            double perturbation = 5.0E-11;

            double upperturbateFreeCoreNutationRate = freeCoreNutationRate + perturbation ;


            double downperturbateFreeCoreNutationRate = freeCoreNutationRate - perturbation ;

            marsRotationModelSettings->setCoreFactorAndFreeCoreNutation( coreFactor , upperturbateFreeCoreNutationRate );

            std::shared_ptr< RotationalEphemeris > upperturbatedMarsRotationModel =
                    createRotationModel( marsRotationModelSettings, "Mars" );

            Eigen::Matrix3d upperturbedRotationMatrix =
                    upperturbatedMarsRotationModel->getRotationToBaseFrame( testTime ).toRotationMatrix( );

            Eigen::Matrix3d upperturbedRotationMatrixDerivative =
                    upperturbatedMarsRotationModel->getDerivativeOfRotationToBaseFrame( testTime );

            marsRotationModelSettings->setCoreFactorAndFreeCoreNutation( coreFactor, downperturbateFreeCoreNutationRate );

            std::shared_ptr< RotationalEphemeris > downperturbatedMarsRotationModel =
                    createRotationModel( marsRotationModelSettings, "Mars" );

            Eigen::Matrix3d downperturbedRotationMatrix =
                    downperturbatedMarsRotationModel->getRotationToBaseFrame( testTime ).toRotationMatrix( );

            Eigen::Matrix3d downperturbedRotationMatrixDerivative =
                    downperturbatedMarsRotationModel->getDerivativeOfRotationToBaseFrame( testTime );

            Eigen::Matrix3d numericalRotationMatrixPartial =
                    ( upperturbedRotationMatrix - downperturbedRotationMatrix ) / ( 2.0 * perturbation );
            Eigen::Matrix3d numericalRotationMatrixDerivativePartial =
                    ( upperturbedRotationMatrixDerivative - downperturbedRotationMatrixDerivative ) /
                    ( 2.0 * perturbation );

            Eigen::Matrix3d matrixDifference = rotationMatrixPartials.at( 0 ) - numericalRotationMatrixPartial;

//            std::cout << "analitical Rotation Matrix Partial"<< std::endl << rotationMatrixPartials.at( 0 ) << std::endl;

//            std::cout << "numerical Rotation Matrix Partial"<< std::endl << numericalRotationMatrixPartial << std::endl;

//            std::cout << "rotation Matrix Partial difference"<< std::endl << matrixDifference << std::endl;

            // Compare analytical and numerical result.

            for( unsigned int i = 0; i < 3; i++ )
            {
                for( unsigned int j = 0; j < 3; j++ )
                {
                    BOOST_CHECK_SMALL( std::fabs( matrixDifference( i, j ) ), 1.0E-3 );
                }
            }

            matrixDifference = rotationMatrixDerivativePartials.at( 0 ) - numericalRotationMatrixDerivativePartial;

//            std::cout << "analitical Rotation Matrix derivative Partial"<< std::endl <<
//                         rotationMatrixDerivativePartials.at( 0 ) << std::endl;

//            std::cout << "numerical Rotation Matrix derivative Partial"<< std::endl <<
//                         numericalRotationMatrixDerivativePartial << std::endl;

//            std::cout << "rotation Matrix derivative Partial difference"<< std::endl << matrixDifference << std::endl;

            //Compare analytical and numerical result.

            for( unsigned int i = 0; i < 3; i++ )
            {
                for( unsigned int j = 0; j < 3; j++ )
                {
                    BOOST_CHECK_SMALL( std::fabs( matrixDifference( i, j ) ), 1.0E-7 );
                }
            }
        }
    }

}
BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat





