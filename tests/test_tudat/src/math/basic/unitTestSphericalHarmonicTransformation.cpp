/*    Copyright (c) 2010-2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *
 *    References
 *
 */

#define BOOST_TEST_MAIN

#include <limits>
#include <string>

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>

#include "Tudat/InputOutput/basicInputOutput.h"

#include "Tudat/Mathematics/BasicMathematics/sphericalHarmonicTransformations.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createGravityField.h"


namespace tudat
{
namespace unit_tests
{

using namespace tudat::basic_mathematics;

BOOST_AUTO_TEST_SUITE( test_spherical_harmonic_transformations )

BOOST_AUTO_TEST_CASE( testAnalyticalTransformations )
{
    Eigen::MatrixXd nominalCosineCoefficients = Eigen::Matrix3d::Zero( );
    Eigen::MatrixXd nominalSineCoefficients = Eigen::Matrix3d::Zero( );

    Eigen::MatrixXd nominalNormalizedCosineCoefficients = Eigen::Matrix3d::Zero( );
    Eigen::MatrixXd nominalNormalizedSineCoefficients = Eigen::Matrix3d::Zero( );
    nominalCosineCoefficients( 0, 0 ) = 1.0;

    Eigen::MatrixXd transformedCosineCoefficients = Eigen::Matrix3d::Zero( );
    Eigen::MatrixXd transformedSineCoefficients = Eigen::Matrix3d::Zero( );

    Eigen::MatrixXd transformedNormalizedCosineCoefficients = Eigen::Matrix3d::Zero( );
    Eigen::MatrixXd transformedNormalizedSineCoefficients = Eigen::Matrix3d::Zero( );


    Eigen::MatrixXd transformedRenormalizedCosineCoefficients = Eigen::Matrix3d::Zero( );
    Eigen::MatrixXd transformedRenormalizedSineCoefficients = Eigen::Matrix3d::Zero( );

    SphericalHarmonicTransformationCache sphericalHarmonicTransformationCache( 2, 2 );

    double anglePsi = 0.0;
    double anglePhi = 0.0;
    double angleTheta = mathematical_constants::PI / 2.0;

    sphericalHarmonicTransformationCache.updateFrom313EulerAngles( anglePhi, angleTheta, anglePsi );

    double perturbationMagnitude = 0.1;

    {
        Eigen::MatrixXd normalizedCosineCoefficients = Eigen::Matrix3d::Zero( );
        Eigen::MatrixXd normalizedSineCoefficients = Eigen::Matrix3d::Zero( );
        nominalCosineCoefficients( 1, 0 ) = perturbationMagnitude;

        sphericalHarmonicTransformationCache.transformCoefficientsAtDegree(
                    nominalCosineCoefficients, nominalSineCoefficients,
                    transformedCosineCoefficients, transformedSineCoefficients, 0 );

        basic_mathematics::convertUnnormalizedToGeodesyNormalizedCoefficients(
                    nominalCosineCoefficients, nominalSineCoefficients,
                    nominalNormalizedCosineCoefficients, nominalNormalizedSineCoefficients );

        sphericalHarmonicTransformationCache.transformCoefficientsAtDegree(
                    nominalNormalizedCosineCoefficients, nominalNormalizedSineCoefficients,
                    transformedNormalizedCosineCoefficients, transformedNormalizedSineCoefficients, 1 );

        basic_mathematics::convertGeodesyNormalizedToUnnormalizedCoefficients(
                    transformedNormalizedCosineCoefficients, transformedNormalizedSineCoefficients,
                    transformedRenormalizedCosineCoefficients, transformedRenormalizedSineCoefficients );

        nominalCosineCoefficients( 1, 0 ) = 0.0;

        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                std::cout<<i<<" "<<j<<std::endl;
                if( i == 0 && j == 0 )
                {
                    BOOST_CHECK_SMALL( std::fabs( transformedCosineCoefficients( i, j ) - 1.0 ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedSineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedCosineCoefficients( i, j ) - 1.0 ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedSineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                }
                else if( i == 1 && j == 1 )
                {
                    BOOST_CHECK_SMALL( std::fabs( transformedCosineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedSineCoefficients( i, j ) + perturbationMagnitude ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedCosineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedSineCoefficients( i, j ) + perturbationMagnitude ), std::numeric_limits< double >::epsilon( ) );
                }
                else
                {
                    BOOST_CHECK_SMALL( std::fabs( transformedCosineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedSineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedCosineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedSineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                }

            }
        }

//        std::cout<<"Transformed: "<<std::endl<<transformedNormalizedCosineCoefficients<<std::endl<<std::endl<<
//                   transformedNormalizedSineCoefficients<<std::endl<<std::endl;
//        std::cout<<"Original: "<<std::endl<<nominalNormalizedCosineCoefficients<<std::endl<<std::endl<<
//                   nominalNormalizedSineCoefficients<<std::endl<<std::endl<<std::endl;
    }

    {
        nominalCosineCoefficients( 1, 1 ) = perturbationMagnitude;
        sphericalHarmonicTransformationCache.transformCoefficientsAtDegree(
                    nominalCosineCoefficients, nominalSineCoefficients,
                    transformedCosineCoefficients, transformedSineCoefficients, 0 );

        basic_mathematics::convertUnnormalizedToGeodesyNormalizedCoefficients(
                    nominalCosineCoefficients, nominalSineCoefficients,
                    nominalNormalizedCosineCoefficients, nominalNormalizedSineCoefficients );

        sphericalHarmonicTransformationCache.transformCoefficientsAtDegree(
                    nominalNormalizedCosineCoefficients, nominalNormalizedSineCoefficients,
                    transformedNormalizedCosineCoefficients, transformedNormalizedSineCoefficients, 1 );

        basic_mathematics::convertGeodesyNormalizedToUnnormalizedCoefficients(
                    transformedNormalizedCosineCoefficients, transformedNormalizedSineCoefficients,
                    transformedRenormalizedCosineCoefficients, transformedRenormalizedSineCoefficients );

        nominalCosineCoefficients( 1, 1 ) = 0.0;

        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                if( i == 0 && j == 0 )
                {
                    BOOST_CHECK_SMALL( std::fabs( transformedCosineCoefficients( i, j ) - 1.0 ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedSineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedCosineCoefficients( i, j ) - 1.0 ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedSineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                }
                else if( i == 1 && j == 1 )
                {
                    BOOST_CHECK_SMALL( std::fabs( transformedCosineCoefficients( i, j ) - perturbationMagnitude ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedSineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedCosineCoefficients( i, j )  - perturbationMagnitude ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedSineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                }
                else
                {
                    BOOST_CHECK_SMALL( std::fabs( transformedCosineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedSineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedCosineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedSineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                }

            }
        }
    }

    {
        nominalSineCoefficients( 1, 1 ) = perturbationMagnitude;
        sphericalHarmonicTransformationCache.transformCoefficientsAtDegree(
                    nominalCosineCoefficients, nominalSineCoefficients,
                    transformedCosineCoefficients, transformedSineCoefficients, 0 );

        basic_mathematics::convertUnnormalizedToGeodesyNormalizedCoefficients(
                    nominalCosineCoefficients, nominalSineCoefficients,
                    nominalNormalizedCosineCoefficients, nominalNormalizedSineCoefficients );

        sphericalHarmonicTransformationCache.transformCoefficientsAtDegree(
                    nominalNormalizedCosineCoefficients, nominalNormalizedSineCoefficients,
                    transformedNormalizedCosineCoefficients, transformedNormalizedSineCoefficients, 1 );

        basic_mathematics::convertGeodesyNormalizedToUnnormalizedCoefficients(
                    transformedNormalizedCosineCoefficients, transformedNormalizedSineCoefficients,
                    transformedRenormalizedCosineCoefficients, transformedRenormalizedSineCoefficients );

//        std::cout<<"Transformed: "<<std::endl<<transformedNormalizedCosineCoefficients<<std::endl<<std::endl<<
//                   transformedNormalizedSineCoefficients<<std::endl<<std::endl;
//        std::cout<<"Original: "<<std::endl<<nominalNormalizedCosineCoefficients<<std::endl<<std::endl<<
//                   nominalNormalizedSineCoefficients<<std::endl<<std::endl<<std::endl;

//        std::cout<<"Transformed: "<<std::endl<<transformedCosineCoefficients<<std::endl<<std::endl<<
//                   transformedSineCoefficients<<std::endl<<std::endl;
//        std::cout<<"Original: "<<std::endl<<nominalCosineCoefficients<<std::endl<<std::endl<<
//                   nominalSineCoefficients<<std::endl<<std::endl<<std::endl;

        nominalSineCoefficients( 1, 1 ) = 0.0;

        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                if( i == 0 && j == 0 )
                {
                    BOOST_CHECK_SMALL( std::fabs( transformedCosineCoefficients( i, j ) - 1.0 ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedSineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedCosineCoefficients( i, j ) - 1.0 ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedSineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                }
                else if( i == 1 && j == 0 )
                {
                    BOOST_CHECK_SMALL( std::fabs( transformedCosineCoefficients( i, j ) - perturbationMagnitude ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedSineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedCosineCoefficients( i, j ) - perturbationMagnitude ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedSineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                }
                else
                {
                    BOOST_CHECK_SMALL( std::fabs( transformedCosineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedSineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedCosineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedSineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                }

            }
        }
    }

    {
        nominalCosineCoefficients( 2, 0 ) = perturbationMagnitude;

        sphericalHarmonicTransformationCache.transformCoefficientsAtDegree(
                    nominalCosineCoefficients, nominalSineCoefficients,
                    transformedCosineCoefficients, transformedSineCoefficients, 0 );

        basic_mathematics::convertUnnormalizedToGeodesyNormalizedCoefficients(
                    nominalCosineCoefficients, nominalSineCoefficients,
                    nominalNormalizedCosineCoefficients, nominalNormalizedSineCoefficients );

        sphericalHarmonicTransformationCache.transformCoefficientsAtDegree(
                    nominalNormalizedCosineCoefficients, nominalNormalizedSineCoefficients,
                    transformedNormalizedCosineCoefficients, transformedNormalizedSineCoefficients, 1 );

        basic_mathematics::convertGeodesyNormalizedToUnnormalizedCoefficients(
                    transformedNormalizedCosineCoefficients, transformedNormalizedSineCoefficients,
                    transformedRenormalizedCosineCoefficients, transformedRenormalizedSineCoefficients );


        nominalSineCoefficients( 2, 0 ) = 0.0;

        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                if( i == 0 && j == 0 )
                {
                    BOOST_CHECK_SMALL( std::fabs( transformedCosineCoefficients( i, j ) - 1.0 ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedSineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedCosineCoefficients( i, j ) - 1.0 ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedSineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                }
                else if( i == 2 && j == 0 )
                {
                    BOOST_CHECK_SMALL( std::fabs( transformedCosineCoefficients( i, j ) + perturbationMagnitude / 2.0 ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedSineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedCosineCoefficients( i, j )  + perturbationMagnitude / 2.0 ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedSineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                }

                else if( i == 2 && j == 2 )
                {
                    BOOST_CHECK_SMALL( std::fabs( transformedCosineCoefficients( i, j ) + perturbationMagnitude / 4.0 ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedSineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedCosineCoefficients( i, j )  + perturbationMagnitude / 4.0 ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedSineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                }
                else
                {
                    BOOST_CHECK_SMALL( std::fabs( transformedCosineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedSineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedCosineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( transformedRenormalizedSineCoefficients( i, j ) ), std::numeric_limits< double >::epsilon( ) );
                }

            }
        }
    }
}

//BOOST_AUTO_TEST_CASE( testReverseSphericalHarmonicTransformation )
//{
//    for( int testBool = 0; testBool <= 1; testBool++ )
//    {
//        std::pair< std::pair< double, double >, std::pair< Eigen::MatrixXd, Eigen::MatrixXd > > marsGravityFieldSettings =
//                tudat::readPdsGravityFieldFile(
//                    tudat::input_output::getDataFilesRootPath( ) + "GravityFields/jgmro_110c_sha.tab.txt", 12, 12 );
//        Eigen::MatrixXd nominalCosineCoefficients = marsGravityFieldSettings.second.first;
//        Eigen::MatrixXd nominalSineCoefficients = marsGravityFieldSettings.second.second;

//        Eigen::MatrixXd transformedCosineCoefficients = Eigen::MatrixXd::Zero( 13, 13 );
//        Eigen::MatrixXd transformedSineCoefficients = Eigen::MatrixXd::Zero( 13, 13 );

//        Eigen::MatrixXd nominalNormalizedCosineCoefficients = Eigen::MatrixXd::Zero( 13, 13 );
//        Eigen::MatrixXd nominalNormalizedSineCoefficients = Eigen::MatrixXd::Zero( 13, 13 );

//        Eigen::MatrixXd transformedNormalizedCosineCoefficients = Eigen::MatrixXd::Zero( 13, 13 );
//        Eigen::MatrixXd transformedNormalizedSineCoefficients = Eigen::MatrixXd::Zero( 13, 13 );

//        Eigen::MatrixXd transformedRenormalizedCosineCoefficients = Eigen::MatrixXd::Zero( 13, 13 );
//        Eigen::MatrixXd transformedRenormalizedSineCoefficients = Eigen::MatrixXd::Zero( 13, 13 );

//        SphericalHarmonicTransformationCache sphericalHarmonicTransformationCache( 12, 12 );

//        if( !testBool )
//        {
//            basic_mathematics::convertUnnormalizedToGeodesyNormalizedCoefficients(
//                        nominalCosineCoefficients, nominalSineCoefficients,
//                        nominalNormalizedCosineCoefficients, nominalNormalizedSineCoefficients );
//        }

//        double anglePsi = 0.4;
//        double anglePhi = 0.3;
//        double angleTheta = -0.2;

//        sphericalHarmonicTransformationCache.update( angleTheta, anglePhi, anglePsi );
//        sphericalHarmonicTransformationCache.transformCoefficientsAtDegree(
//                    nominalCosineCoefficients, nominalSineCoefficients,
//                    transformedCosineCoefficients, transformedSineCoefficients, testBool );

//        if( !testBool )
//        {
//            sphericalHarmonicTransformationCache.transformCoefficientsAtDegree(
//                        nominalNormalizedCosineCoefficients, nominalNormalizedSineCoefficients,
//                        transformedNormalizedCosineCoefficients, transformedNormalizedSineCoefficients, 1 );
//            basic_mathematics::convertGeodesyNormalizedToUnnormalizedCoefficients(
//                        transformedNormalizedCosineCoefficients, transformedNormalizedSineCoefficients,
//                        transformedRenormalizedCosineCoefficients, transformedRenormalizedSineCoefficients );
//        }

//        Eigen::MatrixXd retransformedCosineCoefficients = Eigen::MatrixXd::Zero( 13, 13 );
//        Eigen::MatrixXd retransformedSineCoefficients = Eigen::MatrixXd::Zero( 13, 13 );

//        sphericalHarmonicTransformationCache.update( -angleTheta, -anglePsi, -anglePhi );
//        sphericalHarmonicTransformationCache.transformCoefficientsAtDegree(
//                    transformedCosineCoefficients, transformedSineCoefficients,
//                    retransformedCosineCoefficients, retransformedSineCoefficients, testBool );

//        Eigen::MatrixXd cosineCoefficientError = retransformedCosineCoefficients - nominalCosineCoefficients;
//        Eigen::MatrixXd sineCoefficientError = retransformedSineCoefficients - nominalSineCoefficients;

//        for( unsigned int i = 0; i < 13; i++ )
//        {
//            for( unsigned int j = 0; j < 13; j++ )
//            {
//                std::cout<<i<<" "<<j<<" "<<cosineCoefficientError( i, j )<<" "<<sineCoefficientError( i, j )<<std::endl;
//                BOOST_CHECK_SMALL( std::fabs( cosineCoefficientError( i, j ) ), 1.0E-12 );
//                BOOST_CHECK_SMALL( std::fabs( sineCoefficientError( i, j ) ), 1.0E-12 );

//                if( !testBool )
//                {
//                    BOOST_CHECK_SMALL( std::fabs( transformedCosineCoefficients( i, j ) - transformedRenormalizedCosineCoefficients( i, j ) ), 1.0E-12 );
//                    BOOST_CHECK_SMALL( std::fabs( transformedSineCoefficients( i, j )  - transformedRenormalizedSineCoefficients( i, j ) ), 1.0E-16 );
//                }

//            }
//        }
//    }
//}

//BOOST_AUTO_TEST_CASE( testHFunctionPartials )
//{
//    std::pair< std::pair< double, double >, std::pair< Eigen::MatrixXd, Eigen::MatrixXd > > marsGravityFieldSettings =
//            tudat::readPdsGravityFieldFile(
//                tudat::input_output::getDataFilesRootPath( ) + "GravityFields/jgmro_110c_sha.tab.txt", 12, 12 );
//    //for( int i = 0; i < 13; i++ )
//    //marsGravityFieldSettings.second.first( i, 0 ) = 0.0;

//    int maximumDegree = 12;
//    int maximumOrder = 12;

//    SphericalHarmonicTransformationCache sphericalHarmonicTransformationCache( 12, 12 );
//    sphericalHarmonicTransformationCache.setUpdatePartials( );

//    double angleTheta = 0.3;
//    double anglePhi = -0.8;
//    double anglePsi = 1.1;

//    sphericalHarmonicTransformationCache.update( angleTheta, anglePhi, anglePsi );

//    Eigen::MatrixXd transformedCosineCoefficients = Eigen::MatrixXd( 13, 13 );
//    Eigen::MatrixXd transformedSineCoefficients = Eigen::MatrixXd( 13, 13 );

//    Eigen::MatrixXd upPerturbedCosineCoefficients = Eigen::MatrixXd( 13, 13 );
//    Eigen::MatrixXd upPerturbedSineCoefficients = Eigen::MatrixXd( 13, 13 );

//    Eigen::MatrixXd downPerturbedCosineCoefficients = Eigen::MatrixXd( 13, 13 );
//    Eigen::MatrixXd downPerturbedSineCoefficients = Eigen::MatrixXd( 13, 13 );

//    sphericalHarmonicTransformationCache.transformCoefficientsAtDegree(
//            marsGravityFieldSettings.second.first,
//            marsGravityFieldSettings.second.second,
//            transformedCosineCoefficients,
//            transformedSineCoefficients, 0 );

//    std::vector< Eigen::MatrixXd > currentCosineCoefficientPartials;
//    std::vector< Eigen::MatrixXd > currentSineCoefficientPartials;

//    sphericalHarmonicTransformationCache.getPartialDerivativesOfTransformedCoefficientsWrtEulerAngles(
//                marsGravityFieldSettings.second.first,
//                marsGravityFieldSettings.second.second,
//                currentCosineCoefficientPartials, currentSineCoefficientPartials, 0 );

//    std::map< int, std::map< int, std::map< int, double > > > hFunctionPartials;
//    std::map< int, std::map< int, std::map< int, std::complex< double > > > > eFunctionPartials;

//    sphericalHarmonicTransformationCache.getHFunctionPartials( hFunctionPartials );
//    sphericalHarmonicTransformationCache.getEFunctionPartials( eFunctionPartials );

//    double anglePerturbation = 1.0E-6;

//    sphericalHarmonicTransformationCache.update( angleTheta + anglePerturbation, anglePhi, anglePsi );
//    std::map< int, std::map< int, std::vector< double > > > upPerturbedHFunctions;
//    std::map< int, std::map< int, std::map< int, std::complex< double > > > > upPerturbedEFunctions;

//    sphericalHarmonicTransformationCache.getHFunctions( upPerturbedHFunctions );
//    sphericalHarmonicTransformationCache.getEFunctions( upPerturbedEFunctions );
//    sphericalHarmonicTransformationCache.transformCoefficientsAtDegree(
//            marsGravityFieldSettings.second.first,
//            marsGravityFieldSettings.second.second,
//            upPerturbedCosineCoefficients,
//            upPerturbedSineCoefficients, 0 );


//    sphericalHarmonicTransformationCache.update( angleTheta - anglePerturbation, anglePhi, anglePsi );
//    std::map< int, std::map< int, std::vector< double > > > downPerturbedHFunctions;
//    std::map< int, std::map< int, std::map< int, std::complex< double > > > > downPerturbedEFunctions;

//    sphericalHarmonicTransformationCache.getHFunctions( downPerturbedHFunctions );
//    sphericalHarmonicTransformationCache.getEFunctions( downPerturbedEFunctions );
//    sphericalHarmonicTransformationCache.transformCoefficientsAtDegree(
//            marsGravityFieldSettings.second.first,
//            marsGravityFieldSettings.second.second,
//            downPerturbedCosineCoefficients,
//            downPerturbedSineCoefficients, 0 );

//    for( int l = 0; l <= maximumDegree; l++ )
//    {
//        for( int m = -std::min( l, maximumOrder ); m <= 0; m++ )
//        {
//            for( int k = -l; k <= l; k++ )
//            {
//                BOOST_CHECK_CLOSE_FRACTION(
//                            hFunctionPartials[ k ][ l ][ m ],
//                            ( upPerturbedHFunctions[ k ][ l ][ m + maximumDegree ] - downPerturbedHFunctions[ k ][ l ][ m + maximumDegree ] ) / ( 2.0 * anglePerturbation ), 1.0E-8 );
//                BOOST_CHECK_CLOSE_FRACTION(
//                            eFunctionPartials[ k ][ l ][ m ].real( ),
//                            ( upPerturbedEFunctions[ k ][ l ][ m ].real( ) - downPerturbedEFunctions[ k ][ l ][ m ].real( ) ) / ( 2.0 * anglePerturbation ), 1.0E-8 );
//                BOOST_CHECK_CLOSE_FRACTION(
//                            eFunctionPartials[ k ][ l ][ m ].imag( ),
//                            ( upPerturbedEFunctions[ k ][ l ][ m ].imag( ) - downPerturbedEFunctions[ k ][ l ][ m ].imag( ) ) / ( 2.0 * anglePerturbation ), 1.0E-8 );
//            }
//        }
//    }

//    Eigen::Vector3d nominalEulerAngles = ( Eigen::Vector3d( )<<anglePsi, anglePhi, angleTheta ).finished( );
//    Eigen::Vector3d perturbedEulerAngles;

//    for( unsigned int i = 0; i < 3; i++ )
//    {
//        perturbedEulerAngles = nominalEulerAngles;
//        perturbedEulerAngles( i ) += anglePerturbation;
//        sphericalHarmonicTransformationCache.update( perturbedEulerAngles( 2 ), perturbedEulerAngles( 1 ), perturbedEulerAngles( 0 ) );

//        sphericalHarmonicTransformationCache.transformCoefficientsAtDegree(
//                marsGravityFieldSettings.second.first,
//                marsGravityFieldSettings.second.second,
//                upPerturbedCosineCoefficients,
//                upPerturbedSineCoefficients, 0 );

//        perturbedEulerAngles = nominalEulerAngles;
//        perturbedEulerAngles( i ) -= anglePerturbation;
//        sphericalHarmonicTransformationCache.update( perturbedEulerAngles( 2 ), perturbedEulerAngles( 1 ), perturbedEulerAngles( 0 ) );

//        sphericalHarmonicTransformationCache.update( perturbedEulerAngles( 2 ), perturbedEulerAngles( 1 ), perturbedEulerAngles( 0 ) );
//        sphericalHarmonicTransformationCache.transformCoefficientsAtDegree(
//                marsGravityFieldSettings.second.first,
//                marsGravityFieldSettings.second.second,
//                downPerturbedCosineCoefficients,
//                downPerturbedSineCoefficients, 0 );

//        for( unsigned int j = 0; j <= maximumDegree; j++ )
//        {
//            for( unsigned int k = 0; k <= maximumDegree; k++ )
//            {
//                BOOST_CHECK_CLOSE_FRACTION(
//                            ( upPerturbedCosineCoefficients( j, k ) - downPerturbedCosineCoefficients( j, k ) ) / ( 2.0 * anglePerturbation ),
//                            currentCosineCoefficientPartials[ i ]( j, k ), 1.0E-8 );
//                BOOST_CHECK_CLOSE_FRACTION(
//                            ( upPerturbedSineCoefficients( j, k ) - downPerturbedSineCoefficients( j, k ) ) / ( 2.0 * anglePerturbation ),
//                            currentSineCoefficientPartials[ i ]( j, k ), 1.0E-8 );
//            }
//        }
//    }
//}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat




