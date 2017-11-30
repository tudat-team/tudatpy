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

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat




