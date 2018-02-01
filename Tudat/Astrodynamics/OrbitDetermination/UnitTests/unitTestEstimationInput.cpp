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

#include "Tudat/Basics/testMacros.h"

#include "Tudat/Astrodynamics/OrbitDetermination/UnitTests/orbitDeterminationTestCases.h"
#include "Tudat/Astrodynamics/OrbitDetermination/podProcessing.h"


namespace tudat
{
namespace unit_tests
{
BOOST_AUTO_TEST_SUITE( test_estimation_input_output )

//! This test checks whether the input/output of the estimation (weights, a priori covariance, unscaled covariance) are
//! correctly handed
BOOST_AUTO_TEST_CASE( test_EstimationInputAndOutput )
{
    int simulationType = 0;

    Eigen::VectorXd parameterPerturbation = getDefaultInitialParameterPerturbation( );

    // Define stringent a priori covariance
    Eigen::MatrixXd inverseAPrioriCovariance = 1.0E32 * Eigen::MatrixXd::Identity( 7, 7 );

    // Define moderate a priori covariance
    Eigen::MatrixXd moderateInverseAPriopriCovariance = Eigen::MatrixXd::Zero( 7, 7 );
    for( unsigned int i = 0; i < 7; i++ )
    {
        moderateInverseAPriopriCovariance( i, i ) = 1.0 / ( 1.0E-6 * parameterPerturbation( i ) * parameterPerturbation( i ) );
    }

    // Run estimation with strong a priori covariance
    std::pair< boost::shared_ptr< PodOutput< double > >, Eigen::VectorXd > estimationOutputWithAprioriCovariance =
            executePlanetaryParameterEstimation< double, double >(
                simulationType, parameterPerturbation, inverseAPrioriCovariance );

    // Run estimation with effectively zero covariance
    std::pair< boost::shared_ptr< PodOutput< double > >, Eigen::VectorXd > estimationOutputWithSmallAprioriCovariance =
            executePlanetaryParameterEstimation< double, double >(
                simulationType, parameterPerturbation, 1.0E-64 * inverseAPrioriCovariance );

    // Run estimation with moderate a priori covariance
    std::pair< boost::shared_ptr< PodOutput< double > >, Eigen::VectorXd > estimationOutputWithModerateAprioriCovariance =
            executePlanetaryParameterEstimation< double, double >(
                simulationType, parameterPerturbation,  moderateInverseAPriopriCovariance );

    // Run estimation without a priori covariance
    std::pair< boost::shared_ptr< PodOutput< double > >, Eigen::VectorXd > estimationOutputWithoutAprioriCovariance =
            executePlanetaryParameterEstimation< double, double >(
                simulationType, parameterPerturbation );

    // Run estimation without a priori covariance and increased weights
    double constantWeight = 100.0;
    std::pair< boost::shared_ptr< PodOutput< double > >, Eigen::VectorXd > estimationOutputWithoutAprioriCovarianceAndWeakWeight =
            executePlanetaryParameterEstimation< double, double >(
                simulationType, parameterPerturbation, Eigen::MatrixXd::Zero( 7, 7 ), constantWeight);

    // Retrieve estimation errors and a priori covariances
    Eigen::MatrixXd tightConstraintInverseCovariance  =
            estimationOutputWithAprioriCovariance.first->getUnnormalizedInverseCovarianceMatrix( );
    Eigen::MatrixXd weakConstraintInverseCovariance  =
            estimationOutputWithSmallAprioriCovariance.first->getUnnormalizedInverseCovarianceMatrix( );
    Eigen::MatrixXd moderateConstraintInverseCovariance  =
            estimationOutputWithModerateAprioriCovariance.first->getUnnormalizedInverseCovarianceMatrix( );
    Eigen::MatrixXd noConstraintInverseCovariance  =
            estimationOutputWithoutAprioriCovariance.first->getUnnormalizedInverseCovarianceMatrix( );
    Eigen::MatrixXd noConstraintInverseCovarianceWithWeakWeight  =
            estimationOutputWithoutAprioriCovarianceAndWeakWeight.first->getUnnormalizedInverseCovarianceMatrix( );

    Eigen::VectorXd tightConstraintError  =
            estimationOutputWithAprioriCovariance.second;
    Eigen::VectorXd weakConstraintError  =
            estimationOutputWithSmallAprioriCovariance.second;
    Eigen::VectorXd moderateConstraintError  =
            estimationOutputWithModerateAprioriCovariance.second;
    Eigen::VectorXd noConstraintError  =
            estimationOutputWithoutAprioriCovariance.second;
    Eigen::VectorXd noConstraintWeakWeightError  =
            estimationOutputWithoutAprioriCovarianceAndWeakWeight.second;

    // Check if (effectively) unconstrained solutions converge at expected level
    for( unsigned int i = 0; i < 3; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( weakConstraintError( i ) ), 1.0E-2 );
        BOOST_CHECK_SMALL( std::fabs( weakConstraintError( i + 3 ) ), 1.0E-7 );

        BOOST_CHECK_SMALL( std::fabs( noConstraintError( i ) ), 1.0E-2 );
        BOOST_CHECK_SMALL( std::fabs( noConstraintError( i + 3 ) ), 1.0E-7 );

        BOOST_CHECK_SMALL( std::fabs( noConstraintWeakWeightError( i ) ), 1.0E-2 );
        BOOST_CHECK_SMALL( std::fabs( noConstraintWeakWeightError( i + 3 ) ), 1.0E-7 );
    }

    BOOST_CHECK_SMALL( std::fabs( weakConstraintError( 6 ) ), 500.0 );
    BOOST_CHECK_SMALL( std::fabs( noConstraintError( 6 ) ), 500.0 );
    BOOST_CHECK_SMALL( std::fabs( noConstraintWeakWeightError( 6 ) ), 500.0 );

    std::cout << weakConstraintError( 6 ) << " " << weakConstraintError( 6 ) << " " << noConstraintWeakWeightError( 6 ) << std::endl;
    for( unsigned int i = 0; i < 7; i++ )
    {
        // Check if moderately constrained solution has intermediate accuracy
        BOOST_CHECK_EQUAL( std::fabs( moderateConstraintError( i ) ) > std::fabs( noConstraintError( i ) ), true );
        BOOST_CHECK_EQUAL( std::fabs( moderateConstraintError( i ) ) < std::fabs( tightConstraintError( i ) ), true );

        // Check if very tightly constrained solution has not differed from a priori error
        BOOST_CHECK_CLOSE_FRACTION( tightConstraintError( i ), parameterPerturbation( i ), 1.0E-8 );

        for( unsigned int j = 0; j < 7; j++ )
        {
            // Check if weights are correctly processed into covarince
            BOOST_CHECK_CLOSE_FRACTION( constantWeight * noConstraintInverseCovariance( i, j ),
                                        noConstraintInverseCovarianceWithWeakWeight( i, j ), 1.0E-8 );

            // Check if tight a priori constraints are processed correctly to a posteriori covariance
            if( i == j )
            {
                BOOST_CHECK_CLOSE_FRACTION(
                            tightConstraintInverseCovariance( i, j ), 1.0E32, 1.0E-10 );
            }
            else
            {
                BOOST_CHECK_SMALL( tightConstraintInverseCovariance( i, j ) / tightConstraintInverseCovariance( i, i ), 1.0E-10 );

            }
        }
    }
}

//! Test whether the covariance is correctly computed as a function of time
BOOST_AUTO_TEST_CASE( test_CovarianceAsFunctionOfTime )
{
    std::pair< boost::shared_ptr< PodOutput< double > >, boost::shared_ptr< PodInput< double, double > > > podData;

    // Simulate covariances directly by propagating to different final tomes
    std::map< int, Eigen::MatrixXd > manualCovarianes;
    for( unsigned int i = 1; i < 5; i++ )
    {
        executeEarthOrbiterParameterEstimation< double, double >(
                    podData, 1.0E7, i, 0, false );
        manualCovarianes[ i ] = podData.first->getUnnormalizedCovarianceMatrix( );
    }

    // Use final calculations to compute covariance as a function of time
    std::map< double, Eigen::MatrixXd > automaticCovariances = simulation_setup::calculateCovarianceMatrixAsFunctionOfTime(
                podData.second, podData.first, 86400.0 - 1.0 );

    // Check consistency
    int counter = 1;
    for( std::map< double, Eigen::MatrixXd >::const_iterator covarianceIterator = automaticCovariances.begin( );
         covarianceIterator != automaticCovariances.end( ); covarianceIterator++ )
    {
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( covarianceIterator->second, manualCovarianes.at( counter ), 1.0E-8 );
        counter++;
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}

}



