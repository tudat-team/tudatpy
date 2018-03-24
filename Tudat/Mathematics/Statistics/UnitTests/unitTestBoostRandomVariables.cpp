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

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/random.hpp>

#include "Tudat/Mathematics/Statistics/boostProbabilityDistributions.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_boost_distributions )

//! Test wrappers for continuous probability distributions from boost.
BOOST_AUTO_TEST_CASE( testContinuousBoostDistribution )
{

    // Set independent variables.
    std::vector< double > testVector;
    std::vector< double > positiveTestVector;
    for( unsigned int i = 0; i < 1001; i++ )
    {
        testVector.push_back( -5.0 + static_cast< double >( i ) * 0.01 );
        positiveTestVector.push_back( static_cast< double >( i ) * 0.01 );
    }

    // Set probability vector (used for quantile computation).
    std::vector< double > probabilityVector;
    for( unsigned int i = 1; i < 100; i++ )
    {
        probabilityVector.push_back( static_cast< double >( i ) * 0.01 );
    }

    // Test uniform_distribution
    {
        std::vector< double > parameters;
        parameters.push_back( -1.0 );
        parameters.push_back( 2.5 );

        boost::shared_ptr< statistics::InvertibleContinuousProbabilityDistribution< double > > randomVariable =
                statistics::createBoostRandomVariable(
                    statistics::uniform_boost_distribution, parameters );
        boost::math::uniform_distribution< > manualDistribution( parameters.at( 0 ), parameters.at( 1 ) );

        // Test pdf and cdf.
        for( unsigned int i = 0; i < testVector.size( ); i++ )
        {
            BOOST_CHECK_EQUAL(
                        randomVariable->evaluatePdf( testVector.at( i ) ),
                        boost::math::pdf< double >( manualDistribution, testVector.at( i ) ) );
            BOOST_CHECK_EQUAL(
                        randomVariable->evaluateCdf( testVector.at( i ) ),
                        boost::math::cdf< double >( manualDistribution, testVector.at( i ) ) );
        }

        for( unsigned int i = 0; i < probabilityVector.size( ); i++ )
        {
            BOOST_CHECK_EQUAL(
                        randomVariable->evaluateInverseCdf( probabilityVector.at( i ) ),
                        boost::math::quantile< double >( manualDistribution, probabilityVector.at( i ) ) );
        }
    }

    // Test normal_boost_distribution
    {
        std::vector< double > parameters;
        parameters.push_back( -1.0 );
        parameters.push_back( 2.5 );

        boost::shared_ptr< statistics::InvertibleContinuousProbabilityDistribution< double > > randomVariable =
                statistics::createBoostRandomVariable(
                    statistics::normal_boost_distribution, parameters );
        boost::math::normal_distribution< > manualDistribution( parameters.at( 0 ), parameters.at( 1 ) );

        // Test pdf and cdf.
        for( unsigned int i = 0; i < testVector.size( ); i++ )
        {
            BOOST_CHECK_EQUAL(
                        randomVariable->evaluatePdf( testVector.at( i ) ),
                        boost::math::pdf< double >( manualDistribution, testVector.at( i ) ) );
            BOOST_CHECK_EQUAL(
                        randomVariable->evaluateCdf( testVector.at( i ) ),
                        boost::math::cdf< double >( manualDistribution, testVector.at( i ) ) );
        }

        // Test inverse cdf (quantile).
        for( unsigned int i = 0; i < probabilityVector.size( ); i++ )
        {
            BOOST_CHECK_EQUAL(
                        randomVariable->evaluateInverseCdf( probabilityVector.at( i ) ),
                        boost::math::quantile< double >( manualDistribution, probabilityVector.at( i ) ) );
        }

    }

    // Test exponential_distribution
    {
        std::vector< double > parameters;
        parameters.push_back( 2.5 );

        boost::shared_ptr< statistics::InvertibleContinuousProbabilityDistribution< double > > randomVariable =
                statistics::createBoostRandomVariable(
                    statistics::exponential_boost_distribution, parameters );
        boost::math::exponential_distribution< > manualDistribution( parameters.at( 0 ) );

        // Test pdf and cdf.
        for( unsigned int i = 0; i < testVector.size( ); i++ )
        {
            BOOST_CHECK_EQUAL(
                        randomVariable->evaluatePdf( positiveTestVector.at( i ) ),
                        boost::math::pdf< double >( manualDistribution, positiveTestVector.at( i ) ) );
            BOOST_CHECK_EQUAL(
                        randomVariable->evaluateCdf( positiveTestVector.at( i ) ),
                        boost::math::cdf< double >( manualDistribution, positiveTestVector.at( i ) ) );
        }

        // Test inverse cdf (quantile).
        for( unsigned int i = 0; i < probabilityVector.size( ); i++ )
        {
            BOOST_CHECK_EQUAL(
                        randomVariable->evaluateInverseCdf( probabilityVector.at( i ) ),
                        boost::math::quantile< double >( manualDistribution, probabilityVector.at( i ) ) );
        }

    }

    // Test gamma_distribution
    {
        std::vector< double > parameters;
        parameters.push_back( 1.0 );
        parameters.push_back( 2.5 );

        boost::shared_ptr< statistics::InvertibleContinuousProbabilityDistribution< double > > randomVariable =
                statistics::createBoostRandomVariable(
                    statistics::gamma_boost_distribution, parameters );
        boost::math::gamma_distribution< > manualDistribution( parameters.at( 0 ), parameters.at( 1 ) );

        // Test pdf and cdf.
        for( unsigned int i = 0; i < testVector.size( ); i++ )
        {
            BOOST_CHECK_EQUAL(
                        randomVariable->evaluatePdf( positiveTestVector.at( i ) ),
                        boost::math::pdf< double >( manualDistribution, positiveTestVector.at( i ) ) );
            BOOST_CHECK_EQUAL(
                        randomVariable->evaluateCdf( positiveTestVector.at( i ) ),
                        boost::math::cdf< double >( manualDistribution, positiveTestVector.at( i ) ) );
        }

        // Test inverse cdf (quantile).
        for( unsigned int i = 0; i < probabilityVector.size( ); i++ )
        {
            BOOST_CHECK_EQUAL(
                        randomVariable->evaluateInverseCdf( probabilityVector.at( i ) ),
                        boost::math::quantile< double >( manualDistribution, probabilityVector.at( i ) ) );
        }
    }

    // Test lognormal_distribution
    {
        std::vector< double > parameters;
        parameters.push_back( -1.0 );
        parameters.push_back( 2.5 );

        boost::shared_ptr< statistics::InvertibleContinuousProbabilityDistribution< double > > randomVariable =
                statistics::createBoostRandomVariable(
                    statistics::lognormal_boost_distribution, parameters );
        boost::math::lognormal_distribution< > manualDistribution( parameters.at( 0 ), parameters.at( 1 ) );

        // Test pdf and cdf.
        for( unsigned int i = 0; i < testVector.size( ); i++ )
        {
            BOOST_CHECK_EQUAL(
                        randomVariable->evaluatePdf( positiveTestVector.at( i ) ),
                        boost::math::pdf< double >( manualDistribution, positiveTestVector.at( i ) ) );
            BOOST_CHECK_EQUAL(
                        randomVariable->evaluateCdf( positiveTestVector.at( i ) ),
                        boost::math::cdf< double >( manualDistribution, positiveTestVector.at( i ) ) );
        }

        // Test inverse cdf (quantile).
        for( unsigned int i = 0; i < probabilityVector.size( ); i++ )
        {
            BOOST_CHECK_EQUAL(
                        randomVariable->evaluateInverseCdf( probabilityVector.at( i ) ),
                        boost::math::quantile< double >( manualDistribution, probabilityVector.at( i ) ) );
        }
    }

    // Test beta_boost_distribution
    {
        std::vector< double > parameters;
        parameters.push_back( 1.0 );
        parameters.push_back( 2.5 );

        boost::shared_ptr< statistics::InvertibleContinuousProbabilityDistribution< double > > randomVariable =
                statistics::createBoostRandomVariable(
                    statistics::beta_boost_distribution, parameters );
        boost::math::beta_distribution< > manualDistribution( parameters.at( 0 ), parameters.at( 1 ) );

        // Test pdf and cdf.
        for( unsigned int i = 0; i < probabilityVector.size( ); i++ )
        {
            BOOST_CHECK_EQUAL(
                        randomVariable->evaluatePdf( probabilityVector.at( i ) ),
                        boost::math::pdf< double >( manualDistribution, probabilityVector.at( i ) ) );
            BOOST_CHECK_EQUAL(
                        randomVariable->evaluateCdf( probabilityVector.at( i ) ),
                        boost::math::cdf< double >( manualDistribution, probabilityVector.at( i ) ) );
        }

        // Test inverse cdf (quantile).
        for( unsigned int i = 0; i < probabilityVector.size( ); i++ )
        {
            BOOST_CHECK_EQUAL(
                        randomVariable->evaluateInverseCdf( probabilityVector.at( i ) ),
                        boost::math::quantile< double >( manualDistribution, probabilityVector.at( i ) ) );
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat

