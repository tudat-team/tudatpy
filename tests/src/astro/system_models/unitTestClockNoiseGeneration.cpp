
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <iostream>
#include <ctime>

#include <boost/format.hpp>
#include <boost/test/unit_test.hpp>

#include "tudat/astro/system_models/timingSystem.h"
#include "tudat/math/interpolators/linearInterpolator.h"
#include "tudat/math/statistics/allanVariance.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_clock_noise_generator )

using namespace tudat::statistics;
using namespace tudat::system_models;
//
//// NOTE: test highly similar to power law noise test, modify to be 'Unit' test.
//BOOST_AUTO_TEST_CASE( testInterpolatedPowerLawClockNoise )
//{
//    int counter = 0;
//    double timeIntervalRatio = 0.0;
//    double allanVarianceRatio = 0.0;
//    std::map< double, double >::iterator startIterator;
//    std::map< double, double >::iterator endIterator;
//
//    for( int i = -2; i <=1; i++ )
//    {
//
//        std::map< int, double > allanVariances;
//        allanVariances[ i ] = 2.0;
//
//        int numberOfRuns = 100;
//        double currentValue = 0.0;
//        double summedValue = 0.0;
//
//        double currentExponentValue = 0.0;
//        double summedExponentValue = 0.0;
//
//        std::cout<<i<<std::endl;
//        for( int k = 0; k < numberOfRuns; k++ )
//        {
//            std::pair< std::vector< double >, double > timeDomainNoiseOutput = generateClockNoise(
//                        allanVariances, 0.0, 1024.0, std::pow( 2, 18 ) );
//            std::map< double, double > allanVariance = statistics::calculateAllanVarianceOfTimeDataSet(
//                        timeDomainNoiseOutput.first, timeDomainNoiseOutput.second );
//
//            startIterator = allanVariance.begin( );
//            endIterator = allanVariance.end( );
//            for( int m = 0; m < 4; m++ )
//            {
//                endIterator--;
//            }
//
//            timeIntervalRatio = endIterator->first / startIterator->first;
//            allanVarianceRatio = endIterator->second / startIterator->second;
//
//            currentExponentValue = std::log( allanVarianceRatio ) / std::log( timeIntervalRatio );
//            summedExponentValue += currentExponentValue;
//
//            interpolators::LinearInterpolator< double, double > allanVarianceInterpolator =
//                    interpolators::LinearInterpolator< double, double >( allanVariance );
//            currentValue = allanVarianceInterpolator.interpolate( 1.0 );
//            summedValue += currentValue;
//
//            //output::writeDoubleMapToFile( allanVariance, "allanVariance" + boost::lexical_cast< std::string >( counter ) + ".dat" );
//            //counter++;
//        }
//
//        double expectedAmplitde = 2.0;
//
//        BOOST_CHECK_CLOSE_FRACTION( expectedAmplitde, summedValue / static_cast< double >( numberOfRuns ), 0.1 );
//        BOOST_CHECK_SMALL( std::fabs( summedExponentValue / static_cast< double >( numberOfRuns ) - static_cast< double >( i ) ), 0.1 );
//    }
//}

BOOST_AUTO_TEST_CASE( testCombinedPowerClockNoise )
{
    std::map< int, double > allanVariances;
    allanVariances[ -2 ] = 200.0;
    allanVariances[ -1 ] = 25.0;
    allanVariances[ 0 ] = 0.5;
    allanVariances[ 1 ] = 0.1;

    double currentValue = 0.0;
    double summedValue = 0.0;

    int numberOfRuns = 500;
    for( int i = 0; i < numberOfRuns; i++ )
    {
        std::pair< std::vector< double >, double > timeDomainNoiseOutput = generateClockNoise(
                    allanVariances, 0.0, 8192.0, std::pow( 2, 16 ), 0, static_cast< double >( i ) * time( 0 ) );
        std::map< double, double > allanVariance = statistics::calculateAllanVarianceOfTimeDataSet(
                    timeDomainNoiseOutput.first, timeDomainNoiseOutput.second );

        interpolators::LinearInterpolator< double, double > allanVarianceInterpolator =
                interpolators::LinearInterpolator< double, double >( allanVariance );
        currentValue = allanVarianceInterpolator.interpolate( 16.0 );
        summedValue += currentValue;
    }

    double expectedAmplitde = 0.0;
    for( std::map< int, double >::iterator amplitudeIterator = allanVariances.begin( ); amplitudeIterator != allanVariances.end( );
         amplitudeIterator++ )
    {
        expectedAmplitde += amplitudeIterator->second * std::pow( 16.0, amplitudeIterator->first );
    }

    BOOST_CHECK_CLOSE_FRACTION( expectedAmplitde, summedValue / static_cast< double >( numberOfRuns ), 1.0E-2 );
}

BOOST_AUTO_TEST_SUITE_END( )

}

}


