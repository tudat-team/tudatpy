#include <cmath>
#include <iostream>
#include <iomanip>

#include <boost/make_shared.hpp>

#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/astro/system_models/timingSystem.h"
#include "tudat/math/basic/linearAlgebra.h"
//#include "Mathematics/Statistics/powerLawNoiseGeneration.h"

namespace tudat
{

namespace system_models
{
//
////! Function to convert allan variance amplitudes (time domain) to phase noise amplitudes (frequency domain)
//std::map< int, double > convertAllanVarianceAmplitudesToPhaseNoiseAmplitudes(
//        const std::map< int, double > allanVarianceAmplitudes,
//        const double frequencyDomainCutoffFrequency,
//        const bool isInverseSquareTermFlickerPhaseNoise )
//{
//    std::map< int, double > phaseNoiseAmplitudes;
//
//    // Iterate over all allan variance amplitudes and calculate associated phase noise amplitude.
//    for( std::map< int, double >::const_iterator amplitudeIterator = allanVarianceAmplitudes.begin( );
//         amplitudeIterator != allanVarianceAmplitudes.end( ); amplitudeIterator++ )
//    {
//        // Check if requested amplitude is within implemented range.
//        if( amplitudeIterator->first < -2 || amplitudeIterator->first > 1 )
//        {
//            std::cerr<<"Error, Allan variance amplitude not in correct range, found "<<amplitudeIterator->first<<std::endl;
//        }
//        else
//        {
//            switch( amplitudeIterator->first )
//            {
//            case -2:
//                if( isInverseSquareTermFlickerPhaseNoise )
//                {
//                    phaseNoiseAmplitudes[ -1 ] = amplitudeIterator->second * 1.0 / (
//                                3.0 * ( mathematical_constants::EULER_MASCHERONI_CONSTANT +  2.0 *
//                                        std::log( 2.0 * mathematical_constants::PI * frequencyDomainCutoffFrequency ) ) - std::log( 2 ) );
//                }
//                else
//                {
//                    phaseNoiseAmplitudes[ 0 ] = amplitudeIterator->second / ( 3.0 * frequencyDomainCutoffFrequency );
//                }
//                break;
//            case -1:
//                phaseNoiseAmplitudes[ -2 ] = amplitudeIterator->second / ( 2.0 * mathematical_constants::PI * mathematical_constants::PI );
//                break;
//            case 0:
//                phaseNoiseAmplitudes[ -3 ] = amplitudeIterator->second / ( 8.0 * mathematical_constants::PI * mathematical_constants::PI * std::log( 2.0 ) );
//                break;
//            case 1:
//                phaseNoiseAmplitudes[ -4 ] = 3.0 * amplitudeIterator->second / ( 8.0 * std::pow( mathematical_constants::PI, 4.0 ) );
//                break;
//            default:
//                std::cerr<<"Error, did not recognize Allan variance power of "<<amplitudeIterator->first<<std::endl;
//            }
//        }
//    }
//    return phaseNoiseAmplitudes;
//}
//
////! Function to generate clock noise for a clock with given allan variance behaviour
//std::pair< std::vector< double >, double > generateClockNoise( const std::map< int, double > allanVarianceAmplitudes,
//                                                               const double startTime, const double endTime,
//                                                               const int numberOfTimeSteps, const bool isInverseSquareTermFlickerPhaseNoise,
//                                                               const double seed )
//{
//    using namespace statistics;
//
//    // Calculate time step between subsequent realizations of stochastic process
//    double timeStep = ( endTime - startTime ) / static_cast< double >( numberOfTimeSteps );
//
//    // Calculate maximum noise frequency that can be generated with goven time step
//    double maximumFrequency = 0.5 / timeStep;
//
//    // Calculate number of steps in frequency domain (for real time domain data)
//    int numberOfFrequencySteps = numberOfTimeSteps / 2;
//
//    if( numberOfTimeSteps % 2 == 1 )
//    {
//        numberOfFrequencySteps++;
//    }
//
//    // Convert time domain amplitudes to frequency domain amplitudes.
//    std::map< int, double > phaseNoiseAmplitudes = convertAllanVarianceAmplitudesToPhaseNoiseAmplitudes(
//                allanVarianceAmplitudes, maximumFrequency, isInverseSquareTermFlickerPhaseNoise );
//
//    std::map< double, double > doublePhaseNoiseAmplitudes;
//
//    for( std::map< int, double >::iterator noiseIterator = phaseNoiseAmplitudes.begin( ); noiseIterator != phaseNoiseAmplitudes.end( );
//         noiseIterator++ )
//    {
//        doublePhaseNoiseAmplitudes[ static_cast< double >( noiseIterator->first ) ] = noiseIterator->second;
//
//    }
//    // return clock noise with time step used.
//    return generatePowerLawNoise(
//        maximumFrequency, numberOfFrequencySteps, doublePhaseNoiseAmplitudes, seed );
//}
//
//std::map< double, double > generateClockNoiseMap( const std::map< int, double >& allanVarianceAmplitudes,
//                                                  const double startTime, const double endTime,
//                                                  const double timeStep, const bool isInverseSquareTermFlickerPhaseNoise, const double seed )
//{
//    int numberOfTimeSteps = std::ceil( ( endTime - startTime ) / timeStep );
//
//    std::cout<<"Generating clock noise with :"<<numberOfTimeSteps<<" steps"<<std::endl;
//
//    std::pair< std::vector< double >, double > clockNoise = generateClockNoise(
//                allanVarianceAmplitudes, startTime, endTime, numberOfTimeSteps, isInverseSquareTermFlickerPhaseNoise, seed );
//    std::map< double, double > clockNoiseMap;
//
//    for( unsigned int i = 0; i < clockNoise.first.size( ); i++ )
//    {
//        clockNoiseMap[ startTime + static_cast< double >( i * clockNoise.second ) ] = clockNoise.first[ i ];
//    }
//
//    return clockNoiseMap;
//}
//
//std::function< double( const double ) > getClockNoiseInterpolator(
//        const std::map< int, double > allanVarianceAmplitudes,
//        const double startTime, const double endTime,
//        const double timeStep, const bool isInverseSquareTermFlickerPhaseNoise, const double seed )
//{
//    std::map< double, double > clockNoiseMap = generateClockNoiseMap( allanVarianceAmplitudes, startTime, endTime, timeStep,
//                                                                      isInverseSquareTermFlickerPhaseNoise, seed );
//
//    typedef interpolators::OneDimensionalInterpolator< double, double > LocalInterpolator;
//
//    return std::bind(
//                static_cast< double( LocalInterpolator::* )( const double ) >
//                ( &LocalInterpolator::interpolate ), std::make_shared< interpolators::LinearInterpolatorDouble >( clockNoiseMap ), _1 );
//}

TimingSystem::TimingSystem( const std::vector< Time > arcTimes,
              const std::vector< double > allArcsPolynomialDriftCoefficients,
              const std::function< std::function< double( const double ) >( const double, const double, const double ) >
              clockNoiseGenerationFunction,  const double clockNoiseTimeStep )
{
    for( unsigned int i = 0; i < arcTimes.size( ) - 1; i++ )
    {
        polynomialDriftCoefficients_.push_back( allArcsPolynomialDriftCoefficients );
        clockErrorFunction_.push_back( clockNoiseGenerationFunction( arcTimes[ i ], arcTimes[ i + 1 ], clockNoiseTimeStep ) );
        synchronizationTimes_.push_back( arcTimes[ i ] );
    }

    currentTimeArcLookup_ = std::make_shared< interpolators::HuntingAlgorithmLookupScheme< Time > >( synchronizationTimes_ );
}

TimingSystem::TimingSystem( const std::vector< Time > arcTimes,
              const std::vector< std::vector< double > > polynomialDriftCoefficients,
              const std::function< std::function< double( const double ) >( const double, const double, const double ) >
              clockNoiseGenerationFunction,
              const double clockNoiseTimeStep )
{
    if( arcTimes.size( ) != ( polynomialDriftCoefficients.size( ) + 1 ) )
    {
        throw std::runtime_error(
                "Error when making timing system, inconsistent polynomial error vector size" +
                std::to_string( arcTimes.size( ) ) + ", " + std::to_string( polynomialDriftCoefficients.size( ) ) );
    }

    polynomialDriftCoefficients_ = polynomialDriftCoefficients;

    for( unsigned int i = 0; i < arcTimes.size( ) - 1; i++ )
    {
        clockErrorFunction_.push_back( clockNoiseGenerationFunction( arcTimes[ i ], arcTimes[ i + 1 ], clockNoiseTimeStep ) );
        synchronizationTimes_.push_back( arcTimes[ i ] );
    }

    currentTimeArcLookup_ = std::make_shared< interpolators::HuntingAlgorithmLookupScheme< Time > >( synchronizationTimes_ );
}

Time TimingSystem::generateTimeTag( const Time trueTime )
{
    return trueTime + getCompleteClockError( trueTime );
}

Time TimingSystem::getTimeIntoCurrentArc( const Time trueTime )
{
    return getTimeIntoCurrentArcAndArcIndex( trueTime ).first;
}

std::pair< Time, int > TimingSystem::getTimeIntoCurrentArcAndArcIndex( const Time trueTime )
{
    int currentArcIndex = currentTimeArcLookup_->findNearestLowerNeighbour( trueTime );


    return std::make_pair( trueTime - synchronizationTimes_[ currentArcIndex ], currentArcIndex );
}

double TimingSystem::getCompleteClockError( const Time trueTime )
{
    std::pair< Time, int > timeIntoArcAndIndex = getTimeIntoCurrentArcAndArcIndex( trueTime );

    return getDeterministicClockError( timeIntoArcAndIndex.first, timeIntoArcAndIndex.second ) +
            clockErrorFunction_[ timeIntoArcAndIndex.second ](
                static_cast< double >( trueTime ) );
}

double TimingSystem::getDeterministicClockError( const Time timeIntoArc, const int currentArcIndex )
{
    double deterministicError = 0.0;

    for( unsigned int i = 0; i < polynomialDriftCoefficients_[ currentArcIndex ].size( ); i++ )
    {
        deterministicError += polynomialDriftCoefficients_[ currentArcIndex ][ i ] * std::pow( static_cast< double >( timeIntoArc ), i );
    }

    return deterministicError;
}

double TimingSystem::getDeterministicClockError( const Time trueTime )
{
    std::pair< Time, int > timeIntoArcAndIndex = getTimeIntoCurrentArcAndArcIndex( trueTime );

    return getDeterministicClockError( timeIntoArcAndIndex.first, timeIntoArcAndIndex.second );

}

double TimingSystem::getPolynomialDriftCoefficients( const int power, const int currentArcIndex  )
{
    return polynomialDriftCoefficients_[ currentArcIndex ].at( power );
}

double TimingSystem::getMeanPolynomialDriftCoefficients( const int power )
{
    double meanVariation = 0.0;

    for( unsigned int i = 0; i < polynomialDriftCoefficients_.size( ); i++ )
    {
        meanVariation += polynomialDriftCoefficients_[ i ][ power ];
    }
    return meanVariation / static_cast< double >( polynomialDriftCoefficients_.size( ) );
}

void TimingSystem::setGlobalPolynomialClockCorrections( const std::map< int, double >& newCorrections )
{
    std::map< int, double >::const_iterator correctionIterator;
    for( unsigned int i = 0; i < polynomialDriftCoefficients_.size( ); i++ )
    {
        for( correctionIterator = newCorrections.begin( ); correctionIterator != newCorrections.end( ); correctionIterator++ )
        {
            polynomialDriftCoefficients_[ i ][ correctionIterator->first ] = correctionIterator->second;
        }
    }
}

void TimingSystem::setSingleArcPolynomialClockCorrections( const std::map< int, double >& newCorrections,
                                             const int arcIndex )
{
    for( std::map< int, double >::const_iterator correctionIterator = newCorrections.begin( );
         correctionIterator != newCorrections.end( ); correctionIterator++ )
    {
        polynomialDriftCoefficients_[ arcIndex ][ correctionIterator->first ] = correctionIterator->second;
    }
}

}

}
