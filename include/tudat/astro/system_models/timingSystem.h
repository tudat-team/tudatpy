#ifndef TIMINGSYSTEM_H
#define TIMINGSYSTEM_H

#include <map>
#include <vector>
#include <iostream>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/basics/timeType.h"
#include "tudat/math/interpolators/lookupScheme.h"

namespace tudat
{

namespace system_models
{

//
////! Function to convert allan variance amplitudes (time domain) to phase noise amplitudes (frequency domain)
///*!
// *  Function to convert allan variance amplitudes (time domain) to phase noise amplitudes (frequency domain) using NIST Special Publication 1065.
// *  \param allanVarianceAmplitudes Allan variance amplitudes, with key as power of time and value the mulitplying factor (amplitude).
// *  \param frequencyDomainCutoffFrequency Cut off frequency of noise generation (needed for time powers < -1, since their integration over all
// *  frequencies leads to a divergent solution, i.e. infinite power)
// *  \param isInverseSquareTermFlickerPhaseNoise Boolean determining whether inverse square time term is flicker phase noise (freq.^(-1)) or
// *  white phase modulation (freq.^(0)). Note that for allan variance of flicker phase noise, the logarithmic term will be neglected when using this
// *  function.
// *  \param Map with phase noise amplitudes, powers of frequency as keys and ampliutdes as values.
// */
//std::map< int, double > convertAllanVarianceAmplitudesToPhaseNoiseAmplitudes(
//        const std::map< int, double > allanVarianceAmplitudes,
//        const double frequencyDomainCutoffFrequency,
//        const bool isInverseSquareTermFlickerPhaseNoise = 0 );
//
////! Function to generate clock noise for a clock with given allan variance behaviour
///*!
// *  Function to generate clock noise for a clock with given allan variance behaviour (polynomial with integer time powers -2 >= power <= 1 )
// *  The function uses the convertAllanVarianceAmplitudesToPhaseNoiseAmplitudes function to convert allan variance to phase noise amplitudes in
// *  frequency domain. Subsequently, a discretized phase noise spectrum is generated in frequency domain, which is subsequently Fourier transformed
// *  to time domain to yoield a time series of random errors, which are a realization of the requested allan variance behaviour.
// *  \param allanVarianceAmplitudes Allan variance amplitudes, with key as power of time and value the mulitplying factor (amplitude).
// *  \param startTime Start time of noise time series generation
// *  \param endTime End time of noise time series generation
// *  \param numberOfTimeSteps Number of discrete relaizations of continuous random process to generate, i.e. number of data points in clock
// *  noise vector.
// *  \param isInverseSquareTermFlickerPhaseNoise Boolean determining whether inverse square time term is flicker phase noise (freq.^(-1)) or
// *  white phase modulation (freq.^(0)). Note that for allan variance of flicker phase noise, the logarithmic term will be neglected when using this
// *  function.
// *  \return Pair containing first: vector of clock noise realizations, second: time step between subsequent realizations of clock noise.
// */
//std::pair< std::vector< double >, double > generateClockNoise(
//        const std::map< int, double > allanVarianceAmplitudes, const double startTime, const double endTime,
//        const int numberOfTimeSteps, const bool isInverseSquareTermFlickerPhaseNoise = 0, const double seed = time( 0 ) );
//
////! Function to generate clock noise for a clock with given allan variance behaviour
///*!
// *  Function to generate clock noise for a clock with given allan variance behaviour (polynomial with integer time powers -2 >= power <= 1 )
// *  The function uses the convertAllanVarianceAmplitudesToPhaseNoiseAmplitudes function to convert allan variance to phase noise amplitudes in
// *  frequency domain. Subsequently, a discretized phase noise spectrum is generated in frequency domain, which is subsequently Fourier transformed
// *  to time domain to yoield a time series of random errors, which are a realization of the requested allan variance behaviour.
// *  \param allanVarianceAmplitudes Allan variance amplitudes, with key as power of time and value the mulitplying factor (amplitude).
// *  \param startTime Start time of noise time series generation
// *  \param endTime End time of noise time series generation
// *  \param numberOfTimeSteps Number of discrete relaizations of continuous random process to generate, i.e. number of data points in clock
// *  noise vector.
// *  \param isInverseSquareTermFlickerPhaseNoise Boolean determining whether inverse square time term is flicker phase noise (freq.^(-1)) or
// *  white phase modulation (freq.^(0)). Note that for allan variance of flicker phase noise, the logarithmic term will be neglected when using this
// *  function.
// *  \return Map of clock noise realizations, with associated times as keys and noise as values.
// */
//std::map< double, double > generateClockNoiseMap( const std::map< int, double >& allanVarianceAmplitudes,
//                                                  const double startTime, const double endTime,
//                                                  const double timeStep, const bool isInverseSquareTermFlickerPhaseNoise = 0,
//                                                  const double seed = time( 0 ) );
//
////! Function to generate a function to retrieve clock noise for a clock with given allan variance behaviour
///*!
// *  Function to generate a function to retrieve clock noise for a clock with given allan variance behaviour
// *  (polynomial with integer time powers -2 >= power <= 1 )
// *  The function uses the convertAllanVarianceAmplitudesToPhaseNoiseAmplitudes function to convert allan variance to phase noise amplitudes in
// *  frequency domain. Subsequently, a discretized phase noise spectrum is generated in frequency domain, which is subsequently Fourier transformed
// *  to time domain to yoield a time series of random errors, which are a realization of the requested allan variance behaviour. Subsequently,
// *  an interpolator is created using the discrete noise realizations. Note that the power spectrum is no longer valid at time differences smaller
// *  than the time between two subsequent noise realizations.
// *  \param allanVarianceAmplitudes Allan variance amplitudes, with key as power of time and value the mulitplying factor (amplitude).
// *  \param startTime Start time of noise time series generation
// *  \param endTime End time of noise time series generation
// *  \param numberOfTimeSteps Number of discrete relaizations of continuous random process to generate, i.e. number of data points in clock
// *  noise vector.
// *  \param isInverseSquareTermFlickerPhaseNoise Boolean determining whether inverse square time term is flicker phase noise (freq.^(-1)) or
// *  white phase modulation (freq.^(0)). Note that for allan variance of flicker phase noise, the logarithmic term will be neglected when using this
// *  function.
// *  \return Function that returns clock noise as a function of time.
// */
//std::function< double( const double ) > getClockNoiseInterpolator(
//        const std::map< int, double > allanVarianceAmplitudes,
//        const double startTime, const double endTime,
//        const double timeStep, const bool isInverseSquareTermFlickerPhaseNoise = 0, const double seed = time( 0 ) );

//! Class to represent timing system hardware (such as USO).
class TimingSystem
{
public:
    TimingSystem( const std::vector< Time > arcTimes,
                  const std::vector< double > allArcsPolynomialDriftCoefficients = std::vector< double >( ),
                  const std::function< std::function< double( const double ) >( const double, const double, const double ) >
                  clockNoiseGenerationFunction = nullptr,
                  const double clockNoiseTimeStep = 1.0E-3 );

    TimingSystem( const std::vector< Time > arcTimes,
                  const std::vector< std::vector< double > > polynomialDriftCoefficients,
                  const std::function< std::function< double( const double ) >( const double, const double, const double ) >
                  clockNoiseGenerationFunction = nullptr,
                  const double clockNoiseTimeStep = 1.0E-3 );


    TimingSystem( const std::vector< std::vector< double > > polynomialDriftCoefficients,
                  const std::vector< std::function< double( const double ) > > stochasticClockNoiseFunctions,
                  const std::vector< Time > arcTimes ):
        polynomialDriftCoefficients_( polynomialDriftCoefficients ), clockErrorFunction_( stochasticClockNoiseFunctions )

    {
        for( unsigned int i = 0; i < arcTimes.size( ) - 1; i++ )
        {
            synchronizationTimes_.push_back( arcTimes[ i ] );
        }

        currentTimeArcLookup_ = std::make_shared< interpolators::HuntingAlgorithmLookupScheme< Time > >( synchronizationTimes_ );
    }

    ~TimingSystem( ){ }

    //! Function to generate a time tag of an event.
    /*!
     *  Function to generate a time tag of an event, takes the true time as input and uses the internal clock error models to simulate
     *  a time tag.
     */
    Time generateTimeTag( const Time trueTime );

    //! Computes the time into the current arc (trueTime - arc start time)
    Time getTimeIntoCurrentArc( const Time trueTime );

    //! Computes the time into the current arc (trueTime - arc start time) and the index of the current arc
    std::pair< Time, int > getTimeIntoCurrentArcAndArcIndex( const Time trueTime );

    //! Computes the total clock error at the current (true) time
    double getCompleteClockError( const Time trueTime );

    double getDeterministicClockError( const Time timeIntoArc, const int currentArcIndex );

    //! Computes the deterministic clock error at the current (true) time
    double getDeterministicClockError( const Time trueTime );

    double getPolynomialDriftCoefficients( const int power, const int currentArcIndex );

    double getMeanPolynomialDriftCoefficients( const int power );

    void setGlobalPolynomialClockCorrections( const std::map< int, double >& newCorrections );

    void setSingleArcPolynomialClockCorrections(
            const std::map< int, double >& newCorrections, const int arcIndex );

    int getNumberOfSynchronizationTimes( )
    {
        return synchronizationTimes_.size( );
    }

    std::vector< Time > getSynchronizationTimes( )
    {
        return synchronizationTimes_;
    }

    void resetClockErrorFunction( const std::vector< std::function< double( const double ) > >& clockErrorFunction )
    {
        if( clockErrorFunction.size( ) != clockErrorFunction_.size( ) )
        {
            std::cerr<<"Error when resetting clock error function, size is incompatible"<<std::endl;
        }
        else
        {
            clockErrorFunction_ = clockErrorFunction;
        }
    }

    void resetClockParameters( const std::vector< std::vector< double > > polynomialDriftCoefficients,
                               const std::vector< std::function< double( const double ) > > stochasticClockNoiseFunctions,
                               const std::vector< Time > arcTimes )
    {
        polynomialDriftCoefficients_ = polynomialDriftCoefficients;
        clockErrorFunction_ = stochasticClockNoiseFunctions;

        synchronizationTimes_.clear( );
        for( unsigned int i = 0; i < arcTimes.size( ) - 1; i++ )
        {
            synchronizationTimes_.push_back( arcTimes[ i ] );
        }

        currentTimeArcLookup_ = std::make_shared< interpolators::HuntingAlgorithmLookupScheme< Time > >( synchronizationTimes_ );
    }


protected:

    std::shared_ptr< interpolators::LookUpScheme< Time > > currentTimeArcLookup_;

    std::vector< Time > synchronizationTimes_;

    std::vector< std::vector< double > > polynomialDriftCoefficients_; // In seconds per second^n

    std::vector< std::function< double( const double ) > > clockErrorFunction_;//Input in seconds since synchronization epoch.

};

}

}

#endif // TIMINGSYSTEM_H
