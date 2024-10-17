#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>

#include <boost/make_shared.hpp>

#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/astro/system_models/timingSystem.h"
#include "tudat/math/basic/linearAlgebra.h"
#include "tudat/math/interpolators/linearInterpolator.h"
#include "tudat/math/quadrature/gaussianQuadrature.h"
#if(TUDAT_BUILD_WITH_FFTW3 )
#include "tudat/math/statistics/powerLawNoiseGeneration.h"
#endif

namespace tudat
{

namespace system_models
{
    //! Function to convert an allan variance to polynomial form
/*!
 *  Function to convert an allan variance nodes to polynomial form.
 *  (polynomial with integer time powers -2 >= power <= 1 )
 *  \param allanVarianceNodes Allan deviation nodes, with key as integration time (tau) and value as the avar
 *  \return allanVarianceAmplitudes Allan variance amplitudes, with key as integration time (tau) and value
 *  as a pair containing: the slope of the following node & the mulitplying factor (amplitude).
 */
std::map< int, std::pair<double, double> > convertAllanVarianceNodesToAmplitudes(
    const std::map< int, double > allanVarianceNodes,
    const std::string& varianceType ) {

    std::map< int, std::pair<double, double> > allanVarianceAmplitudes;
    std::vector< double > Mus;

    for( std::map< int, double >::const_iterator iter = allanVarianceNodes.begin( );
         iter != std::prev(allanVarianceNodes.end( )); iter++ ){
        double mu = ( std::log(std::next(iter) ->second ) - std::log(iter -> second ) ) /
                    ( std::log(std::next(iter) ->first )-std::log(iter -> first) );
        if ( ( ( mu < -2 || mu > 2 ) & (varianceType.compare("Allan") || varianceType.compare("Overlapping Allan")) )
             || ( ( mu < -2 || mu > 4 ) & (varianceType.compare("Hadamard") || varianceType.compare("Overlapping Hadamard")) ) )
        { std::cerr<<"Error, Allan variance slope outside boundaries (-2, 1), found "<< mu <<std::endl;}
        else{Mus.push_back( mu );}
    }

    int counter = 0;
    double Bi = 0.;
    for( std::map< int, double >::const_iterator iter = std::next(allanVarianceNodes.begin( ));
         iter != allanVarianceNodes.end( ); iter++ ){
        if(counter == 0){
            Bi = iter->second * std::pow(iter->first, - Mus[counter]);
            //std::cout << counter << " " << iter->second << " " << iter->first << " " << Mus[counter] << std::endl;
        }
        else {
            //std::cout << Bi << std::endl;
            Bi = Bi * std::pow( std::prev(iter)->first, (Mus[counter-1] - Mus[counter]));
            //std::cout << counter << " " << std::prev(iter)->first << " " << Mus[counter -1] << " " << Mus[counter] << std::endl;
        }
        allanVarianceAmplitudes[ std::prev(iter)->first ] = std::make_pair(Mus[counter], Bi);
        //allanVarianceAmplitudes[ Mus[counter] ] = Bi;
        counter++;
        }
    return allanVarianceAmplitudes;
}


//! Function to convert allan variance amplitudes (time domain) to phase noise amplitudes (frequency domain)
std::map< int, double > convertAllanVarianceAmplitudesToPhaseNoiseAmplitudes(
        const std::map< int, double > allanVarianceAmplitudes,
        const double frequencyDomainCutoffFrequency,
        const bool isInverseSquareTermFlickerPhaseNoise )
{
    std::map< int, double > phaseNoiseAmplitudes;

    // Iterate over all allan variance amplitudes and calculate associated phase noise amplitude.
    for( std::map< int, double >::const_iterator amplitudeIterator = allanVarianceAmplitudes.begin( );
         amplitudeIterator != allanVarianceAmplitudes.end( ); amplitudeIterator++ )
    {
        // Check if requested amplitude is within implemented range.
        if( amplitudeIterator->first < -2 || amplitudeIterator->first > 1 )
        {
            std::cerr<<"Error, Allan variance amplitude not in correct range, found "<<amplitudeIterator->first<<std::endl;
        }
        else
        {
            switch( amplitudeIterator->first )
            {
            case -2:
                if( isInverseSquareTermFlickerPhaseNoise )
                {
                    phaseNoiseAmplitudes[ -1 ] = amplitudeIterator->second * 1.0 / (
                                3.0 * ( mathematical_constants::EULER_MASCHERONI_CONSTANT +  2.0 *
                                        std::log( 2.0 * mathematical_constants::PI * frequencyDomainCutoffFrequency ) ) - std::log( 2 ) );
                }
                else
                {
                    phaseNoiseAmplitudes[ 0 ] = amplitudeIterator->second / ( 3.0 * frequencyDomainCutoffFrequency );
                }
                break;
            case -1:
                phaseNoiseAmplitudes[ -2 ] = amplitudeIterator->second / ( 2.0 * mathematical_constants::PI * mathematical_constants::PI );
                break;
            case 0:
                phaseNoiseAmplitudes[ -3 ] = amplitudeIterator->second / ( 8.0 * mathematical_constants::PI * mathematical_constants::PI * std::log( 2.0 ) );
                break;
            case 1:
                phaseNoiseAmplitudes[ -4 ] = 3.0 * amplitudeIterator->second / ( 8.0 * std::pow( mathematical_constants::PI, 4.0 ) );
                break;
            default:
                std::cerr<<"Error, did not recognize Allan variance power of "<<amplitudeIterator->first<<std::endl;
            }
        }
    }
    return phaseNoiseAmplitudes;
}



double genericIntegrationFunction(const double x, const double mu, const double sinExponent, const double factor){
    return factor * std::pow( std::sin( x ), sinExponent) / (std::pow(x, 3 + mu ) );
}
//! Function to convert allan variance amplitudes (time domain) to phase noise amplitudes (frequency domain)
/*!
 *  Function to convert allan variance amplitudes (time domain) to phase noise amplitudes (frequency domain) using DeMarchi et al method.
 *  \param allanVarianceAmplitudes Allan variance amplitudes, with key as power of time and value the mulitplying factor (amplitude).
 *  \param frequencyDomainCutoffFrequency Cut off frequency of noise generation (needed for time powers < -1, since their integration over all
 *  frequencies leads to a divergent solution, i.e. infinite power)
 *  \param isInverseSquareTermFlickerPhaseNoise Boolean determining whether inverse square time term is flicker phase noise (freq.^(-1)) or
 *  white phase modulation (freq.^(0)). Note that for allan variance of flicker phase noise, the logarithmic term will be neglected when using this
 *  function.
 *  \return Pair containing first: Map with phase noise amplitudes, powers of frequency as keys and ampliutdes as values, second: frequency_nodes
 */
std::pair< std::map<double, double>, std::vector<double > > convertAllanVarianceAmplitudesToApproximatedPhaseNoiseSpectrum(
        const std::map< int, std::pair < double, double > > allanVarianceAmplitudes,
        const std::string& varianceType)
{
    using namespace numerical_quadrature;
    //std::map< double, double > PowerSpectralDensity_Sy;
    std::map< double, double > PowerSpectralDensity_Sx;
    std::vector<double> his;
    std::vector<double> alphas;
    std::vector<double> frequencyNodes;
    double sinExponent;
    double factor;

    if (varianceType.compare("Allan") || varianceType.compare("Overlapping Allan")){
        factor = 1;
        sinExponent = 4;}
    else if (varianceType.compare("Hadamard") || varianceType.compare("Overlapping Hadamard")){
        factor = 8;
        sinExponent = 6;}
    else{ throw std::runtime_error( "Unknown varianceType" ); }

    const unsigned int numberOfNodes = 64;
    const double lowerLimit = 0.0;
    const double upperLimit = 10.0; //std::numeric_limits< double >::infinity( );
    for( std::map< int, std::pair <double,double> >::const_iterator iter  = allanVarianceAmplitudes.begin( );
         iter != allanVarianceAmplitudes.end( ); iter++ ) {
    //for( unsigned int i = 0; i != allanVarianceAmplitudes.size(); i++){
        double mu = iter->second.first;
        alphas.push_back( - mu - 1);
        auto integrationFunction = std::bind(&genericIntegrationFunction, std::placeholders::_1, mu, sinExponent, factor);
        GaussianQuadrature<double, double> integrator(integrationFunction, lowerLimit,  upperLimit, numberOfNodes);
        double integral = integrator.getQuadrature();
        double B = iter->second.second;
        double hi = B / ( 2 * integral * std::pow( mathematical_constants::PI, mu) );
        his.push_back(hi);
        // std::cout << " mu " << mu << " bi " <<  B <<  " integ " << integral << " hi " << hi << std::endl;
    }
    for( unsigned int i = 0; i != his.size() - 1; i++ ) {
        double frequencyNode = std::pow( his[i] / his[i+1], 1 / ( alphas[i+1] - alphas[i]) );
        // std::cout << frequencyNode << std::endl;
        frequencyNodes.push_back(frequencyNode);
    }
    std::reverse(frequencyNodes.begin(), frequencyNodes.end());
    frequencyNodes.insert(frequencyNodes.begin(), 0.0);
    frequencyNodes.push_back(std::numeric_limits< double >::infinity( ));
    std::reverse(alphas.begin(), alphas.end());
    std::reverse(his.begin(), his.end());

    for( unsigned int i = 0; i != allanVarianceAmplitudes.size(); i++ ){
        PowerSpectralDensity_Sx[ alphas[i] - 2] = his[i] / (4 * mathematical_constants::PI * mathematical_constants::PI);
    }

    return std::make_pair(PowerSpectralDensity_Sx, frequencyNodes);

}

#if( TUDAT_BUILD_WITH_FFTW3 )
//! Function to generate clock noise for a clock with given allan variance behaviour
std::pair< std::vector< double >, double > generateClockNoise( const std::map< int, double > allanVarianceAmplitudes,
                                                               const double startTime, const double endTime,
                                                               const int numberOfTimeSteps, const bool isInverseSquareTermFlickerPhaseNoise,
                                                               const double seed )
{
    using namespace statistics;

    // Calculate time step between subsequent realizations of stochastic process
    double timeStep = ( endTime - startTime ) / static_cast< double >( numberOfTimeSteps );

    // Calculate maximum noise frequency that can be generated with goven time step
    double maximumFrequency = 0.5 / timeStep;

    // Calculate number of steps in frequency domain (for real time domain data)
    int numberOfFrequencySteps = numberOfTimeSteps / 2;

    if( numberOfTimeSteps % 2 == 1 )
    {
        numberOfFrequencySteps++;
    }

    // Convert time domain amplitudes to frequency domain amplitudes.
    std::map< int, double > phaseNoiseAmplitudes = convertAllanVarianceAmplitudesToPhaseNoiseAmplitudes(
                allanVarianceAmplitudes, maximumFrequency, isInverseSquareTermFlickerPhaseNoise );

    std::map< double, double > doublePhaseNoiseAmplitudes;

    for( std::map< int, double >::iterator noiseIterator = phaseNoiseAmplitudes.begin( ); noiseIterator != phaseNoiseAmplitudes.end( );
         noiseIterator++ )
    {
        doublePhaseNoiseAmplitudes[ static_cast< double >( noiseIterator->first ) ] = noiseIterator->second;

    }
    // return clock noise with time step used.
    return generatePowerLawNoise(
        maximumFrequency, numberOfFrequencySteps, doublePhaseNoiseAmplitudes, seed );
}


//! Function to generate clock noise for a clock with given allan variance behaviour
std::pair< std::vector< double >, double > generateColoredClockNoise( const std::map< int, double > allanVarianceNodes,
                                                               const std::string& varianceType,
                                                               const double startTime, const double endTime,
                                                               const int numberOfTimeSteps,
                                                               const double seed )
{
    using namespace statistics;

    // Calculate time step between subsequent realizations of stochastic process
    double timeStep = ( endTime - startTime ) / static_cast< double >( numberOfTimeSteps );

    // Calculate maximum noise frequency that can be generated with goven time step
    double maximumFrequency = 0.5 / timeStep;

    // Calculate number of steps in frequency domain (for real time domain data)
    int numberOfFrequencySteps = numberOfTimeSteps / 2;

    if( numberOfTimeSteps % 2 == 1 )
    {
        numberOfFrequencySteps++;
    }

    // convert allan variance nodes  to allan variance amplitudes

    std::map< int, std::pair<double, double > > allanVarianceAmplitudes = convertAllanVarianceNodesToAmplitudes(allanVarianceNodes, varianceType);

    // Convert time domain amplitudes to frequency domain amplitudes.
    std::pair< std::map<double, double>, std::vector<double> > dummy = convertAllanVarianceAmplitudesToApproximatedPhaseNoiseSpectrum( allanVarianceAmplitudes, varianceType );
    std::map< double, double > phaseNoiseAmplitudes = dummy.first;
    std::vector<double> frequencyNodes = dummy.second;

    // return clock noise with time step used.
    return generatePowerLawNoise(
            maximumFrequency, numberOfFrequencySteps, phaseNoiseAmplitudes, frequencyNodes, seed );
}


std::map< double, double > generateClockNoiseMap( const std::map< int, double >& allanVarianceAmplitudes,
                                                  const double startTime, const double endTime,
                                                  const double timeStep, const bool isInverseSquareTermFlickerPhaseNoise, const double seed )
{
    int numberOfTimeSteps = std::ceil( ( endTime - startTime ) / timeStep );

    std::cout<<"Generating clock noise with :"<<numberOfTimeSteps<<" steps"<<std::endl;

    std::pair< std::vector< double >, double > clockNoise = generateClockNoise(
                allanVarianceAmplitudes, startTime, endTime, numberOfTimeSteps, isInverseSquareTermFlickerPhaseNoise, seed );
    std::map< double, double > clockNoiseMap;

    for( unsigned int i = 0; i < clockNoise.first.size( ); i++ )
    {
        clockNoiseMap[ startTime + static_cast< double >( i * clockNoise.second ) ] = clockNoise.first[ i ];
    }

    return clockNoiseMap;
}

std::map< double, double > generateColoredClockNoiseMap( const std::map< int, double >& allanVarianceNodes,
                                                         const std::string& varianceType,
                                                  const double startTime, const double endTime,
                                                  const double timeStep, const double seed )
{
    int numberOfTimeSteps = std::ceil( ( endTime - startTime ) / timeStep );

    std::cout<<"Generating clock noise with :"<<numberOfTimeSteps<<" steps"<<std::endl;

    std::pair< std::vector< double >, double > clockNoise = generateColoredClockNoise(
            allanVarianceNodes, varianceType, startTime, endTime, numberOfTimeSteps,  seed );
    std::map< double, double > clockNoiseMap;

    for( unsigned int i = 0; i < clockNoise.first.size( ); i++ )
    {
        clockNoiseMap[ startTime + static_cast< double >( i * clockNoise.second ) ] = clockNoise.first[ i ];
    }

    return clockNoiseMap;
}
std::function< double( const double ) > getClockNoiseInterpolator(
        const std::map< int, double > allanVarianceAmplitudes,
        const double startTime, const double endTime,
        const double timeStep, const bool isInverseSquareTermFlickerPhaseNoise, const double seed )
{
    std::map< double, double > clockNoiseMap = generateClockNoiseMap( allanVarianceAmplitudes, startTime, endTime, timeStep,
                                                                      isInverseSquareTermFlickerPhaseNoise, seed );

    typedef interpolators::OneDimensionalInterpolator< double, double > LocalInterpolator;

    return std::bind(
                static_cast< double( LocalInterpolator::* )( const double ) >
                ( &LocalInterpolator::interpolate ), std::make_shared< interpolators::LinearInterpolatorDouble >( clockNoiseMap ),
                        std::placeholders::_1 );
}

std::function< double( const double ) > getColoredClockNoiseInterpolator(
        const std::map< int, double > allanVarianceNodes,
        const std::string& varianceType,
        const double startTime, const double endTime,
        const double timeStep, const double seed)
        {
            std::map< double, double > clockNoiseMap = generateColoredClockNoiseMap( allanVarianceNodes, varianceType, startTime, endTime, timeStep , seed );

            typedef interpolators::OneDimensionalInterpolator< double, double > LocalInterpolator;

            return std::bind(
                    static_cast< double( LocalInterpolator::* )( const double ) >
                    ( &LocalInterpolator::interpolate ), std::make_shared< interpolators::LinearInterpolatorDouble >( clockNoiseMap ),
                    std::placeholders::_1 );
        }
#endif


TimingSystem::TimingSystem( const std::vector< Time > arcTimes,
              const std::vector< double > allArcsPolynomialDriftCoefficients,
              const std::function< std::function< double( const double ) >( const double, const double, const double ) >
              clockNoiseGenerationFunction,  const double clockNoiseTimeStep )
{
    for( unsigned int i = 0; i < arcTimes.size( ) - 1; i++ )
    {
        polynomialDriftCoefficients_.push_back( allArcsPolynomialDriftCoefficients );
        synchronizationTimes_.push_back( arcTimes[ i ] );
        if( clockNoiseGenerationFunction != nullptr )
        {
            clockErrorFunction_.push_back( clockNoiseGenerationFunction( arcTimes[ i ], arcTimes[ i + 1 ], clockNoiseTimeStep ) );
        }
        else
        {
            clockErrorFunction_.push_back( [=](const double){return 0.0;} );
        }

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
        if( clockNoiseGenerationFunction != nullptr )
        {
            clockErrorFunction_.push_back( clockNoiseGenerationFunction( arcTimes[ i ], arcTimes[ i + 1 ], clockNoiseTimeStep ) );
        }
        else
        {
            clockErrorFunction_.push_back( [=](const double){return 0.0;} );
        }
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
