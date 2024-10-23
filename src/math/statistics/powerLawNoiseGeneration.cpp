#include "tudat/math/interpolators/linearInterpolator.h"
#include "tudat/math/statistics/fastFourierTransform.h"
#include "tudat/math/statistics/powerLawNoiseGeneration.h"


#include <valarray>

namespace tudat
{

namespace statistics
{

double evaluateAmplitudeSpectralDensity( const std::function< double( const double ) > powerSpectralDensityFunction, const double frequency )
{
    return std::sqrt( powerSpectralDensityFunction( frequency ) );
}

double evaluatePolynomialSum( const std::map< double, double >& amplitudesAtUnitFrequencyPerPower, const double frequency )
{
    double result = 0;
    for( std::map< double, double >::const_iterator powerIterator = amplitudesAtUnitFrequencyPerPower.begin( );
         powerIterator != amplitudesAtUnitFrequencyPerPower.end( ); powerIterator++ )
    {
        result += powerIterator->second * std::pow( frequency, powerIterator->first );
    }
    return result;
}

double evaluatePolynomialSquareRootSum( const std::map< double, double >& amplitudesAtUnitFrequencyPerPower, const double frequency )
{
    double result = 0.0;
    for( std::map< double, double >::const_iterator powerIterator = amplitudesAtUnitFrequencyPerPower.begin( );
         powerIterator != amplitudesAtUnitFrequencyPerPower.end( ); powerIterator++ )
    {
        result += powerIterator->second * std::pow( frequency, powerIterator->first );
    }
    return std::sqrt( result );
}

// In the sum, must now take into account frequency nodes
double evaluatePolynomialSquareRootSumBetweenFrequencyNodes( const std::map< double, double >& amplitudesAtUnitFrequencyPerPower, const std::vector<double> frequencyNodes, const double frequency )
{
    const std::vector<double>::const_iterator it = std::find_if(frequencyNodes.begin(), frequencyNodes.end(), [frequency] (double element) { return (frequency < element); });
    const int index = std::distance(frequencyNodes.begin(), it);

    double result = 0.0;
    int counter = 0;
    for( std::map< double, double >::const_iterator powerIterator = amplitudesAtUnitFrequencyPerPower.begin( );
         powerIterator != amplitudesAtUnitFrequencyPerPower.end( ); powerIterator++ )
    {
        //result += powerIterator->second * std::pow( frequency, powerIterator->first );
        if(counter == index-1 ){
            result = powerIterator->second * std::pow( frequency, powerIterator->first );
        }
        counter++;
    }
    return std::sqrt( result );
}


std::function< double( const double ) > createSpectralDensityFunctionFromTabulatedData(
        const std::map< double, double > spectralDensityMap )
{
    typedef interpolators::OneDimensionalInterpolator< double, double > LocalInterpolator;

    return std::bind( static_cast< double( LocalInterpolator::* )( const double ) >
                        ( &LocalInterpolator::interpolate ), std::make_shared< interpolators::LinearInterpolator< double, double > >(
                            spectralDensityMap ), std::placeholders::_1 );
}

std::pair< std::vector< double >, double > generateAmplitudeSpectrumNoise(
        const double maximumFrequency, const int numberOfFrequencySamples,
        const std::function< double( const double ) > ampitudeSpectalDensityFunction,
        const double seed )
{
    if( numberOfFrequencySamples <= 1 )
    {
        std::cerr<<"Error when making power law noise, number of points must be at least 2, is now "<<numberOfFrequencySamples<<std::endl;
    }

    double frequencyStepSize = maximumFrequency / static_cast< double >( numberOfFrequencySamples - 1 );

    // Generate set of discrete frequencies at which stochastic values of spectrum are to be generated.
    std::vector< double > fourierFrequencies;
    fourierFrequencies.reserve( numberOfFrequencySamples );
    for( int i = 0; i < numberOfFrequencySamples; i++ )
    {
        fourierFrequencies[ i ] = static_cast< double >( i ) * frequencyStepSize;
    }

    // Create random number generator (unit gaussian)
    std::vector< double > gaussianNoiseSettings;
    gaussianNoiseSettings.push_back( 0.0 );
    gaussianNoiseSettings.push_back( 1.0 );
    std::shared_ptr< RandomVariableGenerator< double > > gaussianRandomNumberGenerator =
            createBoostContinuousRandomVariableGenerator( normal_boost_distribution, gaussianNoiseSettings, seed );

    // Generate frequency domain noise
    std::vector< std::complex< double > > frequencyDomainNoise;
    frequencyDomainNoise.resize( numberOfFrequencySamples );
    double currentAmplitude;
    double amplitudeScaling = 2.0 * (
                static_cast< double >( numberOfFrequencySamples ) - 1.0 ) * std::sqrt( frequencyStepSize / 2.0 );

    //std::cout<<maximumFrequency * (
    //               static_cast< double >( numberOfFrequencySamples ) - 1.0 )<<" "<<1.0/frequencyStepSize<<" "<<1.0 / maximumFrequency<< std::endl;

    for( int j = 1; j < numberOfFrequencySamples; j++ )
    {
        currentAmplitude = std::sqrt( 0.5 ) * amplitudeScaling * ampitudeSpectalDensityFunction( fourierFrequencies[ j ] );
        frequencyDomainNoise[ j ] = ( std::complex< double >( currentAmplitude * gaussianRandomNumberGenerator->getRandomVariableValue( ),
                                                              currentAmplitude * gaussianRandomNumberGenerator->getRandomVariableValue( ) ) );
    }

    frequencyDomainNoise[ 0 ] = std::complex< double >( 0.0, 0.0 );
    frequencyDomainNoise[ numberOfFrequencySamples - 1 ] = std::complex< double >( frequencyDomainNoise[ numberOfFrequencySamples - 1 ].real( ), 0.0 );

    // Transform frequency domain noise to time domain.
    std::vector< double > timeDomainNoiseVector = fftw_interface::performInverseFftToRealData( frequencyDomainNoise );

    double timeStep = 1.0 / ( 2.0 * maximumFrequency );

    return std::make_pair( timeDomainNoiseVector, timeStep );
}

std::pair< std::vector< double >, double > generateAmplitudeSpectrumNoise(
        const double maximumFrequency, const double frequencyStep, const std::function< double( const double ) > ampitudeSpectalDensityFunction,
        const double seed )
{
    // Calculate number of ppints in frequency domain.
    int numberOfPoints = std::ceil( maximumFrequency / frequencyStep ) + 1;

    return generateAmplitudeSpectrumNoise( maximumFrequency, numberOfPoints, ampitudeSpectalDensityFunction, seed );
}

std::pair< std::vector< double >, double > generatePowerSpectrumNoise(
        const double maximumFrequency, const int numberOfFrequencySamples, const std::function< double( const double ) > powerSpectalDensityFunction,
        const double seed )
{
    return generateAmplitudeSpectrumNoise(
                maximumFrequency, numberOfFrequencySamples, std::bind( &evaluateAmplitudeSpectralDensity, powerSpectalDensityFunction, std::placeholders::_1 ), seed );
}

std::pair< std::vector< double >, double > generatePowerSpectrumNoise(
        const double maximumFrequency, const double frequencyStep, const std::function< double( const double ) > powerSpectalDensityFunction,
        const double seed )
{
    // Calculate number of ppints in frequency domain.
    int numberOfPoints = std::ceil( maximumFrequency / frequencyStep ) + 1;


    return generatePowerSpectrumNoise( maximumFrequency, numberOfPoints, powerSpectalDensityFunction, seed );
}

//! Function to generate power law noise by discretizing in frequency domain.
std::pair< std::vector< double >, double > generatePowerLawNoise(
        const double maximumFrequency, const int numberOfFrequencySamples, const std::map< double, double >& powersAndUnitAmplitudes,
        const double seed )
{
    return generateAmplitudeSpectrumNoise( maximumFrequency, numberOfFrequencySamples,
                                           std::bind( &evaluatePolynomialSquareRootSum, powersAndUnitAmplitudes, std::placeholders::_1 ), seed );
}


std::pair< std::vector< double >, double > generatePowerLawNoise(
        const double maximumFrequency, const int numberOfFrequencySamples, const std::map< double, double >& powersAndUnitAmplitudes,
        const std::vector<double> frequencyNodes,
        const double seed )
{
    return generateAmplitudeSpectrumNoise( maximumFrequency, numberOfFrequencySamples,
                                           std::bind( &evaluatePolynomialSquareRootSumBetweenFrequencyNodes, powersAndUnitAmplitudes, frequencyNodes, std::placeholders::_1 ), seed );
}

//! Function to generate power law noise by discretizing in frequency domain.
std::pair< std::vector< double >, double > generatePowerLawNoise(
        const double maximumFrequency, const double frequencyStep, const std::map< double, double >& powersAndUnitAmplitudes,
        const double seed )
{
    // Calculate number of ppints in frequency domain.
    int numberOfPoints = std::ceil( maximumFrequency / frequencyStep ) + 1;

    return generatePowerLawNoise( maximumFrequency, numberOfPoints, powersAndUnitAmplitudes, seed );
}

std::pair< std::vector< double >, double > generatePowerLawNoise(
        const double maximumFrequency, const int numberOfFrequencySamples, const double power, const double amplitudeAtUnityFrequency,
        const double seed )
{
    std::map< double, double > powerMap;
    powerMap[ power ] = amplitudeAtUnityFrequency;
    return generatePowerLawNoise( maximumFrequency, numberOfFrequencySamples, amplitudeAtUnityFrequency, seed );
}

//! Function to generate power law noise by discretizing in frequency domain.
std::pair< std::vector< double >, double > generatePowerLawNoise(
        const double maximumFrequency, const double frequencyStep, const double power, const double amplitudeAtUnityFrequency,
        const double seed )
{
    std::map< double, double > powerMap;
    powerMap[ power ] = amplitudeAtUnityFrequency;
    return generatePowerLawNoise( maximumFrequency, frequencyStep, powerMap, seed );
}


}

}

