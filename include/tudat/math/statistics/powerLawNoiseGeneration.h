#ifndef POWERLAWNOISEGENERATION_H
#define POWERLAWNOISEGENERATION_H

#include <vector>
#include <map>

#include "tudat/math/statistics/randomVariableGenerator.h"

namespace tudat
{

namespace statistics
{

double evaluateAmplitudeSpectralDensity( const std::function< double( const double ) > powerSpectralDensityFunction, const double frequency );

double evaluatePolynomialSum( const std::map< double, double >& amplitudesAtUnitFrequencyPerPower, const double frequency );

double evaluatePolynomialSquareRootSum( const std::map< double, double >& amplitudesAtUnitFrequencyPerPower, const double frequency );

double evaluatePolynomialSquareRootSumBetweenFrequencyNodes( const std::map< double, double >& amplitudesAtUnitFrequencyPerPower, const std::vector<double> frequencyNodes, const double frequency );

std::function< double( const double ) > createSpectralDensityFunctionFromTabulatedData(
        const std::map< double, double > spectralDensityMap );

std::pair< std::vector< double >, double > generateAmplitudeSpectrumNoise(
        const double maximumFrequency, const int numberOfFrequencySamples, const std::function< double( const double ) > ampitudeSpectalDensityFunction,
        const double seed = statistics::defaultRandomSeedGenerator->getRandomVariableValue( ) );

std::pair< std::vector< double >, double > generateAmplitudeSpectrumNoise(
        const double maximumFrequency, const double frequencyStep, const std::function< double( const double ) > ampitudeSpectalDensityFunction,
        const double seed = statistics::defaultRandomSeedGenerator->getRandomVariableValue( ) );

std::pair< std::vector< double >, double > generatePowerSpectrumNoise(
        const double maximumFrequency, const int numberOfFrequencySamples, const std::function< double( const double ) > powerSpectalDensityFunction,
        const double seed = statistics::defaultRandomSeedGenerator->getRandomVariableValue( ) );

std::pair< std::vector< double >, double > generatePowerSpectrumNoise(
        const double maximumFrequency, const double frequencyStep, const std::function< double( const double ) > powerSpectalDensityFunction,
        const double seed = statistics::defaultRandomSeedGenerator->getRandomVariableValue( ) );

std::pair< std::vector< double >, double > generatePowerLawNoise(
        const double maximumFrequency, const int numberOfFrequencySamples, const std::map< double, double >& powersAndUnitAmplitudes,
        const double seed = statistics::defaultRandomSeedGenerator->getRandomVariableValue( ) );

std::pair< std::vector< double >, double > generatePowerLawNoise(
        const double maximumFrequency, const int numberOfFrequencySamples, const std::map< double, double >& powersAndUnitAmplitudes,
        const std::vector<double> frequencyNodes,
        const double seed = statistics::defaultRandomSeedGenerator->getRandomVariableValue( ) );

//! Function to generate power law noise by discretizing in frequency domain.
std::pair< std::vector< double >, double > generatePowerLawNoise(
        const double maximumFrequency, const double frequencyStep, const std::map< double, double >& powersAndUnitAmplitudes,
        const double seed = statistics::defaultRandomSeedGenerator->getRandomVariableValue( ) );

std::pair< std::vector< double >, double > generatePowerLawNoise(
        const double maximumFrequency, const int numberOfFrequencySamples, const double power, const double amplitudeAtUnityFrequency,
        const double seed = statistics::defaultRandomSeedGenerator->getRandomVariableValue( ) );

//! Function to generate power law noise by discretizing in frequency domain.
std::pair< std::vector< double >, double > generatePowerLawNoise(
        const double maximumFrequency, const double frequencyStep, const double power, const double amplitudeAtUnityFrequency,
        const double seed = statistics::defaultRandomSeedGenerator->getRandomVariableValue( ) );


}

}

#endif // POWERLAWNOISEGENERATION_H
