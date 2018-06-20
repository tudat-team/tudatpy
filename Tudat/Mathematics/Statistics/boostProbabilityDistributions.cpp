/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <boost/random.hpp>
#include <boost/make_shared.hpp>

#include "Tudat/Mathematics/Statistics/boostProbabilityDistributions.h"

namespace tudat
{

namespace statistics
{

//! Function to create a random variable class of BoostContinuousProbabilityDistribution type
std::shared_ptr< InvertibleContinuousProbabilityDistribution< double > > createBoostRandomVariable(
        const ContinuousBoostStatisticalDistributions boostDistribution, const std::vector< double >& parameters )
{
    using namespace boost::math;

    std::shared_ptr< InvertibleContinuousProbabilityDistribution< double > > continuousRandomVariable;

    // Check which distribution type is requested
    switch( boostDistribution )
    {
    case uniform_boost_distribution:
    {
        // Check number of provided parameters
        if ( parameters.size( ) != 2 )
        {
            throw std::runtime_error( "Uniform distribution requires two parameters." );
        }
        else
        {
            // Create uniform distribution. parameters 0: lower bound, 1: upper bound
            uniform_distribution< > uniformDistribution( parameters.at( 0 ), parameters.at( 1 ) );
            continuousRandomVariable =
                    std::make_shared< BoostContinuousProbabilityDistribution< uniform_distribution< > > >(
                        uniformDistribution );
        }
        break;
    }
    case normal_boost_distribution:
    {
        // Check number of provided parameters
        if ( parameters.size( ) != 2 )
        {
            throw std::runtime_error( "Normal distribution requires two parameters." );
        }
        else
        {
            // Create uniform distribution. parameters 0: mean, 1: standard deviation
            normal_distribution< > normalDistribution( parameters.at( 0 ), parameters.at( 1 ) );
            continuousRandomVariable =
                    std::make_shared< BoostContinuousProbabilityDistribution< normal_distribution< > > >(
                        normalDistribution );
        }
        break;
    }
    case exponential_boost_distribution:
    {
        // Check number of provided parameters
        if ( parameters.size( ) != 1 )
        {
            throw std::runtime_error( "Exponential distribution requires one parameter." );
        }
        else
        {
            // Create uniform distribution. parameters 0: lambda parameter
            exponential_distribution< > exponentialDistribution( parameters.at( 0 ) );
            continuousRandomVariable =
                    std::make_shared< BoostContinuousProbabilityDistribution< exponential_distribution< > > >(
                        exponentialDistribution );
        }
        break;
    }
    case gamma_boost_distribution:
    {
        // Check number of provided parameters
        if ( parameters.size( ) != 2 )
        {
            throw std::runtime_error( "Gamma distribution requires two parameters." );
        }
        else
        {
            // Create uniform distribution. parameters 0: shape parameter, parameter 1: scale parameter
            gamma_distribution< > gammaDistribution( parameters.at( 0 ), parameters.at( 1 ) );
            continuousRandomVariable =
                    std::make_shared< BoostContinuousProbabilityDistribution< gamma_distribution< > > >(
                        gammaDistribution );
        }
        break;
    }
    case lognormal_boost_distribution:
    {
        // Check number of provided parameters
        if ( parameters.size( ) != 2 )
        {
            throw std::runtime_error( "Log-normal distribution requires two parameters." );
        }
        else
        {
            // Create uniform distribution. parameters 0: location parameter, parameter 1: scale parameter
            lognormal_distribution< > logNormalDistribution( parameters.at( 0 ), parameters.at( 1 ) );
            continuousRandomVariable =
                    std::make_shared< BoostContinuousProbabilityDistribution< lognormal_distribution< > > >(
                        logNormalDistribution );
        }
        break;
    }
    case beta_boost_distribution:
    {
        // Check number of provided parameters
        if ( parameters.size( ) != 2 )
        {
            throw std::runtime_error( "Beta distribution requires two parameters." );
        }
        else
        {
            // Create uniform distribution. parameters 0: alpha parameter, parameter 1: beta parameter
            beta_distribution< > betaDistribution( parameters.at( 0 ), parameters.at( 1 ) );
            continuousRandomVariable =
                    std::make_shared< BoostContinuousProbabilityDistribution< beta_distribution< > > >(
                        betaDistribution );
        }
        break;
    }
    default:
        throw std::runtime_error( "Boost distribution not recognized." );
    }
    return continuousRandomVariable;
}

}

}
