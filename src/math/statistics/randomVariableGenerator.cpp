/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <boost/make_shared.hpp>
#include <boost/bind/bind.hpp>
using namespace boost::placeholders;


#include "tudat/math/statistics/randomVariableGenerator.h"
#include "tudat/math/statistics/boostProbabilityDistributions.h"

namespace tudat
{

namespace statistics
{

std::function< Eigen::VectorXd( const double ) > getIndependentGaussianNoiseFunction(
        const double standardDeviation,
        const double mean,
        const double seed,
        const int outputSize )
{
    std::function< double( ) > inputFreeNoiseFunction = statistics::createBoostContinuousRandomVariableGeneratorFunction(
                statistics::normal_boost_distribution, { mean, standardDeviation }, seed );
    if( outputSize == 1 )
    {
        return [=](const double){ return ( Eigen::VectorXd( outputSize )<<
                                           inputFreeNoiseFunction( ) ).finished( ); };
    }
    else if( outputSize == 2 )
    {
        return [=](const double){ return ( Eigen::VectorXd( outputSize )<<
                                           inputFreeNoiseFunction( ), inputFreeNoiseFunction( ) ).finished( ); };
    }
    else if( outputSize == 3 )
    {
        return [=](const double){ return ( Eigen::VectorXd( outputSize )<<
                                           inputFreeNoiseFunction( ), inputFreeNoiseFunction( ), inputFreeNoiseFunction( ) ).finished( ); };
    }
    else if( outputSize == 6 )
    {
        return [=](const double){ return ( Eigen::VectorXd( outputSize )<<
                                           inputFreeNoiseFunction( ), inputFreeNoiseFunction( ), inputFreeNoiseFunction( ),
                                           inputFreeNoiseFunction( ), inputFreeNoiseFunction( ), inputFreeNoiseFunction( ) ).finished( ); };
    }
    else
    {
        throw std::runtime_error( "Cannot simulate observation noise of size " + std::to_string( outputSize ) );
    }
}

//! Function to create a random number generating function from a continuous univariate distribution implemented in boost
std::function< double( ) > createBoostContinuousRandomVariableGeneratorFunction(
        const ContinuousBoostStatisticalDistributions boostDistribution,
        const std::vector< double >& parameters,
        const double seed )
{
    return std::bind( &RandomVariableGenerator< double >::getRandomVariableValue,
                        createBoostContinuousRandomVariableGenerator( boostDistribution, parameters, seed ) );
}

//! Function to create a random number generator from a continuous univariate distribution implemented in boost
std::shared_ptr< RandomVariableGenerator< double > > createBoostContinuousRandomVariableGenerator(
        const ContinuousBoostStatisticalDistributions boostDistribution,
        const std::vector< double >& parameters,
        const double seed )
{
    return std::make_shared< ContinuousRandomVariableGenerator >(
                createBoostRandomVariable( boostDistribution, parameters ), seed );
}

}

}

