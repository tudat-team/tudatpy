/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <boost/make_shared.hpp>
#include <boost/bind.hpp>

#include "Tudat/Mathematics/Statistics/randomVariableGenerator.h"
#include "Tudat/Mathematics/Statistics/boostProbabilityDistributions.h"

namespace tudat
{

namespace statistics
{

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

