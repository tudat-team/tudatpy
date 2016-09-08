#include <boost/make_shared.hpp>
#include <boost/bind.hpp>

#include "Tudat/Mathematics/Statistics/randomVariableGenerator.h"
#include "Tudat/Mathematics/Statistics/boostProbabilityDistributions.h"

namespace tudat
{

namespace statistics
{

boost::function< double( ) > createBoostContinuousRandomVariableGeneratorFunction(
        const ContinuousBoostStatisticalDistributions boostDistribution,
        const std::vector< double >& parameters,
        const double seed )
{
    return boost::bind( &RandomVariableGenerator< double >::getRandomVariableValue,
                        createBoostContinuousRandomVariableGenerator( boostDistribution, parameters, seed ) );
}

boost::shared_ptr< RandomVariableGenerator< double > > createBoostContinuousRandomVariableGenerator(
        const ContinuousBoostStatisticalDistributions boostDistribution,
        const std::vector< double >& parameters,
        const double seed )
{
    return boost::make_shared< ContinuousRandomVariableGenerator >(
                createBoostRandomVariable( boostDistribution, parameters ), seed );
}

}

}

