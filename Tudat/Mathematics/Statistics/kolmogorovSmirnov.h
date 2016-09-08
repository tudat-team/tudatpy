#ifndef TUDAT_KOLMOGOROVSMIRNOV_H
#define TUDAT_KOLMOGOROVSMIRNOV_H

#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>

#include "Tudat/Mathematics/Statistics/continuousProbabilityDistributions.h"

namespace tudat
{

namespace statistics
{

boost::function< double( const double ) > createEmpiricalCdfFunction( const std::vector< double >& sortedData );


bool performKolmogorovSmirnovTest( const ContinuousBoostStatisticalDistributions distributionType,
                                   const std::vector< double >& distributionParameters,
                                   const std::vector< double >& dataToTest,
                                   const double toleranceCriterion );

bool performKolmogorovSmirnovTest(
        const boost::shared_ptr< ContinuousProbabilityDistribution< double > > expectedDistribution,
        std::vector< double > dataToTest,
        const double toleranceCriterion );

}

}

#endif // TUDAT_KOLMOGOROVSMIRNOV_H
