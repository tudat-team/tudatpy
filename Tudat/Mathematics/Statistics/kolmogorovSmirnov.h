/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


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


std::pair< bool, double > performKolmogorovSmirnovTest(
        const ContinuousBoostStatisticalDistributions distributionType,
        const std::vector< double >& distributionParameters,
        const std::vector< double >& dataToTest,
        const double toleranceCriterion );

std::pair< bool, double > performKolmogorovSmirnovTest(
        const boost::shared_ptr< ContinuousProbabilityDistribution< double > > expectedDistribution,
        std::vector< double > dataToTest,
        const double toleranceCriterion );

}

}

#endif // TUDAT_KOLMOGOROVSMIRNOV_H
