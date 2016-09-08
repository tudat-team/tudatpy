/*    Copyright (c) 2010-2016, Delft University of Technology
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

#include "Tudat/Mathematics/Statistics/kolmogorovSmirnov.h"
#include "Tudat/Mathematics/Statistics/boostProbabilityDistributions.h"
#include "Tudat/Mathematics/Interpolators/piecewiseConstantInterpolator.h"


namespace tudat
{

namespace statistics
{

boost::function< double( const double ) > createEmpiricalCdfFunction( const std::vector< double >& sortedData )
{

    std::vector< double > cdfValues;
    cdfValues.resize( sortedData.size( ) );

    double numberOfDataPoints = static_cast< double >( cdfValues.size( ) );

    for( unsigned int i = 0; i < sortedData.size( ); i++ )
    {
        cdfValues[ i ] = ( static_cast< double >( i ) + 0.5 ) / numberOfDataPoints;
    }


    typedef interpolators::OneDimensionalInterpolator< double, double > LocalInterpolator;

    boost::shared_ptr< LocalInterpolator > cdfInterpolator =
            boost::make_shared< interpolators::PiecewiseConstantInterpolator< double, double > >( sortedData, cdfValues );

    return boost::bind(
                static_cast< double( LocalInterpolator::* )( const double ) >
                ( &LocalInterpolator::interpolate ), cdfInterpolator, _1 );
}

std::pair< bool, double > performKolmogorovSmirnovTest(
        const tudat::statistics::ContinuousBoostStatisticalDistributions distributionType,
        const std::vector< double >& distributionParameters,
        const std::vector< double >& dataToTest,
        const double toleranceCriterion )
{
    return performKolmogorovSmirnovTest( createBoostRandomVariable( distributionType, distributionParameters ),
                                         dataToTest, toleranceCriterion );
}

std::pair< bool, double > performKolmogorovSmirnovTest(
        const boost::shared_ptr< ContinuousProbabilityDistribution< double > > expectedDistribution,
        std::vector< double > dataToTest,
        const double toleranceCriterion )
{
    boost::function< double( const double ) > empiricalCdfFunction =
            createEmpiricalCdfFunction( dataToTest );

    double maximumDifference = 0.0;

    std::sort( dataToTest.begin( ), dataToTest.end( ) );

    double currentDataPoint, currentDifference;
    for( unsigned int i = 0; i < dataToTest.size( ) - 1; i++ )
    {
        currentDataPoint = dataToTest[ i ] + ( dataToTest[ i + 1 ] - dataToTest[ i ] ) / 2.0;

        currentDifference = empiricalCdfFunction( currentDataPoint ) -
                expectedDistribution->evaluateCdf( currentDataPoint );

        if( std::fabs( currentDifference ) > maximumDifference )
        {
            maximumDifference = std::fabs( currentDifference );
        }
    }

    bool isFailed = 1;

    if( maximumDifference > toleranceCriterion )
    {
        isFailed = 1;
    }
    else
    {
        isFailed = 0;
    }

    return std::make_pair( isFailed, maximumDifference );
}

}

}
