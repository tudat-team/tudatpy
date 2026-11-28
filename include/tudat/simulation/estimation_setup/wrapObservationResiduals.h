/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_WRAPOBSERVATIONRESIDUALS_H
#define TUDAT_WRAPOBSERVATIONRESIDUALS_H

#include <cmath>
#include <limits>
#include <utility>
#include <vector>

#include <Eigen/Core>

#include "tudat/astro/observation_models/observableTypes.h"

namespace tudat
{

namespace simulation_setup
{

template< typename ObservationScalarType >
void wrapObservationResiduals(
        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& residuals,
        const std::pair< int, int > observableResidualStartAndSize,
        const observation_models::ObservableType observableType,
        const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& observedObservationBlock =
                Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >( ),
        const observation_models::ResidualWrappingSettings& residualWrappingSettings = observation_models::ResidualWrappingSettings( ) )
{
    const int residualBlockSize = observableResidualStartAndSize.second;
    if( residualBlockSize == 0 )
    {
        return;
    }

    const int singleObservationSize = observation_models::getObservableSize( observableType );
    const int numberOfObservations = residualBlockSize / singleObservationSize;

    const std::vector< int >& wrappedComponentIndices = observation_models::getResidualWrappingComponentIndices( observableType );
    if( wrappedComponentIndices.empty( ) )
    {
        return;
    }

    const std::vector< observation_models::ResidualWrappingRange > wrappingRanges =
            observation_models::getResidualWrappingRanges( observableType );

    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > currentResidualBlock =
            residuals.block( observableResidualStartAndSize.first, 0, residualBlockSize, 1 );

    // Wrap each periodic component per observation
    for( int currentObservationIndex = 0; currentObservationIndex < numberOfObservations; currentObservationIndex++ )
    {
        for( const int currentObservableComponentIndex : wrappedComponentIndices )
        {
            const int currentResidualComponentIndex = currentObservationIndex * singleObservationSize + currentObservableComponentIndex;
            double residualWrappingPeriod = wrappingRanges[ currentObservableComponentIndex ].period( );
            double residualWrappingCenter = wrappingRanges[ currentObservableComponentIndex ].center( );

            if( observableType == observation_models::angular_position && residualWrappingSettings.normalizeRightAscension &&
                currentObservableComponentIndex == 0 )
            {
                const double rightAscensionNormalizationFactor =
                        std::abs( std::cos( static_cast< double >( observedObservationBlock( currentResidualComponentIndex + 1 ) ) ) );
                if( rightAscensionNormalizationFactor <= 10.0 * std::numeric_limits< double >::epsilon( ) )
                {
                    continue;
                }
                residualWrappingPeriod *= rightAscensionNormalizationFactor;
                residualWrappingCenter *= rightAscensionNormalizationFactor;
            }

            currentResidualBlock( currentResidualComponentIndex, 0 ) = currentResidualBlock( currentResidualComponentIndex, 0 ) -
                    residualWrappingPeriod *
                            std::round( ( currentResidualBlock( currentResidualComponentIndex, 0 ) - residualWrappingCenter ) /
                                        residualWrappingPeriod );
        }
    }

    residuals.block( observableResidualStartAndSize.first, 0, residualBlockSize, 1 ) = currentResidualBlock;
}

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_WRAPOBSERVATIONRESIDUALS_H
