/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "tudat/astro/orbit_determination/metric_partials/metricPartial.h"

namespace tudat
{

namespace orbit_determination
{

namespace partial_derivatives
{

std::vector< Eigen::Matrix< double, 4, 4 > > getMetricPartialWrtGroundStationPosition(
        const std::vector< Eigen::Matrix< double, 4, 4 > >& metricPartialWrtBodyPosition,
        const Eigen::Matrix3d& rotationToLocalFrame )
{
    std::vector< Eigen::Matrix< double, 4, 4 > > partialsWrtGroundStationPosition( 3 );

    for( int i = 0; i < 3; ++i )
    {
        partialsWrtGroundStationPosition[ i ].setZero( );
        for( int j = 0; j < 3; ++j )
        {
            partialsWrtGroundStationPosition[ i ] += metricPartialWrtBodyPosition.at( j ) * rotationToLocalFrame( j, i );
        }
    }

    return partialsWrtGroundStationPosition;
}

std::vector< Eigen::Matrix< double, 4, 4 > > calculateCurrentChristoffelSymbolsFromMetricPartials(
        const Eigen::Matrix< double, 6, 1 >& currentState,
        const double currentTime,
        const std::pair< std::string, std::string >& metricEvaluationPoint,
        const std::shared_ptr< relativity::Metric > metric,
        const std::shared_ptr< MetricPartial > metricPartial )
{
    metric->update( currentState, currentTime, true, false );
    metricPartial->update( );

    const Eigen::Matrix< double, 4, 4 > contravariantMetric = metric->getCurrentContravariantMetric( );

    std::vector< Eigen::Matrix< double, 4, 4 > > metricPartials =
            metricPartial->wrtStateOfIntegratedBody( metricEvaluationPoint, propagators::translational_state );

    std::vector< Eigen::Matrix< double, 4, 4 > > spaceTimeMetricPartials;
    spaceTimeMetricPartials.push_back( metricPartial->wrtScaledTime( ) );
    spaceTimeMetricPartials.insert(
                spaceTimeMetricPartials.end( ), metricPartials.begin( ), metricPartials.end( ) );

    std::vector< Eigen::Matrix< double, 4, 4 > > christoffelSymbols( 4 );

    for( int mu = 0; mu < 4; ++mu )
    {
        christoffelSymbols[ mu ].setZero( );
        const Eigen::Matrix< double, 1, 4 > currentMetricRow = contravariantMetric.block( mu, 0, 1, 4 );

        for( int alpha = 0; alpha < 4; ++alpha )
        {
            for( int beta = 0; beta < 4; ++beta )
            {
                christoffelSymbols[ mu ]( alpha, beta ) += 0.5 *
                        ( currentMetricRow *
                          ( spaceTimeMetricPartials[ beta ].block( 0, alpha, 4, 1 ) +
                            spaceTimeMetricPartials[ alpha ].block( 0, beta, 4, 1 ) ) )( 0, 0 );
            }
        }

        for( int nu = 0; nu < 4; ++nu )
        {
            christoffelSymbols[ mu ] -= 0.5 * contravariantMetric( mu, nu ) * spaceTimeMetricPartials[ nu ];
        }
    }

    return christoffelSymbols;
}

}  // namespace partial_derivatives

}  // namespace orbit_determination

}  // namespace tudat

