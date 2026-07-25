/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "mroDsnObservationModelTestHelpers.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::simulation_setup;

BOOST_AUTO_TEST_SUITE( test_dsn_n_way_range_observation_model )

BOOST_AUTO_TEST_CASE( testMroTrk234DsnNWayRangeModel )
{
    using namespace tudat::unit_tests::mro_dsn_test;

    loadMroSpiceKernels( );

    // The MRO DSN range fixture is generated from mromagr2012_076_0840xmmmv1.tnf over the same
    // one-hour interval as the Doppler fixture, at the Tudat trk234 converter boundary. The CSV
    // stores the exact values that are otherwise passed directly to SingleObservationSet creation.
    std::shared_ptr< observation_models::ObservationCollection< long double, Time > > observedObservationCollection =
            createObservationCollectionFromTrk234Csv( trk234InputsDirectory + "range_single_observation_set_inputs.csv",
                                                      observation_models::dsn_n_way_range );

    std::pair< Time, Time > timeBounds = observedObservationCollection->getTimeBounds( );
    SystemOfBodies bodies = createMroSystemOfBodies( timeBounds.first, timeBounds.second );
    setRampFrequencyInterpolatorsInBodies( bodies );
    applyMroNotebookObservationCollectionPostProcessing( observedObservationCollection, bodies );

    Eigen::VectorXd residuals = simulateAndGetResiduals( observedObservationCollection, bodies, true );
    BOOST_CHECK_EQUAL( residuals.rows( ), 18 );
    BOOST_TEST_MESSAGE( "MRO DSN range residuals: " << residuals.transpose( ) );

    std::vector< double > absoluteResiduals;
    for( Eigen::Index i = 0; i < residuals.rows( ); ++i )
    {
        absoluteResiduals.push_back( std::fabs( residuals( i ) ) );
    }
    std::sort( absoluteResiduals.begin( ), absoluteResiduals.end( ) );
    const double medianAbsoluteResidual =
            0.5 * ( absoluteResiduals.at( absoluteResiduals.size( ) / 2 - 1 ) + absoluteResiduals.at( absoluteResiduals.size( ) / 2 ) );
    const double residualThreshold = 3.0 * medianAbsoluteResidual;

    std::vector< double > filteredResiduals;
    for( Eigen::Index i = 0; i < residuals.rows( ); ++i )
    {
        if( std::fabs( residuals( i ) ) <= residualThreshold )
        {
            filteredResiduals.push_back( residuals( i ) );
        }
    }

    double filteredResidualSum = 0.0;
    double filteredResidualSquaredSum = 0.0;
    for( const double residual : filteredResiduals )
    {
        filteredResidualSum += residual;
        filteredResidualSquaredSum += residual * residual;
    }

    const double filteredMeanResidual = filteredResidualSum / static_cast< double >( filteredResiduals.size( ) );
    const double filteredRmsResidual = std::sqrt( filteredResidualSquaredSum / static_cast< double >( filteredResiduals.size( ) ) );
    BOOST_TEST_MESSAGE( "MRO DSN range median absolute residual: " << medianAbsoluteResidual << " RU" );
    BOOST_TEST_MESSAGE( "MRO DSN range retained observations: " << filteredResiduals.size( ) << " / " << residuals.rows( ) );
    BOOST_TEST_MESSAGE( "MRO DSN range filtered residual mean: " << filteredMeanResidual << " RU" );
    BOOST_TEST_MESSAGE( "MRO DSN range filtered residual RMS: " << filteredRmsResidual << " RU" );

    BOOST_TEST( std::fabs( filteredMeanResidual ) < 10.0 );
    BOOST_TEST( filteredRmsResidual < 10.0 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
