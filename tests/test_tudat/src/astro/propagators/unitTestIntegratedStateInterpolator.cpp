/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include "tudat/math/interpolators/createInterpolator.h"
#include "tudat/math/interpolators/lagrangeInterpolator.h"
#include "tudat/simulation/propagation_setup/propagationProcessingSettings.h"
#include "tudat/simulation/propagation_setup/setNumericallyIntegratedStates.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_integrated_state_interpolator )

BOOST_AUTO_TEST_CASE( testCreateStateInterpolatorOrder )
{
    std::map< double, Eigen::Matrix< double, 6, 1 > > stateMap;
    for( unsigned int i = 0; i < 20; i++ )
    {
        stateMap[ static_cast< double >( i ) ] = Eigen::Matrix< double, 6, 1 >::Constant( static_cast< double >( i ) );
    }

    std::shared_ptr< interpolators::LagrangeInterpolator< double, Eigen::Matrix< double, 6, 1 > > > defaultInterpolator =
            std::dynamic_pointer_cast< interpolators::LagrangeInterpolator< double, Eigen::Matrix< double, 6, 1 > > >(
                    propagators::createStateInterpolator( stateMap ) );
    BOOST_CHECK( defaultInterpolator != nullptr );
    BOOST_CHECK_EQUAL( defaultInterpolator->getNumberOfStages( ), 6 );

    std::shared_ptr< interpolators::LagrangeInterpolator< double, Eigen::Matrix< double, 6, 1 > > > eighthOrderInterpolator =
            std::dynamic_pointer_cast< interpolators::LagrangeInterpolator< double, Eigen::Matrix< double, 6, 1 > > >(
                    propagators::createStateInterpolator( stateMap, interpolators::lagrangeInterpolation( 8 ) ) );
    BOOST_CHECK( eighthOrderInterpolator != nullptr );
    BOOST_CHECK_EQUAL( eighthOrderInterpolator->getNumberOfStages( ), 8 );
}

BOOST_AUTO_TEST_CASE( testProcessingSettingsInterpolatorSettings )
{
    std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings = interpolators::lagrangeInterpolation( 10 );
    propagators::SingleArcPropagatorProcessingSettings processingSettings;
    BOOST_CHECK( processingSettings.getInterpolatorSettings( ) == nullptr );

    processingSettings.setInterpolatorSettings( interpolatorSettings );
    BOOST_CHECK( processingSettings.getInterpolatorSettings( ) == interpolatorSettings );
    BOOST_CHECK_EQUAL(
            std::dynamic_pointer_cast< interpolators::LagrangeInterpolatorSettings >( processingSettings.getInterpolatorSettings( ) )
                    ->getInterpolatorOrder( ),
            10 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
