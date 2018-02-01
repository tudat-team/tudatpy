/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include "Tudat/JsonInterface/UnitTests/unitTestSupport.h"
#include "Tudat/JsonInterface/Mathematics/interpolation.h"

namespace tudat
{

namespace unit_tests
{

#define INPUT( filename ) \
    ( json_interface::inputDirectory( ) / boost::filesystem::path( __FILE__ ).stem( ) / filename ).string( )

BOOST_AUTO_TEST_SUITE( test_json_interpolation )

// Test 1: interpolator types
BOOST_AUTO_TEST_CASE( test_json_interpolation_types )
{
    BOOST_CHECK_EQUAL_ENUM( INPUT( "types" ),
                            interpolators::oneDimensionalInterpolatorTypes,
                            interpolators::unsupportedOneDimensionalInterpolatorTypes );
}

// Test 2: lookup schemes
BOOST_AUTO_TEST_CASE( test_json_interpolation_lookupSchemes )
{
    BOOST_CHECK_EQUAL_ENUM( INPUT( "lookupSchemes" ),
                            interpolators::lookupSchemeTypes,
                            interpolators::unsupportedLookupSchemeTypes );
}

// Test 3: boundary handlings
BOOST_AUTO_TEST_CASE( test_json_interpolation_boundaryHandlings )
{
    BOOST_CHECK_EQUAL_ENUM( INPUT( "boundaryHandlings" ),
                            interpolators::lagrangeInterpolatorBoundaryHandlings,
                            interpolators::unsupportedLagrangeInterpolatorBoundaryHandlings );
}

// Test 4: interpolator
BOOST_AUTO_TEST_CASE( test_json_interpolation_interpolator )
{
    using namespace interpolators;
    using namespace json_interface;

    // Create InterpolatorSettings from JSON file
    const boost::shared_ptr< InterpolatorSettings > fromFileSettings =
            parseJSONFile< boost::shared_ptr< InterpolatorSettings > >( INPUT( "interpolator" ) );

    // Create InterpolatorSettings manually
    const boost::shared_ptr< InterpolatorSettings > manualSettings =
            boost::make_shared< InterpolatorSettings >( piecewise_constant_interpolator, binarySearch, true );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 5: Lagrange interpolator
BOOST_AUTO_TEST_CASE( test_json_interpolation_lagrangeInterpolator )
{
    using namespace interpolators;
    using namespace json_interface;

    // Create InterpolatorSettings from JSON file
    const boost::shared_ptr< InterpolatorSettings > fromFileSettings =
            parseJSONFile< boost::shared_ptr< InterpolatorSettings > >( INPUT( "lagrangeInterpolator" ) );

    // Create InterpolatorSettings manually
    const boost::shared_ptr< InterpolatorSettings > manualSettings =
            boost::make_shared< LagrangeInterpolatorSettings >( 8 );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 6: model interpolation
BOOST_AUTO_TEST_CASE( test_json_interpolation_modelInterpolation )
{
    using namespace interpolators;
    using namespace simulation_setup;
    using namespace json_interface;

    // Create ModelInterpolationSettings from JSON file
    const boost::shared_ptr< ModelInterpolationSettings > fromFileSettings =
            parseJSONFile< boost::shared_ptr< ModelInterpolationSettings > >( INPUT( "modelInterpolation" ) );

    // Create ModelInterpolationSettings manually
    const boost::shared_ptr< ModelInterpolationSettings > manualSettings =
            boost::make_shared< ModelInterpolationSettings >(
                -5.0, 5.0, 0.5, boost::make_shared< InterpolatorSettings >( linear_interpolator ) );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 7: data map
BOOST_AUTO_TEST_CASE( test_json_interpolation_dataMap )
{
    using namespace interpolators;
    using namespace simulation_setup;
    using namespace json_interface;

    // Create DataMapSettings from JSON file
    const boost::shared_ptr< DataMapSettings< double, Eigen::Vector2d > > fromFileSettings =
            parseJSONFile< boost::shared_ptr< DataMapSettings< double, Eigen::Vector2d > > >( INPUT( "dataMap" ) );

    // Create DataMapSettings manually
    std::map< double, Eigen::Vector2d > dataMap;
    dataMap[ 1.0 ] = ( Eigen::Vector2d( ) << 0.5, 1.5 ).finished( );
    dataMap[ 2.0 ] = ( Eigen::Vector2d( ) << 1.5, 2.5 ).finished( );


    const boost::shared_ptr< DataMapSettings< double, Eigen::Vector2d > > manualSettings =
            boost::make_shared< DataMapSettings< double, Eigen::Vector2d > >( dataMap );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 8: from file data map
BOOST_AUTO_TEST_CASE( test_json_interpolation_fromFileDataMap )
{
    using namespace interpolators;
    using namespace simulation_setup;
    using namespace json_interface;

    // Create DataMapSettings from JSON file
    const boost::shared_ptr< DataMapSettings< double, Eigen::Vector2d > > fromFileSettings =
            parseJSONFile< boost::shared_ptr< DataMapSettings< double, Eigen::Vector2d > > >(
                INPUT( "fromFileDataMap" ) );

    // Create DataMapSettings manually
    const boost::shared_ptr< DataMapSettings< double, Eigen::Vector2d > > manualSettings =
            boost::make_shared< FromFileDataMapSettings< Eigen::Vector2d > >( "thrustValues.txt" );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 9: independent-dependent data map
BOOST_AUTO_TEST_CASE( test_json_interpolation_independentDependentDataMap )
{
    using namespace interpolators;
    using namespace simulation_setup;
    using namespace json_interface;

    // Create DataMapSettings from JSON file
    const boost::shared_ptr< DataMapSettings< double, Eigen::Vector2d > > fromFileSettings =
            parseJSONFile< boost::shared_ptr< DataMapSettings< double, Eigen::Vector2d > > >(
                INPUT( "independentDependentDataMap" ) );

    // Create DataMapSettings manually
    std::vector< double > independentVariableValues;
    independentVariableValues.push_back( 1.0 );
    independentVariableValues.push_back( 2.0 );

    std::vector< Eigen::Vector2d > dependentVariableValues;
    dependentVariableValues.push_back( ( Eigen::Vector2d( ) << 0.5, 1.5 ).finished( ) );
    dependentVariableValues.push_back( ( Eigen::Vector2d( ) << 1.5, 2.5 ).finished( ) );

    const boost::shared_ptr< DataMapSettings< double, Eigen::Vector2d > > manualSettings =
            boost::make_shared< IndependentDependentDataMapSettings< double, Eigen::Vector2d > >(
                independentVariableValues, dependentVariableValues );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 10: independent-dependent data map
BOOST_AUTO_TEST_CASE( test_json_interpolation_hermiteDataMap )
{
    using namespace interpolators;
    using namespace simulation_setup;
    using namespace json_interface;

    // Create DataMapSettings from JSON file
    const boost::shared_ptr< DataMapSettings< double, Eigen::Vector2d > > fromFileSettings =
            parseJSONFile< boost::shared_ptr< DataMapSettings< double, Eigen::Vector2d > > >(
                INPUT( "hermiteDataMap" ) );

    // Create DataMapSettings manually
    std::map< double, Eigen::Vector2d > dataMap;
    dataMap[ 1.0 ] = ( Eigen::Vector2d( ) << 0.5, 1.5 ).finished( );
    dataMap[ 2.0 ] = ( Eigen::Vector2d( ) << 1.5, 2.5 ).finished( );

    std::vector< Eigen::Vector2d > dependentVariableFirstDerivativeValues;
    dependentVariableFirstDerivativeValues.push_back( ( Eigen::Vector2d( ) << 1.0, 0.8 ).finished( ) );
    dependentVariableFirstDerivativeValues.push_back( ( Eigen::Vector2d( ) << 0.5, 0.4 ).finished( ) );

    const boost::shared_ptr< DataMapSettings< double, Eigen::Vector2d > > manualSettings =
            boost::make_shared< HermiteDataSettings< double, Eigen::Vector2d > >(
                dataMap, dependentVariableFirstDerivativeValues );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 11: independent-dependent data map
BOOST_AUTO_TEST_CASE( test_json_interpolation_dataInterpolation )
{
    using namespace interpolators;
    using namespace simulation_setup;
    using namespace json_interface;

    // Create DataInterpolationSettings from JSON file
    const boost::shared_ptr< DataInterpolationSettings< double, Eigen::Vector2d > > fromFileSettings =
            parseJSONFile< boost::shared_ptr< DataInterpolationSettings< double, Eigen::Vector2d > > >(
                INPUT( "dataInterpolation" ) );

    // Create DataInterpolationSettings manually
    const boost::shared_ptr< DataInterpolationSettings< double, Eigen::Vector2d > > manualSettings =
            boost::make_shared< DataInterpolationSettings< double, Eigen::Vector2d > >(
                boost::make_shared< FromFileDataMapSettings< Eigen::Vector2d > >( "thrustValues.txt" ),
                boost::make_shared< LagrangeInterpolatorSettings >( 4 ) );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
