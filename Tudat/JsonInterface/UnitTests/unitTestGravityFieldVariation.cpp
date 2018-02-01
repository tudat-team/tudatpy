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
#include "Tudat/JsonInterface/Environment/gravityFieldVariation.h"

namespace tudat
{

namespace unit_tests
{

#define INPUT( filename ) \
    ( json_interface::inputDirectory( ) / boost::filesystem::path( __FILE__ ).stem( ) / filename ).string( )

BOOST_AUTO_TEST_SUITE( test_json_gravityFieldVariation )

// Test 1: gravity field variation types
BOOST_AUTO_TEST_CASE( test_json_gravityFieldVariation_bodyDeformationTypes )
{
    BOOST_CHECK_EQUAL_ENUM( INPUT( "bodyDeformationTypes" ),
                            gravitation::bodyDeformationTypes,
                            gravitation::unsupportedBodyDeformationTypes );
}

// Test 2: basic solid body gravity field variation
BOOST_AUTO_TEST_CASE( test_json_gravityFieldVariation_basicSolidBody )
{
    using namespace simulation_setup;
    using namespace json_interface;

    // Create GravityFieldVariationSettings from JSON file
    const boost::shared_ptr< GravityFieldVariationSettings > fromFileSettings =
            parseJSONFile< boost::shared_ptr< GravityFieldVariationSettings > >( INPUT( "basicSolidBody" ) );

    // Create GravityFieldVariationSettings manually
    const std::vector< std::string > deformingBodies = { "Moon" };
    const std::vector< std::vector< std::complex< double > > > loveNumbers =
    {
        { std::complex< double >( 1.0, 2.0 ), std::complex< double >( 2.0, -1.0 ), std::complex< double >( 0.3, -5.0 ) },
        { std::complex< double >( 0.0, 0.5 ), std::complex< double >( 0.0, 2.0 ), std::complex< double >( 4.0, -2.0 ) },
        { std::complex< double >( -3.0, 0.0 ), std::complex< double >( -5.0, 1.0 ), std::complex< double >( 6.0, 0.5 ) }
    };
    const double referenceRadius = 6.4e6;
    const boost::shared_ptr< GravityFieldVariationSettings > manualSettings =
            boost::make_shared< BasicSolidBodyGravityFieldVariationSettings >( deformingBodies,
                                                                               loveNumbers,
                                                                               referenceRadius );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 3: tabulated field variation
BOOST_AUTO_TEST_CASE( test_json_gravityFieldVariation_tabulated )
{
    using namespace interpolators;
    using namespace simulation_setup;
    using namespace json_interface;

    // Create GravityFieldVariationSettings from JSON file
    const boost::shared_ptr< GravityFieldVariationSettings > fromFileSettings =
            parseJSONFile< boost::shared_ptr< GravityFieldVariationSettings > >( INPUT( "tabulated" ) );

    // Create GravityFieldVariationSettings manually
    const std::map< double, Eigen::MatrixXd > cosineCoefficientCorrections =
    {
        { 0.0, ( Eigen::Matrix2d( ) << 0.0, 1.0, -2.0, 1.0 ).finished( ) },
        { 1.0, ( Eigen::Matrix2d( ) << 3.0, 2.5, -1.0, 0.0 ).finished( ) }
    };
    const std::map< double, Eigen::MatrixXd > sineCoefficientCorrections =
    {
        { 0.0, ( Eigen::Matrix2d( ) << -1.0, 4.0, 0.0, 1.0 ).finished( ) },
        { 1.0, ( Eigen::Matrix2d( ) << 3.0, 0.5, -1.0, -2.0 ).finished( ) }
    };
    const int minimumDegree = 4;
    const int minimumOrder = 2;
    const boost::shared_ptr< InterpolatorSettings > interpolatorSettings =
            boost::make_shared< InterpolatorSettings >( cubic_spline_interpolator );
    const boost::shared_ptr< GravityFieldVariationSettings > manualSettings =
            boost::make_shared< TabulatedGravityFieldVariationSettings >( cosineCoefficientCorrections,
                                                                          sineCoefficientCorrections,
                                                                          minimumDegree,
                                                                          minimumOrder,
                                                                          interpolatorSettings );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}



BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
