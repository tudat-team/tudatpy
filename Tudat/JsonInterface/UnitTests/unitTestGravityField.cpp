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
#include "Tudat/JsonInterface/Environment/gravityField.h"

namespace tudat
{

namespace unit_tests
{

#define INPUT( filename ) \
    ( json_interface::inputDirectory( ) / boost::filesystem::path( __FILE__ ).stem( ) / filename ).string( )

BOOST_AUTO_TEST_SUITE( test_json_gravityField )

// Test 1: gravity field types
BOOST_AUTO_TEST_CASE( test_json_gravityField_types )
{
    BOOST_CHECK_EQUAL_ENUM( INPUT( "types" ),
                            simulation_setup::gravityFieldTypes,
                            simulation_setup::unsupportedGravityFieldTypes );
}

// Test 2: spherical harmonic models
BOOST_AUTO_TEST_CASE( test_json_gravityField_sphericalHarmonicModels )
{
    BOOST_CHECK_EQUAL_ENUM( INPUT( "sphericalHarmonicModels" ),
                            simulation_setup::sphericalHarmonicsModels,
                            simulation_setup::unsupportedSphericalHarmonicsModels );
}

// Test 3: point mass gravity field
BOOST_AUTO_TEST_CASE( test_json_gravityField_pointMass )
{
    using namespace simulation_setup;
    using namespace json_interface;

    // Create GravityFieldSettings from JSON file
    const std::shared_ptr< GravityFieldSettings > fromFileSettings =
            parseJSONFile< std::shared_ptr< GravityFieldSettings > >( INPUT( "pointMass" ) );

    // Create GravityFieldSettings manually
    const double gravitationalParameter = 4.0e14;
    const std::shared_ptr< GravityFieldSettings > manualSettings =
            std::make_shared< CentralGravityFieldSettings >( gravitationalParameter );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 4: point mass Spice gravity field
BOOST_AUTO_TEST_CASE( test_json_gravityField_pointMassSpice )
{
    using namespace simulation_setup;
    using namespace json_interface;

    // Create GravityFieldSettings from JSON file
    const std::shared_ptr< GravityFieldSettings > fromFileSettings =
            parseJSONFile< std::shared_ptr< GravityFieldSettings > >( INPUT( "pointMassSpice" ) );

    // Create GravityFieldSettings manually
    const std::shared_ptr< GravityFieldSettings > manualSettings =
            std::make_shared< GravityFieldSettings >( central_spice );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 5: spherical harmonic gravity field (from named model)
BOOST_AUTO_TEST_CASE( test_json_gravityField_sphericalHarmonic_model )
{
    using namespace simulation_setup;
    using namespace json_interface;

    // Create GravityFieldSettings from JSON file
    const std::shared_ptr< GravityFieldSettings > fromFileSettings =
            parseJSONFile< std::shared_ptr< GravityFieldSettings > >( INPUT( "sphericalHarmonic_model" ) );

    // Create GravityFieldSettings manually
    const std::shared_ptr< GravityFieldSettings > manualSettings =
            std::make_shared< FromFileSphericalHarmonicsGravityFieldSettings >( ggm02c );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 6: spherical harmonic gravity field (from file)
BOOST_AUTO_TEST_CASE( test_json_gravityField_sphericalHarmonic_file )
{
    using namespace simulation_setup;
    using namespace json_interface;

    // Create GravityFieldSettings from JSON file
    const std::shared_ptr< GravityFieldSettings > fromFileSettings =
            parseJSONFile< std::shared_ptr< GravityFieldSettings > >( INPUT( "sphericalHarmonic_file" ) );

    // Create GravityFieldSettings manually
    const std::string file = "sh.txt";
    const std::string associatedReferenceFrame = "IAU_Earth";
    const unsigned int maximumDegree = 2;
    const unsigned int maximumOrder = 1;
    const int gravitationalParameterIndex = 0;
    const int referenceRadiusIndex = 1;
    const std::shared_ptr< GravityFieldSettings > manualSettings =
            std::make_shared< FromFileSphericalHarmonicsGravityFieldSettings >( file,
                                                                               associatedReferenceFrame,
                                                                               maximumDegree,
                                                                               maximumOrder,
                                                                               gravitationalParameterIndex,
                                                                               referenceRadiusIndex );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 7: spherical harmonic gravity field (from file, manual parameters)
BOOST_AUTO_TEST_CASE( test_json_gravityField_sphericalHarmonic_file_manualparam )
{
    using namespace simulation_setup;
    using namespace json_interface;

    // Create GravityFieldSettings from JSON file
    const std::shared_ptr< GravityFieldSettings > fromFileSettings =
            parseJSONFile< std::shared_ptr< GravityFieldSettings > >( INPUT( "sphericalHarmonic_file_manualparam" ) );

    // Create GravityFieldSettings manually
    const std::string file = "sh_manualparam.txt";
    const std::string associatedReferenceFrame = "IAU_Earth";
    const unsigned int maximumDegree = 2;
    const unsigned int maximumOrder = 1;
    const int gravitationalParameterIndex = -1;
    const int referenceRadiusIndex = -1;
    const double gravitationalParameter = 4.0e14;
    const double referenceRadius = 6.4e6;
    const std::shared_ptr< GravityFieldSettings > manualSettings =
            std::make_shared< FromFileSphericalHarmonicsGravityFieldSettings >( file,
                                                                               associatedReferenceFrame,
                                                                               maximumDegree,
                                                                               maximumOrder,
                                                                               gravitationalParameterIndex,
                                                                               referenceRadiusIndex,
                                                                               gravitationalParameter,
                                                                               referenceRadius );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 8: spherical harmonic gravity field (direct)
BOOST_AUTO_TEST_CASE( test_json_gravityField_sphericalHarmonic_direct )
{
    using namespace simulation_setup;
    using namespace json_interface;

    // Create GravityFieldSettings from JSON file
    const std::shared_ptr< GravityFieldSettings > fromFileSettings =
            parseJSONFile< std::shared_ptr< GravityFieldSettings > >( INPUT( "sphericalHarmonic_direct" ) );

    // Create GravityFieldSettings manually
    const double gravitationalParameter = 4.0e14;
    const double referenceRadius = 6.4e6;
    Eigen::MatrixXd cosineCoefficients = Eigen::MatrixXd::Zero( 2, 2 );
    cosineCoefficients( 0, 0 ) = 1.0;
    const Eigen::MatrixXd sisineCoefficients = Eigen::MatrixXd::Zero( 2, 2 );
    const std::string associatedReferenceFrame = "IAU_Earth";
    const std::shared_ptr< GravityFieldSettings > manualSettings =
            std::make_shared< SphericalHarmonicsGravityFieldSettings >( gravitationalParameter,
                                                                          referenceRadius,
                                                                          cosineCoefficients,
                                                                          sisineCoefficients,
                                                                          associatedReferenceFrame );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
