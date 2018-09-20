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
#include "Tudat/JsonInterface/Environment/aerodynamics.h"

namespace tudat
{

namespace unit_tests
{

#define INPUT( filename ) \
    ( json_interface::inputDirectory( ) / boost::filesystem::path( __FILE__ ).stem( ) / filename ).string( )

BOOST_AUTO_TEST_SUITE( test_json_aerodynamics )

// Test 1: aerodynamic coefficients types
BOOST_AUTO_TEST_CASE( test_json_aerodynamics_coefficientsTypes )
{
    BOOST_CHECK_EQUAL_ENUM( INPUT( "coefficientsTypes" ),
                            simulation_setup::aerodynamicCoefficientTypes,
                            simulation_setup::unsupportedAerodynamicCoefficientTypes );
}

// Test 2: aerodynamic variables
BOOST_AUTO_TEST_CASE( test_json_aerodynamics_variables )
{
    BOOST_CHECK_EQUAL_ENUM( INPUT( "variables" ),
                            aerodynamics::aerodynamicVariables,
                            aerodynamics::unsupportedAerodynamicVariables );
}

// Test 3: constant aerodynamics (only drag coefficient)
BOOST_AUTO_TEST_CASE( test_json_aerodynamics_dragCoefficient )
{
    using namespace simulation_setup;
    using namespace json_interface;

    // Create AerodynamicCoefficientSettings from JSON file
    const std::shared_ptr< AerodynamicCoefficientSettings > fromFileSettings =
            parseJSONFile< std::shared_ptr< AerodynamicCoefficientSettings > >( INPUT( "dragCoefficient" ) );

    // Create AerodynamicCoefficientSettings manually
    const double referenceArea = 10.5;
    const double dragCoefficient = 2.2;
    Eigen::Vector3d forceCoefficients = Eigen::Vector3d::Zero( );
    forceCoefficients( 0 ) = dragCoefficient;
    const std::shared_ptr< AerodynamicCoefficientSettings > manualSettings =
            std::make_shared< ConstantAerodynamicCoefficientSettings >( referenceArea,
                                                                          forceCoefficients );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 4: constant aerodynamics (full)
BOOST_AUTO_TEST_CASE( test_json_aerodynamics_constant )
{
    using namespace simulation_setup;
    using namespace json_interface;

    // Create AerodynamicCoefficientSettings from JSON file
    const std::shared_ptr< AerodynamicCoefficientSettings > fromFileSettings =
            parseJSONFile< std::shared_ptr< AerodynamicCoefficientSettings > >( INPUT( "constant" ) );

    // Create AerodynamicCoefficientSettings manually
    const double referenceLength = 5.0;
    const double referenceArea = 10.5;
    const double lateralReferenceLength = 4.0;
    const Eigen::Vector3d momentReferencePoint = ( Eigen::Vector3d( ) << 0.7, 0.8, 0.9 ).finished( );
    const Eigen::Vector3d forceCoefficients = ( Eigen::Vector3d( ) << 1.0, 2.0, 3.0 ).finished( );
    const Eigen::Vector3d momentCoefficients = ( Eigen::Vector3d( ) << 0.0, 1.0e-3, -0.1 ).finished( );
    const bool areCoefficientsInAerodynamicFrame = true;
    const bool areCoefficientsInNegativeAxisDirection = false;
    const std::shared_ptr< AerodynamicCoefficientSettings > manualSettings =
            std::make_shared< ConstantAerodynamicCoefficientSettings >( referenceLength,
                                                                          referenceArea,
                                                                          lateralReferenceLength,
                                                                          momentReferencePoint,
                                                                          forceCoefficients,
                                                                          momentCoefficients,
                                                                          areCoefficientsInAerodynamicFrame,
                                                                          areCoefficientsInNegativeAxisDirection );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 5: tabulated aerodynamics (1 dimension)
BOOST_AUTO_TEST_CASE( test_json_aerodynamics_tabulated1 )
{
    using namespace aerodynamics;
    using namespace interpolators;
    using namespace simulation_setup;
    using namespace json_interface;

    // Create AerodynamicCoefficientSettings from JSON file
    const std::shared_ptr< AerodynamicCoefficientSettings > fromFileSettings =
            parseJSONFile< std::shared_ptr< AerodynamicCoefficientSettings > >( INPUT( "tabulated1" ) );

    // Create AerodynamicCoefficientSettings manually
    const std::vector< double > independentVariables = { 0.0, 1.0, 2.0, 3.0 };
    const std::vector< Eigen::Vector3d > forceCoefficients =
    {
        ( Eigen::Vector3d( ) << 0.7, 0.8, 0.9 ).finished( ),
        ( Eigen::Vector3d( ) << 1.7, 1.8, 1.9 ).finished( ),
        ( Eigen::Vector3d( ) << 2.7, 2.8, 2.9 ).finished( ),
        ( Eigen::Vector3d( ) << 3.7, 3.8, 3.9 ).finished( )
    };
    const std::vector< Eigen::Vector3d > momentCoefficients =
    {
        ( Eigen::Vector3d( ) << 1.0, 2.0, 3.0 ).finished( ),
        ( Eigen::Vector3d( ) << 1.0, 1.0, 1.0 ).finished( ),
        ( Eigen::Vector3d( ) << 2.0, 2.0, 2.0 ).finished( ),
        ( Eigen::Vector3d( ) << 3.0, 3.0, 3.0 ).finished( )
    };
    const double referenceLength = 5.0;
    const double referenceArea = 10.5;
    const double lateralReferenceLength = 4.0;
    const Eigen::Vector3d momentReferencePoint = ( Eigen::Vector3d( ) << 0.7, 0.8, 0.9 ).finished( );
    const AerodynamicCoefficientsIndependentVariables independentVariableName = angle_of_sideslip_dependent;
    const std::shared_ptr< InterpolatorSettings > interpolatorSettings =
            std::make_shared< InterpolatorSettings >( cubic_spline_interpolator );
    const bool areCoefficientsInAerodynamicFrame = false;
    const bool areCoefficientsInNegativeAxisDirection = false;
    const std::shared_ptr< AerodynamicCoefficientSettings > manualSettings =
            std::make_shared< TabulatedAerodynamicCoefficientSettings< 1 > >( independentVariables,
                                                                                forceCoefficients,
                                                                                momentCoefficients,
                                                                                referenceLength,
                                                                                referenceArea,
                                                                                lateralReferenceLength,
                                                                                momentReferencePoint,
                                                                                independentVariableName,
                                                                                areCoefficientsInAerodynamicFrame,
                                                                                areCoefficientsInNegativeAxisDirection,
                                                                                interpolatorSettings );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 6: tabulated aerodynamics (N dimensions)
BOOST_AUTO_TEST_CASE( test_json_aerodynamics_tabulatedN )
{
    using namespace aerodynamics;
    using namespace interpolators;
    using namespace simulation_setup;
    using namespace json_interface;

    // Create AerodynamicCoefficientSettings from JSON file
    const std::shared_ptr< AerodynamicCoefficientSettings > fromFileSettings =
            parseJSONFile< std::shared_ptr< AerodynamicCoefficientSettings > >( INPUT( "tabulatedN" ) );

    // Create AerodynamicCoefficientSettings manually
    const std::map< int, std::string > forceCoefficientsFiles = { { 0, "aurora_CD.txt" }, { 2, "aurora_CL.txt" } };
    const std::map< int, std::string > momentCoefficientsFiles = { { 1, "aurora_Cm.txt" } };
    const double referenceLength = 5.0;
    const double referenceArea = 10.5;
    const double lateralReferenceLength = 4.0;
    const Eigen::Vector3d momentReferencePoint = ( Eigen::Vector3d( ) << 0.7, 0.8, 0.9 ).finished( );
    const std::vector< AerodynamicCoefficientsIndependentVariables > independentVariableNames =
    { mach_number_dependent, angle_of_attack_dependent };
    const bool areCoefficientsInAerodynamicFrame = true;
    const bool areCoefficientsInNegativeAxisDirection = true;
    const std::shared_ptr< AerodynamicCoefficientSettings > manualSettings =
            readTabulatedAerodynamicCoefficientsFromFiles(
                                        forceCoefficientsFiles,
                                        momentCoefficientsFiles,
                                        referenceLength,
                                        referenceArea,
                                        lateralReferenceLength,
                                        momentReferencePoint,
                                        independentVariableNames,
                                        areCoefficientsInAerodynamicFrame,
                                        areCoefficientsInNegativeAxisDirection );

    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
