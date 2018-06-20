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
#include "Tudat/JsonInterface/Propagation/propagator.h"

namespace tudat
{

namespace unit_tests
{

#define INPUT( filename ) \
    ( json_interface::inputDirectory( ) / boost::filesystem::path( __FILE__ ).stem( ) / filename ).string( )

BOOST_AUTO_TEST_SUITE( test_json_propagator )

// Test 1: propagator types
BOOST_AUTO_TEST_CASE( test_json_propagator_types )
{
    BOOST_CHECK_EQUAL_ENUM( INPUT( "types" ),
                            propagators::integratedStateTypes,
                            propagators::unsupportedIntegratedStateTypes );
}

// Test 2: translational propagator types
BOOST_AUTO_TEST_CASE( test_json_propagator_translationalTypes )
{
    BOOST_CHECK_EQUAL_ENUM( INPUT( "translationalTypes" ),
                            propagators::translationalPropagatorTypes,
                            propagators::unsupportedTranslationalPropagatorTypes );
}

// Test 3: translational propagator
BOOST_AUTO_TEST_CASE( test_json_propagator_translational )
{
    using namespace basic_astrodynamics;
    using namespace propagators;
    using namespace simulation_setup;
    using namespace json_interface;

    // Create TranslationalStatePropagatorSettings from JSON file
    const std::shared_ptr< SingleArcPropagatorSettings< double > > fromFileSettings =
            parseJSONFile< std::shared_ptr< SingleArcPropagatorSettings< double > > >( INPUT( "translational" ) );

    // Create TranslationalStatePropagatorSettings manually
    const TranslationalPropagatorType propagatorType = encke;
    const std::vector< std::string > bodiesToPropagate = { "a", "b" };
    const std::vector< std::string > centralBodies = { "A", "B" };
    const std::shared_ptr< PropagationTerminationSettings > terminationSettings;
    const Eigen::VectorXd initialStates = ( Eigen::VectorXd( 12 ) <<
                                            1.0, 2.0, 3.0,  4.0,  5.0,  6.0,
                                            7.0, 8.0, 9.0, 10.0, 11.0, 12.0 ).finished( );
    SelectedAccelerationMap accelerations;
    accelerations[ "a" ][ "A" ] = { std::make_shared< AccelerationSettings >( point_mass_gravity ) };
    accelerations[ "a" ][ "B" ] = { std::make_shared< AccelerationSettings >( point_mass_gravity ) };
    accelerations[ "b" ][ "A" ] = { std::make_shared< AccelerationSettings >( point_mass_gravity ) };
    accelerations[ "b" ][ "B" ] = { std::make_shared< AccelerationSettings >( point_mass_gravity ) };
    const std::shared_ptr< SingleArcPropagatorSettings< double > > manualSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >( centralBodies,
                                                                                  accelerations,
                                                                                  bodiesToPropagate,
                                                                                  initialStates,
                                                                                  terminationSettings,
                                                                                  propagatorType );
    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 4: mass propagator
BOOST_AUTO_TEST_CASE( test_json_propagator_mass )
{
    using namespace propagators;
    using namespace simulation_setup;
    using namespace json_interface;

    // Create MassPropagatorSettings from JSON file
    const std::shared_ptr< SingleArcPropagatorSettings< double > > fromFileSettings =
            parseJSONFile< std::shared_ptr< SingleArcPropagatorSettings< double > > >( INPUT( "mass" ) );

    // Create MassPropagatorSettings manually
    const std::vector< std::string > bodiesToPropagate = { "a", "b" };
    const std::shared_ptr< PropagationTerminationSettings > terminationSettings;
    const Eigen::VectorXd initialStates = ( Eigen::VectorXd( 2 ) << 100.0, 200.0 ).finished( );
    SelectedMassRateModelMap massRateModels;
    massRateModels[ "a" ] = { std::make_shared< FromThrustMassModelSettings >( ) };
    massRateModels[ "b" ] = { std::make_shared< FromThrustMassModelSettings >( ) };
    const std::shared_ptr< SingleArcPropagatorSettings< double > > manualSettings =
            std::make_shared< MassPropagatorSettings< double > >( bodiesToPropagate,
                                                                    massRateModels,
                                                                    initialStates,
                                                                    terminationSettings );
    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

// Test 5: rotational propagator
BOOST_AUTO_TEST_CASE( test_json_propagator_rotational )
{
    using namespace basic_astrodynamics;
    using namespace propagators;
    using namespace simulation_setup;
    using namespace json_interface;

    // Create RotationalStatePropagatorSettings from JSON file
    const std::shared_ptr< SingleArcPropagatorSettings< double > > fromFileSettings =
            parseJSONFile< std::shared_ptr< SingleArcPropagatorSettings< double > > >( INPUT( "rotational" ) );

    // Create RotationalStatePropagatorSettings manually
    const std::vector< std::string > bodiesToPropagate = { "A", "B" };
    const std::shared_ptr< PropagationTerminationSettings > terminationSettings;
    const Eigen::VectorXd initialStates = ( Eigen::VectorXd( 14 ) <<
                                            0.0, 1.0, 2.0,  3.0,  4.0,  5.0,  6.0,
                                            7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0 ).finished( );
    SelectedTorqueMap torques;
    torques[ "A" ][ "B" ] = { std::make_shared< TorqueSettings >( second_order_gravitational_torque ) };
    torques[ "B" ][ "A" ] = { std::make_shared< TorqueSettings >( second_order_gravitational_torque ),
            std::make_shared< TorqueSettings >( aerodynamic_torque ) };
    const std::shared_ptr< SingleArcPropagatorSettings< double > > manualSettings =
            std::make_shared< RotationalStatePropagatorSettings< double > >( torques,
                                                                               bodiesToPropagate,
                                                                               initialStates,
                                                                               terminationSettings );
    // Compare
    BOOST_CHECK_EQUAL_JSON( fromFileSettings, manualSettings );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
