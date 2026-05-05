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

#include <limits>
#include <memory>
#include <sstream>
#include <vector>

#include <boost/test/unit_test.hpp>
#include <Eigen/Core>

#include "tudat/io/serialization.h"

#include "tudat/basics/testMacros.h"
#include "tudat/simulation/propagation_setup/accelerationSettings.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_settings_serialization )

// Test Acceleration settings objects
BOOST_AUTO_TEST_CASE( test_AccelerationSettingsSerialization )
{
    using namespace simulation_setup;
    /*
    Tests:
    - AccelerationSettings
    - RadiationPressureAccelerationSettings
    - SphericalHarmonicAccelerationSettings
    - MutualSphericalHarmonicAccelerationSettings
    - RelativisticAccelerationCorrectionSettings
    - EmpiricalAccelerationSettings
    - YarkovskyAccelerationSettings
    - ThrustAccelerationSettings
    - CustomAccelerationSettings (Should fail)
    - RTGAccelerationSettings
    - DirectTidalDissipationAccelerationSettings
    - MomentumWheelDesaturationAccelerationSettings
    */

    // Create objects
    std::shared_ptr< AccelerationSettings > accelerationSettings =
            std::make_shared< AccelerationSettings >( basic_astrodynamics::spherical_harmonic_gravity );
    std::shared_ptr< RadiationPressureAccelerationSettings > radiationPressureSettings =
            std::make_shared< RadiationPressureAccelerationSettings >( cannonball_target );
    std::shared_ptr< SphericalHarmonicAccelerationSettings > sphericalHarmonicSettings =
            std::make_shared< SphericalHarmonicAccelerationSettings >( 5, 5 );
    std::shared_ptr< MutualSphericalHarmonicAccelerationSettings > mutualSphericalHarmonicSettings =
            std::make_shared< MutualSphericalHarmonicAccelerationSettings >( 5, 5, 10, 10 );
    std::shared_ptr< RelativisticAccelerationCorrectionSettings > relativisticSettings =
            std::make_shared< RelativisticAccelerationCorrectionSettings >( true, false, false );
    std::shared_ptr< EmpiricalAccelerationSettings > empiricalSettings = std::make_shared< EmpiricalAccelerationSettings >(
            Eigen::Vector3d( 1.0, 2.0, 3.0 ), Eigen::Vector3d( 0.1, 0.2, 0.3 ), Eigen::Vector3d( 0.01, 0.02, 0.03 ) );
    std::shared_ptr< YarkovskyAccelerationSettings > yarkovskySettings = std::make_shared< YarkovskyAccelerationSettings >( 0.1 );
    std::shared_ptr< ThrustAccelerationSettings > thrustSettings =
            std::make_shared< ThrustAccelerationSettings >( std::vector< std::string >( { "engine1", "engine2" } ) );
    std::shared_ptr< CustomAccelerationSettings > customSettings =
            std::make_shared< CustomAccelerationSettings >( []( const double ) { return Eigen::Vector3d( 1.0, 2.0, 3.0 ); } );
    std::shared_ptr< RTGAccelerationSettings > rtgSettings =
            std::make_shared< RTGAccelerationSettings >( Eigen::Vector3d( 0.1, 0.0, 0.0 ), 0.1, 0.0 );
    std::shared_ptr< DirectTidalDissipationAccelerationSettings > directTidalSettings =
            std::make_shared< DirectTidalDissipationAccelerationSettings >( 0.1, 0.01, true, false, true );
    std::shared_ptr< MomentumWheelDesaturationAccelerationSettings > momentumWheelSettings =
            std::make_shared< MomentumWheelDesaturationAccelerationSettings >(
                    std::vector< double >( { 100.0, 200.0 } ),
                    std::vector< Eigen::Vector3d >( { Eigen::Vector3d( 1.0, 2.0, 3.0 ), Eigen::Vector3d( 0.1, 0.2, 0.3 ) } ),
                    500.0,
                    50.0 );

    // Create vector of base class pointers to objects so we can loop through
    std::vector< std::shared_ptr< AccelerationSettings > > settingsVector = { accelerationSettings,
                                                                              radiationPressureSettings,
                                                                              sphericalHarmonicSettings,
                                                                              mutualSphericalHarmonicSettings,
                                                                              relativisticSettings,
                                                                              empiricalSettings,
                                                                              yarkovskySettings,
                                                                              thrustSettings,
                                                                                customSettings,
                                                                              rtgSettings,
                                                                              directTidalSettings,
                                                                              momentumWheelSettings };

    // Serialize, de-serialize and check equality using public members. CustomAccelerationSettings is expected to fail.
    for( const auto& settings : settingsVector )
    {
        std::stringstream ss;

        const bool isCustom = std::dynamic_pointer_cast< CustomAccelerationSettings >( settings ) != nullptr;

        try
        {
            {
                cereal::BinaryOutputArchive oarchive( ss );
                oarchive( settings );
            }

            std::shared_ptr< AccelerationSettings > deserializedSettings;

            {
                cereal::BinaryInputArchive iarchive( ss );
                iarchive( deserializedSettings );
            }

            BOOST_CHECK_MESSAGE( !isCustom, "CustomAccelerationSettings unexpectedly serialized/deserialized" );

            if( !isCustom )
            {
                BOOST_CHECK( *settings == *deserializedSettings );
                if( !( *settings == *deserializedSettings ) )
                {
                    std::cerr << "FAILED TYPE: " << static_cast< int >( settings->accelerationType_ ) << std::endl;
                }
            }
        }
        catch( ... )
        {
            BOOST_CHECK_MESSAGE( isCustom, "Serialization failed for non-custom settings" );
        }
    }
}  // namespace

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat