/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <limits>
#include <memory>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <boost/test/included/unit_test.hpp>
#include <Eigen/Core>

#include "tudat/io/serialization.h"

#include "tudat/basics/testMacros.h"
#include "tudat/simulation/propagation_setup/accelerationSettings.h"
#include "tudat/simulation/propagation_setup/propagationOutputSettings.h"
#include "tudat/simulation/propagation_setup/propagationTerminationSettings.h"
#include "tudat/simulation/propagation_setup/torqueSettings.h"
#include "tudat/simulation/propagation_setup/propagationSettings.h"
#include "tudat/simulation/estimation_setup/observationOutputSettings.h"
#include "tudat/astro/observation_models/observationAncillarySettings.h"
#include "tudat/math/root_finders/createRootFinder.h"

// Headers for new serialization tests
#include "tudat/simulation/estimation_setup/observationOutput.h"
#include "tudat/simulation/propagation_setup/propagationTermination.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/astro/observation_models/observableTypes.h"

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
    - FullTwoBodySphericalHarmonicAccelerationSettings
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
    std::shared_ptr< FullTwoBodySphericalHarmonicAccelerationSettings > fullTwoBodySphericalHarmonicSettings =
            std::make_shared< FullTwoBodySphericalHarmonicAccelerationSettings >( 3, 2, 4, 3, 2, 1 );
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
                                                                              fullTwoBodySphericalHarmonicSettings,
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
}

// Test VariableSettings / SaveSettings objects
BOOST_AUTO_TEST_CASE( test_VariableSettingsSerialization )
{
    using namespace propagators;
    using namespace reference_frames;
    using namespace gravitation;
    using namespace basic_astrodynamics;

    // Create objects
    std::shared_ptr< VariableSettings > varSettings = std::make_shared< VariableSettings >( cpuTimeVariable );
    std::shared_ptr< SingleDependentVariableSaveSettings > singleDepVarSettings =
            std::make_shared< SingleDependentVariableSaveSettings >( altitude_dependent_variable, "Earth" );
    std::shared_ptr< SingleAccelerationDependentVariableSaveSettings > singleAccSettings =
            std::make_shared< SingleAccelerationDependentVariableSaveSettings >( point_mass_gravity, "Earth", "Moon", false, -1 );
    std::shared_ptr< SphericalHarmonicAccelerationTermsDependentVariableSaveSettings > sphericalHarmonicTermsSettings =
            std::make_shared< SphericalHarmonicAccelerationTermsDependentVariableSaveSettings >( "Earth", "Moon", 2, 2, -1, false );
    std::shared_ptr< SingleTorqueDependentVariableSaveSettings > singleTorqueSettings =
            std::make_shared< SingleTorqueDependentVariableSaveSettings >(
                    basic_astrodynamics::aerodynamic_torque, "Earth", "Moon", false, -1 );
    std::shared_ptr< IntermediateAerodynamicRotationVariableSaveSettings > intermediateAeroRotSettings =
            std::make_shared< IntermediateAerodynamicRotationVariableSaveSettings >(
                    "Earth", corotating_frame, aerodynamic_frame, "Sun", -1 );
    std::shared_ptr< BodyAerodynamicAngleVariableSaveSettings > bodyAeroAngleSettings =
            std::make_shared< BodyAerodynamicAngleVariableSaveSettings >( "Earth", angle_of_attack, "Sun" );
    std::shared_ptr< LocalWindVelocityDependentVariableSaveSettings > localWindSettings =
            std::make_shared< LocalWindVelocityDependentVariableSaveSettings >( "Earth", "Earth", corotating_frame, -1 );
    std::shared_ptr< ControlSurfaceCoefficientDependentVariableSettings > controlSurfaceSettings =
            std::make_shared< ControlSurfaceCoefficientDependentVariableSettings >(
                    "Earth", "aileron", "Earth", aerodynamic_control_surface_force_coefficients_increment_dependent_variable );
    std::shared_ptr< SingleVariationSphericalHarmonicAccelerationSaveSettings > singleVarSettings =
            std::make_shared< SingleVariationSphericalHarmonicAccelerationSaveSettings >( "Earth", "Moon", basic_solid_body );
    std::shared_ptr< SingleVariationSingleTermSphericalHarmonicAccelerationSaveSettings > singleVarTermSettings =
            std::make_shared< SingleVariationSingleTermSphericalHarmonicAccelerationSaveSettings >(
                    "Earth", "Moon", 2, 2, basic_solid_body );
    std::shared_ptr< AccelerationPartialWrtStateSaveSettings > accPartialSettings =
            std::make_shared< AccelerationPartialWrtStateSaveSettings >( "Earth", "Moon", point_mass_gravity, "Earth" );
    std::shared_ptr< TotalAccelerationPartialWrtStateSaveSettings > totalAccPartialSettings =
            std::make_shared< TotalAccelerationPartialWrtStateSaveSettings >( "Earth", "Earth" );
    std::shared_ptr< MinimumConstellationDistanceDependentVariableSaveSettings > minConstDistSettings =
            std::make_shared< MinimumConstellationDistanceDependentVariableSaveSettings >(
                    "Earth", std::vector< std::string >( { "Moon", "Mars" } ) );
    std::shared_ptr< MinimumConstellationStationDistanceDependentVariableSaveSettings > minConstStationDistSettings =
            std::make_shared< MinimumConstellationStationDistanceDependentVariableSaveSettings >(
                    "Earth", "station1", std::vector< std::string >( { "Moon", "Mars" } ), 10.0 );
    std::shared_ptr< CrossSectionDependentVariableSaveSettings > crossSectionSettings =
            std::make_shared< CrossSectionDependentVariableSaveSettings >(
                    actual_cross_section, "Earth", "Sun", "cannon_ball_radiation_pressure" );
    std::shared_ptr< IlluminatedPanelFractionDependentVariableSaveSettings > illumPanelSettings =
            std::make_shared< IlluminatedPanelFractionDependentVariableSaveSettings >( "Earth", "Sun" );
    std::shared_ptr< TotalGravityFieldVariationSettings > totalGravFieldVarSettings =
            std::make_shared< TotalGravityFieldVariationSettings >( "Earth", 2, 5, 0, 5, true );
    std::shared_ptr< CustomDependentVariableSaveSettings > customDepVarSettings =
            std::make_shared< CustomDependentVariableSaveSettings >( []( ) { return Eigen::VectorXd::Zero( 3 ); }, 3 );

    // Collect in base-class pointer vector
    std::vector< std::shared_ptr< VariableSettings > > settingsVector = { varSettings,
                                                                          singleDepVarSettings,
                                                                          singleAccSettings,
                                                                          sphericalHarmonicTermsSettings,
                                                                          singleTorqueSettings,
                                                                          intermediateAeroRotSettings,
                                                                          bodyAeroAngleSettings,
                                                                          localWindSettings,
                                                                          controlSurfaceSettings,
                                                                          singleVarSettings,
                                                                          singleVarTermSettings,
                                                                          accPartialSettings,
                                                                          totalAccPartialSettings,
                                                                          minConstDistSettings,
                                                                          minConstStationDistSettings,
                                                                          crossSectionSettings,
                                                                          illumPanelSettings,
                                                                          totalGravFieldVarSettings,
                                                                          customDepVarSettings };

    // Round-trip test. CustomDependentVariableSaveSettings is expected to fail.
    for( const auto& settings : settingsVector )
    {
        std::stringstream ss;

        const bool isCustom = std::dynamic_pointer_cast< CustomDependentVariableSaveSettings >( settings ) != nullptr;

        try
        {
            {
                cereal::BinaryOutputArchive oarchive( ss );
                oarchive( settings );
            }

            std::shared_ptr< VariableSettings > deserializedSettings;

            {
                cereal::BinaryInputArchive iarchive( ss );
                iarchive( deserializedSettings );
            }

            BOOST_CHECK_MESSAGE( !isCustom, "CustomDependentVariableSaveSettings unexpectedly serialized/deserialized" );

            if( !isCustom )
            {
                BOOST_CHECK( *settings == *deserializedSettings );
            }
        }
        catch( ... )
        {
            BOOST_CHECK_MESSAGE( isCustom, "Serialization failed for non-custom dependent variable settings" );
        }
    }
}

// Test PropagationTerminationSettings objects
BOOST_AUTO_TEST_CASE( test_PropagationTerminationSettingsSerialization )
{
    using namespace propagators;

    // Create objects
    std::shared_ptr< PropagationTerminationSettings > timeTermSettings = std::make_shared< PropagationTimeTerminationSettings >( 1000.0 );
    std::shared_ptr< PropagationTerminationSettings > cpuTimeTermSettings =
            std::make_shared< PropagationCPUTimeTerminationSettings >( 60.0 );
    std::shared_ptr< PropagationTerminationSettings > depVarTermSettings =
            std::make_shared< PropagationDependentVariableTerminationSettings >(
                    std::make_shared< SingleDependentVariableSaveSettings >( altitude_dependent_variable, "Earth" ), 100.0, false );
    std::shared_ptr< PropagationTerminationSettings > hybridTermSettings = std::make_shared< PropagationHybridTerminationSettings >(
            std::vector< std::shared_ptr< PropagationTerminationSettings > >(
                    { std::make_shared< PropagationTimeTerminationSettings >( 500.0 ),
                      std::make_shared< PropagationCPUTimeTerminationSettings >( 30.0 ) } ),
            true );
    std::shared_ptr< PropagationTerminationSettings > nonSeqTermSettings = std::make_shared< NonSequentialPropagationTerminationSettings >(
            std::make_shared< PropagationTimeTerminationSettings >( 1000.0 ),
            std::make_shared< PropagationTimeTerminationSettings >( -1000.0 ) );
    std::shared_ptr< PropagationCustomTerminationSettings > customTermSettings =
            std::make_shared< PropagationCustomTerminationSettings >( []( const double ) { return true; } );
    std::shared_ptr< PropagationTerminationSettings > baseTermSettings =
            std::make_shared< PropagationTerminationSettings >( time_stopping_condition, false );

    const std::vector< std::pair< std::string, std::shared_ptr< PropagationTerminationSettings > > > settingsVector = {
        { "time", timeTermSettings },     { "CPU time", cpuTimeTermSettings },      { "dependent variable", depVarTermSettings },
        { "hybrid", hybridTermSettings }, { "non-sequential", nonSeqTermSettings }, { "custom", customTermSettings },
        { "base", baseTermSettings }
    };

    for( const auto& namedSettings : settingsVector )
    {
        const std::string& settingsName = namedSettings.first;
        const auto& settings = namedSettings.second;
        std::stringstream ss;

        const bool isCustom = std::dynamic_pointer_cast< PropagationCustomTerminationSettings >( settings ) != nullptr;
        std::string serializationPhase = "binary save";

        try
        {
            {
                cereal::BinaryOutputArchive oarchive( ss );
                oarchive( settings );
            }

            std::shared_ptr< PropagationTerminationSettings > deserializedSettings;

            serializationPhase = "binary load";
            {
                cereal::BinaryInputArchive iarchive( ss );
                iarchive( deserializedSettings );
            }

            BOOST_CHECK_MESSAGE( !isCustom, "PropagationCustomTerminationSettings unexpectedly serialized/deserialized" );

            if( !isCustom )
            {
                BOOST_REQUIRE_MESSAGE( deserializedSettings != nullptr,
                                       "Deserialization returned null for " << settingsName << " termination settings" );
                BOOST_CHECK_MESSAGE( *settings == *deserializedSettings,
                                     "Round-trip equality failed for " << settingsName << " termination settings" );
            }
        }
        catch( const std::exception& exception )
        {
            if( isCustom )
            {
                BOOST_CHECK_MESSAGE( std::string( exception.what( ) ).find( "checkStopCondition_" ) != std::string::npos,
                                     "Custom termination serialization threw the wrong exception: " << exception.what( ) );
            }
            else
            {
                BOOST_ERROR( "Serialization failed during " << serializationPhase << " for non-custom " << settingsName
                                                            << " termination settings: " << exception.what( ) );
            }
        }
        catch( ... )
        {
            BOOST_ERROR( "Serialization failed with a non-standard exception for " << settingsName << " termination settings" );
        }
    }
}

// Test TorqueSettings objects
BOOST_AUTO_TEST_CASE( test_TorqueSettingsSerialization )
{
    using namespace simulation_setup;
    using namespace basic_astrodynamics;

    std::shared_ptr< TorqueSettings > torqueSettings = std::make_shared< TorqueSettings >( aerodynamic_torque );
    std::shared_ptr< SphericalHarmonicTorqueSettings > sphericalHarmonicTorqueSettings =
            std::make_shared< SphericalHarmonicTorqueSettings >( 2, 2 );
    std::shared_ptr< FullTwoBodySphericalHarmonicTorqueSettings > fullTwoBodyTorqueSettings =
            std::make_shared< FullTwoBodySphericalHarmonicTorqueSettings >(
                    std::make_shared< FullTwoBodySphericalHarmonicAccelerationSettings >( 3, 2, 4, 3, 2, 1 ) );
    std::shared_ptr< FourthDegreeFullTwoBodyGravitationalTorqueSettings > fourthDegreeFullTwoBodyTorqueSettings =
            std::make_shared< FourthDegreeFullTwoBodyGravitationalTorqueSettings >( );
    std::shared_ptr< CustomTorqueSettings > customTorqueSettings =
            std::make_shared< CustomTorqueSettings >( []( const double ) { return Eigen::Vector3d( 1.0, 0.0, 0.0 ); } );

    std::vector< std::shared_ptr< TorqueSettings > > settingsVector = { torqueSettings,
                                                                        sphericalHarmonicTorqueSettings,
                                                                        fullTwoBodyTorqueSettings,
                                                                        fourthDegreeFullTwoBodyTorqueSettings,
                                                                        customTorqueSettings };

    for( const auto& settings : settingsVector )
    {
        std::stringstream ss;

        const bool isCustom = std::dynamic_pointer_cast< CustomTorqueSettings >( settings ) != nullptr;

        try
        {
            {
                cereal::BinaryOutputArchive oarchive( ss );
                oarchive( settings );
            }

            std::shared_ptr< TorqueSettings > deserializedSettings;

            {
                cereal::BinaryInputArchive iarchive( ss );
                iarchive( deserializedSettings );
            }

            BOOST_CHECK_MESSAGE( !isCustom, "CustomTorqueSettings unexpectedly serialized/deserialized" );

            if( !isCustom )
            {
                BOOST_CHECK( *settings == *deserializedSettings );
            }
        }
        catch( ... )
        {
            BOOST_CHECK_MESSAGE( isCustom, "Serialization failed for non-custom torque settings" );
        }
    }
}

// Test all PropagatorProcessingSettings derived types through the polymorphic base.
// The hybrid case specifically checks that consistentArcPrintSettings_ survives the round trip.
BOOST_AUTO_TEST_CASE( test_PropagatorProcessingSettingsSerialization )
{
    using namespace propagators;

    auto singleArcSettings = std::make_shared< SingleArcPropagatorProcessingSettings >(
            true, true, 4, 12.0, std::make_shared< PropagationPrintSettings >( true, false, 5.0, 3, true ), true );
    auto multiArcSettings = std::make_shared< MultiArcPropagatorProcessingSettings >(
            std::make_shared< PropagationPrintSettings >( false, true, 7.0, 2, false, true ), true, false, true, true, true );
    auto hybridArcSettings = std::make_shared< HybridArcPropagatorProcessingSettings >(
            std::make_shared< PropagationPrintSettings >( true, true, 9.0, 5, true, true, true ), true, false, true, true );

    // Concrete round trip isolates the hybrid print-setting field from polymorphic registration.
    const std::string serializedHybridArcSettings = serialization::serializeToBinaryString( *hybridArcSettings );
    const auto deserializedHybridArcSettings =
            serialization::deserializeFromBinaryString< HybridArcPropagatorProcessingSettings >( serializedHybridArcSettings );
    BOOST_CHECK( *hybridArcSettings == deserializedHybridArcSettings );

    std::vector< std::shared_ptr< PropagatorProcessingSettings > > settingsVector = { singleArcSettings,
                                                                                      multiArcSettings,
                                                                                      hybridArcSettings };

    for( const auto& settings : settingsVector )
    {
        std::stringstream stream;
        {
            cereal::BinaryOutputArchive archive( stream );
            archive( settings );
        }

        std::shared_ptr< PropagatorProcessingSettings > deserializedSettings;
        {
            cereal::BinaryInputArchive archive( stream );
            archive( deserializedSettings );
        }

        BOOST_REQUIRE( deserializedSettings != nullptr );
        BOOST_CHECK( *settings == *deserializedSettings );
    }
}

// Test ObservationDependentVariableSettings objects
BOOST_AUTO_TEST_CASE( test_ObservationDependentVariableSettingsSerialization )
{
    using namespace simulation_setup;
    using namespace observation_models;

    std::shared_ptr< ObservationDependentVariableSettings > obsDepVarSettings =
            std::make_shared< ObservationDependentVariableSettings >( station_elevation_angle );
    std::shared_ptr< StationAngleObservationDependentVariableSettings > stationAngleSettings =
            std::make_shared< StationAngleObservationDependentVariableSettings >( station_elevation_angle );
    std::shared_ptr< InterlinkObservationDependentVariableSettings > interlinkSettings =
            std::make_shared< InterlinkObservationDependentVariableSettings >( integration_time_dependent_variable );
    std::shared_ptr< AncillaryObservationDependentVariableSettings > ancillarySettings =
            std::make_shared< AncillaryObservationDependentVariableSettings >( integration_time_dependent_variable, one_way_range );
    std::shared_ptr< LightTimeCorrectionComponentsDependentVariableSettings > lightTimeSettings =
            std::make_shared< LightTimeCorrectionComponentsDependentVariableSettings >(
                    transmitter, receiver, LinkEndId( "", "" ), LinkEndId( "", "" ) );

    std::vector< std::shared_ptr< ObservationDependentVariableSettings > > settingsVector = {
        obsDepVarSettings, stationAngleSettings, interlinkSettings, ancillarySettings, lightTimeSettings
    };

    for( const auto& settings : settingsVector )
    {
        std::stringstream ss;

        try
        {
            {
                cereal::BinaryOutputArchive oarchive( ss );
                oarchive( settings );
            }

            std::shared_ptr< ObservationDependentVariableSettings > deserializedSettings;

            {
                cereal::BinaryInputArchive iarchive( ss );
                iarchive( deserializedSettings );
            }

            BOOST_CHECK( *settings == *deserializedSettings );
        }
        catch( std::exception& e )
        {
            BOOST_ERROR( "Serialization failed for observation dependent variable settings: " << e.what( ) );
        }
    }
}

// Test ObservationAncillarySimulationSettings
BOOST_AUTO_TEST_CASE( test_ObservationAncillarySimulationSettingsSerialization )
{
    using namespace observation_models;

    std::shared_ptr< ObservationAncillarySimulationSettings > settings = std::make_shared< ObservationAncillarySimulationSettings >( );

    std::stringstream ss;

    try
    {
        {
            cereal::BinaryOutputArchive oarchive( ss );
            oarchive( settings );
        }

        std::shared_ptr< ObservationAncillarySimulationSettings > deserializedSettings;

        {
            cereal::BinaryInputArchive iarchive( ss );
            iarchive( deserializedSettings );
        }

        BOOST_CHECK( *settings == *deserializedSettings );
    }
    catch( std::exception& e )
    {
        BOOST_ERROR( "Serialization failed for ObservationAncillarySimulationSettings: " << e.what( ) );
    }
}

// Test RootFinderSettings
BOOST_AUTO_TEST_CASE( test_RootFinderSettingsSerialization )
{
    using namespace root_finders;

    std::shared_ptr< RootFinderSettings > settings =
            std::make_shared< RootFinderSettings >( bisection_root_finder, 1.0e-6, 1.0e-6, 1.0e-8, 1000 );

    std::stringstream ss;

    try
    {
        {
            cereal::BinaryOutputArchive oarchive( ss );
            oarchive( settings );
        }

        std::shared_ptr< RootFinderSettings > deserializedSettings;

        {
            cereal::BinaryInputArchive iarchive( ss );
            iarchive( deserializedSettings );
        }

        BOOST_CHECK( *settings == *deserializedSettings );
    }
    catch( std::exception& e )
    {
        BOOST_ERROR( "Serialization failed for RootFinderSettings: " << e.what( ) );
    }
}

// Test ObservationDependentVariableBookkeeping serialization
BOOST_AUTO_TEST_CASE( test_ObservationDependentVariableBookkeepingSerialization )
{
    using namespace observation_models;
    using namespace simulation_setup;

    // Create a minimal bookkeeping with observable type and link ends
    LinkDefinition linkEnds;
    linkEnds[ transmitter ] = LinkEndId( "Earth", "station" );
    linkEnds[ receiver ] = LinkEndId( "Moon", "" );

    auto bookkeeping = std::make_shared< ObservationDependentVariableBookkeeping >( one_way_range, linkEnds );
    auto deferredSetting = std::make_shared< LightTimeCorrectionComponentsDependentVariableSettings >(
            transmitter, receiver, LinkEndId( "Earth", "station" ), LinkEndId( "Moon", "" ) );
    bookkeeping->addDeferredSetting( deferredSetting );

    std::stringstream ss;

    try
    {
        {
            cereal::BinaryOutputArchive oarchive( ss );
            oarchive( bookkeeping );
        }

        std::shared_ptr< ObservationDependentVariableBookkeeping > deserializedBookkeeping;

        {
            cereal::BinaryInputArchive iarchive( ss );
            iarchive( deserializedBookkeeping );
        }

        BOOST_REQUIRE( deserializedBookkeeping != nullptr );
        BOOST_REQUIRE_EQUAL( deserializedBookkeeping->getDeferredSettings( ).size( ), 1 );
        BOOST_CHECK( *bookkeeping == *deserializedBookkeeping );
        BOOST_CHECK( *deferredSetting == *deserializedBookkeeping->getDeferredSettings( ).front( ) );
        BOOST_CHECK_EQUAL( deserializedBookkeeping->getTotalDependentVariableSize( ), 0 );
    }
    catch( std::exception& e )
    {
        BOOST_ERROR( "Serialization failed for ObservationDependentVariableBookkeeping: " << e.what( ) );
    }
}

BOOST_AUTO_TEST_CASE( test_SerializationArchiveValidation )
{
    std::ostringstream provenanceWarnings;
    std::streambuf* originalErrorBuffer = std::cerr.rdbuf( provenanceWarnings.rdbuf( ) );
    serialization::checkProvenance( serialization::currentProvenance( ) );
    const std::string warningsForMatchingProvenance = provenanceWarnings.str( );
    serialization::checkProvenance( "test-version@test-commit@test-date" );
    std::cerr.rdbuf( originalErrorBuffer );

    BOOST_CHECK( warningsForMatchingProvenance.empty( ) );
    BOOST_CHECK_NE( provenanceWarnings.str( ).find( "[Tudat] Warning: loading serialized file" ), std::string::npos );

    const auto makeDimensionArchive = []( const Eigen::Index rows, const Eigen::Index cols ) {
        std::stringstream stream;
        {
            cereal::BinaryOutputArchive archive( stream );
            archive( cereal::make_nvp( "rows", rows ), cereal::make_nvp( "cols", cols ) );
        }
        return stream.str( );
    };

    const auto loadDynamicMatrix = []( const std::string& data ) {
        std::stringstream stream( data );
        cereal::BinaryInputArchive archive( stream );
        Eigen::MatrixXd matrix;
        archive( matrix );
    };

    BOOST_CHECK_THROW( loadDynamicMatrix( makeDimensionArchive( -1, 3 ) ), std::runtime_error );
    BOOST_CHECK_THROW(
            loadDynamicMatrix( makeDimensionArchive( 1, static_cast< Eigen::Index >( cereal::kMaximumSerializedEigenCoefficients + 1 ) ) ),
            std::runtime_error );
    BOOST_CHECK_THROW(
            loadDynamicMatrix( makeDimensionArchive( 0, static_cast< Eigen::Index >( cereal::kMaximumSerializedEigenCoefficients + 1 ) ) ),
            std::runtime_error );

    const std::string fixedSizeData = makeDimensionArchive( 2, 3 );
    BOOST_CHECK_THROW(
            [ &fixedSizeData ]( ) {
                std::stringstream stream( fixedSizeData );
                cereal::BinaryInputArchive archive( stream );
                Eigen::Matrix3d matrix;
                archive( matrix );
            }( ),
            std::runtime_error );

    std::stringstream oversizedProvenance( std::ios::in | std::ios::out | std::ios::binary );
    const std::uint32_t magic = serialization::kBinaryMagic;
    const std::uint32_t oversizedLength = serialization::kMaximumBinaryProvenanceLength + 1;
    oversizedProvenance.write( reinterpret_cast< const char* >( &magic ), sizeof( magic ) );
    oversizedProvenance.write( reinterpret_cast< const char* >( &oversizedLength ), sizeof( oversizedLength ) );
    oversizedProvenance.seekg( 0 );
    BOOST_CHECK_THROW( serialization::readBinaryProvenance( oversizedProvenance ), std::runtime_error );

    std::stringstream truncatedProvenance( std::ios::in | std::ios::out | std::ios::binary );
    const std::uint32_t truncatedLength = 8;
    truncatedProvenance.write( reinterpret_cast< const char* >( &magic ), sizeof( magic ) );
    truncatedProvenance.write( reinterpret_cast< const char* >( &truncatedLength ), sizeof( truncatedLength ) );
    truncatedProvenance.write( "short", 5 );
    truncatedProvenance.seekg( 0 );
    BOOST_CHECK_THROW( serialization::readBinaryProvenance( truncatedProvenance ), std::runtime_error );
}

// Test PropagationTerminationDetailsFromHybridCondition serialization
BOOST_AUTO_TEST_CASE( test_PropagationTerminationDetailsFromHybridConditionSerialization )
{
    using namespace propagators;

    // Create a hybrid termination condition to populate the details
    auto fixedTimeCondition = std::make_shared< FixedTimePropagationTerminationCondition >( 1000.0, true );
    std::vector< std::shared_ptr< PropagationTerminationCondition > > conditions = { fixedTimeCondition };
    auto hybridCondition = std::make_shared< HybridPropagationTerminationCondition >( conditions, true );

    auto details = std::make_shared< PropagationTerminationDetailsFromHybridCondition >( true, hybridCondition );

    // Serialize through the base type pointer to exercise polymorphic tagging
    std::shared_ptr< PropagationTerminationDetails > detailsBase = details;

    std::stringstream ss;

    try
    {
        {
            cereal::BinaryOutputArchive oarchive( ss );
            oarchive( detailsBase );
        }

        std::shared_ptr< PropagationTerminationDetails > deserializedDetails;

        {
            cereal::BinaryInputArchive iarchive( ss );
            iarchive( deserializedDetails );
        }

        BOOST_REQUIRE( deserializedDetails != nullptr );
        BOOST_CHECK( *detailsBase == *deserializedDetails );
    }
    catch( std::exception& e )
    {
        BOOST_ERROR( "Serialization failed for PropagationTerminationDetailsFromHybridCondition: " << e.what( ) );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
