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

#include <iostream>
#include <ctime>

#include <boost/format.hpp>
#include <boost/test/unit_test.hpp>

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/environment_setup/createSystemModel.h"
#include "tudat/astro/system_models/vehicleExteriorPanels.h"
#include "tudat/astro/system_models/selfShadowing.h"
#include "tudat/astro/ephemerides/constantRotationalEphemeris.h"
#include "tudat/simulation/environment_setup/createBodies.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"

namespace tudat
{

using namespace tudat::basic_astrodynamics;
using namespace tudat::simulation_setup;
using namespace tudat::ephemerides;
using namespace tudat::electromagnetism;
using namespace tudat::system_models;
using mathematical_constants::PI;

using namespace tudat::propagators;
using namespace tudat::numerical_integrators;
using namespace tudat::orbital_element_conversions;
using namespace tudat::basic_mathematics;
using namespace tudat::gravitation;
using namespace tudat::estimatable_parameters;

namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_self_shadowing )

BOOST_AUTO_TEST_CASE( testFractionAnalytical )
{
    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Create bodies needed in simulation
    double initialEphemerisTime = 0.0;
    double finalEphemerisTime = 1.1 * 365.25 * 86400.0;
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Sun" );
    SystemOfBodies bodies = createSystemOfBodies( getDefaultBodySettings( bodyNames, initialEphemerisTime, finalEphemerisTime ) );

    std::map< std::string, std::shared_ptr< MaterialProperties > > materialPropertiesMap;
    materialPropertiesMap[ "dummy" ] = materialProperties( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 );
    materialPropertiesMap[ "TO_BE_SHADOWED" ] = materialProperties( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 );
    materialPropertiesMap[ "TO_BE_LIT" ] = materialProperties( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 );

    std::map< std::string, bool > instantaneousReradiation;
    instantaneousReradiation[ "dummy" ] = true;
    instantaneousReradiation[ "TO_BE_SHADOWED" ] = true;
    instantaneousReradiation[ "TO_BE_LIT" ] = true;

    std::vector< std::shared_ptr< BodyPanelSettings > > bodyPanelSettingList =
            bodyPanelSettingsListFromDae( tudat::paths::getTudatTestDataPath( ) + "selfShadowingUnitTest.dae",
                                          Eigen::Vector3d::Zero( ),
                                          materialPropertiesMap,
                                          instantaneousReradiation );

    std::shared_ptr< FullPanelledBodySettings > panelSettings = fullPanelledBodySettings( bodyPanelSettingList );

    bodies.createEmptyBody( "L_SHAPED" );
    // Define constant rotational ephemeris
    Eigen::Vector7d rotationalStateVehicle;
    rotationalStateVehicle.segment( 0, 4 ) =
            linear_algebra::convertQuaternionToVectorFormat( Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ) );
    rotationalStateVehicle.segment( 4, 3 ) = Eigen::Vector3d::Zero( );
    bodies.at( "L_SHAPED" )
            ->setRotationalEphemeris( std::make_shared< ConstantRotationalEphemeris >( rotationalStateVehicle, "ECLIPJ2000" ) );
    addBodyExteriorPanelledShape( panelSettings, "L_SHAPED", bodies );

    std::map< std::string, std::vector< std::shared_ptr< VehicleExteriorPanel > > > sortedBodyPanelMap =
            bodies.at( "L_SHAPED" )->getVehicleSystems( )->getVehicleExteriorPanels( );

    std::vector< std::shared_ptr< VehicleExteriorPanel > > bodyFixedPanels = sortedBodyPanelMap.at( "" );
    std::map< std::string, std::vector< std::string > > sourceToTargetOccultingBodies =
            std::map< std::string, std::vector< std::string > >( );
    const std::map< std::string, int > maximumNumberOfPixelsPerSource = { { "Sun", 1000 } };
    for( const auto& pair: panelSettings->partRotationModelSettings_ )
    {
        std::cout << pair.first << std::endl;
    }

    PaneledRadiationPressureTargetModel targetModel(
            bodyFixedPanels,
            bodyFixedPanels,
            std::map< std::string, std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > > >( ),
            std::map< std::string, std::function< Eigen::Quaterniond( ) > >( ),
            sourceToTargetOccultingBodies,
            maximumNumberOfPixelsPerSource,
            true );

    std::vector< double > angles = { PI / 50, PI / 20, PI / 10, PI / 8, PI / 6, PI / 4 };

    std::vector< int > indexesLit;
    std::vector< int > indexesShadowed;
    for( unsigned int i = 0; i < 20; i++ )
    {
        if( bodyFixedPanels.at( i )->getPanelTypeId( ) == "TO_BE_SHADOWED" )
        {
            indexesShadowed.push_back( i );
        }
        if( bodyFixedPanels.at( i )->getPanelTypeId( ) == "TO_BE_LIT" )
        {
            indexesLit.push_back( i );
        }
    }

    std::map< std::string, std::shared_ptr< SelfShadowing > > mapSSH = targetModel.getSelfShadowingPerSources( );
    Eigen::Vector3d incomingDirection = Eigen::Vector3d::UnitX( );
    bodies.at( "L_SHAPED" )->getVehicleSystems( )->updatePartOrientations( 0.0 );

    for( unsigned int i = 0; i < angles.size( ); i++ )
    {
        incomingDirection( 0 ) = -std::sin( angles[ i ] );
        incomingDirection( 1 ) = 0;
        incomingDirection( 2 ) = -std::cos( angles[ i ] );

        mapSSH[ "Sun" ]->reset( );
        mapSSH[ "Sun" ]->updateIlluminatedPanelFractions( incomingDirection );
        std::vector< double > illuminatedPanelFractions = mapSSH[ "Sun" ]->getIlluminatedPanelFractions( );
        double trueFractionShaded = 1.0 - std::tan( angles[ i ] );
        double trueFractionLit = 1.0;

        double actualFractionShadowed =
                0.5 * ( illuminatedPanelFractions[ indexesShadowed[ 0 ] ] + illuminatedPanelFractions[ indexesShadowed[ 1 ] ] );
        double actualFractionLit = 0.25 *
                ( illuminatedPanelFractions[ indexesLit[ 0 ] ] + illuminatedPanelFractions[ indexesLit[ 1 ] ] +
                  illuminatedPanelFractions[ indexesLit[ 2 ] ] + illuminatedPanelFractions[ indexesLit[ 3 ] ] );

        BOOST_CHECK( actualFractionLit == trueFractionLit );
        BOOST_CHECK( std::abs( actualFractionShadowed - trueFractionShaded ) < 1e-3 );
    }
}

BOOST_AUTO_TEST_CASE( testComputationalEfficiency )
{
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation time settings.
    const double simulationStartEpoch = 0.0;
    const double simulationEndEpoch = 0.5 * tudat::physical_constants::JULIAN_DAY;

    // Define body settings for simulation.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Moon" );

    // Create body objects.
    BodyListSettings bodySettings = getDefaultBodySettings( bodiesToCreate, "SSB", "J2000" );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::map< std::string, std::shared_ptr< MaterialProperties > > materialPropertiesMap;
    materialPropertiesMap[ "dummy" ] = materialProperties( 0.4, 0.4, 0.0, 0.0, 0.0, 0.0 );
    materialPropertiesMap[ "TO_BE_SHADOWED" ] = materialProperties( 0.4, 0.4, 0.0, 0.0, 0.0, 0.0 );
    materialPropertiesMap[ "TO_BE_LIT" ] = materialProperties( 0.4, 0.4, 0.0, 0.0, 0.0, 0.0 );

    std::map< std::string, bool > instantaneousReradiation;
    instantaneousReradiation[ "dummy" ] = true;
    instantaneousReradiation[ "TO_BE_SHADOWED" ] = true;
    instantaneousReradiation[ "TO_BE_LIT" ] = true;

    std::vector< std::shared_ptr< BodyPanelSettings > > bodyPanelSettingList =
            bodyPanelSettingsListFromDae( tudat::paths::getTudatTestDataPath( ) + "selfShadowingUnitTest.dae",
                                          Eigen::Vector3d::Zero( ),
                                          materialPropertiesMap,
                                          instantaneousReradiation );

    std::shared_ptr< FullPanelledBodySettings > panelSettings = fullPanelledBodySettings( bodyPanelSettingList );
    // Create spacecraft object.
    bodies.createEmptyBody( "Vehicle" );
    bodies.at( "Vehicle" )->setConstantBodyMass( 400.0 );
    // Define constant rotational ephemeris
    Eigen::Vector7d rotationalStateVehicle;
    rotationalStateVehicle.segment( 0, 4 ) =
            linear_algebra::convertQuaternionToVectorFormat( Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ) );
    rotationalStateVehicle.segment( 4, 3 ) = Eigen::Vector3d::Zero( );
    bodies.at( "Vehicle" )
            ->setRotationalEphemeris( std::make_shared< ConstantRotationalEphemeris >( rotationalStateVehicle, "ECLIPJ2000" ) );
    addBodyExteriorPanelledShape( panelSettings, "Vehicle", bodies );

    std::map< std::string, std::vector< std::string > > sourceToTargetOccultingBodies =
            std::map< std::string, std::vector< std::string > >( );
    const std::map< std::string, int > maximumNumberOfPixelsPerSource = { { "Sun", 100 } };

    bodies.at( "Vehicle" )
            ->setRadiationPressureTargetModels(
                    { createRadiationPressureTargetModel( std::make_shared< PaneledRadiationPressureTargetModelSettings >(
                                                                  sourceToTargetOccultingBodies, maximumNumberOfPixelsPerSource ),
                                                          "Vehicle",
                                                          bodies ) } );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Define propagation settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfAsterix;
    accelerationsOfAsterix[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );
    accelerationsOfAsterix[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );
    accelerationsOfAsterix[ "Moon" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );
    accelerationsOfAsterix[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::radiation_pressure ) );

    accelerationMap[ "Vehicle" ] = accelerationsOfAsterix;
    bodiesToPropagate.push_back( "Vehicle" );
    centralBodies.push_back( "Earth" );

    basic_astrodynamics::AccelerationMap accelerationModelMap =
            createAccelerationModelsMap( bodies, accelerationMap, bodiesToPropagate, centralBodies );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set Keplerian elements for Asterix.
    Eigen::Vector6d asterixInitialStateInKeplerianElements;
    asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7500.0E3;
    asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.1;
    asterixInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
    asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = unit_conversions::convertDegreesToRadians( 235.7 );
    asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = unit_conversions::convertDegreesToRadians( 23.4 );
    asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );

    double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
    const Eigen::Vector6d asterixInitialState =
            convertKeplerianToCartesianElements( asterixInitialStateInKeplerianElements, earthGravitationalParameter );

    const double fixedStepSize = 10.0;
    std::shared_ptr< IntegratorSettings<> > integratorSettings =
            rungeKuttaFixedStepSettings< double >( fixedStepSize, numerical_integrators::rungeKuttaFehlberg78 );

    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >( centralBodies,
                                                                                accelerationModelMap,
                                                                                bodiesToPropagate,
                                                                                asterixInitialState,
                                                                                simulationStartEpoch,
                                                                                integratorSettings,
                                                                                propagationTimeTerminationSettings( simulationEndEpoch ),
                                                                                cowell );
    SingleArcDynamicsSimulator<> dynamicsSimulator( bodies, propagatorSettings );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
