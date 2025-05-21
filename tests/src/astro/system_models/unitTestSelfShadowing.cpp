
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

    std::map< std::string, std::vector< double > > materialProperties;
    materialProperties[ "dummy" ] = { 0.0, 0.0 };
    materialProperties[ "TO_BE_SHADOWED" ] = { 0.0, 0.0 };
    materialProperties[ "TO_BE_LIT" ] = { 0.0, 0.0 };

    std::map< std::string, bool > instantaneousReradiation;
    instantaneousReradiation[ "dummy" ] = true;
    instantaneousReradiation[ "TO_BE_SHADOWED" ] = true;
    instantaneousReradiation[ "TO_BE_LIT" ] = true;

    std::vector< std::shared_ptr< BodyPanelSettings > > bodyPanelSettingList = bodyPanelSettingsListFromDae( 
        tudat::paths::getTudatTestDataPath( ) + "selfShadowingUnitTest.dae",
        Eigen::Vector3d::Zero( ),
        materialProperties,
        instantaneousReradiation
    );

    std::shared_ptr< FullPanelledBodySettings > panelSettings = fullPanelledBodySettings( bodyPanelSettingList );

    bodies.createEmptyBody( "L_SHAPED" );
    // Define constant rotational ephemeris
    Eigen::Vector7d rotationalStateVehicle;
    rotationalStateVehicle.segment( 0, 4 ) =
            linear_algebra::convertQuaternionToVectorFormat( Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ) );
    rotationalStateVehicle.segment( 4, 3 ) = Eigen::Vector3d::Zero( );
    bodies.at( "L_SHAPED" )
            ->setRotationalEphemeris(
                    std::make_shared< ConstantRotationalEphemeris >( rotationalStateVehicle, "ECLIPJ2000" ) );
    addBodyExteriorPanelledShape( panelSettings,
        "L_SHAPED",
        bodies );
    
    std::map< std::string, std::vector< std::shared_ptr< VehicleExteriorPanel > > > sortedBodyPanelMap =
        bodies.at( "L_SHAPED" )->getVehicleSystems( )->getVehicleExteriorPanels( );

    std::vector< std::shared_ptr< VehicleExteriorPanel > > bodyFixedPanels = sortedBodyPanelMap.at( "" );
    std::map< std::string, std::vector< std::string > > sourceToTargetOccultingBodies = std::map< std::string, std::vector< std::string > >( );
    const std::map< std::string, int > maximumNumberOfPixelsPerSource = {{"Sun", 1000}};

    PaneledRadiationPressureTargetModel targetModel( bodyFixedPanels, 
            std::map< std::string, std::vector< std::shared_ptr< system_models::VehicleExteriorPanel > > >( ),
            std::map< std::string, std::function< Eigen::Quaterniond( ) > >( ),
            sourceToTargetOccultingBodies, 
            maximumNumberOfPixelsPerSource );
    
    std::vector< double > angles = { PI/50, PI/20, PI/10, PI/8, PI/6, PI/4 };

    std::vector< int > indexesLit;
    std::vector< int > indexesShadowed;
    for ( unsigned int i = 0; i<20; i++)
    {
        if ( bodyFixedPanels.at( i )->getPanelTypeId( ) == "TO_BE_SHADOWED" )
        {
            indexesShadowed.push_back( i );
        }
        if ( bodyFixedPanels.at( i )->getPanelTypeId( ) == "TO_BE_LIT" )
        {
            indexesLit.push_back( i );
        }
    }
    
    std::map< std::string, std::shared_ptr< SelfShadowing > > mapSSH = targetModel.getSelfShadowingPerSources( );
    Eigen::Vector3d incomingDirection = Eigen::Vector3d::UnitX( );
    mapSSH[ "Sun" ]->updateIlluminatedPanelFractions( incomingDirection );

    for ( unsigned int i = 0; i<angles.size( ); i++ )
    {
        incomingDirection( 0 ) = -std::sin( angles[ i ] );
        incomingDirection( 1 ) = 0;
        incomingDirection( 2 ) = -std::cos( angles[ i ] );

        mapSSH[ "Sun" ]->updateIlluminatedPanelFractions( incomingDirection );
        std::vector< double > illuminatedPanelFractions = mapSSH[ "Sun" ]->getIlluminatedPanelFractions( );
        double trueFractionShaded = 1.0 - std::tan( angles[ i ] );
        double trueFractionLit = 1.0;

        double actualFractionShadowed = 0.5*( illuminatedPanelFractions[ indexesShadowed[ 0 ] ] + 
            illuminatedPanelFractions[ indexesShadowed[ 1 ] ]);
        double actualFractionLit = 0.25*( illuminatedPanelFractions[ indexesLit[ 0 ] ] +
            illuminatedPanelFractions[ indexesLit[ 1 ] ] +
            illuminatedPanelFractions[ indexesLit[ 2 ] ] +
            illuminatedPanelFractions[ indexesLit[ 3 ] ] );
                    
        BOOST_CHECK( actualFractionLit == trueFractionLit );
        BOOST_CHECK( std::abs( actualFractionShadowed - trueFractionShaded ) < 1e-3);
    }
   
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
