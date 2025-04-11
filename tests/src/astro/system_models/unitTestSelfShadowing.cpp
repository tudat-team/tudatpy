
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

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_self_shadowing )

using namespace tudat::system_models;
using namespace tudat::simulation_setup;
using mathematical_constants::PI;


BOOST_AUTO_TEST_CASE( testFractionAnalytical )
{
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

    SystemOfBodies bodies;
    bodies.createEmptyBody( "L_SHAPED" );
    addBodyExteriorPanelledShape( panelSettings,
        "L_SHAPED",
        bodies );

    std::map< std::string, std::vector< std::shared_ptr< VehicleExteriorPanel > > > sortedBodyPanelMap =
        bodies.at( "L_SHAPED" )->getVehicleSystems( )->getVehicleExteriorPanels( );

    std::vector< std::shared_ptr< VehicleExteriorPanel > > bodyFixedPanels = sortedBodyPanelMap[ "" ];
    std::vector< double > angles = { PI/10, PI/6, PI/4 };
    int maximumNumberOfPixels = 1000;

    std::vector< int > indexesLit;
    std::vector< int > indexesShadowed;
    for ( unsigned int i = 0; i<bodyFixedPanels.size( ); i++)
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

    for ( unsigned int i = 0; i<angles.size( ); i++ )
    {
        Eigen::Vector3d incomingDirection;
        incomingDirection( 0 ) = -std::cos( angles[ i ] );
        incomingDirection( 1 ) = 0;
        incomingDirection( 2 ) = -std::sin( angles[ i ] );

        SelfShadowing selfShadowing( bodyFixedPanels, incomingDirection, maximumNumberOfPixels );
        std::vector< double > illuminatedPanelFractions = selfShadowing.getIlluminatedPanelFractions( );
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
