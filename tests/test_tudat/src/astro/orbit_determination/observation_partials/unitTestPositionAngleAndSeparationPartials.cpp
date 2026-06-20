/*    Copyright (c) 2010-2024, Delft University of Technology
 *    All rights reserved
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
#include <string>
#include <vector>

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/simulation/estimation_setup/createObservationModelFactory.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/constantRotationRate.h"
#include "tudat/simulation/estimation_setup/createObservationPartials.h"
#include "tudat/support/numericalObservationPartial.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/createEphemeris.h"
#include "tudat/support/observationPartialTestFunctions.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::gravitation;
using namespace tudat::ephemerides;
using namespace tudat::observation_models;
using namespace tudat::simulation_setup;
using namespace tudat::spice_interface;
using namespace tudat::observation_partials;
using namespace tudat::estimatable_parameters;

BOOST_AUTO_TEST_SUITE( test_position_angle_and_separation_partials )

BOOST_AUTO_TEST_CASE( testPositionAngleAndSeparationPartials )
{
    std::vector< std::pair< std::string, std::string > > groundStations;
    groundStations.resize( 3 );
    groundStations[ 0 ] = std::make_pair( "Earth", "Graz" );
    groundStations[ 1 ] = std::make_pair( "Mars", "MSL" );
    groundStations[ 2 ] = std::make_pair( "Moon", "" );

    SystemOfBodies bodies = setupEnvironment( groundStations, 1.0E7, 1.2E7, 1.1E7, true );

    LinkDefinition linkEnds;
    linkEnds[ receiver ] = groundStations[ 0 ];
    linkEnds[ transmitter ] = groundStations[ 1 ];
    linkEnds[ transmitter2 ] = groundStations[ 2 ];

    std::vector< std::string > perturbingBodies;
    perturbingBodies.push_back( "Earth" );

    // Test position angle partials
    {
        std::shared_ptr< ObservationModel< 1 > > paModel = ObservationModelCreator< 1, double, double >::createObservationModel(
                std::make_shared< PositionAngleObservationModelSettings >(
                        linkEnds,
                        std::vector< std::shared_ptr< LightTimeCorrectionSettings > >{
                                std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( perturbingBodies ) } ),
                bodies );

        std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet = createEstimatableParameters( bodies, 1.1E7 );

        auto paPartials =
                ObservationPartialCreator< 1, double, double >::createObservationPartials( paModel, bodies, fullEstimatableParameterSet );

        BOOST_CHECK( paPartials.first.size( ) > 0 );
        BOOST_CHECK( paPartials.second != nullptr );

        double observationTime = 1.1E7;
        std::vector< double > linkEndTimes;
        std::vector< Eigen::Vector6d > linkEndStates;
        Eigen::VectorXd observation = paModel->computeObservationsWithLinkEndData( observationTime, receiver, linkEndTimes, linkEndStates );

        paPartials.second->update( linkEndStates, linkEndTimes, receiver, observation );

        for( auto it : paPartials.first )
        {
            auto partialValues = it.second->calculatePartial( linkEndStates, linkEndTimes, receiver );
            BOOST_CHECK( partialValues.size( ) > 0 );
        }
    }

    // Test separation partials
    {
        std::shared_ptr< ObservationModel< 1 > > sepModel = ObservationModelCreator< 1, double, double >::createObservationModel(
                std::make_shared< SeparationObservationModelSettings >(
                        linkEnds,
                        std::vector< std::shared_ptr< LightTimeCorrectionSettings > >{
                                std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( perturbingBodies ) } ),
                bodies );

        std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet = createEstimatableParameters( bodies, 1.1E7 );

        auto sepPartials =
                ObservationPartialCreator< 1, double, double >::createObservationPartials( sepModel, bodies, fullEstimatableParameterSet );

        BOOST_CHECK( sepPartials.first.size( ) > 0 );
        BOOST_CHECK( sepPartials.second != nullptr );

        double observationTime = 1.1E7;
        std::vector< double > linkEndTimes;
        std::vector< Eigen::Vector6d > linkEndStates;
        Eigen::VectorXd observation =
                sepModel->computeObservationsWithLinkEndData( observationTime, receiver, linkEndTimes, linkEndStates );

        sepPartials.second->update( linkEndStates, linkEndTimes, receiver, observation );

        for( auto it : sepPartials.first )
        {
            auto partialValues = it.second->calculatePartial( linkEndStates, linkEndTimes, receiver );
            BOOST_CHECK( partialValues.size( ) > 0 );
        }
    }

    // Test combined position angle and separation partials
    {
        std::shared_ptr< ObservationModel< 2 > > pasModel = ObservationModelCreator< 2, double, double >::createObservationModel(
                std::make_shared< PositionAngleAndSeparationObservationModelSettings >(
                        linkEnds,
                        std::vector< std::shared_ptr< LightTimeCorrectionSettings > >{
                                std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( perturbingBodies ) } ),
                bodies );

        std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet = createEstimatableParameters( bodies, 1.1E7 );

        auto pasPartials =
                ObservationPartialCreator< 2, double, double >::createObservationPartials( pasModel, bodies, fullEstimatableParameterSet );

        BOOST_CHECK( pasPartials.first.size( ) > 0 );
        BOOST_CHECK( pasPartials.second != nullptr );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
