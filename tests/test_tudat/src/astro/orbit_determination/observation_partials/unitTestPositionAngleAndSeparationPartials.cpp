/*    Copyright (c) 2010-2024, Delft University of Technology
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
#include <string>
#include <vector>

#include <boost/test/included/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/simulation/estimation_setup/createObservationModelFactory.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/constantRotationRate.h"
#include "tudat/simulation/estimation_setup/createObservationPartials.h"
#include "tudat/support/numericalObservationPartial.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
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
    // Define and create ground stations.
    std::vector< std::pair< std::string, std::string > > groundStations;
    groundStations.resize( 3 );
    groundStations[ 0 ] = std::make_pair( "Earth", "Graz" );
    groundStations[ 1 ] = std::make_pair( "Mars", "MSL" );
    groundStations[ 2 ] = std::make_pair( "Moon", "" );

    std::vector< std::string > perturbingBodies;
    perturbingBodies.push_back( "Earth" );

    // Test partials with constant ephemerides (allows test of position partials)
    {
        // Create environment
        SystemOfBodies bodies = setupEnvironment( groundStations, 1.0E7, 1.2E7, 1.1E7, true );

        // Set link ends for observation model
        LinkDefinition linkEnds;
        linkEnds[ receiver ] = groundStations[ 0 ];
        linkEnds[ transmitter ] = groundStations[ 1 ];
        linkEnds[ transmitter2 ] = groundStations[ 2 ];

        // Test position angle partials
        {
            std::shared_ptr< ObservationModel< 1 > > positionAngleModel =
                    observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                            std::make_shared< observation_models::PositionAngleObservationModelSettings >(
                                    linkEnds,
                                    std::vector< std::shared_ptr< LightTimeCorrectionSettings > >{
                                            std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( perturbingBodies ) } ),
                            bodies );

            std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet = createEstimatableParameters( bodies, 1.1E7 );

            testObservationPartials< 1 >(
                    positionAngleModel, bodies, fullEstimatableParameterSet, linkEnds, position_angle, 1.0E-4, true, true );
        }

        // Test separation distance partials
        {
            std::shared_ptr< ObservationModel< 1 > > separationModel =
                    observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                            std::make_shared< observation_models::SeparationObservationModelSettings >(
                                    linkEnds,
                                    std::vector< std::shared_ptr< LightTimeCorrectionSettings > >{
                                            std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( perturbingBodies ) } ),
                            bodies );

            std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet = createEstimatableParameters( bodies, 1.1E7 );

            testObservationPartials< 1 >(
                    separationModel, bodies, fullEstimatableParameterSet, linkEnds, separation_distance, 1.0E-4, true, true );
        }

        // Test combined position angle and separation distance partials
        {
            std::shared_ptr< ObservationModel< 2 > > positionAngleAndSeparationModel =
                    observation_models::ObservationModelCreator< 2, double, double >::createObservationModel(
                            std::make_shared< observation_models::PositionAngleAndSeparationObservationModelSettings >(
                                    linkEnds,
                                    std::vector< std::shared_ptr< LightTimeCorrectionSettings > >{
                                            std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( perturbingBodies ) } ),
                            bodies );

            std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet = createEstimatableParameters( bodies, 1.1E7 );

            testObservationPartials< 2 >( positionAngleAndSeparationModel,
                                          bodies,
                                          fullEstimatableParameterSet,
                                          linkEnds,
                                          position_angle_and_separation,
                                          1.0E-4,
                                          true,
                                          true );
        }
    }

    // Test partials with real ephemerides (without test of position partials)
    {
        // Create environment
        SystemOfBodies bodies = setupEnvironment( groundStations, 1.0E7, 1.2E7, 1.1E7, false );

        // Set link ends for observation model
        LinkDefinition linkEnds;
        linkEnds[ receiver ] = groundStations[ 0 ];
        linkEnds[ transmitter ] = groundStations[ 1 ];
        linkEnds[ transmitter2 ] = groundStations[ 2 ];

        // Test position angle partials
        {
            std::shared_ptr< ObservationModel< 1 > > positionAngleModel =
                    observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                            std::make_shared< observation_models::PositionAngleObservationModelSettings >(
                                    linkEnds,
                                    std::vector< std::shared_ptr< LightTimeCorrectionSettings > >{
                                            std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( perturbingBodies ) } ),
                            bodies );

            std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet = createEstimatableParameters( bodies, 1.1E7 );

            testObservationPartials< 1 >(
                    positionAngleModel, bodies, fullEstimatableParameterSet, linkEnds, position_angle, 1.0E-4, false, true );
        }

        // Test separation distance partials
        {
            std::shared_ptr< ObservationModel< 1 > > separationModel =
                    observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                            std::make_shared< observation_models::SeparationObservationModelSettings >(
                                    linkEnds,
                                    std::vector< std::shared_ptr< LightTimeCorrectionSettings > >{
                                            std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( perturbingBodies ) } ),
                            bodies );

            std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet = createEstimatableParameters( bodies, 1.1E7 );

            testObservationPartials< 1 >(
                    separationModel, bodies, fullEstimatableParameterSet, linkEnds, separation_distance, 1.0E-4, false, true );
        }

        // Test combined position angle and separation distance partials
        {
            std::shared_ptr< ObservationModel< 2 > > positionAngleAndSeparationModel =
                    observation_models::ObservationModelCreator< 2, double, double >::createObservationModel(
                            std::make_shared< observation_models::PositionAngleAndSeparationObservationModelSettings >(
                                    linkEnds,
                                    std::vector< std::shared_ptr< LightTimeCorrectionSettings > >{
                                            std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( perturbingBodies ) } ),
                            bodies );

            std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet = createEstimatableParameters( bodies, 1.1E7 );

            testObservationPartials< 2 >( positionAngleAndSeparationModel,
                                          bodies,
                                          fullEstimatableParameterSet,
                                          linkEnds,
                                          position_angle_and_separation,
                                          1.0E-4,
                                          false,
                                          true );
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
