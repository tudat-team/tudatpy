/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <limits>
#include <string>
#include <vector>

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>
#include <boost/lambda/lambda.hpp>

#include "tudat/basics/testMacros.h"

#include "tudat/io/basicInputOutput.h"
#include "tudat/interface/spice/spiceInterface.h"

#include "tudat/simulation/estimation_setup/createObservationModel.h"
#include "tudat/astro/observation_models/oneWayRangeObservationModel.h"
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

BOOST_AUTO_TEST_SUITE( test_N_way_observation_partials)


std::vector< double > getRetransmissionDelays( const double evaluationTime, const int numberOfRetransmitters )
{
    std::vector< double > retransmissionDelays;

        for( int i = 0; i < numberOfRetransmitters; i++ )
        {
            retransmissionDelays.push_back( evaluationTime * 5.0E-17 * static_cast< double >( i + 1 ) );
        }
    return retransmissionDelays;
}

//! Test partial derivatives of one-way range observable, using general test suite of observation partials.
BOOST_AUTO_TEST_CASE( testnWayRangePartials )
{
    // Define and create ground stations.
    std::vector< std::pair< std::string, std::string > > groundStations;
    groundStations.resize( 2 );
    groundStations[ 0 ] = std::make_pair( "Earth", "Graz" );
    groundStations[ 1 ] = std::make_pair( "Mars", "MSL" );


    // Perform test for 2, 3 and 4-way
    for( unsigned int linkNumber = 0; linkNumber < 3; linkNumber++ )
    {
        // Set link ends for observation model
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = groundStations[ 0 ];
        linkEnds[ reflector1 ] = groundStations[ 1 ];
        if( linkNumber > 0 )
        {
            linkEnds[ reflector2 ] = groundStations[ 0 ];
        }
        if( linkNumber > 1 )
        {
            linkEnds[ reflector3 ] = groundStations[ 1 ];
        }

        if( linkNumber % 2 == 0 )
        {
            linkEnds[ receiver ] = groundStations[ 0 ];
        }
        else
        {
            linkEnds[ receiver ] = groundStations[ 1 ];
        }

        // Test partials with constant ephemerides (allows test of position partials)
        {
            // Create environment
            SystemOfBodies bodies = setupEnvironment( groundStations, 1.0E7, 1.2E7, 1.1E7, true );

            // Generate n-way range model
            std::vector< std::string > perturbingBodies;
            perturbingBodies.push_back( "Earth" );
            std::vector< std::shared_ptr< observation_models::ObservationModelSettings > > legObservationModels;
            for( unsigned int i = 0; i < linkNumber + 2; i ++ )
            {
                legObservationModels.push_back(
                            std::make_shared< observation_models::ObservationModelSettings >(
                                one_way_range,
                                getSingleLegLinkEnds( linkEnds, i ),
                                std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                                    perturbingBodies ) ) );
            }

            std::shared_ptr< ObservationModel< 1 > > nWayRangeModel =
                    observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                        std::make_shared< observation_models::NWayRangeObservationSettings >(
                            legObservationModels, std::bind( &getRetransmissionDelays, std::placeholders::_1, linkNumber + 1 ) ), bodies  );

            // Create parameter objects.
            std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet =
                    createEstimatableParameters( bodies, 1.1E7 );

            // Test observation partials
            testObservationPartials< 1 >(
                        nWayRangeModel, bodies, fullEstimatableParameterSet, linkEnds, n_way_range, 2.0E-6, true, true, 1.0,
                        ( Eigen::Vector4d( ) << 10.0, 1.0, 1.0, 10.0 ).finished( ) );
        }

        // Test partials with real ephemerides (without test of position partials)
        {
            // Create environment
            SystemOfBodies bodies = setupEnvironment( groundStations, 1.0E7, 1.2E7, 1.1E7, false );

            // Generate n-way range model
            std::vector< std::string > perturbingBodies;
            perturbingBodies.push_back( "Earth" );
            std::vector< std::shared_ptr< observation_models::ObservationModelSettings > > legObservationModels;
            for( unsigned int i = 0; i < linkNumber + 2; i ++ )
            {
                legObservationModels.push_back(
                            std::make_shared< observation_models::ObservationModelSettings >(
                                one_way_range,
                                getSingleLegLinkEnds( linkEnds, i ),
                                std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                                    perturbingBodies ) ) );
            }
            std::shared_ptr< ObservationModel< 1 > > nWayRangeModel =
                    observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                        std::make_shared< observation_models::NWayRangeObservationSettings >(
                            legObservationModels, std::bind( &getRetransmissionDelays, std::placeholders::_1, linkNumber + 1 ) ), bodies  );

            // Create parameter objects.
            std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet =
                    createEstimatableParameters( bodies, 1.1E7 );

            // Test observation partials
            testObservationPartials< 1 >(
                        nWayRangeModel, bodies, fullEstimatableParameterSet, linkEnds, n_way_range, 2.0E-6, false, true, 1.0,
                        ( Eigen::Vector4d( ) << 10.0, 1.0, 1.0, 20.0 ).finished( ) );
        }
    }
}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat




