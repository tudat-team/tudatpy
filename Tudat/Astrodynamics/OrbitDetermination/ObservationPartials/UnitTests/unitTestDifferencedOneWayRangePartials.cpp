/*    Copyright (c) 2010-2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *
 *    References
 *
 */

//#define BOOST_TEST_MAIN

#include <limits>
#include <string>
#include <vector>

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>
#include <boost/lambda/lambda.hpp>

#include "Tudat/Basics/testMacros.h"

#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/External/SpiceInterface/spiceInterface.h"

#include "Tudat/SimulationSetup/EstimationSetup/createObservationModel.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/constantRotationRate.h"
#include "Tudat/SimulationSetup/EstimationSetup/createObservationPartials.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/UnitTests/numericalObservationPartial.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createGroundStations.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.h"

#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/UnitTests/observationPartialTestFunctions.h"

//namespace tudat
//{
//namespace unit_tests
//{

using namespace tudat::gravitation;
using namespace tudat::ephemerides;
using namespace tudat::observation_models;
using namespace tudat::simulation_setup;
using namespace tudat::spice_interface;
using namespace tudat::observation_partials;
using namespace tudat::estimatable_parameters;
using namespace tudat::unit_tests;
using namespace tudat;

//BOOST_AUTO_TEST_SUITE( test_one_way_observation_partials)

//! Test partial derivatives of one-way differenced range observable, using general test suite of observation partials.
int main( )
//BOOST_AUTO_TEST_CASE( testOneWayRangePartials )
{

    Eigen::VectorXd parameterPerturbationMultipliers =
            ( Eigen::VectorXd( 4 ) << 100.0, 100.0, 10.0, 10.0 ).finished( );

    // Define and create ground stations.
    std::vector< std::pair< std::string, std::string > > groundStations;
    groundStations.resize( 2 );
    groundStations[ 0 ] = std::make_pair( "Earth", "Graz" );
    groundStations[ 1 ] = std::make_pair( "Mars", "MSL" );

    // Test partials with constant ephemerides (allows test of position partials)
    {
        // Create environment
        NamedBodyMap bodyMap = setupEnvironment( groundStations, 1.0E7, 1.2E7, 1.1E7, true );

        // Set link ends for observation model
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = groundStations[ 1 ];
        linkEnds[ receiver ] = groundStations[ 0 ];

        // Generate one-way differenced range model
        boost::shared_ptr< ObservationModel< 1 > > oneWayDifferencedRangeModel =
                observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                    linkEnds, boost::make_shared< observation_models::OneWayDifferencedRangeRateObservationSettings >(
                        boost::lambda::constant( 60.0 ) ), bodyMap  );

        // Create parameter objects.
        boost::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet =
                createEstimatableParameters( bodyMap, 1.1E7 );

        testObservationPartials< 1 >(
                    oneWayDifferencedRangeModel, bodyMap, fullEstimatableParameterSet, linkEnds,
                    one_way_differenced_range, 1.0E-4, true, true, 1000.0, parameterPerturbationMultipliers );
    }

    // Test partials with real ephemerides (without test of position partials)
    {
        // Create environment
        NamedBodyMap bodyMap = setupEnvironment( groundStations, 1.0E7, 1.2E7, 1.1E7, false );

        // Set link ends for observation model
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = groundStations[ 1 ];
        linkEnds[ receiver ] = groundStations[ 0 ];

        // Generate one-way range model
        boost::shared_ptr< ObservationModel< 1 > > oneWayDifferencedRangeModel =
                observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                    linkEnds, boost::make_shared< observation_models::OneWayDifferencedRangeRateObservationSettings >(
                        boost::lambda::constant( 60.0 ) ), bodyMap  );

        // Create parameter objects.
        boost::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet =
                createEstimatableParameters( bodyMap, 1.1E7 );

        testObservationPartials< 1 >(
                    oneWayDifferencedRangeModel, bodyMap, fullEstimatableParameterSet, linkEnds,
                    one_way_differenced_range, 1.0E-4, false, true, 1000.0, parameterPerturbationMultipliers );
    }
}


//BOOST_AUTO_TEST_SUITE_END( )

//} // namespace unit_tests

//} // namespace tudat




