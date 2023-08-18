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
#include <string>
#include <vector>

#include <boost/test/unit_test.hpp>

#include <boost/lambda/lambda.hpp>

#include "tudat/basics/testMacros.h"

#include "tudat/io/basicInputOutput.h"
#include "tudat/interface/spice/spiceInterface.h"

#include "tudat/astro/basic_astro/sphericalBodyShapeModel.h"
#include "tudat/astro/ephemerides/constantEphemeris.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/constantRotationRate.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/constantRotationalOrientation.h"
#include "tudat/simulation/estimation_setup/createCartesianStatePartials.h"
#include "tudat/support/numericalObservationPartial.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/simulation/environment_setup/createBodies.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/estimation_setup/createEstimatableParameters.h"
#include "tudat/support/observationPartialTestFunctions.h"


namespace tudat
{
namespace unit_tests
{

using namespace tudat::basic_astrodynamics;
using namespace tudat::ground_stations;
using namespace tudat::gravitation;
using namespace tudat::ephemerides;
using namespace tudat::observation_models;
using namespace tudat::simulation_setup;
using namespace tudat::spice_interface;
using namespace tudat::orbit_determination;
using namespace tudat::estimatable_parameters;
using namespace tudat::observation_partials;

BOOST_AUTO_TEST_SUITE( test_relative_position_partials)

//! Test partial derivatives of ideal relative positions w.r.t. parameters.
BOOST_AUTO_TEST_CASE( testRelativePositionStatePartials )
{
    // Define and create ground stations.
    std::vector< std::pair< std::string, std::string > > groundStations;
    groundStations.resize( 2 );
    groundStations[ 0 ] = std::make_pair( "Earth", "" );
    groundStations[ 1 ] = std::make_pair( "Mars", "" );

    // Create environment
    SystemOfBodies bodies = setupEnvironment( groundStations, 1.0E7, 1.2E7, 1.1E7, true );

    // Set link ends for observation model
    LinkEnds linkEnds;
    linkEnds[ observer ] = groundStations[ 0 ];
    linkEnds[ observed_body ] = groundStations[ 1 ];


    // Generate one-way range model
    std::shared_ptr< ObservationModel< 3 > > relativePositionModel = observation_models::ObservationModelCreator< 3, double, double >::createObservationModel(
            std::make_shared< observation_models::ObservationModelSettings >(
                    observation_models::relative_position_observable, linkEnds, std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ) ), bodies );

    // Create parameter objects.
    std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > > initialStateParameters;
    initialStateParameters.push_back( std::make_shared< InitialTranslationalStateParameter< double > >(
            "Earth", propagators::getInitialStateOfBody( "Earth", "SSB", bodies, 1.1E7 ) ) );
    initialStateParameters.push_back( std::make_shared< InitialTranslationalStateParameter< double > >(
            "Mars", propagators::getInitialStateOfBody( "Mars", "SSB", bodies, 1.1E7 ) ) );

    std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate = std::make_shared< EstimatableParameterSet< double > >(
            std::vector< std::shared_ptr< EstimatableParameter< double > > >( ), std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > >( ),
                    initialStateParameters );
    printEstimatableParameterEntries( parametersToEstimate );

    // Create observation partials.
    std::pair< std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< 3 > > >,
            std::shared_ptr< PositionPartialScaling > > analyticalPartialSet =
            ObservationPartialCreator< 3, double, double >::createObservationPartials(
                    relativePositionModel, bodies, parametersToEstimate );

    std::shared_ptr< PositionPartialScaling > positionPartialScaler = analyticalPartialSet.second;

    // Iterate over link ends, compute and test partials for observable referenced at each link end.
    for ( LinkEnds::const_iterator linkEndIterator = linkEnds.begin( ) ; linkEndIterator != linkEnds.end( ) ; linkEndIterator++ )
    {
            // Evaluate nominal observation values
            std::vector<Eigen::Vector6d> vectorOfStates;
            std::vector<double> vectorOfTimes;
            double observationTime = 1.1E7;
            Eigen::VectorXd currentObservation = relativePositionModel->computeObservationsWithLinkEndData(
                    observationTime, linkEndIterator->first, vectorOfTimes, vectorOfStates, nullptr );

            // Calculate analytical observation partials.
            if (positionPartialScaler != NULL) {
                positionPartialScaler->update(vectorOfStates, vectorOfTimes, static_cast< LinkEndType >( linkEndIterator->first ),
                                              currentObservation);
            }


        std::vector< std::vector< std::pair< Eigen::Matrix< double, 3, Eigen::Dynamic >, double > > > partialList;
        for( typename std::map< std::pair< int, int >,
                std::shared_ptr< ObservationPartial< 3 > > >::const_iterator partialIterator =
                analyticalPartialSet.first.begin( ); partialIterator != analyticalPartialSet.first.end( ); partialIterator++ )
        {
            partialList.push_back( partialIterator->second->calculatePartial( vectorOfStates, vectorOfTimes, linkEndIterator->first, currentObservation ) );
        }

        Eigen::MatrixXd expectedPartialsWrtObserverState = Eigen::MatrixXd::Zero( 3, 6 );
        expectedPartialsWrtObserverState.block( 0, 0, 3, 3 ) = - Eigen::Matrix3d::Identity( );

        Eigen::MatrixXd expectedPartialsWrtObservedState = Eigen::MatrixXd::Zero( 3, 6 );
        expectedPartialsWrtObservedState.block( 0, 0, 3, 3 ) = Eigen::Matrix3d::Identity( );


        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                partialList.at( 0 ).at( 0 ).first, expectedPartialsWrtObserverState, std::numeric_limits< double >::epsilon( ) );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                partialList.at( 1 ).at( 0 ).first, expectedPartialsWrtObservedState, std::numeric_limits< double >::epsilon( ) );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat




