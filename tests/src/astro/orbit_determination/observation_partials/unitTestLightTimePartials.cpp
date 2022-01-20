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

#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/astro/observation_models/lightTimeSolution.h"
#include "tudat/simulation/estimation_setup/createObservationModel.h"
#include "tudat/simulation/estimation_setup/createEstimatableParameters.h"
#include "tudat/simulation/estimation_setup/createLightTimeCorrectionPartials.h"
#include "tudat/astro/orbit_determination/observation_partials/firstOrderRelativisticPartial.h"
#include "tudat/support/observationPartialTestFunctions.h"

namespace tudat
{
namespace unit_tests
{

//Using declarations.
using namespace tudat::gravitation;
using namespace tudat::ephemerides;
using namespace tudat::observation_models;
using namespace tudat::simulation_setup;
using namespace tudat::spice_interface;
using namespace tudat::observation_partials;
using namespace tudat::orbit_determination;
using namespace tudat::observation_models;
using namespace tudat::estimatable_parameters;
using namespace tudat::unit_tests;

BOOST_AUTO_TEST_SUITE( test_light_time_partials )

BOOST_AUTO_TEST_CASE( testOneWayRangePartialsWrtLightTimeParameters )
{
    {
        // Define and create ground stations.
        std::vector< std::pair< std::string, std::string > > groundStations;
        groundStations.resize( 2 );
        groundStations[ 0 ] = std::make_pair( "Earth", "Graz" );
        groundStations[ 1 ] = std::make_pair( "Mars", "MSL" );

        // Create environment
        SystemOfBodies bodies = setupEnvironment( groundStations );

        double ephemerisEvaluationTime = basic_astrodynamics::calculateJulianDaySinceEpoch< double >(
                    boost::gregorian::date( 2002, 8, 10 ), 0.0 ) * physical_constants::JULIAN_DAY;

        std::dynamic_pointer_cast< ConstantEphemeris >( bodies.at( "Earth" )->getEphemeris( ) )->updateConstantState(
                    getBodyCartesianStateAtEpoch( "Earth", "SSB", "ECLIPJ2000", "NONE", ephemerisEvaluationTime ) );
        std::dynamic_pointer_cast< ConstantEphemeris >( bodies.at( "Mars" )->getEphemeris( ) )->updateConstantState(
                    getBodyCartesianStateAtEpoch( "Mars", "SSB", "ECLIPJ2000", "NONE", ephemerisEvaluationTime ) );


        // Set link ends for observation model
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = groundStations[ 1 ];
        linkEnds[ receiver ] = groundStations[ 0 ];

        // Generate one-way range model
        std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrections;
        std::vector< std::string > relativisticPerturbingBodies = { "Sun" };
        lightTimeCorrections.push_back( std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                                            relativisticPerturbingBodies ) );
        std::shared_ptr< ObservationModelSettings > observationSettings = std::make_shared<
                ObservationModelSettings >( one_way_range, linkEnds, lightTimeCorrections );
        std::shared_ptr< OneWayRangeObservationModel< double, double > > oneWayRangeModel =
                std::dynamic_pointer_cast< OneWayRangeObservationModel< double, double > >(
                    observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                        observationSettings, bodies ) );


        // Create parameters for which partials are to be computed
        std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > > parameterSettings;
        parameterSettings.push_back( std::make_shared< estimatable_parameters::EstimatableParameterSettings >(
                                         "Sun", gravitational_parameter ) );
        parameterSettings.push_back( std::make_shared< estimatable_parameters::EstimatableParameterSettings >(
                                         "global_metric", ppn_parameter_gamma ) );
        std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
                createParametersToEstimate( parameterSettings, bodies );

        // Create partial objects.
        std::pair< SingleLinkObservationPartialList, std::shared_ptr< PositionPartialScaling > > partialList =
                createOneWayRangePartials( linkEnds, bodies, parametersToEstimate,
                                           oneWayRangeModel->getLightTimeCalculator( )->getLightTimeCorrection( ) );
        std::shared_ptr< PositionPartialScaling > positionPartialScaler = partialList.second;

        // Compute current observation and link end times/states
        std::vector< double > linkEndTimes;
        std::vector< Eigen::Matrix< double, 6, 1 > > linkEndStates;
        double testTime = 1.1E7;
        Eigen::VectorXd currentPositionObservation =
                oneWayRangeModel->computeObservationsWithLinkEndData( testTime, transmitter, linkEndTimes, linkEndStates );

        // Update position partial scaler for current observation
        positionPartialScaler->update( linkEndStates, linkEndTimes, transmitter, currentPositionObservation );


        // Define numerical partial settings
        std::vector< double > perturbations = { 1.0E16, 1.0E8 };
        std::vector< double > tolerances = { 1.0E-4, 10E-4 };

        // Compute numerical partials for each parameter and compare to analytical result.
        std::function< double( const double ) > observationFunction = std::bind(
                    &ObservationModel< 1, double, double >::computeObservationEntry, oneWayRangeModel, std::placeholders::_1, transmitter, 0 );
        for( SingleLinkObservationPartialList::iterator partialIterator = partialList.first.begin( ); partialIterator != partialList.first.end( );
             partialIterator++ )
        {
            // Compute total analytical partial
            std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > calculatedPartial =
                    partialIterator->second->calculatePartial( linkEndStates, linkEndTimes, transmitter );
            Eigen::Matrix< double, 1, Eigen::Dynamic > totalPartial = Eigen::Matrix< double, 1, Eigen::Dynamic >::Zero( 1, 1 );
            for( unsigned int j = 0; j < calculatedPartial.size( ); j++ )
            {
                totalPartial += calculatedPartial.at( j ).first;
            }

            // Compute numerical partial and compare to analytical result.
            BOOST_CHECK_CLOSE_FRACTION(
                        calculateNumericalObservationParameterPartial(
                            parametersToEstimate->getDoubleParameters( )[ partialIterator->first.first ],
                        perturbations.at( partialIterator->first.first ), observationFunction, testTime ).x( ), totalPartial.x( ),
                    tolerances.at( partialIterator->first.first ) );
        }
    }

    {

        // Define and create ground stations.
        std::vector< std::pair< std::string, std::string > > groundStations;
        groundStations.resize( 2 );
        groundStations[ 0 ] = std::make_pair( "Earth", "Graz" );
        groundStations[ 1 ] = std::make_pair( "Mars", "MSL" );

        // Create environment
        SystemOfBodies bodies = setupEnvironment( groundStations, 1.0E7, 1.2E7, 1.65E7 );


        // Set link ends for observation model
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = groundStations[ 1 ];
        linkEnds[ receiver ] = groundStations[ 0 ];

        // Generate one-way range model
        std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrections;
        std::vector< std::string > perturbingBodyList = { "Earth", "Sun" };
        lightTimeCorrections.push_back( std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                                            perturbingBodyList ) );
        std::shared_ptr< ObservationModelSettings > observationSettings = std::make_shared<
                ObservationModelSettings >( one_way_range, linkEnds, lightTimeCorrections );
        std::shared_ptr< ObservationModel< 1 > > oneWayRangeModel =
                observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                    observationSettings, bodies  );

        std::map< LinkEnds, std::shared_ptr< ObservationModel< 1 > > > oneWayRangeModelMap;
        oneWayRangeModelMap[ linkEnds ] = oneWayRangeModel;

        // Create parameter objects.
        std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Sun", gravitational_parameter ) );
        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Earth", gravitational_parameter ) );
        parameterNames.push_back( std::make_shared< estimatable_parameters::EstimatableParameterSettings >(
                                      "global_metric", ppn_parameter_gamma ) );
        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Mars", gravitational_parameter ) );
        std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
                createParametersToEstimate< double >( parameterNames, bodies );
        std::vector< std::shared_ptr< EstimatableParameter< double > > > doubleParameterVector =
                parametersToEstimate->getEstimatedDoubleParameters( );

        // Create observation partials.
        std::shared_ptr< ObservationPartialCreator< 1, double, double > > observationPartialCreator =
                std::make_shared< ObservationPartialCreator< 1, double, double > >( );
        std::pair< std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< 1 > > >,
                std::shared_ptr< PositionPartialScaling > > fullAnalyticalPartialSet =
                observationPartialCreator->createObservationPartials(
                    one_way_range, oneWayRangeModelMap, bodies, parametersToEstimate ).begin( )->second;

        std::shared_ptr< PositionPartialScaling > positionPartialScaler = fullAnalyticalPartialSet.second;

        // Compute partials for each refernce link end.
        for( LinkEnds::const_iterator linkEndIterator = linkEnds.begin( ); linkEndIterator != linkEnds.end( );
             linkEndIterator++ )
        {
            // Evaluate nominal observation values
            std::vector< Eigen::Vector6d > vectorOfStates;
            std::vector< double > vectorOfTimes;
            double observationTime = 1.1E7;
            Eigen::VectorXd currentRangeObservation = oneWayRangeModel->computeObservationsWithLinkEndData(
                        observationTime, linkEndIterator->first, vectorOfTimes, vectorOfStates );

            // Calculate analytical observation partials.
            positionPartialScaler->update( vectorOfStates, vectorOfTimes, static_cast< LinkEndType >( linkEndIterator->first ),
                                           currentRangeObservation );
            typedef std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > ObservationPartialReturnType;
            std::vector< ObservationPartialReturnType > analyticalObservationPartials =
                    calculateAnalyticalPartials(
                        fullAnalyticalPartialSet.first, vectorOfStates, vectorOfTimes, linkEndIterator->first );

            // Test evaliuuation time of partials
            for( unsigned int i = 0; i < analyticalObservationPartials.size( ); i++ )
            {
                for( unsigned int j = 0; j < analyticalObservationPartials.at( i ).size( ); j++ )
                {
                    BOOST_CHECK_CLOSE_FRACTION( analyticalObservationPartials.at( i ).at( j ).second,
                                                ( vectorOfTimes.at( 0 ) + vectorOfTimes.at( 1 ) ) / 2.0,
                                                std::numeric_limits< double >::epsilon( ) );
                }

            }

            // Settings for body state partials
            std::function< Eigen::VectorXd( const double ) > observationFunction = std::bind(
                        &ObservationModel< 1, double, double >::computeObservations, oneWayRangeModel, std::placeholders::_1,
                        linkEndIterator->first );

            // Settings for parameter partial functions.
            std::vector< double > parameterPerturbations = { 1.0E19, 1.0E16, 1.0E15, 1.0E8 };
            std::vector< std::function< void( ) > > updateFunctionList;
            updateFunctionList.push_back( emptyVoidFunction );
            updateFunctionList.push_back( emptyVoidFunction );
            updateFunctionList.push_back( emptyVoidFunction );
            updateFunctionList.push_back( emptyVoidFunction );

            // Calculate and test analytical against numerical partials.
            std::vector< Eigen::VectorXd > numericalPartialsWrtDoubleParameters = calculateNumericalPartialsWrtDoubleParameters(
                        doubleParameterVector, updateFunctionList, parameterPerturbations, observationFunction, observationTime );

            // Compare analytical and numerical partials
            for( unsigned int i = 0; i < analyticalObservationPartials.size( ); i++ )
            {

                double currentParameterPartial = 0.0;
                for( unsigned int j = 0; j < analyticalObservationPartials.at( i ).size( ); j++ )
                {
                    currentParameterPartial += analyticalObservationPartials.at( i ).at( j ).first.x( );

                }

                BOOST_CHECK_CLOSE_FRACTION( currentParameterPartial, numericalPartialsWrtDoubleParameters.at( i ).x( ), 1.0E-4 );
            }

            BOOST_CHECK_EQUAL( numericalPartialsWrtDoubleParameters[ 3 ].x( ), 0.0 );

        }
    }
}
BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
