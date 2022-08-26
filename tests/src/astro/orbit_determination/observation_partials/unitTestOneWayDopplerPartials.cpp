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
#include "tudat/astro/observation_models/oneWayDopplerObservationModel.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/constantRotationRate.h"
#include "tudat/simulation/estimation_setup/createObservationPartials.h"
#include "tudat/support/numericalObservationPartial.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/math/basic/numericalDerivative.h"
#include "tudat/support/observationPartialTestFunctions.h"

namespace tudat
{
namespace unit_tests
{


BOOST_AUTO_TEST_SUITE( test_one_way_observation_partials)

Eigen::Vector3d computeUnitVectorToReceiverFromTransmitterState(
        const Eigen::Vector3d receiverPosition,
        const std::function< Eigen::Vector6d( const double ) > transmitterStateFunction,
        const double evaluationTime )
{
    return ( receiverPosition - transmitterStateFunction( evaluationTime ).segment( 0, 3 ) ).normalized( );
}


Eigen::Vector3d computeUnitVectorToReceiverFromReceiverState(
        const std::function< Eigen::Vector6d( const double ) > receiverStateFunction,
        const Eigen::Vector3d transmitterPosition,
        const double evaluationTime )
{
    return ( receiverStateFunction( evaluationTime ).segment( 0, 3 ) - transmitterPosition ).normalized( );
}

Eigen::VectorXd getProperTimeRateInVectorForm(
        std::shared_ptr< DopplerProperTimeRateInterface > properTimeRateCalculator,
        const std::vector< double >& linkEndTimes,
        const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
        const LinkEndType linkEndAssociatedWithTime )
{
    return ( Eigen::Vector1d( ) << properTimeRateCalculator->getOberverProperTimeDeviation(
                 linkEndTimes, linkEndStates ) ).finished( );
}

//! Test partial derivatives of one-way doppler observable, using general test suite of observation partials.
BOOST_AUTO_TEST_CASE( testOneWayDopplerPartials )
{

    using namespace tudat::gravitation;
    using namespace tudat::gravitation;
    using namespace tudat::ephemerides;
    using namespace tudat::observation_models;
    using namespace tudat::simulation_setup;
    using namespace tudat::spice_interface;
    using namespace tudat::observation_partials;
    using namespace tudat::estimatable_parameters;

    // Define and create ground stations.
    std::vector< std::pair< std::string, std::string > > groundStations;
    groundStations.resize( 2 );
    groundStations[ 0 ] = std::make_pair( "Earth", "Graz" );
    groundStations[ 1 ] = std::make_pair( "Mars", "MSL" );


    // Test ancilliary functions
    {
        double nominalEvaluationTime = 1.1E7;

        // Create environment
        SystemOfBodies bodies = setupEnvironment( groundStations, 1.0E7, 1.2E7, 1.1E7, false );

        // Set link ends for observation model
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = groundStations[ 1 ];
        linkEnds[ receiver ] = groundStations[ 0 ];

        // Create transmitter/receriver state functions
        std::function< Eigen::Vector6d( const double ) > transmitterStateFunction =
                getLinkEndCompleteEphemerisFunction< double, double >( linkEnds[ transmitter ], bodies );
        std::function< Eigen::Vector6d( const double ) > receiverStateFunction =
                getLinkEndCompleteEphemerisFunction< double, double >( linkEnds[ receiver ], bodies );

        // Define (independent!) transmission/reception times
        double transmissionTime = nominalEvaluationTime;
        double receptionTime = nominalEvaluationTime + 1.0E3;

        // Compute associated states
        Eigen::Vector6d nominalTransmitterState = transmitterStateFunction( transmissionTime );
        Eigen::Vector6d nominalReceiverState = receiverStateFunction( receptionTime );
        Eigen::Vector3d nominalVectorToReceiver = ( nominalReceiverState - nominalTransmitterState ).segment( 0, 3 );

        double timePerturbation = 100.0;

        // Partials for fixed receiver
        {
            // Compute numerical derivative of transmitter state for acceleration)
            Eigen::Vector6d numericalStateDerivative = numerical_derivatives::computeCentralDifferenceFromFunction(
                        transmitterStateFunction, transmissionTime, timePerturbation, numerical_derivatives::order8 );

            // Compute unit vector derivative numerically
            std::function< Eigen::Vector3d( const double ) > unitVectorFunction =
                    std::bind( &computeUnitVectorToReceiverFromTransmitterState,
                                 nominalReceiverState.segment( 0, 3 ), transmitterStateFunction, std::placeholders::_1 );
            Eigen::Vector3d numericalUnitVectorDerivative = numerical_derivatives::computeCentralDifferenceFromFunction(
                        unitVectorFunction, transmissionTime, timePerturbation, numerical_derivatives::order8 );

            // Compute projected velocoty vector derivative numerically
            std::function< double( const double) > projectedVelocityFunction =
                    std::bind( &calculateLineOfSightVelocityAsCFractionFromTransmitterStateFunction< double, double >,
                                 nominalReceiverState.segment( 0, 3 ), transmitterStateFunction, std::placeholders::_1 );
            double numericalProjectedVelocityDerivative =
                    numerical_derivatives::computeCentralDifferenceFromFunction(
                        projectedVelocityFunction, transmissionTime, timePerturbation, numerical_derivatives::order8 );

            // Compute analytical partial derivatives
            Eigen::Vector3d analyticalUnitVectorDerivative =
                    -computePartialOfUnitVectorWrtLinkEndTime(
                        nominalVectorToReceiver, nominalVectorToReceiver.normalized( ),
                        nominalVectorToReceiver.norm( ), nominalTransmitterState.segment( 3, 3 ) );
            double analyticalProjectedVelocityDerivative = computePartialOfProjectedLinkEndVelocityWrtAssociatedTime(
                        nominalVectorToReceiver,
                        nominalTransmitterState.segment( 3, 3 ),
                        nominalTransmitterState.segment( 3, 3 ),
                        numericalStateDerivative.segment( 3, 3 ), false );


            for( unsigned int i = 0; i < 3; i++ )
            {
                BOOST_CHECK_SMALL( std::fabs( analyticalUnitVectorDerivative( i ) - numericalUnitVectorDerivative( i ) ), 1.0E-16 );

            }
            BOOST_CHECK_SMALL( std::fabs( analyticalProjectedVelocityDerivative / physical_constants::SPEED_OF_LIGHT -
                                          numericalProjectedVelocityDerivative ), 1.0E-21 );
        }


        // Partials for fixed transmitter
        {
            // Compute numerical derivative of receiver state for acceleration)
            Eigen::Vector6d numericalStateDerivative = numerical_derivatives::computeCentralDifferenceFromFunction(
                        receiverStateFunction, receptionTime, timePerturbation, numerical_derivatives::order8 );

            // Compute unit vector derivative numerically
            std::function< Eigen::Vector3d( const double ) > unitVectorFunction =
                    std::bind( &computeUnitVectorToReceiverFromReceiverState,
                                 receiverStateFunction, nominalTransmitterState.segment( 0, 3 ), std::placeholders::_1 );
            Eigen::Vector3d numericalUnitVectorDerivative = numerical_derivatives::computeCentralDifferenceFromFunction(
                        unitVectorFunction, receptionTime, timePerturbation, numerical_derivatives::order8 );

            // Compute projected velocoty vector derivative numerically
            std::function< double( const double) > projectedVelocityFunction =
                    std::bind( &calculateLineOfSightVelocityAsCFractionFromReceiverStateFunction< double, double >,
                                 receiverStateFunction, nominalTransmitterState.segment( 0, 3 ), std::placeholders::_1 );
            double numericalProjectedVelocityDerivative =
                    numerical_derivatives::computeCentralDifferenceFromFunction(
                        projectedVelocityFunction, receptionTime, timePerturbation, numerical_derivatives::order8 );

            // Compute analytical partial derivatives
            Eigen::Vector3d analyticalUnitVectorDerivative =
                    computePartialOfUnitVectorWrtLinkEndTime(
                        nominalVectorToReceiver, nominalVectorToReceiver.normalized( ),
                        nominalVectorToReceiver.norm( ), nominalReceiverState.segment( 3, 3 ) );
            double analyticalProjectedVelocityDerivative = computePartialOfProjectedLinkEndVelocityWrtAssociatedTime(
                        nominalVectorToReceiver,
                        nominalReceiverState.segment( 3, 3 ),
                        nominalReceiverState.segment( 3, 3 ),\
                        numericalStateDerivative.segment( 3, 3 ), true );

            for( unsigned int i = 0; i < 3; i++ )
            {
                BOOST_CHECK_SMALL( std::fabs( analyticalUnitVectorDerivative( i ) - numericalUnitVectorDerivative( i ) ), 1.0E-17 );

            }
            BOOST_CHECK_SMALL( std::fabs( analyticalProjectedVelocityDerivative / physical_constants::SPEED_OF_LIGHT -
                                          numericalProjectedVelocityDerivative ), 1.5E-22 );
        }

    }

    // Test partials with constant ephemerides (allows test of position partials)
    {
        // Create environment
        SystemOfBodies bodies = setupEnvironment( groundStations, 1.0E7, 1.2E7, 1.1E7, true );

        // Set link ends for observation model
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = groundStations[ 1 ];
        linkEnds[ receiver ] = groundStations[ 0 ];

        for( unsigned int estimationCase  = 0; estimationCase  < 3; estimationCase ++ )
        {
            std::cout << "Case " << estimationCase << std::endl;
            // Generate one-way doppler model
            std::shared_ptr< ObservationModel< 1 > > oneWayDopplerModel;
            std::vector< std::string > perturbingBodies;
            perturbingBodies.push_back( "Earth" );
            if( estimationCase  == 0 )
            {
                oneWayDopplerModel =
                        observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                            std::make_shared< observation_models::ObservationModelSettings >(
                                observation_models::one_way_doppler, linkEnds,
                                std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                                    perturbingBodies ) ), bodies  );
            }
            else
            {
                oneWayDopplerModel =
                        observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                            std::make_shared< OneWayDopplerObservationSettings >
                            (  linkEnds, std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                                   perturbingBodies ),
                               std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Mars" ),
                               std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Earth" ) ), bodies  );
            }

            // Create parameter objects.
            std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet;
            Eigen::VectorXd parameterPerturbationMultipliers = Eigen::Vector4d::Constant( 1.0 );
            if( estimationCase < 2 )
            {
                fullEstimatableParameterSet = createEstimatableParameters( bodies, 1.1E7 );
            }
            else
            {
                fullEstimatableParameterSet = createEstimatableParameters( bodies, 1.1E7, true );
                parameterPerturbationMultipliers( 2 ) = 1.0E-4;
            }

            testObservationPartials< 1 >(
                        oneWayDopplerModel, bodies, fullEstimatableParameterSet, linkEnds, one_way_doppler, 1.0E-5,
                        true, true, 10.0, parameterPerturbationMultipliers );
            std::cout << "Case " << estimationCase << std::endl;

        }
    }

    // Test partials with real ephemerides (without test of position partials)
    {
        // Create environment
        SystemOfBodies bodies = setupEnvironment( groundStations, 1.0E7, 1.2E7, 1.1E7, false );

        // Set link ends for observation model
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = groundStations[ 1 ];
        linkEnds[ receiver ] = groundStations[ 0 ];

        for( unsigned int estimationCase  = 0; estimationCase  < 3; estimationCase ++ )
        {
            std::cout << "Rates: " << estimationCase << std::endl;
            // Generate one-way doppler model
            std::shared_ptr< ObservationModel< 1 > > oneWayDopplerModel;
            std::vector< std::string > perturbingBodies;
            perturbingBodies.push_back( "Earth" );
            if( estimationCase  == 0 )
            {
                oneWayDopplerModel =
                        observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                            std::make_shared< observation_models::ObservationModelSettings >(
                                observation_models::one_way_doppler, linkEnds,
                                std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                                    perturbingBodies ) ), bodies  );
            }
            else
            {
                oneWayDopplerModel =
                        observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                            std::make_shared< OneWayDopplerObservationSettings >
                            (  linkEnds, std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                                   perturbingBodies ),
                               std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Mars" ),
                               std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Earth" ) ), bodies  );
            }
            // Create parameter objects.
            std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet;
            Eigen::VectorXd parameterPerturbationMultipliers = Eigen::Vector4d::Constant( 1.0 );
            if( estimationCase < 2 )
            {
                fullEstimatableParameterSet = createEstimatableParameters( bodies, 1.1E7 );
            }
            else
            {
                fullEstimatableParameterSet = createEstimatableParameters( bodies, 1.1E7, true );
                parameterPerturbationMultipliers( 2 ) = 1.0E-4;
            }

            testObservationPartials< 1 >(
                        oneWayDopplerModel, bodies, fullEstimatableParameterSet, linkEnds, one_way_doppler, 1.0E-4, false, true,
                        1.0, parameterPerturbationMultipliers );
        }
    }


    // Test partials with constant ephemerides (allows test of position partials)
    {
        // Create environment
        SystemOfBodies bodies = setupEnvironment( groundStations, 1.0E7, 1.2E7, 1.1E7, true, 1000000.0 );

        // Set link ends for observation model (Mars to Earth)
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = groundStations[ 1 ];
        linkEnds[ receiver ] = groundStations[ 0 ];


        // Create one-way doppler model
        std::shared_ptr< OneWayDopplerObservationModel< > > oneWayDopplerModel =
                std::dynamic_pointer_cast< OneWayDopplerObservationModel< > >(
                    observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                        std::make_shared< OneWayDopplerObservationSettings >
                        (  linkEnds, std::shared_ptr< LightTimeCorrectionSettings >( ),
                           std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Earth" ),
                           std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Mars" ) ), bodies ) );


        // Extract proper time calculators
        std::shared_ptr< DopplerProperTimeRateInterface > receiverProperTimeRateCalculator =
                oneWayDopplerModel->getReceiverProperTimeRateCalculator( );
        std::shared_ptr< DopplerProperTimeRateInterface > transmitterProperTimeRateCalculator =
                oneWayDopplerModel->getTransmitterProperTimeRateCalculator( );

        // Create parameter objects.
        std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet =
                createEstimatableParameters( bodies, 1.1E7 );

        // Create partials for Doppler with proper time rates
        std::map< LinkEnds, std::shared_ptr< ObservationModel< 1 > > > observationModelList;
        observationModelList[ linkEnds ] = oneWayDopplerModel;
        std::map< LinkEnds, std::pair< SingleLinkObservationPartialList, std::shared_ptr< PositionPartialScaling > > > dopplerPartials =
                createOneWayDopplerPartials( observationModelList, bodies, fullEstimatableParameterSet );

        // Retrieve  scaling objects and partials with proper time
        std::shared_ptr< OneWayDopplerScaling > partialScalingObject =
                std::dynamic_pointer_cast< OneWayDopplerScaling >( dopplerPartials.begin( )->second.second );

        std::shared_ptr< OneWayDopplerProperTimeComponentScaling > transmitterProperTimePartials =
                partialScalingObject->getTransmitterProperTimePartials( );
        std::shared_ptr< OneWayDopplerProperTimeComponentScaling > receiverProperTimePartials =
                partialScalingObject->getReceiverProperTimePartials( );

        std::shared_ptr< OneWayDopplerPartial > earthStatePartial =
                std::dynamic_pointer_cast< OneWayDopplerPartial >(
                    ( dopplerPartials.begin( )->second.first ).begin( )->second );
        std::shared_ptr< OneWayDopplerPartial > marsStatePartial =
                std::dynamic_pointer_cast< OneWayDopplerPartial >(
                    ( ++( ( dopplerPartials.begin( )->second.first ).begin( ) ) )->second );

        // Compute nominal observation with proper time
        double observationTime = 1.1E7;
        std::vector< double > linkEndTimes;
        std::vector< Eigen::Vector6d > linkEndStates;
        LinkEndType referenceLinkEnd = transmitter;
        Eigen::VectorXd nominalObservable = oneWayDopplerModel->computeIdealObservationsWithLinkEndData(
                    observationTime, referenceLinkEnd, linkEndTimes, linkEndStates );

        // Compute partials with proper time.
        partialScalingObject->update(
                    linkEndStates, linkEndTimes, referenceLinkEnd, nominalObservable );
        std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > earthStatePartialOutput =
                earthStatePartial->calculatePartial( linkEndStates, linkEndTimes, referenceLinkEnd, nominalObservable );
        std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > marsStatePartialOutput =
                marsStatePartial->calculatePartial( linkEndStates, linkEndTimes, referenceLinkEnd, nominalObservable );

        // Compute numerical proper time rate partials and compare to analytical results
        {
            std::function< Eigen::VectorXd( const double ) > transmitterProperTimeRateFunction =
                    std::bind( &getProperTimeRateInVectorForm,
                                 transmitterProperTimeRateCalculator,
                                 linkEndTimes, linkEndStates, referenceLinkEnd );
            Eigen::Matrix< double, Eigen::Dynamic, 3 > numericalTransmitterProperTimePartialsWrtMarsPosition =
                    calculatePartialWrtConstantBodyState(
                        "Earth", bodies, Eigen::Vector3d::Constant( 1000.0E3 ), transmitterProperTimeRateFunction, 1.1E7, 1 );
            Eigen::Matrix< double, Eigen::Dynamic, 3 > numericalTransmitterProperTimePartialsWrtEarthPosition =
                    calculatePartialWrtConstantBodyState(
                        "Mars", bodies, Eigen::Vector3d::Constant( 1000.0E3 ), transmitterProperTimeRateFunction, 1.1E7, 1 );
            Eigen::Matrix< double, Eigen::Dynamic, 3 > numericalTransmitterProperTimePartialsWrtMarsVelocity =
                    calculatePartialWrtConstantBodyVelocity(
                        "Earth", bodies, Eigen::Vector3d::Constant( 1.0E0 ), transmitterProperTimeRateFunction, 1.1E7, 1 );
            Eigen::Matrix< double, Eigen::Dynamic, 3 > numericalTransmitterProperTimePartialsWrtEarthVelocity =
                    calculatePartialWrtConstantBodyVelocity(
                        "Mars", bodies, Eigen::Vector3d::Constant( 1.0E0 ), transmitterProperTimeRateFunction, 1.1E7, 1 );

            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        ( transmitterProperTimePartials->getPositionScalingFactor( transmitter ) ),
                        numericalTransmitterProperTimePartialsWrtMarsPosition, 1.0E-6 );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        ( transmitterProperTimePartials->getPositionScalingFactor( receiver ) ),
                        numericalTransmitterProperTimePartialsWrtEarthPosition, 1.0E-6 );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        ( transmitterProperTimePartials->getVelocityScalingFactor( transmitter ) ),
                        numericalTransmitterProperTimePartialsWrtMarsVelocity, 1.0E-6 );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        ( transmitterProperTimePartials->getVelocityScalingFactor( receiver ) ),
                        numericalTransmitterProperTimePartialsWrtEarthVelocity, 1.0E-6 );

            std::function< Eigen::VectorXd( const double ) > receiverProperTimeRateFunction =
                    std::bind( &getProperTimeRateInVectorForm,
                                 receiverProperTimeRateCalculator,
                                 linkEndTimes, linkEndStates, referenceLinkEnd );
            Eigen::Matrix< double, Eigen::Dynamic, 3 > numericalReceiverProperTimePartialsWrtMarsPosition =
                    calculatePartialWrtConstantBodyState(
                        "Earth", bodies, Eigen::Vector3d::Constant( 10000.0 ), receiverProperTimeRateFunction, 1.1E7, 1 );
            Eigen::Matrix< double, Eigen::Dynamic, 3 > numericalReceiverProperTimePartialsWrtEarthPosition =
                    calculatePartialWrtConstantBodyState(
                        "Mars", bodies, Eigen::Vector3d::Constant( 10000.0 ), receiverProperTimeRateFunction, 1.1E7, 1 );
            Eigen::Matrix< double, Eigen::Dynamic, 3 > numericalReceiverProperTimePartialsWrtMarsVelocity =
                    calculatePartialWrtConstantBodyVelocity(
                        "Earth", bodies, Eigen::Vector3d::Constant( 1000.0 ), receiverProperTimeRateFunction, 1.1E7, 1 );
            Eigen::Matrix< double, Eigen::Dynamic, 3 > numericalReceiverProperTimePartialsWrtEarthVelocity =
                    calculatePartialWrtConstantBodyVelocity(
                        "Mars", bodies, Eigen::Vector3d::Constant( 1000.0 ), receiverProperTimeRateFunction, 1.1E7, 1 );

            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        ( receiverProperTimePartials->getPositionScalingFactor( receiver ) ),
                        numericalReceiverProperTimePartialsWrtEarthPosition, 1.0E-6 );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        ( receiverProperTimePartials->getPositionScalingFactor( transmitter ) ),
                        numericalReceiverProperTimePartialsWrtMarsPosition, 1.0E-6 );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        ( receiverProperTimePartials->getVelocityScalingFactor( transmitter ) ),
                        numericalReceiverProperTimePartialsWrtMarsVelocity, 1.0E-6 );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        ( receiverProperTimePartials->getVelocityScalingFactor( receiver ) ),
                        numericalReceiverProperTimePartialsWrtEarthVelocity, 1.0E-6 );
        }


        // Create one-way doppler model without proper time rates
        std::shared_ptr< OneWayDopplerObservationModel< > > oneWayDopplerModelWithoutProperTime =
                std::dynamic_pointer_cast< OneWayDopplerObservationModel< > >(
                    observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                        std::make_shared< ObservationModelSettings >
                        (  one_way_doppler, linkEnds, std::shared_ptr< LightTimeCorrectionSettings >( ) ), bodies ) );

        // Create partials for Doppler without proper time rates
        observationModelList.clear( );
        observationModelList[ linkEnds ] = oneWayDopplerModelWithoutProperTime;
        std::map< LinkEnds, std::pair< SingleLinkObservationPartialList, std::shared_ptr< PositionPartialScaling > > > dopplerPartialsWithoutProperTime =
                createOneWayDopplerPartials( observationModelList, bodies, fullEstimatableParameterSet );

        // Retrieve partial object without proper time
        std::shared_ptr< OneWayDopplerScaling > partialScalingObjectWithoutProperTime =
                std::dynamic_pointer_cast< OneWayDopplerScaling >( dopplerPartialsWithoutProperTime.begin( )->second.second );
        std::shared_ptr< OneWayDopplerPartial > earthStatePartialWithoutProperTime =
                std::dynamic_pointer_cast< OneWayDopplerPartial >(
                    ( dopplerPartialsWithoutProperTime.begin( )->second.first ).begin( )->second );
        std::shared_ptr< OneWayDopplerPartial > marsStatePartialWithoutProperTime =
                std::dynamic_pointer_cast< OneWayDopplerPartial >(
                    ( ++( ( dopplerPartialsWithoutProperTime.begin( )->second.first ).begin( ) ) )->second );

        // Compute nominal observation without proper time
        std::vector< double > linkEndTimesWithoutProperTime;
        std::vector< Eigen::Vector6d > linkEndStatesWithoutProperTime;
        Eigen::VectorXd nominalObservableWithoutProperTime = oneWayDopplerModelWithoutProperTime->computeIdealObservationsWithLinkEndData(
                    observationTime, referenceLinkEnd, linkEndTimesWithoutProperTime, linkEndStatesWithoutProperTime );

        // Compute partials with proper time.
        partialScalingObjectWithoutProperTime->update(
                    linkEndStatesWithoutProperTime, linkEndTimesWithoutProperTime,
                    referenceLinkEnd, nominalObservableWithoutProperTime );
        std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > earthStatePartialOutputWithoutProperTime =
                earthStatePartialWithoutProperTime->calculatePartial(
                    linkEndStates, linkEndTimes, referenceLinkEnd, nominalObservable );
        std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > marsStatePartialOutputWithoutProperTime =
                marsStatePartialWithoutProperTime->calculatePartial(
                    linkEndStates, linkEndTimes, referenceLinkEnd, nominalObservable );

        Eigen::MatrixXd partialWrtEarthState = earthStatePartialOutput.at( 0 ).first;
        Eigen::MatrixXd partialWrtEarthStateWithoutProperTime = earthStatePartialOutputWithoutProperTime.at( 0 ).first;

        Eigen::MatrixXd partialWrtMarsState = marsStatePartialOutput.at( 0 ).first;
        Eigen::MatrixXd partialWrtMarsStateWithoutProperTime = marsStatePartialOutputWithoutProperTime.at( 0 ).first;

        Eigen::MatrixXd properTimePartialWrtMarsPosition = transmitterProperTimePartials->getPositionScalingFactor( transmitter );
        Eigen::MatrixXd properTimePartialWrtEarthPosition = receiverProperTimePartials->getPositionScalingFactor( receiver );

        Eigen::MatrixXd properTimePartialWrtMarsVelocity = transmitterProperTimePartials->getVelocityScalingFactor( transmitter );
        Eigen::MatrixXd properTimePartialWrtEarthVelocity = receiverProperTimePartials->getVelocityScalingFactor( receiver );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    ( ( partialWrtMarsState - partialWrtMarsStateWithoutProperTime ).block( 0, 0, 1, 3 ) ),
                    properTimePartialWrtMarsPosition, 1.0E-9 );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    ( -( partialWrtEarthState - partialWrtEarthStateWithoutProperTime ).block( 0, 0, 1, 3 ) ),
                    properTimePartialWrtEarthPosition, 1.0E-9 );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    ( ( partialWrtMarsState - partialWrtMarsStateWithoutProperTime ).block( 0, 3, 1, 3 ) ),
                    properTimePartialWrtMarsVelocity, 1.0E-8 );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    ( -( partialWrtEarthState - partialWrtEarthStateWithoutProperTime ).block( 0, 3, 1, 3 ) ),
                    properTimePartialWrtEarthVelocity, 1.0E-8 );

    }
}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat




