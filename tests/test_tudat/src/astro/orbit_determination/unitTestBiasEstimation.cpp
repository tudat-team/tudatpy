/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <limits>

#include <boost/test/included/unit_test.hpp>

#include "tudat/simulation/estimation_setup/executeEarthOrbiterBiasEstimationTestCase.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_estimation_from_positions )

//! This test checks whether observation biases are correctly estimated, using a variety of different settings
//! for types of observables/biases.
BOOST_AUTO_TEST_CASE( test_EstimationFromPosition )
{
    for( int estimateRangeBiases = 0; estimateRangeBiases < 2; estimateRangeBiases++ )
    {
        for( int estimateTwoWayBiases = 0; estimateTwoWayBiases < 2; estimateTwoWayBiases++ )
        {
            for( int useSingleBiasModel = 0; useSingleBiasModel < 2; useSingleBiasModel++ )
            {
                for( int estimateAbsoluteBiases = 0; estimateAbsoluteBiases < 2; estimateAbsoluteBiases++ )
                {
                    for( int estimateMultiArcBiases = 0; estimateMultiArcBiases < 2; estimateMultiArcBiases++ )
                    {
                        for( int estimateTimeBiases = 0; estimateTimeBiases < 2; estimateTimeBiases++ )
                        {
                            if( !( ( static_cast< bool >( estimateTwoWayBiases ) == true ) &&
                                   ( static_cast< bool >( estimateRangeBiases ) == false ) ) )
                            {
                                std::cout << "=========== Running Case: " << estimateRangeBiases << " " << estimateTwoWayBiases << " "
                                          << useSingleBiasModel << " " << estimateAbsoluteBiases << " " << estimateMultiArcBiases << " "
                                          << estimateTimeBiases << std::endl;

                                // Simulate estimated parameter error.
                                Eigen::VectorXd totalError = executeEarthOrbiterBiasEstimation< double, double >( estimateRangeBiases,
                                                                                                                  estimateTwoWayBiases,
                                                                                                                  useSingleBiasModel,
                                                                                                                  estimateAbsoluteBiases,
                                                                                                                  false,
                                                                                                                  estimateMultiArcBiases,
                                                                                                                  estimateTimeBiases )
                                                                     .first;

                                for( unsigned int j = 0; j < 3; j++ )
                                {
                                    BOOST_CHECK_SMALL( std::fabs( totalError( j ) ), 1.0E-5 );
                                    BOOST_CHECK_SMALL( std::fabs( totalError( j + 3 ) ), 1.0E-8 );
                                }

                                for( unsigned int j = 6; j < totalError.rows( ); j++ )
                                {
                                    if( estimateTimeBiases )
                                    {
                                        if( estimateRangeBiases || estimateMultiArcBiases )
                                        {
                                            BOOST_CHECK_SMALL( std::fabs( totalError( j ) ), 2.0E-11 );
                                        }
                                        else
                                        {
                                            BOOST_CHECK_SMALL( std::fabs( totalError( j ) ), 2.0E-14 );
                                        }
                                    }
                                    else
                                    {
                                        if( estimateAbsoluteBiases )
                                        {
                                            if( estimateRangeBiases )
                                            {
                                                BOOST_CHECK_SMALL( std::fabs( totalError( j ) ), 2.0E-7 );
                                            }
                                            else
                                            {
                                                BOOST_CHECK_SMALL( std::fabs( totalError( j ) ), 1.0E-10 );
                                            }
                                        }
                                        else if( !estimateMultiArcBiases )
                                        {
                                            BOOST_CHECK_SMALL( std::fabs( totalError( j ) ), 2.0E-14 );
                                        }
                                        else
                                        {
                                            BOOST_CHECK_SMALL( std::fabs( totalError( j ) ), 1.0E-13 );
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    BOOST_CHECK_EQUAL( executeEarthOrbiterBiasEstimation( true, false, true, true, true, false ).second, true );
}

//! Add matching model and simulation settings, and record angular links for the selection checks.
void addSharedBiasTestObservationModel( const std::vector< double >& observationTimes,
                                        std::vector< std::shared_ptr< ObservationModelSettings > >& models,
                                        std::vector< std::shared_ptr< ObservationSimulationSettings< double > > >& simulations,
                                        std::vector< LinkEnds >& angularLinks,
                                        const ObservableType type,
                                        const LinkEnds& links,
                                        const std::shared_ptr< ObservationBiasSettings >& bias )
{
    models.push_back( std::make_shared< ObservationModelSettings >( type, links, nullptr, bias ) );
    simulations.push_back( std::make_shared< TabulatedObservationSimulationSettings< double > >(
            type, links, observationTimes, type == position_observable ? observed_body : receiver ) );
    if( type == angular_position )
    {
        angularLinks.push_back( links );
    }
}

//! Create a separate angular-bias setting for each selected model, initialized to the common truth values.
std::shared_ptr< ObservationBiasSettings > createSharedAngularBiasTestSettings( const ObservationBiasTypes biasType,
                                                                                const Eigen::VectorXd& trueAngularBias,
                                                                                const std::vector< double >& arcTimes,
                                                                                const LinkEndType timeLinkEnd )
{
    if( biasType == arc_wise_constant_absolute_bias )
    {
        return arcWiseAbsoluteBias( arcTimes, { trueAngularBias.head( 2 ), trueAngularBias.tail( 2 ) }, timeLinkEnd );
    }
    return biasType == constant_absolute_bias ? constantAbsoluteBias( trueAngularBias ) : constantRelativeBias( trueAngularBias );
}

//! Build one angular model per target after the observer in the body-name list, using the supplied biases.
std::vector< std::shared_ptr< ObservationModelSettings > > createSharedBiasTestObservationModels(
        const std::vector< std::string >& names,
        const LinkEndId& observerId,
        const std::vector< std::shared_ptr< ObservationBiasSettings > >& biases )
{
    std::vector< std::shared_ptr< ObservationModelSettings > > models;
    for( unsigned int i = 0; i < biases.size( ); ++i )
    {
        const LinkEnds links = { { transmitter, LinkEndId( names.at( i + 1 ), "" ) }, { receiver, observerId } };
        models.push_back( std::make_shared< ObservationModelSettings >( angular_position, links, nullptr, biases.at( i ) ) );
    }
    return models;
}

//! Exercise shared-bias validation through the normal parameter and estimator creation interfaces.
std::shared_ptr< OrbitDeterminationManager< double, double > > createSharedBiasTestEstimator(
        const SystemOfBodies& bodies,
        const std::vector< std::string >& names,
        const LinkEndId& observerId,
        const std::shared_ptr< EstimatableParameterSettings >& settings,
        const std::vector< std::shared_ptr< ObservationBiasSettings > >& biases,
        const std::vector< std::shared_ptr< EstimatableParameterSettings > >& additionalSettings = {} )
{
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterSettings = { settings };
    parameterSettings.insert( parameterSettings.end( ), additionalSettings.begin( ), additionalSettings.end( ) );
    auto parameters = createParametersToEstimate< double, double >( parameterSettings, bodies );
    const std::shared_ptr< propagators::PropagatorSettings< double > > noPropagator;
    return std::make_shared< OrbitDeterminationManager< double, double > >(
            bodies, parameters, createSharedBiasTestObservationModels( names, observerId, biases ), noPropagator );
}

//! Create Earth stations, a nearly circular observer orbit, and three higher target orbits under Earth gravity.
SystemOfBodies createSharedBiasTestBodies( )
{
    const double earthGravitationalParameter = spice_interface::getBodyGravitationalParameter( "Earth" );
    BodyListSettings bodySettings = getDefaultBodySettings( { "Earth" }, "Earth", "J2000" );
    bodySettings.at( "Earth" )->gravityFieldSettings = centralGravitySettings( earthGravitationalParameter );
    bodySettings.at( "Earth" )->groundStationSettings = {
        groundStationSettings( "Station1", Eigen::Vector3d( 0.0, 0.35, 0.0 ), coordinate_conversions::geodetic_position ),
        groundStationSettings( "Station2", Eigen::Vector3d( 0.0, -0.55, 2.0 ), coordinate_conversions::geodetic_position )
    };

    // The observer follows a prescribed two-body orbit with a 10,000 km semi-major axis and eccentricity 0.01.
    Eigen::Vector6d observerElements;
    observerElements << 1.0E7, 0.01, unit_conversions::convertDegreesToRadians( 25.0 ), unit_conversions::convertDegreesToRadians( 35.0 ),
            unit_conversions::convertDegreesToRadians( 15.0 ), unit_conversions::convertDegreesToRadians( 80.0 );
    bodySettings.addSettings( "ObserverSatellite" );
    bodySettings.at( "ObserverSatellite" )->ephemerisSettings =
            keplerEphemerisSettings( observerElements, 0.0, earthGravitationalParameter, "Earth", "J2000" );

    // The targets have distinct 20,000, 24,000, and 28,000 km orbits, all above the observer's orbit.
    const std::vector< std::string > targets = { "A", "B", "C" };
    for( unsigned int i = 0; i < targets.size( ); ++i )
    {
        Eigen::Vector6d targetElements;
        targetElements << 2.0E7 + i * 4.0E6, 0.03 + i * 0.01, unit_conversions::convertDegreesToRadians( 30.0 + i * 15.0 ),
                unit_conversions::convertDegreesToRadians( 40.0 + i * 20.0 ), unit_conversions::convertDegreesToRadians( 20.0 + i * 70.0 ),
                unit_conversions::convertDegreesToRadians( 10.0 + i * 90.0 );
        bodySettings.addSettings( targets.at( i ) );
        bodySettings.at( targets.at( i ) )->ephemerisSettings =
                keplerEphemerisSettings( targetElements, 0.0, earthGravitationalParameter, "Earth", "J2000" );
    }
    return createSystemOfBodies( bodySettings );
}

//! Exercise the public setup, propagation, simulation, and estimation interfaces together.
//! Estimate three Earth-orbiting target states and shared angular/range biases, selecting three of eight angular links.
//! Repeat for a spacecraft receiver and an Earth station to test both body and station-specific selection.
//! The other links differ in receiver body, reference point, or role, so incorrect selection changes the fitted data.
//! Cover constant absolute, constant relative, and arc-wise absolute biases, including both arc reference-time roles.
//! After estimation, check the shared partials against the theoretical values and active arc blocks.
BOOST_AUTO_TEST_CASE( test_SharedObservationBiasEstimation )
{
    spice_interface::loadStandardSpiceKernels( );

    using namespace observation_models;
    using namespace estimatable_parameters;
    using namespace simulation_setup;

    // For each receiver, the last two bias cases use different link-end event times to choose the active arc.
    const std::vector< ObservationBiasTypes > biasTypes = {
        constant_absolute_bias, constant_relative_bias, arc_wise_constant_absolute_bias, arc_wise_constant_absolute_bias
    };
    const std::vector< LinkEndId > receiverIds = { LinkEndId( "ObserverSatellite", "" ), LinkEndId( "Earth", "Station1" ) };
    const double relativeCaseTimeBias = 15.0;
    for( unsigned int testCase = 0; testCase < receiverIds.size( ) * biasTypes.size( ); ++testCase )
    {
        BOOST_TEST_CONTEXT( "shared bias case " << testCase )
        {
            const unsigned int biasCase = testCase % biasTypes.size( );
            const unsigned int receiverCase = testCase / biasTypes.size( );
            const auto biasType = biasTypes.at( biasCase );
            const LinkEndId selectedReceiver = receiverIds.at( receiverCase );
            const LinkEndId otherReceiver = receiverIds.at( 1 - receiverCase );
            const bool arcWise = biasType == arc_wise_constant_absolute_bias;
            const LinkEndType timeLinkEnd = biasCase == 3 ? transmitter : receiver;
            const std::vector< double > arcTimes = arcWise ? std::vector< double >{ 0.0, 3600.0 } : std::vector< double >{};
            Eigen::VectorXd trueAngularBias( arcWise ? 4 : 2 );
            trueAngularBias.head< 2 >( ) = Eigen::Vector2d( 2.0E-5, -3.0E-5 );
            if( arcWise )
            {
                trueAngularBias.tail< 2 >( ) = Eigen::Vector2d( -4.0E-5, 5.0E-5 );
            }
            if( biasType == constant_relative_bias )
            {
                trueAngularBias *= 100.0;
            }

            SystemOfBodies bodies = createSharedBiasTestBodies( );
            const std::vector< std::string > targets = { "A", "B", "C" };
            const std::vector< std::string > centralBodies( targets.size( ), "Earth" );
            Eigen::VectorXd initialState( 18 );
            SelectedAccelerationMap accelerationSettings;

            // Start from the Keplerian orbit, then let numerical propagation supply each target's ephemeris.
            for( unsigned int i = 0; i < targets.size( ); ++i )
            {
                initialState.segment< 6 >( 6 * i ) = bodies.at( targets.at( i ) )->getStateInBaseFrameFromEphemeris( 0.0 );
                bodies.at( targets.at( i ) )
                        ->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris<> >(
                                std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector6d > >( ),
                                "Earth",
                                "J2000" ) );

                // Earth point-mass gravity is the only acceleration acting on each propagated target.
                accelerationSettings[ targets.at( i ) ][ "Earth" ] = { pointMassGravityAcceleration( ) };
            }
            const auto accelerations = createAccelerationModelsMap( bodies, accelerationSettings, targets, centralBodies );
            auto propagator = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >(
                    centralBodies,
                    accelerations,
                    targets,
                    initialState,
                    0.0,
                    numerical_integrators::rungeKuttaFixedStepSettings( 60.0, numerical_integrators::rungeKuttaFehlberg78 ),
                    std::make_shared< propagators::PropagationTimeTerminationSettings >( 7200.0 ) );

            std::vector< double > observationTimes;
            for( double time = 600.0; time <= 6600.0; time += 300.0 )
            {
                observationTimes.push_back( time );
            }
            std::vector< std::shared_ptr< ObservationModelSettings > > models;
            std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > simulations;
            std::vector< LinkEnds > angularLinks;

            // Fixed RA/Dec offsets (radians) on angular links outside the shared receiver selection.
            // These biases are not estimated: shared-parameter updates must leave their values unchanged.
            const Eigen::Vector2d unsharedAngularBias( 7.0E-6, 9.0E-6 );
            for( const auto& target : targets )
            {
                const LinkEnds selected = { { transmitter, LinkEndId( target, "" ) }, { receiver, selectedReceiver } };

                // Also exercise finding a bias within a combined model.
                const auto secondaryBias = biasType == constant_relative_bias ? constantAbsoluteBias( Eigen::Vector2d( 4.0E-6, -2.0E-6 ) )
                                                                              : constantRelativeBias( Eigen::Vector2d::Zero( ) );
                std::vector< std::shared_ptr< ObservationBiasSettings > > combinedBiases = {
                    createSharedAngularBiasTestSettings( biasType, trueAngularBias, arcTimes, timeLinkEnd ), secondaryBias
                };

                // Shared relative partials must recompute at the biased event time, not the nominal observation time.
                if( biasType == constant_relative_bias )
                {
                    combinedBiases.push_back( constantTimeBias( relativeCaseTimeBias, receiver ) );
                }
                addSharedBiasTestObservationModel( observationTimes,
                                                   models,
                                                   simulations,
                                                   angularLinks,
                                                   angular_position,
                                                   selected,
                                                   multipleObservationBiasSettings( combinedBiases ) );
                addSharedBiasTestObservationModel( observationTimes,
                                                   models,
                                                   simulations,
                                                   angularLinks,
                                                   one_way_range,
                                                   selected,
                                                   constantAbsoluteBias( Eigen::VectorXd::Constant( 1, 12.0 ) ) );
                addSharedBiasTestObservationModel( observationTimes,
                                                   models,
                                                   simulations,
                                                   angularLinks,
                                                   angular_position,
                                                   { { transmitter, LinkEndId( target, "" ) }, { receiver, otherReceiver } },
                                                   constantAbsoluteBias( unsharedAngularBias ) );

                // Independent position observations make all three target states identifiable.
                addSharedBiasTestObservationModel( observationTimes,
                                                   models,
                                                   simulations,
                                                   angularLinks,
                                                   position_observable,
                                                   { { observed_body, LinkEndId( target, "" ) } },
                                                   nullptr );
            }
            addSharedBiasTestObservationModel( observationTimes,
                                               models,
                                               simulations,
                                               angularLinks,
                                               angular_position,
                                               { { transmitter, LinkEndId( "A", "" ) }, { receiver, LinkEndId( "Earth", "Station2" ) } },
                                               constantAbsoluteBias( unsharedAngularBias ) );
            addSharedBiasTestObservationModel( observationTimes,
                                               models,
                                               simulations,
                                               angularLinks,
                                               angular_position,
                                               { { transmitter, selectedReceiver }, { receiver, LinkEndId( "A", "" ) } },
                                               constantAbsoluteBias( unsharedAngularBias ) );

            // Ensure the fixture contains five angular links that must be excluded as well as the three selected links.
            BOOST_REQUIRE_EQUAL( angularLinks.size( ), 8 );

            std::vector< std::shared_ptr< EstimatableParameterSettings > > settings;
            for( int i = 0; i < 3; ++i )
            {
                settings.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
                        targets.at( i ), initialState.segment< 6 >( 6 * i ), "Earth" ) );
            }
            settings.push_back( sharedObservationBias( biasType,
                                                       angular_position,
                                                       receiver,
                                                       selectedReceiver,
                                                       arcTimes,
                                                       biasCase == 3 ? transmitter : unidentified_link_end ) );
            settings.push_back( sharedObservationBias( constant_absolute_bias, one_way_range, receiver, selectedReceiver ) );
            auto parameters = createParametersToEstimate< double, double >( settings, bodies, propagator );
            auto angularParameter =
                    std::dynamic_pointer_cast< SharedObservationBiasParameter >( parameters->getVectorParameters( ).at( 18 ) );
            auto rangeParameter = std::dynamic_pointer_cast< SharedObservationBiasParameter >(
                    parameters->getVectorParameters( ).at( 18 + trueAngularBias.size( ) ) );

            // Both factory settings must produce shared parameters, with separate vectors for angular and range biases.
            BOOST_REQUIRE( angularParameter != nullptr );
            BOOST_REQUIRE( rangeParameter != nullptr );

            OrbitDeterminationManager< double, double > manager( bodies, parameters, models, propagator );

            // Binding must find exactly three models for each observable and allocate one vector per shared parameter.
            BOOST_CHECK_EQUAL( angularParameter->getMembers( ).size( ), 3 );
            BOOST_CHECK_EQUAL( rangeParameter->getMembers( ).size( ), 3 );
            BOOST_CHECK_EQUAL( parameters->getParameterSetSize( ), 18 + trueAngularBias.size( ) + 1 );
            const Eigen::VectorXd truth = parameters->getFullParameterValues< double >( );
            auto observations = simulateObservations< double, double >( simulations, manager.getObservationSimulators( ), bodies );

            // Perturb both the target states and shared biases before estimating them jointly.
            Eigen::VectorXd initialGuess = truth;
            for( int i = 0; i < 3; ++i )
            {
                initialGuess.segment< 3 >( 6 * i ) += Eigen::Vector3d( 10.0, -20.0, 15.0 );
                initialGuess.segment< 3 >( 6 * i + 3 ) += Eigen::Vector3d( 0.01, -0.02, 0.03 );
            }
            initialGuess.segment( 18, trueAngularBias.size( ) ).array( ) += 1.0E-4;
            initialGuess.tail( 1 )( 0 ) += 5.0;
            manager.resetParameterEstimate( initialGuess );
            auto input = std::make_shared< EstimationInput< double, double > >( observations );
            input->defineEstimationSettings( true, true, true, false, true, false );
            auto output = manager.estimateParameters( input );

            // Estimation must complete its inversion successfully.
            BOOST_REQUIRE( !output->exceptionDuringInversion_ );

            // Recover all three target states and both shared biases from perturbed guesses, leaving small residuals.
            BOOST_CHECK_SMALL( ( output->parameterEstimate_.head( 18 ) - truth.head( 18 ) ).norm( ), 1.0E-5 );
            BOOST_CHECK_SMALL( ( output->parameterEstimate_.segment( 18, trueAngularBias.size( ) ) - trueAngularBias ).norm( ), 1.0E-11 );
            BOOST_CHECK_SMALL( std::abs( output->parameterEstimate_.tail( 1 )( 0 ) - 12.0 ), 1.0E-6 );
            BOOST_CHECK_SMALL( output->residuals_.cwiseAbs( ).maxCoeff( ), 1.0E-5 );

            // Each underlying selected model must receive the final common angular-bias value.
            for( const auto& member : angularParameter->getMembers( ) )
            {
                BOOST_CHECK_SMALL( ( member.second->getParameterValue( ) - trueAngularBias ).norm( ), 1.0E-11 );
            }

            // Check the fitted angular models against theoretical partials before, at, and after the arc boundary.
            const auto angularManager = manager.getObservationManagers( ).at( angular_position );
            const auto simulator = std::dynamic_pointer_cast< ObservationSimulator< 2 > >( angularManager->getObservationSimulator( ) );
            BOOST_REQUIRE( simulator != nullptr );
            const std::vector< double > checkTimes = { 1800.0, 3600.0, 5400.0 };
            for( const auto& links : angularLinks )
            {
                Eigen::VectorXd values;
                Eigen::MatrixXd partials;
                angularManager->computeObservationsWithPartials( checkTimes, links, receiver, nullptr, values, partials );
                Eigen::MatrixXd expectedPartials = Eigen::MatrixXd::Zero( 2 * checkTimes.size( ), trueAngularBias.size( ) + 1 );

                // Absolute partials form an identity block in the active arc; other arc blocks stay zero.
                // At reception time 3600, transmission is still in the first arc.
                if( links.at( receiver ) == selectedReceiver )
                {
                    for( unsigned int i = 0; i < checkTimes.size( ); ++i )
                    {
                        const int column =
                                arcWise && ( timeLinkEnd == receiver ? checkTimes.at( i ) >= 3600.0 : checkTimes.at( i ) > 3600.0 ) ? 2 : 0;
                        expectedPartials.block< 2, 2 >( 2 * i, column ).setIdentity( );

                        // Relative partials equal the ideal RA/Dec at the time shifted by the fixed receiver time bias.
                        if( biasType == constant_relative_bias )
                        {
                            expectedPartials.block< 2, 2 >( 2 * i, 0 ) =
                                    simulator->getObservationModel( links )
                                            ->computeIdealObservations( checkTimes.at( i ) - relativeCaseTimeBias, receiver )
                                            .asDiagonal( );
                        }
                    }
                }
                else
                {
                    const auto bias = std::dynamic_pointer_cast< ConstantObservationBias< 2 > >(
                            simulator->getObservationModel( links )->getObservationBiasCalculator( ) );

                    // Unselected links keep their fixed biases and have zero shared-bias partials.
                    BOOST_REQUIRE( bias != nullptr );
                    BOOST_CHECK_SMALL( ( bias->getConstantObservationBias( ) - unsharedAngularBias ).norm( ), 1.0E-30 );
                }

                // Compare all shared columns, including off-diagonal zeros and the zero range-bias column.
                BOOST_CHECK_SMALL( ( partials.rightCols( expectedPartials.cols( ) ) - expectedPartials ).norm( ), 1.0E-14 );
            }

            // Range observations have unit partials for their shared absolute bias and zero angular-bias partials.
            Eigen::VectorXd rangeValues;
            Eigen::MatrixXd rangePartials;
            manager.getObservationManagers( )
                    .at( one_way_range )
                    ->computeObservationsWithPartials( checkTimes, angularLinks.front( ), receiver, nullptr, rangeValues, rangePartials );
            Eigen::MatrixXd expectedRangePartials = Eigen::MatrixXd::Zero( checkTimes.size( ), trueAngularBias.size( ) + 1 );
            expectedRangePartials.rightCols( 1 ).setOnes( );
            BOOST_CHECK_SMALL( ( rangePartials.rightCols( expectedRangePartials.cols( ) ) - expectedRangePartials ).norm( ), 1.0E-15 );
        }
    }
}

//! Check invalid shared-bias settings, model compatibility, and values assigned before binding.
//! Reject unsupported types, invalid arc definitions, inconsistent parameter objects, and ambiguous or missing model matches.
//! These checks prevent silently fitting the wrong links or choosing an arbitrary model's initial value.
//! Verify that a pre-binding assignment initializes every member once, and that new estimators use their new model settings.
BOOST_AUTO_TEST_CASE( test_SharedObservationBiasValidation )
{
    spice_interface::loadStandardSpiceKernels( );

    using namespace observation_models;
    using namespace estimatable_parameters;
    using namespace simulation_setup;

    const LinkEndId observerId( "ObserverSatellite", "" );
    const auto constantSettings = sharedObservationBias( constant_absolute_bias, angular_position, receiver, observerId );

    // Exposing all enum values must not enable any of the eight unsupported shared-bias types.
    for( const auto type : { multiple_observation_biases,
                             arc_wise_constant_relative_bias,
                             constant_time_drift_bias,
                             arc_wise_time_drift_bias,
                             constant_time_bias,
                             arc_wise_time_bias,
                             clock_induced_bias,
                             two_way_range_time_scale_bias } )
    {
        BOOST_CHECK_THROW( sharedObservationBias( type, angular_position, receiver, observerId ), std::runtime_error );
    }

    // Arc-wise biases require arc boundaries; constant biases must not accept them.
    BOOST_CHECK_THROW( sharedObservationBias( arc_wise_constant_absolute_bias, angular_position, receiver, observerId ),
                       std::runtime_error );
    BOOST_CHECK_THROW( sharedObservationBias( constant_absolute_bias, angular_position, receiver, observerId, { 0.0 } ),
                       std::runtime_error );

    // Arc boundaries must be ordered, distinct, and finite so their parameter blocks are unambiguous.
    BOOST_CHECK_THROW( sharedObservationBias( arc_wise_constant_absolute_bias, angular_position, receiver, observerId, { 2.0, 1.0 } ),
                       std::runtime_error );
    BOOST_CHECK_THROW( sharedObservationBias( arc_wise_constant_absolute_bias, angular_position, receiver, observerId, { 0.0, 0.0 } ),
                       std::runtime_error );
    BOOST_CHECK_THROW( sharedObservationBias( arc_wise_constant_absolute_bias, angular_position, receiver, observerId, { TUDAT_NAN } ),
                       std::runtime_error );

    // Prescribed Earth orbits suffice here: these checks validate model binding without fitting target states.
    const std::vector< std::string > names = { "ObserverSatellite", "A", "B", "C" };
    SystemOfBodies bodies = createSharedBiasTestBodies( );
    const std::shared_ptr< propagators::PropagatorSettings< double > > noPropagator;
    const auto zeroBias = constantAbsoluteBias( Eigen::Vector2d::Zero( ) );
    const auto otherBias = constantAbsoluteBias( Eigen::Vector2d::Ones( ) );

    // An inconsistent parameter identifier must produce an exception rather than a null-pointer dereference.
    const LinkEnds testLinks = { { transmitter, LinkEndId( "A", "" ) }, { receiver, observerId } };
    const auto inconsistentParameter =
            std::make_shared< SingleArcObservationBiasParameter >( shared_observation_bias, nullptr, nullptr, testLinks, angular_position );
    BOOST_CHECK_THROW( doesObservationParameterMatchObservable( inconsistentParameter, angular_position ), std::runtime_error );
    BOOST_CHECK_THROW( observation_partials::createObservationPartialWrtLinkProperty< 2 >(
                               testLinks, angular_position, inconsistentParameter, bodies ),
                       std::runtime_error );

    // Every selected model must contain the requested bias component.
    BOOST_CHECK_THROW( createSharedBiasTestEstimator( bodies, names, observerId, constantSettings, { zeroBias, nullptr, zeroBias } ),
                       std::runtime_error );

    // Different initial model values must be rejected unless a common value was explicitly assigned.
    BOOST_CHECK_THROW( createSharedBiasTestEstimator( bodies, names, observerId, constantSettings, { zeroBias, otherBias, zeroBias } ),
                       std::runtime_error );

    // Two components of the requested type in one model would make the shared selection ambiguous.
    BOOST_CHECK_THROW(
            createSharedBiasTestEstimator(
                    bodies, names, observerId, constantSettings, { multipleObservationBiasSettings( { zeroBias, otherBias } ) } ),
            std::runtime_error );

    // Nested combined biases must still contain exactly one component of the selected type.
    const auto nestedBias = multipleObservationBiasSettings(
            { constantRelativeBias( Eigen::Vector2d::Zero( ) ), multipleObservationBiasSettings( { zeroBias } ) } );
    BOOST_CHECK_NO_THROW(
            createSharedBiasTestEstimator( bodies, names, observerId, constantSettings, { nestedBias, zeroBias, zeroBias } ) );
    BOOST_CHECK_THROW(
            createSharedBiasTestEstimator(
                    bodies, names, observerId, constantSettings, { multipleObservationBiasSettings( { zeroBias, nestedBias } ) } ),
            std::runtime_error );

    // Selecting ObserverSatellite as transmitter must not accidentally match links on which it is the receiver.
    BOOST_CHECK_THROW(
            createSharedBiasTestEstimator( bodies,
                                           names,
                                           observerId,
                                           sharedObservationBias( constant_absolute_bias, angular_position, transmitter, observerId ),
                                           { zeroBias, zeroBias, zeroBias } ),
            std::runtime_error );

    // A shared range parameter must fail when the available models are all angular observations.
    BOOST_CHECK_THROW( createSharedBiasTestEstimator( bodies,
                                                      names,
                                                      observerId,
                                                      sharedObservationBias( constant_absolute_bias, one_way_range, receiver, observerId ),
                                                      { zeroBias, zeroBias, zeroBias } ),
                       std::runtime_error );
    const auto arcSettings =
            sharedObservationBias( arc_wise_constant_absolute_bias, angular_position, receiver, observerId, { 0.0, 10.0 } );
    const std::vector< Eigen::VectorXd > arcValues = { Eigen::Vector2d::Zero( ), Eigen::Vector2d::Ones( ) };
    const auto goodArc = arcWiseAbsoluteBias( { 0.0, 10.0 }, arcValues, receiver );

    // Equal vector sizes are insufficient: all members must use the same boundaries and arc reference-time role.
    BOOST_CHECK_THROW(
            createSharedBiasTestEstimator(
                    bodies, names, observerId, arcSettings, { goodArc, arcWiseAbsoluteBias( { 0.0, 11.0 }, arcValues, receiver ) } ),
            std::runtime_error );
    BOOST_CHECK_THROW(
            createSharedBiasTestEstimator(
                    bodies, names, observerId, arcSettings, { goodArc, arcWiseAbsoluteBias( { 0.0, 10.0 }, arcValues, transmitter ) } ),
            std::runtime_error );

    // A deliberate pre-closure assignment initializes every member, even if their settings differ.
    auto parameters = createParametersToEstimate< double, double >( { constantSettings }, bodies );
    const auto shared = std::dynamic_pointer_cast< SharedObservationBiasParameter >( parameters->getVectorParameters( ).at( 0 ) );

    // Verify the factory type before testing assignments made while the member list is still empty.
    BOOST_REQUIRE( shared != nullptr );
    const Eigen::Vector2d desired( 1.0E-4, -2.0E-4 );
    parameters->resetParameterValues( Eigen::VectorXd( desired ) );
    OrbitDeterminationManager< double, double > manager(
            bodies,
            parameters,
            createSharedBiasTestObservationModels( names, observerId, { zeroBias, otherBias, zeroBias } ),
            noPropagator );

    // Diagnostics and parameter identifiers should contain readable link-end and bias names.
    BOOST_CHECK( shared->getParameterDescription( ).find( "receiver" ) != std::string::npos );
    BOOST_CHECK( shared->getParameterDescription( ).find( "constant_absolute_bias" ) != std::string::npos );
    BOOST_CHECK( shared->getParameterName( ).second.second.find( "receiver" ) != std::string::npos );
    BOOST_CHECK( shared->getParameterName( ).second.second.find( "constant_absolute_bias" ) != std::string::npos );

    // The explicit value must override differing model initial values and be readable from every member.
    BOOST_CHECK_EQUAL( shared->getMembers( ).size( ), 3 );
    BOOST_CHECK_SMALL( ( shared->getParameterValue( ) - desired ).norm( ), 1.0E-30 );
    for( const auto& member : shared->getMembers( ) )
    {
        BOOST_CHECK_SMALL( ( member.second->getParameterValue( ) - desired ).norm( ), 1.0E-30 );
    }

    // Parameter updates must preserve the angular observable's two-component vector size.
    BOOST_CHECK_THROW( shared->setParameterValue( Eigen::Vector3d::Zero( ) ), std::runtime_error );

    // Detect an out-of-band change to one member instead of returning an arbitrary shared value.
    shared->getMembers( ).begin( )->second->setParameterValue( Eigen::Vector2d::Zero( ) );
    BOOST_CHECK_THROW( shared->getParameterValue( ), std::runtime_error );
    shared->setParameterValue( desired );
    auto simulator = std::dynamic_pointer_cast< ObservationSimulator< 2 > >(
            manager.getObservationManagers( ).at( angular_position )->getObservationSimulator( ) );
    performObservationParameterEstimationClosure( simulator, parameters );

    // Rebinding the same models rebuilds three members and reads their existing common value.
    BOOST_CHECK_EQUAL( shared->getMembers( ).size( ), 3 );
    BOOST_CHECK_SMALL( ( shared->getParameterValue( ) - desired ).norm( ), 1.0E-30 );

    // Reusing the parameter set must read the new models, not the old pre-binding assignment or an estimator update.
    OrbitDeterminationManager< double, double > secondManager(
            bodies,
            parameters,
            createSharedBiasTestObservationModels( names, observerId, { otherBias, otherBias, otherBias } ),
            noPropagator );
    BOOST_CHECK_SMALL( ( shared->getParameterValue( ) - Eigen::Vector2d::Ones( ) ).norm( ), 1.0E-30 );
    parameters->resetParameterValues( Eigen::VectorXd( desired ) );
    OrbitDeterminationManager< double, double > thirdManager(
            bodies,
            parameters,
            createSharedBiasTestObservationModels( names, observerId, { zeroBias, zeroBias, zeroBias } ),
            noPropagator );
    BOOST_CHECK_SMALL( shared->getParameterValue( ).norm( ), 1.0E-30 );

    // A reused parameter set must still reject inconsistent initial values in the replacement models.
    BOOST_CHECK_THROW( ( OrbitDeterminationManager< double, double >(
                               bodies,
                               parameters,
                               createSharedBiasTestObservationModels( names, observerId, { zeroBias, otherBias, zeroBias } ),
                               noPropagator ) ),
                       std::runtime_error );
}

//! Recognize ownership conflicts specifically, so unrelated setup errors cannot satisfy the regression tests.
bool isObservationBiasOwnershipError( const std::runtime_error& error )
{
    return std::string( error.what( ) ).find( "Only one parameter may own a bias component" ) != std::string::npos;
}

//! Reject ordinary/shared and shared/shared ownership conflicts for all three supported bias types in either order.
//! Distinct bias components on one link and disjoint link selections must remain valid, including nested combined biases.
BOOST_AUTO_TEST_CASE( test_SharedObservationBiasOwnership )
{
    spice_interface::loadStandardSpiceKernels( );
    const auto bodies = createSharedBiasTestBodies( );
    const std::vector< std::string > names = { "ObserverSatellite", "A", "B", "C" };
    const LinkEndId observerId( "ObserverSatellite", "" );
    const LinkEnds links = { { transmitter, LinkEndId( "A", "" ) }, { receiver, observerId } };
    for( const auto type : { constant_absolute_bias, constant_relative_bias, arc_wise_constant_absolute_bias } )
    {
        const bool arcWise = type == arc_wise_constant_absolute_bias;
        const std::vector< double > arcTimes = arcWise ? std::vector< double >{ 0.0, 10.0 } : std::vector< double >{};
        const auto timeRole = arcWise ? receiver : unidentified_link_end;
        const auto sharedSettings = sharedObservationBias( type, angular_position, receiver, observerId, arcTimes, timeRole );
        const auto ordinarySettings = arcWise ? arcwiseObservationBias( links, angular_position, arcTimes )
                                              : ( type == constant_absolute_bias ? observationBias( links, angular_position )
                                                                                 : relativeObservationBias( links, angular_position ) );
        const auto bias = createSharedAngularBiasTestSettings( type, Eigen::VectorXd::Zero( arcWise ? 4 : 2 ), arcTimes, receiver );
        const auto nestedBias = multipleObservationBiasSettings( { multipleObservationBiasSettings( { bias } ) } );

        // Both parameter orders must report ownership conflicts, even when the component is nested.
        BOOST_CHECK_EXCEPTION( createSharedBiasTestEstimator(
                                       bodies, names, observerId, sharedSettings, { nestedBias, bias, bias }, { ordinarySettings } ),
                               std::runtime_error,
                               isObservationBiasOwnershipError );
        BOOST_CHECK_EXCEPTION( createSharedBiasTestEstimator(
                                       bodies, names, observerId, ordinarySettings, { nestedBias, bias, bias }, { sharedSettings } ),
                               std::runtime_error,
                               isObservationBiasOwnershipError );

        // Shared receiver and transmitter selectors overlap on A's link and cannot both own that component.
        const auto transmitterSettings =
                sharedObservationBias( type, angular_position, transmitter, LinkEndId( "A", "" ), arcTimes, timeRole );
        BOOST_CHECK_EXCEPTION( createSharedBiasTestEstimator(
                                       bodies, names, observerId, sharedSettings, { nestedBias, bias, bias }, { transmitterSettings } ),
                               std::runtime_error,
                               isObservationBiasOwnershipError );
        BOOST_CHECK_EXCEPTION( createSharedBiasTestEstimator(
                                       bodies, names, observerId, transmitterSettings, { nestedBias, bias, bias }, { sharedSettings } ),
                               std::runtime_error,
                               isObservationBiasOwnershipError );

        // A separate bias component on the same link is a distinct estimation parameter.
        const auto otherSettings = type == constant_relative_bias ? observationBias( links, angular_position )
                                                                  : relativeObservationBias( links, angular_position );
        const auto otherBias = type == constant_relative_bias ? constantAbsoluteBias( Eigen::Vector2d::Zero( ) )
                                                              : constantRelativeBias( Eigen::Vector2d::Zero( ) );
        BOOST_CHECK_NO_THROW( createSharedBiasTestEstimator( bodies,
                                                             names,
                                                             observerId,
                                                             sharedSettings,
                                                             { multipleObservationBiasSettings( { nestedBias, otherBias } ), bias, bias },
                                                             { otherSettings } ) );

        // Shared parameters of the same type may independently own disjoint transmitter selections.
        const auto secondTransmitterSettings =
                sharedObservationBias( type, angular_position, transmitter, LinkEndId( "B", "" ), arcTimes, timeRole );
        BOOST_CHECK_NO_THROW( createSharedBiasTestEstimator(
                bodies, names, observerId, transmitterSettings, { bias, bias, bias }, { secondTransmitterSettings } ) );
    }
}

//! A vector observable's arc-wise relative partial is a matrix with a diagonal block in the active arc.
//! Check both arcs, the boundary, times before the first arc, and reset of the inactive block when revisiting an earlier time.
BOOST_AUTO_TEST_CASE( test_ArcWiseRelativeBiasPartialShape )
{
    const LinkEnds links = { { transmitter, LinkEndId( "A", "" ) }, { receiver, LinkEndId( "ObserverSatellite", "" ) } };
    const auto lookup = std::make_shared< interpolators::HuntingAlgorithmLookupScheme< double > >(
            std::vector< double >{ 0.0, 10.0, std::numeric_limits< double >::max( ) } );
    observation_partials::ObservationPartialWrtArcWiseRelativeBias< 2 > partial( angular_position, links, lookup, 1, 2 );
    const Eigen::Vector2d observable( 0.4, -0.2 );
    for( const double time : { -1.0, 5.0, 10.0, 15.0, 5.0 } )
    {
        const auto result = partial.calculatePartial( {}, { time - 1.0, time }, receiver, nullptr, observable ).front( );
        Eigen::Matrix< double, 2, 4 > expected = Eigen::Matrix< double, 2, 4 >::Zero( );
        if( time >= 0.0 )
        {
            expected.block< 2, 2 >( 0, time < 10.0 ? 0 : 2 ) = observable.asDiagonal( );
        }

        // Two observable rows, four bias columns, and no cross-component or inactive-arc contributions.
        BOOST_REQUIRE_EQUAL( result.first.cols( ), 4 );
        BOOST_CHECK_SMALL( ( result.first - expected ).norm( ), 1.0E-15 );
        BOOST_CHECK_EQUAL( result.second, time );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
