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

//! Exercise the public setup, propagation, simulation, and estimation interfaces together.
BOOST_AUTO_TEST_CASE( test_SharedObservationBiasEstimation )
{
    using namespace observation_models;
    using namespace estimatable_parameters;
    using namespace simulation_setup;

    // The last two cases use identical arcs but different event times for choosing an arc.
    const std::vector< ObservationBiasTypes > biasTypes = {
        constant_absolute_bias, constant_relative_bias, arc_wise_constant_absolute_bias, arc_wise_constant_absolute_bias
    };
    for( unsigned int testCase = 0; testCase < biasTypes.size( ); ++testCase )
    {
        BOOST_TEST_CONTEXT( "shared bias case " << testCase )
        {
            const auto biasType = biasTypes.at( testCase );
            const bool arcWise = biasType == arc_wise_constant_absolute_bias;
            const LinkEndType timeLinkEnd = testCase == 3 ? transmitter : receiver;
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

            BodyListSettings bodySettings( "SSB", "J2000" );
            for( const auto& name : { "Euclid", "OtherObserver" } )
            {
                bodySettings.addSettings( name );
                Eigen::Vector6d state = Eigen::Vector6d::Zero( );
                if( std::string( name ) == "OtherObserver" )
                {
                    state.head< 3 >( ) = Eigen::Vector3d( 2.0E7, -4.0E7, 1.0E7 );
                }
                bodySettings.at( name )->ephemerisSettings = constantEphemerisSettings( state, "SSB", "J2000" );
            }
            bodySettings.at( "Euclid" )->shapeModelSettings = sphericalBodyShapeSettings( 1.0 );
            bodySettings.at( "Euclid" )->rotationModelSettings =
                    constantRotationModelSettings( "J2000", "Euclid_fixed", Eigen::Quaterniond::Identity( ) );
            bodySettings.at( "Euclid" )
                    ->groundStationSettings.push_back( groundStationSettings( "camera", Eigen::Vector3d( 1.0, 0.0, 0.0 ) ) );
            SystemOfBodies bodies = createSystemOfBodies( bodySettings );
            std::vector< std::string > targets = { "A", "B", "C" };
            Eigen::VectorXd initialState( 18 );
            basic_astrodynamics::AccelerationMap accelerations;
            for( int i = 0; i < 3; ++i )
            {
                bodies.createEmptyBody( targets.at( i ) );
                bodies.at( targets.at( i ) )
                        ->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris<> >(
                                std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector6d > >( ),
                                "SSB",
                                "J2000" ) );
                initialState.segment< 6 >( 6 * i ) << 8.0E7 + i * 2.0E7, 5.0E7 - i * 1.0E7, 3.0E7 + i * 1.0E7, 100.0 + i * 20.0, -50.0,
                        30.0;
                accelerations[ targets.at( i ) ] = {};
            }
            auto propagator = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >(
                    std::vector< std::string >( 3, "SSB" ),
                    accelerations,
                    targets,
                    initialState,
                    0.0,
                    numerical_integrators::rungeKutta4Settings( 60.0 ),
                    std::make_shared< propagators::PropagationTimeTerminationSettings >( 7200.0 ) );

            std::vector< double > observationTimes;
            for( double time = 600.0; time <= 6600.0; time += 300.0 )
            {
                observationTimes.push_back( time );
            }
            std::vector< std::shared_ptr< ObservationModelSettings > > models;
            std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > simulations;
            std::vector< LinkEnds > angularLinks;
            auto addModel = [ & ]( ObservableType type, const LinkEnds& links, std::shared_ptr< ObservationBiasSettings > bias ) {
                models.push_back( std::make_shared< ObservationModelSettings >( type, links, nullptr, bias ) );
                simulations.push_back( std::make_shared< TabulatedObservationSimulationSettings< double > >(
                        type, links, observationTimes, type == position_observable ? observed_body : receiver ) );
                if( type == angular_position )
                {
                    angularLinks.push_back( links );
                }
            };
            auto makeAngularBias = [ & ]( ) -> std::shared_ptr< ObservationBiasSettings > {
                if( arcWise )
                {
                    return arcWiseAbsoluteBias( arcTimes, { trueAngularBias.head( 2 ), trueAngularBias.tail( 2 ) }, timeLinkEnd );
                }
                return biasType == constant_absolute_bias ? constantAbsoluteBias( trueAngularBias )
                                                          : constantRelativeBias( trueAngularBias );
            };
            const Eigen::Vector2d excludedBias( 7.0E-6, 9.0E-6 );
            for( const auto& target : targets )
            {
                const LinkEnds selected = { { transmitter, LinkEndId( target, "" ) }, { receiver, LinkEndId( "Euclid", "" ) } };
                // Also exercise finding a bias within a combined model.
                const auto secondaryBias = biasType == constant_relative_bias ? constantAbsoluteBias( Eigen::Vector2d( 4.0E-6, -2.0E-6 ) )
                                                                              : constantRelativeBias( Eigen::Vector2d::Zero( ) );
                std::vector< std::shared_ptr< ObservationBiasSettings > > combinedBiases = { makeAngularBias( ), secondaryBias };
                if( biasType == constant_relative_bias )
                {
                    // Shared relative partials must recompute at the biased event time, not the nominal observation time.
                    combinedBiases.push_back( constantTimeBias( 15.0, receiver ) );
                }
                addModel( angular_position, selected, multipleObservationBiasSettings( combinedBiases ) );
                addModel( one_way_range, selected, constantAbsoluteBias( Eigen::VectorXd::Constant( 1, 12.0 ) ) );
                addModel( angular_position,
                          { { transmitter, LinkEndId( target, "" ) }, { receiver, LinkEndId( "OtherObserver", "" ) } },
                          constantAbsoluteBias( excludedBias ) );
                // Independent position observations make all three target states identifiable.
                addModel( position_observable, { { observed_body, LinkEndId( target, "" ) } }, nullptr );
            }
            addModel( angular_position,
                      { { transmitter, LinkEndId( "A", "" ) }, { receiver, LinkEndId( "Euclid", "camera" ) } },
                      constantAbsoluteBias( excludedBias ) );
            addModel( angular_position,
                      { { transmitter, LinkEndId( "Euclid", "" ) }, { receiver, LinkEndId( "A", "" ) } },
                      constantAbsoluteBias( excludedBias ) );
            BOOST_REQUIRE_EQUAL( angularLinks.size( ), 8 );

            std::vector< std::shared_ptr< EstimatableParameterSettings > > settings;
            for( int i = 0; i < 3; ++i )
            {
                settings.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
                        targets.at( i ), initialState.segment< 6 >( 6 * i ), "SSB" ) );
            }
            settings.push_back( sharedObservationBias( biasType,
                                                       angular_position,
                                                       receiver,
                                                       LinkEndId( "Euclid", "" ),
                                                       arcTimes,
                                                       testCase == 3 ? transmitter : unidentified_link_end ) );
            settings.push_back( sharedObservationBias( constant_absolute_bias, one_way_range, receiver, LinkEndId( "Euclid", "" ) ) );
            auto parameters = createParametersToEstimate< double, double >( settings, bodies, propagator );
            auto angularParameter =
                    std::dynamic_pointer_cast< SharedObservationBiasParameter >( parameters->getVectorParameters( ).at( 18 ) );
            auto rangeParameter = std::dynamic_pointer_cast< SharedObservationBiasParameter >(
                    parameters->getVectorParameters( ).at( 18 + trueAngularBias.size( ) ) );
            BOOST_REQUIRE( angularParameter != nullptr );
            BOOST_REQUIRE( rangeParameter != nullptr );
            OrbitDeterminationManager< double, double > manager( bodies, parameters, models, propagator );
            BOOST_CHECK_EQUAL( angularParameter->getMembers( ).size( ), 3 );
            BOOST_CHECK_EQUAL( rangeParameter->getMembers( ).size( ), 3 );
            BOOST_CHECK_EQUAL( parameters->getParameterSetSize( ), 18 + trueAngularBias.size( ) + 1 );
            const Eigen::VectorXd truth = parameters->getFullParameterValues< double >( );
            auto observations = simulateObservations< double, double >( simulations, manager.getObservationSimulators( ), bodies );

            // Check full manager partials against finite differences on every angular link.
            const auto angularManager = manager.getObservationManagers( ).at( angular_position );
            const std::vector< double > checkTimes = { 1800.0, 3599.0, 3600.0, 3601.0, 5400.0 };
            for( const auto& links : angularLinks )
            {
                Eigen::VectorXd values;
                Eigen::MatrixXd partials;
                angularManager->computeObservationsWithPartials( checkTimes, links, receiver, nullptr, values, partials );
                const bool selected = angularParameter->matches( links, angular_position );
                if( !selected )
                {
                    BOOST_CHECK_SMALL( partials.middleCols( 18, trueAngularBias.size( ) ).norm( ), 1.0E-30 );
                }
                for( int column = 0; column < trueAngularBias.size( ); ++column )
                {
                    const double step = 1.0E-6;
                    Eigen::VectorXd perturbed = trueAngularBias;
                    perturbed( column ) += step;
                    angularParameter->setParameterValue( perturbed );
                    Eigen::VectorXd plus, minus;
                    Eigen::MatrixXd unused;
                    angularManager->computeObservationsWithPartials( checkTimes, links, receiver, nullptr, plus, unused, true, false );
                    perturbed( column ) -= 2.0 * step;
                    angularParameter->setParameterValue( perturbed );
                    angularManager->computeObservationsWithPartials( checkTimes, links, receiver, nullptr, minus, unused, true, false );
                    angularParameter->setParameterValue( trueAngularBias );
                    BOOST_CHECK_SMALL( ( ( plus - minus ) / ( 2.0 * step ) - partials.col( 18 + column ) ).norm( ), 1.0E-8 );
                }
                if( arcWise && selected )
                {
                    const int activeColumn = timeLinkEnd == receiver ? 20 : 18;
                    BOOST_CHECK_SMALL( ( partials.block< 2, 2 >( 4, activeColumn ) - Eigen::Matrix2d::Identity( ) ).norm( ), 1.0E-15 );
                }
            }
            Eigen::VectorXd rangeValues;
            Eigen::MatrixXd rangePartials;
            manager.getObservationManagers( )
                    .at( one_way_range )
                    ->computeObservationsWithPartials( checkTimes, angularLinks.front( ), receiver, nullptr, rangeValues, rangePartials );
            BOOST_CHECK_SMALL( rangePartials.middleCols( 18, trueAngularBias.size( ) ).norm( ), 1.0E-30 );
            BOOST_CHECK_SMALL( ( rangePartials.col( 18 + trueAngularBias.size( ) ).array( ) - 1.0 ).matrix( ).norm( ), 1.0E-15 );

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
            BOOST_REQUIRE( !output->exceptionDuringInversion_ );
            BOOST_CHECK_SMALL( ( output->parameterEstimate_.head( 18 ) - truth.head( 18 ) ).norm( ), 1.0E-5 );
            BOOST_CHECK_SMALL( ( output->parameterEstimate_.segment( 18, trueAngularBias.size( ) ) - trueAngularBias ).norm( ), 1.0E-11 );
            BOOST_CHECK_SMALL( std::abs( output->parameterEstimate_.tail( 1 )( 0 ) - 12.0 ), 1.0E-6 );
            BOOST_CHECK_SMALL( output->residuals_.cwiseAbs( ).maxCoeff( ), 1.0E-5 );
            for( const auto& member : angularParameter->getMembers( ) )
            {
                BOOST_CHECK_SMALL( ( member.second->getParameterValue( ) - trueAngularBias ).norm( ), 1.0E-11 );
            }
            auto simulator = std::dynamic_pointer_cast< ObservationSimulator< 2 > >( angularManager->getObservationSimulator( ) );
            for( const auto& links : angularLinks )
            {
                if( !angularParameter->matches( links, angular_position ) )
                {
                    const auto bias = std::dynamic_pointer_cast< ConstantObservationBias< 2 > >(
                            simulator->getObservationModel( links )->getObservationBiasCalculator( ) );
                    BOOST_REQUIRE( bias != nullptr );
                    BOOST_CHECK_SMALL( ( bias->getConstantObservationBias( ) - excludedBias ).norm( ), 1.0E-30 );
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( test_SharedObservationBiasValidation )
{
    using namespace observation_models;
    using namespace estimatable_parameters;
    using namespace simulation_setup;

    const LinkEndId observerId( "Euclid", "" );
    const auto constantSettings = sharedObservationBias( constant_absolute_bias, angular_position, receiver, observerId );
    BOOST_CHECK_THROW( sharedObservationBias( constant_time_bias, angular_position, receiver, observerId ), std::runtime_error );
    BOOST_CHECK_THROW( sharedObservationBias( arc_wise_constant_relative_bias, angular_position, receiver, observerId ),
                       std::runtime_error );
    BOOST_CHECK_THROW( sharedObservationBias( arc_wise_constant_absolute_bias, angular_position, receiver, observerId ),
                       std::runtime_error );
    BOOST_CHECK_THROW( sharedObservationBias( constant_absolute_bias, angular_position, receiver, observerId, { 0.0 } ),
                       std::runtime_error );
    BOOST_CHECK_THROW( sharedObservationBias( arc_wise_constant_absolute_bias, angular_position, receiver, observerId, { 2.0, 1.0 } ),
                       std::runtime_error );
    BOOST_CHECK_THROW( sharedObservationBias( arc_wise_constant_absolute_bias, angular_position, receiver, observerId, { 0.0, 0.0 } ),
                       std::runtime_error );
    BOOST_CHECK_THROW( sharedObservationBias( arc_wise_constant_absolute_bias, angular_position, receiver, observerId, { TUDAT_NAN } ),
                       std::runtime_error );

    BodyListSettings bodySettings( "SSB", "J2000" );
    const std::vector< std::string > names = { "Euclid", "A", "B", "C" };
    for( unsigned int i = 0; i < names.size( ); ++i )
    {
        bodySettings.addSettings( names.at( i ) );
        Eigen::Vector6d state = Eigen::Vector6d::Zero( );
        state.head< 3 >( ) = Eigen::Vector3d( i * 1.0E8, i * 2.0E8, i * 3.0E8 );
        bodySettings.at( names.at( i ) )->ephemerisSettings = constantEphemerisSettings( state, "SSB", "J2000" );
    }
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    auto makeModels = [ & ]( const std::vector< std::shared_ptr< ObservationBiasSettings > >& biases ) {
        std::vector< std::shared_ptr< ObservationModelSettings > > models;
        for( unsigned int i = 0; i < biases.size( ); ++i )
        {
            LinkEnds links = { { transmitter, LinkEndId( names.at( i + 1 ), "" ) }, { receiver, observerId } };
            models.push_back( std::make_shared< ObservationModelSettings >( angular_position, links, nullptr, biases.at( i ) ) );
        }
        return models;
    };
    const std::shared_ptr< propagators::PropagatorSettings< double > > noPropagator;
    auto makeManager = [ & ]( const std::shared_ptr< EstimatableParameterSettings >& settings,
                              const std::vector< std::shared_ptr< ObservationBiasSettings > >& biases ) {
        auto parameters = createParametersToEstimate< double, double >( { settings }, bodies );
        return std::make_shared< OrbitDeterminationManager< double, double > >( bodies, parameters, makeModels( biases ), noPropagator );
    };
    const auto zeroBias = constantAbsoluteBias( Eigen::Vector2d::Zero( ) );
    const auto otherBias = constantAbsoluteBias( Eigen::Vector2d::Ones( ) );
    BOOST_CHECK_THROW( makeManager( constantSettings, { zeroBias, nullptr, zeroBias } ), std::runtime_error );
    BOOST_CHECK_THROW( makeManager( constantSettings, { zeroBias, otherBias, zeroBias } ), std::runtime_error );
    BOOST_CHECK_THROW( makeManager( constantSettings, { multipleObservationBiasSettings( { zeroBias, otherBias } ) } ),
                       std::runtime_error );
    BOOST_CHECK_THROW( makeManager( sharedObservationBias( constant_absolute_bias, angular_position, transmitter, observerId ),
                                    { zeroBias, zeroBias, zeroBias } ),
                       std::runtime_error );
    BOOST_CHECK_THROW( makeManager( sharedObservationBias( constant_absolute_bias, one_way_range, receiver, observerId ),
                                    { zeroBias, zeroBias, zeroBias } ),
                       std::runtime_error );
    const auto arcSettings =
            sharedObservationBias( arc_wise_constant_absolute_bias, angular_position, receiver, observerId, { 0.0, 10.0 } );
    const std::vector< Eigen::VectorXd > arcValues = { Eigen::Vector2d::Zero( ), Eigen::Vector2d::Ones( ) };
    const auto goodArc = arcWiseAbsoluteBias( { 0.0, 10.0 }, arcValues, receiver );
    BOOST_CHECK_THROW( makeManager( arcSettings, { goodArc, arcWiseAbsoluteBias( { 0.0, 11.0 }, arcValues, receiver ) } ),
                       std::runtime_error );
    BOOST_CHECK_THROW( makeManager( arcSettings, { goodArc, arcWiseAbsoluteBias( { 0.0, 10.0 }, arcValues, transmitter ) } ),
                       std::runtime_error );

    // A deliberate pre-closure assignment initializes every member, even if their settings differ.
    auto parameters = createParametersToEstimate< double, double >( { constantSettings }, bodies );
    const auto shared = std::dynamic_pointer_cast< SharedObservationBiasParameter >( parameters->getVectorParameters( ).at( 0 ) );
    BOOST_REQUIRE( shared != nullptr );
    const Eigen::Vector2d desired( 1.0E-4, -2.0E-4 );
    parameters->resetParameterValues( Eigen::VectorXd( desired ) );
    OrbitDeterminationManager< double, double > manager(
            bodies, parameters, makeModels( { zeroBias, otherBias, zeroBias } ), noPropagator );
    BOOST_CHECK_EQUAL( shared->getMembers( ).size( ), 3 );
    BOOST_CHECK_SMALL( ( shared->getParameterValue( ) - desired ).norm( ), 1.0E-30 );
    for( const auto& member : shared->getMembers( ) )
    {
        BOOST_CHECK_SMALL( ( member.second->getParameterValue( ) - desired ).norm( ), 1.0E-30 );
    }
    BOOST_CHECK_THROW( shared->setParameterValue( Eigen::Vector3d::Zero( ) ), std::runtime_error );
    shared->getMembers( ).begin( )->second->setParameterValue( Eigen::Vector2d::Zero( ) );
    BOOST_CHECK_THROW( shared->getParameterValue( ), std::runtime_error );
    shared->setParameterValue( desired );
    auto simulator = std::dynamic_pointer_cast< ObservationSimulator< 2 > >(
            manager.getObservationManagers( ).at( angular_position )->getObservationSimulator( ) );
    performObservationParameterEstimationClosure( simulator, parameters );
    BOOST_CHECK_EQUAL( shared->getMembers( ).size( ), 3 );
    BOOST_CHECK_SMALL( ( shared->getParameterValue( ) - desired ).norm( ), 1.0E-30 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
