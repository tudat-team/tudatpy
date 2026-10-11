#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include <limits>

#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/astro/ephemerides/customEphemeris.h"
#include "tudat/astro/gravitation/gravityFieldModel.h"
#include "tudat/simulation/estimation_setup/createEstimatableParametersFactory.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"
#include "tudat/simulation/propagation_setup/createAccelerationModels.h"

namespace tudat
{
namespace unit_tests
{
using namespace simulation_setup;
using namespace propagators;
using namespace observation_models;
using namespace estimatable_parameters;
using namespace numerical_integrators;
using LongState = Eigen::Matrix< long double, 6, 1 >;

namespace
{
const std::string centerName = "Jupiter_barycenter";
const double sunMu = 1.32712440018E20;
const double jupiterMu = 1.26686534E17;
const double day = 86400.0;
// Some platforms represent long double with double precision. Preserve the
// physical checks without requiring sub-ulp accuracy at heliocentric distances.
const double positionRoundoff =
        std::max( 2.0E-7, static_cast< double >( 2.0L * std::numeric_limits< long double >::epsilon( ) * 7.785E11L ) );
const double velocityRoundoff =
        std::max( 2.0E-12, static_cast< double >( 2.0L * std::numeric_limits< long double >::epsilon( ) * 3.0E4L ) );
const double offsetAccelerationRoundoff =
        std::max( 1.0E-13, static_cast< double >( 8.0L * std::numeric_limits< long double >::epsilon( ) * 3.0E4L / 10.0L ) );
const std::vector< std::string > members = { "Jupiter", "Io", "Europa", "Ganymede", "Callisto" };
const std::vector< double > moonMu = { 5.9599E12, 3.2027E12, 9.8878E12, 7.1793E12 };
const std::vector< double > moonRadius = { 4.217E8, 6.711E8, 1.0704E9, 1.8827E9 };

// Analytic, velocity-consistent ephemerides keep this regression independent of kernel downloads.
LongState circularState( const long double time,
                         const long double radius,
                         const long double mu,
                         const long double phase,
                         const long double inclination )
{
    const long double n = std::sqrt( mu / ( radius * radius * radius ) );
    const long double angle = n * time + phase;
    LongState state;
    state << radius * std::cos( angle ), radius * std::sin( angle ) * std::cos( inclination ),
            radius * std::sin( angle ) * std::sin( inclination ), -radius * n * std::sin( angle ),
            radius * n * std::cos( angle ) * std::cos( inclination ), radius * n * std::cos( angle ) * std::sin( inclination );
    return state;
}

LongState moonState( const int moon, const long double time )
{
    return circularState( time, moonRadius.at( moon ), jupiterMu, 0.7L * moon, 0.01L * ( moon + 1 ) );
}

LongState nominalCenter( const long double time )
{
    return circularState( time, 7.785E11L, sunMu, 0.4L, 0.02L );
}

SystemOfBodies makeEnvironment( const bool hybrid = false )
{
    SystemOfBodies bodies;
    for( const auto& name : { "Sun", "Earth", "Jupiter", "Io", "Europa", "Ganymede", "Callisto", "Juno" } )
    {
        bodies.createEmptyBody( name, false );
    }
    bodies.at( "Sun" )->setEphemeris( std::make_shared< ephemerides::CustomEphemeris< long double, long double > >(
            []( const long double ) -> LongState { return LongState::Zero( ); } ) );
    bodies.at( "Sun" )->setGravityFieldModel( std::make_shared< gravitation::GravityFieldModel >( sunMu ) );
    bodies.at( "Earth" )->setEphemeris( std::make_shared< ephemerides::CustomEphemeris< long double, long double > >(
            []( const long double time ) -> LongState { return circularState( time, 1.496E11L, sunMu, 1.3L, 0.0L ); } ) );
    bodies.at( "Jupiter" )
            ->setEphemeris( std::make_shared< ephemerides::CustomEphemeris< long double, long double > >(
                    []( const long double time ) -> LongState {
                        LongState weightedOffset = LongState::Zero( );
                        long double totalMu = jupiterMu;
                        for( unsigned int i = 0; i < moonMu.size( ); ++i )
                        {
                            weightedOffset += moonMu[ i ] * moonState( i, time );
                            totalMu += moonMu[ i ];
                        }
                        return nominalCenter( time ) - weightedOffset / totalMu;
                    } ) );
    bodies.at( "Jupiter" )->setGravityFieldModel( std::make_shared< gravitation::GravityFieldModel >( jupiterMu ) );
    for( unsigned int i = 0; i < moonMu.size( ); ++i )
    {
        bodies.at( members[ i + 1 ] )
                ->setEphemeris( std::make_shared< ephemerides::CustomEphemeris< long double, long double > >(
                        [ i ]( const long double time ) -> LongState { return moonState( i, time ); }, "Jupiter" ) );
        bodies.at( members[ i + 1 ] )->setGravityFieldModel( std::make_shared< gravitation::GravityFieldModel >( moonMu[ i ] ) );
    }
    addEmptyTabulatedEphemeris< double, double >( bodies, "Juno", "Jupiter", hybrid );
    bodies.processBodyFrameDefinitions< long double, Time >( );
    bodies.addBarycenter( centerName, members );
    return bodies;
}

Eigen::Vector6d flybyInitialState( const int arc )
{
    // Juno-like eccentric, nearly polar Jupiter orbit, starting before perijove.
    Eigen::Vector6d elements;
    elements << 1.5E9, 0.95, 1.55, 0.4 + 0.15 * arc, 0.2 + 0.3 * arc, -2.0;
    return orbital_element_conversions::convertKeplerianToCartesianElements( elements, jupiterMu );
}

std::shared_ptr< TranslationalStatePropagatorSettings<> > translationSettings( const std::vector< std::string >& centralBodies,
                                                                               const basic_astrodynamics::AccelerationMap& accelerations,
                                                                               const std::vector< std::string >& propagated,
                                                                               const Eigen::VectorXd& initialState,
                                                                               const double start,
                                                                               const double end,
                                                                               const double step )
{
    return std::make_shared< TranslationalStatePropagatorSettings<> >( centralBodies,
                                                                       accelerations,
                                                                       propagated,
                                                                       initialState,
                                                                       start,
                                                                       rungeKuttaFixedStepSettings( step, rungeKuttaFehlberg78 ),
                                                                       std::make_shared< PropagationTimeTerminationSettings >( end ) );
}

void runJunoEstimation( const bool hybrid )
{
    auto bodies = makeEnvironment( hybrid );
    const int numberOfArcs = hybrid ? 3 : 1;
    std::vector< double > starts, ends, observationTimes;
    for( int arc = 0; arc < numberOfArcs; ++arc )
    {
        const double epoch = arc * 90.0 * day;  // Three months are defined as 90 Julian days in this fixture.
        starts.push_back( epoch - 10800.0 );
        ends.push_back( epoch + 10800.0 );
        for( int sample = 0; sample <= 30; ++sample )
        {
            observationTimes.push_back( epoch - 7200.0 + 480.0 * sample );
        }
    }

    SelectedAccelerationMap selected;
    selected[ centerName ][ "Sun" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );
    auto barycenterAccelerations = createAccelerationModelsMap( bodies, selected, { centerName }, { "SSB" } );
    selected.clear( );
    for( const auto& source : members )
    {
        selected[ "Juno" ][ source ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );
    }
    selected[ "Juno" ][ "Sun" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );
    auto junoAccelerations = createAccelerationModelsMap( bodies, selected, { "Juno" }, { "Jupiter" } );
    auto accelerations = barycenterAccelerations;
    accelerations.insert( junoAccelerations.begin( ), junoAccelerations.end( ) );

    std::shared_ptr< PropagatorSettings<> > propagation;
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterSettings;
    if( hybrid )
    {
        const double start = starts.front( ) - day;
        auto single = translationSettings( { "SSB" },
                                           barycenterAccelerations,
                                           { centerName },
                                           nominalCenter( start ).cast< double >( ),
                                           start,
                                           ends.back( ) + day,
                                           7200.0 );
        std::vector< std::shared_ptr< SingleArcPropagatorSettings<> > > arcs;
        for( int arc = 0; arc < numberOfArcs; ++arc )
        {
            arcs.push_back( translationSettings(
                    { "Jupiter" }, junoAccelerations, { "Juno" }, flybyInitialState( arc ), starts[ arc ], ends[ arc ], 30.0 ) );
        }
        auto settings = std::make_shared< HybridArcPropagatorSettings<> >(
                single, std::make_shared< MultiArcPropagatorSettings<> >( arcs, false ) );
        parameterSettings = getInitialHybridArcParameterSettings<>( settings, bodies, starts );
        propagation = settings;
    }
    else
    {
        Eigen::VectorXd initial( 12 );
        initial << nominalCenter( starts.front( ) ).cast< double >( ), flybyInitialState( 0 );
        auto settings = translationSettings(
                { "SSB", "Jupiter" }, accelerations, { centerName, "Juno" }, initial, starts.front( ), ends.front( ), 30.0 );
        parameterSettings = getInitialStateParameterSettings< double >( settings, bodies );
        propagation = settings;
    }
    auto parameters = createParametersToEstimate<>( parameterSettings, bodies );
    const Eigen::VectorXd truth = parameters->getFullParameterValues< double >( );
    BOOST_REQUIRE_EQUAL( truth.size( ), 6 * ( 1 + numberOfArcs ) );

    if( !hybrid )
    {
        const double time = starts.front( );
        auto center = bodies.getBarycenter( centerName );
        center->beginStateUpdate( Time( time ) );
        center->setCurrentState< double >( truth.head< 6 >( ) );
        bodies.at( "Sun" )->setState( Eigen::Vector6d::Zero( ) );
        for( const auto& member : members )
        {
            bodies.at( member )->setState( bodies.at( member )->getEphemeris( )->getCartesianState( time ) );
        }
        bodies.at( "Juno" )->setState( bodies.at( "Jupiter" )->getState( ) + truth.tail< 6 >( ) );
        Eigen::Vector3d relativeAcceleration = Eigen::Vector3d::Zero( );
        for( const auto& source : junoAccelerations.at( "Juno" ) )
        {
            for( const auto& model : source.second )
            {
                model->resetCurrentTime( );
                model->updateMembers( time );
                relativeAcceleration += model->getAcceleration( );
            }
        }
        Eigen::Vector3d directAcceleration = Eigen::Vector3d::Zero( );
        for( const auto& source : { "Sun", "Jupiter", "Io", "Europa", "Ganymede", "Callisto" } )
        {
            const Eigen::Vector3d displacement = bodies.at( "Juno" )->getPosition( ) - bodies.at( source )->getPosition( );
            directAcceleration -= bodies.at( source )->getGravitationalParameter( ) * displacement / std::pow( displacement.norm( ), 3 );
        }
        const Eigen::Vector3d positionB = truth.head< 3 >( );
        const Eigen::Vector3d accelerationB = -sunMu * positionB / std::pow( positionB.norm( ), 3 );
        BOOST_CHECK_SMALL(
                ( relativeAcceleration + accelerationB + center->getMemberOffsetAcceleration( "Jupiter", time ) - directAcceleration )
                        .norm( ),
                1.0E-13 );

        auto partials = createAccelerationPartialsMap( accelerations, bodies, parameters );
        VariationalEquations equations( { { translational_state, partials } },
                                        parameters,
                                        { { translational_state, 0 } },
                                        -1,
                                        {},
                                        bodies.getTranslationalStateAliases( ) );
        equations.updatePartials< double >( time, { { translational_state, truth } } );
        const Eigen::MatrixXd identity = Eigen::MatrixXd::Identity( 12, 12 );
        Eigen::MatrixXd jacobian = Eigen::MatrixXd::Zero( 12, 12 );
        equations.getBodyInitialStatePartialMatrix< double >( identity, jacobian.block( 0, 0, 12, 12 ) );
        const Eigen::Matrix3d expected = acceleration_partials::calculatePartialOfPointMassGravityWrtPositionOfAcceleratedBody(
                                                 bodies.at( "Juno" )->getPosition( ), sunMu ) -
                acceleration_partials::calculatePartialOfPointMassGravityWrtPositionOfAcceleratedBody( positionB, sunMu );
        // All Jupiter/moon contributions cancel. Only the external-field difference remains.
        BOOST_CHECK_SMALL( ( jacobian.block< 3, 3 >( 9, 0 ) - expected ).norm( ), 1.0E-20 );
        BOOST_CHECK_SMALL( jacobian.block( 6, 3, 3, 3 ).norm( ), 1.0E-15 );
        center->clearCurrentState( );
    }

    LinkEnds positionLink, rangeLink, memberLink;
    positionLink[ observed_body ] = LinkEndId( "Juno" );
    rangeLink[ transmitter ] = LinkEndId( "Earth" );
    rangeLink[ receiver ] = LinkEndId( "Juno" );
    memberLink[ observed_body ] = LinkEndId( "Ganymede" );
    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettings = {
        std::make_shared< ObservationModelSettings >( position_observable, positionLink ),
        std::make_shared< ObservationModelSettings >( one_way_range, rangeLink ),
        std::make_shared< ObservationModelSettings >( position_observable, memberLink )
    };
    OrbitDeterminationManager<> manager( bodies, parameters, observationSettings, propagation );

    // Check the defining identity after installing propagated ephemerides, also in
    // the long gap between spacecraft arcs, where only the barycenter is propagated.
    auto center = bodies.getBarycenter( centerName );
    for( const double time : { observationTimes.front( ), hybrid ? 45.0 * day : 0.0, observationTimes.back( ) } )
    {
        const LongState stateB = center->getEphemeris( )->getCartesianLongState( time );
        LongState weightedOffset = LongState::Zero( );
        for( unsigned int i = 0; i < members.size( ); ++i )
        {
            weightedOffset +=
                    center->getWeights( )[ i ] * ( bodies.at( members[ i ] )->getEphemeris( )->getCartesianLongState( time ) - stateB );
        }
        BOOST_CHECK_SMALL( static_cast< double >( weightedOffset.head< 3 >( ).norm( ) ), positionRoundoff );
        BOOST_CHECK_SMALL( static_cast< double >( weightedOffset.tail< 3 >( ).norm( ) ), velocityRoundoff );
    }

    Eigen::VectorXd position, range;
    Eigen::MatrixXd positionPartials, rangePartials;
    manager.computePartialsAndObservations< 3 >(
            positionLink, nullptr, position_observable, observed_body, observationTimes, position, positionPartials );
    manager.computePartialsAndObservations< 1 >( rangeLink, nullptr, one_way_range, receiver, observationTimes, range, rangePartials );
    Eigen::VectorXd memberPosition;
    Eigen::MatrixXd memberPartials;
    manager.computePartialsAndObservations< 3 >(
            memberLink, nullptr, position_observable, observed_body, observationTimes, memberPosition, memberPartials );
    BOOST_CHECK_SMALL( memberPartials.rightCols( truth.size( ) - 6 ).norm( ), 1.0E-15 );
    BOOST_CHECK( positionPartials.leftCols( 6 ).norm( ) > 1.0 );
    BOOST_CHECK( rangePartials.leftCols( 6 ).norm( ) > 1.0 );

    // End-to-end finite differences cover every barycenter and local-arc initial component,
    // including reconstruction of member states and the light-time observation chain rule.
    for( int column = 0; column < truth.size( ); ++column )
    {
        // Larger common-state perturbations resolve weak out-of-plane range
        // sensitivities above double-precision roundoff at interplanetary distances.
        const double step = column % 6 < 3 ? ( column < 6 ? 1.0E5 : 1.0E3 ) : ( column < 6 ? 0.1 : 0.01 );
        Eigen::VectorXd plusPosition, minusPosition, plusRange, minusRange, plusMember, minusMember;
        Eigen::MatrixXd ignored;
        Eigen::VectorXd perturbed = truth;
        perturbed[ column ] += step;
        manager.resetParameterEstimate( perturbed, false );
        manager.computePartialsAndObservations< 3 >(
                positionLink, nullptr, position_observable, observed_body, observationTimes, plusPosition, ignored );
        manager.computePartialsAndObservations< 1 >( rangeLink, nullptr, one_way_range, receiver, observationTimes, plusRange, ignored );
        if( column < 6 )
        {
            manager.computePartialsAndObservations< 3 >(
                    memberLink, nullptr, position_observable, observed_body, observationTimes, plusMember, ignored );
        }
        perturbed[ column ] -= 2.0 * step;
        manager.resetParameterEstimate( perturbed, false );
        manager.computePartialsAndObservations< 3 >(
                positionLink, nullptr, position_observable, observed_body, observationTimes, minusPosition, ignored );
        manager.computePartialsAndObservations< 1 >( rangeLink, nullptr, one_way_range, receiver, observationTimes, minusRange, ignored );
        if( column < 6 )
        {
            manager.computePartialsAndObservations< 3 >(
                    memberLink, nullptr, position_observable, observed_body, observationTimes, minusMember, ignored );
            const Eigen::VectorXd numericalMember = ( plusMember - minusMember ) / ( 2.0 * step );
            const double memberError = ( numericalMember - memberPartials.col( column ) ).norm( ) / memberPartials.col( column ).norm( );
            BOOST_CHECK_MESSAGE( memberError < 2.0E-5,
                                 "Member partial column " << column << ": relative error " << memberError << ", numerical norm "
                                                          << numericalMember.norm( ) << ", analytic norm "
                                                          << memberPartials.col( column ).norm( ) );
        }
        const Eigen::VectorXd numericalPosition = ( plusPosition - minusPosition ) / ( 2.0 * step );
        const Eigen::VectorXd numericalRange = ( plusRange - minusRange ) / ( 2.0 * step );
        BOOST_CHECK_MESSAGE(
                ( numericalPosition - positionPartials.col( column ) ).norm( ) / positionPartials.col( column ).norm( ) < 2.0E-5,
                "Cartesian partial column " << column );
        const double rangeError = ( numericalRange - rangePartials.col( column ) ).norm( ) / rangePartials.col( column ).norm( );
        BOOST_TEST_MESSAGE( "Range finite-difference relative error, column " << column << ": " << rangeError );
        BOOST_CHECK_MESSAGE( rangeError < 2.0E-5,
                             "Range partial column " << column << ": relative error " << rangeError << ", numerical norm "
                                                     << numericalRange.norm( ) << ", analytic norm "
                                                     << rangePartials.col( column ).norm( ) );
    }
    if( hybrid )
    {
        for( int arc = 0; arc < numberOfArcs; ++arc )
        {
            BOOST_CHECK( positionPartials.block( arc * 93, 0, 93, 6 ).norm( ) > 1.0 );
            for( int otherArc = 0; otherArc < numberOfArcs; ++otherArc )
            {
                if( arc != otherArc )
                {
                    BOOST_CHECK_SMALL( positionPartials.block( arc * 93, 6 + otherArc * 6, 93, 6 ).norm( ), 1.0E-15 );
                }
            }
        }
    }

    manager.resetParameterEstimate( truth );
    std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > simulationSettings = {
        std::make_shared< TabulatedObservationSimulationSettings<> >( position_observable, positionLink, observationTimes, observed_body ),
        std::make_shared< TabulatedObservationSimulationSettings<> >( one_way_range, rangeLink, observationTimes, receiver )
    };
    auto observations = simulateObservations<>( simulationSettings, manager.getObservationSimulators( ), bodies );
    Eigen::VectorXd initialEstimate = truth;
    for( int i = 0; i < truth.size( ); ++i )
    {
        initialEstimate[ i ] += ( i % 6 < 3 ? 100.0 : 0.001 ) * ( i % 2 ? -1.0 : 1.0 );
    }
    parameters->resetParameterValues( initialEstimate );
    auto input = std::make_shared< EstimationInput<> >( observations );
    input->defineEstimationSettings( true, true, true, false );
    input->setConvergenceChecker( std::make_shared< EstimationConvergenceChecker >( 5 ) );
    const auto output = manager.estimateParameters( input );
    const Eigen::VectorXd error = output->parameterEstimate_ - truth;
    for( int i = 0; i < error.size( ); ++i )
    {
        BOOST_CHECK_SMALL( error[ i ], i % 6 < 3 ? 2.0 : 1.0E-3 );
    }
    BOOST_CHECK_SMALL( output->residuals_.norm( ) / std::sqrt( output->residuals_.size( ) ), 0.01 );
}
}  // namespace

BOOST_AUTO_TEST_SUITE( test_barycenter_estimation )

BOOST_AUTO_TEST_CASE( registrationAndPrescribedStates )
{
    auto bodies = makeEnvironment( );
    auto center = bodies.getBarycenter( centerName );
    BOOST_CHECK( !bodies.doesBodyExist( centerName ) );
    BOOST_CHECK( bodies.doesBarycenterExist( centerName ) );
    BOOST_CHECK_EQUAL( bodies.getNumberOfBodies( ), 8 );
    BOOST_CHECK_THROW( bodies.at( centerName ), std::runtime_error );
    BOOST_CHECK_THROW( bodies.addBarycenter( "overlap", { "Jupiter", "Io" } ), std::runtime_error );
    BOOST_CHECK_THROW( bodies.createEmptyBody( "Jupiter" ), std::runtime_error );
    BOOST_CHECK_THROW( bodies.addBarycenter( "invalid", { "Sun" }, { -1.0 } ), std::runtime_error );
    const auto weights = center->getWeights( );
    bodies.at( "Jupiter" )->getGravityFieldModel( )->resetGravitationalParameter( 1.01 * jupiterMu );
    BOOST_CHECK( center->getWeights( ) == weights );

    // A member origin must resolve from this stage's integrated vector, even when
    // the spacecraft precedes its barycenter and the environment holds an old state.
    auto origins = createCentralBodyData< long double, Time >( { "Jupiter", "SSB" }, { "Juno", centerName }, bodies );
    Eigen::Matrix< long double, Eigen::Dynamic, 1 > reversedState( 12 );
    reversedState << flybyInitialState( 0 ).cast< long double >( ), nominalCenter( 0.0 );
    reversedState[ 6 ] += 1.0E6L;
    center->beginStateUpdate( Time( 0.0 ) );
    center->setCurrentState( nominalCenter( 0.0 ) );
    std::vector< LongState > originStates;
    origins->getReferenceFrameOriginInertialStates( reversedState, Time( 0.0 ), originStates );
    BOOST_CHECK_SMALL(
            static_cast< double >(
                    ( originStates.at( 0 ) - reversedState.tail< 6 >( ) - center->getMemberOffset( "Jupiter", Time( 0.0 ) ) ).norm( ) ),
            positionRoundoff );
    center->setCurrentState< long double >( reversedState.tail< 6 >( ) );
    const LongState prescribed = bodies.getTranslationalStateFromEphemeris< long double, Time >( "Jupiter", Time( 0.0 ) );
    BOOST_CHECK_SMALL(
            static_cast< double >( ( prescribed - nominalCenter( 0.0 ) - center->getMemberOffset( "Jupiter", Time( 0.0 ) ) ).norm( ) ),
            positionRoundoff );
    center->clearCurrentState( );

    for( const double time : { 0.0, 90.0 * day, 180.0 * day } )
    {
        LongState weightedOffset = LongState::Zero( );
        for( unsigned int i = 0; i < members.size( ); ++i )
        {
            weightedOffset += weights[ i ] * center->getMemberOffset( members[ i ], Time( time ) );
        }
        BOOST_CHECK_SMALL( static_cast< double >( weightedOffset.head< 3 >( ).norm( ) ), positionRoundoff );
        BOOST_CHECK_SMALL( static_cast< double >( weightedOffset.tail< 3 >( ).norm( ) ), velocityRoundoff );
        LongState expectedOffsetAcceleration = LongState::Zero( );
        for( unsigned int moon = 0; moon < moonMu.size( ); ++moon )
        {
            expectedOffsetAcceleration += weights[ moon + 1 ] * jupiterMu / std::pow( moonRadius[ moon ], 3 ) * moonState( moon, time );
        }
        BOOST_CHECK_SMALL(
                ( center->getMemberOffsetAcceleration( "Jupiter", time ) - expectedOffsetAcceleration.head< 3 >( ).cast< double >( ) )
                        .norm( ),
                offsetAccelerationRoundoff );
        LongState displaced = nominalCenter( time );
        displaced.head< 3 >( ) += Eigen::Vector3d( 1000.0, -2000.0, 3000.0 ).cast< long double >( );
        displaced.tail< 3 >( ) += Eigen::Vector3d( 1.0, -2.0, 3.0 ).cast< long double >( );
        center->beginStateUpdate( Time( time ) );
        center->setCurrentState( displaced );
        for( const auto& member : members )
        {
            const LongState actual =
                    bodies.at( member )->getEphemeris( )->getTemplatedStateFromEphemeris< long double, Time >( Time( time ) );
            const LongState expected = displaced + center->getMemberOffset( member, Time( time ) );
            BOOST_CHECK_SMALL( static_cast< double >( ( actual - expected ).norm( ) ), 1.0E-12 );
            bodies.at( member )->recomputeStateOnNextCall( );
            const LongState globalState = bodies.at( member )->getStateInBaseFrameFromEphemeris< long double, Time >( Time( time ) );
            BOOST_CHECK_SMALL( static_cast< double >( ( globalState - expected ).norm( ) ), 1.0E-12 );
        }
        // A repeated epoch with a new integration stage must use the new state.
        displaced[ 0 ] += 10.0L;
        center->beginStateUpdate( Time( time ) );
        center->setCurrentState( displaced );
        BOOST_CHECK_SMALL( static_cast< double >( ( bodies.at( "Io" )->getEphemeris( )->getCartesianLongState( time ) - displaced -
                                                    center->getMemberOffset( "Io", Time( time ) ) )
                                                          .norm( ) ),
                           positionRoundoff );
        center->clearCurrentState( );
    }
    SelectedAccelerationMap forbidden;
    forbidden[ "Juno" ][ centerName ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );
    BOOST_CHECK_THROW( createAccelerationModelsMap( bodies, forbidden, { "Juno" }, { "SSB" } ), std::runtime_error );
}

BOOST_AUTO_TEST_CASE( junoSingleArc )
{
    runJunoEstimation( false );
}
BOOST_AUTO_TEST_CASE( junoThreeFlybysHybridArc )
{
    runJunoEstimation( true );
}

BOOST_AUTO_TEST_SUITE_END( )
}  // namespace unit_tests
}  // namespace tudat
