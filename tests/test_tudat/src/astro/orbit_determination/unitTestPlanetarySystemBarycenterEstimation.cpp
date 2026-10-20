/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Tests the estimation of the state of a propagated body (here a planetary system barycenter) w.r.t. which the
 *    ephemerides of other bodies (planet and moons) are defined, together with the state of a spacecraft propagated w.r.t.
 *    such a body (here the planet). The setup is a simplified, synthetic, Juno-like case: the Jupiter system barycenter is
 *    propagated in a single arc, and the spacecraft in separate arcs around a number of perijoves.
 */

#define BOOST_TEST_MAIN

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <string>
#include <vector>

#include <boost/test/included/unit_test.hpp>

#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/astro/ephemerides/customEphemeris.h"
#include "tudat/astro/ephemerides/tabulatedEphemeris.h"
#include "tudat/astro/gravitation/gravityFieldModel.h"
#include "tudat/math/interpolators/createInterpolator.h"
#include "tudat/simulation/estimation_setup/createEstimatableParametersFactory.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"
#include "tudat/simulation/estimation_setup/singleArcVariationalEquationsSolver.h"
#include "tudat/simulation/propagation_setup/createAccelerationModels.h"
#include "tudat/simulation/propagation_setup/dynamicsSimulator.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::simulation_setup;
using namespace tudat::propagators;
using namespace tudat::observation_models;
using namespace tudat::estimatable_parameters;
using namespace tudat::numerical_integrators;
using namespace tudat::basic_astrodynamics;

namespace
{

const std::string barycenterName = "Jupiter_barycenter";
const double sunGravitationalParameter = 1.32712440018E20;
const double jupiterGravitationalParameter = 1.26686534E17;
const std::vector< std::string > moonNames = { "Io", "Europa", "Ganymede", "Callisto" };
const std::vector< double > moonGravitationalParameters = { 5.9599E12, 3.2027E12, 9.8878E12, 7.1793E12 };
const std::vector< double > moonOrbitRadii = { 4.217E8, 6.711E8, 1.0704E9, 1.8827E9 };

// Time between subsequent perijoves, and duration of each spacecraft arc
const double perijovePeriod = 53.0 * 86400.0;
const double arcDuration = 6.0 * 3600.0;

// Analytical circular orbit, used for synthetic (velocity-consistent) ephemerides
Eigen::Vector6d getCircularOrbitState( const double time,
                                       const double radius,
                                       const double gravitationalParameter,
                                       const double phase,
                                       const double inclination )
{
    const double meanMotion = std::sqrt( gravitationalParameter / ( radius * radius * radius ) );
    const double angle = meanMotion * time + phase;
    Eigen::Vector6d state;
    state << radius * std::cos( angle ), radius * std::sin( angle ) * std::cos( inclination ),
            radius * std::sin( angle ) * std::sin( inclination ), -radius * meanMotion * std::sin( angle ),
            radius * meanMotion * std::cos( angle ) * std::cos( inclination ),
            radius * meanMotion * std::cos( angle ) * std::sin( inclination );
    return state;
}

Eigen::Vector6d getNominalBarycenterState( const double time )
{
    return getCircularOrbitState( time, 7.785E11, sunGravitationalParameter, 0.4, 0.02 );
}

Eigen::Vector6d getMoonStateWrtJupiter( const int moonIndex, const double time )
{
    return getCircularOrbitState(
            time, moonOrbitRadii.at( moonIndex ), jupiterGravitationalParameter, 0.7 * moonIndex, 0.01 * ( moonIndex + 1 ) );
}

// State of Jupiter w.r.t. system barycenter, consistent with the moon states
Eigen::Vector6d getJupiterStateWrtBarycenter( const double time )
{
    Eigen::Vector6d weightedMoonStates = Eigen::Vector6d::Zero( );
    double totalGravitationalParameter = jupiterGravitationalParameter;
    for( unsigned int i = 0; i < moonNames.size( ); i++ )
    {
        weightedMoonStates += moonGravitationalParameters.at( i ) * getMoonStateWrtJupiter( i, time );
        totalGravitationalParameter += moonGravitationalParameters.at( i );
    }
    return -weightedMoonStates / totalGravitationalParameter;
}

// Juno-like eccentric, near-polar Jupiter orbit, starting about two hours before perijove
Eigen::Vector6d getSpacecraftInitialState( const int arcIndex )
{
    Eigen::Vector6d keplerianElements;
    keplerianElements << 1.5E9, 0.95, 1.55, 0.4 + 0.15 * arcIndex, 0.2 + 0.3 * arcIndex, -2.0;
    return orbital_element_conversions::convertKeplerianToCartesianElements( keplerianElements, jupiterGravitationalParameter );
}

// Create environment in which Jupiter and the moons have ephemerides w.r.t. the Jupiter system barycenter, which is a
// regular body. Io and Europa are defined w.r.t. Jupiter (and indirectly w.r.t. the barycenter), Ganymede and Callisto directly
// w.r.t. the barycenter.
SystemOfBodies createJovianSystemBodies( const bool spacecraftIsMultiArc )
{
    SystemOfBodies bodies = SystemOfBodies( "SSB", "ECLIPJ2000" );
    for( std::string bodyName :
         std::vector< std::string >( { "Sun", "Earth", barycenterName, "Jupiter", "Io", "Europa", "Ganymede", "Callisto", "Juno" } ) )
    {
        bodies.createEmptyBody( bodyName, false );
    }

    bodies.at( "Sun" )->setEphemeris( std::make_shared< ephemerides::CustomEphemeris<> >(
            []( const double ) -> Eigen::Vector6d { return Eigen::Vector6d::Zero( ); }, "SSB", "ECLIPJ2000" ) );
    bodies.at( "Sun" )->setGravityFieldModel( std::make_shared< gravitation::GravityFieldModel >( sunGravitationalParameter ) );
    bodies.at( "Earth" )->setEphemeris( std::make_shared< ephemerides::CustomEphemeris<> >(
            []( const double time ) -> Eigen::Vector6d {
                return getCircularOrbitState( time, 1.496E11, sunGravitationalParameter, 1.3, 0.0 );
            },
            "SSB",
            "ECLIPJ2000" ) );

    // Barycenter: tabulated (propagated) ephemeris w.r.t. SSB, nominal values from analytical orbit
    std::map< double, Eigen::Vector6d > nominalBarycenterStates;
    for( double time = -5.0 * 86400.0; time < 3.0 * perijovePeriod + 5.0 * 86400.0; time += 3600.0 )
    {
        nominalBarycenterStates[ time ] = getNominalBarycenterState( time );
    }
    bodies.at( barycenterName )
            ->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris<> >(
                    interpolators::createOneDimensionalInterpolator( nominalBarycenterStates,
                                                                     std::make_shared< interpolators::LagrangeInterpolatorSettings >( 8 ) ),
                    "SSB",
                    "ECLIPJ2000" ) );

    bodies.at( "Jupiter" )
            ->setEphemeris( std::make_shared< ephemerides::CustomEphemeris<> >(
                    []( const double time ) -> Eigen::Vector6d { return getJupiterStateWrtBarycenter( time ); },
                    barycenterName,
                    "ECLIPJ2000" ) );
    bodies.at( "Jupiter" )->setGravityFieldModel( std::make_shared< gravitation::GravityFieldModel >( jupiterGravitationalParameter ) );
    for( unsigned int i = 0; i < moonNames.size( ); i++ )
    {
        if( i < 2 )
        {
            bodies.at( moonNames.at( i ) )
                    ->setEphemeris( std::make_shared< ephemerides::CustomEphemeris<> >(
                            [ = ]( const double time ) -> Eigen::Vector6d { return getMoonStateWrtJupiter( i, time ); },
                            "Jupiter",
                            "ECLIPJ2000" ) );
        }
        else
        {
            bodies.at( moonNames.at( i ) )
                    ->setEphemeris( std::make_shared< ephemerides::CustomEphemeris<> >(
                            [ = ]( const double time ) -> Eigen::Vector6d {
                                return getMoonStateWrtJupiter( i, time ) + getJupiterStateWrtBarycenter( time );
                            },
                            barycenterName,
                            "ECLIPJ2000" ) );
        }
        bodies.at( moonNames.at( i ) )
                ->setGravityFieldModel( std::make_shared< gravitation::GravityFieldModel >( moonGravitationalParameters.at( i ) ) );
    }

    addEmptyTabulatedEphemeris< double, double >( bodies, "Juno", "Jupiter", spacecraftIsMultiArc );
    bodies.processBodyFrameDefinitions( );
    return bodies;
}

// Barycenter: any acceleration may be used; here solar point mass gravity and relativistic correction
AccelerationMap createBarycenterAccelerations( const SystemOfBodies& bodies )
{
    SelectedAccelerationMap accelerationSettings;
    accelerationSettings[ barycenterName ][ "Sun" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationSettings[ barycenterName ][ "Sun" ].push_back( std::make_shared< RelativisticAccelerationCorrectionSettings >( ) );
    return createAccelerationModelsMap( bodies, accelerationSettings, { barycenterName }, { "SSB" } );
}

// Spacecraft, propagated w.r.t. Jupiter: regular (third-body) accelerations, including velocity-dependent relativistic correction
AccelerationMap createSpacecraftAccelerations( const SystemOfBodies& bodies )
{
    SelectedAccelerationMap accelerationSettings;
    accelerationSettings[ "Juno" ][ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationSettings[ "Juno" ][ "Jupiter" ].push_back( std::make_shared< RelativisticAccelerationCorrectionSettings >( ) );
    for( std::string moonName : moonNames )
    {
        accelerationSettings[ "Juno" ][ moonName ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    }
    accelerationSettings[ "Juno" ][ "Sun" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    return createAccelerationModelsMap( bodies, accelerationSettings, { "Juno" }, { "Jupiter" } );
}

std::shared_ptr< TranslationalStatePropagatorSettings< double > > createPropagatorSettings(
        const std::vector< std::string >& centralBodies,
        const AccelerationMap& accelerations,
        const std::vector< std::string >& bodiesToPropagate,
        const Eigen::VectorXd& initialState,
        const double startTime,
        const double endTime,
        const double timeStep,
        const std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > >& dependentVariables =
                std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > >( ) )
{
    return std::make_shared< TranslationalStatePropagatorSettings< double > >(
            centralBodies,
            accelerations,
            bodiesToPropagate,
            initialState,
            startTime,
            rungeKuttaFixedStepSettings< double >( timeStep, rungeKuttaFehlberg78 ),
            std::make_shared< PropagationTimeTerminationSettings >( endTime ),
            cowell,
            dependentVariables );
}

}  // namespace

BOOST_AUTO_TEST_SUITE( test_planetary_system_barycenter_estimation )

// Check that, during propagation, bodies with an ephemeris w.r.t. a propagated body use its current propagated state, and that
// the state transition matrix (including the dependency of the spacecraft state on the barycenter state through Jupiter) is
// consistent with numerical differentiation.
BOOST_AUTO_TEST_CASE( testSingleArcPropagationAndStateTransitionMatrix )
{
    const double startTime = -2.0 * 3600.0;
    const double endTime = startTime + arcDuration;
    const double timeStep = 30.0;

    std::function< std::shared_ptr< TranslationalStatePropagatorSettings< double > >(
            const SystemOfBodies&, const Eigen::VectorXd&, const bool ) >
            createSettings = [ = ]( const SystemOfBodies& bodies, const Eigen::VectorXd& initialState, const bool saveDependentVariables ) {
                AccelerationMap accelerations = createBarycenterAccelerations( bodies );
                AccelerationMap spacecraftAccelerations = createSpacecraftAccelerations( bodies );
                accelerations.insert( spacecraftAccelerations.begin( ), spacecraftAccelerations.end( ) );

                std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
                if( saveDependentVariables )
                {
                    dependentVariables.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                            relative_position_dependent_variable, "Juno", "Jupiter" ) );
                    dependentVariables.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                            relative_position_dependent_variable, barycenterName, "Sun" ) );
                    dependentVariables.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                            relative_position_dependent_variable, "Jupiter", "Sun" ) );
                    dependentVariables.push_back(
                            std::make_shared< SingleDependentVariableSaveSettings >( relative_position_dependent_variable, "Io", "Sun" ) );
                }
                return createPropagatorSettings( { "SSB", "Jupiter" },
                                                 accelerations,
                                                 { barycenterName, "Juno" },
                                                 initialState,
                                                 startTime,
                                                 endTime,
                                                 timeStep,
                                                 dependentVariables );
            };

    Eigen::VectorXd nominalInitialState = Eigen::VectorXd::Zero( 12 );
    nominalInitialState << getNominalBarycenterState( startTime ), getSpacecraftInitialState( 0 );

    // Propagate with barycenter initial state far (10000 km) from its ephemeris: the environment must use the propagated state.
    {
        SystemOfBodies bodies = createJovianSystemBodies( false );
        Eigen::VectorXd perturbedInitialState = nominalInitialState;
        perturbedInitialState.segment( 0, 3 ) += Eigen::Vector3d( 1.0E7, -2.0E7, 0.5E7 );

        SingleArcDynamicsSimulator< double, double > dynamicsSimulator( bodies, createSettings( bodies, perturbedInitialState, true ) );
        std::map< double, Eigen::VectorXd > stateHistory = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
        std::map< double, Eigen::VectorXd > dependentVariableHistory = dynamicsSimulator.getDependentVariableHistory( );

        BOOST_CHECK( dependentVariableHistory.size( ) > 100 );
        for( auto variableIterator : dependentVariableHistory )
        {
            const double time = variableIterator.first;
            const Eigen::VectorXd& dependentVariables = variableIterator.second;
            const Eigen::Vector3d propagatedBarycenterPosition = stateHistory.at( time ).segment( 0, 3 );

            // Spacecraft position w.r.t. Jupiter in environment must be equal to propagated relative position
            BOOST_CHECK_SMALL( ( dependentVariables.segment( 0, 3 ) - stateHistory.at( time ).segment( 6, 3 ) ).norm( ), 1.0E-3 );

            // Barycenter in environment must have propagated state (not be reset from its ephemeris), and Jupiter and Io must be
            // shifted along with it.
            BOOST_CHECK_SMALL( ( dependentVariables.segment( 3, 3 ) - propagatedBarycenterPosition ).norm( ), 1.0E-3 );
            BOOST_CHECK_SMALL( ( dependentVariables.segment( 6, 3 ) - propagatedBarycenterPosition -
                                 getJupiterStateWrtBarycenter( time ).segment( 0, 3 ) )
                                       .norm( ),
                               1.0E-3 );
            BOOST_CHECK_SMALL( ( dependentVariables.segment( 9, 3 ) - propagatedBarycenterPosition -
                                 ( getJupiterStateWrtBarycenter( time ) + getMoonStateWrtJupiter( 0, time ) ).segment( 0, 3 ) )
                                       .norm( ),
                               1.0E-3 );
        }
    }

    // Compare state transition matrix of absolute states (spacecraft state shifted by barycenter state) to numerical derivatives
    {
        SystemOfBodies bodies = createJovianSystemBodies( false );
        std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                createSettings( bodies, nominalInitialState, false );
        std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate =
                createParametersToEstimate< double >( getInitialStateParameterSettings< double >( propagatorSettings, bodies ), bodies );

        SingleArcVariationalEquationsSolver< double, double > variationalEquationsSolver(
                bodies, propagatorSettings, parametersToEstimate );

        // Evaluate at integration epoch before the final time (no interpolation error)
        const double evaluationTime = startTime + 600.0 * timeStep;
        Eigen::MatrixXd stateTransitionMatrix =
                variationalEquationsSolver.getStateTransitionMatrixInterface( )->getCombinedStateTransitionAndSensitivityMatrix(
                        evaluationTime );
        BOOST_CHECK_EQUAL( stateTransitionMatrix.rows( ), 12 );
        BOOST_CHECK_EQUAL( stateTransitionMatrix.cols( ), 12 );

        // The interface provides the partials of the absolute (w.r.t. SSB) spacecraft state, which is its propagated state w.r.t.
        // Jupiter, plus the barycenter state, plus the (fixed) state of Jupiter w.r.t. the barycenter. Subtract the barycenter
        // rows to retrieve the partials of the propagated state, which are compared to numerical derivatives.
        Eigen::MatrixXd propagatedStateTransitionMatrix = stateTransitionMatrix;
        propagatedStateTransitionMatrix.block( 6, 0, 6, 12 ) -= stateTransitionMatrix.block( 0, 0, 6, 12 );

        std::function< Eigen::VectorXd( const Eigen::VectorXd& ) > getPropagatedStates = [ & ]( const Eigen::VectorXd& initialState ) {
            SingleArcDynamicsSimulator< double, double > dynamicsSimulator( bodies, createSettings( bodies, initialState, false ) );
            return Eigen::VectorXd( dynamicsSimulator.getEquationsOfMotionNumericalSolution( ).at( evaluationTime ) );
        };

        for( int column = 0; column < 12; column++ )
        {
            const double step = ( column % 6 < 3 ) ? ( column < 6 ? 1.0E3 : 1.0 ) : ( column < 6 ? 1.0E-2 : 1.0E-3 );
            Eigen::VectorXd perturbedInitialState = nominalInitialState;
            perturbedInitialState( column ) += step;
            Eigen::VectorXd upperState = getPropagatedStates( perturbedInitialState );
            perturbedInitialState( column ) -= 2.0 * step;
            Eigen::VectorXd lowerState = getPropagatedStates( perturbedInitialState );
            Eigen::VectorXd numericalColumn = ( upperState - lowerState ) / ( 2.0 * step );

            // Check position and velocity blocks of both bodies separately
            for( int block = 0; block < 4; block++ )
            {
                const Eigen::Vector3d analyticalBlock = propagatedStateTransitionMatrix.block( 3 * block, column, 3, 1 );
                const Eigen::Vector3d numericalBlock = numericalColumn.segment( 3 * block, 3 );
                BOOST_CHECK_MESSAGE( ( analyticalBlock - numericalBlock ).norm( ) <=
                                             1.0E-5 * numericalBlock.norm( ) + ( block % 2 == 0 ? 1.0E-5 : 1.0E-8 ),
                                     "STM column " << column << ", block " << block << ": analytical " << analyticalBlock.transpose( )
                                                   << ", numerical " << numericalBlock.transpose( ) );
            }
        }

        // Spacecraft absolute position depends (near-)identically on barycenter initial position, while its position w.r.t.
        // Jupiter is (near-)independent of it: the local Jovian gravity does not depend on the common translation
        BOOST_CHECK_SMALL( ( stateTransitionMatrix.block( 6, 0, 3, 3 ) - Eigen::Matrix3d::Identity( ) ).norm( ), 1.0E-2 );
        BOOST_CHECK_SMALL( propagatedStateTransitionMatrix.block( 6, 0, 6, 6 ).norm( ), 1.0E-6 );
    }
}

// Hybrid-arc estimation: single arc for the Jupiter system barycenter over three perijoves, separate spacecraft arcs at each
// perijove (propagated w.r.t. Jupiter). Observation partials are compared to numerical derivatives, and the estimation must
// recover the true initial states.
BOOST_AUTO_TEST_CASE( testHybridArcEstimationOverThreePerijoves )
{
    const int numberOfArcs = 3;
    SystemOfBodies bodies = createJovianSystemBodies( true );

    // Define arcs and observation times
    std::vector< double > arcStartTimes, arcEndTimes, observationTimes;
    for( int arc = 0; arc < numberOfArcs; arc++ )
    {
        arcStartTimes.push_back( arc * perijovePeriod );
        arcEndTimes.push_back( arc * perijovePeriod + arcDuration );
        for( double time = arcStartTimes.back( ) + 1800.0; time <= arcEndTimes.back( ) - 1800.0; time += 600.0 )
        {
            observationTimes.push_back( time );
        }
    }
    const int numberOfObservationsPerArc = static_cast< int >( observationTimes.size( ) ) / numberOfArcs;

    // Create propagation settings
    const double singleArcStartTime = arcStartTimes.front( ) - 3600.0;
    const double singleArcEndTime = arcEndTimes.back( ) + 3600.0;
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > singleArcSettings =
            createPropagatorSettings( { "SSB" },
                                      createBarycenterAccelerations( bodies ),
                                      { barycenterName },
                                      getNominalBarycenterState( singleArcStartTime ),
                                      singleArcStartTime,
                                      singleArcEndTime,
                                      3600.0 );

    AccelerationMap spacecraftAccelerations = createSpacecraftAccelerations( bodies );
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > arcSettings;
    for( int arc = 0; arc < numberOfArcs; arc++ )
    {
        arcSettings.push_back( createPropagatorSettings( { "Jupiter" },
                                                         spacecraftAccelerations,
                                                         { "Juno" },
                                                         getSpacecraftInitialState( arc ),
                                                         arcStartTimes.at( arc ),
                                                         arcEndTimes.at( arc ),
                                                         30.0 ) );
    }
    std::shared_ptr< HybridArcPropagatorSettings< double > > propagatorSettings = std::make_shared< HybridArcPropagatorSettings< double > >(
            singleArcSettings, std::make_shared< MultiArcPropagatorSettings< double > >( arcSettings ) );

    // Estimate barycenter state and spacecraft state in each arc
    std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate = createParametersToEstimate< double >(
            getInitialHybridArcParameterSettings< double >( propagatorSettings, bodies, arcStartTimes ), bodies );
    const Eigen::VectorXd truthParameters = parametersToEstimate->getFullParameterValues< double >( );
    const int numberOfParameters = truthParameters.rows( );
    BOOST_CHECK_EQUAL( numberOfParameters, 6 * ( 1 + numberOfArcs ) );

    // Define observation models: spacecraft position, one-way range and Doppler from Earth, and position of Io (ephemeris
    // w.r.t. Jupiter, w.r.t. barycenter) and Ganymede (ephemeris directly w.r.t. barycenter)
    LinkEnds spacecraftLinkEnds, earthSpacecraftLinkEnds, ioLinkEnds, ganymedeLinkEnds;
    spacecraftLinkEnds[ observed_body ] = LinkEndId( "Juno", "" );
    earthSpacecraftLinkEnds[ transmitter ] = LinkEndId( "Earth", "" );
    earthSpacecraftLinkEnds[ receiver ] = LinkEndId( "Juno", "" );
    ioLinkEnds[ observed_body ] = LinkEndId( "Io", "" );
    ganymedeLinkEnds[ observed_body ] = LinkEndId( "Ganymede", "" );

    struct ObservationTypeSettings {
        std::string name;
        ObservableType observable;
        LinkEnds linkEnds;
        LinkEndType referenceLinkEnd;
        double relativeTolerance;
    };
    const std::vector< ObservationTypeSettings > observationTypes = {
        { "Juno position", position_observable, spacecraftLinkEnds, observed_body, 2.0E-5 },
        { "range", one_way_range, earthSpacecraftLinkEnds, receiver, 2.0E-5 },
        { "Doppler", one_way_doppler, earthSpacecraftLinkEnds, receiver, 1.0E-4 },
        { "Io position", position_observable, ioLinkEnds, observed_body, 2.0E-5 },
        { "Ganymede position", position_observable, ganymedeLinkEnds, observed_body, 2.0E-5 }
    };

    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettings;
    for( const ObservationTypeSettings& type : observationTypes )
    {
        observationSettings.push_back( std::make_shared< ObservationModelSettings >( type.observable, type.linkEnds ) );
    }

    OrbitDeterminationManager< double, double > orbitDeterminationManager(
            bodies, parametersToEstimate, observationSettings, propagatorSettings );

    // Compute observations and partials of all types (concatenated), and start row of each type
    std::vector< int > typeStartIndices = { 0 };
    std::function< void( Eigen::VectorXd&, Eigen::MatrixXd& ) > computeObservations = [ & ]( Eigen::VectorXd& observations,
                                                                                             Eigen::MatrixXd& partials ) {
        std::vector< Eigen::VectorXd > observationsPerType( observationTypes.size( ) );
        std::vector< Eigen::MatrixXd > partialsPerType( observationTypes.size( ) );
        typeStartIndices = { 0 };
        for( unsigned int i = 0; i < observationTypes.size( ); i++ )
        {
            const ObservationTypeSettings& type = observationTypes.at( i );
            if( getObservableSize( type.observable ) == 3 )
            {
                orbitDeterminationManager.computePartialsAndObservations< 3 >( type.linkEnds,
                                                                               nullptr,
                                                                               type.observable,
                                                                               type.referenceLinkEnd,
                                                                               observationTimes,
                                                                               observationsPerType[ i ],
                                                                               partialsPerType[ i ] );
            }
            else
            {
                orbitDeterminationManager.computePartialsAndObservations< 1 >( type.linkEnds,
                                                                               nullptr,
                                                                               type.observable,
                                                                               type.referenceLinkEnd,
                                                                               observationTimes,
                                                                               observationsPerType[ i ],
                                                                               partialsPerType[ i ] );
            }
            typeStartIndices.push_back( typeStartIndices.back( ) + observationsPerType[ i ].rows( ) );
        }

        observations.resize( typeStartIndices.back( ) );
        partials.resize( typeStartIndices.back( ), numberOfParameters );
        for( unsigned int i = 0; i < observationTypes.size( ); i++ )
        {
            observations.segment( typeStartIndices[ i ], observationsPerType[ i ].rows( ) ) = observationsPerType[ i ];
            partials.block( typeStartIndices[ i ], 0, partialsPerType[ i ].rows( ), numberOfParameters ) = partialsPerType[ i ];
        }
    };

    Eigen::VectorXd nominalObservations;
    Eigen::MatrixXd analyticalPartials;
    computeObservations( nominalObservations, analyticalPartials );

    // Numerical noise floor per type, relative to largest partial (for weak dependencies, e.g. Doppler on barycenter position)
    std::vector< double > typeAbsoluteTolerances;
    for( unsigned int type = 0; type < observationTypes.size( ); type++ )
    {
        double maximumColumnNorm = 0.0;
        for( int column = 0; column < numberOfParameters; column++ )
        {
            maximumColumnNorm = std::max(
                    maximumColumnNorm,
                    analyticalPartials.block( typeStartIndices[ type ], column, typeStartIndices[ type + 1 ] - typeStartIndices[ type ], 1 )
                            .norm( ) );
        }
        typeAbsoluteTolerances.push_back( 1.0E-8 * maximumColumnNorm );
    }

    // Partials of moon positions: identity w.r.t. barycenter initial position (shifted by propagation), none w.r.t. spacecraft
    for( unsigned int type = 3; type < 5; type++ )
    {
        BOOST_CHECK( analyticalPartials.block( typeStartIndices[ type ], 0, typeStartIndices[ type + 1 ] - typeStartIndices[ type ], 6 )
                             .norm( ) > 1.0 );
        BOOST_CHECK_SMALL( analyticalPartials
                                   .block( typeStartIndices[ type ],
                                           6,
                                           typeStartIndices[ type + 1 ] - typeStartIndices[ type ],
                                           numberOfParameters - 6 )
                                   .norm( ),
                           1.0E-15 );
    }

    // Spacecraft observations in each arc depend on barycenter and spacecraft state in that arc only
    for( int arc = 0; arc < numberOfArcs; arc++ )
    {
        const int startRow = 3 * arc * numberOfObservationsPerArc;
        const int numberOfRows = 3 * numberOfObservationsPerArc;
        BOOST_CHECK( analyticalPartials.block( startRow, 0, numberOfRows, 6 ).norm( ) > 1.0 );
        BOOST_CHECK( analyticalPartials.block( startRow, 6 * ( arc + 1 ), numberOfRows, 6 ).norm( ) > 1.0 );
        for( int otherArc = 0; otherArc < numberOfArcs; otherArc++ )
        {
            if( otherArc != arc )
            {
                BOOST_CHECK_SMALL( analyticalPartials.block( startRow, 6 * ( otherArc + 1 ), numberOfRows, 6 ).norm( ), 1.0E-15 );
            }
        }
    }

    // Compare observation partials to central differences, for all initial state parameters
    for( int column = 0; column < numberOfParameters; column++ )
    {
        const double step = ( column % 6 < 3 ) ? ( column < 6 ? 1.0E5 : 1.0E3 ) : ( column < 6 ? 0.1 : 0.01 );
        Eigen::VectorXd perturbedParameters = truthParameters;
        Eigen::VectorXd upperObservations, lowerObservations;
        Eigen::MatrixXd unusedPartials;

        perturbedParameters( column ) += step;
        orbitDeterminationManager.resetParameterEstimate( perturbedParameters, false );
        computeObservations( upperObservations, unusedPartials );

        perturbedParameters( column ) -= 2.0 * step;
        orbitDeterminationManager.resetParameterEstimate( perturbedParameters, false );
        computeObservations( lowerObservations, unusedPartials );

        const Eigen::VectorXd numericalPartials = ( upperObservations - lowerObservations ) / ( 2.0 * step );
        for( unsigned int type = 0; type < observationTypes.size( ); type++ )
        {
            const int numberOfRows = typeStartIndices[ type + 1 ] - typeStartIndices[ type ];
            const Eigen::VectorXd numericalBlock = numericalPartials.segment( typeStartIndices[ type ], numberOfRows );
            const Eigen::VectorXd analyticalBlock = analyticalPartials.block( typeStartIndices[ type ], column, numberOfRows, 1 );

            // Moon positions do not depend on spacecraft states
            if( type >= 3 && column >= 6 )
            {
                continue;
            }
            const double relativeError = ( numericalBlock - analyticalBlock ).norm( ) / analyticalBlock.norm( );
            BOOST_CHECK_MESSAGE(
                    ( numericalBlock - analyticalBlock ).norm( ) <
                            observationTypes.at( type ).relativeTolerance * analyticalBlock.norm( ) + typeAbsoluteTolerances.at( type ),
                    observationTypes.at( type ).name << " partial, parameter " << column << ": relative error " << relativeError
                                                     << ", analytical norm " << analyticalBlock.norm( ) << ", numerical norm "
                                                     << numericalBlock.norm( ) );
        }
    }
    orbitDeterminationManager.resetParameterEstimate( truthParameters, false );

    // Simulate (noise-free) spacecraft observations with true parameters
    std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > simulationSettings;
    for( unsigned int i = 0; i < 3; i++ )
    {
        simulationSettings.push_back(
                std::make_shared< TabulatedObservationSimulationSettings<> >( observationTypes.at( i ).observable,
                                                                              observationTypes.at( i ).linkEnds,
                                                                              observationTimes,
                                                                              observationTypes.at( i ).referenceLinkEnd ) );
    }
    std::shared_ptr< ObservationCollection< double, double > > observations =
            simulateObservations< double, double >( simulationSettings, orbitDeterminationManager.getObservationSimulators( ), bodies );

    std::map< std::shared_ptr< ObservationCollectionParser >, double > weightPerObservable;
    weightPerObservable[ observationParser( position_observable ) ] = 1.0;
    weightPerObservable[ observationParser( one_way_range ) ] = 1.0;
    weightPerObservable[ observationParser( one_way_doppler ) ] = 1.0E16;
    observations->setConstantWeightPerObservable( weightPerObservable );

    // Perturb initial states, and estimate
    Eigen::VectorXd initialParameterEstimate = truthParameters;
    for( int i = 0; i < numberOfParameters; i++ )
    {
        initialParameterEstimate( i ) += ( i % 6 < 3 ? 100.0 : 1.0E-3 ) * ( i % 2 == 0 ? 1.0 : -1.0 );
    }
    parametersToEstimate->resetParameterValues( initialParameterEstimate );

    std::shared_ptr< EstimationInput< double, double > > estimationInput =
            std::make_shared< EstimationInput< double, double > >( observations );
    estimationInput->setConvergenceChecker( std::make_shared< EstimationConvergenceChecker >( 6 ) );
    auto estimationOutput = orbitDeterminationManager.estimateParameters( estimationInput );

    const Eigen::VectorXd estimationError = estimationOutput->parameterEstimate_ - truthParameters;
    for( int i = 0; i < numberOfParameters; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( estimationError( i ) ), ( i % 6 < 3 ) ? 1.0E-1 : 1.0E-6 );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
