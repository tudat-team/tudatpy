/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Tests the propagation and estimation of a body (here a planetary system barycenter) that is the ephemeris origin of
 *    other bodies, in a simplified Jupiter flyby setup: Jupiter and Ganymede have an ephemeris w.r.t. the propagated barycenter,
 *    Io w.r.t. Jupiter (so indirectly w.r.t. the barycenter), and the spacecraft is propagated w.r.t. Jupiter.
 */

#define BOOST_TEST_MAIN

#include <cmath>
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
const double saturnGravitationalParameter = 3.7931187E16;
const double jupiterRadius = 71492.0E3;

// Four flybys, one month apart, with a data arc from two hours before to two hours after each flyby. All arcs are propagated
// forward and backward from periapsis
const int numberOfArcs = 4;
const double flybySeparation = 30.0 * 86400.0;
const double arcDuration = 4.0 * 3600.0;

// Analytical circular orbit, used for (velocity-consistent) synthetic ephemerides
Eigen::Vector6d getCircularOrbitState( const double time, const double radius, const double gravitationalParameter, const double phase )
{
    const double meanMotion = std::sqrt( gravitationalParameter / ( radius * radius * radius ) );
    const double angle = meanMotion * time + phase;
    Eigen::Vector6d state;
    state << radius * std::cos( angle ), radius * std::sin( angle ), 0.01 * radius * std::sin( angle ),
            -radius * meanMotion * std::sin( angle ), radius * meanMotion * std::cos( angle ),
            0.01 * radius * meanMotion * std::cos( angle );
    return state;
}

Eigen::Vector6d getNominalBarycenterState( const double time )
{
    return getCircularOrbitState( time, 7.785E11, sunGravitationalParameter, 0.4 );
}

// Offset of Jupiter w.r.t. barycenter (about 10000 km, varying with a period of about a week): strongly exaggerated reflex
// motion, opposite to the motion of Ganymede
Eigen::Vector6d getJupiterStateWrtBarycenter( const double time )
{
    return -1.0E-2 * getCircularOrbitState( time, 1.0704E9, jupiterGravitationalParameter, 2.0 );
}

Eigen::Vector6d getIoStateWrtJupiter( const double time )
{
    return getCircularOrbitState( time, 4.217E8, jupiterGravitationalParameter, 0.0 );
}

Eigen::Vector6d getGanymedeStateWrtBarycenter( const double time )
{
    return getCircularOrbitState( time, 1.0704E9, jupiterGravitationalParameter, 2.0 );
}

double getPeriapsisTime( const int arcIndex )
{
    return arcIndex * flybySeparation;
}

// Observation times in a single arc, starting and ending one cadence step from the arc boundaries
std::vector< double > getArcObservationTimes( const int arcIndex, const double cadence )
{
    std::vector< double > observationTimes;
    const double arcStartTime = getPeriapsisTime( arcIndex ) - arcDuration / 2.0;
    for( double time = arcStartTime + cadence; time < arcStartTime + arcDuration; time += cadence )
    {
        observationTimes.push_back( time );
    }
    return observationTimes;
}

// Periapsis state of flyby orbit w.r.t. Jupiter (1000 km periapsis altitude, eccentricity 0.9)
Eigen::Vector6d getSpacecraftPeriapsisState( const int arcIndex )
{
    const double eccentricity = 0.9;
    Eigen::Vector6d keplerianElements;
    keplerianElements << ( jupiterRadius + 1000.0E3 ) / ( 1.0 - eccentricity ), eccentricity, 1.55, 0.4 + 0.15 * arcIndex,
            0.2 + 0.3 * arcIndex, 0.0;
    return orbital_element_conversions::convertKeplerianToCartesianElements( keplerianElements, jupiterGravitationalParameter );
}

std::shared_ptr< ephemerides::Ephemeris > createCustomEphemeris( const std::function< Eigen::Vector6d( const double ) > stateFunction,
                                                                 const std::string& origin )
{
    return std::make_shared< ephemerides::CustomEphemeris<> >( stateFunction, origin, "ECLIPJ2000" );
}

SystemOfBodies createBodies( const bool spacecraftIsMultiArc )
{
    SystemOfBodies bodies = SystemOfBodies( "SSB", "ECLIPJ2000" );
    for( std::string bodyName : { "Sun", "Earth", "Saturn", "Jupiter_barycenter", "Jupiter", "Io", "Ganymede", "Juno" } )
    {
        bodies.createEmptyBody( bodyName, false );
    }
    bodies.at( "Sun" )->setEphemeris( createCustomEphemeris( []( const double ) { return Eigen::Vector6d::Zero( ).eval( ); }, "SSB" ) );
    bodies.at( "Earth" )->setEphemeris( createCustomEphemeris(
            []( const double time ) { return getCircularOrbitState( time, 1.496E11, sunGravitationalParameter, 1.3 ); }, "SSB" ) );
    bodies.at( "Saturn" )->setEphemeris( createCustomEphemeris(
            []( const double time ) { return getCircularOrbitState( time, 1.4335E12, sunGravitationalParameter, 1.0 ); }, "SSB" ) );
    bodies.at( "Jupiter" )->setEphemeris( createCustomEphemeris( getJupiterStateWrtBarycenter, barycenterName ) );
    bodies.at( "Io" )->setEphemeris( createCustomEphemeris( getIoStateWrtJupiter, "Jupiter" ) );
    bodies.at( "Ganymede" )->setEphemeris( createCustomEphemeris( getGanymedeStateWrtBarycenter, barycenterName ) );

    bodies.at( "Sun" )->setGravityFieldModel( std::make_shared< gravitation::GravityFieldModel >( sunGravitationalParameter ) );
    bodies.at( "Saturn" )->setGravityFieldModel( std::make_shared< gravitation::GravityFieldModel >( saturnGravitationalParameter ) );
    bodies.at( "Jupiter" )->setGravityFieldModel( std::make_shared< gravitation::GravityFieldModel >( jupiterGravitationalParameter ) );
    bodies.at( "Io" )->setGravityFieldModel( std::make_shared< gravitation::GravityFieldModel >( 5.9599E12 ) );
    bodies.at( "Ganymede" )->setGravityFieldModel( std::make_shared< gravitation::GravityFieldModel >( 9.8878E12 ) );

    // Tabulated (propagated) barycenter ephemeris w.r.t. SSB, with nominal states
    std::map< double, Eigen::Vector6d > nominalBarycenterStates;
    for( double time = -5.0 * 86400.0; time < numberOfArcs * flybySeparation + 5.0 * 86400.0; time += 3600.0 )
    {
        nominalBarycenterStates[ time ] = getNominalBarycenterState( time );
    }
    bodies.at( barycenterName )
            ->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris<> >(
                    interpolators::createOneDimensionalInterpolator( nominalBarycenterStates,
                                                                     std::make_shared< interpolators::LagrangeInterpolatorSettings >( 8 ) ),
                    "SSB",
                    "ECLIPJ2000" ) );
    addEmptyTabulatedEphemeris< double, double >( bodies, "Juno", "Jupiter", spacecraftIsMultiArc );
    bodies.processBodyFrameDefinitions( );
    return bodies;
}

// Propagator settings, propagating forward and backward from the initial time, to the end and start time
std::shared_ptr< TranslationalStatePropagatorSettings< double > > createSingleArcPropagatorSettings(
        const SystemOfBodies& bodies,
        const std::vector< std::string >& bodiesToPropagate,
        const Eigen::VectorXd& initialState,
        const double initialTime,
        const double startTime,
        const double endTime,
        const double timeStep,
        const std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > >& dependentVariables =
                std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > >( ) )
{
    std::vector< std::string > centralBodies;
    SelectedAccelerationMap accelerationSettings;
    for( const std::string& bodyName : bodiesToPropagate )
    {
        if( bodyName == "Juno" )
        {
            // Spacecraft propagated w.r.t. Jupiter, perturbed by Jupiter, moons and Sun
            centralBodies.push_back( "Jupiter" );
            for( std::string perturbingBody : { "Jupiter", "Io", "Ganymede", "Sun" } )
            {
                accelerationSettings[ bodyName ][ perturbingBody ] = { pointMassGravityAcceleration( ) };
            }
        }
        else if( bodyName == barycenterName )
        {
            // Barycenter propagated w.r.t. Sun, perturbed by Sun and Saturn
            centralBodies.push_back( "Sun" );
            accelerationSettings[ bodyName ][ "Sun" ] = { pointMassGravityAcceleration( ) };
            accelerationSettings[ bodyName ][ "Saturn" ] = { pointMassGravityAcceleration( ) };
        }
    }
    return std::make_shared< TranslationalStatePropagatorSettings< double > >(
            centralBodies,
            createAccelerationModelsMap( bodies, accelerationSettings, bodiesToPropagate, centralBodies ),
            bodiesToPropagate,
            initialState,
            initialTime,
            rungeKuttaFixedStepSettings< double >( timeStep, rungeKuttaFehlberg78 ),
            std::make_shared< NonSequentialPropagationTerminationSettings >(
                    std::make_shared< PropagationTimeTerminationSettings >( endTime ),
                    std::make_shared< PropagationTimeTerminationSettings >( startTime ) ),
            cowell,
            dependentVariables );
}

// Barycenter propagated in a single arc from the first periapsis (covering all spacecraft arcs), spacecraft in one arc per flyby,
// from periapsis
std::shared_ptr< HybridArcPropagatorSettings< double > > createHybridArcPropagatorSettings( const SystemOfBodies& bodies )
{
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > arcSettings;
    for( int arc = 0; arc < numberOfArcs; arc++ )
    {
        const double periapsisTime = getPeriapsisTime( arc );
        arcSettings.push_back( createSingleArcPropagatorSettings( bodies,
                                                                  { "Juno" },
                                                                  getSpacecraftPeriapsisState( arc ),
                                                                  periapsisTime,
                                                                  periapsisTime - arcDuration / 2.0,
                                                                  periapsisTime + arcDuration / 2.0,
                                                                  10.0 ) );
    }

    const double singleArcInitialTime = getPeriapsisTime( 0 );
    return std::make_shared< HybridArcPropagatorSettings< double > >(
            createSingleArcPropagatorSettings( bodies,
                                               { barycenterName },
                                               getNominalBarycenterState( singleArcInitialTime ),
                                               singleArcInitialTime,
                                               singleArcInitialTime - arcDuration / 2.0 - 3600.0,
                                               getPeriapsisTime( numberOfArcs - 1 ) + arcDuration / 2.0 + 3600.0,
                                               3600.0 ),
            std::make_shared< MultiArcPropagatorSettings< double > >( arcSettings ) );
}

std::vector< double > getObservationTimesOfAllArcs( const double cadence )
{
    std::vector< double > observationTimes;
    for( int arc = 0; arc < numberOfArcs; arc++ )
    {
        std::vector< double > arcObservationTimes = getArcObservationTimes( arc, cadence );
        observationTimes.insert( observationTimes.end( ), arcObservationTimes.begin( ), arcObservationTimes.end( ) );
    }
    return observationTimes;
}

// Relative position observations of barycenter w.r.t. Sun and spacecraft w.r.t. Jupiter
void addRelativePositionObservations( const std::vector< double >& observationTimes,
                                      std::vector< std::shared_ptr< ObservationModelSettings > >& observationSettings,
                                      std::vector< std::shared_ptr< ObservationSimulationSettings< double > > >& simulationSettings )
{
    for( std::pair< std::string, std::string > observedBodyAndObserver :
         { std::make_pair( barycenterName, std::string( "Sun" ) ), std::make_pair( std::string( "Juno" ), std::string( "Jupiter" ) ) } )
    {
        LinkEnds linkEnds;
        linkEnds[ observed_body ] = LinkEndId( observedBodyAndObserver.first, "" );
        linkEnds[ observer ] = LinkEndId( observedBodyAndObserver.second, "" );
        observationSettings.push_back( relativePositionObservableSettings( linkEnds ) );
        simulationSettings.push_back( std::make_shared< TabulatedObservationSimulationSettings< double > >(
                relative_position_observable, linkEnds, observationTimes, observed_body ) );
    }
}

// Compare partials of spacecraft, Io and Ganymede position, and of relative positions of barycenter w.r.t. Sun and spacecraft
// w.r.t. Jupiter, w.r.t. all estimated initial states (barycenter first) to central differences.
void checkObservationPartials( const SystemOfBodies& bodies,
                               const std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate,
                               const std::shared_ptr< PropagatorSettings< double > > propagatorSettings,
                               const std::vector< double >& observationTimes )
{
    // Define observation models and observations to simulate
    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettings;
    std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > simulationSettings;
    for( std::string observedBody : { "Juno", "Io", "Ganymede" } )
    {
        LinkEnds positionLinkEnds;
        positionLinkEnds[ observed_body ] = LinkEndId( observedBody, "" );
        observationSettings.push_back( std::make_shared< ObservationModelSettings >( position_observable, positionLinkEnds ) );
        simulationSettings.push_back( std::make_shared< TabulatedObservationSimulationSettings< double > >(
                position_observable, positionLinkEnds, observationTimes, observed_body ) );
    }
    addRelativePositionObservations( observationTimes, observationSettings, simulationSettings );

    OrbitDeterminationManager< double, double > orbitDeterminationManager(
            bodies, parametersToEstimate, observationSettings, propagatorSettings );

    // Compute analytical partials
    std::shared_ptr< ObservationCollection< double, double > > observations =
            simulateObservations< double, double >( simulationSettings, orbitDeterminationManager.getObservationSimulators( ), bodies );
    std::shared_ptr< CovarianceAnalysisInput< double, double > > covarianceInput =
            std::make_shared< CovarianceAnalysisInput< double, double > >( observations );
    covarianceInput->defineCovarianceSettings( true, true, true, false );
    Eigen::MatrixXd analyticalPartials = orbitDeterminationManager.computeCovariance( covarianceInput )->getUnnormalizedDesignMatrix( );

    // Compute numerical partials
    const Eigen::VectorXd nominalParameters = parametersToEstimate->getFullParameterValues< double >( );
    Eigen::MatrixXd numericalPartials = Eigen::MatrixXd::Zero( analyticalPartials.rows( ), analyticalPartials.cols( ) );
    for( int parameterIndex = 0; parameterIndex < nominalParameters.rows( ); parameterIndex++ )
    {
        // Position and velocity steps, larger for the barycenter (first six parameters)
        const bool isPositionParameter = ( parameterIndex % 6 < 3 );
        double step = isPositionParameter ? 1.0E3 : 1.0;
        if( parameterIndex < 6 )
        {
            step = isPositionParameter ? 1.0E5 : 10.0;
        }

        Eigen::VectorXd perturbedParameters = nominalParameters;
        perturbedParameters( parameterIndex ) += step;
        orbitDeterminationManager.resetParameterEstimate( perturbedParameters, false );
        Eigen::VectorXd upperturbedObservations =
                simulateObservations< double, double >( simulationSettings, orbitDeterminationManager.getObservationSimulators( ), bodies )
                        ->getConcatenatedObservations( );

        perturbedParameters( parameterIndex ) -= 2.0 * step;
        orbitDeterminationManager.resetParameterEstimate( perturbedParameters, false );
        Eigen::VectorXd downperturbedObservations =
                simulateObservations< double, double >( simulationSettings, orbitDeterminationManager.getObservationSimulators( ), bodies )
                        ->getConcatenatedObservations( );

        numericalPartials.col( parameterIndex ) = ( upperturbedObservations - downperturbedObservations ) / ( 2.0 * step );
    }
    orbitDeterminationManager.resetParameterEstimate( nominalParameters, false );

    // Compare partials per observable type, with tolerance relative to the partial itself
    for( auto observableIterator : observations->getObservationTypeStartAndSize( ) )
    {
        const int startIndex = observableIterator.second.first;
        const int numberOfObservations = observableIterator.second.second;
        const Eigen::MatrixXd currentAnalyticalPartials = analyticalPartials.middleRows( startIndex, numberOfObservations );
        const Eigen::MatrixXd currentNumericalPartials = numericalPartials.middleRows( startIndex, numberOfObservations );
        for( int parameterIndex = 0; parameterIndex < nominalParameters.rows( ); parameterIndex++ )
        {
            const double analyticalNorm = currentAnalyticalPartials.col( parameterIndex ).norm( );
            const double partialError =
                    ( currentNumericalPartials.col( parameterIndex ) - currentAnalyticalPartials.col( parameterIndex ) ).norm( );
            BOOST_CHECK_MESSAGE( partialError < 1.0E-5 * analyticalNorm,
                                 "Observable " << getObservableName( observableIterator.first ) << ", parameter " << parameterIndex
                                               << ": analytical norm " << analyticalNorm << ", error " << partialError );
        }
    }
}

}  // namespace

BOOST_AUTO_TEST_SUITE( test_planetary_system_barycenter_estimation )

// Check that, during propagation, bodies with an ephemeris (directly or indirectly) w.r.t. the propagated barycenter, and the
// spacecraft propagated w.r.t. such a body, use its propagated state (10000 km from its ephemeris), and not its ephemeris.
BOOST_AUTO_TEST_CASE( testPropagatedEphemerisOriginInEnvironment )
{
    SystemOfBodies bodies = createBodies( false );
    const double initialTime = getPeriapsisTime( 0 );

    Eigen::VectorXd initialState = Eigen::VectorXd( 12 );
    initialState << getNominalBarycenterState( initialTime ), getSpacecraftPeriapsisState( 0 );
    initialState.segment( 0, 3 ) += Eigen::Vector3d( 1.0E7, -2.0E7, 0.5E7 );

    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
    for( std::string bodyName : { "Jupiter_barycenter", "Jupiter", "Io", "Ganymede", "Juno" } )
    {
        dependentVariables.push_back( relativePositionDependentVariable( bodyName, "Sun" ) );
    }
    SingleArcDynamicsSimulator< double, double > dynamicsSimulator(
            bodies,
            createSingleArcPropagatorSettings( bodies,
                                               { barycenterName, "Juno" },
                                               initialState,
                                               initialTime,
                                               initialTime - arcDuration / 2.0,
                                               initialTime + arcDuration / 2.0,
                                               30.0,
                                               dependentVariables ) );
    std::map< double, Eigen::VectorXd > stateHistory = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > dependentVariableHistory = dynamicsSimulator.getDependentVariableHistory( );

    BOOST_CHECK( dependentVariableHistory.size( ) > 100 );
    for( auto variableIterator : dependentVariableHistory )
    {
        const double time = variableIterator.first;
        const Eigen::Vector3d barycenterPosition = stateHistory.at( time ).segment( 0, 3 );
        const Eigen::Vector3d jupiterPosition = barycenterPosition + getJupiterStateWrtBarycenter( time ).segment( 0, 3 );

        Eigen::VectorXd expectedPositions = Eigen::VectorXd( 15 );
        expectedPositions << barycenterPosition, jupiterPosition, jupiterPosition + getIoStateWrtJupiter( time ).segment( 0, 3 ),
                barycenterPosition + getGanymedeStateWrtBarycenter( time ).segment( 0, 3 ),
                jupiterPosition + stateHistory.at( time ).segment( 6, 3 );
        for( int i = 0; i < 5; i++ )
        {
            BOOST_CHECK_SMALL( ( variableIterator.second.segment( 3 * i, 3 ) - expectedPositions.segment( 3 * i, 3 ) ).norm( ), 1.0E-3 );
        }
    }
}

// Single-arc estimation of barycenter and spacecraft state
BOOST_AUTO_TEST_CASE( testSingleArcObservationPartials )
{
    SystemOfBodies bodies = createBodies( false );
    const double initialTime = getPeriapsisTime( 0 );

    Eigen::VectorXd initialState = Eigen::VectorXd( 12 );
    initialState << getNominalBarycenterState( initialTime ), getSpacecraftPeriapsisState( 0 );
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            createSingleArcPropagatorSettings( bodies,
                                               { barycenterName, "Juno" },
                                               initialState,
                                               initialTime,
                                               initialTime - arcDuration / 2.0,
                                               initialTime + arcDuration / 2.0,
                                               10.0 );
    std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double >( getInitialStateParameterSettings< double >( propagatorSettings, bodies ), bodies );

    checkObservationPartials( bodies, parametersToEstimate, propagatorSettings, getArcObservationTimes( 0, 1200.0 ) );
}

// Hybrid-arc estimation: barycenter in a single arc, spacecraft in one arc per flyby
BOOST_AUTO_TEST_CASE( testHybridArcObservationPartials )
{
    SystemOfBodies bodies = createBodies( true );

    std::shared_ptr< HybridArcPropagatorSettings< double > > propagatorSettings = createHybridArcPropagatorSettings( bodies );
    std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate = createParametersToEstimate< double >(
            getInitialHybridArcParameterSettings< double >( propagatorSettings, bodies ), bodies );
    BOOST_CHECK_EQUAL( parametersToEstimate->getParameterSetSize( ), 6 * ( numberOfArcs + 1 ) );

    checkObservationPartials( bodies, parametersToEstimate, propagatorSettings, getObservationTimesOfAllArcs( 1200.0 ) );
}

// Full estimation of single-arc barycenter state (w.r.t. Sun), multi-arc spacecraft state (w.r.t. Jupiter) and Jupiter
// gravitational parameter from noise-free relative position observations of the barycenter w.r.t. the Sun and the spacecraft
// w.r.t. Jupiter, every 60 s during each flyby
BOOST_AUTO_TEST_CASE( testHybridArcEstimation )
{
    SystemOfBodies bodies = createBodies( true );
    std::shared_ptr< HybridArcPropagatorSettings< double > > propagatorSettings = createHybridArcPropagatorSettings( bodies );

    // Create parameters to estimate
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterSettings =
            getInitialHybridArcParameterSettings< double >( propagatorSettings, bodies );
    parameterSettings.push_back( gravitationalParameter( "Jupiter" ) );
    std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double >( parameterSettings, bodies );

    // Create observation models and simulate observations
    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettings;
    std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > simulationSettings;
    addRelativePositionObservations( getObservationTimesOfAllArcs( 60.0 ), observationSettings, simulationSettings );
    OrbitDeterminationManager< double, double > orbitDeterminationManager(
            bodies, parametersToEstimate, observationSettings, propagatorSettings );
    std::shared_ptr< ObservationCollection< double, double > > observations =
            simulateObservations< double, double >( simulationSettings, orbitDeterminationManager.getObservationSimulators( ), bodies );

    // Perturb barycenter state (first six parameters), spacecraft states, and Jupiter gravitational parameter (last parameter)
    const Eigen::VectorXd truthParameters = parametersToEstimate->getFullParameterValues< double >( );
    Eigen::VectorXd perturbedParameters = truthParameters;
    perturbedParameters.segment( 0, 3 ) += Eigen::Vector3d::Constant( 1.0E3 );
    perturbedParameters.segment( 3, 3 ) += Eigen::Vector3d::Constant( 1.0E-3 );
    for( int arc = 0; arc < numberOfArcs; arc++ )
    {
        perturbedParameters.segment( 6 * ( arc + 1 ), 3 ) += Eigen::Vector3d::Constant( 100.0 );
        perturbedParameters.segment( 6 * ( arc + 1 ) + 3, 3 ) += Eigen::Vector3d::Constant( 1.0E-2 );
    }
    perturbedParameters( perturbedParameters.rows( ) - 1 ) *= 1.0 + 1.0E-6;
    parametersToEstimate->resetParameterValues( perturbedParameters );

    // Estimate parameters, and check that the initial errors are removed
    std::shared_ptr< EstimationInput< double, double > > estimationInput =
            std::make_shared< EstimationInput< double, double > >( observations );
    estimationInput->setConvergenceChecker( std::make_shared< EstimationConvergenceChecker >( 5 ) );
    std::shared_ptr< EstimationOutput< double, double > > estimationOutput =
            orbitDeterminationManager.estimateParameters( estimationInput );

    const Eigen::VectorXd initialError = perturbedParameters - truthParameters;
    const Eigen::VectorXd estimationError = estimationOutput->parameterEstimate_ - truthParameters;
    for( int i = 0; i < estimationError.rows( ); i++ )
    {
        BOOST_CHECK_MESSAGE( std::fabs( estimationError( i ) ) < 1.0E-4 * std::fabs( initialError( i ) ),
                             "Parameter " << i << ": initial error " << initialError( i ) << ", final error " << estimationError( i ) );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
