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
 *    other bodies, in a simplified Juno-like setup: Jupiter and Ganymede have an ephemeris w.r.t. the propagated barycenter,
 *    Io w.r.t. Jupiter (so indirectly w.r.t. the barycenter), and the spacecraft is propagated w.r.t. Jupiter.
 */

#define BOOST_TEST_MAIN

#include <cmath>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include <boost/test/included/unit_test.hpp>

#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/astro/ephemerides/customEphemeris.h"
#include "tudat/astro/ephemerides/tabulatedEphemeris.h"
#include "tudat/astro/gravitation/gravityFieldModel.h"
#include "tudat/math/interpolators/createInterpolator.h"
#include "tudat/simulation/estimation_setup/createEstimatableParametersFactory.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"
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
const double perijovePeriod = 53.0 * 86400.0;
const double arcDuration = 6.0 * 3600.0;

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

// Reflex motion of Jupiter w.r.t. barycenter, opposite to (and scaled from) the motion of Ganymede
Eigen::Vector6d getJupiterStateWrtBarycenter( const double time )
{
    return -1.0E-4 * getCircularOrbitState( time, 1.0704E9, jupiterGravitationalParameter, 2.0 );
}

Eigen::Vector6d getIoStateWrtJupiter( const double time )
{
    return getCircularOrbitState( time, 4.217E8, jupiterGravitationalParameter, 0.0 );
}

Eigen::Vector6d getGanymedeStateWrtBarycenter( const double time )
{
    return getCircularOrbitState( time, 1.0704E9, jupiterGravitationalParameter, 2.0 );
}

// Juno-like eccentric, near-polar orbit w.r.t. Jupiter, starting about two hours before perijove
Eigen::Vector6d getSpacecraftInitialState( const int arcIndex )
{
    Eigen::Vector6d keplerianElements;
    keplerianElements << 1.5E9, 0.95, 1.55, 0.4 + 0.15 * arcIndex, 0.2 + 0.3 * arcIndex, -2.0;
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
    for( std::string bodyName : { "Sun", "Earth", "Jupiter_barycenter", "Jupiter", "Io", "Ganymede", "Juno" } )
    {
        bodies.createEmptyBody( bodyName, false );
    }
    bodies.at( "Sun" )->setEphemeris( createCustomEphemeris( []( const double ) { return Eigen::Vector6d::Zero( ).eval( ); }, "SSB" ) );
    bodies.at( "Earth" )->setEphemeris( createCustomEphemeris(
            []( const double time ) { return getCircularOrbitState( time, 1.496E11, sunGravitationalParameter, 1.3 ); }, "SSB" ) );
    bodies.at( "Jupiter" )->setEphemeris( createCustomEphemeris( getJupiterStateWrtBarycenter, barycenterName ) );
    bodies.at( "Io" )->setEphemeris( createCustomEphemeris( getIoStateWrtJupiter, "Jupiter" ) );
    bodies.at( "Ganymede" )->setEphemeris( createCustomEphemeris( getGanymedeStateWrtBarycenter, barycenterName ) );

    bodies.at( "Sun" )->setGravityFieldModel( std::make_shared< gravitation::GravityFieldModel >( sunGravitationalParameter ) );
    bodies.at( "Jupiter" )->setGravityFieldModel( std::make_shared< gravitation::GravityFieldModel >( jupiterGravitationalParameter ) );
    bodies.at( "Io" )->setGravityFieldModel( std::make_shared< gravitation::GravityFieldModel >( 5.9599E12 ) );
    bodies.at( "Ganymede" )->setGravityFieldModel( std::make_shared< gravitation::GravityFieldModel >( 9.8878E12 ) );

    // Tabulated (propagated) barycenter ephemeris w.r.t. SSB, with nominal states
    std::map< double, Eigen::Vector6d > nominalBarycenterStates;
    for( double time = -5.0 * 86400.0; time < 2.0 * perijovePeriod; time += 3600.0 )
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

// Barycenter is propagated w.r.t. the SSB (perturbed by the Sun), spacecraft w.r.t. Jupiter (perturbed by Jupiter, moons and Sun)
std::shared_ptr< TranslationalStatePropagatorSettings< double > > createPropagatorSettings(
        const SystemOfBodies& bodies,
        const std::vector< std::string >& bodiesToPropagate,
        const Eigen::VectorXd& initialState,
        const double startTime,
        const double endTime,
        const double timeStep,
        const std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > >& dependentVariables =
                std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > >( ) )
{
    SelectedAccelerationMap accelerationSettings;
    accelerationSettings[ barycenterName ][ "Sun" ] = { pointMassGravityAcceleration( ) };
    accelerationSettings[ "Juno" ][ "Jupiter" ] = { pointMassGravityAcceleration( ),
                                                    std::make_shared< RelativisticAccelerationCorrectionSettings >( ) };
    for( std::string perturbingBody : { "Io", "Ganymede", "Sun" } )
    {
        accelerationSettings[ "Juno" ][ perturbingBody ] = { pointMassGravityAcceleration( ) };
    }

    std::vector< std::string > centralBodies;
    SelectedAccelerationMap accelerationsToUse;
    for( const std::string& bodyName : bodiesToPropagate )
    {
        centralBodies.push_back( bodyName == barycenterName ? "SSB" : "Jupiter" );
        accelerationsToUse[ bodyName ] = accelerationSettings.at( bodyName );
    }
    return std::make_shared< TranslationalStatePropagatorSettings< double > >(
            centralBodies,
            createAccelerationModelsMap( bodies, accelerationsToUse, bodiesToPropagate, centralBodies ),
            bodiesToPropagate,
            initialState,
            startTime,
            rungeKuttaFixedStepSettings< double >( timeStep, rungeKuttaFehlberg78 ),
            std::make_shared< PropagationTimeTerminationSettings >( endTime ),
            cowell,
            dependentVariables );
}

// Compare partials of spacecraft position and Doppler, and of Io and Ganymede position, w.r.t. all estimated initial states
// (barycenter first) to central differences.
void checkObservationPartials( const SystemOfBodies& bodies,
                               const std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate,
                               const std::shared_ptr< PropagatorSettings< double > > propagatorSettings,
                               const std::vector< double >& observationTimes )
{
    LinkEnds earthToSpacecraft;
    earthToSpacecraft[ transmitter ] = LinkEndId( "Earth", "" );
    earthToSpacecraft[ receiver ] = LinkEndId( "Juno", "" );
    const std::vector< std::tuple< ObservableType, LinkEnds, LinkEndType > > observationTypes = {
        { position_observable, LinkEnds( { { observed_body, LinkEndId( "Juno", "" ) } } ), observed_body },
        { one_way_doppler, earthToSpacecraft, receiver },
        { position_observable, LinkEnds( { { observed_body, LinkEndId( "Io", "" ) } } ), observed_body },
        { position_observable, LinkEnds( { { observed_body, LinkEndId( "Ganymede", "" ) } } ), observed_body }
    };

    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettings;
    for( auto observationType : observationTypes )
    {
        observationSettings.push_back(
                std::make_shared< ObservationModelSettings >( std::get< 0 >( observationType ), std::get< 1 >( observationType ) ) );
    }
    OrbitDeterminationManager< double, double > orbitDeterminationManager(
            bodies, parametersToEstimate, observationSettings, propagatorSettings );

    auto computeObservationsAndPartials = [ & ]( const int typeIndex, Eigen::VectorXd& observations, Eigen::MatrixXd& partials ) {
        const auto& type = observationTypes.at( typeIndex );
        if( getObservableSize( std::get< 0 >( type ) ) == 3 )
        {
            orbitDeterminationManager.computePartialsAndObservations< 3 >( std::get< 1 >( type ),
                                                                           nullptr,
                                                                           std::get< 0 >( type ),
                                                                           std::get< 2 >( type ),
                                                                           observationTimes,
                                                                           observations,
                                                                           partials );
        }
        else
        {
            orbitDeterminationManager.computePartialsAndObservations< 1 >( std::get< 1 >( type ),
                                                                           nullptr,
                                                                           std::get< 0 >( type ),
                                                                           std::get< 2 >( type ),
                                                                           observationTimes,
                                                                           observations,
                                                                           partials );
        }
    };

    const Eigen::VectorXd nominalParameters = parametersToEstimate->getFullParameterValues< double >( );
    for( unsigned int type = 0; type < observationTypes.size( ); type++ )
    {
        Eigen::VectorXd observations, upperObservations, lowerObservations;
        Eigen::MatrixXd analyticalPartials, unusedPartials;
        orbitDeterminationManager.resetParameterEstimate( nominalParameters, false );
        computeObservationsAndPartials( type, observations, analyticalPartials );

        for( int column = 0; column < nominalParameters.rows( ); column++ )
        {
            // Position and velocity steps; larger for the barycenter (first six parameters)
            const double step = ( column % 6 < 3 ) ? ( column < 6 ? 1.0E5 : 1.0E3 ) : ( column < 6 ? 0.1 : 0.01 );
            Eigen::VectorXd perturbedParameters = nominalParameters;
            perturbedParameters( column ) += step;
            orbitDeterminationManager.resetParameterEstimate( perturbedParameters, false );
            computeObservationsAndPartials( type, upperObservations, unusedPartials );
            perturbedParameters( column ) -= 2.0 * step;
            orbitDeterminationManager.resetParameterEstimate( perturbedParameters, false );
            computeObservationsAndPartials( type, lowerObservations, unusedPartials );

            // Tolerance relative to partial itself, and to largest partial of this type (for weak dependencies)
            const Eigen::VectorXd numericalPartial = ( upperObservations - lowerObservations ) / ( 2.0 * step );
            const Eigen::VectorXd analyticalPartial = analyticalPartials.col( column );
            const double tolerance = 1.0E-4 * analyticalPartial.norm( ) + 1.0E-8 * analyticalPartials.colwise( ).norm( ).maxCoeff( );
            BOOST_CHECK_MESSAGE( ( numericalPartial - analyticalPartial ).norm( ) < tolerance,
                                 "Observable " << type << ", parameter " << column << ": analytical norm " << analyticalPartial.norm( )
                                               << ", numerical norm " << numericalPartial.norm( ) );
        }
    }
    orbitDeterminationManager.resetParameterEstimate( nominalParameters, false );
}

}  // namespace

BOOST_AUTO_TEST_SUITE( test_planetary_system_barycenter_estimation )

// Check that, during propagation, bodies with an ephemeris (directly or indirectly) w.r.t. the propagated barycenter, and the
// spacecraft propagated w.r.t. such a body, use its propagated state (10000 km from its ephemeris), and not its ephemeris.
BOOST_AUTO_TEST_CASE( testPropagatedEphemerisOriginInEnvironment )
{
    SystemOfBodies bodies = createBodies( false );
    const double startTime = -2.0 * 3600.0;

    Eigen::VectorXd initialState = Eigen::VectorXd( 12 );
    initialState << getNominalBarycenterState( startTime ), getSpacecraftInitialState( 0 );
    initialState.segment( 0, 3 ) += Eigen::Vector3d( 1.0E7, -2.0E7, 0.5E7 );

    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
    for( std::string bodyName : { "Jupiter_barycenter", "Jupiter", "Io", "Ganymede", "Juno" } )
    {
        dependentVariables.push_back( relativePositionDependentVariable( bodyName, "Sun" ) );
    }
    SingleArcDynamicsSimulator< double, double > dynamicsSimulator(
            bodies,
            createPropagatorSettings(
                    bodies, { barycenterName, "Juno" }, initialState, startTime, startTime + arcDuration, 30.0, dependentVariables ) );
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
    const double startTime = -2.0 * 3600.0;

    Eigen::VectorXd initialState = Eigen::VectorXd( 12 );
    initialState << getNominalBarycenterState( startTime ), getSpacecraftInitialState( 0 );
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            createPropagatorSettings( bodies, { barycenterName, "Juno" }, initialState, startTime, startTime + arcDuration, 30.0 );
    std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double >( getInitialStateParameterSettings< double >( propagatorSettings, bodies ), bodies );

    std::vector< double > observationTimes;
    for( double time = startTime + 1800.0; time < startTime + arcDuration - 1800.0; time += 1200.0 )
    {
        observationTimes.push_back( time );
    }
    checkObservationPartials( bodies, parametersToEstimate, propagatorSettings, observationTimes );
}

// Hybrid-arc estimation: barycenter in a single arc, spacecraft in two arcs (one per perijove)
BOOST_AUTO_TEST_CASE( testHybridArcObservationPartials )
{
    SystemOfBodies bodies = createBodies( true );

    std::vector< double > arcStartTimes, observationTimes;
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > arcSettings;
    for( int arc = 0; arc < 2; arc++ )
    {
        arcStartTimes.push_back( arc * perijovePeriod );
        arcSettings.push_back( createPropagatorSettings(
                bodies, { "Juno" }, getSpacecraftInitialState( arc ), arcStartTimes.back( ), arcStartTimes.back( ) + arcDuration, 30.0 ) );
        for( double time = arcStartTimes.back( ) + 1800.0; time < arcStartTimes.back( ) + arcDuration - 1800.0; time += 1200.0 )
        {
            observationTimes.push_back( time );
        }
    }
    const double singleArcStartTime = arcStartTimes.front( ) - 3600.0;
    std::shared_ptr< HybridArcPropagatorSettings< double > > propagatorSettings = std::make_shared< HybridArcPropagatorSettings< double > >(
            createPropagatorSettings( bodies,
                                      { barycenterName },
                                      getNominalBarycenterState( singleArcStartTime ),
                                      singleArcStartTime,
                                      arcStartTimes.back( ) + arcDuration + 3600.0,
                                      3600.0 ),
            std::make_shared< MultiArcPropagatorSettings< double > >( arcSettings ) );
    std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate = createParametersToEstimate< double >(
            getInitialHybridArcParameterSettings< double >( propagatorSettings, bodies, arcStartTimes ), bodies );
    BOOST_CHECK_EQUAL( parametersToEstimate->getParameterSetSize( ), 18 );

    checkObservationPartials( bodies, parametersToEstimate, propagatorSettings, observationTimes );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
