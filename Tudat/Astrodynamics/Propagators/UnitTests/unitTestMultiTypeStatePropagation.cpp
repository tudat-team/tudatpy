/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <boost/bind.hpp>
#include <boost/make_shared.hpp>
#include <memory>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "Tudat/Basics/testMacros.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Mathematics/NumericalIntegrators/rungeKutta4Integrator.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"
#include "Tudat/Basics/basicTypedefs.h"
#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/SimulationSetup/PropagationSetup/dynamicsSimulator.h"
#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createBodies.h"
#include "Tudat/SimulationSetup/EstimationSetup/createNumericalSimulator.h"



namespace tudat
{
namespace unit_tests
{

using namespace mathematical_constants;

BOOST_AUTO_TEST_SUITE( test_hybrid_state_derivative_model )

std::map< double, Eigen::VectorXd > propagateKeplerOrbitAndMassState(
        const int simulationCase )
{
    using namespace simulation_setup;
    using namespace propagators;
    using namespace numerical_integrators;
    using namespace orbital_element_conversions;
    using namespace basic_mathematics;
    using namespace gravitation;
    using namespace numerical_integrators;
    using namespace unit_conversions;

    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation end epoch.
    const double simulationEndEpoch = tudat::physical_constants::JULIAN_DAY;

    // Set numerical integration fixed step size.
    const double fixedStepSize = 60.0;

    // Define body settings for simulation.
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings;
    bodySettings[ "Earth" ] = std::make_shared< BodySettings >( );
    bodySettings[ "Earth" ]->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >(
                Eigen::Vector6d::Zero( ), "SSB", "J2000" );
    bodySettings[ "Earth" ]->gravityFieldSettings = std::make_shared< GravityFieldSettings >( central_spice );

    // Create Earth object
    NamedBodyMap bodyMap = createBodies( bodySettings );

    // Create spacecraft object.
    bodyMap[ "Asterix" ] = std::make_shared< simulation_setup::Body >( );
    bodyMap[ "Asterix" ]->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                                            std::shared_ptr< interpolators::OneDimensionalInterpolator
                                                < double, Eigen::Vector6d > >( ), "Earth", "J2000" ) );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Define propagation settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfAsterix;
    accelerationsOfAsterix[ "Earth" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::central_gravity ) );
    accelerationMap[ "Asterix" ] = accelerationsOfAsterix;
    bodiesToPropagate.push_back( "Asterix" );
    centralBodies.push_back( "Earth" );

    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );


    // Set Keplerian elements for Asterix.
    Eigen::Vector6d asterixInitialStateInKeplerianElements;
    asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7500.0E3;
    asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.1;
    asterixInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 85.3 );
    asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
            = convertDegreesToRadians( 235.7 );
    asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
            = convertDegreesToRadians( 23.4 );
    asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 139.87 );

    // Convert Asterix state from Keplerian elements to Cartesian elements.
    double earthGravitationalParameter = bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
    Eigen::VectorXd systemInitialState = convertKeplerianToCartesianElements(
                asterixInitialStateInKeplerianElements,
                earthGravitationalParameter );


    std::shared_ptr< SingleArcPropagatorSettings< double > > translationalPropagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState,
              std::make_shared< PropagationTimeTerminationSettings >( simulationEndEpoch ) );

    // Create mass rate model and mass propagation settings
    std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;
    massRateModels[ "Vehicle" ] = std::make_shared< basic_astrodynamics::CustomMassRateModel >(
                [ ]( const double ){ return -0.01; } );
    Eigen::VectorXd initialMass = Eigen::VectorXd( 1 );
    initialMass( 0 ) = 500.0;
    std::shared_ptr< SingleArcPropagatorSettings< double > > massPropagatorSettings =
            std::make_shared< MassPropagatorSettings< double > >(
                std::vector< std::string >{ "Asterix" }, massRateModels, initialMass,
                std::make_shared< PropagationTimeTerminationSettings >( simulationEndEpoch ) );

    // Create total propagator settings, depending on current case.
    std::shared_ptr< SingleArcPropagatorSettings< double > > propagatorSettings;
    if( ( simulationCase  % 3 ) == 0 )
    {
        propagatorSettings = translationalPropagatorSettings;
    }
    else if( ( simulationCase  % 3 ) == 1 )
    {
        propagatorSettings = massPropagatorSettings;
    }
    else if( ( simulationCase  % 3 ) == 2 )
    {
        std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > >  propagatorSettingsList;
        propagatorSettingsList.push_back( translationalPropagatorSettings );
        propagatorSettingsList.push_back( massPropagatorSettings );

        propagatorSettings = std::make_shared< MultiTypePropagatorSettings< double > >(
                    propagatorSettingsList,
                    std::make_shared< PropagationTimeTerminationSettings >( simulationEndEpoch ) );
    }



    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4, 0.0, fixedStepSize );

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodyMap, integratorSettings, propagatorSettings, true, false, true );

    // Return propagated dynamics (if simulationCase < 3) or interpolated dynamics (else)
    if( simulationCase < 3 )
    {
        return dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    }
    else
    {
        std::map< double, Eigen::VectorXd > returnMap;
        double currentTime = 4.0 * fixedStepSize;

        double interpolationTimeStep = 32.1;

        // Interpolate dynamics for range of times.
        while( currentTime < simulationEndEpoch - 4.0 * fixedStepSize )
        {

            // Interpolate propagated dynamics type
            if( simulationCase == 3 )
            {
                returnMap[ currentTime ] = bodyMap[ "Asterix" ]->getEphemeris( )->getCartesianState( currentTime );
            }
            else if( simulationCase == 4 )
            {
                returnMap[ currentTime ] = Eigen::VectorXd::Zero( 1 );
                returnMap[ currentTime ]( 0 ) = bodyMap[ "Asterix" ]->getBodyMassFunction( )( currentTime );
            }
            else if( simulationCase == 5 )
            {
                returnMap[ currentTime ] = Eigen::VectorXd::Zero( 7 );
                returnMap[ currentTime ].segment( 0, 6 ) = bodyMap[ "Asterix" ]->getEphemeris( )->getCartesianState( currentTime );
                returnMap[ currentTime ]( 6 ) = bodyMap[ "Asterix" ]->getBodyMassFunction( )( currentTime );

            }

            // Increment time (non-resonant with integration time step).
            currentTime += interpolationTimeStep;
        }

        return returnMap;
    }

}
//! Test if conversion from Keplerian elements to Cartesian elements is working correctly.
BOOST_AUTO_TEST_CASE( testHybridStateDerivativeModel )
{
    // Compare separate and multitype (independent) propagation directly from propagation (useCase = 0) and
    // interpolated from reset dynamics (useCase = 1)
    for( int useCase = 0; useCase < 2; useCase++ )
    {
        int simulationCaseToAdd = 3;
        if( useCase == 0 )
        {
             simulationCaseToAdd = 0;
        }

        // Propagate dynamics for translational, mass and combined state.
        std::map< double, Eigen::VectorXd >  translationalState = propagateKeplerOrbitAndMassState( 0 + simulationCaseToAdd );
        std::map< double, Eigen::VectorXd >  massState = propagateKeplerOrbitAndMassState( 1 + simulationCaseToAdd );
        std::map< double, Eigen::VectorXd >  combinedState = propagateKeplerOrbitAndMassState( 2  + simulationCaseToAdd );

        std::map< double, Eigen::VectorXd >::const_iterator stateIterator = translationalState.begin( );
        std::map< double, Eigen::VectorXd >::const_iterator massIterator = massState.begin( );
        std::map< double, Eigen::VectorXd >::const_iterator combinedIterator = combinedState.begin( );

        // Compare separate and multitype dynamics of each type.
        for( unsigned int i = 0; i < translationalState.size( ); i++ )
        {
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( stateIterator->second, combinedIterator->second.segment( 0, 6 ),
                                               std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( massIterator->second, combinedIterator->second.segment( 6, 1 ),
                                               std::numeric_limits< double >::epsilon( ) );
            stateIterator++;
            massIterator++;
            combinedIterator++;
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}

}
