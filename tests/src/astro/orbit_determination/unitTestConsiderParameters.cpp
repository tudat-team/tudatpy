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

#include <string>

#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/math/basic/linearAlgebra.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/orbitalElementConversions.h"

#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/math/integrators/rungeKuttaCoefficients.h"
#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/astro/basic_astro/keplerPropagator.h"
#include "tudat/simulation/simulation.h"
#include "tudat/simulation/propagation_setup/dynamicsSimulator.h"
#include "tudat/simulation/estimation.h"

namespace tudat {

namespace unit_tests {

//! Using declarations.
using namespace interpolators;
using namespace numerical_integrators;
using namespace spice_interface;
using namespace ephemerides;
using namespace simulation_setup;
using namespace basic_astrodynamics;
using namespace orbital_element_conversions;
using namespace propagators;
using namespace tudat::estimatable_parameters;

BOOST_AUTO_TEST_SUITE( test_consider_parameters )

BOOST_AUTO_TEST_CASE( testConsiderParameters )
{
    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Jupiter" );
    bodyNames.push_back( "Io" );
    bodyNames.push_back( "Europa" );
    bodyNames.push_back( "Ganymede" );
    bodyNames.push_back( "Callisto" );

    // Specify arc-wise start and end times
    double initialEpoch = 0.0;
    double finalEpoch = 2.0 * 86400.0;
    std::vector< double > arcStartTimes = { 0.0, 12.0 * 3600.0, 24.0 * 3600.0 };
    std::vector< double > arcEndTimes = { 12.0 * 3600.0, 24.0 * 3600.0, 36.0 * 3600.0 };
    std::vector< double > midArcTimes;
    for ( unsigned int i = 0 ; i < arcStartTimes.size( ) ; i++ )
    {
        midArcTimes.push_back( ( arcStartTimes.at( i ) + arcEndTimes.at( i ) ) / 2.0 );
    }

    unsigned int nbArcs = arcStartTimes.size( );

    std::string globalFrameOrientation = "ECLIPJ2000";
    std::string globalFrameOrigin = "Jupiter";

    // Create bodies.
    BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, initialEpoch - 86400.0, finalEpoch + 86400.0 );
    bodySettings.at( "Io" )->ephemerisSettings = std::make_shared< InterpolatedSpiceEphemerisSettings >(
            initialEpoch, finalEpoch, 3600.0, globalFrameOrigin, globalFrameOrientation );
    bodySettings.at( "Europa" )->ephemerisSettings = std::make_shared< InterpolatedSpiceEphemerisSettings >(
            initialEpoch, finalEpoch, 3600.0, globalFrameOrigin, globalFrameOrientation );
    bodySettings.at( "Ganymede" )->ephemerisSettings = std::make_shared< InterpolatedSpiceEphemerisSettings >(
            initialEpoch, finalEpoch, 3600.0, globalFrameOrigin, globalFrameOrientation );
    bodySettings.at( "Callisto" )->ephemerisSettings = std::make_shared< InterpolatedSpiceEphemerisSettings >(
            initialEpoch, finalEpoch, 3600.0, globalFrameOrigin, globalFrameOrientation );

    bodySettings.at( "Io" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );
    bodySettings.at( "Europa" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );
    bodySettings.at( "Ganymede" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );
    bodySettings.at( "Callisto" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    //! Define bodies to propagate.
    std::vector< std::string > bodiesToPropagate = { "Io", "Europa", "Ganymede", "Callisto" };
    std::vector< std::string > centralBodies = { "Jupiter", "Jupiter", "Jupiter", "Jupiter" };

    // Set accelerations.
    SelectedAccelerationMap accelerationSettings;
    for ( unsigned int i = 0 ; i < bodiesToPropagate.size( ) ; i++ )
    {
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfSatellite;
        accelerationsOfSatellite[ "Jupiter" ].push_back( std::make_shared< MutualSphericalHarmonicAccelerationSettings >( 8, 0, 2, 2  )  );
        for ( unsigned int j = 0 ; j < bodiesToPropagate.size( ) ; j++ )
        {
            if ( i != j )
            {
                accelerationsOfSatellite[ bodiesToPropagate[ j ] ].push_back( std::make_shared< MutualSphericalHarmonicAccelerationSettings >( 2, 2, 2, 2, 8, 0 ) );
            }
        }
        accelerationSettings[ bodiesToPropagate[ i ] ] = accelerationsOfSatellite;
    }

    basic_astrodynamics::AccelerationMap accelerationsMap = createAccelerationModelsMap(
            bodies, accelerationSettings, bodiesToPropagate, centralBodies );

    // Define integrator settings
    double timeStep = 60.0;
    std::shared_ptr< numerical_integrators::IntegratorSettings< > > integratorSettings =
            std::make_shared< numerical_integrators::RungeKuttaVariableStepSizeSettingsScalarTolerances< double > >(
                    TUDAT_NAN, timeStep, CoefficientSets::rungeKuttaFehlberg78, timeStep, timeStep, 1.0e3, 1.0e3 );

    // Define arc-wise initial states
    std::vector< Eigen::VectorXd > initialStatesMoons;
    Eigen::VectorXd initialStatesIo, initialStatesEuropa, initialStatesGanymede, initialStatesCallisto;
    initialStatesIo.resize( 6 * nbArcs );
    initialStatesEuropa.resize( 6 * nbArcs );
    initialStatesGanymede.resize( 6 * nbArcs );
    initialStatesCallisto.resize( 6 * nbArcs );
    for ( unsigned int i = 0 ; i < nbArcs ; i++ )
    {
        initialStatesMoons.push_back( propagators::getInitialStatesOfBodies( bodiesToPropagate, centralBodies, bodies, midArcTimes.at( i ) ) );
        initialStatesIo.segment( i * 6, 6 ) = initialStatesMoons[ i ].segment( 0, 6 );
        initialStatesEuropa.segment( i * 6, 6 ) = initialStatesMoons[ i ].segment( 6, 6 );
        initialStatesGanymede.segment( i * 6, 6 ) = initialStatesMoons[ i ].segment( 12, 6 );
        initialStatesCallisto.segment( i * 6, 6 ) = initialStatesMoons[ i ].segment( 18, 6 );
    }
    std::map< std::string, Eigen::VectorXd > initialStatesPerBody;
    initialStatesPerBody[ "Io" ] = initialStatesIo;
    initialStatesPerBody[ "Europa" ] = initialStatesEuropa;
    initialStatesPerBody[ "Ganymede" ] = initialStatesGanymede;
    initialStatesPerBody[ "Callisto" ] = initialStatesCallisto;

    // Define propagator settings.
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagationSettingsList;
    for ( unsigned int i = 0 ; i < nbArcs ; i++ )
    {
        std::shared_ptr< NonSequentialPropagationTerminationSettings > terminationSettings = std::make_shared< NonSequentialPropagationTerminationSettings >(
                std::make_shared< PropagationTimeTerminationSettings >( arcEndTimes.at( i ) ),
                std::make_shared< PropagationTimeTerminationSettings >( arcStartTimes.at( i ) ) );
        std::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< > >(
                centralBodies, accelerationsMap, bodiesToPropagate, initialStatesMoons.at( i ), midArcTimes.at( i ), integratorSettings, terminationSettings );
        propagationSettingsList.push_back( propagatorSettings );
    }
    std::shared_ptr< MultiArcPropagatorSettings< > > propagatorSettings = std::make_shared< MultiArcPropagatorSettings< > >( propagationSettingsList );


    // Define parameters to estimate for non-sequential propagation / estimation
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    for( unsigned int i = 0; i < bodiesToPropagate.size( ); i++ )
    {
        parameterNames.push_back( std::make_shared< ArcWiseInitialTranslationalStateEstimatableParameterSettings< double > >(
                bodiesToPropagate.at( i ), initialStatesPerBody.at( bodiesToPropagate.at( i ) ), midArcTimes, centralBodies.at( i ) ) );

        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( bodiesToPropagate.at( i ), gravitational_parameter ) );
    }

    std::vector< std::shared_ptr< EstimatableParameterSettings > > considerParameterNames;
    for( unsigned int i = 0; i < bodiesToPropagate.size( ); i++ )
    {
        considerParameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                2, 0,  2, 2, bodiesToPropagate.at( i ), spherical_harmonics_cosine_coefficient_block ) );
        considerParameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                2, 1,  2, 2, bodiesToPropagate.at( i ), spherical_harmonics_sine_coefficient_block ) );
    }
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parameters =
            createParametersToEstimate< double >( parameterNames, bodies, propagatorSettings, considerParameterNames );
    printEstimatableParameterEntries( parameters );

    std::cout << "size parameters to estimate: " << parameters->getFullParameterValues< double >( ).size( ) << "\n\n";

    // Define links and observations.
    std::vector< observation_models::LinkEnds > linkEndsList;
    std::vector< std::shared_ptr< observation_models::ObservationModelSettings > > observationSettingsList;
    linkEndsList.resize( bodiesToPropagate.size( ) );
    for( unsigned int i = 0; i < bodiesToPropagate.size( ) ; i++ )
    {
        linkEndsList[ i ][ observation_models::observed_body ] = observation_models::LinkEndId( bodiesToPropagate.at( i ), "" );
        observationSettingsList.push_back( std::make_shared< observation_models::ObservationModelSettings >(
                observation_models::position_observable, linkEndsList[ i ] ) );
    }

    // Define observation times
    std::vector< double > observationTimes;
    for ( unsigned int i = 0 ; i < nbArcs ; i++ )
    {
        for ( double time = arcStartTimes.at( i ) + 3600.0 ; time < arcEndTimes.at( i ) - 3600.0 ; time += 3.0 * 3600.0 )
        {
            observationTimes.push_back( time );
        }
    }

    // Define observation settings
    std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > measurementInput;
    for( unsigned int i = 0; i < bodiesToPropagate.size( ); i++ )
    {
        measurementInput.push_back( std::make_shared< TabulatedObservationSimulationSettings< > >(
                observation_models::position_observable, linkEndsList[ i ], observationTimes, observation_models::observed_body ) );
    }

    // Create orbit determination manager
    OrbitDeterminationManager< > orbitDeterminationManager = OrbitDeterminationManager< >(
            bodies, parameters, observationSettingsList, propagatorSettings );

    // Simulate observations
    std::shared_ptr< observation_models::ObservationCollection< > > observationsAndTimes = simulateObservations< >(
            measurementInput, orbitDeterminationManager.getObservationSimulators( ), bodies  );

    // Define estimation input
    std::shared_ptr< EstimationInput< double, double  > > estimationInput =
            std::make_shared< EstimationInput< double, double > >( observationsAndTimes, Eigen::MatrixXd::Zero( 0, 0 ), std::make_shared< EstimationConvergenceChecker >( 1 ) );

    // Perform estimation
    std::shared_ptr< EstimationOutput< double, double > > estimationOutput = orbitDeterminationManager.estimateParameters( estimationInput );

    // Retrieve covariance matrix
    Eigen::MatrixXd covariance = estimationOutput->getUnnormalizedCovarianceMatrix( );
    std::cout << "size cov: " << covariance.rows( ) << "\n\n";


}

BOOST_AUTO_TEST_SUITE_END( )

}

}