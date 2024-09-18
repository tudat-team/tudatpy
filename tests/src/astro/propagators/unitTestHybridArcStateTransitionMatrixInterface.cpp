/* git    Copyright (c) 2010-2019, Delft University of Technology
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

#include <string>
#include <thread>

#include <boost/test/unit_test.hpp>
#include "tudat/basics/testMacros.h"

#include "tudat/simulation/estimation.h"

#include "tudat/math/basic/linearAlgebra.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/unitConversions.h"

#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/math/integrators/rungeKuttaCoefficients.h"
#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/io/basicInputOutput.h"

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/estimation_setup/variationalEquationsSolver.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/createBodies.h"
#include "tudat/simulation/estimation_setup/createNumericalSimulator.h"
#include "tudat/simulation/estimation_setup/createEstimatableParameters.h"

#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameterSet.h"

namespace tudat
{

namespace unit_tests
{

// Using declarations.
using namespace tudat;
using namespace tudat::estimatable_parameters;
using namespace tudat::orbit_determination;
using namespace tudat::interpolators;
using namespace tudat::numerical_integrators;
using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;
using namespace tudat::basic_astrodynamics;
using namespace tudat::orbital_element_conversions;
using namespace tudat::ephemerides;
using namespace tudat::propagators;
using namespace tudat::observation_models;

BOOST_AUTO_TEST_SUITE( test_hybrid_arc_state_transition_matrix_interface )

BOOST_AUTO_TEST_CASE( testHybridArcStateTransitionMatrixInterface )
{
    std::cout.precision( 20 );

    //Load spice kernels
    spice_interface::loadStandardSpiceKernels( );

    // Simulation parameters
    double initialEpoch = 34.8 * physical_constants::JULIAN_YEAR;
    double finalEpoch = 34.9 * physical_constants::JULIAN_YEAR;

    std::vector< double > arcStartTime = { 34.8 * physical_constants::JULIAN_YEAR };
    std::vector< double > arcEndTime = { 34.8 * physical_constants::JULIAN_YEAR + 86400.0 };
    double testEpoch = arcStartTime.at( 0 ) + 1.0 * 3600.0;

    std::shared_ptr< IntegratorSettings< > > singleArcIntegratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettingsScalarTolerances< double > >(
            initialEpoch, 60.0, CoefficientSets::rungeKuttaFehlberg78, 60.0, 60.0, 1.0e3, 1.0e3 );

    std::shared_ptr< IntegratorSettings< > > multiArcIntegratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettingsScalarTolerances< double > >(
            TUDAT_NAN, 60.0, CoefficientSets::rungeKuttaFehlberg78, 60.0, 60.0, 1.0e3, 1.0e3 );


    // test case 0 : Ganymede as multi-arc body ; test case 1 : Ganymede as single-arc body
    for ( unsigned int testCase = 0 ; testCase < 2 ; testCase++ )
    {
        // Define bodies settings for simulation
        std::vector<std::string> bodiesToCreate = { "Ganymede", "Jupiter", "Sun" };
        BodyListSettings bodySettings = getDefaultBodySettings( bodiesToCreate, initialEpoch, finalEpoch, "Sun", "J2000" );

        // Create multi-arc ephemerides for Ganymede if necessary
        if ( testCase == 0 )
        {
            bodySettings.at( "Ganymede" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );
        }

//        bodySettings.addSettings( "-28" );
//        bodySettings.at( "-28" )->ephemerisSettings = std::make_shared< DirectSpiceEphemerisSettings >( "Jupiter", "J2000" );
//        bodySettings.at( "-28" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );

        Eigen::VectorXd unperturbedStateJupiter, unperturbedStateGanymede, unperturbedStateJuice, perturbedStateJupiter, perturbedStateGanymede, perturbedStateJuice;
        Eigen::MatrixXd stateTransitionSensitivityMatrix;
        for ( unsigned int perturbed = 0 ; perturbed < 2 ; perturbed++ )
        {
            // Create systemOfBodies
            SystemOfBodies bodies = createSystemOfBodies( bodySettings );
            bodies.createEmptyBody( "JUICE" );
            bodies.at( "JUICE" )->setEphemeris( std::make_shared< ephemerides::MultiArcEphemeris >(
                    std::map< double, std::shared_ptr< ephemerides::Ephemeris > >( ), "Ganymede", "J2000" ) );

            // Define acceleration
            SelectedAccelerationMap accelerationSettings;
            accelerationSettings[ "Jupiter" ][ "Sun" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );
            accelerationSettings[ "Ganymede" ][ "Jupiter" ].push_back( std::make_shared< MutualSphericalHarmonicAccelerationSettings >( 2, 2, 2, 2 ) );
            accelerationSettings[ "JUICE" ][ "Ganymede" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );

            // Define single-arc and multi-arc propagation settings
            std::vector< std::string > singleArcBodiesToPropagate, singleArcCentralBodies, multiArcBodiesToPropagate, multiArcCentralBodies;
            SelectedAccelerationMap singleArcAccelerationSettings, multiArcAccelerationSettings;
            singleArcAccelerationSettings[ "Jupiter" ] = accelerationSettings.at( "Jupiter" );
            multiArcAccelerationSettings[ "JUICE" ] = accelerationSettings.at( "JUICE" );
            if ( testCase == 0 )
            {
                singleArcBodiesToPropagate = { "Jupiter" };
                singleArcCentralBodies = { "Sun" };
                multiArcBodiesToPropagate = { "Ganymede", "JUICE" };
                multiArcCentralBodies = { "Jupiter", "Ganymede" };

                multiArcAccelerationSettings[ "Ganymede" ] = accelerationSettings.at( "Ganymede" );

            }
            else
            {
                singleArcBodiesToPropagate = { "Jupiter", "Ganymede" };
                singleArcCentralBodies = { "Sun", "Jupiter" };
                multiArcBodiesToPropagate = { "JUICE" };
                multiArcCentralBodies = { "Ganymede" };

                singleArcAccelerationSettings[ "Ganymede" ] = accelerationSettings.at( "Ganymede" );
            }

            // Define initial states
            Eigen::Vector6d initialStateJupiter = spice_interface::getBodyCartesianStateAtEpoch( "Jupiter", "Sun", "J2000", "None", initialEpoch );
            Eigen::Vector6d initialStateGanymede;
            if ( testCase == 0 )
            {
                initialStateGanymede = spice_interface::getBodyCartesianStateAtEpoch( "Ganymede", "Jupiter", "J2000", "None", initialEpoch );
            }
            else
            {
                initialStateGanymede = spice_interface::getBodyCartesianStateAtEpoch( "Ganymede", "Jupiter", "J2000", "None", arcStartTime.at( 0 ) );
            }
            Eigen::Vector6d keplerianInitialStateJuice = ( Eigen::Vector6d( ) << 500.0e3, 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( );
            Eigen::Vector6d initialStateJuice = convertKeplerianToCartesianElements( keplerianInitialStateJuice, bodies.at( "Ganymede" ) ->getGravitationalParameter( ) );

            Eigen::VectorXd singleArcInitialStates, multiArcInitialStates;
            if ( testCase == 0 )
            {
                singleArcInitialStates = initialStateJupiter;
                multiArcInitialStates = Eigen::VectorXd::Zero( 12 );
                multiArcInitialStates.segment( 0, 6 ) = initialStateGanymede;
                multiArcInitialStates.segment( 6, 6 ) = initialStateJuice;
            }
            else
            {
                singleArcInitialStates = Eigen::VectorXd::Zero( 12 );
                singleArcInitialStates.segment( 0, 6 ) = initialStateJupiter;
                singleArcInitialStates.segment( 6, 6 ) = initialStateGanymede;
                multiArcInitialStates = initialStateJuice;
            }

            Eigen::VectorXd fullInitialState = Eigen::VectorXd::Zero( 18 );
            fullInitialState.segment( 0, 6 ) = initialStateJupiter;
            fullInitialState.segment( 6, 6 ) = initialStateGanymede;
            fullInitialState.segment( 12, 6 ) = initialStateJuice;


            // Create single-arc propagator settings
            basic_astrodynamics::AccelerationMap singleArcAccelerations = createAccelerationModelsMap(
                    bodies, singleArcAccelerationSettings, singleArcBodiesToPropagate, singleArcCentralBodies );

            std::shared_ptr< TranslationalStatePropagatorSettings< > > singleArcPropagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< > >(
                    singleArcCentralBodies, singleArcAccelerations, singleArcBodiesToPropagate, singleArcInitialStates, initialEpoch,
                    singleArcIntegratorSettings, std::make_shared< PropagationTimeTerminationSettings >( finalEpoch ) );

            // Create multi-arc propagator settings
            basic_astrodynamics::AccelerationMap multiArcAccelerations = createAccelerationModelsMap(
                    bodies, multiArcAccelerationSettings, multiArcBodiesToPropagate, multiArcCentralBodies );

            std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsList;
            propagatorSettingsList.push_back( std::make_shared< TranslationalStatePropagatorSettings< double > >(
                    multiArcCentralBodies, multiArcAccelerations, multiArcBodiesToPropagate, multiArcInitialStates, arcStartTime.at( 0 ),
                    multiArcIntegratorSettings, std::make_shared< PropagationTimeTerminationSettings >( arcEndTime.at( 0 ) ) ) );


            if ( perturbed == 0 )
            {
                std::shared_ptr< HybridArcPropagatorSettings< > > hybridArcPropagatorSettings = std::make_shared< HybridArcPropagatorSettings< > >(
                        singleArcPropagatorSettings, std::make_shared< MultiArcPropagatorSettings< > >( propagatorSettingsList ) );
                HybridArcDynamicsSimulator< > dynamicsSimulator = HybridArcDynamicsSimulator< >( bodies, hybridArcPropagatorSettings );

                // Create state history interpolators
                std::map< double, Eigen::VectorXd > singleArcStateHistory = dynamicsSimulator.getSingleArcDynamicsSimulator( )->getEquationsOfMotionNumericalSolution( );
                std::map< double, Eigen::VectorXd > multiArcStateHistory = dynamicsSimulator.getMultiArcDynamicsSimulator( )->getEquationsOfMotionNumericalSolution( ).at( 0 );

                std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > singleArcInterpolator =
                        std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd > >(
                                utilities::createVectorFromMapKeys< Eigen::VectorXd, double >( singleArcStateHistory ),
                                utilities::createVectorFromMapValues< Eigen::VectorXd, double >( singleArcStateHistory ), 8 );

                std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > multiArcInterpolator =
                        std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd > >(
                                utilities::createVectorFromMapKeys< Eigen::VectorXd, double >( multiArcStateHistory ),
                                utilities::createVectorFromMapValues< Eigen::VectorXd, double >( multiArcStateHistory ), 8 );

                // Retrieve unperturbed states at test epoch
                unperturbedStateJupiter = singleArcInterpolator->interpolate( testEpoch ).segment( 0, 6 );
                if ( testCase == 0 )
                {
                    unperturbedStateGanymede = multiArcInterpolator->interpolate( testEpoch ).segment( 0, 6 ) + unperturbedStateJupiter;
                    unperturbedStateJuice = multiArcInterpolator->interpolate( testEpoch ).segment( 6, 6 ) + unperturbedStateGanymede;
                }
                else
                {
                    unperturbedStateGanymede = singleArcInterpolator->interpolate( testEpoch ).segment( 6, 6 ) + unperturbedStateJupiter;
                    unperturbedStateJuice = multiArcInterpolator->interpolate( testEpoch ).segment( 0, 6 ) + unperturbedStateGanymede;
                }

                std::cout << "unperturbedStateJupiter: " << unperturbedStateJupiter.transpose( ) << "\n\n";
                std::cout << "unperturbedStateGanymede: " << unperturbedStateGanymede.transpose( ) << "\n\n";
                std::cout << "unperturbedStateJuice: " << unperturbedStateJuice.transpose( ) << "\n\n";


                // Create parameters to estimate
                std::vector< std::shared_ptr< EstimatableParameterSettings > > parametersNames;
                parametersNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
                        "Jupiter", initialStateJupiter, "Sun", "J2000" ) );
                if ( testCase == 0 )
                {
                    parametersNames.push_back( std::make_shared< ArcWiseInitialTranslationalStateEstimatableParameterSettings< double > >(
                            "Ganymede", initialStateGanymede, arcStartTime, "Jupiter", "J2000" ) );
                }
                else
                {
                    parametersNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
                            "Ganymede", initialStateGanymede, "Jupiter", "J2000" ) );
                }
                parametersNames.push_back( std::make_shared< ArcWiseInitialTranslationalStateEstimatableParameterSettings< double > >(
                        "JUICE", initialStateJuice, arcStartTime, "Ganymede", "J2000" ) );
                parametersNames.push_back( std::make_shared< EstimatableParameterSettings >( "Jupiter", gravitational_parameter ) );

                std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
                        createParametersToEstimate< double >( parametersNames, bodies, hybridArcPropagatorSettings );
                printEstimatableParameterEntries( parametersToEstimate );

                // Dummy observations list
                LinkEnds dummyLinkEnd;
                dummyLinkEnd[ observed_body ] = std::make_pair< std::string, std::string >( "Jupiter" , ""  );
                std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;
                observationSettingsList.push_back( std::make_shared< ObservationModelSettings >( position_observable, dummyLinkEnd, std::shared_ptr< LightTimeCorrectionSettings >( ) ) );

                // Create OD manager to later retrieve state transition and sensitivity matrix interface
                OrbitDeterminationManager< > orbitDeterminationManager = OrbitDeterminationManager< >(
                        bodies, parametersToEstimate, observationSettingsList, hybridArcPropagatorSettings, true );

                std::shared_ptr< HybridArcVariationalEquationsSolver< > > variationalEquationsSolver =
                        std::dynamic_pointer_cast< HybridArcVariationalEquationsSolver< > >( orbitDeterminationManager.getVariationalEquationsSolver( ) );

                // Retrieve state transition matrix at test epoch
                stateTransitionSensitivityMatrix = orbitDeterminationManager.getStateTransitionAndSensitivityMatrixInterface( )->
                        getFullCombinedStateTransitionAndSensitivityMatrix( testEpoch );
            }
            else
            {
                // perturbed parameter = 0 -> Jupiter's state
                // perturbed parameter = 1 -> Ganymede's state
                // perturbed parameter = 2 -> JUICE's state
                // perturbed parameter = 3 -> Jupiter's gravitational parameter

                double muJupiter = bodies.at( "Jupiter" )->getGravitationalParameter( );
                for ( unsigned int perturbedParameter = 0 ; perturbedParameter < 4 ; perturbedParameter++ )
                {
                    Eigen::Vector6d perturbation = Eigen::Vector6d::Zero( );
                    Eigen::VectorXd fullPerturbedParameters = Eigen::VectorXd::Zero( 19 );

                    if ( perturbedParameter < 3 )
                    {
                        if ( perturbedParameter == 0 || perturbedParameter == 1 )
                        {
                            perturbation = ( Eigen::Vector6d( ) << 10.0, 10.0, 10.0, 0.1, 0.1, 0.1 ).finished( );
                        }
                        else
                        {
                            perturbation = ( Eigen::Vector6d( ) << 1.0, 1.0, 1.0, 0.001, 0.001, 0.001 ).finished( );
                        }

                        Eigen::VectorXd singleArcPerturbedInitialState = singleArcInitialStates;
                        Eigen::VectorXd multiArcPerturbedInitialState = multiArcInitialStates;
                        if ( perturbedParameter == 0 ) // perturbing Jupiter's state
                        {
                            singleArcPerturbedInitialState.segment( 0, 6 ) += perturbation;
                            singleArcPropagatorSettings->resetInitialStates( singleArcPerturbedInitialState );
                            fullPerturbedParameters.segment( 0, 6 ) = perturbation;
                        }
                        else if ( perturbedParameter == 1 ) // perturbing Ganymede's state
                        {
                            if ( testCase == 0 )
                            {
                                multiArcPerturbedInitialState.segment( 0, 6 ) += perturbation;
                                propagatorSettingsList.at( 0 )->resetInitialStates( multiArcPerturbedInitialState );
                            }
                            else
                            {
                                singleArcPerturbedInitialState.segment( 6, 6 ) += perturbation;
                                singleArcPropagatorSettings->resetInitialStates( singleArcPerturbedInitialState );
                            }
                            fullPerturbedParameters.segment( 6, 6 ) = perturbation;
                        }
                        else if ( perturbedParameter == 2 ) // perturbing JUICE's state
                        {
                            if ( testCase == 0 )
                            {
                                multiArcPerturbedInitialState.segment( 6, 6 ) += perturbation;
                                propagatorSettingsList.at( 0 )->resetInitialStates( multiArcPerturbedInitialState );
                            }
                            else
                            {
                                multiArcPerturbedInitialState.segment( 0, 6 ) += perturbation;
                                propagatorSettingsList.at( 0 )->resetInitialStates( multiArcPerturbedInitialState );
                            }
                            fullPerturbedParameters.segment( 12, 6 ) = perturbation;
                        }
                        else if ( perturbedParameter == 3 ) // Perturbing Jupiter's gravitational parameter
                        {
                            bodies.at( "Jupiter" )->getGravityFieldModel( )->resetGravitationalParameter( muJupiter * 1.0001 );
                            fullPerturbedParameters[ 18 ] = 0.0001 * muJupiter;
                        }
                    }

                    // Create new propagator settings
                    std::shared_ptr< HybridArcPropagatorSettings< > > perturbedHybridArcPropagatorSettings = std::make_shared< HybridArcPropagatorSettings< > >(
                            singleArcPropagatorSettings, std::make_shared< MultiArcPropagatorSettings< > >( propagatorSettingsList ) );
                    HybridArcDynamicsSimulator< > perturbedDynamicsSimulator = HybridArcDynamicsSimulator< >( bodies, perturbedHybridArcPropagatorSettings );

                    // Create state history interpolators
                    std::map< double, Eigen::VectorXd > perturbedSingleArcStateHistory =
                            perturbedDynamicsSimulator.getSingleArcDynamicsSimulator( )->getEquationsOfMotionNumericalSolution( );
                    std::map< double, Eigen::VectorXd > perturbedMultiArcStateHistory =
                            perturbedDynamicsSimulator.getMultiArcDynamicsSimulator( )->getEquationsOfMotionNumericalSolution( ).at( 0 );

                    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > perturbedSingleArcInterpolator =
                            std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd > >(
                                    utilities::createVectorFromMapKeys< Eigen::VectorXd, double >( perturbedSingleArcStateHistory ),
                                    utilities::createVectorFromMapValues< Eigen::VectorXd, double >( perturbedSingleArcStateHistory ), 8 );

                    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > perturbedMultiArcInterpolator =
                            std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd > >(
                                    utilities::createVectorFromMapKeys< Eigen::VectorXd, double >( perturbedMultiArcStateHistory ),
                                    utilities::createVectorFromMapValues< Eigen::VectorXd, double >( perturbedMultiArcStateHistory ), 8 );


                    // Retrieve perturbed states at test epoch
                    perturbedStateJupiter = perturbedSingleArcInterpolator->interpolate( testEpoch ).segment( 0, 6 );
                    if ( testCase == 0 )
                    {
                        perturbedStateGanymede = perturbedMultiArcInterpolator->interpolate( testEpoch ).segment( 0, 6 ) + perturbedStateJupiter;
                        perturbedStateJuice = perturbedMultiArcInterpolator->interpolate( testEpoch ).segment( 6, 6 ) + perturbedStateGanymede;
                    }
                    else
                    {
                        perturbedStateGanymede = perturbedSingleArcInterpolator->interpolate( testEpoch ).segment( 6, 6 ) + perturbedStateJupiter;
                        perturbedStateJuice = perturbedMultiArcInterpolator->interpolate( testEpoch ).segment( 0, 6 ) + perturbedStateGanymede;
                    }


                    Eigen::VectorXd predictedChangeState = stateTransitionSensitivityMatrix * fullPerturbedParameters;
                    std::cout << "expected change: " << predictedChangeState.transpose( ) << "\n\n";

                    Eigen::VectorXd changeStateJupiter = perturbedStateJupiter - unperturbedStateJupiter;
                    Eigen::VectorXd changeStateGanymede = perturbedStateGanymede - unperturbedStateGanymede;
                    Eigen::VectorXd changeStateJuice = perturbedStateJuice - unperturbedStateJuice;

                    std::cout << "changeStateJupiter: " << changeStateJupiter.transpose( ) << "\n\n";
                    std::cout << "changeStateGanymede: " << changeStateGanymede.transpose( ) << "\n\n";
                    std::cout << "changeStateJuice: " << changeStateJuice.transpose( ) << "\n\n";

                    if ( perturbedParameter < 3 )
                    {
                        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( predictedChangeState.segment( 0, 6 ), changeStateJupiter, 5.0E-4 );
                    }
                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( predictedChangeState.segment( 6, 6 ), changeStateGanymede, 5.0E-4 );
                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( predictedChangeState.segment( 12, 6 ), changeStateJuice, 5.0E-4 );

                    // Reset original parameters values
                    singleArcPropagatorSettings->resetInitialStates( singleArcInitialStates );
                    propagatorSettingsList.at( 0 )->resetInitialStates( multiArcInitialStates );
                    bodies.at( "Jupiter" )->getGravityFieldModel( )->resetGravitationalParameter( muJupiter );

                }
            }
        }
    }
}


BOOST_AUTO_TEST_SUITE_END( )

}

}
