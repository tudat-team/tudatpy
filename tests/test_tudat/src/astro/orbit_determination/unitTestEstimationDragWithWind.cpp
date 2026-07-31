/*    Copyright (c) 2010-2019, Delft University of Technology
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

#include <cmath>
#include <iostream>

#include <boost/test/unit_test.hpp>

#include "tudat/simulation/environment_setup/createAtmosphereModel.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/estimation_setup/createEstimatableParametersFactory.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"
#include "tudat/simulation/estimation_setup/simulatePseudoObservations.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"
#include "tudat/simulation/propagation_setup/singleArcDynamicsSimulator.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat;
using namespace tudat::observation_models;
using namespace tudat::orbit_determination;
using namespace tudat::estimatable_parameters;
using namespace tudat::simulation_setup;
using namespace tudat::basic_astrodynamics;
using namespace tudat::propagators;

BOOST_AUTO_TEST_SUITE( test_estimation_drag_with_wind )

// Regression test for the NaN bug in variational-equations propagation when a wind model is
// active together with aerodynamic drag and Cd is estimated. The bug surfaces in Python as:
//
//     RuntimeError: Error: Lagrange interpolator cannot identify zero entry.
//
// thrown from lagrangeInterpolator.h:131 — the integrated state map fed to the interpolator
// contains NaN entries because AerodynamicAccelerationPartial::update() produces NaN partials
// when wind is present.
//
// This test mirrors the structure of unitTestEstimationDragScaling.cpp::test_EstimationDragScaling
// but adds a wind model. It uses ExponentialAtmosphere (no coma — that's covered by a separate
// test) so it isolates the aerodynamic-partial + wind interaction.
BOOST_AUTO_TEST_CASE( test_VariationalEquationsWithWindModel )
{
    spice_interface::loadStandardSpiceKernels( );

    spice_interface::loadSpiceKernelInTudat( paths::getTudatTestDataPath( ) +
                                             "/dsn_n_way_doppler_observation_model/mgs_map1_ipng_mgs95j.bsp" );

    // Sweep wind reference frame: corotating and vertical (both are used in production setups).
    for( int frameIndex = 0; frameIndex < 2; ++frameIndex )
    {
        const reference_frames::AerodynamicsReferenceFrames windFrame =
                ( frameIndex == 0 ) ? reference_frames::corotating_frame : reference_frames::vertical_frame;
        BOOST_TEST_MESSAGE( "Sub-case: wind frame = " << ( frameIndex == 0 ? "corotating" : "vertical" ) );

        const double initialTime = DateTime( 1999, 3, 10, 0, 0, 0.0 ).epoch< double >( );
        const double finalTime = initialTime + 3600 * 6;

        std::vector< std::string > bodyNames = { "Mars", "Sun" };
        BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, "Mars" );
        bodySettings.addSettings( "MGS" );
        bodySettings.at( "MGS" )->ephemerisSettings =
                std::make_shared< InterpolatedSpiceEphemerisSettings >( initialTime - 3600.0, finalTime + 3600.0, 30.0, "Mars" );
        bodySettings.at( "MGS" )->radiationPressureTargetModelSettings = cannonballRadiationPressureTargetModelSettings( 10.0, 1.2 );

        // Mars atmosphere with a constant wind. The wind is non-zero so the aerodynamic partial
        // exercises the same code path as the failing Rosetta cases.
        const Eigen::Vector3d constantWindInFrame( 100.0, -50.0, 200.0 );
        const auto windFunction = [ constantWindInFrame ]( const double, const double, const double, const double ) {
            return constantWindInFrame;
        };
        bodySettings.at( "Mars" )->atmosphereSettings = std::make_shared< ExponentialAtmosphereSettings >( aerodynamics::mars );
        bodySettings.at( "Mars" )->atmosphereSettings->setWindSettings(
                std::make_shared< CustomWindModelSettings >( windFunction, windFrame ) );

        SystemOfBodies bodies = createSystemOfBodies( bodySettings );
        bodies.at( "MGS" )->setConstantBodyMass( 2.0 );

        // Aerodynamic coefficients.
        const Eigen::Vector3d forceCoefficients( 1.0, 0.0, 0.0 );
        std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
                std::make_shared< ConstantAerodynamicCoefficientSettings >(
                        2000.0, forceCoefficients, aerodynamics::negative_aerodynamic_frame_coefficients );
        bodies.at( "MGS" )->setAerodynamicCoefficientInterface(
                createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "MGS", bodies ) );

        // Acceleration setup: aerodynamic + point-mass gravity (drag is what exercises the bug).
        SelectedAccelerationMap accelerationMap;
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfSpacecraft;
        accelerationsOfSpacecraft[ "Mars" ].push_back( std::make_shared< AccelerationSettings >( aerodynamic ) );
        accelerationsOfSpacecraft[ "Mars" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
        accelerationMap[ "MGS" ] = accelerationsOfSpacecraft;

        std::vector< std::string > bodiesToEstimate = { "MGS" };
        std::vector< std::string > centralBodies = { "Mars" };
        AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodies, accelerationMap, bodiesToEstimate, centralBodies );

        // Estimate initial state + Cd — the exact parameter combination that triggers the bug.
        std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > > additionalParameterNames;
        additionalParameterNames.push_back( estimatable_parameters::constantDragCoefficient( "MGS" ) );

        const Eigen::Matrix< double, Eigen::Dynamic, 1 > initialState =
                getInitialStatesOfBodies( bodiesToEstimate, centralBodies, bodies, initialTime );

        std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;
        dependentVariablesList.push_back( singleAccelerationDependentVariable( basic_astrodynamics::aerodynamic, "MGS", "Mars" ) );

        std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >(
                        centralBodies,
                        accelerationModelMap,
                        bodiesToEstimate,
                        initialState,
                        initialTime,
                        numerical_integrators::rungeKuttaFixedStepSettings( 120.0, numerical_integrators::rungeKuttaFehlberg78 ),
                        std::make_shared< PropagationTimeTerminationSettings >( finalTime ),
                        cowell,
                        dependentVariablesList );

        std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
                getInitialStateParameterSettings< double, double >( propagatorSettings, bodies );
        parameterNames.insert( parameterNames.end( ), additionalParameterNames.begin( ), additionalParameterNames.end( ) );

        std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
                createParametersToEstimate< double, double >( parameterNames, bodies, propagatorSettings );

        // Pseudo-observations — needed by OrbitDeterminationManager constructor.
        std::pair< std::vector< std::shared_ptr< observation_models::ObservationModelSettings > >,
                   std::shared_ptr< observation_models::ObservationCollection< double > > >
                observationCollectionAndModelSettings =
                        simulatePseudoObservations( bodies, bodiesToEstimate, centralBodies, initialTime, finalTime, 120.0 );
        const std::vector< std::shared_ptr< observation_models::ObservationModelSettings > > observationModelSettingsList =
                observationCollectionAndModelSettings.first;

        // === The actual bug surface. ===
        // The OrbitDeterminationManager constructor propagates the variational equations together
        // with the dynamics. If AerodynamicAccelerationPartial::update() produces NaN partials,
        // the integrated variational state has NaN entries; the subsequent Lagrange interpolator
        // construction (state ephemeris) throws with the message the user sees in Python.
        std::shared_ptr< OrbitDeterminationManager< double > > orbitDeterminationManager;
        BOOST_REQUIRE_NO_THROW( orbitDeterminationManager = std::make_shared< OrbitDeterminationManager< double > >(
                                        bodies, parametersToEstimate, observationModelSettingsList, propagatorSettings ) );
        BOOST_REQUIRE( orbitDeterminationManager != nullptr );

        // Probe the state-transition + sensitivity interface at a few epochs. If the
        // interpolators inside it were built from NaN-laden matrices, query results are NaN —
        // catches the case where the upstream catch block silently swallowed the Lagrange
        // exception and left the interface populated with null interpolators that produce NaN
        // when queried.
        auto stateTransitionInterface = orbitDeterminationManager->getStateTransitionAndSensitivityMatrixInterface( );
        BOOST_REQUIRE( stateTransitionInterface != nullptr );

        const std::vector< double > probeTimes = { initialTime + 600.0, initialTime + 1800.0, initialTime + 3600.0, initialTime + 5400.0 };
        for( const double t : probeTimes )
        {
            Eigen::MatrixXd combined;
            BOOST_REQUIRE_NO_THROW( combined = stateTransitionInterface->getCombinedStateTransitionAndSensitivityMatrix( t ) );
            for( int row = 0; row < combined.rows( ); ++row )
            {
                for( int col = 0; col < combined.cols( ); ++col )
                {
                    BOOST_CHECK_MESSAGE( std::isfinite( combined( row, col ) ),
                                         "combined state transition + sensitivity matrix has non-finite entry at t="
                                                 << t << " (" << row << "," << col << ") = " << combined( row, col ) );
                }
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
