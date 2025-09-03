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

#include <iostream>
#include <ctime>

#include <boost/format.hpp>
#include <boost/test/unit_test.hpp>

#include "tudat/simulation/estimation_setup/fitOrbitToEphemeris.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"
#include "tudat/simulation/estimation_setup/observationSimulationSettings.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"

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

BOOST_AUTO_TEST_SUITE( test_estimation_drag_scaling )

BOOST_AUTO_TEST_CASE( test_EstimationDragScaling )
{
    
    std::vector< double > residuals;
    std::vector< double > parameterEstimateList;
    spice_interface::loadStandardSpiceKernels( );

    spice_interface::loadSpiceKernelInTudat( paths::getTudatTestDataPath( ) +
                                             "/dsn_n_way_doppler_observation_model/mgs_map1_ipng_mgs95j.bsp" );
    for( unsigned int i = 0; i < 5; i++ )
    {
        double initialTime = DateTime( 1999, 3, 10, 0, 0, 0.0 ).epoch< double >( );
        double finalTime = initialTime + 3600 * 6;

        std::vector< std::string > bodyNames;
        bodyNames.push_back( "Mars" );
        bodyNames.push_back( "Sun" );

        BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, "Mars" );
        bodySettings.addSettings( "MGS" );
        bodySettings.at( "MGS" )->ephemerisSettings =
            std::make_shared< InterpolatedSpiceEphemerisSettings >( initialTime - 3600.0, finalTime + 3600.0, 30.0, "Mars" );
        bodySettings.at( "MGS" )->radiationPressureTargetModelSettings = cannonballRadiationPressureTargetModelSettings(
                10.0, 1.2 );
        // Create bodies needed in simulation
        SystemOfBodies bodies = createSystemOfBodies( bodySettings );
        bodies.at( "Mars" )->setAtmosphereModel( 
                createAtmosphereModel( std::make_shared< ExponentialAtmosphereSettings >( aerodynamics::mars ),
                                       "Mars" ) );
        bodies.at( "MGS" )->setConstantBodyMass( 2.0 );

        Eigen::Vector3d forceCoefficients( 1.0, 0.0, 0.0 );

        std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
                std::make_shared< ConstantAerodynamicCoefficientSettings >(
                        2000.0, forceCoefficients, aerodynamics::negative_aerodynamic_frame_coefficients );
        bodies.at( "MGS" )
                ->setAerodynamicCoefficientInterface(
                        createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "MGS", bodies ) );

        // Set accelerations between bodies that are to be taken into account.
        SelectedAccelerationMap accelerationMap;
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfSpacecraft;
        accelerationsOfSpacecraft[ "Mars" ].push_back( std::make_shared< AccelerationSettings >( aerodynamic ) );
        accelerationsOfSpacecraft[ "Mars" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
        if( i == 4 )
        {
            accelerationsOfSpacecraft[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( radiation_pressure ) );
        }

        accelerationMap[ "MGS" ] = accelerationsOfSpacecraft;

        // Set bodies for which initial state is to be estimated and integrated.
        std::vector< std::string > bodiesToEstimate = { "MGS" };
        std::vector< std::string > centralBodies = { "Mars" };

        AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodies, accelerationMap, bodiesToEstimate, centralBodies );

        std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > > additionalParameterNames;

        if( i == 0 )
        {
            additionalParameterNames.push_back( estimatable_parameters::constantDragCoefficient( "MGS" ) );
        }
        else if( i == 1 )
        {
            additionalParameterNames.push_back( estimatable_parameters::dragComponentScaling( "MGS" ) );
        }
        else if( i == 2 )
        {
            additionalParameterNames.push_back( estimatable_parameters::fullAccelerationScaling(
                    "MGS", "Mars", basic_astrodynamics::aerodynamic ) );
        }
        else if( i == 3 || i == 4 )
        {
            additionalParameterNames.push_back( estimatable_parameters::areaToMassScaling(
                    "MGS" ) );
        }

        Eigen::Matrix< double, Eigen::Dynamic, 1 > initialState =
            getInitialStatesOfBodies( bodiesToEstimate, centralBodies, bodies, initialTime );

        std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;
        dependentVariablesList.push_back( singleAccelerationDependentVariable(
                basic_astrodynamics::aerodynamic, "MGS", "Mars" ) );
        if( i == 4 )
        {
            dependentVariablesList.push_back(
                    singleAccelerationDependentVariable( basic_astrodynamics::radiation_pressure, "MGS", "Sun" ) );
        }
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

        SingleArcDynamicsSimulator< > dynamicsSimulatorOriginal( bodies, propagatorSettings );

        std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
                getInitialStateParameterSettings< double, double >(
                    propagatorSettings,
                    bodies );
        parameterNames.insert( parameterNames.end( ), additionalParameterNames.begin( ), additionalParameterNames.end( ) );

        std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
                createParametersToEstimate< double, double >( parameterNames, bodies, propagatorSettings );
        printEstimatableParameterEntries( parametersToEstimate );

        std::pair< std::vector< std::shared_ptr< observation_models::ObservationModelSettings > >,
                   std::shared_ptr< observation_models::ObservationCollection< double > > >
                observationCollectionAndModelSettings = simulatePseudoObservations(
                        bodies, bodiesToEstimate, centralBodies, initialTime, finalTime, 120.0 );
        std::shared_ptr< observation_models::ObservationCollection< double > > observationCollection =
                observationCollectionAndModelSettings.second;

        std::vector< std::shared_ptr< observation_models::ObservationModelSettings > > observationModelSettingsList =
                observationCollectionAndModelSettings.first;

        OrbitDeterminationManager< double > orbitDeterminationManager =
                OrbitDeterminationManager< double >(
                        bodies, parametersToEstimate, observationModelSettingsList, propagatorSettings );

        std::shared_ptr< EstimationInput< double > > estimationInput =
                std::make_shared< EstimationInput< double > >( observationCollection );
        estimationInput->defineEstimationSettings( 0, 1, 0, 1, 1, 1 );

        std::shared_ptr< EstimationOutput<> > estimationOutput = orbitDeterminationManager.estimateParameters( estimationInput );

        residuals.push_back( estimationOutput->residualStandardDeviation_ );
        parameterEstimateList.push_back( estimationOutput->parameterHistory_.at(  estimationOutput->parameterHistory_.size( ) -1 )( 6 ) );

        double finalParameterValue = estimationOutput->parameterHistory_.at(  estimationOutput->parameterHistory_.size( ) -1 )( 6 );
        if( i == 4 )
        {
            Eigen::VectorXd currentParameters =
                    estimationOutput->parameterHistory_.at(  estimationOutput->parameterHistory_.size( ) -1 );
            currentParameters( 6 ) = 1.0;
            propagatorSettings->resetInitialStates( currentParameters.segment( 0, 6 ) );
            parametersToEstimate->resetParameterValues( currentParameters );
            SingleArcDynamicsSimulator< > dynamicsSimulatorRecomputed( bodies, propagatorSettings );


            int numberOfIterations = estimationOutput->getSimulationResults( ).size( );

            auto originalDependentVariableHistory =
                    dynamicsSimulatorOriginal.getSingleArcPropagationResults( )->getDependentVariableHistory( );
            auto recomputedDependentVariableHistory =
                    dynamicsSimulatorRecomputed.getSingleArcPropagationResults( )->getDependentVariableHistory( );
            auto estimatedDependentVariableHistory =
                    std::dynamic_pointer_cast< SingleArcVariationalSimulationResults< double, double > >(
                            estimationOutput->getSimulationResults( ).at( numberOfIterations - 1 ) )->getDynamicsResults( )
                            ->getDependentVariableHistory( );

            Eigen::VectorXd originalAccelerationRatios =
                    ( originalDependentVariableHistory.begin( )->second ).cwiseQuotient(
                            estimatedDependentVariableHistory.begin( )->second );
            Eigen::VectorXd recomputedAccelerationRatios =
                    ( recomputedDependentVariableHistory.begin( )->second ).cwiseQuotient(
                        estimatedDependentVariableHistory.begin( )->second );

            for( int i = 0; i < 6; i++ )
            {
                BOOST_CHECK_CLOSE_FRACTION( recomputedAccelerationRatios( i ), 1.0 / finalParameterValue, 1.0E-9 );
                if( i < 3 )
                {
                    BOOST_CHECK_CLOSE_FRACTION( originalAccelerationRatios( i ), 1.0 / finalParameterValue, 0.2 );
                }
                else
                {
                    BOOST_CHECK_CLOSE_FRACTION( originalAccelerationRatios( i ), 1.0 / finalParameterValue, 1.0E-6 );

                }
            }


        }
    }
    BOOST_CHECK( std::abs( residuals.at( 0 ) - residuals.at( 1 ) ) < 1e-8 );
    BOOST_CHECK_CLOSE_FRACTION( parameterEstimateList.at( 0 ), parameterEstimateList.at( 1 ), 1e-8 );
    BOOST_CHECK( std::abs( residuals.at( 0 ) - residuals.at( 2 ) ) < 1e-8 );
    BOOST_CHECK_CLOSE_FRACTION( parameterEstimateList.at( 0 ), parameterEstimateList.at( 2 ), 1e-8 );
    BOOST_CHECK( std::abs( residuals.at( 0 ) - residuals.at( 3 ) ) < 1e-8 );
    BOOST_CHECK_CLOSE_FRACTION( parameterEstimateList.at( 0 ), parameterEstimateList.at( 3 ), 1e-8 );
}
BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat