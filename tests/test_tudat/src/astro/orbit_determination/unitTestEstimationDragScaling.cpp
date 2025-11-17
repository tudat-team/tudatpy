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
#include "tudat/simulation/propagation_setup/setNumericallyIntegratedStates.h"
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

        std::cout << "i: " << i << std::endl;

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


// unit test for testing of recovery of original arc-wise aerodynamic parameters (in conjunction with states parameters)
BOOST_AUTO_TEST_CASE( test_EstimationArcwiseDragScaling )
{
    std::vector< double > residuals;
    std::vector< double > parameterEstimateList;
    spice_interface::loadStandardSpiceKernels( );

    spice_interface::loadSpiceKernelInTudat( paths::getTudatTestDataPath( ) +
                                             "/dsn_n_way_doppler_observation_model/mgs_map1_ipng_mgs95j.bsp" );

    std::cout << "start "  << std::endl;

    for( unsigned int i = 0; i < 3; i++ )
    {
        std::cout << "i: " << i << std::endl;

        double initialTime = DateTime( 1999, 3, 10, 0, 0, 0.0 ).epoch< double >( );
        double finalTime = initialTime + 3600 * 10;

        std::cout << "arcStartTimes " << std::endl;

        std::vector<double> arcStartTimes = {initialTime, initialTime+4*3600, initialTime+7.5*3600};

        std::vector< std::string > bodyNames;
        bodyNames.push_back( "Mars" );
        bodyNames.push_back( "Sun" );

        BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, "Mars" );

        std::map<double, Eigen::Vector6d> emptyMap;

        bodySettings.addSettings( "MGS" );
        bodySettings.at( "MGS" )->ephemerisSettings =
            std::make_shared< InterpolatedSpiceEphemerisSettings >( initialTime - 3600.0, finalTime + 3600.0, 30.0, "Mars" );
        //bodySettings.at("MGS")->ephemerisSettings = tabulatedEphemerisSettings(emptyMap, "Mars", "ECLIPJ2000");
        bodySettings.at( "MGS" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );

        bodySettings.at( "MGS" )->radiationPressureTargetModelSettings = cannonballRadiationPressureTargetModelSettings(
                10.0, 1.2 );


        // Create bodies needed in simulation
        SystemOfBodies bodies = createSystemOfBodies( bodySettings );

        //bodies.createEmptyBody( "MGS" );
        //addEmptyTabulatedEphemeris<double, double>(bodies, "MGS", "", true);
        bodies.at( "MGS" )->setConstantBodyMass( 2.0 );


        Eigen::Vector3d forceCoefficients( 1.0, 0.0, 1.0 );
        std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
                std::make_shared< ConstantAerodynamicCoefficientSettings >(
                        2000.0, forceCoefficients, aerodynamics::negative_aerodynamic_frame_coefficients );
        bodies.at( "MGS" )
                ->setAerodynamicCoefficientInterface(
                        createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "MGS", bodies ) );


        bodies.at( "Mars" )->setAtmosphereModel(
                createAtmosphereModel( std::make_shared< ExponentialAtmosphereSettings >( aerodynamics::mars ),
                                       "Mars" ) );


        // Set accelerations between bodies that are to be taken into account.
        SelectedAccelerationMap accelerationMap;
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfSpacecraft;
        accelerationsOfSpacecraft[ "Mars" ].push_back( std::make_shared< AccelerationSettings >( aerodynamic ) );
        accelerationsOfSpacecraft[ "Mars" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
        accelerationsOfSpacecraft[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( radiation_pressure ) );

        accelerationMap[ "MGS" ] = accelerationsOfSpacecraft;

        // Set bodies for which initial state is to be estimated and integrated.
        std::vector< std::string > bodiesToEstimate = { "MGS" };
        std::vector< std::string > centralBodies = { "Mars" };

        AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodies, accelerationMap, bodiesToEstimate, centralBodies );

        std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > > additionalParameterNames;

        if( i == 0 )
        {
            additionalParameterNames.push_back( estimatable_parameters::arcwiseDragCoefficient( "MGS",  arcStartTimes) );
        }
        if( i == 1 )
        {
            additionalParameterNames.push_back( estimatable_parameters::arcwiseDragComponentScaling( "MGS", arcStartTimes ) );
        }
        if( i == 2 )
        {
            additionalParameterNames.push_back( estimatable_parameters::arcwiseDragComponentScaling( "MGS", arcStartTimes ) );
            additionalParameterNames.push_back( estimatable_parameters::arcwiseLiftComponentScaling( "MGS", arcStartTimes ) );
        }

        std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;
        dependentVariablesList.push_back( singleAccelerationDependentVariable(
                basic_astrodynamics::aerodynamic, "MGS", "Mars" ) );

        std::vector< double > integrationArcStartTimes = arcStartTimes;
        std::vector< double > integrationArcEndTimes(arcStartTimes.begin() + 1, arcStartTimes.end());
        for (int i = 0; i < integrationArcEndTimes.size( ); ++i) {
            integrationArcEndTimes.at(i) = integrationArcEndTimes.at(i) - 200*(i+1);
        }
        integrationArcEndTimes.push_back( finalTime );

        std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsList;
        for( unsigned int j = 0; j < integrationArcStartTimes.size( ); j++ )
        {
            std::cout << "Querying initial states " << std::endl;
            Eigen::Matrix< double, Eigen::Dynamic, 1 > currentInitialState =
                spice_interface::getBodyCartesianStateAtEpoch( "MGS", "Mars", "J2000", "None", integrationArcStartTimes.at(j));
            std::cout << "Got: " << currentInitialState << std::endl;

            propagatorSettingsList.push_back(
                std::make_shared< TranslationalStatePropagatorSettings< double > >(
                    centralBodies,
                    accelerationModelMap,
                    bodiesToEstimate,
                    currentInitialState,
                    integrationArcStartTimes.at( j ),
                    numerical_integrators::rungeKuttaFixedStepSettings( 30.0, numerical_integrators::rungeKuttaFehlberg78 ),
                    propagationTimeTerminationSettings( integrationArcEndTimes.at( j ) ),
                    cowell,
                    dependentVariablesList ) );

        }

        std::cout << "MultiArcPropagatorSettings" << std::endl;
        std::shared_ptr< MultiArcPropagatorSettings< double > > propagatorSettings =
                std::make_shared< MultiArcPropagatorSettings< double > >( propagatorSettingsList,
        false, std::make_shared< MultiArcPropagatorProcessingSettings >(false, true));


        std::cout << "************************************* RUNNING TEST *************************" << std::endl;
        std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
                getInitialMultiArcParameterSettings< double >( propagatorSettings, bodies, integrationArcStartTimes );

        parameterNames.insert( parameterNames.end( ), additionalParameterNames.begin( ), additionalParameterNames.end( ) );

        std::cout << "EstimatableParameterSet" << std::endl;
        std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
                createParametersToEstimate< double, double >( parameterNames, bodies, propagatorSettings );
        printEstimatableParameterEntries( parametersToEstimate );

        std::vector< double > observationTimes;
        double dataPointInterval = 120.;
        for( int k = 0; k < integrationArcStartTimes.size( ); k++)
        {
            double currentTime = integrationArcStartTimes.at( k ) + 300.0;
            while( currentTime < integrationArcEndTimes.at( k ) - 300.0 )

            {
                observationTimes.push_back( currentTime );
                currentTime += static_cast< double >( dataPointInterval );
            }
        }

        std::cout << "Integrate Equations of Motion" << std::endl;
        // Create simulation object (but do not propagate dynamics).
        MultiArcDynamicsSimulator<> dynamicsSimulator( bodies, propagatorSettings, true );

        std::cout << bodies.getBody( "MGS" )->getEphemeris(  ) << std::endl;
        std::cout << bodies.getBody( "MGS" )->getEphemeris(  )->getCartesianState( arcStartTimes.at(0) ) << std::endl;
        std::cout << bodies.getBody( "MGS" )->getEphemeris(  )->getCartesianState( arcStartTimes.at(0) + 10.) << std::endl;


        std::pair< std::vector< std::shared_ptr< observation_models::ObservationModelSettings > >,
                   std::shared_ptr< observation_models::ObservationCollection< double > > >
                observationCollectionAndModelSettings = simulatePseudoObservations(
                        bodies, bodiesToEstimate, centralBodies, observationTimes );
        std::shared_ptr< observation_models::ObservationCollection< double > > observationCollection =
                observationCollectionAndModelSettings.second;

        std::vector< std::shared_ptr< observation_models::ObservationModelSettings > > observationModelSettingsList =
                observationCollectionAndModelSettings.first;

        std::cout << "getFullParameterValues: " << std::endl;
        Eigen::VectorXd truthParameters = parametersToEstimate->getFullParameterValues< double >( );
        std::cout << "truthParameters: " << truthParameters << std::endl;


        // resetParameterValues
        Eigen::VectorXd perturbation0(3);
        perturbation0 << 0.1, -0.1, 0.0;
        Eigen::VectorXd perturbation1(3);
        perturbation1 << 0.1, -0.1, 0.0;
        Eigen::VectorXd perturbation2(6);
        perturbation2 << 0.1, -0.1, 0.0, -0.6, 0.2, 0.4;

        if (i == 0)
        {
            Eigen::VectorXd paramVector = parametersToEstimate->getFullParameterValues<double>(  );
            paramVector.tail(perturbation0.size()) += perturbation0;
            parametersToEstimate->resetParameterValues( paramVector );
        }
        if (i == 1)
        {
            Eigen::VectorXd paramVector = parametersToEstimate->getFullParameterValues<double>(  );
            paramVector.tail(perturbation1.size()) += perturbation1;
            parametersToEstimate->resetParameterValues( paramVector );
        }
        if (i == 2)
        {

            Eigen::VectorXd paramVector = parametersToEstimate->getFullParameterValues<double>(  );

            for( unsigned int m = 0; m < 3; m++ )
            {
                paramVector[(m*6)] += 10;
                paramVector[(m*6)+1] -= 10;
                paramVector[(m*6)+2] += 10;

                paramVector[(m*6)+3] -= 0.1;
                paramVector[(m*6)+4] += 0.1;
                paramVector[(m*6)+5] += 0.1;
            }
            paramVector.tail(perturbation2.size()) += perturbation2;
            parametersToEstimate->resetParameterValues( paramVector );
        }


        // parametersToEstimate->getFullParameterValues< double >( )
        OrbitDeterminationManager< double > orbitDeterminationManager =
                OrbitDeterminationManager< double >(
                        bodies, parametersToEstimate, observationModelSettingsList, propagatorSettings );

        std::shared_ptr< EstimationInput< double > > estimationInput =
                std::make_shared< EstimationInput< double > >( observationCollection );
        estimationInput->defineEstimationSettings( 0, 1, 0, 1, 1, 1 );

        std::shared_ptr< EstimationOutput<> > estimationOutput = orbitDeterminationManager.estimateParameters( estimationInput );

        residuals.push_back( estimationOutput->residualStandardDeviation_ );
        parameterEstimateList.push_back( estimationOutput->parameterHistory_.at(  estimationOutput->parameterHistory_.size( ) -1 )( 6 ) );

        // Check if parameters are correctly estimated
        Eigen::VectorXd estimatedParametervalues = estimationOutput->parameterEstimate_;

        // Initial states (3 arcs).
        for( unsigned int l = 0; l < 3; l++ )
        {
            BOOST_CHECK_SMALL( std::fabs( truthParameters( l ) - estimationOutput->parameterEstimate_( l ) ), 0.1 );
            BOOST_CHECK_SMALL( std::fabs( truthParameters( l + 6 ) - estimationOutput->parameterEstimate_( l + 6 ) ), 0.1 );
            BOOST_CHECK_SMALL( std::fabs( truthParameters( l + 12 ) - estimationOutput->parameterEstimate_( l + 12 ) ), 0.1 );

            BOOST_CHECK_SMALL( std::fabs( truthParameters( l + 3 ) - estimationOutput->parameterEstimate_( l + 3 ) ), 1.0E-6 );
            BOOST_CHECK_SMALL( std::fabs( truthParameters( l + 9 ) - estimationOutput->parameterEstimate_( l + 9 ) ), 1.0E-6 );
            BOOST_CHECK_SMALL( std::fabs( truthParameters( l + 15 ) - estimationOutput->parameterEstimate_( l + 15 ) ), 1.0E-6 );

        }

        // arcwise aerodynamic coefficients.
        BOOST_CHECK_SMALL( std::fabs( truthParameters( 18 ) - estimationOutput->parameterEstimate_( 18 ) ), 1.0e-4 );
        BOOST_CHECK_SMALL( std::fabs( truthParameters( 19 ) - estimationOutput->parameterEstimate_( 19 ) ), 1.0e-4 );
        BOOST_CHECK_SMALL( std::fabs( truthParameters( 20 ) - estimationOutput->parameterEstimate_( 20 ) ), 1.0e-4 );

        // more arcwise aerodynamic coefficients.
        if (truthParameters.size() == 24)
        {
            BOOST_CHECK_SMALL( std::fabs( truthParameters( 21 ) - estimationOutput->parameterEstimate_( 21 ) ), 1.0e-4 );
            BOOST_CHECK_SMALL( std::fabs( truthParameters( 22 ) - estimationOutput->parameterEstimate_( 22 ) ), 1.0e-4 );
            BOOST_CHECK_SMALL( std::fabs( truthParameters( 23 ) - estimationOutput->parameterEstimate_( 23 ) ), 1.0e-4 );

        }

    }
}



/* Fitting to MGS spice ephemeris, not useful as stand-alone unit test...
BOOST_AUTO_TEST_CASE( test_EstimationArcwiseDragScaling )
{
    std::vector< double > residuals;
    std::vector< double > parameterEstimateList;
    spice_interface::loadStandardSpiceKernels( );

    spice_interface::loadSpiceKernelInTudat( paths::getTudatTestDataPath( ) +
                                             "/dsn_n_way_doppler_observation_model/mgs_map1_ipng_mgs95j.bsp" );

    std::cout << "start "  << std::endl;

    for( unsigned int i = 0; i < 3; i++ )
    {
        std::cout << "i: " << i << std::endl;

        double initialTime = DateTime( 1999, 3, 10, 0, 0, 0.0 ).epoch< double >( );
        double finalTime = initialTime + 3600 * 10;

        std::vector<double> arcStartTimes = {initialTime, initialTime+4*3600, initialTime+7.5*3600};

        std::vector< std::string > bodyNames;
        bodyNames.push_back( "Mars" );
        bodyNames.push_back( "Sun" );

        BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, "Mars" );
        bodySettings.addSettings( "MGS" );
        bodySettings.at( "MGS" )->ephemerisSettings =
            std::make_shared< InterpolatedSpiceEphemerisSettings >( initialTime - 3600.0, finalTime + 3600.0, 30.0, "Mars" );
        bodySettings.at( "MGS" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );

        bodySettings.at( "MGS" )->radiationPressureTargetModelSettings = cannonballRadiationPressureTargetModelSettings(
                10.0, 1.2 );
        // Create bodies needed in simulation
        SystemOfBodies bodies = createSystemOfBodies( bodySettings );
        bodies.at( "Mars" )->setAtmosphereModel(
                createAtmosphereModel( std::make_shared< ExponentialAtmosphereSettings >( aerodynamics::mars ),
                                       "Mars" ) );
        bodies.at( "MGS" )->setConstantBodyMass( 2.0 );

        Eigen::Vector3d forceCoefficients( 1.0, 0.0, 1.0 );

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
        accelerationsOfSpacecraft[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( radiation_pressure ) );

        accelerationMap[ "MGS" ] = accelerationsOfSpacecraft;

        // Set bodies for which initial state is to be estimated and integrated.
        std::vector< std::string > bodiesToEstimate = { "MGS" };
        std::vector< std::string > centralBodies = { "Mars" };

        AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodies, accelerationMap, bodiesToEstimate, centralBodies );

        std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > > additionalParameterNames;

        if( i == 0 )
        {
            additionalParameterNames.push_back( estimatable_parameters::arcwiseDragCoefficient( "MGS",  arcStartTimes) );
        }
        if( i == 1 )
        {
            additionalParameterNames.push_back( estimatable_parameters::arcwiseDragComponentScaling( "MGS", arcStartTimes ) );
        }
        if( i == 2 )
        {
            additionalParameterNames.push_back( estimatable_parameters::arcwiseDragComponentScaling( "MGS", arcStartTimes ) );
            additionalParameterNames.push_back( estimatable_parameters::arcwiseLiftComponentScaling( "MGS", arcStartTimes ) );
        }

        std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;
        dependentVariablesList.push_back( singleAccelerationDependentVariable(
                basic_astrodynamics::aerodynamic, "MGS", "Mars" ) );

        std::vector< double > integrationArcStartTimes = arcStartTimes;
        std::vector< double > integrationArcEndTimes(arcStartTimes.begin() + 1, arcStartTimes.end());
        for (int i = 0; i < integrationArcEndTimes.size( ); ++i) {
            integrationArcEndTimes.at(i) = integrationArcEndTimes.at(i) - 200*(i+1);
        }
        integrationArcEndTimes.push_back( finalTime );

        std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsList;
        for( unsigned int j = 0; j < integrationArcStartTimes.size( ); j++ )
        {
            Eigen::Matrix< double, Eigen::Dynamic, 1 > currentInitialState =
                getInitialStatesOfBodies( bodiesToEstimate, centralBodies, bodies, integrationArcStartTimes.at( j ) );

            propagatorSettingsList.push_back(
                std::make_shared< TranslationalStatePropagatorSettings< double > >(
                    centralBodies,
                    accelerationModelMap,
                    bodiesToEstimate,
                    currentInitialState,
                    integrationArcStartTimes.at( j ),
                    numerical_integrators::rungeKuttaFixedStepSettings( 30.0, numerical_integrators::rungeKuttaFehlberg78 ),
                    propagationTimeTerminationSettings( integrationArcEndTimes.at( j ) ),
                    cowell,
                    dependentVariablesList ) );

        }
        std::cout << "MultiArcPropagatorSettings" << std::endl;
        std::shared_ptr< MultiArcPropagatorSettings< double > > propagatorSettings =
                std::make_shared< MultiArcPropagatorSettings< double > >( propagatorSettingsList );

        std::cout << "************************************* RUNNING TEST *************************" << std::endl;
        std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
                getInitialMultiArcParameterSettings< double >( propagatorSettings, bodies, integrationArcStartTimes );

        parameterNames.insert( parameterNames.end( ), additionalParameterNames.begin( ), additionalParameterNames.end( ) );

        std::cout << "EstimatableParameterSet" << std::endl;
        std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
                createParametersToEstimate< double, double >( parameterNames, bodies, propagatorSettings );
        printEstimatableParameterEntries( parametersToEstimate );

        std::vector< double > observationTimes;
        double dataPointInterval = 120.;
        for( int k = 0; k < integrationArcStartTimes.size( ); k++)
        {
            double currentTime = integrationArcStartTimes.at( k ) + 300.0;
            while( currentTime < integrationArcEndTimes.at( k ) - 300.0 )

            {
                observationTimes.push_back( currentTime );
                currentTime += static_cast< double >( dataPointInterval );
            }
        }

        std::pair< std::vector< std::shared_ptr< observation_models::ObservationModelSettings > >,
                   std::shared_ptr< observation_models::ObservationCollection< double > > >
                observationCollectionAndModelSettings = simulatePseudoObservations(
                        bodies, bodiesToEstimate, centralBodies, observationTimes );
        std::shared_ptr< observation_models::ObservationCollection< double > > observationCollection =
                observationCollectionAndModelSettings.second;

        std::vector< std::shared_ptr< observation_models::ObservationModelSettings > > observationModelSettingsList =
                observationCollectionAndModelSettings.first;

        Eigen::VectorXd truthParameters = parametersToEstimate->getFullParameterValues< double >( );
        std::cout << "truthParameters: " << truthParameters << std::endl;


        // parametersToEstimate->getFullParameterValues< double >( )
        OrbitDeterminationManager< double > orbitDeterminationManager =
                OrbitDeterminationManager< double >(
                        bodies, parametersToEstimate, observationModelSettingsList, propagatorSettings );

        std::shared_ptr< EstimationInput< double > > estimationInput =
                std::make_shared< EstimationInput< double > >( observationCollection );
        estimationInput->defineEstimationSettings( 0, 1, 0, 1, 1, 1 );

        std::shared_ptr< EstimationOutput<> > estimationOutput = orbitDeterminationManager.estimateParameters( estimationInput );

        residuals.push_back( estimationOutput->residualStandardDeviation_ );
        parameterEstimateList.push_back( estimationOutput->parameterHistory_.at(  estimationOutput->parameterHistory_.size( ) -1 )( 6 ) );


        Eigen::VectorXd estimationError = estimationOutput->parameterEstimate_ - truthParameters;
        std::cout << "estimation error: " << ( estimationError ).transpose( ) << std::endl;

        Eigen::MatrixXd correlationMatrix = estimationOutput->getCorrelationMatrix( );
        std::cout << "correlationMatrix: " << correlationMatrix << std::endl;

        // Check if parameters are correctly estimated
        Eigen::VectorXd estimatedParametervalues = estimationOutput->parameterEstimate_;

        // Initial states (3 arcs).
        for( unsigned int l = 0; l < 3; l++ )
        {
            BOOST_CHECK_SMALL( std::fabs( truthParameters( l ) - estimationOutput->parameterEstimate_( l ) ), 0.1 );
            BOOST_CHECK_SMALL( std::fabs( truthParameters( l + 6 ) - estimationOutput->parameterEstimate_( l + 6 ) ), 0.1 );
            BOOST_CHECK_SMALL( std::fabs( truthParameters( l + 12 ) - estimationOutput->parameterEstimate_( l + 12 ) ), 0.1 );

            BOOST_CHECK_SMALL( std::fabs( truthParameters( l + 3 ) - estimationOutput->parameterEstimate_( l + 3 ) ), 1.0E-6 );
            BOOST_CHECK_SMALL( std::fabs( truthParameters( l + 9 ) - estimationOutput->parameterEstimate_( l + 9 ) ), 1.0E-6 );
            BOOST_CHECK_SMALL( std::fabs( truthParameters( l + 15 ) - estimationOutput->parameterEstimate_( l + 15 ) ), 1.0E-6 );

        }

        // arcwise aerodynamic coefficients.
        BOOST_CHECK_SMALL( std::fabs( truthParameters( 18 ) - estimationOutput->parameterEstimate_( 18 ) ), 1.0e-4 );
        BOOST_CHECK_SMALL( std::fabs( truthParameters( 19 ) - estimationOutput->parameterEstimate_( 19 ) ), 1.0e-4 );
        BOOST_CHECK_SMALL( std::fabs( truthParameters( 20 ) - estimationOutput->parameterEstimate_( 20 ) ), 1.0e-4 );

        // more arcwise aerodynamic coefficients.
        if (truthParameters.size() == 24)
        {
            BOOST_CHECK_SMALL( std::fabs( truthParameters( 21 ) - estimationOutput->parameterEstimate_( 21 ) ), 1.0e-4 );
            BOOST_CHECK_SMALL( std::fabs( truthParameters( 22 ) - estimationOutput->parameterEstimate_( 22 ) ), 1.0e-4 );
            BOOST_CHECK_SMALL( std::fabs( truthParameters( 23 ) - estimationOutput->parameterEstimate_( 23 ) ), 1.0e-4 );

        }

    }
}

*/

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat