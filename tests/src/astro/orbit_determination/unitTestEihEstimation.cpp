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

#include <boost/test/unit_test.hpp>

#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/simulation/simulation.h"
#include "tudat/simulation/estimation.h"

namespace tudat
{

namespace unit_tests
{

using namespace tudat;
using namespace tudat::simulation_setup;
using namespace tudat::propagators;
using namespace tudat::numerical_integrators;
using namespace tudat::orbital_element_conversions;
using namespace tudat::basic_mathematics;
using namespace tudat::basic_astrodynamics;
using namespace tudat::unit_conversions;
using namespace tudat::estimatable_parameters;
using namespace tudat::orbit_determination;
using namespace tudat::observation_models;


BOOST_AUTO_TEST_SUITE( test_eih_estimation )

BOOST_AUTO_TEST_CASE( testEihEstimation )
{
    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation end epoch.
    const double simulationStartEpoch = 0.0 * tudat::physical_constants::JULIAN_YEAR;
    const double simulationEndEpoch = 25.0 * 365.0 * tudat::physical_constants::JULIAN_DAY;
    const double stepSize = 2.0 * 86400.0;


    std::vector< Eigen::Vector6d > rmsDifferences;
    std::vector< Eigen::Vector6d > maximumDifferences;

    std::vector<std::string> bodiesToPropagate = { "Venus" };
    std::vector<std::string> centralBodies = std::vector<std::string>( bodiesToPropagate.size( ), "SSB" );

    std::vector< Eigen::VectorXd > parameterCorrections;
    for( unsigned int test = 0; test < 2; test++ )
    {
        // Create body objects.
        std::vector<std::string> bodiesToCreate =
        { "Sun", "Mercury", "Venus", "Earth", "Moon", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune" };
        BodyListSettings bodySettings =
            getDefaultBodySettings( bodiesToCreate );
        for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
        {
            if( bodiesToCreate.at( i ) != "Sun" && bodiesToCreate.at( i ) != "Mercury" && bodiesToCreate.at( i ) != "Venus" && bodiesToCreate.at( i ) != "Earth" && bodiesToCreate.at( i ) != "Moon"  )
            {
                bodySettings.at( bodiesToCreate.at( i ) )->ephemerisSettings = directSpiceEphemerisSettings(
                    "SSB", "ECLIPJ2000", bodiesToCreate.at( i ) + "_Barycenter" );
            }
            bodySettings.at( bodiesToCreate.at( i ) )->gravityFieldSettings = centralGravityFromSpiceSettings( );
        }

        for( unsigned int i = 0; i < bodiesToPropagate.size( ); i++ )
        {
            bodySettings.at( bodiesToPropagate.at( i ) )->ephemerisSettings = tabulatedEphemerisSettings(
                bodySettings.at( bodiesToPropagate.at( i ) )->ephemerisSettings, simulationStartEpoch - 10.0 * stepSize, simulationEndEpoch + 10.0 * stepSize, stepSize );
        }

        std::vector<std::string> minorBodies = { "Ceres", "Vesta", "Pallas", "2000010", "2000704", "2000052" };
        for( unsigned int i = 0; i < minorBodies.size( ); i++ )
        {
            bodySettings.addSettings( minorBodies.at( i ) );
            bodySettings.at( minorBodies.at( i ) )->ephemerisSettings = directSpiceEphemerisSettings( "SSB", "ECLIPJ2000", false, false, false );
            bodySettings.at( minorBodies.at( i ) )->gravityFieldSettings = centralGravityFromSpiceSettings( );
        }


        // Create Earth object
        SystemOfBodies bodies = createSystemOfBodies( bodySettings );

        auto accelerationType = ( test == 0 ? einstein_infeld_hoffmann_acceleration : point_mass_gravity );
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Define propagator settings variables.
        SelectedAccelerationMap accelerationMap;


        // Define propagation settings.
        for ( unsigned int i = 0; i < bodiesToPropagate.size( ); i++ )
        {
            std::map<std::string, std::vector<std::shared_ptr<AccelerationSettings> > > accelerationsOfCurrentBody;
            for ( unsigned int j = 0; j < bodiesToCreate.size( ); j++ )
            {
                if ( bodiesToPropagate.at( i ) != bodiesToCreate.at( j ))
                {
                    accelerationsOfCurrentBody[ bodiesToCreate.at( j ) ].push_back(
                        std::make_shared<AccelerationSettings>(
                            accelerationType ));
                }
            }

            for ( unsigned int j = 0; j < minorBodies.size( ); j++ )
            {
                accelerationsOfCurrentBody[ minorBodies.at( j )].push_back(
                    std::make_shared<AccelerationSettings>(
                        point_mass_gravity ) );
            }

            accelerationMap[ bodiesToPropagate.at( i ) ] = accelerationsOfCurrentBody;
        }

        // Create acceleration models and propagation settings.
        basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
            bodies, accelerationMap, bodiesToPropagate, centralBodies );



        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        Eigen::VectorXd systemInitialState = getInitialStatesOfBodies(
            bodiesToPropagate, centralBodies, bodies, simulationStartEpoch );

        std::shared_ptr<IntegratorSettings<> >
            integratorSettings = std::make_shared<RungeKuttaFixedStepSizeSettings<> >( stepSize, CoefficientSets::rungeKutta87DormandPrince );

        std::shared_ptr<TranslationalStatePropagatorSettings<double> > propagatorSettings =
            std::make_shared<TranslationalStatePropagatorSettings<double> >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, simulationStartEpoch,
                  integratorSettings,
                  std::make_shared<PropagationTimeTerminationSettings>( simulationEndEpoch ),
                  cowell );
//        propagatorSettings->getOutputSettings( )->setResultsSaveFrequencyInSteps( 20 );
//        propagatorSettings->getPrintSettings( )->setResultsPrintFrequencyInSeconds( )

        std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
            getInitialStateParameterSettings< double >( propagatorSettings, bodies );
        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Sun", gravitational_parameter ) );
        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Mercury", gravitational_parameter ) );
//        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Venus", gravitational_parameter ) );
        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Earth", gravitational_parameter ) );
        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Mars", gravitational_parameter ) );
        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Jupiter", gravitational_parameter ) );
        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Saturn", gravitational_parameter ) );
        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Uranus", gravitational_parameter ) );
        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Neptune", gravitational_parameter ) );


        
        // Create initial state estimation objects
        std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double >( parameterNames, bodies );
        
        // Add current body to list of observed bodies
        std::vector< std::shared_ptr< ObservationModelSettings > >  observationSettingsList;
        std::vector< LinkEnds > linkEndsList;

        for( unsigned int i = 0; i < bodiesToPropagate.size( ); i++ )
        {
            LinkEnds observationLinkEnds;
            observationLinkEnds[ observed_body ] = std::pair<std::string, std::string>( std::make_pair( bodiesToPropagate.at( i ), "" ));
            linkEndsList.push_back( observationLinkEnds );
            observationSettingsList.push_back( std::make_shared<ObservationModelSettings>(
                position_observable, observationLinkEnds ));
        }


        std::vector< double > baseTimeList;
        double currentTime = simulationStartEpoch + 10.0 * stepSize;
        while( currentTime < simulationEndEpoch - 10.0 * 86400 )
        {
            baseTimeList.push_back( currentTime );
            currentTime += stepSize;
        }
        std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > measurementSimulationInput;
        for( unsigned int i = 0; i < bodiesToPropagate.size( ); i++ )
        {
            measurementSimulationInput.push_back(
                std::make_shared<TabulatedObservationSimulationSettings<> >(
                    position_observable, linkEndsList.at( i ), baseTimeList, observed_body ));
        }

        // Simulate ideal observations
        std::shared_ptr< ObservationCollection< > > observationsAndTimes = simulateObservations< double, double >(
            measurementSimulationInput, createObservationSimulators(
                observationSettingsList, bodies ), bodies );

        // Create orbit determination object.
        OrbitDeterminationManager< double, double > orbitDeterminationManager =
            OrbitDeterminationManager< double, double >(
                bodies, parametersToEstimate, observationSettingsList,
                propagatorSettings );
        
        // Define estimation input
        std::shared_ptr< EstimationInput< double, double > > estimationInput =
            std::make_shared< EstimationInput< double, double > >(
                observationsAndTimes );
//        estimationInput->defineEstimationSettings( true, true, true, true, true );
        estimationInput->setConvergenceChecker(
            std::make_shared< EstimationConvergenceChecker >( 2 ) );

    

        // Fit nominal dynamics to pertrubed dynamical model
        Eigen::VectorXd trueParameters = parametersToEstimate->getFullParameterValues< double >( );
        std::shared_ptr< EstimationOutput< double > > estimationOutput = orbitDeterminationManager.estimateParameters(
            estimationInput );

        if( test == 0 )
        {
            BOOST_CHECK_EQUAL( estimationOutput->residualStandardDeviation_ < 15.0, true );
        }
        else
        {
            BOOST_CHECK_EQUAL( estimationOutput->residualStandardDeviation_ > 150.0, true );
        }

        Eigen::VectorXd estimatedParameters = parametersToEstimate->getFullParameterValues< double >( );
        Eigen::VectorXd parameterCorrection = estimatedParameters - trueParameters;
        std::cout<<std::endl<<"PARAMETER CORRECTION  ************"<<std::endl<<
        parameterCorrection.transpose( )<<std::endl<<std::endl;

        parameterCorrections.push_back( parameterCorrection );
//        Eigen::VectorXd bestResiduals = estimationOutput->residuals_;
//        input_output::writeMatrixToFile( bestResiduals, "residuals_eih_" + std::to_string( test ) + ".dat", 16 );
    }

    BOOST_CHECK_SMALL( parameterCorrections.at( 0 ).segment( 0, 3 ).norm( ) /
                       parameterCorrections.at( 1 ).segment( 0, 3 ).norm( ), 0.025 );

    BOOST_CHECK_SMALL( parameterCorrections.at( 0 ).segment( 3, 3 ).norm( ) /
                       parameterCorrections.at( 1 ).segment( 3, 3 ).norm( ), 0.025 );

    BOOST_CHECK_SMALL( parameterCorrections.at( 0 )( 6 ) /
                       parameterCorrections.at( 1 )( 6 ), 0.00025 );
}

BOOST_AUTO_TEST_SUITE_END( )

}

}


