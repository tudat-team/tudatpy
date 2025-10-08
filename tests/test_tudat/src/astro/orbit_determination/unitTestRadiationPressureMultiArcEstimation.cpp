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

#include <string>
#include <thread>
#include <iomanip>

#include <limits>

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/simulation/estimation.h"
#include "tudat/astro/ephemerides/constantRotationalEphemeris.h"

namespace tudat
{

namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_radiation_pressure_estimation )

// Using declarations.
using namespace tudat;
using namespace tudat::observation_models;
using namespace tudat::orbit_determination;
using namespace tudat::estimatable_parameters;
using namespace tudat::interpolators;
using namespace tudat::numerical_integrators;
using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;
using namespace tudat::orbital_element_conversions;
using namespace tudat::ephemerides;
using namespace tudat::propagators;
using namespace tudat::basic_astrodynamics;
using namespace tudat::coordinate_conversions;
using namespace tudat::ground_stations;
using namespace tudat::observation_models;
using namespace tudat::ephemerides;
using namespace tudat::electromagnetism;
using namespace tudat;

BOOST_AUTO_TEST_CASE( test_RadiationPressureMultiArcVariationalEquations )
{
    // Define initial epoch and analysis period
    std::string initialTimeString = "2012-04-06 01:44:04.686";
    double initialTime = basic_astrodynamics::timeFromIsoString< double >( initialTimeString );
    double arcLength = 600.0;
    double integratorStep = 60.0;
    unsigned int numberOfArcs = 3;

    // Load spice kernels.
    std::string kernelsPath = paths::getSpiceKernelPath( );
    spice_interface::loadStandardSpiceKernels( );

    // Test for using same acceleration model for each arc, or separate models per arc
    for( int createSeparateAccelerations = 0; createSeparateAccelerations < 2; createSeparateAccelerations++ )
    {
        // Test for separate source, and double source
        /*
         *  DIFFERENT ACCELERATION CASES:
         *  CASE 0: Sun isotropic source on Vehicle cannonball target
         *  CASE 1: Moon panelled source on Vehicle cannonball target
         *  CASE 2: Moon panelled source on Vehicle cannonball target AND Sun isotropic source on Vehicle cannonball target
         *  CASE 3: Sun isotropic source on Vehicle panelled target
         *  CASE 4: Moon panelled source on Vehicle panelled target
         *  CASE 5: Moon panelled source on Vehicle panelled target AND Sun isotropic source on Vehicle panelled target
         *  CASE 6: Moon panelled source on Vehicle panelled target AND Sun isotropic source on Vehicle cannonball target
         *  CASE 7: Moon panelled source on Vehicle cannonball target AND Sun isotropic source on Vehicle panelled target
         *
         *  DIFFERENT PARAMETER CASES:
         *  CASE 0: Constant radiation pressure coefficient
         *  CASE 1: Perpendicular scaling factor, Sun source
         *  CASE 2: Parallel scaling factor, Sun source
         *  CASE 3: Perpendicular scaling factor, Moon source
         *  CASE 4: Parallel scaling factor, Moon source
         *  CASE 5: Bus specular reflectivity
         *  CASE 6: Solar panel specular reflectivity
         *  CASE 7: Bus diffuse reflectivity
         *  CASE 8: Solar panel diffuse reflectivity
         *
         *  DIFFERENT PARAMETER CASES (SECONDARY):
         *  CASE 5: Arc-wise radiation pressure coefficient (same arcs as state)
         *  CASE 6: Arc-wise radiation pressure coefficient (shorter arcs than state)
         *  CASE 7: Arc-wise radiation pressure coefficient (longer arcs than state)
         */
        for( unsigned int test = 0; test < 8; test++ )
        {
            // Create settings for default bodies
            std::vector< std::string > bodiesToCreate = { "Earth", "Sun", "Mercury", "Venus", "Mars", "Jupiter", "Saturn", "Moon" };
            std::string globalFrameOrigin = "Moon";
            std::string globalFrameOrientation = "J2000";
            BodyListSettings bodySettings = getDefaultBodySettings( bodiesToCreate, globalFrameOrigin, globalFrameOrientation );

            // Add Moon radiation pressure models
            std::vector< std::shared_ptr< PanelRadiosityModelSettings > > panelRadiosityModels;
            panelRadiosityModels.push_back( angleBasedThermalPanelRadiosityModelSettings( 95.0, 385.0, 0.95, "Sun" ) );
            panelRadiosityModels.push_back(
                    albedoPanelRadiosityModelSettings( SphericalHarmonicsSurfacePropertyDistributionModel::albedo_dlam1, "Sun" ) );
            std::map< std::string, std::vector< std::string > > originalSourceToSourceOccultingBodies;
            originalSourceToSourceOccultingBodies[ "Sun" ].push_back( "Earth" );
            bodySettings.at( "Moon" )->radiationSourceModelSettings = extendedRadiationSourceModelSettingsWithOccultationMap(
                    panelRadiosityModels, { 4, 8 }, originalSourceToSourceOccultingBodies );

            // Add spacecraft settings
            std::string spacecraftName = "GRAIL-A";
            std::string spacecraftCentralBody = "Moon";
            bodySettings.addSettings( spacecraftName );
            bodySettings.at( spacecraftName )->constantMass = 150.0;

            // Create vehicle radiation pressure target settings
            std::map< std::string, std::vector< std::string > > sourceToTargetOccultingBodies;
            sourceToTargetOccultingBodies[ "Sun" ].push_back( "Moon" );
            double referenceAreaRadiation = 5.0E5;
            double radiationPressureCoefficient = 1.5;
            std::vector< std::shared_ptr< RadiationPressureTargetModelSettings > > radiationPressureTargetModelSettings;
            radiationPressureTargetModelSettings.push_back(
                    paneledRadiationPressureTargetModelSettingsWithOccultationMap( sourceToTargetOccultingBodies ) );
            radiationPressureTargetModelSettings.push_back( cannonballRadiationPressureTargetModelSettingsWithOccultationMap(
                    referenceAreaRadiation, radiationPressureCoefficient, sourceToTargetOccultingBodies ) );
            bodySettings.at( spacecraftName )->radiationPressureTargetModelSettings =
                    std::make_shared< MultiRadiationPressureTargetModelSettings >( radiationPressureTargetModelSettings );

            // Create spacecraft panel settings
            bodySettings.get( spacecraftName )->bodyExteriorPanelSettings_ = bodyWingPanelledGeometry( 2.0 * 100.0 * std::sqrt( 10.0 ),
                                                                                                       1.0 * 100.0 * std::sqrt( 10.0 ),
                                                                                                       4.0 * 100.0 * std::sqrt( 10.0 ),
                                                                                                       1.0E5,
                                                                                                       0.1,
                                                                                                       0.2,
                                                                                                       0.2,
                                                                                                       0.3,
                                                                                                       false,
                                                                                                       false );

            // Create dummy spacecraft rotation settings
            bodySettings.get( spacecraftName )->rotationModelSettings =
                    constantRotationModelSettings( "J2000", spacecraftName + "_Fixed", Eigen::Matrix3d::Identity( ) );

            // Create bodies
            SystemOfBodies bodies = createSystemOfBodies< long double, Time >( bodySettings );

            // Set accelerations between bodies that are to be taken into account.
            SelectedAccelerationMap accelerationMap;
            std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfSpacecraft;
            std::vector< int > parameterIndicesToUse;
            if( test == 0 || test == 2 || test == 6 )
            {
                accelerationsOfSpacecraft[ "Sun" ].push_back( std::make_shared< RadiationPressureAccelerationSettings >( cannonball_target ) );
                parameterIndicesToUse.push_back( 0 );
                parameterIndicesToUse.push_back( 1 );
            }
            if( test == 1 || test == 2 || test == 7 )
            {
                accelerationsOfSpacecraft[ "Moon" ].push_back( std::make_shared< RadiationPressureAccelerationSettings >( cannonball_target ) );
                parameterIndicesToUse.push_back( 0 );
                parameterIndicesToUse.push_back( 3 );
            }
            if( test == 3 || test == 5 || test == 7 )
            {
                accelerationsOfSpacecraft[ "Sun" ].push_back( std::make_shared< RadiationPressureAccelerationSettings >( paneled_target ) );
                parameterIndicesToUse.push_back( 1 );
                parameterIndicesToUse.push_back( 2 );
                parameterIndicesToUse.push_back( 5 );
                parameterIndicesToUse.push_back( 6 );
                parameterIndicesToUse.push_back( 7 );
                parameterIndicesToUse.push_back( 8 );
            }
            if( test == 4 || test == 5 || test == 6 )
            {
                accelerationsOfSpacecraft[ "Moon" ].push_back( std::make_shared< RadiationPressureAccelerationSettings >( paneled_target ) );
                parameterIndicesToUse.push_back( 3 );
                parameterIndicesToUse.push_back( 4 );
                parameterIndicesToUse.push_back( 5 );
                parameterIndicesToUse.push_back( 6 );
                parameterIndicesToUse.push_back( 7 );
                parameterIndicesToUse.push_back( 8 );
            }
            utilities::removeDuplicates( parameterIndicesToUse );
            accelerationMap[ "GRAIL-A" ] = accelerationsOfSpacecraft;

            // Set bodies for which initial state is to be estimated and integrated.
            std::vector< std::string > bodiesToEstimate = { "GRAIL-A" };
            std::vector< std::string > centralBodies = { "Moon" };
            AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodies, accelerationMap, bodiesToEstimate, centralBodies );

            // Reasonable initial states from Webgeocalc with semi-arbitrary initial times
            std::vector< Eigen::VectorXd > initialStates;
            initialStates.push_back(
                    ( Eigen::VectorXd( 6 ) << -618712.81427, 557495.52023, 1597758.84330, 530.76325, -1397.57332, 687.35537 ).finished( ) );
            initialStates.push_back(
                    ( Eigen::VectorXd( 6 ) << -270294.33519, 1160199.58017, -1335000.99409, -738.74766, 1046.55086, 1050.72524 ).finished( ) );
            initialStates.push_back(
                    ( Eigen::VectorXd( 6 ) << 836709.02648, -1521295.66730, -444616.38949, 94.55259, 515.43686, -1569.18467 ).finished( ) );
            initialStates.push_back(
                    ( Eigen::VectorXd( 6 ) << -437317.69399, 129623.10246, 1736785.90701, 661.17560, -1486.15366, 281.93859 ).finished( ) );

            std::vector< std::shared_ptr< SingleArcPropagatorSettings< long double, Time > > > propagatorSettingsList;
            for( unsigned int arcIndex = 0; arcIndex < numberOfArcs; arcIndex++ )
            {
                std::shared_ptr< TranslationalStatePropagatorSettings< long double, Time > > propagatorSettingsArc =
                        std::make_shared< TranslationalStatePropagatorSettings< long double, Time > >(
                                centralBodies,
                                ( createSeparateAccelerations == 0 ) ? accelerationModelMap : createAccelerationModelsMap( bodies, accelerationMap, bodiesToEstimate, centralBodies ),
                                bodiesToEstimate,
                                initialStates.at( arcIndex % initialStates.size( ) ).template cast< long double >( ),
                                Time( initialTime + arcIndex * arcLength ),
                                numerical_integrators::rungeKuttaFixedStepSettings< Time >( integratorStep,
                                                                                            numerical_integrators::rungeKuttaFehlberg78 ),
                                std::make_shared< PropagationTimeTerminationSettings >( initialTime + ( arcIndex + 1 ) * arcLength ) );
                propagatorSettingsList.push_back( propagatorSettingsArc );
            }

            // Define multi-arc propagator settings
            std::shared_ptr< MultiArcPropagatorSettings< long double, Time > > propagatorSettings =
                    std::make_shared< MultiArcPropagatorSettings< long double, Time > >( propagatorSettingsList );

            // Define arc durations over shorter and longer arcs thatn states
            std::vector< double > initialStateArcStartTimes = propagatorSettings->getArcStartTimes( );
            std::vector< double > shortArcParameterStartTimes;
            for( unsigned int i = 0; i < initialStateArcStartTimes.size( ); i = i + 2 )
            {
                shortArcParameterStartTimes.push_back( initialStateArcStartTimes.at( i ) );
                shortArcParameterStartTimes.push_back( initialStateArcStartTimes.at( i ) + 2.0 / 3.0 * arcLength );
                shortArcParameterStartTimes.push_back( initialStateArcStartTimes.at( i ) + 4.0 / 3.0 * arcLength );
            }

            std::vector< double > longArcParameterStartTimes;
            for( unsigned int i = 0; i < initialStateArcStartTimes.size( ); i = i + 3 )
            {
                longArcParameterStartTimes.push_back( initialStateArcStartTimes.at( i ) );
                longArcParameterStartTimes.push_back( initialStateArcStartTimes.at( i ) + 1.5 * arcLength );
            }

            // Define full list of estimated parameter cases
            std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > > fullParameterNames;

            fullParameterNames.push_back( estimatable_parameters::radiationPressureCoefficient( "GRAIL-A" ) );
            fullParameterNames.push_back( estimatable_parameters::radiationPressureTargetDirectionScaling( "GRAIL-A", "Sun" ) );
            fullParameterNames.push_back( estimatable_parameters::radiationPressureTargetPerpendicularDirectionScaling( "GRAIL-A", "Sun" ) );
            fullParameterNames.push_back( estimatable_parameters::radiationPressureTargetDirectionScaling( "GRAIL-A", "Moon" ) );
            fullParameterNames.push_back( estimatable_parameters::radiationPressureTargetPerpendicularDirectionScaling( "GRAIL-A", "Moon" ) );
            fullParameterNames.push_back( estimatable_parameters::specularReflectivity( "GRAIL-A", "Bus" ) );
            fullParameterNames.push_back( estimatable_parameters::specularReflectivity( "GRAIL-A", "SolarPanel" ) );
            fullParameterNames.push_back( estimatable_parameters::diffuseReflectivity( "GRAIL-A", "Bus" ) );
            fullParameterNames.push_back( estimatable_parameters::diffuseReflectivity( "GRAIL-A", "SolarPanel" ) );

            // Create parameters for initial states
            std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
                    getInitialStateParameterSettings< long double, Time >( propagatorSettings, bodies );

            // Add specific parameters for current run
            for( unsigned int i = 0; i < parameterIndicesToUse.size( ); i++ )
            {
                parameterNames.push_back( fullParameterNames.at( parameterIndicesToUse.at( i ) ) );
            }

            std::shared_ptr< estimatable_parameters::EstimatableParameterSet< long double > > parametersToEstimate =
                    createParametersToEstimate< long double, Time >( parameterNames, bodies, propagatorSettings );

            // Propagate variational equations
            MultiArcVariationalEquationsSolver< long double, Time > variationalEquationsSolver(
                    bodies, propagatorSettings, parametersToEstimate, true );
            std::shared_ptr< CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionMatrixInterface =
                    variationalEquationsSolver.getStateTransitionMatrixInterface( );
            auto variationalPropagationResults = variationalEquationsSolver.getMultiArcVariationalPropagationResults( );
            std::vector< std::shared_ptr< SingleArcVariationalSimulationResults< long double, Time > > > singleArcResults = variationalPropagationResults->getSingleArcResults( );

            std::vector< std::map< double, Eigen::MatrixXd > > stateTransitionResults;
            std::vector< std::map< double, Eigen::MatrixXd > > sensitivityResults;
            for( unsigned int arc = 0; arc < numberOfArcs; arc++ )
            {
                stateTransitionResults.push_back( singleArcResults.at( arc )->getStateTransitionSolution( ) );
                sensitivityResults.push_back( singleArcResults.at( arc )->getSensitivitySolution( ) );
            }

            // Iterate over all parameters
            printEstimatableParameterEntries( parametersToEstimate );
            auto nominalParameters = parametersToEstimate->getFullParameterValues< long double >( );
            for( unsigned int parameterIndex = 0; parameterIndex < static_cast< unsigned int >( parametersToEstimate->getParameterSetSize( ) );
                 parameterIndex++ )
            {
                // Parameter perturbations and tolerances determined empirically to be acceptable
                // Non-linearity of effect, weakness of Moon radiation, and drift away from Moon, make test tolerance difficult to get more tight
                double parameterPerturbation = TUDAT_NAN;
                if( parameterIndex >= 6 * numberOfArcs )
                {
                    parameterPerturbation = 0.001;
                }
                else
                {
                    if( ( parameterIndex % 6 ) < 3 )
                    {
                        parameterPerturbation = 1.0;
                    }
                    else
                    {
                        parameterPerturbation = 1.0E-3;
                    }
                }
                int scalingIndex = 2;
                if( test == 0 || test == 3 )
                {
                    scalingIndex = 4;
                }
                parameterPerturbation *= std::pow( 10.0, scalingIndex );

                // Set test tolerance
                double stateTolerance = 1.0E-2;
                double parameterTolerance = 1.0E-4;
                if( test == 2 || test == 7 )
                {
                    stateTolerance = 1.0E-3;
                }
                else if( test == 0 || test == 3 )
                {
                    stateTolerance = 1.0E-6;
                    parameterTolerance = 1.0E-10;
                }

                // Propagate with up-perturbed parameters
                auto perturbedParameters = nominalParameters;
                perturbedParameters( parameterIndex ) += parameterPerturbation;
                parametersToEstimate->resetParameterValues( perturbedParameters );
                propagatorSettings->resetInitialStates( perturbedParameters.segment( 0, 6 * numberOfArcs ) );
                MultiArcDynamicsSimulator< long double, Time > upperturbedDynamics( bodies, propagatorSettings );
                std::vector< std::map< Time, Eigen::Matrix< long double, Eigen::Dynamic, 1 > > > upperturbedResults =
                        upperturbedDynamics.getEquationsOfMotionNumericalSolution( );

                // Propagate with down-perturbed parameters
                perturbedParameters = nominalParameters;
                perturbedParameters( parameterIndex ) -= parameterPerturbation;
                parametersToEstimate->resetParameterValues( perturbedParameters );
                propagatorSettings->resetInitialStates( perturbedParameters.segment( 0, 6 * numberOfArcs ) );
                MultiArcDynamicsSimulator< long double, Time > downperturbedDynamics( bodies, propagatorSettings );
                std::vector< std::map< Time, Eigen::Matrix< long double, Eigen::Dynamic, 1 > > > downperturbedResults =
                        downperturbedDynamics.getEquationsOfMotionNumericalSolution( );

                // Reset
                parametersToEstimate->resetParameterValues( nominalParameters );
                propagatorSettings->resetInitialStates( nominalParameters.segment( 0, 6 * numberOfArcs ) );

                // Iterate over all arcs
                //            std::cout<<std::endl;
                for( unsigned int arc = 0; arc < numberOfArcs; arc++ )
                {
                    // Retrieve current arc variational equation results
                    std::map< double, Eigen::MatrixXd > currentArcStateTransitionResults = stateTransitionResults.at( arc );
                    std::map< double, Eigen::MatrixXd > currentArcSensitivityResults = sensitivityResults.at( arc );

                    // Retrieve current arc up- and down-perturbed state results
                    std::map< Time, Eigen::Matrix< long double, Eigen::Dynamic, 1 > > currentArcUpperturbedState = upperturbedResults.at( arc );
                    std::map< Time, Eigen::Matrix< long double, Eigen::Dynamic, 1 > > currentArcDownperturbedState =
                            downperturbedResults.at( arc );

                    // Create iterator over state transition or sensitivity matrix AND RETRIEVE FINAL VALUE
                    bool useStateTransition = parameterIndex < ( 6 * numberOfArcs );
                    auto mapToIterate = useStateTransition ? currentArcStateTransitionResults : currentArcSensitivityResults;
                    auto matrixIterator = mapToIterate.rbegin( );

                    // Define column to test
                    int matrixColumn = 0;
                    if( parameterIndex > 6 * numberOfArcs )
                    {
                        matrixColumn = parameterIndex - 6 * numberOfArcs;
                    }
                    else
                    {
                        matrixColumn = parameterIndex % 6;
                    }

                    // Compare values
                    //                std::cout<<std::setprecision( 12 );
                    //                std::cout<<"INITIAL UP  :   "<<currentArcUpperturbedState.begin( )->second<<std::endl;
                    //                std::cout<<"INITIAL DOWN:   "<<currentArcUpperturbedState.begin( )->second<<std::endl;
                    //
                    //                std::cout<<"UP:   "<<currentArcUpperturbedState.at( matrixIterator->first )<<std::endl;
                    //                std::cout<<"DOWN: "<<currentArcDownperturbedState.at( matrixIterator->first )<<std::endl;
                    //                std::cout<<std::setprecision( 6 );
                    Eigen::VectorXd analyticalValue = matrixIterator->second.block( 0, matrixColumn, 6, 1 );
                    Eigen::VectorXd numericalValue = ( ( currentArcUpperturbedState.at( matrixIterator->first ) -
                                                         currentArcDownperturbedState.at( matrixIterator->first ) )
                                                               .transpose( ) /
                                                       ( 2.0 * parameterPerturbation ) )
                                                             .cast< double >( );
                    // Remove initial unity value for testing resolution
                    if( ( ( arc == parameterIndex / 6 ) && useStateTransition ) || !useStateTransition )
                    {
                        if( parameterIndex < 6 * numberOfArcs )
                        {
                            analyticalValue( matrixColumn ) -= 1.0;
                            numericalValue( matrixColumn ) -= 1.0;
                        }
                        Eigen::VectorXd errorLevels = ( ( numericalValue - analyticalValue ).cwiseQuotient( analyticalValue ) );

                        if( parameterIndex < 6 * numberOfArcs )
                        {
                            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( numericalValue, analyticalValue, stateTolerance );
                        }
                        else
                        {
                            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( numericalValue, analyticalValue, parameterTolerance );
                        }
                        std::cout << " CASE: " << test << " PARAMETER: " << parameterIndex << " ARC: " << arc
                                  << " MATRIX COLUMN: " << matrixColumn <<" CREATE "<<createSeparateAccelerations
                                  << " ERROR:                                     " << errorLevels.array( ).abs( ).maxCoeff( )
                                  << std::endl;
                        //                    {
                        //                        std::cout << "ANALYTICAL: " << analyticalValue.transpose( ) << std::endl;
                        //                        std::cout << "NUMERICAL : " << numericalValue.transpose( ) << std::endl;
                        //                        std::cout << "RATIO     : "
                        //                                  << ( ( numericalValue - analyticalValue ).cwiseQuotient( analyticalValue ) ).transpose( ) << std::endl
                        //                                  << std::endl;
                        //                        std::cout<<"Current matrix "<<matrixIterator->second<<std::endl<<std::endl;
                        //
                        //                    }
                    }
                }
            }
        }
    }
}


BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
