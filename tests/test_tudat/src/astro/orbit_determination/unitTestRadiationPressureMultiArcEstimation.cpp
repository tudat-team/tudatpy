/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

//#define BOOST_TEST_DYN_LINK
//#define BOOST_TEST_MAIN

#include <string>
#include <thread>

#include <limits>

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/simulation/estimation.h"
#include "tudat/astro/ephemerides/constantRotationalEphemeris.h"

//namespace tudat
//{
//
//namespace unit_tests
//{

//BOOST_AUTO_TEST_SUITE( test_radiation_pressure_estimation )

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

int main( )
{
    // Define initial epoch and analysis period
    std::string initialTimeString = "2012-04-06 01:44:04.686";
    double initialTime = basic_astrodynamics::timeFromIsoString< double >( initialTimeString );
    double arcLength = 600.0;
    double integratorStep = 60.0;
    int numberOfArcs = 3;

    // Load spice kernels.
    std::string kernelsPath = paths::getSpiceKernelPath( );
    spice_interface::loadStandardSpiceKernels( );

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
    bodySettings.at( "Moon" )->radiationSourceModelSettings =
            extendedRadiationSourceModelSettingsWithOccultationMap( panelRadiosityModels, { 4, 8 }, originalSourceToSourceOccultingBodies );

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
    radiationPressureTargetModelSettings.push_back( 
                    cannonballRadiationPressureTargetModelSettingsWithOccultationMap(
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
        std::vector< int > parameterIndicesToUse;

        // Create bodies
        SystemOfBodies bodies = createSystemOfBodies< long double, Time >( bodySettings );

        // Set accelerations between bodies that are to be taken into account.
        SelectedAccelerationMap accelerationMap;
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfSpacecraft;
        if( test == 0 || test == 2 || test == 6 )
        {
            accelerationsOfSpacecraft[ "Sun" ].push_back( std::make_shared< RadiationPressureAccelerationSettings >( cannonball_target ) );
            parameterIndicesToUse.push_back( 0 );
            parameterIndicesToUse.push_back( 2 );
        }
        if( test == 1 || test == 2 || test == 7 )
        {
            accelerationsOfSpacecraft[ "Moon" ].push_back( std::make_shared< RadiationPressureAccelerationSettings >( cannonball_target ) );
            parameterIndicesToUse.push_back( 0 );
            parameterIndicesToUse.push_back( 4 );
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
        std::vector<Eigen::VectorXd> initialStates;
        initialStates.push_back((Eigen::VectorXd(6) << -618712.81427, 557495.52023, 1597758.84330, 530.76325, -1397.57332, 687.35537).finished());
        initialStates.push_back((Eigen::VectorXd(6) << -270294.33519, 1160199.58017, -1335000.99409, -738.74766, 1046.55086, 1050.72524).finished());
        initialStates.push_back((Eigen::VectorXd(6) << 836709.02648, -1521295.66730, -444616.38949, 94.55259, 515.43686, -1569.18467).finished());
        initialStates.push_back((Eigen::VectorXd(6) << -437317.69399, 129623.10246, 1736785.90701, 661.17560, -1486.15366, 281.93859).finished());


        std::vector< std::shared_ptr< SingleArcPropagatorSettings< long double, Time > > > propagatorSettingsList;
        for( unsigned int arcIndex = 0; arcIndex < numberOfArcs; arcIndex++ )
        {
            std::shared_ptr< TranslationalStatePropagatorSettings< long double, Time > > propagatorSettingsArc =
                    std::make_shared< TranslationalStatePropagatorSettings< long double, Time > >(
                            centralBodies,
                            accelerationModelMap,
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
        for( int i = 0; i < initialStateArcStartTimes.size( ); i = i + 2 )
        {
            shortArcParameterStartTimes.push_back( initialStateArcStartTimes.at( i ) );
            shortArcParameterStartTimes.push_back( initialStateArcStartTimes.at( i ) + 2.0 / 3.0 * arcLength );
            shortArcParameterStartTimes.push_back( initialStateArcStartTimes.at( i ) + 4.0 / 3.0 * arcLength );
        }

        std::vector< double > longArcParameterStartTimes;
        for( int i = 0; i < initialStateArcStartTimes.size( ); i = i + 3 )
        {
            longArcParameterStartTimes.push_back( initialStateArcStartTimes.at( i ) );
            longArcParameterStartTimes.push_back( initialStateArcStartTimes.at( i ) + 1.5 * arcLength );
        }

        // Define full list of estimated parameter cases
        std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > > fullParameterNames;
        fullParameterNames.push_back( estimatable_parameters::radiationPressureCoefficient( "GRAIL-A" ) );
        fullParameterNames.push_back( estimatable_parameters::radiationPressureTargetPerpendicularDirectionScaling( "GRAIL-A", "Sun" ) );
        fullParameterNames.push_back( estimatable_parameters::radiationPressureTargetDirectionScaling( "GRAIL-A", "Sun" ) );
        fullParameterNames.push_back( estimatable_parameters::radiationPressureTargetPerpendicularDirectionScaling( "GRAIL-A", "Moon" ) );
        fullParameterNames.push_back( estimatable_parameters::radiationPressureTargetDirectionScaling( "GRAIL-A", "Moon" ) );
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
        std::vector< std::vector< std::map< double, Eigen::MatrixXd > > > directVariationalResults =
                variationalEquationsSolver.getNumericalVariationalEquationsSolution( );
        std::shared_ptr< CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionMatrixInterface =
                variationalEquationsSolver.getStateTransitionMatrixInterface( );
//        auto variationalPropagationResults = variationalEquationsSolver.getMultiArcVariationalPropagationResults( );
//        std::vector< std::shared_ptr< SingleArcVariationalSimulationResults< long double, Time > > > singleArcResults = variationalPropagationResults->getSingleArcResults( );

        std::vector< std::map< double, Eigen::MatrixXd > > stateTransitionResults;
        std::vector< std::map< double, Eigen::MatrixXd > > sensitivityResults;
        for( int arc = 0; arc < numberOfArcs; arc++ )
        {
            stateTransitionResults.push_back( directVariationalResults.at( arc ).at( 0 ) );
            sensitivityResults.push_back( directVariationalResults.at( arc ).at( 1 ) );
//            std::cout<<"Arc sizes direct "<<directVariationalResults.at( arc ).at( 0 ).size( )<<" "<<directVariationalResults.at( arc ).at( 1 ).size( )<<std::endl;
//            std::cout<<"Arc sizes object "<<stateTransitionResults.at( arc ).size( )<<" "<<sensitivityResults.at( arc ).size( )<<std::endl;
        }


        // Iterate over all parameters
        auto nominalParameters = parametersToEstimate->getFullParameterValues< long double >( );
        for( unsigned int parameterIndex = 6 * numberOfArcs; parameterIndex < static_cast< unsigned int >( parametersToEstimate->getParameterSetSize( ) );
             parameterIndex++ )
        {
            std::cout <<"CASE: "<< test << " PARAMETER: " << parameterIndex << std::endl;

            // Parameter perturbations and tolerances determined empirically to be acceptable
            int scalingIndex = 4;
            double toleranceStates = 1E-3;
            double toleranceParameter = 1E-12;
            if( parameterIndex > 6 * numberOfArcs )
            {
                scalingIndex = 1;
                toleranceParameter = 1.0E-7;
            }

            // Perturb parameters
            auto perturbedParameters = nominalParameters;
            double parameterPerturbation = TUDAT_NAN;
            double tolerance = TUDAT_NAN;
            
            if( parameterIndex > 6 * numberOfArcs )
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
            parameterPerturbation *= std::pow( 10.0, scalingIndex );

            // Propagate with up-perturbed parameters
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
            propagatorSettings->resetInitialStates( perturbedParameters.segment( 0, 6 * numberOfArcs  ) );
            MultiArcDynamicsSimulator< long double, Time > downperturbedDynamics( bodies, propagatorSettings );
            std::vector< std::map< Time, Eigen::Matrix< long double, Eigen::Dynamic, 1 > > > downperturbedResults =
                    downperturbedDynamics.getEquationsOfMotionNumericalSolution( );

            // Reset
            parametersToEstimate->resetParameterValues( nominalParameters );
            propagatorSettings->resetInitialStates( nominalParameters.segment( 0, 6 * numberOfArcs  ) );


            // Iterate over all arcs
            printEstimatableParameterEntries( parametersToEstimate );
            std::cout<<std::endl;
            for( int arc = 0; arc < numberOfArcs; arc++ )
            {
                std::cout <<"CASE: "<< test << " PARAMETER: " << parameterIndex << " ARC: "<< arc << std::endl;

                // Retrieve current arc variational equation results
                std::map< double, Eigen::MatrixXd > currentArcStateTransitionResults = stateTransitionResults.at( arc );
                std::map< double, Eigen::MatrixXd > currentArcSensitivityResults = sensitivityResults.at( arc );

                // Retrieve current arc up- and down-perturbed state results
                std::map< Time, Eigen::Matrix< long double, Eigen::Dynamic, 1 > > currentArcUpperturbedState = upperturbedResults.at( arc );
                std::map< Time, Eigen::Matrix< long double, Eigen::Dynamic, 1 > > currentArcDownperturbedState = downperturbedResults.at( arc );

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
                Eigen::VectorXd analyticalValue = matrixIterator->second.block( 0, matrixColumn, 6, 1 );
                Eigen::VectorXd numericalValue = ( ( currentArcUpperturbedState.at( matrixIterator->first ) - currentArcDownperturbedState.at( matrixIterator->first ) ).transpose( ) /
                                                   ( 2.0 * parameterPerturbation ) )
                                                         .cast< double >( );
//                Eigen::MatrixXd fullCombinedVariationalMatrix = stateTransitionMatrixInterface->getFullCombinedStateTransitionAndSensitivityMatrix( matrixIterator->first );

                // Remove initial unity value for testing resolution
                if( parameterIndex < 6 * numberOfArcs  )
                {
                    analyticalValue( matrixColumn ) -= 1.0 * arcLength;
                    numericalValue( matrixColumn ) -= 1.0 * arcLength;
                }
                if( ( arc == parameterIndex / 6 ) && useStateTransition || !useStateTransition )
                {
                    std::cout << "ANALYTICAL: " << analyticalValue.transpose( ) << std::endl;
                    std::cout << "NUMERICAL : " << numericalValue.transpose( ) << std::endl;
                    std::cout << "RATIO     : " << ( ( numericalValue - analyticalValue ).cwiseQuotient( analyticalValue ) ).transpose( ) << std::endl << std::endl;

                }
            }
        }
    }
}

//int main( )
////BOOST_AUTO_TEST_CASE( test_PanelledRadiationPressureEstimation )
//{
//    spice_interface::loadStandardSpiceKernels( );
//
//    for( int testCase = 0; testCase < 1; testCase++ )
//    {
//        // Define bodies in simulation
//        std::vector< std::string > bodyNames;
//        bodyNames.push_back( "Earth" );
//        bodyNames.push_back( "Sun" );
//
//        // Specify initial time
//        double initialEphemerisTime = double( 1.0E7 );
//        double finalEphemerisTime = initialEphemerisTime + 0.5 * 86400.0;
//
//        // Create bodies needed in simulation
//        BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, initialEphemerisTime - 3600.0, finalEphemerisTime + 3600.0 );
//        setSimpleRotationSettingsFromSpice( bodySettings, "Earth", initialEphemerisTime );
//
//        // Add Earth radiation pressure models
//        bodySettings.at( "Earth" )->radiationSourceModelSettings =
//                                        extendedRadiationSourceModelSettings({
//                                            albedoPanelRadiosityModelSettings(KnockeTypeSurfacePropertyDistributionModel::albedo_knocke, "Sun"),
//                                            delayedThermalPanelRadiosityModelSettings(KnockeTypeSurfacePropertyDistributionModel::emissivity_knocke, "Sun")
//                                        }, {3, 4});
//
//
//        bodySettings.addSettings( "Vehicle" );
//        bodySettings.at( "Vehicle" )->ephemerisSettings = constantEphemerisSettings( Eigen::Vector6d::Zero( ), "Earth", "ECLIPJ2000" );
//        bodySettings.at( "Vehicle" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );
//
//        SystemOfBodies bodies = createSystemOfBodies( bodySettings );
//
//        // Create spacecraft object.
//        bodies.at( "Vehicle" )->setConstantBodyMass( 400.0 );
//
//        Eigen::Vector7d rotationalStateVehicle;
//        rotationalStateVehicle.segment( 0, 4 ) =
//                linear_algebra::convertQuaternionToVectorFormat( Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ) );
//        rotationalStateVehicle.segment( 4, 3 ) = Eigen::Vector3d::Zero( );
//        bodies.at( "Vehicle" )
//                ->setRotationalEphemeris(
//                        std::make_shared< ConstantRotationalEphemeris >( rotationalStateVehicle, "ECLIPJ2000", "VehicleFixed" ) );
//        {
////            std::shared_ptr< FullPanelledBodySettings > panelSettings;
////            if( testCase < 2 )
////            {
////                std::vector< std::shared_ptr< BodyPanelSettings > > panelSettingsList;
////
////                double specular1 = 0.35, specular2 = 0.35, diffuse1 = 0.25, diffuse2 = 0.25;
////                if( testCase == 1 )
////                {
////                    specular1 += 0.1;
////                    specular2 -= 0.1;
////                    diffuse1 += 0.05;
////                    diffuse2 -= 0.05;
////                }
////                panelSettingsList.push_back( bodyPanelSettings( frameFixedPanelGeometry( Eigen::Vector3d::UnitX( ), 9.9 ),
////                                                                specularDiffuseBodyPanelReflectionLawSettings( specular1, diffuse1, false ),
////                                                                "Bus" ) );
////                panelSettingsList.push_back( bodyPanelSettings( frameFixedPanelGeometry( Eigen::Vector3d::UnitY( ), 9.9 ),
////                                                                specularDiffuseBodyPanelReflectionLawSettings( specular2, diffuse2, false ),
////                                                                "Bus" ) );
////                panelSettings = fullPanelledBodySettings( panelSettingsList );
////            }
////            else
////            {
////                panelSettings = bodyWingPanelledGeometry( 2., 3., 4., 0., 0.35, 0.25, 0.35, 0.25, false, false );
////            }
////            addBodyExteriorPanelledShape( panelSettings, "Vehicle", bodies );
//            bodies.at( "Vehicle" )
//                    ->setRadiationPressureTargetModels( { createRadiationPressureTargetModel(
//                            std::make_shared< CannonballRadiationPressureTargetModelSettings >( 10.0, 1.5 ), "Vehicle", bodies ) } );
//        }
//
//        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//        ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
//        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//        // Set accelerations on Vehicle that are to be taken into account.
//        SelectedAccelerationMap accelerationSettingsList;
//        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
//
//        accelerationsOfVehicle[ "Earth" ] = { sphericalHarmonicAcceleration( 2, 2 ), radiationPressureAcceleration( ) };
//        accelerationsOfVehicle[ "Sun" ] = { pointMassGravityAcceleration( ), radiationPressureAcceleration( ) };
//
//        accelerationSettingsList[ "Vehicle" ] = accelerationsOfVehicle;
//
//        // Set bodies for which initial state is to be estimated and integrated.
//        std::vector< std::string > bodiesToIntegrate = { "Vehicle" };
//        std::vector< std::string > centralBodies = { "Earth" };
//
//        // Create acceleration models
//        AccelerationMap accelerationModelMap =
//                createAccelerationModelsMap( bodies, accelerationSettingsList, bodiesToIntegrate, centralBodies );
//
//        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//        ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
//        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//        // Set Keplerian elements for Vehicle.
//        Eigen::Vector6d vehicleInitialStateInKeplerianElements;
//        vehicleInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7200.0E3;
//        vehicleInitialStateInKeplerianElements( eccentricityIndex ) = 0.05;
//        vehicleInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
//        vehicleInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = unit_conversions::convertDegreesToRadians( 235.7 );
//        vehicleInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = unit_conversions::convertDegreesToRadians( 23.4 );
//        vehicleInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );
//
//        // Set initial state.
//        double earthGravitationalParameter = getBodyGravitationalParameter( bodies, "Earth" );
//        Eigen::Matrix< double, 6, 1 > systemInitialState =
//                convertKeplerianToCartesianElements( vehicleInitialStateInKeplerianElements, earthGravitationalParameter );
//
//        vehicleInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7300.0E3;
//        vehicleInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 95.3 );
//        Eigen::Matrix< double, 6, 1 > systemInitialState2 =
//                convertKeplerianToCartesianElements( vehicleInitialStateInKeplerianElements, earthGravitationalParameter );
//
//        // Create integrator settings
//        std::shared_ptr< IntegratorSettings< double > > integratorSettings = rungeKuttaFixedStepSettings( 40.0, rungeKuttaFehlberg78 );
//
//        // Create propagator settings
//        std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettingsArc1 =
//                std::make_shared< TranslationalStatePropagatorSettings< double > >(
//                        centralBodies,
//                        accelerationModelMap,
//                        bodiesToIntegrate,
//                        systemInitialState,
//                        initialEphemerisTime,
//                        integratorSettings,
//                        propagationTimeTerminationSettings( finalEphemerisTime ) );
//
//        // Create propagator settings
//        double arcLength = 0.5 * 86400;
//        std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettingsArc2 =
//                std::make_shared< TranslationalStatePropagatorSettings< double > >(
//                        centralBodies,
//                        accelerationModelMap,
//                        bodiesToIntegrate,
//                        systemInitialState2,
//                        finalEphemerisTime,
//                        integratorSettings,
//                        propagationTimeTerminationSettings( finalEphemerisTime + arcLength ) );
//
//        std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsList;
//        propagatorSettingsList.push_back( propagatorSettingsArc1 );
//        propagatorSettingsList.push_back( propagatorSettingsArc2 );
//
//        std::shared_ptr< MultiArcPropagatorSettings< double > > propagatorSettings =
//                std::make_shared< MultiArcPropagatorSettings< double > >( propagatorSettingsList );
//
//
//        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//        ///////////////////////    DEFINE PARAMETERS THAT ARE TO BE ESTIMATED      ////////////////////////////////////////////
//        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//        // Define list of parameters to estimate.
//        std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
//        parameterNames = getInitialStateParameterSettings< double >( propagatorSettings, bodies );
////        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Vehicle", radiation_pressure_coefficient ) );
//        parameterNames.push_back( radiationPressureTargetDirectionScaling( "Vehicle", "Earth" ) );
//        parameterNames.push_back( radiationPressureTargetDirectionScaling( "Vehicle", "Earth" ) );
//        parameterNames.push_back( radiationPressureTargetPerpendicularDirectionScaling( "Vehicle", "Sun" ) );
//        parameterNames.push_back( radiationPressureTargetPerpendicularDirectionScaling( "Vehicle", "Sun" ) );
//
//
//        // Create parameters
//        std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
//                createParametersToEstimate< double, double >( parameterNames, bodies, propagatorSettings );
//
//        // Print identifiers and indices of parameters to terminal.
//        printEstimatableParameterEntries( parametersToEstimate );
//
//        std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;
//        LinkEnds linkEnds;
//        linkEnds[ observed_body ] = LinkEndId( "Vehicle", "" );
//        observationSettingsList.push_back( std::make_shared< ObservationModelSettings >( position_observable, linkEnds ) );
//
//        // Create orbit determination object.
//        OrbitDeterminationManager< double, double > orbitDeterminationManager =
//                OrbitDeterminationManager< double, double >( bodies, parametersToEstimate, observationSettingsList, propagatorSettings );
//
//        // Compute list of observation times.
//        std::vector< double > baseTimeList;
//        double observationTime = initialEphemerisTime + 1000.0;
//        double observationInterval = 5.0;
//        while( observationTime < finalEphemerisTime +  arcLength - 1000.0 )
//        {
//            baseTimeList.push_back( observationTime );
//            observationTime += observationInterval;
//        }
//
//        std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > measurementSimulationInput;
//        measurementSimulationInput.push_back( std::make_shared< TabulatedObservationSimulationSettings<> >(
//                position_observable, linkEnds, baseTimeList, observed_body ) );
//
//        // Simulate observations.
//        std::shared_ptr< ObservationCollection<> > observationsAndTimes = simulateObservations< double, double >(
//                measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );
//
//        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//        //////////////////    PERTURB PARAMETER VECTOR AND ESTIMATE PARAMETERS     ///////////        //
///////////////////////////////////
//        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//        // Perturb parameter estimate.
//        Eigen::Matrix< double, Eigen::Dynamic, 1 > initialParameterEstimate =
//                parametersToEstimate->template getFullParameterValues< double >( );
//        Eigen::Matrix< double, Eigen::Dynamic, 1 > truthParameters = initialParameterEstimate;
//        Eigen::Matrix< double, Eigen::Dynamic, 1 > parameterPerturbation =
//                Eigen::Matrix< double, Eigen::Dynamic, 1 >::Zero( truthParameters.rows( ) );
//
//        // Perturbe initial state estimate.
//        parameterPerturbation.segment( 0, 3 ) = Eigen::Vector3d::Constant( 1.0 );
//        parameterPerturbation.segment( 3, 3 ) = Eigen::Vector3d::Constant( 1.E-3 );
//        parameterPerturbation.segment( 6, 3 ) = Eigen::Vector3d::Constant( 1.0 );
//        parameterPerturbation.segment( 9, 3 ) = Eigen::Vector3d::Constant( 1.E-3 );
//        for( int i = 12; i < truthParameters.rows( ); i++ )
//        {
//            parameterPerturbation.segment( i, 1 ) = Eigen::Vector1d::Constant( 0.1 );
//        }
////        parameterPerturbation.segment( 13, 1 ) = Eigen::Vector1d::Constant( 0.1 );
//
//        initialParameterEstimate += parameterPerturbation;
//        parametersToEstimate->resetParameterValues( initialParameterEstimate );
//        Eigen::Matrix< double, Eigen::Dynamic, 1 > testValues = parametersToEstimate->template getFullParameterValues< double >( );
////        std::cout << "TEST A:" << truthParameters.transpose( ) << std::endl;
////        std::cout << "TEST B:" << testValues.transpose( ) << std::endl;
//
//        // Define estimation input
//        std::shared_ptr< EstimationInput< double, double > > estimationInput =
//                std::make_shared< EstimationInput< double, double > >( observationsAndTimes );
//
//        estimationInput->defineEstimationSettings( true, true, true, true, true );
//        estimationInput->setConvergenceChecker( std::make_shared< EstimationConvergenceChecker >( 4 ) );
//
//        // Perform estimation
//        std::shared_ptr< EstimationOutput< double > > estimationOutput = orbitDeterminationManager.estimateParameters( estimationInput );
//
//        double rmsResidual = linear_algebra::getVectorEntryRootMeanSquare( observationsAndTimes->getConcatenatedResiduals( ) );
//        BOOST_CHECK_SMALL( rmsResidual, 1.0E-4 );
//        Eigen::VectorXd parameterEstimate = estimationOutput->parameterEstimate_;
//        std::cout << "parameter estimate: " << ( parameterEstimate ).transpose( ) << std::endl;
//
//        Eigen::VectorXd estimationError = parameterEstimate - truthParameters;
////        for( int i = 0; i < 3; i++ )
////        {
////            BOOST_CHECK_SMALL( std::fabs( estimationError( i ) ), 2.0E-4 );
////            BOOST_CHECK_SMALL( std::fabs( estimationError( i + 3 ) ), 2.0E-7 );
////        }
////        BOOST_CHECK_SMALL( std::fabs( estimationError( 6 ) ), 2.0E-4 );
////        BOOST_CHECK_SMALL( std::fabs( estimationError( 6 ) ), 2.0E-4 );
//
//        Eigen::Matrix< double, Eigen::Dynamic, 1 > newTestValues = parametersToEstimate->template getFullParameterValues< double >( );
////        std::cout << "TEST C:" << newTestValues.transpose( ) << std::endl;
//
//        std::cout << "estimation error: " << ( estimationError ).transpose( ) << std::endl;
//        std::cout << "truth value: " << ( truthParameters ).transpose( ) << std::endl;
//
//        input_output::writeMatrixToFile( estimationOutput->getUnnormalizedDesignMatrix( ),
//                                         "multiArcPartials.dat", 16,
//                                         "/home/dominic/Downloads/" );
//
//
//    }
//}

//BOOST_AUTO_TEST_SUITE_END( )
//
//}  // namespace unit_tests
//
//}  // namespace tudat
