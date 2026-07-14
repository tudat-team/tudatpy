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
#include <algorithm>
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/createGravityField.h"
#include "tudat/simulation/environment_setup/createGravityFieldVariations.h"
#include <thread>

#include <limits>

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/basics/utilities.h"
#include "tudat/simulation/estimation_setup/createEstimatableParametersFactory.h"
#include "tudat/astro/propagators/propagateCovariance.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/simulation/estimation_setup/createObservationModelFactory.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"
#include "tudat/simulation/estimation_setup/podProcessing.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/directTidalTimeLag.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/tidalLoveNumber.h"
#include "tudat/astro/gravitation/basicSolidBodyTideGravityFieldVariations.h"

namespace tudat
{
namespace unit_tests
{
BOOST_AUTO_TEST_SUITE( test_tidal_propery_estimation )

// Using declarations.
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

//! Unit test to check if tidal time lag parameters are estimated correctly
BOOST_AUTO_TEST_CASE( test_DissipationParameterEstimation )
{
    // Load spice kernels
    spice_interface::loadStandardSpiceKernels( );

    // Define list of bodies
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Jupiter" );
    bodyNames.push_back( "Io" );
    bodyNames.push_back( "Europa" );

    std::vector< std::string > satelliteNames;
    satelliteNames.push_back( "Io" );
    satelliteNames.push_back( "Europa" );

    // Specify initial time
    double initialEphemerisTime = 0.0;
    double finalEphemerisTime = 0.25 * physical_constants::JULIAN_YEAR;
    double fixedStepSize = 450.0;
    double buffer = 10.0 * fixedStepSize;

    // Get default settings
    BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    Eigen::MatrixXd cosineCoefficients = Eigen::MatrixXd::Zero( 3, 3 );

    // Update gravity field settings with sh field, required by tidal models (reference radius)
    cosineCoefficients( 0, 0 ) = 1.0;
    Eigen::MatrixXd sineCoefficients = Eigen::MatrixXd::Zero( 3, 3 );
    for( unsigned int i = 0; i < bodyNames.size( ); i++ )
    {
        if( bodyNames.at( i ) != "Sun" )
        {
            bodySettings.at( bodyNames.at( i ) )->gravityFieldSettings =
                    std::make_shared< SphericalHarmonicsGravityFieldSettings >( getBodyGravitationalParameter( bodyNames.at( i ) ),
                                                                                getAverageRadius( bodyNames.at( i ) ),
                                                                                cosineCoefficients,
                                                                                sineCoefficients,
                                                                                "IAU_" + bodyNames.at( i ) );
        }
    }

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // Test estimation of Jupiter dissipation without and with frequency dependence,
    // First estimating the tidal time lag (tests 0 and 1) and then the tidal quality factor Q (tests 2 and 3)
    for( unsigned int test = 0; test < 4; test++ )
    {
        // Define accelerations acting on Io and Europa
        double satelliteLoveNumber = 1.0, satelliteTimeLag = 1000.0;
        double jupiterLoveNumber = 1.0, jupiterTimeLag = 100.0;

        double jupiterRotationRate = 2.0 * mathematical_constants::PI / ( 9.8 * 3600.0 );

        double ioPeriod = 1.769 * 86400.0;
        double synodicPeriodIo =
                2.0 * mathematical_constants::PI / ( 2.0 * std::fabs( jupiterRotationRate - 2.0 * mathematical_constants::PI / ioPeriod ) );
        double jupiterInvQ = std::tan( jupiterTimeLag / synodicPeriodIo * 2.0 * mathematical_constants::PI );
        double ioInvQ = std::tan( satelliteTimeLag / ioPeriod * 2.0 * mathematical_constants::PI );
        double europaPeriod = 3.551 * 86400.0;
        double synodicPeriodEuropa = 2.0 * mathematical_constants::PI /
                ( 2.0 * std::fabs( jupiterRotationRate - 2.0 * mathematical_constants::PI / europaPeriod ) );
        double europaInvQ = std::tan( satelliteTimeLag / europaPeriod * 2.0 * mathematical_constants::PI );

        SelectedAccelerationMap accelerationMap;
        for( unsigned int i = 0; i < satelliteNames.size( ); i++ )
        {
            std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfSatellite;
            accelerationsOfSatellite[ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
            accelerationsOfSatellite[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );

            if( test < 2 )
            {
                accelerationsOfSatellite[ "Jupiter" ].push_back( std::make_shared< DirectTidalDissipationAccelerationSettings >(
                        satelliteLoveNumber, satelliteTimeLag, false, false ) );
                accelerationsOfSatellite[ "Jupiter" ].push_back(
                        std::make_shared< DirectTidalDissipationAccelerationSettings >( jupiterLoveNumber, jupiterTimeLag, false, true ) );
            }
            else
            {
                if( satelliteNames.at( i ) == "Io" )
                {
                    accelerationsOfSatellite[ "Jupiter" ].push_back( std::make_shared< DirectTidalDissipationAccelerationSettings >(
                            satelliteLoveNumber, ioInvQ, ioPeriod, false, false ) );
                    accelerationsOfSatellite[ "Jupiter" ].push_back( std::make_shared< DirectTidalDissipationAccelerationSettings >(
                            jupiterLoveNumber, jupiterInvQ, synodicPeriodIo, false, true ) );
                }
                else if( satelliteNames.at( i ) == "Europa" )
                {
                    accelerationsOfSatellite[ "Jupiter" ].push_back( std::make_shared< DirectTidalDissipationAccelerationSettings >(
                            satelliteLoveNumber, europaInvQ, europaPeriod, false, false ) );
                    accelerationsOfSatellite[ "Jupiter" ].push_back( std::make_shared< DirectTidalDissipationAccelerationSettings >(
                            jupiterLoveNumber, jupiterInvQ, synodicPeriodEuropa, false, true ) );
                }
            }

            for( unsigned int j = 0; j < satelliteNames.size( ); j++ )
            {
                if( i != j )
                {
                    accelerationsOfSatellite[ satelliteNames.at( j ) ].push_back(
                            std::make_shared< AccelerationSettings >( point_mass_gravity ) );
                }
            }
            accelerationMap[ satelliteNames.at( i ) ] = accelerationsOfSatellite;
        }

        // Set bodies for which initial state is to be estimated and integrated.
        std::vector< std::string > bodiesToEstimate;
        bodiesToEstimate.push_back( "Io" );
        bodiesToEstimate.push_back( "Europa" );
        std::vector< std::string > centralBodies;
        centralBodies.push_back( "Jupiter" );
        centralBodies.push_back( "Jupiter" );
        AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodies, accelerationMap, bodiesToEstimate, centralBodies );

        std::shared_ptr< PropagatorSettings< double > > propagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >(
                        centralBodies,
                        accelerationModelMap,
                        bodiesToEstimate,
                        getInitialStatesOfBodies( bodiesToEstimate, centralBodies, bodies, initialEphemerisTime ),
                        finalEphemerisTime );

        // Set parameters that are to be estimated.
        std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
                getInitialStateParameterSettings< double >( propagatorSettings, bodies );

        if( test < 2 )
        {
            parameterNames.push_back(
                    std::make_shared< tudat::estimatable_parameters::DirectTidalTimeLagEstimatableParameterSettings >( "Io", "" ) );
            parameterNames.push_back(
                    std::make_shared< tudat::estimatable_parameters::DirectTidalTimeLagEstimatableParameterSettings >( "Europa", "" ) );
        }
        else
        {
            parameterNames.push_back(
                    std::make_shared< tudat::estimatable_parameters::InverseTidalQualityFactorEstimatableParameterSettings >( "Io", "" ) );
            parameterNames.push_back(
                    std::make_shared< tudat::estimatable_parameters::InverseTidalQualityFactorEstimatableParameterSettings >( "Europa",
                                                                                                                              "" ) );
        }

        if( test == 0 )
        {
            parameterNames.push_back(
                    std::make_shared< tudat::estimatable_parameters::DirectTidalTimeLagEstimatableParameterSettings >( "Jupiter", "" ) );
        }
        else if( test == 1 )
        {
            parameterNames.push_back(
                    std::make_shared< tudat::estimatable_parameters::DirectTidalTimeLagEstimatableParameterSettings >( "Jupiter", "Io" ) );
            parameterNames.push_back( std::make_shared< tudat::estimatable_parameters::DirectTidalTimeLagEstimatableParameterSettings >(
                    "Jupiter", "Europa" ) );
        }
        else if( test == 2 )
        {
            parameterNames.push_back(
                    std::make_shared< tudat::estimatable_parameters::InverseTidalQualityFactorEstimatableParameterSettings >( "Jupiter",
                                                                                                                              "" ) );
        }
        else if( test == 3 )
        {
            parameterNames.push_back(
                    std::make_shared< tudat::estimatable_parameters::InverseTidalQualityFactorEstimatableParameterSettings >( "Jupiter",
                                                                                                                              "Io" ) );
            parameterNames.push_back(
                    std::make_shared< tudat::estimatable_parameters::InverseTidalQualityFactorEstimatableParameterSettings >( "Jupiter",
                                                                                                                              "Europa" ) );
        }
        std::shared_ptr< tudat::estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
                createParametersToEstimate< double, double >( parameterNames, bodies, propagatorSettings );

        // Define links in simulation.
        std::vector< LinkDefinition > linkEnds;
        linkEnds.resize( 2 );
        linkEnds[ 0 ][ observed_body ] = std::make_pair< std::string, std::string >( "Io", "" );
        linkEnds[ 1 ][ observed_body ] = std::make_pair< std::string, std::string >( "Europa", "" );
        std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;
        observationSettingsList.push_back( std::make_shared< ObservationModelSettings >( position_observable, linkEnds[ 0 ] ) );
        observationSettingsList.push_back( std::make_shared< ObservationModelSettings >( position_observable, linkEnds[ 1 ] ) );

        // Define integrator and propagator settings.
        std::shared_ptr< IntegratorSettings<> > integratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettings<> >(
                0.0, fixedStepSize, rungeKuttaFehlberg78, fixedStepSize, fixedStepSize, 1.0, 1.0 );

        // Create orbit determination object.
        OrbitDeterminationManager< double, double > orbitDeterminationManager = OrbitDeterminationManager< double, double >(
                bodies, parametersToEstimate, observationSettingsList, integratorSettings, propagatorSettings );

        Eigen::VectorXd initialParameterEstimate = parametersToEstimate->template getFullParameterValues< double >( );

        // Test if tidal time lag parameters are correctly linked to acceleration models
        std::vector< std::shared_ptr< EstimatableParameter< double > > > dissipationEstimatedParameters =
                parametersToEstimate->getEstimatedDoubleParameters( );
        if( test == 0 || test == 2 )
        {
            BOOST_CHECK_EQUAL( dissipationEstimatedParameters.size( ), 3 );
        }
        else
        {
            BOOST_CHECK_EQUAL( dissipationEstimatedParameters.size( ), 4 );
        }
        for( unsigned int i = 0; i < dissipationEstimatedParameters.size( ); i++ )
        {
            if( test < 2 )
            {
                std::shared_ptr< DirectTidalTimeLag > currentTidalTimeLagParameter =
                        std::dynamic_pointer_cast< DirectTidalTimeLag >( dissipationEstimatedParameters.at( i ) );
                BOOST_CHECK_EQUAL( ( currentTidalTimeLagParameter == nullptr ), false );
                int numberOfAccelerationModels = currentTidalTimeLagParameter->getTidalAccelerationModels( ).size( );

                if( ( i == 2 ) && ( test == 0 ) )
                {
                    BOOST_CHECK_EQUAL( numberOfAccelerationModels, 2 );
                }
                else
                {
                    BOOST_CHECK_EQUAL( numberOfAccelerationModels, 1 );
                }
            }
            else
            {
                std::shared_ptr< InverseTidalQualityFactor > currentQualityFactorParameter =
                        std::dynamic_pointer_cast< InverseTidalQualityFactor >( dissipationEstimatedParameters.at( i ) );
                BOOST_CHECK_EQUAL( ( currentQualityFactorParameter == nullptr ), false );
                int numberOfAccelerationModels = currentQualityFactorParameter->getTidalAccelerationModels( ).size( );

                if( ( i == 2 ) && ( test == 2 ) )
                {
                    BOOST_CHECK_EQUAL( numberOfAccelerationModels, 2 );
                }
                else
                {
                    BOOST_CHECK_EQUAL( numberOfAccelerationModels, 1 );
                }
            }
        }

        // Define observation simulation times
        std::vector< double > observationTimes;
        double observationTime = initialEphemerisTime + 10.0 * fixedStepSize;
        while( observationTime < finalEphemerisTime - 10.0 * fixedStepSize )
        {
            observationTimes.push_back( observationTime );
            observationTime += 6.0 * 3600.0;
            ;
        }

        std::vector< std::shared_ptr< ObservationSimulationSettings<> > > measurementSimulationInput;
        measurementSimulationInput.push_back( std::make_shared< TabulatedObservationSimulationSettings<> >(
                position_observable, linkEnds[ 0 ], observationTimes, observed_body ) );
        measurementSimulationInput.push_back( std::make_shared< TabulatedObservationSimulationSettings<> >(
                position_observable, linkEnds[ 1 ], observationTimes, observed_body ) );

        // Simulate observations
        std::shared_ptr< ObservationCollection<> > observationsAndTimes = simulateObservations< double, double >(
                measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );

        // Perturb parameter values
        Eigen::VectorXd truthParameters = initialParameterEstimate;
        for( unsigned int i = 0; i < bodiesToEstimate.size( ); i++ )
        {
            initialParameterEstimate[ 0 + 6 * i ] += 1.0E0;
            initialParameterEstimate[ 1 + 6 * i ] += 1.0E0;
            initialParameterEstimate[ 2 + 6 * i ] += 1.0E0;
            initialParameterEstimate[ 3 + 6 * i ] += 1.0E-5;
            initialParameterEstimate[ 4 + 6 * i ] += 1.0E-5;
            initialParameterEstimate[ 5 + 6 * i ] += 1.0E-5;
        }
        for( unsigned int i = 6 * bodiesToEstimate.size( ); i < static_cast< unsigned int >( initialParameterEstimate.rows( ) ); i++ )
        {
            initialParameterEstimate[ i ] *= 10.0;
        }
        parametersToEstimate->resetParameterValues( initialParameterEstimate );

        // Estimate initial states and tidal parameters
        std::shared_ptr< EstimationInput< double, double > > estimationInput =
                std::make_shared< EstimationInput< double, double > >( observationsAndTimes );
        estimationInput->setConvergenceChecker( std::make_shared< EstimationConvergenceChecker >( 3 ) );
        estimationInput->applyFinalParameterCorrection_ = true;
        estimationInput->defineEstimationSettings( 1, 1, 1, 0, 1, 0 );
        std::shared_ptr< EstimationOutput< double > > estimationOutput = orbitDeterminationManager.estimateParameters( estimationInput );

        // Check if parameters are correctly estimated
        Eigen::VectorXd estimatedParametervalues =
                estimationOutput->parameterHistory_.at( estimationOutput->parameterHistory_.size( ) - 1 );

        for( unsigned int i = 0; i < 3; i++ )
        {
            BOOST_CHECK_SMALL( std::fabs( truthParameters( i ) - estimatedParametervalues( i ) ), 0.1 );
            BOOST_CHECK_SMALL( std::fabs( truthParameters( i + 6 ) - estimatedParametervalues( i + 6 ) ), 0.1 );
            BOOST_CHECK_SMALL( std::fabs( truthParameters( i + 3 ) - estimatedParametervalues( i + 3 ) ), 1.0E-6 );
            BOOST_CHECK_SMALL( std::fabs( truthParameters( i + 9 ) - estimatedParametervalues( i + 9 ) ), 1.0E-6 );
        }
        if( test < 2 )
        {
            BOOST_CHECK_SMALL( std::fabs( truthParameters( 12 ) - estimatedParametervalues( 12 ) ), 0.1 );
            BOOST_CHECK_SMALL( std::fabs( truthParameters( 13 ) - estimatedParametervalues( 13 ) ), 10.0 );
            BOOST_CHECK_SMALL( std::fabs( truthParameters( 14 ) - estimatedParametervalues( 14 ) ), 1.0E-3 );
            if( test == 1 )
            {
                BOOST_CHECK_SMALL( std::fabs( truthParameters( 15 ) - estimatedParametervalues( 15 ) ), 1.0E-2 );
            }
        }
        else
        {
            BOOST_CHECK_SMALL( std::fabs( truthParameters( 12 ) - estimatedParametervalues( 12 ) ), 2.0e-6 );
            BOOST_CHECK_SMALL( std::fabs( truthParameters( 13 ) - estimatedParametervalues( 13 ) ), 1.0e-4 );
            BOOST_CHECK_SMALL( std::fabs( truthParameters( 14 ) - estimatedParametervalues( 14 ) ), 1.0E-7 );
            if( test == 1 )
            {
                BOOST_CHECK_SMALL( std::fabs( truthParameters( 15 ) - estimatedParametervalues( 15 ) ), 1.0E-6 );
            }
        }

        std::cout << "Parameter error: " << ( truthParameters - estimatedParametervalues ).transpose( ) << std::endl;
    }
}

//! Function to get tidal deformation model for Earth
std::vector< std::shared_ptr< GravityFieldVariationSettings > > getEarthGravityFieldVariationSettings( )
{
    std::vector< std::shared_ptr< GravityFieldVariationSettings > > gravityFieldVariations;

    std::vector< std::string > deformingBodies;
    deformingBodies.push_back( "Moon" );

    std::map< int, std::vector< std::complex< double > > > loveNumbers;

    std::vector< std::complex< double > > degreeTwoLoveNumbers_;
    degreeTwoLoveNumbers_.push_back( std::complex< double >( 0.29525, -0.00087 ) );
    degreeTwoLoveNumbers_.push_back( std::complex< double >( 0.29525, -0.00087 ) );
    degreeTwoLoveNumbers_.push_back( std::complex< double >( 0.29525, -0.00087 ) );
    std::vector< std::complex< double > > degreeThreeLoveNumbers_;
    degreeThreeLoveNumbers_.push_back( std::complex< double >( 0.093, 0.0 ) );
    degreeThreeLoveNumbers_.push_back( std::complex< double >( 0.093, 0.0 ) );
    degreeThreeLoveNumbers_.push_back( std::complex< double >( 0.093, 0.0 ) );
    degreeThreeLoveNumbers_.push_back( std::complex< double >( 0.093, 0.0 ) );
    loveNumbers[ 2 ] = degreeTwoLoveNumbers_;
    loveNumbers[ 3 ] = degreeThreeLoveNumbers_;

    std::shared_ptr< GravityFieldVariationSettings > moonGravityFieldVariation =
            std::make_shared< BasicSolidBodyGravityFieldVariationSettings >( deformingBodies, loveNumbers );
    gravityFieldVariations.push_back( moonGravityFieldVariation );

    deformingBodies[ 0 ] = "Sun";
    std::shared_ptr< GravityFieldVariationSettings > sunSingleGravityFieldVariation =
            std::make_shared< BasicSolidBodyGravityFieldVariationSettings >( deformingBodies, loveNumbers );
    gravityFieldVariations.push_back( sunSingleGravityFieldVariation );

    return gravityFieldVariations;
}

//! Unit test to check if real/complex Love numbers at different degrees/orders/forcing frequencies are being correctly estimated.
BOOST_AUTO_TEST_CASE( test_LoveNumberEstimationFromOrbiterData )
{
    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Define bodies in simulation
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Moon" );
    bodyNames.push_back( "Mars" );

    // Specify initial time
    int numberOfDaysOfData = 5;
    double initialEphemerisTime = 1.0E7;
    double finalEphemerisTime = initialEphemerisTime + numberOfDaysOfData * 86400.0;

    // Create bodies needed in simulation
    BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, "Earth", "ECLIPJ2000" );
    bodySettings.at( "Earth" )->gravityFieldVariationSettings = getEarthGravityFieldVariationSettings( );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    bodies.createEmptyBody( "Vehicle" );

    // Set accelerations on Vehicle that are to be taken into account.
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 3, 3 ) );
    accelerationsOfVehicle[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );
    accelerationsOfVehicle[ "Moon" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );
    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;

    // Set bodies for which initial state is to be estimated and integrated.
    std::vector< std::string > bodiesToIntegrate;
    std::vector< std::string > centralBodies;
    bodiesToIntegrate.push_back( "Vehicle" );
    centralBodies.push_back( "Earth" );

    // Create acceleration models
    AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodies, accelerationMap, bodiesToIntegrate, centralBodies );

    // Set Keplerian and Cartesian elements for spacecraft.
    Eigen::Vector6d asterixInitialStateInKeplerianElements;
    asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 6600.0E3;
    asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.05;
    asterixInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
    asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = unit_conversions::convertDegreesToRadians( 235.7 );
    asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = unit_conversions::convertDegreesToRadians( 23.4 );
    asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );
    double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
    Eigen::Vector6d systemInitialState =
            convertKeplerianToCartesianElements( asterixInitialStateInKeplerianElements, earthGravitationalParameter );

    // Create propagator settings
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >(
                    centralBodies, accelerationModelMap, bodiesToIntegrate, systemInitialState, double( finalEphemerisTime ), cowell );

    // Create integrator settings
    std::shared_ptr< IntegratorSettings< double > > integratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettings< double > >(
            initialEphemerisTime, 60.0, CoefficientSets::rungeKuttaFehlberg78, 60.0, 60.0, 1.0, 1.0 );

    // Define link ends to use
    LinkDefinition linkEnds;
    linkEnds[ observed_body ] = std::make_pair< std::string, std::string >( "Vehicle", "" );

    // Define parameters to estimate (vehicle initial state and various Love number combinations)
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
            getInitialStateParameterSettings< double >( propagatorSettings, bodies );
    parameterNames.push_back( std::make_shared< SingleDegreeVariableTidalLoveNumberEstimatableParameterSettings >(
            "Earth", 2, std::vector< int >{ 2 }, "Moon", true ) );
    parameterNames.push_back( std::make_shared< SingleDegreeVariableTidalLoveNumberEstimatableParameterSettings >(
            "Earth", 2, std::vector< int >{ 2 }, "Sun", true ) );
    parameterNames.push_back( std::make_shared< SingleDegreeVariableTidalLoveNumberEstimatableParameterSettings >(
            "Earth", 2, std::vector< int >{ 0, 1 }, "Moon", false ) );
    parameterNames.push_back( std::make_shared< FullDegreeTidalLoveNumberEstimatableParameterSettings >( "Earth", 3, "Moon", true ) );

    // Create parameters
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double, double >( parameterNames, bodies );
    printEstimatableParameterEntries( parametersToEstimate );

    // Define observation settings
    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;
    observationSettingsList.push_back( std::make_shared< ObservationModelSettings >( position_observable, linkEnds ) );
    // Create orbit determination object.
    OrbitDeterminationManager< double, double > orbitDeterminationManager = OrbitDeterminationManager< double, double >(
            bodies, parametersToEstimate, observationSettingsList, integratorSettings, propagatorSettings );

    // define observation simulation times
    std::vector< double > baseTimeList;
    double currentObservationTime = initialEphemerisTime + 1000.0;
    double observationInterval = 60.0;
    while( currentObservationTime < finalEphemerisTime - 1000.0 )
    {
        baseTimeList.push_back( currentObservationTime );
        currentObservationTime += observationInterval;
    }

    std::vector< std::shared_ptr< ObservationSimulationSettings<> > > measurementSimulationInput;
    measurementSimulationInput.push_back(
            std::make_shared< TabulatedObservationSimulationSettings<> >( position_observable, linkEnds, baseTimeList, observed_body ) );

    // Simulate observations
    std::shared_ptr< ObservationCollection<> > observationsAndTimes = simulateObservations< double, double >(
            measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );

    // Perturb parameter estimate
    Eigen::VectorXd initialParameterEstimate = parametersToEstimate->template getFullParameterValues< double >( );
    Eigen::VectorXd truthParameters = initialParameterEstimate;
    Eigen::VectorXd parameterPerturbation = Eigen::VectorXd::Zero( truthParameters.rows( ) );
    parameterPerturbation.segment( 0, 3 ) = Eigen::Vector3d::Constant( 1.0 );
    parameterPerturbation.segment( 3, 3 ) = Eigen::Vector3d::Constant( 1.E-3 );

    for( int i = 6; i < parameterPerturbation.rows( ); i++ )
    {
        parameterPerturbation( i ) += 0.001;
    }
    initialParameterEstimate += parameterPerturbation;
    parametersToEstimate->resetParameterValues( initialParameterEstimate );

    // Define estimation input
    std::shared_ptr< EstimationInput< double, double > > estimationInput =
            std::make_shared< EstimationInput< double, double > >( observationsAndTimes );
    estimationInput->setConvergenceChecker( std::make_shared< EstimationConvergenceChecker >( 4 ) );

    // Perform estimation
    std::shared_ptr< EstimationOutput< double > > estimationOutput = orbitDeterminationManager.estimateParameters( estimationInput );

    // Check estimation results
    Eigen::VectorXd estimationError = estimationOutput->parameterEstimate_ - truthParameters;
    for( int i = 0; i < 3; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( estimationError( i ) ), 1.0E-4 );
        BOOST_CHECK_SMALL( std::fabs( estimationError( i + 3 ) ), 5.0E-8 );
    }

    for( int i = 0; i < estimationError.rows( ) - 6; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( estimationError( i + 6 ) ), 5.0E-6 );
    }

    std::cout << "Parameter error" << ( estimationError ).transpose( ) << std::endl;
}

//! Helper to create separate Moon/Sun basic solid-body tide models with optional distinct Love numbers.
std::vector< std::shared_ptr< GravityFieldVariationSettings > > getSeparateMoonSunEarthTideSettings(
        const std::complex< double > moonDegreeTwoLoveNumber,
        const std::complex< double > sunDegreeTwoLoveNumber )
{
    std::vector< std::shared_ptr< GravityFieldVariationSettings > > gravityFieldVariations;

    std::map< int, std::vector< std::complex< double > > > moonLoveNumbers;
    moonLoveNumbers[ 2 ] = std::vector< std::complex< double > >( 3, moonDegreeTwoLoveNumber );
    moonLoveNumbers[ 3 ] = std::vector< std::complex< double > >( 4, std::complex< double >( 0.093, 0.0 ) );
    gravityFieldVariations.push_back(
            std::make_shared< BasicSolidBodyGravityFieldVariationSettings >( std::vector< std::string >{ "Moon" }, moonLoveNumbers ) );

    std::map< int, std::vector< std::complex< double > > > sunLoveNumbers;
    sunLoveNumbers[ 2 ] = std::vector< std::complex< double > >( 3, sunDegreeTwoLoveNumber );
    sunLoveNumbers[ 3 ] = std::vector< std::complex< double > >( 4, std::complex< double >( 0.093, 0.0 ) );
    gravityFieldVariations.push_back(
            std::make_shared< BasicSolidBodyGravityFieldVariationSettings >( std::vector< std::string >{ "Sun" }, sunLoveNumbers ) );

    return gravityFieldVariations;
}

//! Unit test for multi-model full-degree Love-number parameter selection and reset semantics (Issue #734).
BOOST_AUTO_TEST_CASE( test_MultiModelFullDegreeLoveNumberParameterSetup )
{
    spice_interface::loadStandardSpiceKernels( );

    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Moon" );
    BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, "Earth", "ECLIPJ2000" );
    bodySettings.at( "Earth" )->gravityFieldVariationSettings =
            getSeparateMoonSunEarthTideSettings( std::complex< double >( 0.301, -0.001 ), std::complex< double >( 0.295, -0.002 ) );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    auto getBasicTideModels = [ &bodies ]( ) -> std::vector< std::shared_ptr< gravitation::BasicSolidBodyTideGravityFieldVariations > > {
        return utilities::dynamicCastSVectorToTVector< gravitation::GravityFieldVariations,
                                                       gravitation::BasicSolidBodyTideGravityFieldVariations >(
                bodies.at( "Earth" )->getGravityFieldVariationSet( )->getDirectTidalGravityFieldVariations( ) );
    };

    // Test A / G: explicit multi-model creation, both body-list orders, and combined-model compatibility.
    {
        std::vector< std::string > deformingBodiesMoonSun;
        deformingBodiesMoonSun.push_back( "Moon" );
        deformingBodiesMoonSun.push_back( "Sun" );
        std::vector< std::string > deformingBodiesSunMoon;
        deformingBodiesSunMoon.push_back( "Sun" );
        deformingBodiesSunMoon.push_back( "Moon" );

        std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNamesMoonSun;
        parameterNamesMoonSun.push_back(
                std::make_shared< FullDegreeTidalLoveNumberEstimatableParameterSettings >( "Earth", 2, deformingBodiesMoonSun, false ) );
        std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNamesSunMoon;
        parameterNamesSunMoon.push_back(
                std::make_shared< FullDegreeTidalLoveNumberEstimatableParameterSettings >( "Earth", 2, deformingBodiesSunMoon, false ) );

        std::shared_ptr< EstimatableParameterSet< double > > parameterSetMoonSun =
                createParametersToEstimate< double, double >( parameterNamesMoonSun, bodies );
        std::shared_ptr< EstimatableParameterSet< double > > parameterSetSunMoon =
                createParametersToEstimate< double, double >( parameterNamesSunMoon, bodies );

        std::shared_ptr< FullDegreeTidalLoveNumber > loveNumberMoonSun =
                std::dynamic_pointer_cast< FullDegreeTidalLoveNumber >( parameterSetMoonSun->getVectorParameters( ).begin( )->second );
        std::shared_ptr< FullDegreeTidalLoveNumber > loveNumberSunMoon =
                std::dynamic_pointer_cast< FullDegreeTidalLoveNumber >( parameterSetSunMoon->getVectorParameters( ).begin( )->second );
        BOOST_REQUIRE( loveNumberMoonSun != nullptr );
        BOOST_REQUIRE( loveNumberSunMoon != nullptr );
        BOOST_CHECK_EQUAL( loveNumberMoonSun->getParameterSize( ), 1 );
        BOOST_CHECK_EQUAL( loveNumberMoonSun->getGravityFieldVariationModels( ).size( ), 2 );
        BOOST_CHECK_EQUAL( loveNumberSunMoon->getGravityFieldVariationModels( ).size( ), 2 );
        const std::vector< std::string > selectedBodies = loveNumberMoonSun->getDeformingBodies( );
        BOOST_CHECK_EQUAL( selectedBodies.size( ), 2 );
        BOOST_CHECK( std::find( selectedBodies.begin( ), selectedBodies.end( ), "Moon" ) != selectedBodies.end( ) );
        BOOST_CHECK( std::find( selectedBodies.begin( ), selectedBodies.end( ), "Sun" ) != selectedBodies.end( ) );
        BOOST_CHECK( loveNumberMoonSun->getGravityFieldVariationModels( ) == loveNumberSunMoon->getGravityFieldVariationModels( ) );
    }

    // Existing combined model with both deforming bodies remains valid.
    {
        BodyListSettings combinedBodySettings = getDefaultBodySettings( bodyNames, "Earth", "ECLIPJ2000" );
        std::map< int, std::vector< std::complex< double > > > combinedLoveNumbers;
        combinedLoveNumbers[ 2 ] = std::vector< std::complex< double > >( 3, std::complex< double >( 0.3, -0.001 ) );
        std::vector< std::string > combinedDeformingBodies;
        combinedDeformingBodies.push_back( "Moon" );
        combinedDeformingBodies.push_back( "Sun" );
        std::vector< std::shared_ptr< GravityFieldVariationSettings > > combinedVariations;
        combinedVariations.push_back(
                std::make_shared< BasicSolidBodyGravityFieldVariationSettings >( combinedDeformingBodies, combinedLoveNumbers ) );
        combinedBodySettings.at( "Earth" )->gravityFieldVariationSettings = combinedVariations;
        SystemOfBodies combinedBodies = createSystemOfBodies( combinedBodySettings );

        std::vector< std::shared_ptr< EstimatableParameterSettings > > combinedParameterNames;
        combinedParameterNames.push_back(
                std::make_shared< FullDegreeTidalLoveNumberEstimatableParameterSettings >( "Earth", 2, combinedDeformingBodies, false ) );
        bool combinedCreationSucceeded = true;
        try
        {
            createParametersToEstimate< double, double >( combinedParameterNames, combinedBodies );
        }
        catch( const std::runtime_error& )
        {
            combinedCreationSucceeded = false;
        }
        BOOST_CHECK( combinedCreationSucceeded );
    }

    // Test B: real-only reset updates the shared real part and preserves each order's imaginary part.
    {
        std::vector< std::shared_ptr< gravitation::BasicSolidBodyTideGravityFieldVariations > > tideModels = getBasicTideModels( );
        BOOST_REQUIRE_EQUAL( tideModels.size( ), 2 );

        // Deliberately distinct per-order imaginary parts so a mean-imag rewrite would fail.
        std::vector< std::complex< double > > moonDegreeTwoLoveNumbers = { std::complex< double >( 0.301, -0.001 ),
                                                                           std::complex< double >( 0.301, -0.003 ),
                                                                           std::complex< double >( 0.301, -0.005 ) };
        std::vector< std::complex< double > > sunDegreeTwoLoveNumbers = { std::complex< double >( 0.295, -0.002 ),
                                                                          std::complex< double >( 0.295, -0.004 ),
                                                                          std::complex< double >( 0.295, -0.006 ) };
        tideModels.at( 0 )->resetLoveNumbersOfDegree( moonDegreeTwoLoveNumbers, 2 );
        tideModels.at( 1 )->resetLoveNumbersOfDegree( sunDegreeTwoLoveNumbers, 2 );

        std::vector< std::string > deformingBodiesMoonSun;
        deformingBodiesMoonSun.push_back( "Moon" );
        deformingBodiesMoonSun.push_back( "Sun" );
        std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
        parameterNames.push_back(
                std::make_shared< FullDegreeTidalLoveNumberEstimatableParameterSettings >( "Earth", 2, deformingBodiesMoonSun, false ) );
        std::shared_ptr< EstimatableParameterSet< double > > parameterSet =
                createParametersToEstimate< double, double >( parameterNames, bodies );
        std::shared_ptr< FullDegreeTidalLoveNumber > loveNumberParameter =
                std::dynamic_pointer_cast< FullDegreeTidalLoveNumber >( parameterSet->getVectorParameters( ).begin( )->second );
        BOOST_REQUIRE( loveNumberParameter != nullptr );
        BOOST_CHECK_EQUAL( loveNumberParameter->getParameterSize( ), 1 );

        bool wrongSizeRejected = false;
        try
        {
            loveNumberParameter->setParameterValue( Eigen::VectorXd::Zero( 2 ) );
        }
        catch( const std::runtime_error& )
        {
            wrongSizeRejected = true;
        }
        BOOST_CHECK( wrongSizeRejected );

        Eigen::VectorXd resetValue = Eigen::VectorXd::Zero( 1 );
        resetValue( 0 ) = 0.33;
        loveNumberParameter->setParameterValue( resetValue );

        for( int order = 0; order <= 2; order++ )
        {
            BOOST_CHECK_CLOSE_FRACTION( tideModels.at( 0 )->getLoveNumbersOfDegree( 2 ).at( order ).real( ), 0.33, 1.0E-14 );
            BOOST_CHECK_CLOSE_FRACTION( tideModels.at( 1 )->getLoveNumbersOfDegree( 2 ).at( order ).real( ), 0.33, 1.0E-14 );
            BOOST_CHECK_CLOSE_FRACTION( tideModels.at( 0 )->getLoveNumbersOfDegree( 2 ).at( order ).imag( ),
                                        moonDegreeTwoLoveNumbers.at( order ).imag( ),
                                        1.0E-14 );
            BOOST_CHECK_CLOSE_FRACTION( tideModels.at( 1 )->getLoveNumbersOfDegree( 2 ).at( order ).imag( ),
                                        sunDegreeTwoLoveNumbers.at( order ).imag( ),
                                        1.0E-14 );
        }
        BOOST_CHECK( moonDegreeTwoLoveNumbers.at( 0 ).imag( ) != sunDegreeTwoLoveNumbers.at( 0 ).imag( ) );
        BOOST_CHECK( moonDegreeTwoLoveNumbers.at( 0 ).imag( ) != moonDegreeTwoLoveNumbers.at( 1 ).imag( ) );
    }

    // Complex shared reset updates both components in every selected model.
    {
        bodySettings.at( "Earth" )->gravityFieldVariationSettings =
                getSeparateMoonSunEarthTideSettings( std::complex< double >( 0.301, -0.001 ), std::complex< double >( 0.295, -0.002 ) );
        bodies = createSystemOfBodies( bodySettings );

        std::vector< std::string > deformingBodiesMoonSun;
        deformingBodiesMoonSun.push_back( "Moon" );
        deformingBodiesMoonSun.push_back( "Sun" );
        std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
        parameterNames.push_back(
                std::make_shared< FullDegreeTidalLoveNumberEstimatableParameterSettings >( "Earth", 2, deformingBodiesMoonSun, true ) );
        std::shared_ptr< EstimatableParameterSet< double > > parameterSet =
                createParametersToEstimate< double, double >( parameterNames, bodies );
        std::shared_ptr< FullDegreeTidalLoveNumber > loveNumberParameter =
                std::dynamic_pointer_cast< FullDegreeTidalLoveNumber >( parameterSet->getVectorParameters( ).begin( )->second );
        BOOST_REQUIRE( loveNumberParameter != nullptr );

        Eigen::VectorXd resetValue = Eigen::VectorXd::Zero( 2 );
        resetValue( 0 ) = 0.28;
        resetValue( 1 ) = -0.004;
        loveNumberParameter->setParameterValue( resetValue );

        std::vector< std::shared_ptr< gravitation::BasicSolidBodyTideGravityFieldVariations > > tideModels = getBasicTideModels( );
        for( unsigned int modelIndex = 0; modelIndex < tideModels.size( ); modelIndex++ )
        {
            for( int order = 0; order <= 2; order++ )
            {
                BOOST_CHECK_CLOSE_FRACTION( tideModels.at( modelIndex )->getLoveNumbersOfDegree( 2 ).at( order ).real( ), 0.28, 1.0E-14 );
                BOOST_CHECK_CLOSE_FRACTION( tideModels.at( modelIndex )->getLoveNumbersOfDegree( 2 ).at( order ).imag( ), -0.004, 1.0E-14 );
            }
        }
    }

    // Test C: empty deforming-body list selects all compatible basic tide models.
    {
        bodySettings.at( "Earth" )->gravityFieldVariationSettings =
                getSeparateMoonSunEarthTideSettings( std::complex< double >( 0.301, -0.001 ), std::complex< double >( 0.295, -0.002 ) );
        bodies = createSystemOfBodies( bodySettings );

        std::vector< std::string > emptyDeformingBodies;
        std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
        parameterNames.push_back(
                std::make_shared< FullDegreeTidalLoveNumberEstimatableParameterSettings >( "Earth", 2, emptyDeformingBodies, false ) );
        std::shared_ptr< EstimatableParameterSet< double > > parameterSet =
                createParametersToEstimate< double, double >( parameterNames, bodies );
        std::shared_ptr< FullDegreeTidalLoveNumber > loveNumberParameter =
                std::dynamic_pointer_cast< FullDegreeTidalLoveNumber >( parameterSet->getVectorParameters( ).begin( )->second );
        BOOST_REQUIRE( loveNumberParameter != nullptr );
        BOOST_CHECK_EQUAL( loveNumberParameter->getGravityFieldVariationModels( ).size( ), 2 );

        Eigen::VectorXd resetValue = Eigen::VectorXd::Zero( 1 );
        resetValue( 0 ) = 0.27;
        loveNumberParameter->setParameterValue( resetValue );
        std::vector< std::shared_ptr< gravitation::BasicSolidBodyTideGravityFieldVariations > > tideModels = getBasicTideModels( );
        for( unsigned int modelIndex = 0; modelIndex < tideModels.size( ); modelIndex++ )
        {
            for( int order = 0; order <= 2; order++ )
            {
                BOOST_CHECK_CLOSE_FRACTION( tideModels.at( modelIndex )->getLoveNumbersOfDegree( 2 ).at( order ).real( ), 0.27, 1.0E-14 );
            }
        }
    }

    // Test D: subset selection only updates the requested model.
    {
        bodySettings.at( "Earth" )->gravityFieldVariationSettings =
                getSeparateMoonSunEarthTideSettings( std::complex< double >( 0.301, -0.001 ), std::complex< double >( 0.295, -0.002 ) );
        bodies = createSystemOfBodies( bodySettings );

        std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
        parameterNames.push_back( std::make_shared< FullDegreeTidalLoveNumberEstimatableParameterSettings >( "Earth", 2, "Moon", false ) );
        std::shared_ptr< EstimatableParameterSet< double > > parameterSet =
                createParametersToEstimate< double, double >( parameterNames, bodies );
        std::shared_ptr< FullDegreeTidalLoveNumber > loveNumberParameter =
                std::dynamic_pointer_cast< FullDegreeTidalLoveNumber >( parameterSet->getVectorParameters( ).begin( )->second );
        BOOST_REQUIRE( loveNumberParameter != nullptr );

        std::vector< std::shared_ptr< gravitation::BasicSolidBodyTideGravityFieldVariations > > tideModels = getBasicTideModels( );
        const double sunRealBefore = tideModels.at( 1 )->getLoveNumbersOfDegree( 2 ).at( 0 ).real( );
        Eigen::VectorXd resetValue = Eigen::VectorXd::Zero( 1 );
        resetValue( 0 ) = 0.35;
        loveNumberParameter->setParameterValue( resetValue );

        for( int order = 0; order <= 2; order++ )
        {
            BOOST_CHECK_CLOSE_FRACTION( tideModels.at( 0 )->getLoveNumbersOfDegree( 2 ).at( order ).real( ), 0.35, 1.0E-14 );
            BOOST_CHECK_CLOSE_FRACTION( tideModels.at( 1 )->getLoveNumbersOfDegree( 2 ).at( order ).real( ), sunRealBefore, 1.0E-14 );
        }
    }

    // Test E: invalid selections throw actionable errors.
    {
        bodySettings.at( "Earth" )->gravityFieldVariationSettings =
                getSeparateMoonSunEarthTideSettings( std::complex< double >( 0.301, -0.001 ), std::complex< double >( 0.295, -0.002 ) );
        bodies = createSystemOfBodies( bodySettings );

        bool missingBodyCaught = false;
        try
        {
            std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
            parameterNames.push_back(
                    std::make_shared< FullDegreeTidalLoveNumberEstimatableParameterSettings >( "Earth", 2, "Mars", false ) );
            createParametersToEstimate< double, double >( parameterNames, bodies );
        }
        catch( const std::runtime_error& error )
        {
            missingBodyCaught = true;
            BOOST_CHECK( std::string( error.what( ) ).find( "Mars" ) != std::string::npos );
        }
        BOOST_CHECK( missingBodyCaught );

        bool unavailableDegreeCaught = false;
        try
        {
            std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
            parameterNames.push_back(
                    std::make_shared< FullDegreeTidalLoveNumberEstimatableParameterSettings >( "Earth", 4, "Moon", false ) );
            createParametersToEstimate< double, double >( parameterNames, bodies );
        }
        catch( const std::runtime_error& )
        {
            unavailableDegreeCaught = true;
        }
        BOOST_CHECK( unavailableDegreeCaught );

        BodyListSettings combinedBodySettings = getDefaultBodySettings( bodyNames, "Earth", "ECLIPJ2000" );
        std::map< int, std::vector< std::complex< double > > > combinedLoveNumbers;
        combinedLoveNumbers[ 2 ] = std::vector< std::complex< double > >( 3, std::complex< double >( 0.3, 0.0 ) );
        std::vector< std::string > combinedDeformingBodies;
        combinedDeformingBodies.push_back( "Moon" );
        combinedDeformingBodies.push_back( "Sun" );
        std::vector< std::shared_ptr< GravityFieldVariationSettings > > combinedVariations;
        combinedVariations.push_back(
                std::make_shared< BasicSolidBodyGravityFieldVariationSettings >( combinedDeformingBodies, combinedLoveNumbers ) );
        combinedBodySettings.at( "Earth" )->gravityFieldVariationSettings = combinedVariations;
        SystemOfBodies combinedBodies = createSystemOfBodies( combinedBodySettings );

        bool partialOverlapCaught = false;
        try
        {
            std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
            parameterNames.push_back(
                    std::make_shared< FullDegreeTidalLoveNumberEstimatableParameterSettings >( "Earth", 2, "Moon", false ) );
            createParametersToEstimate< double, double >( parameterNames, combinedBodies );
        }
        catch( const std::runtime_error& error )
        {
            partialOverlapCaught = true;
            BOOST_CHECK( std::string( error.what( ) ).find( "partially overlap" ) != std::string::npos );
        }
        BOOST_CHECK( partialOverlapCaught );
    }

    // Degree-aware selection: independent Moon degree-2 / degree-3 models remain selectable.
    {
        BodyListSettings degreeSeparatedBodySettings = getDefaultBodySettings( bodyNames, "Earth", "ECLIPJ2000" );
        std::vector< std::shared_ptr< GravityFieldVariationSettings > > degreeSeparatedVariations;

        std::map< int, std::vector< std::complex< double > > > moonDegreeTwoLoveNumbers;
        moonDegreeTwoLoveNumbers[ 2 ] = std::vector< std::complex< double > >( 3, std::complex< double >( 0.301, -0.001 ) );
        degreeSeparatedVariations.push_back( std::make_shared< BasicSolidBodyGravityFieldVariationSettings >(
                std::vector< std::string >{ "Moon" }, moonDegreeTwoLoveNumbers ) );

        std::map< int, std::vector< std::complex< double > > > moonDegreeThreeLoveNumbers;
        moonDegreeThreeLoveNumbers[ 3 ] = std::vector< std::complex< double > >( 4, std::complex< double >( 0.093, 0.0 ) );
        degreeSeparatedVariations.push_back( std::make_shared< BasicSolidBodyGravityFieldVariationSettings >(
                std::vector< std::string >{ "Moon" }, moonDegreeThreeLoveNumbers ) );

        std::map< int, std::vector< std::complex< double > > > sunDegreeTwoLoveNumbers;
        sunDegreeTwoLoveNumbers[ 2 ] = std::vector< std::complex< double > >( 3, std::complex< double >( 0.295, -0.002 ) );
        degreeSeparatedVariations.push_back( std::make_shared< BasicSolidBodyGravityFieldVariationSettings >(
                std::vector< std::string >{ "Sun" }, sunDegreeTwoLoveNumbers ) );

        degreeSeparatedBodySettings.at( "Earth" )->gravityFieldVariationSettings = degreeSeparatedVariations;
        SystemOfBodies degreeSeparatedBodies = createSystemOfBodies( degreeSeparatedBodySettings );

        std::vector< std::shared_ptr< EstimatableParameterSettings > > degreeTwoMoonParameters;
        degreeTwoMoonParameters.push_back(
                std::make_shared< FullDegreeTidalLoveNumberEstimatableParameterSettings >( "Earth", 2, "Moon", false ) );
        std::shared_ptr< EstimatableParameterSet< double > > degreeTwoMoonParameterSet =
                createParametersToEstimate< double, double >( degreeTwoMoonParameters, degreeSeparatedBodies );
        std::shared_ptr< FullDegreeTidalLoveNumber > degreeTwoMoonLoveNumber = std::dynamic_pointer_cast< FullDegreeTidalLoveNumber >(
                degreeTwoMoonParameterSet->getVectorParameters( ).begin( )->second );
        BOOST_REQUIRE( degreeTwoMoonLoveNumber != nullptr );
        BOOST_CHECK_EQUAL( degreeTwoMoonLoveNumber->getGravityFieldVariationModels( ).size( ), 1 );
        BOOST_CHECK_EQUAL( degreeTwoMoonLoveNumber->getDeformingBodies( ).size( ), 1 );
        BOOST_CHECK_EQUAL( degreeTwoMoonLoveNumber->getDeformingBodies( ).at( 0 ), "Moon" );

        std::vector< std::string > emptyDeformingBodies;
        std::vector< std::shared_ptr< EstimatableParameterSettings > > degreeTwoEmptyParameters;
        degreeTwoEmptyParameters.push_back(
                std::make_shared< FullDegreeTidalLoveNumberEstimatableParameterSettings >( "Earth", 2, emptyDeformingBodies, false ) );
        std::shared_ptr< EstimatableParameterSet< double > > degreeTwoEmptyParameterSet =
                createParametersToEstimate< double, double >( degreeTwoEmptyParameters, degreeSeparatedBodies );
        std::shared_ptr< FullDegreeTidalLoveNumber > degreeTwoEmptyLoveNumber = std::dynamic_pointer_cast< FullDegreeTidalLoveNumber >(
                degreeTwoEmptyParameterSet->getVectorParameters( ).begin( )->second );
        BOOST_REQUIRE( degreeTwoEmptyLoveNumber != nullptr );
        BOOST_CHECK_EQUAL( degreeTwoEmptyLoveNumber->getGravityFieldVariationModels( ).size( ), 2 );

        std::vector< std::shared_ptr< EstimatableParameterSettings > > degreeThreeMoonParameters;
        degreeThreeMoonParameters.push_back(
                std::make_shared< FullDegreeTidalLoveNumberEstimatableParameterSettings >( "Earth", 3, "Moon", false ) );
        std::shared_ptr< EstimatableParameterSet< double > > degreeThreeMoonParameterSet =
                createParametersToEstimate< double, double >( degreeThreeMoonParameters, degreeSeparatedBodies );
        std::shared_ptr< FullDegreeTidalLoveNumber > degreeThreeMoonLoveNumber = std::dynamic_pointer_cast< FullDegreeTidalLoveNumber >(
                degreeThreeMoonParameterSet->getVectorParameters( ).begin( )->second );
        BOOST_REQUIRE( degreeThreeMoonLoveNumber != nullptr );
        BOOST_CHECK_EQUAL( degreeThreeMoonLoveNumber->getGravityFieldVariationModels( ).size( ), 1 );
    }

    // Empty list must reject duplicate body coverage for the same requested degree.
    {
        BodyListSettings ambiguousBodySettings = getDefaultBodySettings( bodyNames, "Earth", "ECLIPJ2000" );
        std::vector< std::shared_ptr< GravityFieldVariationSettings > > ambiguousVariations;
        std::map< int, std::vector< std::complex< double > > > firstMoonLoveNumbers;
        firstMoonLoveNumbers[ 2 ] = std::vector< std::complex< double > >( 3, std::complex< double >( 0.301, 0.0 ) );
        ambiguousVariations.push_back( std::make_shared< BasicSolidBodyGravityFieldVariationSettings >(
                std::vector< std::string >{ "Moon" }, firstMoonLoveNumbers ) );
        std::map< int, std::vector< std::complex< double > > > secondMoonLoveNumbers;
        secondMoonLoveNumbers[ 2 ] = std::vector< std::complex< double > >( 3, std::complex< double >( 0.295, 0.0 ) );
        ambiguousVariations.push_back( std::make_shared< BasicSolidBodyGravityFieldVariationSettings >(
                std::vector< std::string >{ "Moon" }, secondMoonLoveNumbers ) );
        ambiguousBodySettings.at( "Earth" )->gravityFieldVariationSettings = ambiguousVariations;
        SystemOfBodies ambiguousBodies = createSystemOfBodies( ambiguousBodySettings );

        bool emptyListAmbiguityCaught = false;
        try
        {
            std::vector< std::string > emptyDeformingBodies;
            std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
            parameterNames.push_back(
                    std::make_shared< FullDegreeTidalLoveNumberEstimatableParameterSettings >( "Earth", 2, emptyDeformingBodies, false ) );
            createParametersToEstimate< double, double >( parameterNames, ambiguousBodies );
        }
        catch( const std::runtime_error& error )
        {
            emptyListAmbiguityCaught = true;
            BOOST_CHECK( std::string( error.what( ) ).find( "Moon" ) != std::string::npos );
            BOOST_CHECK( std::string( error.what( ) ).find( "ambiguous" ) != std::string::npos );
        }
        BOOST_CHECK( emptyListAmbiguityCaught );

        bool explicitListAmbiguityCaught = false;
        try
        {
            std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
            parameterNames.push_back(
                    std::make_shared< FullDegreeTidalLoveNumberEstimatableParameterSettings >( "Earth", 2, "Moon", false ) );
            createParametersToEstimate< double, double >( parameterNames, ambiguousBodies );
        }
        catch( const std::runtime_error& error )
        {
            explicitListAmbiguityCaught = true;
            BOOST_CHECK( std::string( error.what( ) ).find( "Moon" ) != std::string::npos );
            BOOST_CHECK( std::string( error.what( ) ).find( "ambiguous" ) != std::string::npos );
        }
        BOOST_CHECK( explicitListAmbiguityCaught );
    }

    // Order-varying Love numbers also synchronize across multiple selected models.
    {
        bodySettings.at( "Earth" )->gravityFieldVariationSettings =
                getSeparateMoonSunEarthTideSettings( std::complex< double >( 0.301, -0.001 ), std::complex< double >( 0.295, -0.002 ) );
        bodies = createSystemOfBodies( bodySettings );

        std::vector< std::string > deformingBodiesMoonSun;
        deformingBodiesMoonSun.push_back( "Moon" );
        deformingBodiesMoonSun.push_back( "Sun" );
        std::vector< int > estimatedOrders;
        estimatedOrders.push_back( 1 );
        std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
        parameterNames.push_back( std::make_shared< SingleDegreeVariableTidalLoveNumberEstimatableParameterSettings >(
                "Earth", 2, estimatedOrders, deformingBodiesMoonSun, false ) );
        std::shared_ptr< EstimatableParameterSet< double > > parameterSet =
                createParametersToEstimate< double, double >( parameterNames, bodies );
        std::shared_ptr< SingleDegreeVariableTidalLoveNumber > loveNumberParameter =
                std::dynamic_pointer_cast< SingleDegreeVariableTidalLoveNumber >( parameterSet->getVectorParameters( ).begin( )->second );
        BOOST_REQUIRE( loveNumberParameter != nullptr );
        BOOST_CHECK_EQUAL( loveNumberParameter->getGravityFieldVariationModels( ).size( ), 2 );
        BOOST_CHECK_EQUAL( loveNumberParameter->getParameterSize( ), 1 );

        std::vector< std::shared_ptr< gravitation::BasicSolidBodyTideGravityFieldVariations > > tideModels = getBasicTideModels( );
        BOOST_REQUIRE_EQUAL( tideModels.size( ), 2 );
        const double moonOrderZeroBefore = tideModels.at( 0 )->getLoveNumbersOfDegree( 2 ).at( 0 ).real( );
        const double moonOrderTwoBefore = tideModels.at( 0 )->getLoveNumbersOfDegree( 2 ).at( 2 ).real( );
        const double sunOrderZeroBefore = tideModels.at( 1 )->getLoveNumbersOfDegree( 2 ).at( 0 ).real( );
        const double sunOrderTwoBefore = tideModels.at( 1 )->getLoveNumbersOfDegree( 2 ).at( 2 ).real( );

        Eigen::VectorXd resetValue = Eigen::VectorXd::Zero( 1 );
        resetValue( 0 ) = 0.42;
        loveNumberParameter->setParameterValue( resetValue );

        BOOST_CHECK_CLOSE_FRACTION( tideModels.at( 0 )->getLoveNumbersOfDegree( 2 ).at( 1 ).real( ), 0.42, 1.0E-14 );
        BOOST_CHECK_CLOSE_FRACTION( tideModels.at( 1 )->getLoveNumbersOfDegree( 2 ).at( 1 ).real( ), 0.42, 1.0E-14 );
        BOOST_CHECK_CLOSE_FRACTION( tideModels.at( 0 )->getLoveNumbersOfDegree( 2 ).at( 0 ).real( ), moonOrderZeroBefore, 1.0E-14 );
        BOOST_CHECK_CLOSE_FRACTION( tideModels.at( 0 )->getLoveNumbersOfDegree( 2 ).at( 2 ).real( ), moonOrderTwoBefore, 1.0E-14 );
        BOOST_CHECK_CLOSE_FRACTION( tideModels.at( 1 )->getLoveNumbersOfDegree( 2 ).at( 0 ).real( ), sunOrderZeroBefore, 1.0E-14 );
        BOOST_CHECK_CLOSE_FRACTION( tideModels.at( 1 )->getLoveNumbersOfDegree( 2 ).at( 2 ).real( ), sunOrderTwoBefore, 1.0E-14 );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
