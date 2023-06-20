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

#include "tudat/simulation/simulation.h"
#include "tudat/simulation/estimation.h"
#include <boost/test/unit_test.hpp>
#include "tudat/basics/testMacros.h"

namespace tudat
{
namespace unit_tests
{

//! Using declarations.
using namespace interpolators;
using namespace numerical_integrators;
using namespace spice_interface;
using namespace ephemerides;
using namespace simulation_setup;
using namespace basic_astrodynamics;
using namespace orbital_element_conversions;
using namespace propagators;
using namespace estimatable_parameters;

//! Test suite for astro functions.
BOOST_AUTO_TEST_SUITE( test_non_sequential_variational_equations )

BOOST_AUTO_TEST_CASE( testNonSequentialSingleArcVariationalEquations )
{
    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Jupiter" );
    bodyNames.push_back( "Io" );
    bodyNames.push_back( "Europa" );
    bodyNames.push_back( "Ganymede" );
    bodyNames.push_back( "Callisto" );

    // Specify initial time
    double initialEpoch = 0.0;
    double finalEpoch = 86400.0;
    double midArcEpoch = ( initialEpoch + finalEpoch ) / 2.0;

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
    double initialTimeIntegration = midArcEpoch;
    double timeStep = 3600.0;
    std::shared_ptr< numerical_integrators::IntegratorSettings< > > forwardIntegratorSettings =
            std::make_shared< numerical_integrators::RungeKuttaVariableStepSizeSettingsScalarTolerances< double > >(
                    initialTimeIntegration, timeStep, CoefficientSets::rungeKuttaFehlberg78,
                    timeStep, timeStep, 1.0e3, 1.0e3 );

    std::shared_ptr< numerical_integrators::IntegratorSettings< > > backwardIntegratorSettings =
            std::make_shared< numerical_integrators::RungeKuttaVariableStepSizeSettingsScalarTolerances< double > >(
                    initialTimeIntegration, - timeStep, CoefficientSets::rungeKuttaFehlberg78,
                    - timeStep, - timeStep, 1.0e3, 1.0e3 );

    // Define initial states
    Eigen::VectorXd midArcStatesMoons = propagators::getInitialStatesOfBodies( bodiesToPropagate, centralBodies, bodies, midArcEpoch );

    // Define dependent variables to save
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
    dependentVariables.push_back( std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
            basic_astrodynamics::mutual_spherical_harmonic_gravity, "Io", "Jupiter" ) );
    for ( unsigned int k = 0 ; k < bodiesToPropagate.size( ) ; k++ )
    {
        if ( bodiesToPropagate.at( k ) != "Io" )
        {
            dependentVariables.push_back( std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    basic_astrodynamics::mutual_spherical_harmonic_gravity, "Io", bodiesToPropagate.at( k ) ) );
        }
    }


    // Define propagator settings.
    std::shared_ptr< TranslationalStatePropagatorSettings< > > forwardPropagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< > >(
            centralBodies, accelerationsMap, bodiesToPropagate, midArcStatesMoons, finalEpoch, cowell, dependentVariables );
    std::shared_ptr< TranslationalStatePropagatorSettings< > > backwardPropagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< > >(
            centralBodies, accelerationsMap, bodiesToPropagate, midArcStatesMoons, initialEpoch, cowell, dependentVariables );

    // Define parameters to estimate
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    for( unsigned int i = 0; i < bodiesToPropagate.size( ); i++ )
    {
        parameterNames.push_back(
                std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
                        bodiesToPropagate.at( i ), midArcStatesMoons.segment( i * 6, 6 ), centralBodies.at( i ) ) );

        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >(
                bodiesToPropagate.at( i ), gravitational_parameter ) );

        parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                2, 0,  2, 2, bodiesToPropagate.at( i ), spherical_harmonics_cosine_coefficient_block ) );
        parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                2, 1,  2, 2, bodiesToPropagate.at( i ), spherical_harmonics_sine_coefficient_block ) );
    }

    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > forwardPropagationParametersToEstimate =
            createParametersToEstimate< double >( parameterNames, bodies, forwardPropagatorSettings );
    printEstimatableParameterEntries( forwardPropagationParametersToEstimate );

    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > backwardPropagationParametersToEstimate =
            createParametersToEstimate< double >( parameterNames, bodies, backwardPropagatorSettings );

    // Propagate variational equations (forward leg)
    bool setIntegratedResult = false;
    bool integrateConcurrently = true;
    SingleArcVariationalEquationsSolver< > forwardVariationalEquationsSolver = SingleArcVariationalEquationsSolver< >(
            bodies, forwardIntegratorSettings, forwardPropagatorSettings, forwardPropagationParametersToEstimate, integrateConcurrently,
            forwardIntegratorSettings, false, true, setIntegratedResult );

    // Propagate variational equations (backward leg)
    SingleArcVariationalEquationsSolver< > backwardVariationalEquationsSolver = SingleArcVariationalEquationsSolver< >(
            bodies, backwardIntegratorSettings, backwardPropagatorSettings, backwardPropagationParametersToEstimate, integrateConcurrently,
            backwardIntegratorSettings, false, true, setIntegratedResult );

    // Save propagations outputs
    std::map< double, Eigen::MatrixXd > forwardStmHistory = forwardVariationalEquationsSolver.getStateTransitionMatrixSolution( );
    std::map< double, Eigen::MatrixXd > forwardSemHistory = forwardVariationalEquationsSolver.getSensitivityMatrixSolution( );
    std::map< double, Eigen::MatrixXd > backwardStmHistory = backwardVariationalEquationsSolver.getStateTransitionMatrixSolution( );
    std::map< double, Eigen::MatrixXd > backwardSemHistory = backwardVariationalEquationsSolver.getSensitivityMatrixSolution( );

    //! Create settings for non-sequential propagation
    std::shared_ptr< NonSequentialPropagationTerminationSettings > terminationSettings = std::make_shared< NonSequentialPropagationTerminationSettings >(
            std::make_shared< PropagationTimeTerminationSettings >( finalEpoch ), std::make_shared< PropagationTimeTerminationSettings >( initialEpoch ) );
    std::shared_ptr< TranslationalStatePropagatorSettings< > > nonsequentialPropagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< > >(
            centralBodies, accelerationsMap, bodiesToPropagate, midArcStatesMoons, terminationSettings, cowell, dependentVariables );


    // Propagate variational equations (non-sequentially)
    SingleArcVariationalEquationsSolver< > nonsequentialVariationalEquationsSolver = SingleArcVariationalEquationsSolver< >(
            bodies, forwardIntegratorSettings, nonsequentialPropagatorSettings, forwardPropagationParametersToEstimate, integrateConcurrently,
            forwardIntegratorSettings, false, true, setIntegratedResult );

    std::map< double, Eigen::MatrixXd > nonsequentialStmHistory = nonsequentialVariationalEquationsSolver.getStateTransitionMatrixSolution( );
    std::map< double, Eigen::MatrixXd > nonsequentialSemHistory = nonsequentialVariationalEquationsSolver.getSensitivityMatrixSolution( );

    // Check consistency maps size
    BOOST_CHECK_EQUAL( nonsequentialStmHistory.size( ), forwardStmHistory.size( ) + backwardStmHistory.size( ) - 1 );
    BOOST_CHECK_EQUAL( nonsequentialSemHistory.size( ), forwardSemHistory.size( ) + backwardStmHistory.size( ) - 1 );

    // Check state transition and sensitivity matrices at initial time.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( backwardStmHistory.begin( )->second, nonsequentialStmHistory.begin( )->second, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( backwardSemHistory.begin( )->second, nonsequentialSemHistory.begin( )->second, std::numeric_limits< double >::epsilon( ) );

    // Check state transition and sensitivity matrices at final time.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( forwardStmHistory.rbegin( )->second, nonsequentialStmHistory.rbegin( )->second, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( forwardSemHistory.rbegin( )->second, nonsequentialSemHistory.rbegin( )->second, std::numeric_limits< double >::epsilon( ) );


}

BOOST_AUTO_TEST_CASE( testNonSequentialMultiArcVariationalEquations )
{
    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Jupiter" );
    bodyNames.push_back( "Io" );
    bodyNames.push_back( "Europa" );
    bodyNames.push_back( "Ganymede" );
    bodyNames.push_back( "Callisto" );

    // Specify initial times
    double initialEpoch = 0.0;
    double finalEpoch = 86400.0;
    std::vector< double > arcStartTimes = { 0.0, 8.0 * 3600.0, 16.0 * 3600.0 };
    std::vector< double > arcEndTimes = { 7.0 * 3600.0, 15.0 * 3600.0, 86400.0 };
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
            initialEpoch - 86400.0, finalEpoch + 86400.0, 3600.0, globalFrameOrigin, globalFrameOrientation );
    bodySettings.at( "Europa" )->ephemerisSettings = std::make_shared< InterpolatedSpiceEphemerisSettings >(
            initialEpoch - 86400.0, finalEpoch + 86400.0, 3600.0, globalFrameOrigin, globalFrameOrientation );
    bodySettings.at( "Ganymede" )->ephemerisSettings = std::make_shared< InterpolatedSpiceEphemerisSettings >(
            initialEpoch - 86400.0, finalEpoch + 86400.0, 3600.0, globalFrameOrigin, globalFrameOrientation );
    bodySettings.at( "Callisto" )->ephemerisSettings = std::make_shared< InterpolatedSpiceEphemerisSettings >(
            initialEpoch - 86400.0, finalEpoch + 86400.0, 3600.0, globalFrameOrigin, globalFrameOrientation );

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
    double timeStep = 2700.0;
    std::shared_ptr< numerical_integrators::IntegratorSettings< > > forwardIntegratorSettings =
            std::make_shared< numerical_integrators::RungeKuttaVariableStepSizeSettingsScalarTolerances< double > >(
                    TUDAT_NAN, timeStep, CoefficientSets::rungeKuttaFehlberg78, timeStep, timeStep, 1.0e3, 1.0e3 );

    std::shared_ptr< numerical_integrators::IntegratorSettings< > > backwardIntegratorSettings =
            std::make_shared< numerical_integrators::RungeKuttaVariableStepSizeSettingsScalarTolerances< double > >(
                    TUDAT_NAN, - timeStep, CoefficientSets::rungeKuttaFehlberg78, - timeStep, - timeStep, 1.0e3, 1.0e3 );

    // Define arc-wise initial states
    std::vector< Eigen::VectorXd > midArcStatesMoons;
    Eigen::VectorXd arcWiseMidStatesIo, arcWiseMidStatesEuropa, arcWiseMidStatesGanymede, arcWiseMidStatesCallisto;
    arcWiseMidStatesIo.resize( 6 * nbArcs );
    arcWiseMidStatesEuropa.resize( 6 * nbArcs );
    arcWiseMidStatesGanymede.resize( 6 * nbArcs );
    arcWiseMidStatesCallisto.resize( 6 * nbArcs );
    for ( unsigned int i = 0 ; i < nbArcs ; i++ )
    {
        midArcStatesMoons.push_back( propagators::getInitialStatesOfBodies( bodiesToPropagate, centralBodies, bodies, midArcTimes.at( i ) ) );
        arcWiseMidStatesIo.segment( i * 6, 6 ) = midArcStatesMoons[ i ].segment( 0, 6 );
        arcWiseMidStatesEuropa.segment( i * 6, 6 ) = midArcStatesMoons[ i ].segment( 6, 6 );
        arcWiseMidStatesGanymede.segment( i * 6, 6 ) = midArcStatesMoons[ i ].segment( 12, 6 );
        arcWiseMidStatesCallisto.segment( i * 6, 6 ) = midArcStatesMoons[ i ].segment( 18, 6 );
    }
    std::map< std::string, Eigen::VectorXd > concatenatedArcWiseStatesPerBody;
    concatenatedArcWiseStatesPerBody[ "Io" ] =  arcWiseMidStatesIo;
    concatenatedArcWiseStatesPerBody[ "Europa" ] =  arcWiseMidStatesEuropa;
    concatenatedArcWiseStatesPerBody[ "Ganymede" ] =  arcWiseMidStatesGanymede;
    concatenatedArcWiseStatesPerBody[ "Callisto" ] =  arcWiseMidStatesCallisto;

    // Define dependent variables to save
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
    dependentVariables.push_back( std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
            basic_astrodynamics::mutual_spherical_harmonic_gravity, "Io", "Jupiter" ) );
    for ( unsigned int k = 0 ; k < bodiesToPropagate.size( ) ; k++ )
    {
        if ( bodiesToPropagate.at( k ) != "Io" )
        {
            dependentVariables.push_back( std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    basic_astrodynamics::mutual_spherical_harmonic_gravity, "Io", bodiesToPropagate.at( k ) ) );
        }
    }


    // Define propagator settings.
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > forwardPropagationSettingsList, backwardPropagationSettingsList, nonsequentialPropagationSettingsList;
    for ( unsigned int i = 0 ; i < nbArcs ; i++ )
    {
        forwardPropagationSettingsList.push_back( std::make_shared< TranslationalStatePropagatorSettings< double > >(
                centralBodies, accelerationsMap, bodiesToPropagate, midArcStatesMoons.at( i ), midArcTimes.at( i ), forwardIntegratorSettings,
                std::make_shared< PropagationTimeTerminationSettings >( arcEndTimes.at( i ) ), cowell, dependentVariables ) );
        backwardPropagationSettingsList.push_back( std::make_shared< TranslationalStatePropagatorSettings< double > >(
                centralBodies, accelerationsMap, bodiesToPropagate, midArcStatesMoons.at( i ), midArcTimes.at( i ), backwardIntegratorSettings,
                std::make_shared< PropagationTimeTerminationSettings >( arcStartTimes.at( i ) ), cowell, dependentVariables ) );

        std::shared_ptr< NonSequentialPropagationTerminationSettings > terminationSettings = std::make_shared< NonSequentialPropagationTerminationSettings >(
                std::make_shared< PropagationTimeTerminationSettings >( arcEndTimes.at( i ) ),
                std::make_shared< PropagationTimeTerminationSettings >( arcStartTimes.at( i ) ) );
        std::shared_ptr< TranslationalStatePropagatorSettings< > > nonsequentialPropagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< > >(
                centralBodies, accelerationsMap, bodiesToPropagate, midArcStatesMoons.at( i ), midArcTimes.at( i ), forwardIntegratorSettings,
                terminationSettings, cowell, dependentVariables );
        nonsequentialPropagationSettingsList.push_back( nonsequentialPropagatorSettings );
    }
    std::shared_ptr< MultiArcPropagatorSettings< > > forwardPropagatorSettings = std::make_shared< MultiArcPropagatorSettings< > >( forwardPropagationSettingsList );
    std::shared_ptr< MultiArcPropagatorSettings< > > backwardPropagatorSettings = std::make_shared< MultiArcPropagatorSettings< > >( backwardPropagationSettingsList );

    std::shared_ptr< MultiArcPropagatorSettings< > > nonSequentialPropagatorSettings = std::make_shared< MultiArcPropagatorSettings< > >( nonsequentialPropagationSettingsList );


    // Define parameters to estimate
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    for( unsigned int i = 0; i < bodiesToPropagate.size( ); i++ )
    {
        parameterNames.push_back(
                std::make_shared< ArcWiseInitialTranslationalStateEstimatableParameterSettings< double > >(
                        bodiesToPropagate.at( i ), concatenatedArcWiseStatesPerBody.at( bodiesToPropagate.at( i ) ), midArcTimes, centralBodies.at( i ) ) );

        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >(
                bodiesToPropagate.at( i ), gravitational_parameter ) );

        parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                2, 0, 2, 2, bodiesToPropagate.at( i ), spherical_harmonics_cosine_coefficient_block ) );
        parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                2, 1, 2, 2, bodiesToPropagate.at( i ), spherical_harmonics_sine_coefficient_block ) );
    }

    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > forwardPropagationParametersToEstimate =
            createParametersToEstimate< double >( parameterNames, bodies, forwardPropagatorSettings );
    printEstimatableParameterEntries( forwardPropagationParametersToEstimate );

    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > backwardPropagationParametersToEstimate =
            createParametersToEstimate< double >( parameterNames, bodies, backwardPropagatorSettings );

    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > nonSequentialPropagationParametersToEstimate =
            createParametersToEstimate< double >( parameterNames, bodies, nonSequentialPropagatorSettings );

    // Propagate variational equations (forward leg)
    MultiArcVariationalEquationsSolver< > forwardVariationalEquationsSolver = MultiArcVariationalEquationsSolver< >(
            bodies, forwardPropagatorSettings, forwardPropagationParametersToEstimate, true );


    // Propagate variational equations (backward leg)
    MultiArcVariationalEquationsSolver< > backwardVariationalEquationsSolver = MultiArcVariationalEquationsSolver< >(
            bodies, backwardPropagatorSettings, backwardPropagationParametersToEstimate, true );

    // Propagate variational equations non-sequentially from mid-arc.
    MultiArcVariationalEquationsSolver< > nonSequentialVariationalEquationsSolver = MultiArcVariationalEquationsSolver< >(
            bodies, nonSequentialPropagatorSettings, nonSequentialPropagationParametersToEstimate, true );


    std::vector< std::vector< std::map< double, Eigen::MatrixXd > > > forwardVariationalHistory = forwardVariationalEquationsSolver.getNumericalVariationalEquationsSolution( );
    std::vector< std::vector< std::map< double, Eigen::MatrixXd > > > backwardVariationalHistory = backwardVariationalEquationsSolver.getNumericalVariationalEquationsSolution( );
    std::vector< std::vector< std::map< double, Eigen::MatrixXd > > > nonSequentialVariationalHistory = nonSequentialVariationalEquationsSolver.getNumericalVariationalEquationsSolution( );

    for ( unsigned int i = 0 ; i < arcStartTimes.size( ) ; i++ )
    {
        std::map< double, Eigen::MatrixXd > stmForward = forwardVariationalHistory.at( i )[ 0 ];
        std::map< double, Eigen::MatrixXd > stmBackward = backwardVariationalHistory.at( i )[ 0 ];
        std::map< double, Eigen::MatrixXd > stmNonSequential = nonSequentialVariationalHistory.at( i )[ 0 ];

        std::map< double, Eigen::MatrixXd > semForward = forwardVariationalHistory.at( i )[ 1 ];
        std::map< double, Eigen::MatrixXd > semBackward = backwardVariationalHistory.at( i )[ 1 ];
        std::map< double, Eigen::MatrixXd > semNonSequential = nonSequentialVariationalHistory.at( i )[ 1 ];

        BOOST_CHECK_EQUAL( stmNonSequential.size( ), stmForward.size( ) + stmBackward.size( ) - 1 );
        BOOST_CHECK_EQUAL( semNonSequential.size( ), semForward.size( ) + semBackward.size( ) - 1 );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( stmNonSequential.begin( )->second, stmBackward.begin( )->second,
                                           std::numeric_limits< double >::epsilon( ) );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( stmNonSequential.rbegin( )->second, stmForward.rbegin( )->second,
                                           std::numeric_limits< double >::epsilon( ) );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( semNonSequential.begin( )->second, semBackward.begin( )->second,
                                           std::numeric_limits< double >::epsilon( ) );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( semNonSequential.rbegin( )->second, semForward.rbegin( )->second,
                                           std::numeric_limits< double >::epsilon( ) );
    }

}

BOOST_AUTO_TEST_CASE( testNonSequentialHybridArcVariationalEquations )
{
    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Jupiter" );
    bodyNames.push_back( "Saturn" );
    bodyNames.push_back( "Io" );
    bodyNames.push_back( "Europa" );
    bodyNames.push_back( "Ganymede" );
    bodyNames.push_back( "Callisto" );

    // Specify initial times
    double initialEpoch = 0.0;
    double finalEpoch = 86400.0;
    double midSingleArc = ( finalEpoch + initialEpoch ) / 2.0;
    std::vector< double > arcStartTimes = { 0.0, 8.0 * 3600.0, 16.0 * 3600.0 };
    std::vector< double > arcEndTimes = { 7.0 * 3600.0, 15.0 * 3600.0, 86400.0 };

    std::vector< double > midArcTimes;
    for ( unsigned int i = 0 ; i < arcStartTimes.size( ) ; i++ )
    {
        midArcTimes.push_back( ( arcStartTimes.at( i ) + arcEndTimes.at( i ) ) / 2.0 );
    }

    unsigned int nbArcs = arcStartTimes.size( );

    std::string globalFrameOrientation = "ECLIPJ2000";
    std::string globalFrameOrigin = "Sun";

    // Create bodies.
    BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, initialEpoch - 86400.0, finalEpoch + 86400.0 );
    bodySettings.at( "Jupiter" )->ephemerisSettings = std::make_shared< InterpolatedSpiceEphemerisSettings >(
            initialEpoch - 86400.0, finalEpoch + 86400.0, 3600.0, globalFrameOrigin, globalFrameOrientation );
    bodySettings.at( "Io" )->ephemerisSettings = std::make_shared< InterpolatedSpiceEphemerisSettings >(
            initialEpoch - 86400.0, finalEpoch + 86400.0, 3600.0, globalFrameOrigin, globalFrameOrientation );
    bodySettings.at( "Europa" )->ephemerisSettings = std::make_shared< InterpolatedSpiceEphemerisSettings >(
            initialEpoch - 86400.0, finalEpoch + 86400.0, 3600.0, globalFrameOrigin, globalFrameOrientation );
    bodySettings.at( "Ganymede" )->ephemerisSettings = std::make_shared< InterpolatedSpiceEphemerisSettings >(
            initialEpoch - 86400.0, finalEpoch + 86400.0, 3600.0, globalFrameOrigin, globalFrameOrientation );
    bodySettings.at( "Callisto" )->ephemerisSettings = std::make_shared< InterpolatedSpiceEphemerisSettings >(
            initialEpoch - 86400.0, finalEpoch + 86400.0, 3600.0, globalFrameOrigin, globalFrameOrientation );

    bodySettings.at( "Io" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );
    bodySettings.at( "Europa" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );
    bodySettings.at( "Ganymede" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );
    bodySettings.at( "Callisto" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );


    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    //! Define bodies to propagate.
    std::vector< std::string > singleArcBodiesToPropagate = { "Jupiter" };
    std::vector< std::string > singleArcCentralBodies = { "Sun" };

    std::vector< std::string > multiArcBodiesToPropagate = { "Io", "Europa", "Ganymede", "Callisto" };
    std::vector< std::string > multiArcCentralBodies = { "Jupiter", "Jupiter", "Jupiter", "Jupiter" };

    // Set multi-arc accelerations.
    SelectedAccelerationMap multiArcAccelerationSettings;
    for ( unsigned int i = 0 ; i < multiArcBodiesToPropagate.size( ) ; i++ )
    {
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfSatellite;
        accelerationsOfSatellite[ "Jupiter" ].push_back( std::make_shared< MutualSphericalHarmonicAccelerationSettings >( 8, 0, 2, 2  )  );
        for ( unsigned int j = 0 ; j < multiArcBodiesToPropagate.size( ) ; j++ )
        {
            if ( i != j )
            {
                accelerationsOfSatellite[ multiArcBodiesToPropagate[ j ] ].push_back( std::make_shared< MutualSphericalHarmonicAccelerationSettings >( 2, 2, 2, 2, 8, 0 ) );
            }
        }
        multiArcAccelerationSettings[ multiArcBodiesToPropagate[ i ] ] = accelerationsOfSatellite;
    }

    basic_astrodynamics::AccelerationMap multiArcAccelerationsMap = createAccelerationModelsMap(
            bodies, multiArcAccelerationSettings, multiArcBodiesToPropagate, multiArcCentralBodies );

    // Set single-arc accelerations.
    SelectedAccelerationMap singleArcAccelerationSettings;
    singleArcAccelerationSettings[ "Jupiter" ][ "Sun" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );
    singleArcAccelerationSettings[ "Jupiter" ][ "Saturn" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );
    basic_astrodynamics::AccelerationMap singleArcAccelerationsMap = createAccelerationModelsMap(
            bodies, singleArcAccelerationSettings, singleArcBodiesToPropagate, singleArcCentralBodies );

    // Define integrator settings
    double timeStep = 3600.0;
    std::shared_ptr< numerical_integrators::IntegratorSettings< > > forwardIntegratorSettings =
            std::make_shared< numerical_integrators::RungeKuttaVariableStepSizeSettingsScalarTolerances< double > >(
                    midSingleArc, timeStep, CoefficientSets::rungeKuttaFehlberg78, timeStep, timeStep, 1.0e3, 1.0e3 );

    std::shared_ptr< numerical_integrators::IntegratorSettings< > > backwardIntegratorSettings =
            std::make_shared< numerical_integrators::RungeKuttaVariableStepSizeSettingsScalarTolerances< double > >(
                    midSingleArc, - timeStep, CoefficientSets::rungeKuttaFehlberg78, - timeStep, - timeStep, 1.0e3, 1.0e3 );

    // Define arc-wise initial states
    std::vector< Eigen::VectorXd > midArcStatesMoons;
    Eigen::VectorXd arcWiseMidStatesIo, arcWiseMidStatesEuropa, arcWiseMidStatesGanymede, arcWiseMidStatesCallisto;
    arcWiseMidStatesIo.resize( 6 * nbArcs );
    arcWiseMidStatesEuropa.resize( 6 * nbArcs );
    arcWiseMidStatesGanymede.resize( 6 * nbArcs );
    arcWiseMidStatesCallisto.resize( 6 * nbArcs );
    for ( unsigned int i = 0 ; i < nbArcs ; i++ )
    {
        midArcStatesMoons.push_back( propagators::getInitialStatesOfBodies( multiArcBodiesToPropagate, multiArcCentralBodies, bodies, midArcTimes.at( i ) ) );
        arcWiseMidStatesIo.segment( i * 6, 6 ) = midArcStatesMoons[ i ].segment( 0, 6 );
        arcWiseMidStatesEuropa.segment( i * 6, 6 ) = midArcStatesMoons[ i ].segment( 6, 6 );
        arcWiseMidStatesGanymede.segment( i * 6, 6 ) = midArcStatesMoons[ i ].segment( 12, 6 );
        arcWiseMidStatesCallisto.segment( i * 6, 6 ) = midArcStatesMoons[ i ].segment( 18, 6 );
    }
    std::map< std::string, Eigen::VectorXd > concatenatedArcWiseStatesPerBody;
    concatenatedArcWiseStatesPerBody[ "Io" ] =  arcWiseMidStatesIo;
    concatenatedArcWiseStatesPerBody[ "Europa" ] =  arcWiseMidStatesEuropa;
    concatenatedArcWiseStatesPerBody[ "Ganymede" ] =  arcWiseMidStatesGanymede;
    concatenatedArcWiseStatesPerBody[ "Callisto" ] =  arcWiseMidStatesCallisto;

    // Define single-arc initial states
    Eigen::VectorXd midArcStatesJupiter = propagators::getInitialStatesOfBodies( singleArcBodiesToPropagate, singleArcCentralBodies, bodies, midSingleArc );

    // Define multi-arc propagator settings lists.
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > forwardPropagationSettingsList, backwardPropagationSettingsList, nonSequentialPropagationSettingsList;
    for ( unsigned int i = 0 ; i < nbArcs ; i++ )
    {
        forwardPropagationSettingsList.push_back( std::make_shared< TranslationalStatePropagatorSettings< double > >(
                multiArcCentralBodies, multiArcAccelerationsMap, multiArcBodiesToPropagate, midArcStatesMoons.at( i ), midArcTimes.at( i ), forwardIntegratorSettings,
                std::make_shared< PropagationTimeTerminationSettings >( arcEndTimes.at( i ) ) ) );
        backwardPropagationSettingsList.push_back( std::make_shared< TranslationalStatePropagatorSettings< double > >(
                multiArcCentralBodies, multiArcAccelerationsMap, multiArcBodiesToPropagate, midArcStatesMoons.at( i ), midArcTimes.at( i ), backwardIntegratorSettings,
                std::make_shared< PropagationTimeTerminationSettings >( arcStartTimes.at( i ) ) ) );

        std::shared_ptr< NonSequentialPropagationTerminationSettings > terminationSettings = std::make_shared< NonSequentialPropagationTerminationSettings >(
                std::make_shared< PropagationTimeTerminationSettings >( arcEndTimes.at( i ) ),
                std::make_shared< PropagationTimeTerminationSettings >( arcStartTimes.at( i ) ) );
        std::shared_ptr< TranslationalStatePropagatorSettings< > > nonSequentialPropagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< > >(
                multiArcCentralBodies, multiArcAccelerationsMap, multiArcBodiesToPropagate, midArcStatesMoons.at( i ),
                midArcTimes.at( i ), forwardIntegratorSettings, terminationSettings );
        nonSequentialPropagationSettingsList.push_back( nonSequentialPropagatorSettings );
    }

    // Create single-arc propagator settings
    std::shared_ptr< TranslationalStatePropagatorSettings< > > forwardSingleArcPropagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< > >(
            singleArcCentralBodies, singleArcAccelerationsMap, singleArcBodiesToPropagate, midArcStatesJupiter, midSingleArc, forwardIntegratorSettings,
            std::make_shared< PropagationTimeTerminationSettings >( finalEpoch ) );
    std::shared_ptr< TranslationalStatePropagatorSettings< > > backwardSingleArcPropagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< > >(
            singleArcCentralBodies, singleArcAccelerationsMap, singleArcBodiesToPropagate, midArcStatesJupiter, midSingleArc, backwardIntegratorSettings,
            std::make_shared< PropagationTimeTerminationSettings >( initialEpoch ) );

    std::shared_ptr< NonSequentialPropagationTerminationSettings > terminationSettings = std::make_shared< NonSequentialPropagationTerminationSettings >(
            std::make_shared< PropagationTimeTerminationSettings >( finalEpoch ),
            std::make_shared< PropagationTimeTerminationSettings >( initialEpoch ) );
    std::shared_ptr< TranslationalStatePropagatorSettings< > > nonSequentialSingleArcPropagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< > >(
            singleArcCentralBodies, singleArcAccelerationsMap, singleArcBodiesToPropagate, midArcStatesJupiter, midSingleArc, forwardIntegratorSettings, terminationSettings );

    // Create hybrid-arc propagator settings
    std::shared_ptr< HybridArcPropagatorSettings< > > forwardPropagatorSettings = std::make_shared< HybridArcPropagatorSettings< > >(
            forwardSingleArcPropagatorSettings, std::make_shared< MultiArcPropagatorSettings< > >( forwardPropagationSettingsList ) );
    std::shared_ptr< HybridArcPropagatorSettings< > > backwardPropagatorSettings = std::make_shared< HybridArcPropagatorSettings< > >(
            backwardSingleArcPropagatorSettings, std::make_shared< MultiArcPropagatorSettings< > >( backwardPropagationSettingsList ) );

    std::shared_ptr< HybridArcPropagatorSettings< > > nonSequentialPropagatorSettings = std::make_shared< HybridArcPropagatorSettings< > >(
            nonSequentialSingleArcPropagatorSettings, std::make_shared< MultiArcPropagatorSettings< > >( nonSequentialPropagationSettingsList ) );

        // Define parameters to estimate
        std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
        parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
                singleArcBodiesToPropagate.at( 0 ), midArcStatesJupiter, singleArcCentralBodies.at( 0 ) ) );

        for( unsigned int i = 0; i < multiArcBodiesToPropagate.size( ); i++ )
        {
            parameterNames.push_back( std::make_shared< ArcWiseInitialTranslationalStateEstimatableParameterSettings< double > >(
                    multiArcBodiesToPropagate.at( i ), concatenatedArcWiseStatesPerBody.at( multiArcBodiesToPropagate.at( i ) ), midArcTimes, multiArcCentralBodies.at( i ) ) );

            parameterNames.push_back( std::make_shared< EstimatableParameterSettings >(
                    multiArcBodiesToPropagate.at( i ), gravitational_parameter ) );

            parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                    2, 0, 2, 2, multiArcBodiesToPropagate.at( i ), spherical_harmonics_cosine_coefficient_block ) );
            parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                    2, 1, 2, 2, multiArcBodiesToPropagate.at( i ), spherical_harmonics_sine_coefficient_block ) );
        }

        std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > forwardPropagationParametersToEstimate =
                createParametersToEstimate< double >( parameterNames, bodies, forwardPropagatorSettings );
        printEstimatableParameterEntries( forwardPropagationParametersToEstimate );

        std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > backwardPropagationParametersToEstimate =
                createParametersToEstimate< double >( parameterNames, bodies, backwardPropagatorSettings );

        std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > nonSequentialPropagationParametersToEstimate =
                createParametersToEstimate< double >( parameterNames, bodies, nonSequentialPropagatorSettings );

        // Propagate variational equations (forward leg)
        HybridArcVariationalEquationsSolver< > forwardVariationalEquationsSolver = HybridArcVariationalEquationsSolver< >(
                bodies, forwardPropagatorSettings, forwardPropagationParametersToEstimate, true );

        // Propagate variational equations (backward leg)
        HybridArcVariationalEquationsSolver< > backwardVariationalEquationsSolver = HybridArcVariationalEquationsSolver< >(
                bodies, backwardPropagatorSettings, backwardPropagationParametersToEstimate, true );

        // Propagate variational equations non-sequentially from mid-arc.
        HybridArcVariationalEquationsSolver< > nonSequentialVariationalEquationsSolver = HybridArcVariationalEquationsSolver< >(
                bodies, nonSequentialPropagatorSettings, nonSequentialPropagationParametersToEstimate, true );

        std::vector< std::map< double, Eigen::MatrixXd > > forwardSingleArcVariationalHistory =
                forwardVariationalEquationsSolver.getSingleArcSolver( )->getNumericalVariationalEquationsSolution( );
        std::vector< std::map< double, Eigen::MatrixXd > > backwardSingleArcVariationalHistory =
                backwardVariationalEquationsSolver.getSingleArcSolver( )->getNumericalVariationalEquationsSolution( );
        std::vector< std::map< double, Eigen::MatrixXd > > nonSequentialSingleArcVariationalHistory =
                nonSequentialVariationalEquationsSolver.getSingleArcSolver( )->getNumericalVariationalEquationsSolution( );

        std::vector< std::vector< std::map< double, Eigen::MatrixXd > > > forwardMultiArcVariationalHistory =
                forwardVariationalEquationsSolver.getMultiArcSolver( )->getNumericalVariationalEquationsSolution( );
        std::vector< std::vector< std::map< double, Eigen::MatrixXd > > > backwardMultiArcVariationalHistory =
                backwardVariationalEquationsSolver.getMultiArcSolver( )->getNumericalVariationalEquationsSolution( );
        std::vector< std::vector< std::map< double, Eigen::MatrixXd > > > nonSequentialMultiArcVariationalHistory =
                nonSequentialVariationalEquationsSolver.getMultiArcSolver( )->getNumericalVariationalEquationsSolution( );

        BOOST_CHECK_EQUAL( nonSequentialSingleArcVariationalHistory[ 0 ].size( ),
                           forwardSingleArcVariationalHistory[ 0 ].size( ) + backwardSingleArcVariationalHistory[ 0 ].size( ) - 1 );
        BOOST_CHECK_EQUAL( nonSequentialSingleArcVariationalHistory[ 1 ].size( ),
                           forwardSingleArcVariationalHistory[ 1 ].size( ) + backwardSingleArcVariationalHistory[ 1 ].size( ) - 1 );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( nonSequentialSingleArcVariationalHistory[ 0 ].begin( )->second, backwardSingleArcVariationalHistory[ 0 ].begin( )->second, std::numeric_limits< double >::epsilon( ) );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( nonSequentialSingleArcVariationalHistory[ 1 ].begin( )->second, backwardSingleArcVariationalHistory[ 1 ].begin( )->second, std::numeric_limits< double >::epsilon( ) );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( nonSequentialSingleArcVariationalHistory[ 0 ].rbegin( )->second, forwardSingleArcVariationalHistory[ 0 ].rbegin( )->second, std::numeric_limits< double >::epsilon( ) );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( nonSequentialSingleArcVariationalHistory[ 1 ].rbegin( )->second, forwardSingleArcVariationalHistory[ 1 ].rbegin( )->second, std::numeric_limits< double >::epsilon( ) );

        for ( unsigned int arc = 0 ; arc < arcStartTimes.size( ) ; arc++ )
        {
            BOOST_CHECK_EQUAL( nonSequentialMultiArcVariationalHistory.at( arc )[ 0 ].size( ),
                               forwardMultiArcVariationalHistory.at( arc )[ 0 ].size( ) + backwardMultiArcVariationalHistory.at( arc )[ 0 ].size( ) - 1 );
            BOOST_CHECK_EQUAL( nonSequentialMultiArcVariationalHistory.at( arc )[ 1 ].size( ),
                               forwardMultiArcVariationalHistory.at( arc )[ 1 ].size( ) + backwardMultiArcVariationalHistory.at( arc )[ 1 ].size( ) - 1 );

            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( nonSequentialMultiArcVariationalHistory.at( arc )[ 0 ].begin( )->second,
                                               backwardMultiArcVariationalHistory.at( arc )[ 0 ].begin( )->second, std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( nonSequentialMultiArcVariationalHistory.at( arc )[ 1 ].begin( )->second,
                                               backwardMultiArcVariationalHistory.at( arc )[ 1 ].begin( )->second, std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( nonSequentialMultiArcVariationalHistory.at( arc )[ 0 ].rbegin( )->second,
                                               forwardMultiArcVariationalHistory.at( arc )[ 0 ].rbegin( )->second, std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( nonSequentialMultiArcVariationalHistory.at( arc )[ 1 ].rbegin( )->second,
                                               forwardMultiArcVariationalHistory.at( arc )[ 1 ].rbegin( )->second, std::numeric_limits< double >::epsilon( ) );
        }

}

BOOST_AUTO_TEST_SUITE_END( )

}

}