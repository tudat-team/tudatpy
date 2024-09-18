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
#include "tudat/simulation/propagation_setup/dependentVariablesInterface.h"

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

BOOST_AUTO_TEST_SUITE( test_dependent_variables_interface )

BOOST_AUTO_TEST_CASE( testSingleArcDependentVariablesInterface )
{

    // Set simulation time settings.
    const double simulationStartEpoch = 0.0;
    const double simulationEndEpoch = tudat::physical_constants::JULIAN_DAY;

    // Load spice kernel.
    spice_interface::loadStandardSpiceKernels( );

    // Define body settings for simulation.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Mars" );

    // Create system of bodies
    BodyListSettings bodySettings = getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 300.0, simulationEndEpoch + 300.0 );
    std::shared_ptr< BodySettings > phobosSettings = std::make_shared< BodySettings >( );
    phobosSettings->ephemerisSettings = getDefaultEphemerisSettings( "Phobos" );
    phobosSettings->gravityFieldSettings = getDefaultGravityFieldSettings( "Phobos", simulationStartEpoch - 300.0, simulationEndEpoch + 300.0 );
    bodySettings.addSettings( phobosSettings, "Phobos" );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // Create spacecraft object.
    bodies.createEmptyBody( "AlienSpaceship" );
    bodies.at( "AlienSpaceship" )->setConstantBodyMass( 400.0 );

    // Create radiation pressure settings
    double referenceAreaRadiation = 4.0;
    double radiationPressureCoefficient = 1.2;
    std::vector< std::string > occultingBodies;
    occultingBodies.push_back( "Mars" );
    std::shared_ptr< RadiationPressureInterfaceSettings > radiationPressureSettings = std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
            "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

    // Create and set radiation pressure settings
    bodies.at( "AlienSpaceship" )->setRadiationPressureInterface(
            "Sun", createRadiationPressureInterface( radiationPressureSettings, "AlienSpaceship", bodies ) );


    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Define propagation settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerations;
    accelerations[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 5, 5 ) );
    accelerations[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );
    accelerations[ "Mars" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );
    accelerations[ "Phobos" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );
    accelerations[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::cannon_ball_radiation_pressure ) );

    accelerationMap[ "AlienSpaceship" ] = accelerations;
    bodiesToPropagate.push_back( "AlienSpaceship" );
    centralBodies.push_back( "Phobos" );

    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
            bodies, accelerationMap, bodiesToPropagate, centralBodies );

    // Set Keplerian elements for spacecraft.
    Eigen::Vector6d initialKeplerianElements;
    initialKeplerianElements( semiMajorAxisIndex ) = 100.0E3;
    initialKeplerianElements( eccentricityIndex ) = 0.0;
    initialKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 26.04 );
    initialKeplerianElements( argumentOfPeriapsisIndex ) = unit_conversions::convertDegreesToRadians( 0.0 );
    initialKeplerianElements( longitudeOfAscendingNodeIndex ) = unit_conversions::convertDegreesToRadians( 0.0 );
    initialKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 0.0 );

    double phobosGravitationalParameter = bodies.at( "Phobos" )->getGravityFieldModel( )->getGravitationalParameter( );
    const Eigen::Vector6d initialState = convertKeplerianToCartesianElements(
            initialKeplerianElements, phobosGravitationalParameter );

    // Define list of dependent variables to save.
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;

    dependentVariables.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
            total_acceleration_dependent_variable, "AlienSpaceship" ) );
    dependentVariables.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
            keplerian_state_dependent_variable, "AlienSpaceship", "Phobos" ) );

    const double fixedStepSize = 10.0;
    std::shared_ptr< IntegratorSettings< > > integratorSettings = rungeKutta4Settings( fixedStepSize );

    std::shared_ptr< PropagationTimeTerminationSettings > terminationSettings =
            std::make_shared< propagators::PropagationTimeTerminationSettings >( simulationEndEpoch );
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >(
            centralBodies, accelerationModelMap, bodiesToPropagate, initialState, simulationStartEpoch,
            integratorSettings, terminationSettings, cowell, dependentVariables );
    propagatorSettings->getOutputSettings()->setUpdateDependentVariableInterpolator( true );

    // Define single-arc dynamics simulator
    SingleArcDynamicsSimulator< > simulator = SingleArcDynamicsSimulator< >(
            bodies, propagatorSettings );

    // Retrieve dependent variables history.
    std::map< double, Eigen::VectorXd > dependentVariablesHistory = simulator.getDependentVariableHistory( );

    // Retrieve dependent variables interface.
    std::shared_ptr< SingleArcDependentVariablesInterface< > > dependentVariablesInterface =
            std::dynamic_pointer_cast< SingleArcDependentVariablesInterface< > >( simulator.getDependentVariablesInterface( ) );

    // Create dependent variables interpolator.
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > dependentVariablesInterpolator
            = std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd > >(
                    utilities::createVectorFromMapKeys< Eigen::VectorXd, double >( dependentVariablesHistory ),
                    utilities::createVectorFromMapValues< Eigen::VectorXd, double >( dependentVariablesHistory ), 4 );

    std::vector< double > testEpochs;
    testEpochs.push_back( simulationStartEpoch );
    testEpochs.push_back( ( simulationEndEpoch -  simulationStartEpoch ) / 4.0 );
    testEpochs.push_back( ( simulationEndEpoch -  simulationStartEpoch ) / 2.0 );
    testEpochs.push_back( 3.0 * ( simulationEndEpoch -  simulationStartEpoch ) / 4.0 );
    testEpochs.push_back( simulationEndEpoch );

    // Check consistency between interpolator results and interface results.
    for ( unsigned int i = 0 ; i < testEpochs.size( ) ; i++ )
    {
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( dependentVariablesInterpolator->interpolate( testEpochs[ i ] ),
                                           dependentVariablesInterface->getDependentVariables( testEpochs[ i ] ),
                                           std::numeric_limits< double >::epsilon( ) );
    }


    std::map< double, Eigen::VectorXd > totalAccelerationHistory;
    for ( std::map< double, Eigen::VectorXd >::iterator itr = dependentVariablesHistory.begin( ) ; itr != dependentVariablesHistory.end( ) ; itr++ )
    {
        totalAccelerationHistory[ itr->first ] = itr->second.segment( 0, 3 );
    }

    // Create total acceleration history interpolator.
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > totalAccelerationInterpolator
            = std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd > >(
                    utilities::createVectorFromMapKeys< Eigen::VectorXd, double >( totalAccelerationHistory ),
                    utilities::createVectorFromMapValues< Eigen::VectorXd, double >( totalAccelerationHistory ), 4 );

    // Total acceleration dependent variable settings.
    std::shared_ptr< SingleDependentVariableSaveSettings > totalAccelerationDependentVariable
            = std::make_shared< SingleDependentVariableSaveSettings >( total_acceleration_dependent_variable, "AlienSpaceship" );

    // Check consistency between interpolator results and interface results, for a single dependent variable.
    for ( unsigned int i = 0 ; i < testEpochs.size( ) ; i++ )
    {
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( totalAccelerationInterpolator->interpolate( testEpochs[ i ] ),
                                           dependentVariablesInterface->getSingleDependentVariable( totalAccelerationDependentVariable, testEpochs[ i ] ),
                                           std::numeric_limits< double >::epsilon( ) );
    }
}


BOOST_AUTO_TEST_CASE( testMultiArcDependentVariablesInterface )
{

    // Set simulation time settings.
    const double simulationStartEpoch = 0.0;
    const double simulationEndEpoch = 3.0 * tudat::physical_constants::JULIAN_DAY;

    // Load spice kernel.
    spice_interface::loadStandardSpiceKernels( );

    // Define body settings for simulation.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Mars" );

    // Create system of bodies
    BodyListSettings bodySettings = getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 300.0, simulationEndEpoch + 300.0 );
    std::shared_ptr< BodySettings > phobosSettings = std::make_shared< BodySettings >( );
    phobosSettings->ephemerisSettings = getDefaultEphemerisSettings( "Phobos" );
    phobosSettings->gravityFieldSettings = getDefaultGravityFieldSettings( "Phobos", simulationStartEpoch - 300.0, simulationEndEpoch + 300.0 );
    bodySettings.addSettings( phobosSettings, "Phobos" );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // Create spacecraft object.
    bodies.createEmptyBody( "AlienSpaceship" );
    bodies.at( "AlienSpaceship" )->setConstantBodyMass( 400.0 );
    bodies.at( "AlienSpaceship" )->setEphemeris( std::make_shared< MultiArcEphemeris >(
            std::map< double, std::shared_ptr< Ephemeris > >( ), "Phobos", "ECLIPJ2000" ) );

    // Create radiation pressure settings
    double referenceAreaRadiation = 4.0;
    double radiationPressureCoefficient = 1.2;
    std::vector< std::string > occultingBodies;
    occultingBodies.push_back( "Mars" );
    std::shared_ptr< RadiationPressureInterfaceSettings > radiationPressureSettings = std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
            "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

    // Create and set radiation pressure settings
    bodies.at( "AlienSpaceship" )->setRadiationPressureInterface(
            "Sun", createRadiationPressureInterface( radiationPressureSettings, "AlienSpaceship", bodies ) );


    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Define propagation settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerations;
    accelerations[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 5, 5 ) );
    accelerations[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );
    accelerations[ "Mars" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );
    accelerations[ "Phobos" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );
    accelerations[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::cannon_ball_radiation_pressure ) );

    accelerationMap[ "AlienSpaceship" ] = accelerations;
    bodiesToPropagate.push_back( "AlienSpaceship" );
    centralBodies.push_back( "Phobos" );

    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
            bodies, accelerationMap, bodiesToPropagate, centralBodies );

    // Set Keplerian elements for spacecraft.
    Eigen::Vector6d initialKeplerianElements;
    initialKeplerianElements( semiMajorAxisIndex ) = 100.0E3;
    initialKeplerianElements( eccentricityIndex ) = 0.0;
    initialKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 26.04 );
    initialKeplerianElements( argumentOfPeriapsisIndex ) = unit_conversions::convertDegreesToRadians( 0.0 );
    initialKeplerianElements( longitudeOfAscendingNodeIndex ) = unit_conversions::convertDegreesToRadians( 0.0 );
    initialKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 0.0 );

    std::vector< Eigen::Vector6d > arcWiseKeplerianStates;
    arcWiseKeplerianStates.push_back( initialKeplerianElements );
    initialKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 60.0 );
    arcWiseKeplerianStates.push_back( initialKeplerianElements );
    initialKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 120.0 );
    arcWiseKeplerianStates.push_back( initialKeplerianElements );

    unsigned int nbArcs = arcWiseKeplerianStates.size( );

    std::vector< double > arcStartTimes, arcEndTimes;
    std::vector< Eigen::Vector6d > arcWiseStates;
    double phobosGravitationalParameter = bodies.at( "Phobos" )->getGravityFieldModel( )->getGravitationalParameter( );
    for ( unsigned int i = 0 ; i < nbArcs ; i++ )
    {
        arcWiseStates.push_back( convertKeplerianToCartesianElements( arcWiseKeplerianStates.at( i ), phobosGravitationalParameter ) );
        arcStartTimes.push_back( simulationStartEpoch + i * physical_constants::JULIAN_DAY );
        arcEndTimes.push_back( simulationStartEpoch + ( i + 1 ) * physical_constants::JULIAN_DAY - 3600.0 );
    }

    // Define list of dependent variables to save.
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
    dependentVariables.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
            total_acceleration_dependent_variable, "AlienSpaceship" ) );
    dependentVariables.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
            keplerian_state_dependent_variable, "AlienSpaceship", "Phobos" ) );

    const double fixedStepSize = 10.0;
    std::shared_ptr< IntegratorSettings< > > integratorSettings = rungeKutta4Settings( fixedStepSize );

    std::vector< std::shared_ptr< SingleArcPropagatorSettings< > > > propagatorSettingsList;
    for ( unsigned int i = 0 ; i < nbArcs ; i++ )
    {
        propagatorSettingsList.push_back( std::make_shared< TranslationalStatePropagatorSettings< double > >(
                centralBodies, accelerationModelMap, bodiesToPropagate, arcWiseStates.at( i ), arcStartTimes.at( i ), integratorSettings,
                std::make_shared< propagators::PropagationTimeTerminationSettings >( arcEndTimes.at( i ) ), cowell, dependentVariables ) );
    }


    std::shared_ptr< MultiArcPropagatorSettings< > > multiArcPropagatorSettings =
            std::make_shared<MultiArcPropagatorSettings<> >( propagatorSettingsList );
    multiArcPropagatorSettings->getOutputSettings( )->setUpdateDependentVariableInterpolator( true );

    // Define single-arc dynamics simulator
    MultiArcDynamicsSimulator< > simulator = MultiArcDynamicsSimulator< >(
            bodies, multiArcPropagatorSettings );


    // Retrieve dependent variables history.
    std::vector< std::map< double, Eigen::VectorXd > > dependentVariablesHistory = simulator.getDependentVariableHistory( );

    // Retrieve dependent variables interface.
    std::shared_ptr< MultiArcDependentVariablesInterface< > > dependentVariablesInterface =
            std::dynamic_pointer_cast< MultiArcDependentVariablesInterface< > >( simulator.getDependentVariablesInterface( ) );

    // Create dependent variables interpolator.
    for ( unsigned int i = 0 ; i < nbArcs ; i++ )
    {
        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > dependentVariablesInterpolator
                = std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd > >(
                        utilities::createVectorFromMapKeys< Eigen::VectorXd, double >( dependentVariablesHistory.at( i ) ),
                        utilities::createVectorFromMapValues< Eigen::VectorXd, double >( dependentVariablesHistory.at( i ) ), 4 );

        std::vector< double > testEpochs;
        testEpochs.push_back( arcStartTimes.at( i ) + 10.0 );
        testEpochs.push_back( arcStartTimes.at( i ) + ( arcEndTimes.at( i ) -  arcStartTimes.at( i ) ) / 4.0 );
        testEpochs.push_back( arcStartTimes.at( i ) + ( arcEndTimes.at( i ) -  arcStartTimes.at( i ) ) / 2.0 );
        testEpochs.push_back( arcStartTimes.at( i ) + 3.0 * ( arcEndTimes.at( i ) -  arcStartTimes.at( i ) ) / 4.0 );
        testEpochs.push_back( arcEndTimes.at( i ) - 10.0 );

        // Check consistency between interpolator results and interface results.
        for ( unsigned int j = 0 ; j < testEpochs.size( ) ; j++ )
        {
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( dependentVariablesInterpolator->interpolate( testEpochs[ j ] ),
                                               dependentVariablesInterface->getDependentVariables( testEpochs[ j ] ),
                                               std::numeric_limits< double >::epsilon( ) );
        }


        std::map< double, Eigen::VectorXd > totalAccelerationHistory;
        for ( auto itr : dependentVariablesHistory.at( i ) )
        {
            totalAccelerationHistory[ itr.first ] = itr.second.segment( 0, 3 );
        }

        // Create total acceleration history interpolator.
        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > totalAccelerationInterpolator
                = std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd > >(
                        utilities::createVectorFromMapKeys< Eigen::VectorXd, double >( totalAccelerationHistory ),
                        utilities::createVectorFromMapValues< Eigen::VectorXd, double >( totalAccelerationHistory ), 4 );

        // Total acceleration dependent variable settings.
        std::shared_ptr< SingleDependentVariableSaveSettings > totalAccelerationDependentVariable
                = std::make_shared< SingleDependentVariableSaveSettings >( total_acceleration_dependent_variable, "AlienSpaceship" );

        // Check consistency between interpolator results and interface results, for a single dependent variable.
        for ( unsigned int j = 0 ; j < testEpochs.size( ) ; j++ )
        {
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( totalAccelerationInterpolator->interpolate( testEpochs[ j ] ),
                                               dependentVariablesInterface->getSingleDependentVariable( totalAccelerationDependentVariable, testEpochs[ j ] ),
                                               std::numeric_limits< double >::epsilon( ) );
        }
    }
}

BOOST_AUTO_TEST_CASE( testHybridArcDependentVariablesInterface )
{

    // Set simulation time settings.
    const double simulationStartEpoch = 0.0;
    const double simulationEndEpoch = 3.0 * tudat::physical_constants::JULIAN_DAY;

    // Load spice kernel.
    spice_interface::loadStandardSpiceKernels( );

    // Define body settings for simulation.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Mars" );

    // Create system of bodies
    BodyListSettings bodySettings = getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 300.0, simulationEndEpoch + 300.0 );
    std::shared_ptr< BodySettings > phobosSettings = std::make_shared< BodySettings >( );
    phobosSettings->ephemerisSettings = getDefaultEphemerisSettings( "Phobos" );
    phobosSettings->gravityFieldSettings = getDefaultGravityFieldSettings( "Phobos", simulationStartEpoch - 300.0, simulationEndEpoch + 300.0 );
    bodySettings.addSettings( phobosSettings, "Phobos" );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // Create spacecraft object.
    bodies.createEmptyBody( "AlienSpaceship" );
    bodies.at( "AlienSpaceship" )->setConstantBodyMass( 400.0 );
    bodies.at( "AlienSpaceship" )->setEphemeris( std::make_shared< MultiArcEphemeris >(
            std::map< double, std::shared_ptr< Ephemeris > >( ), "Mars", "ECLIPJ2000" ) );

    // Create radiation pressure settings
    double referenceAreaRadiation = 4.0;
    double radiationPressureCoefficient = 1.2;
    std::vector< std::string > occultingBodies;
    occultingBodies.push_back( "Mars" );
    std::shared_ptr< RadiationPressureInterfaceSettings > radiationPressureSettings = std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
            "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

    // Create and set radiation pressure settings
    bodies.at( "AlienSpaceship" )->setRadiationPressureInterface(
            "Sun", createRadiationPressureInterface( radiationPressureSettings, "AlienSpaceship", bodies ) );

    // Define propagator settings variables.
    SelectedAccelerationMap phobosAccelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Define propagation settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > phobosAccelerations;
    phobosAccelerations[ "Mars" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );

    phobosAccelerationMap[ "Phobos" ] = phobosAccelerations;
    bodiesToPropagate.push_back( "Phobos" );
    centralBodies.push_back( "Mars" );

    basic_astrodynamics::AccelerationMap phobosAccelerationModelMap = createAccelerationModelsMap(
            bodies, phobosAccelerationMap, bodiesToPropagate, centralBodies );

    // Define propagator settings variables.
    SelectedAccelerationMap spacecraftAccelerationMap;
    std::vector< std::string > multiArcBodiesToPropagate;
    std::vector< std::string > multiArcCentralBodies;

    // Define propagation settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > spacecraftAccelerations;
    spacecraftAccelerations[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 5, 5 ) );
    spacecraftAccelerations[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );
    spacecraftAccelerations[ "Mars" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );
    spacecraftAccelerations[ "Phobos" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );
    spacecraftAccelerations[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::cannon_ball_radiation_pressure ) );

    spacecraftAccelerationMap[ "AlienSpaceship" ] = spacecraftAccelerations;
    multiArcBodiesToPropagate.push_back( "AlienSpaceship" );
    multiArcCentralBodies.push_back( "Mars" );

    basic_astrodynamics::AccelerationMap spacecraftAccelerationModelMap = createAccelerationModelsMap(
            bodies, spacecraftAccelerationMap, multiArcBodiesToPropagate, multiArcCentralBodies );

    // Set initial state for Phobos
    Eigen::Vector6d initialStatePhobos = bodies.at( "Phobos" )->getEphemeris( )->getCartesianState( simulationStartEpoch );

    // Set Keplerian elements for spacecraft.
    Eigen::Vector6d initialKeplerianElements;
    initialKeplerianElements( semiMajorAxisIndex ) = 1000.0E3;
    initialKeplerianElements( eccentricityIndex ) = 0.0;
    initialKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 26.04 );
    initialKeplerianElements( argumentOfPeriapsisIndex ) = unit_conversions::convertDegreesToRadians( 0.0 );
    initialKeplerianElements( longitudeOfAscendingNodeIndex ) = unit_conversions::convertDegreesToRadians( 0.0 );
    initialKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 0.0 );

    std::vector< Eigen::Vector6d > arcWiseKeplerianStates;
    arcWiseKeplerianStates.push_back( initialKeplerianElements );
    initialKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 60.0 );
    arcWiseKeplerianStates.push_back( initialKeplerianElements );
    initialKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 120.0 );
    arcWiseKeplerianStates.push_back( initialKeplerianElements );

    unsigned int nbArcs = arcWiseKeplerianStates.size( );

    std::vector< double > arcStartTimes, arcEndTimes;
    std::vector< Eigen::Vector6d > arcWiseStates;
    double marsGravitationalParameter = bodies.at( "Mars" )->getGravityFieldModel( )->getGravitationalParameter( );
    for ( unsigned int i = 0 ; i < nbArcs ; i++ )
    {
        arcWiseStates.push_back( convertKeplerianToCartesianElements( arcWiseKeplerianStates.at( i ), marsGravitationalParameter ) );
        arcStartTimes.push_back( simulationStartEpoch + i * physical_constants::JULIAN_DAY );
        arcEndTimes.push_back( simulationStartEpoch + ( i + 1 ) * physical_constants::JULIAN_DAY - 3600.0 );
    }

    // Define list of dependent variables to save.
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > multiArcDependentVariables;
    multiArcDependentVariables.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
            total_acceleration_dependent_variable, "AlienSpaceship" ) );
    multiArcDependentVariables.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
            keplerian_state_dependent_variable, "AlienSpaceship", "Mars" ) );

    const double fixedStepSize = 10.0;
    std::shared_ptr< IntegratorSettings< > > integratorSettings = rungeKutta4Settings( fixedStepSize );

    std::vector< std::shared_ptr< SingleArcPropagatorSettings< > > > propagatorSettingsList;
    for ( unsigned int i = 0 ; i < nbArcs ; i++ )
    {
        propagatorSettingsList.push_back( std::make_shared< TranslationalStatePropagatorSettings< double > >(
                multiArcCentralBodies, spacecraftAccelerationModelMap, multiArcBodiesToPropagate,
                arcWiseStates.at( i ), arcStartTimes.at( i ), integratorSettings,
                std::make_shared< propagators::PropagationTimeTerminationSettings >( arcEndTimes.at( i ) ), cowell, multiArcDependentVariables ) );
    }

    std::shared_ptr< MultiArcPropagatorSettings< > > multiArcPropagatorSettings =
            std::make_shared<MultiArcPropagatorSettings<> >( propagatorSettingsList );

    std::shared_ptr< TranslationalStatePropagatorSettings< double > > singleArcPropagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< double > >(
            centralBodies, phobosAccelerationModelMap, bodiesToPropagate, initialStatePhobos, simulationStartEpoch, integratorSettings,
            std::make_shared< propagators::PropagationTimeTerminationSettings >( simulationEndEpoch ), cowell );

    std::shared_ptr< HybridArcPropagatorSettings< > > hybridArcPropagatorSettings = std::make_shared< HybridArcPropagatorSettings< > >(
            singleArcPropagatorSettings, multiArcPropagatorSettings );
    hybridArcPropagatorSettings->getOutputSettings( )->setUpdateDependentVariableInterpolator( true );


    // Define single-arc dynamics simulator
    HybridArcDynamicsSimulator< > simulator = HybridArcDynamicsSimulator< >( bodies, hybridArcPropagatorSettings );


    // Retrieve dependent variables history.
    std::map< double, Eigen::VectorXd > singleArcDependentVariablesHistory = simulator.getSingleArcDynamicsSimulator( )->getDependentVariableHistory( );
    std::vector< std::map< double, Eigen::VectorXd > > multiArcDependentVariablesHistory = simulator.getMultiArcDynamicsSimulator( )->getDependentVariableHistory( );

    // Retrieve dependent variables interface.
    std::shared_ptr< HybridArcDependentVariablesInterface< > > dependentVariablesInterface =
            std::dynamic_pointer_cast< HybridArcDependentVariablesInterface< > >( simulator.getDependentVariablesInterface( ) );

    std::cout<<( dependentVariablesInterface == nullptr )<<std::endl;
    // Create dependent variables interpolator.
    for ( unsigned int i = 0 ; i < nbArcs ; i++ )
    {
        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > dependentVariablesInterpolator
                = std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd > >(
                        utilities::createVectorFromMapKeys< Eigen::VectorXd, double >( multiArcDependentVariablesHistory.at( i ) ),
                        utilities::createVectorFromMapValues< Eigen::VectorXd, double >( multiArcDependentVariablesHistory.at( i ) ), 4 );

        std::vector< double > testEpochs;
        testEpochs.push_back( arcStartTimes.at( i ) + 10.0 );
        testEpochs.push_back( arcStartTimes.at( i ) + ( arcEndTimes.at( i ) -  arcStartTimes.at( i ) ) / 4.0 );
        testEpochs.push_back( arcStartTimes.at( i ) + ( arcEndTimes.at( i ) -  arcStartTimes.at( i ) ) / 2.0 );
        testEpochs.push_back( arcStartTimes.at( i ) + 3.0 * ( arcEndTimes.at( i ) -  arcStartTimes.at( i ) ) / 4.0 );
        testEpochs.push_back( arcEndTimes.at( i ) - 10.0 );

        bool exceptionCaught = false;
        try
        {
            dependentVariablesInterface->getDependentVariables( testEpochs[ 0 ] );
        }
        catch( std::runtime_error const& )
        {
            exceptionCaught = true;
        }
        BOOST_CHECK_EQUAL( exceptionCaught, true );


        std::map< double, Eigen::VectorXd > totalAccelerationHistory;
        for ( auto itr : multiArcDependentVariablesHistory.at( i ) )
        {
            totalAccelerationHistory[ itr.first ] = itr.second.segment( 0, 3 );
        }

        // Create total acceleration history interpolator.
        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > totalAccelerationInterpolator
                = std::make_shared< interpolators::LagrangeInterpolator< double, Eigen::VectorXd > >(
                        utilities::createVectorFromMapKeys< Eigen::VectorXd, double >( totalAccelerationHistory ),
                        utilities::createVectorFromMapValues< Eigen::VectorXd, double >( totalAccelerationHistory ), 4 );

//        // Total acceleration dependent variable settings.
//        std::shared_ptr< SingleDependentVariableSaveSettings > totalAccelerationDependentVariable
//                = std::make_shared< SingleDependentVariableSaveSettings >( total_acceleration_dependent_variable, "AlienSpaceship" );
//
//
//        // Check consistency between interpolator results and interface results, for a single dependent variable.
//        for ( unsigned int j = 0 ; j < testEpochs.size( ) ; j++ )
//        {
//            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( totalAccelerationInterpolator->interpolate( testEpochs[ j ] ),
//                                               dependentVariablesInterface->getSingleDependentVariable( totalAccelerationDependentVariable, testEpochs[ j ] ),
//                                               ( 10.0 * std::numeric_limits< double >::epsilon( ) ) );
//        }
    }

}

BOOST_AUTO_TEST_SUITE_END( )

}

}
