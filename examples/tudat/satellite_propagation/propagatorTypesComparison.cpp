/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References:
 *      Mooij, E. "Orbit-State Model Selection for Solar-Sailing Mission Optimization."
 *          AIAA/AAS astro Specialist Conference. 2012.
 */

#include <tudat/simulation/simulation.h>
#include "tudat/io/applicationOutput.h"

#include "tudat/math/statistics/basicStatistics.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/basics/utilities.h"

//! Execute propagation of orbit of Satellite around the Earth.
int main( )
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    using namespace tudat;
    using namespace tudat::aerodynamics;
    using namespace tudat::basic_astrodynamics;
    using namespace tudat::basic_mathematics;
    using namespace tudat::ephemerides;
    using namespace tudat::estimatable_parameters;
    using namespace tudat::gravitation;
    using namespace tudat::input_output;
    using namespace tudat::numerical_integrators;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::propagators;
    using namespace tudat::simulation_setup;
    using namespace tudat::unit_conversions;

    using namespace tudat_applications;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            DEFINE TEST CASES             //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Predefine variables
    bool keplerOrbit = false; // toggle use of accelerations (true == no accelerations)
    double simulationDuration = 10.0 * physical_constants::JULIAN_DAY;
    double integrationRelativeTolerance = 1.0e-12;
    double integrationAbsoluteTolerance = 1.0e-12;
    double integrationReferenceTolerance = 1.0e-15;
    double integrationConstantTimeStepSize = 5.0;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Load Spice kernels
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation time settings
    const double simulationStartEpoch = 7.0 * physical_constants::JULIAN_YEAR + 30.0 * 6.0 * physical_constants::JULIAN_DAY;
    const double simulationEndEpoch = simulationDuration + simulationStartEpoch;

    // Define body settings for simulation
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Moon" );

    // Create body objects
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 1.0e3, simulationEndEpoch + 1.0e3 );
    for ( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
    {
        bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
        bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    }
    bodySettings[ "Earth" ]->gravityFieldSettings = std::make_shared< FromFileSphericalHarmonicsGravityFieldSettings >( ggm02s );
    bodySettings[ "Earth" ]->atmosphereSettings = std::make_shared< ExponentialAtmosphereSettings >( aerodynamics::earth );
    SystemOfBodies bodies = createBodies( bodySettings );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create spacecraft object
    bodies[ "Satellite" ] = std::make_shared< Body >( );
    const double satelliteMass = 1000.0;
    bodies[ "Satellite" ]->setConstantBodyMass( satelliteMass );

    // Set constant aerodynamic drag coefficient
    const double referenceAreaAerodynamic = 37.5;
    const Eigen::Vector3d aerodynamicCoefficients = 2.2 * Eigen::Vector3d::UnitX( ); // only drag coefficient
    std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
            std::make_shared< ConstantAerodynamicCoefficientSettings >( referenceAreaAerodynamic, aerodynamicCoefficients, true, true );

    bodies[ "Satellite" ]->setAerodynamicCoefficientInterface(
                createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Satellite" ) );

    // Create radiation pressure settings
    double referenceAreaRadiation = referenceAreaAerodynamic;
    double radiationPressureCoefficient = 1.25;
    std::vector< std::string > occultingBodies;
    occultingBodies.push_back( "Earth" );
    std::shared_ptr< RadiationPressureInterfaceSettings > SatelliteRadiationPressureSettings =
            std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

    // Create and set radiation pressure settings
    bodies[ "Satellite" ]->setRadiationPressureInterface( "Sun", createRadiationPressureInterface(
                                                               SatelliteRadiationPressureSettings, "Satellite", bodies ) );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodies, "SSB", "J2000" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define propagator settings variables
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Switch between Kepler orbit or perturbed environment
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfSatellite;
    if ( keplerOrbit )
    {
        // Only central gravity
        accelerationsOfSatellite[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
    }
    else
    {
        // Define spherical harmonics, third bodies, solar radiation and aerodynamic forces
        accelerationsOfSatellite[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 4, 4 ) );
        for ( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
        {
            if ( bodiesToCreate.at( i ) != "Earth" )
            {
                accelerationsOfSatellite[ bodiesToCreate.at( i ) ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
            }
        }
        accelerationsOfSatellite[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( cannon_ball_radiation_pressure ) );
        accelerationsOfSatellite[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( aerodynamic ) );
    }

    // Add acceleration information
    accelerationMap[ "Satellite" ] = accelerationsOfSatellite;
    bodiesToPropagate.push_back( "Satellite" );
    centralBodies.push_back( "Earth" );

    AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodies, accelerationMap, bodiesToPropagate, centralBodies );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE INITIAL CONDITIONS              ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set Keplerian initial conditions
    Eigen::Vector6d satelliteInitialStateInKeplerianElements;
    satelliteInitialStateInKeplerianElements( semiMajorAxisIndex ) = 6778136.0;
    satelliteInitialStateInKeplerianElements( eccentricityIndex ) = 0.003;
    satelliteInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 93.0 );
    satelliteInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 158.7 );
    satelliteInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 23.4 );
    satelliteInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 0.0 );

    // Convert to Cartesian elements
    double mainGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
    const Eigen::Vector6d satelliteInitialState = convertKeplerianToCartesianElements(
                satelliteInitialStateInKeplerianElements, mainGravitationalParameter );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             LOOP OVER PROPAGATORS                  ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Vector of function evaluations and times
    std::vector< unsigned int > numberOfFunctionEvaluations;
    std::vector< double > totalPropagationTime;

    // Loop over propagators
    std::vector< string > nameAdditionPropagator = { "_cowell", "_encke", "_kepl", "_equi", "_usm7", "_usm6", "_usmem", "_ref" };
    std::vector< string > nameAdditionIntegrator = { "_var", "_const" };
    for ( unsigned int propagatorType = 0; propagatorType < 8; propagatorType++ )
    {
        // Progress
        std::cout << std::endl << "Propagator: " << propagatorType + 1 << std::endl;

        // Loop over integrators
        for ( unsigned int integratorType = 0; integratorType < 2; integratorType++ )
        {
            // Progress
            std::cout << "Integrator: " << integratorType + 1 << std::endl;

            ///////////////////////     CREATE SIMULATION SETTINGS          ////////////////////////////////////////////

            // Propagator settings
            std::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettings;
            if ( propagatorType == 7 )
            {
                // Reference trajectory
                propagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< > >(
                            centralBodies, accelerationModelMap, bodiesToPropagate, satelliteInitialState,
                            simulationEndEpoch, cowell );
            }
            else
            {
                // Propagator dependent on loop
                propagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< > >(
                            centralBodies, accelerationModelMap, bodiesToPropagate, satelliteInitialState,
                            simulationEndEpoch, static_cast< TranslationalPropagatorType >( propagatorType ) );
            }

            // Integrator settings
            std::shared_ptr< IntegratorSettings< > > integratorSettings;
            if ( propagatorType == 7 )
            {
                // Reference trajectory
                integratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettingsScalarTolerances< > >(
                            simulationStartEpoch, 100.0, rungeKuttaFehlberg78, 1.0e-5, 1.0e5,
                            integrationReferenceTolerance, integrationReferenceTolerance );
            }
            else
            {
                // Integrator dependent on loop
                if ( integratorType == 0 )
                {
                    integratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettingsScalarTolerances< > >(
                                simulationStartEpoch, 100.0, rungeKuttaFehlberg56, 1.0e-5, 1.0e5,
                                integrationRelativeTolerance, integrationAbsoluteTolerance );
                }
                else if ( integratorType == 1 )
                {
                    integratorSettings = std::make_shared< IntegratorSettings< > >(
                                rungeKutta4, simulationStartEpoch, integrationConstantTimeStepSize );
                }
            }

            ///////////////////////     PROPAGATE ORBIT                     ////////////////////////////////////////////

            // Simulate orbit and output computation time
            SingleArcDynamicsSimulator< > dynamicsSimulator( bodies, integratorSettings, propagatorSettings,
                                                             true, false, false, false );

            // Retrieve results
            std::map< double, Eigen::VectorXd > cartesianIntegrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
            if ( propagatorType != 7 )
            {
                if ( integratorType == 0 )
                {
                    // Store number of function evaluations for variable step-size
                    numberOfFunctionEvaluations.push_back( dynamicsSimulator.getCumulativeNumberOfFunctionEvaluations( ).rbegin( )->second );
                    std::cout << "Total Number of Function Evaluations: " << numberOfFunctionEvaluations.back( ) << std::endl;
                }
                else
                {
                    // Store propagation time for constant step-size
                    totalPropagationTime.push_back( dynamicsSimulator.getCumulativeComputationTimeHistory( ).rbegin( )->second );
                    std::cout << "Total propagation Time: " << totalPropagationTime.back( ) << std::endl;
                }
            }
            else
            {
                std::cout << "Total Number of Function Evaluations: " <<
                          dynamicsSimulator.getCumulativeNumberOfFunctionEvaluations( ).rbegin( )->second << std::endl;
                std::cout << "Total propagation Time: " <<
                          dynamicsSimulator.getCumulativeComputationTimeHistory( ).rbegin( )->second << std::endl;
            }

            ///////////////////////     PROVIDE OUTPUT TO FILES             ////////////////////////////////////////////

            // Write perturbed satellite propagation history to file
            writeDataMapToTextFile( cartesianIntegrationResult, "cartesian" + nameAdditionPropagator[ propagatorType ] +
                                    nameAdditionIntegrator[ integratorType ] + ".dat", getOutputPath( "PropagatorTypesComparison/" ) );

            // Break loop if reference propagator
            if ( propagatorType == 7 )
            {
                break;
            }
        }
    }

    // Write function evaluations and times to file
    writeMatrixToFile( utilities::convertStlVectorToEigenVector( numberOfFunctionEvaluations ),
                       "functionEvaluations.dat", 16, getOutputPath( "PropagatorTypesComparison/" ) );
    writeMatrixToFile( utilities::convertStlVectorToEigenVector( totalPropagationTime ),
                       "propagationTime.dat", 16, getOutputPath( "PropagatorTypesComparison/" ) );

    // Final statement
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed
    return EXIT_SUCCESS;
}
