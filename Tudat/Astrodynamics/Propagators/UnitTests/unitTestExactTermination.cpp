/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <limits>
#include <string>
#include <thread>

#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"

#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/Mathematics/NumericalIntegrators/rungeKuttaCoefficients.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/SimulationSetup/PropagationSetup/dynamicsSimulator.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createBodies.h"
#include "Tudat/SimulationSetup/PropagationSetup/createNumericalSimulator.h"
#include "Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_exact_termination )

// Test Encke propagator for point mass, and spherical harmonics central body.
BOOST_AUTO_TEST_CASE( testEnckePopagatorForSphericalHarmonicCentralBodies )
{
    for( unsigned int integratorCase = 1; integratorCase < 2; integratorCase++ )
    {
        for( unsigned int simulationCase = 3; simulationCase < 4; simulationCase++ )
        {
            for( unsigned int direction = 0; direction < 1; direction++ )
            {
                std::cout<<integratorCase<<" "<<direction<<" "<<simulationCase<<std::endl;
                using namespace tudat;
                using namespace simulation_setup;
                using namespace propagators;
                using namespace numerical_integrators;
                using namespace orbital_element_conversions;
                using namespace basic_mathematics;
                using namespace gravitation;
                using namespace numerical_integrators;

                // Load Spice kernels.
                spice_interface::loadStandardSpiceKernels( );

                // Set simulation time settings.
                double simulationStartEpoch;
                double simulationEndEpoch;

                double directionMultiplier = 1.0;
                if( direction == 0 )
                {
                    simulationStartEpoch = 0.0;
                    simulationEndEpoch = 0.2 * tudat::physical_constants::JULIAN_DAY;
                }
                else
                {
                    simulationStartEpoch = 0.2 * tudat::physical_constants::JULIAN_DAY;
                    simulationEndEpoch = 0.0;
                    directionMultiplier = -1.0;
                }


                // Define body settings for simulation.
                std::vector< std::string > bodiesToCreate;
                bodiesToCreate.push_back( "Sun" );
                bodiesToCreate.push_back( "Earth" );
                bodiesToCreate.push_back( "Moon" );

                // Create body objects.
                std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings;
                if( direction == 0 )
                {
                    bodySettings =
                            getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 300.0, simulationEndEpoch + 300.0 );
                }
                else
                {
                    bodySettings =
                            getDefaultBodySettings( bodiesToCreate, simulationEndEpoch - 300.0, simulationStartEpoch + 300.0 );
                }
                NamedBodyMap bodyMap = createBodies( bodySettings );

                // Create spacecraft object.
                bodyMap[ "Vehicle" ] = boost::make_shared< simulation_setup::Body >( );
                bodyMap[ "Vehicle" ]->setConstantBodyMass( 400.0 );
                bodyMap[ "Vehicle" ]->setEphemeris( boost::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                                                        boost::shared_ptr< interpolators::OneDimensionalInterpolator
                                                        < double, Eigen::Vector6d  > >( ), "Earth", "ECLIPJ2000" ) );


                // Finalize body creation.
                setGlobalFrameBodyEphemerides( bodyMap, "Earth", "ECLIPJ2000" );

                // Define propagator settings variables.
                SelectedAccelerationMap accelerationMap;
                std::vector< std::string > bodiesToPropagate;
                std::vector< std::string > centralBodies;

                // Define propagation settings.
                std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;

                {
                    accelerationsOfVehicle[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >(
                                                                     basic_astrodynamics::central_gravity ) );
                    accelerationsOfVehicle[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >(
                                                                   basic_astrodynamics::central_gravity ) );
                    accelerationsOfVehicle[ "Moon" ].push_back( boost::make_shared< AccelerationSettings >(
                                                                    basic_astrodynamics::central_gravity ) );
                }

                accelerationMap[  "Vehicle" ] = accelerationsOfVehicle;
                bodiesToPropagate.push_back( "Vehicle" );
                centralBodies.push_back( "Earth" );
                basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                            bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

                // Set Keplerian elements for Vehicle.
                Eigen::Vector6d vehicleInitialStateInKeplerianElements;
                vehicleInitialStateInKeplerianElements( semiMajorAxisIndex ) = 8000.0E3;
                vehicleInitialStateInKeplerianElements( eccentricityIndex ) = 0.1;
                vehicleInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
                vehicleInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
                        = unit_conversions::convertDegreesToRadians( 235.7 );
                vehicleInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
                        = unit_conversions::convertDegreesToRadians( 23.4 );
                vehicleInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );

                double earthGravitationalParameter = bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
                const Eigen::Vector6d vehicleInitialState = convertKeplerianToCartesianElements(
                            vehicleInitialStateInKeplerianElements, earthGravitationalParameter );

                // Define propagator settings (Cowell)
                boost::shared_ptr< PropagationTerminationSettings > terminationSettings;
                std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
                dependentVariables.push_back(
                            boost::make_shared< SingleDependentVariableSaveSettings >( relative_distance_dependent_variable,
                                                                                       "Vehicle", "Earth" ) );
                double finalTestTime;
                double secondFinalTestTime;

                if( direction == 0 )
                {
                    finalTestTime = 322.5;
                    secondFinalTestTime = 501.0;
                }
                else
                {
                    finalTestTime = 11737.5;
                    secondFinalTestTime = 11701.0;
                }
                if( simulationCase == 0 )
                {
                    terminationSettings = boost::make_shared< PropagationTimeTerminationSettings >(
                                simulationEndEpoch - directionMultiplier * 4.5, true );
                }
                else if( simulationCase == 1 )
                {
                    terminationSettings = boost::make_shared< PropagationDependentVariableTerminationSettings >(
                                dependentVariables.at( 0 ), 8.7E6, false, true,
                                boost::make_shared< root_finders::RootFinderSettings >(
                                    root_finders::bisection_root_finder, 1.0E-6, 100 ) );
                }
                else if( simulationCase == 2 )
                {
                    std::vector< boost::shared_ptr< PropagationTerminationSettings > > terminationSettingsList;
                    terminationSettingsList.push_back(
                                boost::make_shared< PropagationTimeTerminationSettings >( finalTestTime, true ) );
                    terminationSettingsList.push_back(
                                boost::make_shared< PropagationDependentVariableTerminationSettings >(
                                    dependentVariables.at( 0 ), 8.7E6, false, true,
                                    boost::make_shared< root_finders::RootFinderSettings >(
                                        root_finders::bisection_root_finder, 1.0E-6, 100 ) ) );
                    terminationSettings = boost::make_shared< PropagationHybridTerminationSettings >(
                                terminationSettingsList, true );
                }
                else if( simulationCase == 3 )
                {
                    std::vector< boost::shared_ptr< PropagationTerminationSettings > > terminationSettingsList;
                    terminationSettingsList.push_back(
                                boost::make_shared< PropagationTimeTerminationSettings >( finalTestTime, true ) );
                    terminationSettingsList.push_back(
                                boost::make_shared< PropagationDependentVariableTerminationSettings >(
                                    dependentVariables.at( 0 ), 8.7E6, false, true,
                                    boost::make_shared< root_finders::RootFinderSettings >(
                                        root_finders::bisection_root_finder, 1.0E-6, 100 ) ) );
                    terminationSettings = boost::make_shared< PropagationHybridTerminationSettings >(
                                terminationSettingsList, false );
                }
                else if( simulationCase == 4 )
                {
                    std::vector< boost::shared_ptr< PropagationTerminationSettings > > terminationSettingsList;
                    terminationSettingsList.push_back(
                                boost::make_shared< PropagationTimeTerminationSettings >( secondFinalTestTime, true ) );
                    terminationSettingsList.push_back(
                                boost::make_shared< PropagationDependentVariableTerminationSettings >(
                                    dependentVariables.at( 0 ), 8.7E6, false, true,
                                    boost::make_shared< root_finders::RootFinderSettings >(
                                        root_finders::bisection_root_finder, 1.0E-6, 100 ) ) );
                    terminationSettings = boost::make_shared< PropagationHybridTerminationSettings >(
                                terminationSettingsList, true );
                }


                boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                        boost::make_shared< TranslationalStatePropagatorSettings< double > >
                        ( centralBodies, accelerationModelMap, bodiesToPropagate, vehicleInitialState, terminationSettings, cowell,
                          boost::make_shared< DependentVariableSaveSettings >( dependentVariables ) );

                // Define integrator settings.
                const double fixedStepSize = 5.0;
                boost::shared_ptr< IntegratorSettings< > > integratorSettings;
                if( integratorCase == 0 )
                {
                    integratorSettings = boost::make_shared< IntegratorSettings< > >
                            ( rungeKutta4, simulationStartEpoch, directionMultiplier * fixedStepSize );

                }
                else
                {
                    integratorSettings = boost::make_shared< RungeKuttaVariableStepSizeSettings< double > >
                            ( rungeKuttaVariableStepSize, simulationStartEpoch, directionMultiplier * fixedStepSize,
                              RungeKuttaCoefficients::CoefficientSets::rungeKuttaFehlberg45,
                              1.0E-3, 1.0E3, 1.0E-12, 1.0E-12 );
                }

                // Propagate orbit with Cowell method
                SingleArcDynamicsSimulator< double > dynamicsSimulator(
                            bodyMap, integratorSettings, propagatorSettings, true, false, false );
                std::map< double, Eigen::VectorXd > stateHistory = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
                std::map< double, Eigen::VectorXd > dependentVariableHistory = dynamicsSimulator.getDependentVariableHistory( );

                if( simulationCase == 0 )
                {
                    if( direction == 0 )
                    {
                        BOOST_CHECK_SMALL( std::fabs( stateHistory.rbegin( )->first -
                                                      ( simulationEndEpoch - 4.5 ) ), 1.0E-10 );
                    }
                    else
                    {
                        BOOST_CHECK_SMALL( std::fabs( stateHistory.begin( )->first -
                                                      ( simulationEndEpoch + 4.5 ) ), 1.0E-10 );
                    }
                }
                else if( simulationCase == 1 )
                {
                    if( direction == 0 )
                    {
                        BOOST_CHECK_SMALL( std::fabs( dependentVariableHistory.rbegin( )->second( 0 ) - 8.7E6 ), 0.01 );
                    }
                    else
                    {
                        BOOST_CHECK_SMALL( std::fabs( dependentVariableHistory.begin( )->second( 0 ) - 8.7E6 ), 0.01 );
                    }
                }
                else if( simulationCase == 2 )
                {
                    if( direction == 0 )
                    {
                        BOOST_CHECK_SMALL( std::fabs( stateHistory.rbegin( )->first - finalTestTime ), 0.01 );
                    }
                    else
                    {
                        BOOST_CHECK_SMALL( std::fabs( dependentVariableHistory.begin( )->second( 0 ) - 8.7E6 ), 0.01 );
                    }
                }
                else if( simulationCase == 3 )
                {
                    if( direction == 0 )
                    {
                        BOOST_CHECK_SMALL( std::fabs( dependentVariableHistory.rbegin( )->second( 0 ) - 8.7E6 ), 0.01 );
                    }
                    else
                    {
                        BOOST_CHECK_SMALL( std::fabs( stateHistory.begin( )->first - finalTestTime ), 0.01 );
                    }
                }
                else if( simulationCase == 4 )
                {
                    if( direction == 0 )
                    {
                        BOOST_CHECK_SMALL( std::fabs( stateHistory.rbegin( )->first - secondFinalTestTime ), 0.01 );
                    }
                    else
                    {
                        BOOST_CHECK_SMALL( std::fabs( stateHistory.begin( )->first - secondFinalTestTime ), 0.01 );
                    }
                }

                //                std::cout<<"First time/altitude: "<<std::setprecision( 12 )<<dependentVariableHistory.begin( )->first<<" "<<
                //                           dependentVariableHistory.begin( )->second<<std::endl;
                //                std::cout<<"Final time/altitude: "<<std::setprecision( 12 )<<dependentVariableHistory.rbegin( )->first<<" "<<
                //                           dependentVariableHistory.rbegin( )->second<<std::endl;

                //                for( std::map< double, Eigen::VectorXd >::const_iterator it = dependentVariableHistory.begin( );
                //                     it != dependentVariableHistory.end( ); it++ )
                //                {
                //                    std::cout<<"Time/altitude: "<<it->first<<" "<<it->second<<std::endl;
                //                }
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}

}


