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




#include <memory>

#include "tudat/astro/aerodynamics/testApolloCapsuleCoefficients.h"
#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/astro/basic_astro/geodeticCoordinateConversions.h"
#include "tudat/astro/gravitation/gravityFieldVariations.h"
#include "tudat/basics/testMacros.h"
#include "tudat/simulation/estimation.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/basic_astro/sphericalStateConversions.h"
#include "tudat/astro/basic_astro/stateVectorIndices.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/estimation_setup/createNumericalSimulator.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/io/basicInputOutput.h"
#include <limits>
#include <string>
#include "tudat/astro/basic_astro/celestialBodyConstants.h"

#include <Eigen/Core>

namespace tudat
{

namespace unit_tests
{

using namespace tudat;
using namespace ephemerides;
using namespace interpolators;
using namespace numerical_integrators;
using namespace spice_interface;
using namespace simulation_setup;
using namespace basic_astrodynamics;
using namespace orbital_element_conversions;
using namespace unit_conversions;
using namespace propagators;
using namespace aerodynamics;
using namespace basic_mathematics;
using namespace input_output;
using namespace estimatable_parameters;

BOOST_AUTO_TEST_SUITE( test_dependent_variable_output )

Eigen::VectorXd customDependentVariable1( const simulation_setup::SystemOfBodies& bodies )
{
    return bodies.at( "Apollo" )->getPosition( ) / 2.0;
}

Eigen::VectorXd customDependentVariable2( const simulation_setup::SystemOfBodies& bodies )
{
    return ( Eigen::VectorXd( 1 ) << bodies.at( "Apollo" )->getPosition( ).norm( ) / 2.0 ).finished( );
}

//! Propagate entry of Apollo capsule, and save a list of dependent variables during entry. The saved dependent variables
//! are compared against theoretical/manual values in this test.
BOOST_AUTO_TEST_CASE( testDependentVariableOutput )
{
    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation start epoch.
    const double simulationStartEpoch = 0.0;

    // Set simulation end epoch.
    const double simulationEndEpoch = 3300.0;

    // Set numerical integration fixed step size.
    const double fixedStepSize = 1.0;


    // Set Keplerian elements for Capsule.
    Eigen::Vector6d apolloInitialStateInKeplerianElements;
    apolloInitialStateInKeplerianElements( semiMajorAxisIndex ) = spice_interface::getAverageRadius( "Earth" ) + 120.0E3;
    apolloInitialStateInKeplerianElements( eccentricityIndex ) = 0.005;
    apolloInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
    apolloInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
            = unit_conversions::convertDegreesToRadians( 235.7 );
    apolloInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
            = unit_conversions::convertDegreesToRadians( 23.4 );
    apolloInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );

    // Convert apollo state from Keplerian elements to Cartesian elements.
    const double earthGravitationalParameter = getBodyGravitationalParameter( "Earth" );
    const Eigen::Vector6d apolloInitialState = convertKeplerianToCartesianElements(
                apolloInitialStateInKeplerianElements,
                getBodyGravitationalParameter( "Earth" ) );

#if TUDAT_BUILD_WITH_NRLMSISE
    unsigned int maximumTestCase = 3;
#else
    unsigned int maximumTestCase = 2;
#endif
    for( unsigned int testCase = 0; testCase < maximumTestCase; testCase++ )
    {
        std::map< double, Eigen::Vector3d > cowellAcceleration;
        for( unsigned int propagatorType = 0; propagatorType < 2; propagatorType++ )
        {
            // Define simulation body settings.
            BodyListSettings bodySettings =
                    getDefaultBodySettings( { "Earth", "Moon" }, simulationStartEpoch - 10.0 * fixedStepSize,
                                            simulationEndEpoch + 10.0 * fixedStepSize, "Earth", "ECLIPJ2000" );
            bodySettings.at( "Earth" )->gravityFieldSettings =
                    std::make_shared< simulation_setup::GravityFieldSettings >( central_spice );

            if( testCase >= 2 )
            {
                bodySettings.at( "Earth" )->atmosphereSettings =
                        std::make_shared< simulation_setup::AtmosphereSettings >( nrlmsise00 );
            }

            bool isOblateSpheroidUsed = 0;
            double oblateSpheroidEquatorialRadius = 6378.0E3;
            double oblateSpheroidFlattening = 1.0 / 300.0;

            if( testCase% 2 == 0 )
            {
                isOblateSpheroidUsed = 1;
                bodySettings.at( "Earth" )->shapeModelSettings =
                        std::make_shared< simulation_setup::OblateSphericalBodyShapeSettings >(
                            oblateSpheroidEquatorialRadius, oblateSpheroidFlattening );
            }

            // Create Earth object
            simulation_setup::SystemOfBodies bodies = simulation_setup::createSystemOfBodies( bodySettings );

            // Create vehicle objects.
            bodies.createEmptyBody( "Apollo" );

            // Create vehicle aerodynamic coefficients
            bodies.at( "Apollo" )->setAerodynamicCoefficientInterface(
                        unit_tests::getApolloCoefficientInterface( ) );
            bodies.at( "Apollo" )->setConstantBodyMass( 5.0E3 );
            bodies.at( "Apollo" )->setRotationalEphemeris(
                        createRotationModel(
                            std::make_shared< PitchTrimRotationSettings >( "Earth", "ECLIPJ2000", "VehicleFixed" ),
                            "Apollo", bodies ) );

            std::shared_ptr< system_models::VehicleSystems > vehicleSystems =
                    std::make_shared< system_models::VehicleSystems >( );

            double noseRadius = 0.7;
            double wallEmissivity = 0.7;
            vehicleSystems->setNoseRadius( noseRadius );
            vehicleSystems->setWallEmissivity( wallEmissivity );

            bodies.at( "Apollo" )->setVehicleSystems( vehicleSystems );

            // Define propagator settings variables.
            SelectedAccelerationMap accelerationMap;
            std::vector< std::string > bodiesToPropagate;
            std::vector< std::string > centralBodies;

            // Define acceleration model settings.
            std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfApollo;
            accelerationsOfApollo[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
            accelerationsOfApollo[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( aerodynamic ) );
            accelerationsOfApollo[ "Moon" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
            accelerationMap[ "Apollo" ] = accelerationsOfApollo;

            bodiesToPropagate.push_back( "Apollo" );
            centralBodies.push_back( "Earth" );

            // Set initial state
            Eigen::Vector6d systemInitialState = apolloInitialState;

            // Define list of dependent variables to save.
            std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
            dependentVariables.push_back(
                        std::make_shared< SingleDependentVariableSaveSettings >( mach_number_dependent_variable, "Apollo" ) );
            dependentVariables.push_back(
                        std::make_shared< SingleDependentVariableSaveSettings >( altitude_dependent_variable,
                                                                                 "Apollo", "Earth" ) );
            dependentVariables.push_back(
                        std::make_shared< SingleDependentVariableSaveSettings >( relative_distance_dependent_variable,
                                                                                 "Apollo", "Earth" ) );
            dependentVariables.push_back(
                        std::make_shared< SingleDependentVariableSaveSettings >( relative_speed_dependent_variable,
                                                                                 "Apollo", "Earth" ) );
            dependentVariables.push_back(
                        std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                            point_mass_gravity, "Apollo", "Earth", 1 ) );

            dependentVariables.push_back(
                        std::make_shared< SingleDependentVariableSaveSettings >( total_aerodynamic_g_load_variable,
                                                                                 "Apollo", "Earth" ) );
            dependentVariables.push_back(
                        std::make_shared< SingleDependentVariableSaveSettings >( stagnation_point_heat_flux_dependent_variable,
                                                                                 "Apollo", "Earth" ) );

            dependentVariables.push_back(
                        std::make_shared< SingleDependentVariableSaveSettings >( local_temperature_dependent_variable,
                                                                                 "Apollo", "Earth" ) );
            dependentVariables.push_back(
                        std::make_shared< SingleDependentVariableSaveSettings >( geodetic_latitude_dependent_variable,
                                                                                 "Apollo", "Earth" ) );

            dependentVariables.push_back(
                        std::make_shared< SingleDependentVariableSaveSettings >( local_density_dependent_variable,
                                                                                 "Apollo", "Earth" ) );

            dependentVariables.push_back(
                        std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                            "Apollo", reference_frames::latitude_angle ) );
            dependentVariables.push_back(
                        std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                            "Apollo", reference_frames::longitude_angle ) );

            dependentVariables.push_back(
                        std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                            "Apollo", reference_frames::angle_of_attack ) );
            dependentVariables.push_back(
                        std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                            "Apollo", reference_frames::angle_of_sideslip ) );
            dependentVariables.push_back(
                        std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                            "Apollo", reference_frames::bank_angle ) );




            dependentVariables.push_back(
                        std::make_shared< SingleDependentVariableSaveSettings >( relative_position_dependent_variable,
                                                                                 "Apollo", "Earth" ) );
            dependentVariables.push_back(
                        std::make_shared< SingleDependentVariableSaveSettings >( relative_velocity_dependent_variable,
                                                                                 "Apollo", "Earth" ) );
            dependentVariables.push_back(
                        std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                            point_mass_gravity, "Apollo", "Earth", 0 ) );
            dependentVariables.push_back(
                        std::make_shared< SingleDependentVariableSaveSettings >(
                            total_acceleration_dependent_variable, "Apollo" ) );
            dependentVariables.push_back(
                        std::make_shared< SingleDependentVariableSaveSettings >(
                            aerodynamic_moment_coefficients_dependent_variable, "Apollo" ) );
            dependentVariables.push_back(
                        std::make_shared< SingleDependentVariableSaveSettings >(
                            aerodynamic_force_coefficients_dependent_variable, "Apollo" ) );
            dependentVariables.push_back(
                        std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                            aerodynamic, "Apollo", "Earth", 0 ) );
            dependentVariables.push_back(
                        std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                            point_mass_gravity, "Apollo", "Moon", 0 ) );
            dependentVariables.push_back(
                        std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                            third_body_point_mass_gravity, "Apollo", "Moon", 0 ) );
            dependentVariables.push_back(
                        std::make_shared< SingleDependentVariableSaveSettings >(
                            keplerian_state_dependent_variable,  "Apollo", "Earth" ) );
            dependentVariables.push_back(
                        std::make_shared< SingleDependentVariableSaveSettings >(
                            modified_equinocial_state_dependent_variable,  "Apollo", "Earth" ) );
            dependentVariables.push_back(
                        std::make_shared< SingleDependentVariableSaveSettings >(
                            body_fixed_relative_cartesian_position,  "Apollo", "Earth" ) );


            // Create acceleration models and propagation settings.
            basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                        bodies, accelerationMap, bodiesToPropagate, centralBodies );

            std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings;
            if( propagatorType == 0 )
            {
                propagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< double > >
                        ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState,
                          std::make_shared< propagators::PropagationTimeTerminationSettings >( 3200.0 ), cowell,
                          dependentVariables );
            }
            else
            {
                propagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< double > >
                        ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState,
                          std::make_shared< propagators::PropagationTimeTerminationSettings >( 3200.0 ), gauss_modified_equinoctial,
                          dependentVariables );
            }

            std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesToAdd;
            dependentVariablesToAdd.push_back(
                        std::make_shared< SingleDependentVariableSaveSettings >(
                            body_fixed_relative_spherical_position,  "Apollo", "Earth" ) );
            dependentVariablesToAdd.push_back(
                        std::make_shared< SingleDependentVariableSaveSettings >(
                            rsw_to_inertial_frame_rotation_dependent_variable,  "Apollo", "Earth" ) );
            dependentVariablesToAdd.push_back(
                        std::make_shared< CustomDependentVariableSaveSettings >(
                            std::bind( &customDependentVariable1, bodies ), 3 ) );
            dependentVariablesToAdd.push_back(
                        std::make_shared< CustomDependentVariableSaveSettings >(
                            std::bind( &customDependentVariable2, bodies ), 1 ) );
            dependentVariablesToAdd.push_back(
                        std::make_shared< SingleDependentVariableSaveSettings >(
                            gravity_field_potential_dependent_variable,  "Apollo", "Earth" ) );
            dependentVariablesToAdd.push_back(
                        std::make_shared< SingleDependentVariableSaveSettings >(
                            gravity_field_potential_dependent_variable,  "Apollo", "Moon" ) );
            dependentVariablesToAdd.push_back(
                        std::make_shared< SingleDependentVariableSaveSettings >(
                            body_fixed_relative_cartesian_position,  "Apollo", "Moon" ) );

            addDepedentVariableSettings< double, double >( dependentVariablesToAdd, propagatorSettings );


            std::shared_ptr< IntegratorSettings< > > integratorSettings =
                    std::make_shared< IntegratorSettings< > >
                    ( rungeKutta4, simulationStartEpoch, fixedStepSize );

            // Create simulation object and propagate dynamics.
            SingleArcDynamicsSimulator< > dynamicsSimulator(
                        bodies, integratorSettings, propagatorSettings, true, false, false );

            std::cout<<"Propagation finished "<<std::endl;

            // Retrieve numerical solutions for state and dependent variables
            std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > numericalSolution =
                    dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
            std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > rawNumericalSolution =
                    dynamicsSimulator.getEquationsOfMotionNumericalSolutionRaw( );
            std::map< double, Eigen::VectorXd > dependentVariableSolution =
                    dynamicsSimulator.getDependentVariableHistory( );

            std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > originalDependentVariables =
                dynamicsSimulator.getSingleArcPropagationResults( )->getOriginalDependentVariableSettings( );

            std::map <std::pair<int, int>, std::string> dependentVariableId =
                dynamicsSimulator.getSingleArcPropagationResults( )->getDependentVariableId( );

            std::map< std::pair< int, int >, std::shared_ptr< SingleDependentVariableSaveSettings > > orderedDependentVariableSettings =
            dynamicsSimulator.getSingleArcPropagationResults( )-> getOrderedDependentVariableSettings( );

            BOOST_CHECK_EQUAL(
                originalDependentVariables.size( ), propagatorSettings->getDependentVariablesToSave( ).size( ) );
            for( unsigned int k = 0; k < originalDependentVariables.size( ); k++ )
            {
                BOOST_CHECK_EQUAL( originalDependentVariables.at( k ), propagatorSettings->getDependentVariablesToSave( ).at( k ) );
            }


            // Iterate over results for dependent variables, and check against manually retrieved values.
            Eigen::Vector6d currentStateDerivative;
            Eigen::Vector3d manualCentralGravity;
            std::shared_ptr< ephemerides::RotationalEphemeris > earthRotationModel =
                    bodies.at( "Earth" )->getRotationalEphemeris( );
            std::shared_ptr< aerodynamics::AtmosphereModel > earthAtmosphereModel =
                    bodies.at( "Earth" )->getAtmosphereModel( );
            std::shared_ptr< gravitation::GravityFieldModel > earthGravityModel =
                    bodies.at( "Earth" )->getGravityFieldModel( );
            std::shared_ptr< gravitation::GravityFieldModel > moonGravityModel =
                    bodies.at( "Moon" )->getGravityFieldModel( );

            std::shared_ptr< aerodynamics::AtmosphericFlightConditions > apolloFlightConditions =
                    std::dynamic_pointer_cast< aerodynamics::AtmosphericFlightConditions >(
                        bodies.at( "Apollo" )->getFlightConditions( ) );
            std::shared_ptr< aerodynamics::AerodynamicCoefficientInterface > apolloCoefficientInterface =
                    bodies.at( "Apollo" )->getAerodynamicCoefficientInterface( );
            bodies.at( "Apollo" )->setIsBodyInPropagation( true );

            std::cout<<"Number of points: "<<dependentVariableSolution.size( )<<std::endl;
            for( std::map< double, Eigen::VectorXd >::iterator variableIterator = dependentVariableSolution.begin( );
                 variableIterator != dependentVariableSolution.end( ); variableIterator++ )
            {
                Eigen::VectorXd currentDependentVariables = variableIterator->second;

                double machNumber = currentDependentVariables( 0 );
                double altitude = currentDependentVariables( 1 );
                double relativeDistance = currentDependentVariables( 2 );
                double relativeSpeed= currentDependentVariables( 3 );
                double gravitationalAccelerationNorm = currentDependentVariables( 4 );
                double gLoad = currentDependentVariables( 5 );
                double stagnationPointHeatFlux = currentDependentVariables( 6 );
                double freestreamTemperature = currentDependentVariables( 7 );
                double geodeticLatitude = currentDependentVariables( 8 );
                double freestreamDensity = currentDependentVariables( 9 );
                double latitude = currentDependentVariables( 10 );
                double longitude = currentDependentVariables( 11 );
                double angleOfAttack = currentDependentVariables( 12 );
                double sideslipAngle = currentDependentVariables( 13 );
                double bankAngle = currentDependentVariables( 14 );

                Eigen::Vector3d relativePosition = currentDependentVariables.segment( 15, 3 );
                Eigen::Vector3d computedBodyFixedPosition =
                        earthRotationModel->getRotationToTargetFrame( variableIterator->first ) * relativePosition;
                Eigen::Vector3d computedSphericalBodyFixedPosition =
                        coordinate_conversions::convertCartesianToSpherical( computedBodyFixedPosition );

                Eigen::Vector3d relativeVelocity = currentDependentVariables.segment( 18, 3 );
                Eigen::Vector3d gravitationalAcceleration = currentDependentVariables.segment( 21, 3 );
                Eigen::Vector3d totalAcceleration = currentDependentVariables.segment( 24, 3 );
                Eigen::Vector3d momentCoefficients = currentDependentVariables.segment( 27, 3 );
                Eigen::Vector3d forceCoefficients = currentDependentVariables.segment( 30, 3 );
                Eigen::Vector3d aerodynamicAcceleration = currentDependentVariables.segment( 33, 3 );
                Eigen::Vector3d moonAcceleration1 = currentDependentVariables.segment( 36, 3 );
                Eigen::Vector3d moonAcceleration2 = currentDependentVariables.segment( 39, 3 );

                Eigen::Vector6d keplerElements =  currentDependentVariables.segment( 42, 6 );
                Eigen::Vector6d modifiedEquinoctialElements =  currentDependentVariables.segment( 48, 6 );
                Eigen::Vector3d bodyFixedCartesianPosition = currentDependentVariables.segment( 54, 3 );
                Eigen::Vector3d bodyFixedSphericalPosition = currentDependentVariables.segment( 57, 3 );
                Eigen::Matrix3d rswToInertialRotationMatrix =
                        propagators::getMatrixFromVectorRotationRepresentation( currentDependentVariables.segment( 60, 9 ) );
                Eigen::Vector3d customVariable1 = currentDependentVariables.segment( 69, 3 );
                Eigen::Vector1d customVariable2 = currentDependentVariables.segment( 72, 1 );

                double earthGravitationalPotential = variableIterator->second( 73 );
                double moonGravitationalPotential = variableIterator->second( 74 );
                Eigen::Vector3d moonBodyFixedCartesianPosition = variableIterator->second.segment( 75, 3 );

                currentStateDerivative = dynamicsSimulator.getDynamicsStateDerivative( )->computeStateDerivative(
                            variableIterator->first, rawNumericalSolution.at( variableIterator->first ) );

                // Manually compute central gravity.
                manualCentralGravity =
                        -bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( ) *
                        relativePosition /
                        std::pow( relativePosition.norm( ), 3 );

                // Check output time consistency
                BOOST_CHECK_EQUAL( numericalSolution.count( variableIterator->first ), 1 );

                // Check relative position and velocity against state
                for( unsigned int i = 0; i < 3; i++ )
                {
                    BOOST_CHECK_SMALL(
                                std::fabs( numericalSolution.at( variableIterator->first )( i ) -
                                           relativePosition( i ) ), 2.0E-5 );
                    BOOST_CHECK_SMALL(
                                std::fabs( numericalSolution.at( variableIterator->first )( 3 + i ) -
                                           relativeVelocity( i ) ), 5.0E-11 );
                }

                // Check central gravity acceleration
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                            manualCentralGravity.segment( 0, 3 ),
                            gravitationalAcceleration, ( 6.0 * std::numeric_limits< double >::epsilon( ) ) );

                if( propagatorType == 0 )
                {
                    // Check total acceleration (tolerance is not epsilon due to numerical root finding for trim)
                    for( unsigned int i = 0; i < 3; i++ )
                    {
                        BOOST_CHECK_SMALL(
                                    std::fabs( currentStateDerivative( 3 + i ) - totalAcceleration( i ) ), 1.0E-13 );
                    }
                    cowellAcceleration[ variableIterator->first ] = totalAcceleration;
                }
                else if( variableIterator->first < 100.0 )
                {
                    for( unsigned int i = 0; i < 3; i++ )
                    {
                        BOOST_CHECK_SMALL(
                                    std::fabs( cowellAcceleration[ variableIterator->first ]( i ) - totalAcceleration( i ) ), 1.0E-8 );
                    }
                }

                // Check relative position and velocity norm.
                BOOST_CHECK_SMALL(
                            std::fabs( ( numericalSolution.at( variableIterator->first ).segment( 0, 3 ) ).norm( ) -
                                       relativeDistance ), 2.0E-5 );
                BOOST_CHECK_SMALL(
                            std::fabs( ( numericalSolution.at( variableIterator->first ).segment( 3, 3 ) ).norm( ) -
                                       relativeSpeed ), 2.0E-11 );

                // Check central gravity acceleration norm
                BOOST_CHECK_CLOSE_FRACTION(
                            manualCentralGravity.norm( ),
                            gravitationalAccelerationNorm, 6.0 * std::numeric_limits< double >::epsilon( ) );

                // Check Mach number
                BOOST_CHECK_CLOSE_FRACTION(
                            apolloFlightConditions->getCurrentAirspeed( ) /
                            apolloFlightConditions->getCurrentSpeedOfSound( ),
                            machNumber , std::numeric_limits< double >::epsilon( ) );

                // Check altitude.
                BOOST_CHECK_CLOSE_FRACTION(
                            apolloFlightConditions->getCurrentAltitude( ),
                            altitude, std::numeric_limits< double >::epsilon( ) );

                // Check latitude and longitude
                BOOST_CHECK_SMALL(
                            std::fabs( mathematical_constants::PI / 2.0 - computedSphericalBodyFixedPosition( 1 ) - latitude ),
                            6.0 * std::numeric_limits< double >::epsilon( ) );
                BOOST_CHECK_SMALL(
                            std::fabs( computedSphericalBodyFixedPosition( 2 ) - longitude ),
                            8.0 * std::numeric_limits< double >::epsilon( ) );

                // Check geodetic latitude.
                if( !isOblateSpheroidUsed )
                {
                    BOOST_CHECK_SMALL(
                                std::fabs( latitude - geodeticLatitude ),
                                6.0 * std::numeric_limits< double >::epsilon( ) );
                }
                else
                {
                    double computedGeodeticLatitude = coordinate_conversions::calculateGeodeticLatitude(
                                computedBodyFixedPosition, oblateSpheroidEquatorialRadius, oblateSpheroidFlattening, 1.0E-4 );
                    BOOST_CHECK_SMALL(
                                std::fabs( computedGeodeticLatitude - geodeticLatitude ),
                                6.0 * std::numeric_limits< double >::epsilon( ) );
                }
                BOOST_CHECK_SMALL(
                            std::fabs( geodeticLatitude - apolloFlightConditions->getCurrentGeodeticLatitude( ) ),
                            6.0 * std::numeric_limits< double >::epsilon( ) );

                // Check sideslip anf bank angle
                BOOST_CHECK_SMALL(
                            std::fabs( sideslipAngle ),
                            6.0 * std::numeric_limits< double >::epsilon( ) );
                BOOST_CHECK_SMALL(
                            std::fabs( bankAngle ),
                            6.0 * std::numeric_limits< double >::epsilon( ) );

                // Check g-load
                BOOST_CHECK_CLOSE_FRACTION(
                            gLoad, aerodynamicAcceleration.norm( ) / physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION,
                            6.0 * std::numeric_limits< double >::epsilon( ) );

                // Check density and temperature
                BOOST_CHECK_CLOSE_FRACTION(
                            freestreamDensity,
                            earthAtmosphereModel->getDensity( altitude, longitude, latitude, variableIterator->first ),
                            6.0 * std::numeric_limits< double >::epsilon( ) );
                BOOST_CHECK_CLOSE_FRACTION(
                            freestreamDensity,
                            apolloFlightConditions->getCurrentDensity( ),
                            6.0 * std::numeric_limits< double >::epsilon( ) );
                BOOST_CHECK_CLOSE_FRACTION(
                            freestreamTemperature,
                            earthAtmosphereModel->getTemperature( altitude, longitude, latitude, variableIterator->first ),
                            6.0 * std::numeric_limits< double >::epsilon( ) );
                BOOST_CHECK_CLOSE_FRACTION(
                            freestreamTemperature,
                            apolloFlightConditions->getCurrentFreestreamTemperature( ),
                            6.0 * std::numeric_limits< double >::epsilon( ) );

                // Check aerodynamic force coefficients
                std::vector< double > aerodynamicCoefficientInput;
                aerodynamicCoefficientInput.push_back( machNumber );
                aerodynamicCoefficientInput.push_back( angleOfAttack );
                aerodynamicCoefficientInput.push_back( sideslipAngle );
                apolloCoefficientInterface->updateCurrentCoefficients( aerodynamicCoefficientInput );
                Eigen::Vector3d computedForceCoefficients = apolloCoefficientInterface->getCurrentForceCoefficients( );

                BOOST_CHECK_SMALL( std::fabs( computedForceCoefficients( 1 ) ), 1.0E-14 );
                for( unsigned int i = 0; i < 3; i ++ )
                {
                    BOOST_CHECK_SMALL(
                                std::fabs( computedForceCoefficients( i ) - forceCoefficients( i ) ),
                                6.0 * std::numeric_limits< double >::epsilon( ) );
                }

                // Check heat flux
                double expectedHeatFlux = aerodynamics::computeEquilibriumFayRiddellHeatFlux(
                            freestreamDensity, apolloFlightConditions->getCurrentAirspeed( ),
                            freestreamTemperature, machNumber, noseRadius, wallEmissivity );
                BOOST_CHECK_CLOSE_FRACTION(
                            expectedHeatFlux, stagnationPointHeatFlux,
                            6.0 * std::numeric_limits< double >::epsilon( ) );

                // Check trimmed condition (y-term)/symmetric vehicle shape (x- and z-term).
                BOOST_CHECK_SMALL(
                            std::fabs( momentCoefficients( 0 ) ), 1.0E-14 );
                BOOST_CHECK_SMALL(
                            std::fabs( momentCoefficients( 1 ) ), 1.0E-10 );
                BOOST_CHECK_SMALL(
                            std::fabs( momentCoefficients( 2 ) ), 1.0E-14 );

                // Check if third-body gravity saving is done correctly
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( moonAcceleration1, moonAcceleration2,
                                                   ( 2.0 * std::numeric_limits< double >::epsilon( ) ) );

                Eigen::Vector6d expectedKeplerElements =
                        tudat::orbital_element_conversions::convertCartesianToKeplerianElements(
                            Eigen::Vector6d( numericalSolution.at( variableIterator->first ) ), earthGravitationalParameter );

                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedKeplerElements, keplerElements,
                                                   ( 10.0 * std::numeric_limits< double >::epsilon( ) ) );

                Eigen::Vector6d expectedModifiedEquinoctialElements =
                        tudat::orbital_element_conversions::convertCartesianToModifiedEquinoctialElements(
                            Eigen::Vector6d( numericalSolution.at( variableIterator->first ) ), earthGravitationalParameter );

                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedModifiedEquinoctialElements, modifiedEquinoctialElements,
                                                   ( 10.0 * std::numeric_limits< double >::epsilon( ) ) );

                for( unsigned int i = 0; i < 3; i ++ )
                {
                    BOOST_CHECK_SMALL(
                                std::fabs( computedBodyFixedPosition( i ) - bodyFixedCartesianPosition( i ) ),
                                1.0E-8 );
                }
                BOOST_CHECK_SMALL(
                            std::fabs( computedSphericalBodyFixedPosition( 0 ) - bodyFixedSphericalPosition( 0 ) ),
                            1.0E-8 );
                BOOST_CHECK_SMALL(
                            std::fabs( tudat::mathematical_constants::PI / 2.0 - computedSphericalBodyFixedPosition( 1 ) -
                                       bodyFixedSphericalPosition( 1 ) ),
                            10.0 * std::numeric_limits< double >::epsilon( ) );
                BOOST_CHECK_SMALL(
                            std::fabs( computedSphericalBodyFixedPosition( 2 ) - bodyFixedSphericalPosition( 2 ) ),
                            10.0 * std::numeric_limits< double >::epsilon( ) );
                Eigen::Matrix3d computedRswRotationMatrix =
                        tudat::reference_frames::getRswSatelliteCenteredToInertialFrameRotationMatrix( numericalSolution.at( variableIterator->first ) );

                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( rswToInertialRotationMatrix, computedRswRotationMatrix,
                                                   ( 10.0 * std::numeric_limits< double >::epsilon( ) ) );
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(  customVariable1, ( relativePosition / 2.0 ),
                                                    ( 10.0 * std::numeric_limits< double >::epsilon( ) ) );
                BOOST_CHECK_CLOSE_FRACTION(  customVariable2( 0 ), ( relativePosition.norm( ) / 2.0 ),
                                             ( 10.0 * std::numeric_limits< double >::epsilon( ) ) );


                // Check central body gravitational potential
                BOOST_CHECK_CLOSE_FRACTION(
                            earthGravitationalPotential,
                            earthGravityModel->getGravitationalPotential( computedBodyFixedPosition ),
                            6.0 * std::numeric_limits< double >::epsilon( ) );

                // Check 3rd body gravitational potential - not sure why the tolerance needs to be so large
                BOOST_CHECK_CLOSE_FRACTION(
                            moonGravitationalPotential,
                            moonGravityModel->getGravitationalPotential( moonBodyFixedCartesianPosition ),
                            1e-8 );
            }
        }
    }
}

//! Function to test whether separate spherical harmonic acceleration contributions are correctly saved.
BOOST_AUTO_TEST_CASE( testSphericalHarmonicDependentVariableOutput )
{
    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Create body objects.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    BodyListSettings bodySettings = getDefaultBodySettings( bodiesToCreate, "Earth", "ECLIPJ2000" );

    // Create Body objects
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    bodies.createEmptyBody( "Asterix" );

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    bodiesToPropagate.push_back( "Asterix" );
    centralBodies.push_back( "Earth" );

    // Define propagation settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfAsterix;
    accelerationsOfAsterix[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >(
                                                     6, 6 ) );
    accelerationMap[ "Asterix" ] = accelerationsOfAsterix;

    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodies, accelerationMap, bodiesToPropagate, centralBodies );

    // Set Keplerian elements for Asterix.
    Eigen::Vector6d asterixInitialStateInKeplerianElements;
    asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7500.0E3;
    asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.1;
    asterixInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 85.3 );
    asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
            = convertDegreesToRadians( 235.7 );
    asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
            = convertDegreesToRadians( 23.4 );
    asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 139.87 );

    // Convert Asterix state from Keplerian elements to Cartesian elements.
    double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
    Eigen::VectorXd systemInitialState = convertKeplerianToCartesianElements(
                asterixInitialStateInKeplerianElements,
                earthGravitationalParameter );

    double simulationEndEpoch = 10.0;
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, simulationEndEpoch, cowell );

    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;

    dependentVariables.push_back(
                std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    spherical_harmonic_gravity, "Asterix", "Earth" ) );

    std::vector< std::pair< int, int > > separateTermsToSave;
    separateTermsToSave.push_back( std::make_pair( 1, 0 ) );
    separateTermsToSave.push_back( std::make_pair( 3, 0 ) );
    separateTermsToSave.push_back( std::make_pair( 3, 2 ) );
    separateTermsToSave.push_back( std::make_pair( 3, 1 ) );
    dependentVariables.push_back(
                std::make_shared< SphericalHarmonicAccelerationTermsDependentVariableSaveSettings >(
                    "Asterix", "Earth", separateTermsToSave ) );

    dependentVariables.push_back(
                std::make_shared< SphericalHarmonicAccelerationTermsDependentVariableSaveSettings >(
                    "Asterix", "Earth", 6, 6 ) );
    std::vector< std::pair< int, int > > singleTermToSave;
    singleTermToSave.push_back( std::make_pair( 6, 6 ) );
    dependentVariables.push_back( sphericalHarmonicAccelerationTermsNormDependentVariable(
                                      "Asterix", "Earth", singleTermToSave ) );

    addDepedentVariableSettings< double, double >( dependentVariables, propagatorSettings );

    // Create numerical integrator.
    double simulationStartEpoch = 0.0;
    const double fixedStepSize = 10.0;
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4, simulationStartEpoch, fixedStepSize );

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodies, integratorSettings, propagatorSettings );

    std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > depdendentVariableResult = dynamicsSimulator.getDependentVariableHistory( );

    for( std::map< double, Eigen::VectorXd >::iterator variableIterator = depdendentVariableResult.begin( );
         variableIterator != depdendentVariableResult.end( ); variableIterator++ )
    {
        Eigen::VectorXd currentDependentVariables = variableIterator->second;

        Eigen::Vector3d manualAccelerationSum = Eigen::Vector3d::Zero( );
        for( unsigned int i = 5; i < 33; i++ )
        {
            manualAccelerationSum += currentDependentVariables.segment( i * 3, 3 );
        }

        BOOST_CHECK_SMALL(
                    std::fabs( ( currentDependentVariables.segment( 96, 3 ) ).norm( ) - currentDependentVariables( 99 ) ),
                    10.0 * ( currentDependentVariables.segment( 96, 3 ) ).norm( ) * std::numeric_limits< double >::epsilon( ) );
        for( unsigned int i = 0; i < 3; i ++ )
        {
            BOOST_CHECK_SMALL(
                        std::fabs( manualAccelerationSum( i ) - currentDependentVariables( i ) ),
                        10.0 * manualAccelerationSum.norm( ) * std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_SMALL(
                        std::fabs( currentDependentVariables( 3 + i ) - currentDependentVariables( 15 + 3 * 1 + i ) ),
                        10.0 * ( currentDependentVariables.segment( 3, 3 ) ).norm( ) * std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_SMALL(
                        std::fabs( currentDependentVariables( 6 + i ) - currentDependentVariables( 15 + 3 * 6 + i ) ),
                        10.0 * ( currentDependentVariables.segment( 6, 3 ) ).norm( ) * std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_SMALL(
                        std::fabs( currentDependentVariables( 9 + i ) - currentDependentVariables( 15 + 3 * 8 + i ) ),
                        10.0 * ( currentDependentVariables.segment( 9, 3 ) ).norm( ) * std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_SMALL(
                        std::fabs( currentDependentVariables( 12 + i ) - currentDependentVariables( 15 + 3 * 7 + i ) ),
                        10.0 * ( currentDependentVariables.segment( 12, 3 ) ).norm( ) * std::numeric_limits< double >::epsilon( ) );

        }
    }
}

Eigen::VectorXd getCustomDependentVariable(
        const SystemOfBodies& bodies )
{
    Eigen::Vector6d sunState = bodies.at( "Sun" )->getState( );
    Eigen::Vector6d moonState = bodies.at( "Moon" )->getState( );
    return sunState.cwiseQuotient( moonState );

}

BOOST_AUTO_TEST_CASE( testDependentVariableEnvironmentUpdate )
{
    std::string kernelsPath = paths::getSpiceKernelPath( );
    spice_interface::loadStandardSpiceKernels( );

    // Define bodies in simulation.
    unsigned int totalNumberOfBodies = 4;
    std::vector< std::string > bodyNames;
    bodyNames.resize( totalNumberOfBodies );
    bodyNames[ 0 ] = "Moon";
    bodyNames[ 1 ] = "Earth";
    bodyNames[ 2 ] = "Venus";
    bodyNames[ 3 ] = "Sun";

    // Create bodies needed in simulation
    BodyListSettings bodySettings =
            getDefaultBodySettings( bodyNames );
    bodySettings.at("Sun")->radiationSourceModelSettings =
            isotropicPointRadiationSourceModelSettings(
                    irradianceBasedLuminosityModelSettings(1371, physical_constants::ASTRONOMICAL_UNIT));
    auto earthRadius = spice_interface::getAverageRadius("Earth");
    bodySettings.at("Earth")->radiationPressureTargetModelSettings =
            cannonballRadiationPressureTargetModelSettings(
                    4 * mathematical_constants::PI * earthRadius * earthRadius, 1.5);
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );


    SelectedAccelerationMap accelerationMap;
    accelerationMap[ "Earth" ][ "Moon" ].push_back(
                std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationMap[ "Earth" ][ "Sun" ].push_back(
                std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationMap[ "Earth" ][ "Sun" ].push_back(
                std::make_shared< AccelerationSettings >( radiation_pressure ) );

    std::vector< std::string > bodiesToPropagate;
    bodiesToPropagate.push_back( "Earth" );

    std::vector< std::string > centralBodies;
    centralBodies.push_back( "SSB" );

    // Create acceleration models and propagation settings.
    AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodies, accelerationMap, bodiesToPropagate, centralBodies );

    // Specify initial time
    double initialEphemerisTime = 0.0;
    double finalEphemerisTime = 28.0 * 86400.0;

    // Get initial state vector as input to integration.
    Eigen::VectorXd systemInitialState = getInitialStatesOfBodies(
                bodiesToPropagate, centralBodies, bodies, initialEphemerisTime );

    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;

    dependentVariables.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    relative_position_dependent_variable, "Sun", "Venus" ) );
    dependentVariables.push_back(
                std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    "Moon", reference_frames::latitude_angle, "Earth" ) );
    dependentVariables.push_back(
                std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    "Moon", reference_frames::longitude_angle, "Earth" ) );
    dependentVariables.push_back(
                std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    "Moon", reference_frames::heading_angle, "Earth" ) );
    dependentVariables.push_back(
                std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    "Moon", reference_frames::flight_path_angle, "Earth" ) );
    // Earth-Sun distance is not tested, only used to calculate expected received irradiance
    dependentVariables.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    relative_distance_dependent_variable, "Earth", "Sun" ) );
    dependentVariables.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    received_irradiance, "Earth", "Sun" ) );
    dependentVariables.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    received_fraction, "Earth", "Sun" ) );
    dependentVariables.push_back(
                std::make_shared< CustomDependentVariableSaveSettings >(
                    [=]( ){ return getCustomDependentVariable( bodies ); }, 6 ) );


    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, finalEphemerisTime, cowell,
              dependentVariables  );

    // Define numerical integrator settings.
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4, initialEphemerisTime, 3600.0 );


    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodies, integratorSettings, propagatorSettings, true, false, false );

    std::map< double, Eigen::VectorXd > stateResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > depdendentVariableResult = dynamicsSimulator.getDependentVariableHistory( );

    for( std::map< double, Eigen::VectorXd >::iterator variableIterator = depdendentVariableResult.begin( );
         variableIterator != depdendentVariableResult.end( ); variableIterator++ )
    {
        Eigen::VectorXd currentDependentVariables = variableIterator->second;

        Eigen::Vector3d expectedRelativePosition =
                tudat::spice_interface::getBodyCartesianPositionAtEpoch(
                    "Sun", "Venus", "ECLIPJ2000", "None", variableIterator->first );
        Eigen::Vector3d computedRelativePosition = currentDependentVariables.segment( 0, 3 );


        for( unsigned int i = 0; i < 3; i ++ )
        {
            BOOST_CHECK_SMALL(
                        std::fabs( expectedRelativePosition( i ) - computedRelativePosition( i ) ), 1.0E-4 );
        }

        Eigen::Vector6d moonRelativeCartesianState = tudat::spice_interface::getBodyCartesianStateAtEpoch(
                    "Moon", "SSB", "ECLIPJ2000", "None", variableIterator->first ) -
                stateResult.at( variableIterator->first );
        Eigen::Vector6d moonRelativeEarthFixedCartesianState = ephemerides::transformStateToTargetFrame(
                    moonRelativeCartesianState, variableIterator->first, bodies.at( "Earth" )->getRotationalEphemeris( ) );
        Eigen::Vector3d moonSphericalPosition = tudat::coordinate_conversions::convertCartesianToSpherical< double >(
                    Eigen::Vector3d( moonRelativeEarthFixedCartesianState.segment( 0, 3 ) ) );

        BOOST_CHECK_SMALL(
                    std::fabs( currentDependentVariables( 3 ) - (
                                   tudat::mathematical_constants::PI / 2.0 - moonSphericalPosition( 1 ) ) ), 1.0E-14 );
        BOOST_CHECK_SMALL(
                    std::fabs( currentDependentVariables( 4 ) - moonSphericalPosition( 2 ) ), 1.0E-14 );

        Eigen::Vector6d moonRelativeSphericalState =
                tudat::orbital_element_conversions::convertCartesianToSphericalOrbitalState(
                    moonRelativeEarthFixedCartesianState );

        BOOST_CHECK_SMALL(
                    std::fabs( currentDependentVariables( 3 ) - moonRelativeSphericalState(
                                   orbital_element_conversions::latitudeIndex ) ), 1.0E-14 );
        BOOST_CHECK_SMALL(
                    std::fabs( currentDependentVariables( 4 ) - moonRelativeSphericalState(
                                   orbital_element_conversions::longitudeIndex ) ), 1.0E-14 );
        BOOST_CHECK_SMALL(
                    std::fabs( currentDependentVariables( 5 ) - moonRelativeSphericalState(
                                   orbital_element_conversions::headingAngleIndex ) ), 1.0E-14 );
        BOOST_CHECK_SMALL(
                    std::fabs( currentDependentVariables( 6 ) - moonRelativeSphericalState(
                                   orbital_element_conversions::flightPathIndex ) ), 1.0E-14 );

        double earthSunDistance = currentDependentVariables ( 7 );
        double expectedReceivedIrradiance =
                1371 * physical_constants::ASTRONOMICAL_UNIT * physical_constants::ASTRONOMICAL_UNIT
                / (earthSunDistance * earthSunDistance);
        double computedReceivedIrradiance = currentDependentVariables ( 8 );
        BOOST_CHECK_SMALL(std::fabs( expectedReceivedIrradiance - computedReceivedIrradiance ), 1.0E-10 );

        double expectedReceivedFraction = 1.0;
        double computedReceivedFraction = currentDependentVariables ( 9 );
        BOOST_CHECK_SMALL(std::fabs( expectedReceivedFraction - computedReceivedFraction ), 1.0E-10 );

        Eigen::Vector6d sunState = bodies.at( "Sun" )->getStateInBaseFrameFromEphemeris( variableIterator->first );
        Eigen::Vector6d moonState = bodies.at( "Moon" )->getStateInBaseFrameFromEphemeris(variableIterator->first  );
        Eigen::Vector6d customDependentVariable = sunState.cwiseQuotient( moonState );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( customDependentVariable, ( currentDependentVariables.segment( 10, 6 ) ), 1.0E-14 );
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
            std::make_shared< BasicSolidBodyGravityFieldVariationSettings >(
                deformingBodies, loveNumbers );
    gravityFieldVariations.push_back( moonGravityFieldVariation );

    deformingBodies[ 0 ] = "Sun";
    std::shared_ptr< GravityFieldVariationSettings > sunSingleGravityFieldVariation =
            std::make_shared< BasicSolidBodyGravityFieldVariationSettings >(
                deformingBodies, loveNumbers );
    gravityFieldVariations.push_back( sunSingleGravityFieldVariation );

    return gravityFieldVariations;
}

//! Unit test to check if acceleration contributions due to gravity field variations are being correctly stored
BOOST_AUTO_TEST_CASE( test_GravityFieldVariationAccelerationSaving )
{

    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Define bodies in simulation
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Moon" );
    bodyNames.push_back( "Mars" );

    // Specify initial time
    double initialEphemerisTime = 1.0E7;
    double finalEphemerisTime = initialEphemerisTime + 3000.0;

    // Create bodies needed in simulation
    BodyListSettings bodySettings =
            getDefaultBodySettings( bodyNames,  "Earth", "ECLIPJ2000" );
    bodySettings.at( "Earth" )->gravityFieldVariationSettings = getEarthGravityFieldVariationSettings( );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    bodies.createEmptyBody( "Vehicle" );

    // Set accelerations on Vehicle that are to be taken into account.
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 3, 3 ) );
    accelerationsOfVehicle[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                   basic_astrodynamics::point_mass_gravity ) );
    accelerationsOfVehicle[ "Moon" ].push_back( std::make_shared< AccelerationSettings >(
                                                    basic_astrodynamics::point_mass_gravity ) );
    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;

    // Set bodies for which initial state is to be estimated and integrated.
    std::vector< std::string > bodiesToIntegrate;
    std::vector< std::string > centralBodies;
    bodiesToIntegrate.push_back( "Vehicle" );
    centralBodies.push_back( "Earth" );

    // Create acceleration models
    AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodies, accelerationMap, bodiesToIntegrate, centralBodies );

    // Set Keplerian and Cartesian elements for spacecraft.
    Eigen::Vector6d asterixInitialStateInKeplerianElements;
    asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 6600.0E3;
    asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.05;
    asterixInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
    asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
            = unit_conversions::convertDegreesToRadians( 235.7 );
    asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
            = unit_conversions::convertDegreesToRadians( 23.4 );
    asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );
    double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
    Eigen::Vector6d systemInitialState = convertKeplerianToCartesianElements(
                asterixInitialStateInKeplerianElements, earthGravitationalParameter );

    // Create propagator settings
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
    dependentVariables.push_back(
                std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    spherical_harmonic_gravity, "Vehicle", "Earth" ) );
    dependentVariables.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    total_gravity_field_variation_acceleration, "Vehicle", "Earth" ) );
    dependentVariables.push_back(
                std::make_shared< SingleVariationSphericalHarmonicAccelerationSaveSettings >(
                    "Vehicle", "Earth", gravitation::basic_solid_body, "Sun" ) );
    dependentVariables.push_back(
                std::make_shared< SingleVariationSphericalHarmonicAccelerationSaveSettings >(
                    "Vehicle", "Earth", gravitation::basic_solid_body, "Moon" ) );
    dependentVariables.push_back(
                std::make_shared< SingleVariationSingleTermSphericalHarmonicAccelerationSaveSettings >(
                    "Vehicle", "Earth", 3, 3, gravitation::basic_solid_body, "Moon" ) );
    std::vector< std::pair< int, int > > componentIndices = { { 1, 0 }, { 2, 0 }, { 2, 2 }, { 3, 1 } };
    dependentVariables.push_back(
                std::make_shared< SingleVariationSingleTermSphericalHarmonicAccelerationSaveSettings >(
                    "Vehicle", "Earth", componentIndices, gravitation::basic_solid_body, "Moon" ) );
    dependentVariables.push_back(
                std::make_shared< TotalGravityFieldVariationSettings >(
                    "Earth", 1, 3, 0, 3, true ) );
    dependentVariables.push_back(
                std::make_shared< TotalGravityFieldVariationSettings >(
                    "Earth", 1, 3, 1, 3, false ) );

    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToIntegrate, systemInitialState,
              double( finalEphemerisTime ), cowell,
              dependentVariables );

    // Create integrator settings
    std::shared_ptr< IntegratorSettings< double > > integratorSettings =
            std::make_shared< RungeKuttaVariableStepSizeSettings< double > >
            ( double( initialEphemerisTime ), 300.0,
              CoefficientSets::rungeKuttaFehlberg78,
              300.0, 300.0, 1.0, 1.0 );

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodies, integratorSettings, propagatorSettings, true, false, false );

    // Retrieve numerical solutions for state and dependent variables
    std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > numericalSolution =
            dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > dependentVariableSolution =
            dynamicsSimulator.getDependentVariableHistory( );

    std::shared_ptr< gravitation::SphericalHarmonicsGravitationalAccelerationModel >
            sphericalHarmonicAcceleration =
            std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravitationalAccelerationModel >(
                getAccelerationBetweenBodies(
                    "Vehicle", "Earth", dynamicsSimulator.getDynamicsStateDerivative( )->getStateDerivativeModels( ),
                    basic_astrodynamics::spherical_harmonic_gravity ).at( 0 ) );
    std::shared_ptr< gravitation::TimeDependentSphericalHarmonicsGravityField > earthGravityField =
            std::dynamic_pointer_cast< gravitation::TimeDependentSphericalHarmonicsGravityField >(
                bodies.at( "Earth" )->getGravityFieldModel( ) );


    // Iterate over results for dependent variables, and check against manually retrieved values.
    Eigen::Vector6d currentStateDerivative;
    for( std::map< double, Eigen::VectorXd >::iterator variableIterator = dependentVariableSolution.begin( );
         variableIterator != dependentVariableSolution.end( ); variableIterator++ )
    {
        currentStateDerivative = dynamicsSimulator.getDynamicsStateDerivative( )->computeStateDerivative(
                    variableIterator->first, numericalSolution.at( variableIterator->first ) );

        Eigen::Vector3d currentTotalAcceleration = dependentVariableSolution.at( variableIterator->first ).segment( 0, 3 );
        Eigen::Vector3d currentTotalTidalCorrection = dependentVariableSolution.at( variableIterator->first ).segment( 3, 3 );
        Eigen::Vector3d currentTotalSunTidalCorrection = dependentVariableSolution.at( variableIterator->first ).segment( 6, 3 );
        Eigen::Vector3d currentTotalMoonTidalCorrection = dependentVariableSolution.at( variableIterator->first ).segment( 9, 3 );
        Eigen::VectorXd perTermMoonTidalCorrections = dependentVariableSolution.at( variableIterator->first ).segment( 12, 30 );
        Eigen::VectorXd perTermMoonTidalSelectedCorrections = dependentVariableSolution.at( variableIterator->first ).segment( 42, 12 );
        Eigen::VectorXd sphericalHarmonicCosineCoefficientCorrection = dependentVariableSolution.at( variableIterator->first ).segment( 54, 9 );
        Eigen::VectorXd sphericalHarmonicSineCoefficientCorrection = dependentVariableSolution.at( variableIterator->first ).segment( 63, 6 );

        Eigen::Vector3d reconstructedTotalMoonTidalCorrection = Eigen::Vector3d::Zero( );
        for( int j = 0; j < 10; j++ )
        {
            reconstructedTotalMoonTidalCorrection += perTermMoonTidalCorrections.segment( j * 3, 3 );
        }
        Eigen::Vector3d computedTotalCoefficientAcceleration =
                sphericalHarmonicAcceleration->getAccelerationWithAlternativeCoefficients(
                    earthGravityField->getCosineCoefficients( ).block( 0, 0, 4, 4 ),
                    earthGravityField->getSineCoefficients( ).block( 0, 0, 4, 4 ) );
        Eigen::Vector3d computedNominalCoefficientAcceleration =
                sphericalHarmonicAcceleration->getAccelerationWithAlternativeCoefficients(
                    earthGravityField->getNominalCosineCoefficients( ).block( 0, 0, 4, 4 ),
                    earthGravityField->getNominalSineCoefficients( ).block( 0, 0, 4, 4 ) );
        Eigen::Vector3d computedTotalTidalCorrection =
                computedTotalCoefficientAcceleration - computedNominalCoefficientAcceleration;

        Eigen::MatrixXd computedCosineCorrection = earthGravityField->getTotalCosineCoefficientCorrection( 3, 3 );
        Eigen::MatrixXd computedSineCorrection = earthGravityField->getTotalSineCoefficientCorrection( 3, 3 );

        for( int i = 0; i < 3; i++ )
        {
            BOOST_CHECK_SMALL( std::fabs( currentTotalSunTidalCorrection( i ) + currentTotalMoonTidalCorrection( i ) -
                                          currentTotalTidalCorrection( i ) ), 1.0E-14 );

            BOOST_CHECK_SMALL( std::fabs( currentTotalAcceleration( i ) - currentTotalTidalCorrection( i ) -
                                          computedNominalCoefficientAcceleration( i ) ), 1.0E-14 );
            BOOST_CHECK_SMALL( std::fabs( currentTotalTidalCorrection( i ) - computedTotalTidalCorrection( i ) ), 1.0E-14 );
            BOOST_CHECK_SMALL( std::fabs( reconstructedTotalMoonTidalCorrection( i ) - currentTotalMoonTidalCorrection( i ) ),
                               1.0E-14 );
            BOOST_CHECK_SMALL( std::fabs( perTermMoonTidalCorrections( i + 3 ) - perTermMoonTidalSelectedCorrections( i ) ),
                               1.0E-14 );
            BOOST_CHECK_SMALL( std::fabs( perTermMoonTidalCorrections( i + 9 ) - perTermMoonTidalSelectedCorrections( i + 3 ) ),
                               1.0E-14 );
            BOOST_CHECK_SMALL( std::fabs( perTermMoonTidalCorrections( i + 15 ) - perTermMoonTidalSelectedCorrections( i + 6 ) ),
                               1.0E-14 );
            BOOST_CHECK_SMALL( std::fabs( perTermMoonTidalCorrections( i + 21 ) - perTermMoonTidalSelectedCorrections( i + 9 ) ),
                               1.0E-14 );
        }

        BOOST_CHECK_SMALL( std::fabs( computedCosineCorrection( 1, 0 ) - sphericalHarmonicCosineCoefficientCorrection( 0 ) ), 1.0E-25 );
        BOOST_CHECK_SMALL( std::fabs( computedCosineCorrection( 1, 1 ) - sphericalHarmonicCosineCoefficientCorrection( 1 ) ), 1.0E-25 );
        BOOST_CHECK_SMALL( std::fabs( computedCosineCorrection( 2, 0 ) - sphericalHarmonicCosineCoefficientCorrection( 2 ) ), 1.0E-25 );
        BOOST_CHECK_SMALL( std::fabs( computedCosineCorrection( 2, 1 ) - sphericalHarmonicCosineCoefficientCorrection( 3 ) ), 1.0E-25 );
        BOOST_CHECK_SMALL( std::fabs( computedCosineCorrection( 2, 2 ) - sphericalHarmonicCosineCoefficientCorrection( 4 ) ), 1.0E-25 );
        BOOST_CHECK_SMALL( std::fabs( computedCosineCorrection( 3, 0 ) - sphericalHarmonicCosineCoefficientCorrection( 5 ) ), 1.0E-25 );
        BOOST_CHECK_SMALL( std::fabs( computedCosineCorrection( 3, 1 ) - sphericalHarmonicCosineCoefficientCorrection( 6 ) ), 1.0E-25 );
        BOOST_CHECK_SMALL( std::fabs( computedCosineCorrection( 3, 2 ) - sphericalHarmonicCosineCoefficientCorrection( 7 ) ), 1.0E-25 );
        BOOST_CHECK_SMALL( std::fabs( computedCosineCorrection( 3, 3 ) - sphericalHarmonicCosineCoefficientCorrection( 8 ) ), 1.0E-25 );

        BOOST_CHECK_SMALL( std::fabs( computedSineCorrection( 1, 1 ) - sphericalHarmonicSineCoefficientCorrection( 0 ) ), 1.0E-25 );
        BOOST_CHECK_SMALL( std::fabs( computedSineCorrection( 2, 1 ) - sphericalHarmonicSineCoefficientCorrection( 1 ) ), 1.0E-25 );
        BOOST_CHECK_SMALL( std::fabs( computedSineCorrection( 2, 2 ) - sphericalHarmonicSineCoefficientCorrection( 2 ) ), 1.0E-25 );
        BOOST_CHECK_SMALL( std::fabs( computedSineCorrection( 3, 1 ) - sphericalHarmonicSineCoefficientCorrection( 3 ) ), 1.0E-25 );
        BOOST_CHECK_SMALL( std::fabs( computedSineCorrection( 3, 2 ) - sphericalHarmonicSineCoefficientCorrection( 4 ) ), 1.0E-25 );
        BOOST_CHECK_SMALL( std::fabs( computedSineCorrection( 3, 3 ) - sphericalHarmonicSineCoefficientCorrection( 5 ) ), 1.0E-25 );


    }
}

//! Unit test to check if acceleration contributions due to gravity field variations are being correctly stored
BOOST_AUTO_TEST_CASE( test_AccelerationPartialSaving )
{

    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Define bodies in simulation
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Moon" );
    bodyNames.push_back( "Mars" );

    // Specify initial time
    double initialEphemerisTime = 1.0E7;
    double finalEphemerisTime = initialEphemerisTime + 300.0;

    // Create bodies needed in simulation
    BodyListSettings bodySettings =
            getDefaultBodySettings( bodyNames, "Earth", "ECLIPJ2000" );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    bodies.createEmptyBody( "Vehicle" );

    // Set accelerations on Vehicle that are to be taken into account.
    for( int test = 0; test < 5; test++ )
    {
        SelectedAccelerationMap accelerationMap;
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;

        if( test == 0 || test == 3 || test == 4 )
        {
            accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 3, 3 ) );
        }

        if( test > 0 )
        {
            accelerationsOfVehicle[ "Moon" ].push_back( std::make_shared< AccelerationSettings >(
                                                            basic_astrodynamics::point_mass_gravity ) );
        }

        if( test > 2 )
        {
            accelerationsOfVehicle[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                           basic_astrodynamics::point_mass_gravity ) );
        }

        accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;

        // Set bodies for which initial state is to be estimated and integrated.
        std::vector< std::string > bodiesToIntegrate;
        std::vector< std::string > centralBodies;
        bodiesToIntegrate.push_back( "Vehicle" );
        centralBodies.push_back( "Earth" );

        // Create acceleration models
        AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodies, accelerationMap, bodiesToIntegrate, centralBodies );

        // Set Keplerian and Cartesian elements for spacecraft.
        Eigen::Vector6d asterixInitialStateInKeplerianElements;
        asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 9000.0E3;
        asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.05;
        asterixInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
        asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
                = unit_conversions::convertDegreesToRadians( 235.7 );
        asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
                = unit_conversions::convertDegreesToRadians( 23.4 );
        asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );
        double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
        Eigen::Vector6d systemInitialState = convertKeplerianToCartesianElements(
                    asterixInitialStateInKeplerianElements, earthGravitationalParameter );

        // Create propagator settings
        std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
        if( test == 0 || test == 3 || test == 4 )
        {
            dependentVariables.push_back(
                        std::make_shared< AccelerationPartialWrtStateSaveSettings >(
                            "Vehicle", "Earth", spherical_harmonic_gravity, "Vehicle" ) );
        }

        if( test == 0 )
        {
            dependentVariables.push_back(
                        std::make_shared< AccelerationPartialWrtStateSaveSettings >(
                            "Vehicle", "Earth", spherical_harmonic_gravity, "Earth" ) );
            dependentVariables.push_back(
                        std::make_shared< AccelerationPartialWrtStateSaveSettings >(
                            "Vehicle", "Earth", spherical_harmonic_gravity, "Moon" ) );
        }

        if( test == 1 || test == 3 )
        {
            dependentVariables.push_back(
                        std::make_shared< AccelerationPartialWrtStateSaveSettings >(
                            "Vehicle", "Moon", point_mass_gravity, "Vehicle" ) );
            dependentVariables.push_back(
                        std::make_shared< AccelerationPartialWrtStateSaveSettings >(
                            "Vehicle", "Moon", point_mass_gravity, "Moon" ) );
        }
        else if( test == 2 || test == 4 )
        {
            dependentVariables.push_back(
                        std::make_shared< AccelerationPartialWrtStateSaveSettings >(
                            "Vehicle", "Moon", third_body_point_mass_gravity, "Vehicle" ) );
            dependentVariables.push_back(
                        std::make_shared< AccelerationPartialWrtStateSaveSettings >(
                            "Vehicle", "Moon", third_body_point_mass_gravity, "Moon" ) );
        }

        if( test == 3 )
        {
            dependentVariables.push_back(
                        std::make_shared< AccelerationPartialWrtStateSaveSettings >(
                            "Vehicle", "Sun", point_mass_gravity, "Vehicle" ) );
        }
        else if( test == 4 )
        {
            dependentVariables.push_back(
                        std::make_shared< AccelerationPartialWrtStateSaveSettings >(
                            "Vehicle", "Sun", third_body_point_mass_gravity, "Vehicle" ) );
        }
        dependentVariables.push_back(
                std::make_shared< TotalAccelerationPartialWrtStateSaveSettings >( "Vehicle", "Vehicle" ) );
        dependentVariables.push_back(
                std::make_shared< TotalAccelerationPartialWrtStateSaveSettings >( "Vehicle", "Moon" ) );

        std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToIntegrate, systemInitialState,
                  double( finalEphemerisTime ), cowell,
                  dependentVariables );

        // Create integrator settings
        std::shared_ptr< IntegratorSettings< double > > integratorSettings =
                std::make_shared< RungeKuttaVariableStepSizeSettings< double > >
                ( double( initialEphemerisTime ), 0.1,
                  CoefficientSets::rungeKuttaFehlberg78,
                  0.1, 0.1, 1.0, 1.0 );


        // Define list of parameters to estimate.
        std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
                getInitialStateParameterSettings< double >( propagatorSettings, bodies );

        // Create parameters
        std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
                createParametersToEstimate( parameterNames, bodies );

        // Print identifiers and indices of parameters to terminal.
        printEstimatableParameterEntries( parametersToEstimate );

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             PROPAGATE ORBIT AND VARIATIONAL EQUATIONS         /////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


        // Create simulation object and propagate dynamics.
        SingleArcVariationalEquationsSolver< > variationalEquationsSimulator(
                    bodies, integratorSettings, propagatorSettings, parametersToEstimate, true,
                    std::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ), false, true );


        // Retrieve numerical solutions for state and dependent variables
        std::map< double, Eigen::MatrixXd > numericalSolution =
                variationalEquationsSimulator.getStateTransitionMatrixSolution( );
        std::map< double, Eigen::VectorXd > dependentVariableSolution =
                variationalEquationsSimulator.getDynamicsSimulator( )->getDependentVariableHistory( );

        // Iterate over results for dependent variables, and check against manually retrieved values.
        Eigen::Vector6d currentStateDerivative;

        auto variableIteratorBack = dependentVariableSolution.begin( );
        auto variableIteratorMid =  dependentVariableSolution.begin( );
        std::advance( variableIteratorMid, 1 );
        auto variableIteratorForward = dependentVariableSolution.begin( );
        std::advance( variableIteratorForward, 2 );

        for( unsigned int i = 0; i < dependentVariableSolution.size( ) - 2; i++ )
        {

            Eigen::MatrixXd currentPartial;
            if( test < 3 )
            {
                getOutputVectorInMatrixRepresentation( variableIteratorMid->second.segment( 0, 18 ), currentPartial, 3, 6 );
            }
            else
            {
                currentPartial.setZero( 3, 6 );
                for( int j = 0; j < 3; j++ )
                {
                    Eigen::MatrixXd currentPartialComponent;
                    getOutputVectorInMatrixRepresentation(
                                variableIteratorMid->second.segment( j * 18, 18 ), currentPartialComponent, 3, 6 );
                    currentPartial += currentPartialComponent;
                }
            }

            Eigen::MatrixXd currentNumericalPartial =
                    ( numericalSolution.at( variableIteratorForward->first ) - numericalSolution.at( variableIteratorBack->first ) ) /
                    ( variableIteratorForward->first - variableIteratorBack->first );

            currentPartial *= numericalSolution.at( variableIteratorMid->first );

            for( int j = 0; j < 3; j++ )
            {
                for( int k = 0; k < 3; k++ )
                {
                    BOOST_CHECK_SMALL( std::fabs( currentNumericalPartial( j + 3, k ) - currentPartial( j, k ) ), 1.0E-12 );
                    BOOST_CHECK_SMALL( std::fabs( currentNumericalPartial( j + 3, k + 3 ) - currentPartial( j, k + 3 ) ), 1.0E-9 );

                }
            }

            if( test == 0 )
            {
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                            variableIteratorMid->second.segment( 0, 18 ), (-variableIteratorMid->second.segment( 18, 18 ) ),
                            std::numeric_limits< double >::epsilon( ) );
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                            variableIteratorMid->second.segment( 36, 18 ), ( Eigen::VectorXd::Zero( 18 ) ),
                            std::numeric_limits< double >::epsilon( ) );
                // Check consistency of partial of total acceleration of Vehicle w.r.t. Vehicle's translational state.
                Eigen::VectorXd computedTotalAccelerationPartials = variableIteratorBack->second.segment( 0, 18 );
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                            computedTotalAccelerationPartials, variableIteratorBack->second.segment( 54, 18 ),
                            std::numeric_limits< double >::epsilon( ) );
            }
            else if ( test == 1 || test == 2 )
            {
                // Check consistency of partial of total acceleration of Vehicle w.r.t. Vehicle's translational state.
                Eigen::VectorXd computedTotalAccelerationPartials = variableIteratorBack->second.segment( 0, 18 );
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                            computedTotalAccelerationPartials, variableIteratorBack->second.segment( 36, 18 ),
                            std::numeric_limits< double >::epsilon( ) );

                // Check consistency of partial of total acceleration of Vehicle w.r.t. Moon's translational state.
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                            variableIteratorBack->second.segment( 18, 18 ), variableIteratorBack->second.segment( 54, 18 ),
                            std::numeric_limits< double >::epsilon( ) );
            }
            else if ( test == 3 || test == 4 )
            {
                // Check consistency of partial of total acceleration of Vehicle w.r.t. Vehicle's translational state.
                Eigen::VectorXd computedTotalAccelerationPartials = Eigen::VectorXd::Zero( 18 );
                computedTotalAccelerationPartials += variableIteratorBack->second.segment( 0, 18 );
                computedTotalAccelerationPartials += variableIteratorBack->second.segment( 18, 18 );
                computedTotalAccelerationPartials += variableIteratorBack->second.segment( 54, 18 );
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                            computedTotalAccelerationPartials, variableIteratorBack->second.segment( 72, 18 ),
                            std::numeric_limits< double >::epsilon( ) );

                // Check consistency of partial of total acceleration of Vehicle w.r.t. Moon's translational state.
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                            variableIteratorBack->second.segment( 36, 18 ), variableIteratorBack->second.segment( 90, 18 ),
                            std::numeric_limits< double >::epsilon( ) );
            }









            variableIteratorBack++;
            variableIteratorMid++;
            variableIteratorForward++;
        }
    }
}


// Check if gravitational potential and laplacian are being saved correctly for spherical harmonics and polyhedron models
BOOST_AUTO_TEST_CASE( test_GravitationalPotentialAndLaplacianSaving )
{
    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    for ( unsigned int gravityModelsId: {0, 1} )
    {
        // Create body objects.
        std::vector< std::string > bodiesToCreate;
        bodiesToCreate.push_back( "Earth" );
        bodiesToCreate.push_back( "Moon" );
        BodyListSettings bodySettings = getDefaultBodySettings( bodiesToCreate, "Earth", "ECLIPJ2000" );

        // Use polyhedron
        if ( gravityModelsId == 1 )
        {
            // Define cuboid polyhedron dimensions
            const double w = 3000e3; // width
            const double h = 3000e3; // height
            const double l = 3000e3; // length

            // Define cuboid
            Eigen::MatrixXd verticesCoordinates(8,3);
            verticesCoordinates <<
                0.0, 0.0, 0.0,
                l, 0.0, 0.0,
                0.0, w, 0.0,
                l, w, 0.0,
                0.0, 0.0, h,
                l, 0.0, h,
                0.0, w, h,
                l, w, h;
            Eigen::MatrixXi verticesDefiningEachFacet(12,3);
            verticesDefiningEachFacet <<
                2, 1, 0,
                1, 2, 3,
                4, 2, 0,
                2, 4, 6,
                1, 4, 0,
                4, 1, 5,
                6, 5, 7,
                5, 6, 4,
                3, 6, 7,
                6, 3, 2,
                5, 3, 7,
                3, 5, 1;

            bodySettings.at( "Earth" )->gravityFieldSettings = polyhedronGravitySettingsFromMu(
                celestial_body_constants::EARTH_GRAVITATIONAL_PARAMETER, verticesCoordinates,
                verticesDefiningEachFacet, "IAU_Earth");
            bodySettings.at( "Moon" )->gravityFieldSettings = polyhedronGravitySettingsFromMu(
                celestial_body_constants::MOON_GRAVITATIONAL_PARAMETER, verticesCoordinates,
                verticesDefiningEachFacet, "IAU_Moon");
        }

        // Create Body objects
        SystemOfBodies bodies = createSystemOfBodies( bodySettings );
        bodies.createEmptyBody( "Asterix" );

        // Define propagator settings variables.
        SelectedAccelerationMap accelerationMap;
        std::vector< std::string > bodiesToPropagate;
        std::vector< std::string > centralBodies;

        bodiesToPropagate.push_back( "Asterix" );
        centralBodies.push_back( "Earth" );

        // Define propagation settings.
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfAsterix;
        if ( gravityModelsId == 0 )
        {
            accelerationsOfAsterix[ "Earth" ].push_back(
                    std::make_shared< SphericalHarmonicAccelerationSettings >( 6, 6 ) );
            accelerationsOfAsterix[ "Moon" ].push_back(
                    std::make_shared< SphericalHarmonicAccelerationSettings >( 2, 2 ) );
        }
        else
        {
            accelerationsOfAsterix[ "Earth" ].push_back( polyhedronAcceleration( ) );
            accelerationsOfAsterix[ "Moon" ].push_back( polyhedronAcceleration( ) );
        }
        accelerationMap[ "Asterix" ] = accelerationsOfAsterix;

        // Create acceleration models and propagation settings.
        basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodies, accelerationMap, bodiesToPropagate, centralBodies );

        // Set Keplerian elements for Asterix.
        Eigen::Vector6d asterixInitialStateInKeplerianElements;
        asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7500.0E3;
        asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.1;
        asterixInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 85.3 );
        asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
                = convertDegreesToRadians( 235.7 );
        asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
                = convertDegreesToRadians( 23.4 );
        asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 139.87 );

        // Convert Asterix state from Keplerian elements to Cartesian elements.
        double earthGravitationalParameter = bodies.at(
                "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
        Eigen::VectorXd systemInitialState = convertKeplerianToCartesianElements(
                asterixInitialStateInKeplerianElements,
                earthGravitationalParameter );

        double simulationEndEpoch = 10.0;
        std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >
                        ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState,
                          simulationEndEpoch, cowell );

        std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;

        dependentVariables.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                        gravity_field_potential_dependent_variable, "Asterix", "Earth" ) );
        dependentVariables.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                        gravity_field_potential_dependent_variable, "Asterix", "Moon" ) );
        dependentVariables.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                        body_fixed_relative_cartesian_position, "Asterix", "Earth" ) );
        dependentVariables.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                        body_fixed_relative_cartesian_position, "Asterix", "Moon" ) );
//        dependentVariables.push_back(
//                std::make_shared< SingleDependentVariableSaveSettings >(
//                        gravity_field_laplacian_of_potential_dependent_variable, "Asterix", "Earth" ) );
        if ( gravityModelsId == 1 )
        {
            dependentVariables.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                        gravity_field_laplacian_of_potential_dependent_variable, "Asterix", "Earth" ) );
            dependentVariables.push_back(
                    std::make_shared< SingleDependentVariableSaveSettings >(
                            gravity_field_laplacian_of_potential_dependent_variable, "Asterix", "Moon" ) );
        }

        addDepedentVariableSettings< double, double >( dependentVariables, propagatorSettings );

        // Create numerical integrator.
        double simulationStartEpoch = 0.0;
        const double fixedStepSize = 10.0;
        std::shared_ptr< IntegratorSettings< > > integratorSettings =
                std::make_shared< IntegratorSettings< > >
                        ( rungeKutta4, simulationStartEpoch, fixedStepSize );

        // Create simulation object and propagate dynamics.
        SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodies, integratorSettings, propagatorSettings );

        std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
        std::map< double, Eigen::VectorXd > depdendentVariableResult = dynamicsSimulator.getDependentVariableHistory( );

        // Get gravity models
        std::shared_ptr< gravitation::GravityFieldModel > earthGravityModel =
                bodies.at( "Earth" )->getGravityFieldModel( );
        std::shared_ptr< gravitation::GravityFieldModel > moonGravityModel =
                bodies.at( "Moon" )->getGravityFieldModel( );

        for ( std::map< double, Eigen::VectorXd >::iterator variableIterator = depdendentVariableResult.begin( );
              variableIterator != depdendentVariableResult.end( ); variableIterator++ )
        {
            double earthGravitationalPotential = variableIterator->second( 0 );
            double moonGravitationalPotential = variableIterator->second( 1 );
            Eigen::Vector3d earthBodyFixedCartesianPosition = variableIterator->second.segment( 2, 3 );
            Eigen::Vector3d moonBodyFixedCartesianPosition = variableIterator->second.segment( 5, 3 );

            // Spherical harmonics
            if ( gravityModelsId == 0 )
            {
                // Check central body gravitational potential - not sure why the tolerance needs to be so large
                /*BOOST_CHECK_CLOSE_FRACTION(
                        earthGravitationalPotential,
                        earthGravityModel->getGravitationalPotential( earthBodyFixedCartesianPosition ),
                        1e-7 );*/
                BOOST_CHECK( std::fabs( earthGravitationalPotential - earthGravityModel->getGravitationalPotential( earthBodyFixedCartesianPosition ) )
                            < std::fabs( std::min( earthGravitationalPotential, earthGravityModel->getGravitationalPotential( earthBodyFixedCartesianPosition ) ) * 1e-7 ) );

                // Check 3rd body gravitational potential - not sure why the tolerance needs to be so large
                /*BOOST_CHECK_CLOSE_FRACTION(
                        moonGravitationalPotential,
                        moonGravityModel->getGravitationalPotential( moonBodyFixedCartesianPosition ),
                        1e-11 );*/
                BOOST_CHECK( std::fabs( moonGravitationalPotential - moonGravityModel->getGravitationalPotential( moonBodyFixedCartesianPosition ) )
                            < std::fabs( std::min( moonGravitationalPotential, moonGravityModel->getGravitationalPotential( moonBodyFixedCartesianPosition ) ) * 1e-11 ) );
            }
            // Polyhedron
            else
            {
                // Check central body gravitational potential
                /*BOOST_CHECK_CLOSE_FRACTION(
                        earthGravitationalPotential,
                        earthGravityModel->getGravitationalPotential( earthBodyFixedCartesianPosition ),
                        1e-15 );*/
                BOOST_CHECK( std::fabs( earthGravitationalPotential - earthGravityModel->getGravitationalPotential( earthBodyFixedCartesianPosition ) )
                            < std::fabs( std::min( earthGravitationalPotential, earthGravityModel->getGravitationalPotential( earthBodyFixedCartesianPosition ) ) * 1e-15 ) );

                // Check 3rd body gravitational potential
                /*BOOST_CHECK_CLOSE_FRACTION(
                        moonGravitationalPotential,
                        moonGravityModel->getGravitationalPotential( moonBodyFixedCartesianPosition ),
                        1e-15 );*/
                BOOST_CHECK( std::fabs( moonGravitationalPotential - moonGravityModel->getGravitationalPotential( moonBodyFixedCartesianPosition ) )
                            < std::fabs( std::min( moonGravitationalPotential, moonGravityModel->getGravitationalPotential( moonBodyFixedCartesianPosition ) ) * 5.0e-12 ) );

                double earthGravitationalLaplacianOfPotential = variableIterator->second( 8 );
                double moonGravitationalLaplacianOfPotential = variableIterator->second( 9 );

                // Check central body gravitational potential: adding 1 because value is very close to 0
                /*BOOST_CHECK_CLOSE_FRACTION(
                        earthGravitationalLaplacianOfPotential + 1,
                        earthGravityModel->getLaplacianOfPotential( earthBodyFixedCartesianPosition ) + 1,
                        1e-15 );*/
                BOOST_CHECK_SMALL( std::fabs( earthGravitationalLaplacianOfPotential ), 1.0E-21 );

                // Check 3rd body gravitational potential: adding 1 because value is very close to 0
                /*BOOST_CHECK_CLOSE_FRACTION(
                        moonGravitationalLaplacianOfPotential + 1,
                        moonGravityModel->getLaplacianOfPotential( moonBodyFixedCartesianPosition ) + 1,
                        1e-15 );*/
                BOOST_CHECK_SMALL( std::fabs( moonGravitationalLaplacianOfPotential ), 1.0E-21 );
            }
        }
    }
}

std::pair< int, double > getClosestSatelliteDistance(
        const SystemOfBodies& bodies,
        const std::string body,
        const std::vector< std::string >& bodyList,
        const double time )
{
    std::vector< double > distances;
    distances.resize( bodyList.size( ) );
    for( unsigned int i = 0; i < distances.size( ); i++ )
    {
        distances[ i ] = ( bodies.at( body )->getEphemeris( )->getCartesianPosition( time ) -
                bodies.at( bodyList.at( i ) )->getEphemeris( )->getCartesianPosition( time ) ).norm( );
    }
    int minimumDistanceIndex = std::distance(std::begin(distances), std::min_element(std::begin(distances), std::end(distances)));

    return std::make_pair( minimumDistanceIndex, distances[ minimumDistanceIndex ] );
}

std::tuple< int, double, double > getClosestStationSatelliteDistance(
        const SystemOfBodies& bodies,
        const std::shared_ptr< ground_stations::GroundStation > groundStation,
        const std::vector< std::string >& bodyList,
        const double time )
{
    std::vector< double > distances;
    std::vector< double > elevationAngles;
    std::vector< int > indices;


    Eigen::Vector6d stationState = getLinkEndCompleteEphemerisFunction(
            observation_models::LinkEndId( "Earth", groundStation->getStationId( ) ), bodies )( time );

    for( unsigned int i = 0; i < bodyList.size( ); i++ )
    {
        Eigen::Vector3d relativePosition =
                bodies.at( bodyList.at( i ) )->getEphemeris( )->getCartesianPosition( time ) -
                stationState.segment( 0, 3 );
        double elevationAngle = groundStation->getPointingAnglesCalculator( )->calculateElevationAngleFromInertialVector(
                    relativePosition, time );
        if( elevationAngle > 0.0 )
        {
            elevationAngles.push_back( elevationAngle );
            distances.push_back( relativePosition.norm( ) );
            indices.push_back( i );
        }
    }

    if( elevationAngles.size( ) > 0 )
    {
        int minimumDistanceIndex = std::distance(std::begin(distances), std::min_element(std::begin(distances), std::end(distances)));
        return std::make_tuple( indices.at( minimumDistanceIndex ), distances.at( minimumDistanceIndex ), elevationAngles.at( minimumDistanceIndex ) );
    }
    else
    {
        return std::make_tuple( -1, TUDAT_NAN, TUDAT_NAN );
    }

}

BOOST_AUTO_TEST_CASE( test_ConstellationVariables )
{
    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation end epoch.
    const double simulationStartEpoch = 0.0 * tudat::physical_constants::JULIAN_DAY;
    const double simulationEndEpoch = 7.0 * tudat::physical_constants::JULIAN_DAY;

    // Set numerical integration fixed step size.
    const double fixedStepSize = 120.0;

    // Define body settings for simulation.
    BodyListSettings bodySettings =getDefaultBodySettings(
        {"Earth"}, "Earth", "J2000" );

    // Create Earth object
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // Create spacecraft object.
    bodies.createEmptyBody( "Satellite1" );
    bodies.createEmptyBody( "Satellite2" );
    bodies.createEmptyBody( "Satellite3" );
    bodies.createEmptyBody( "Satellite4" );
    bodies.createEmptyBody( "Satellite5" );

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Define propagation settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfSatellite1;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfSatellite2;
    accelerationsOfSatellite1[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >(
                                                     2, 2 ) );
//    accelerationsOfSatellite1[ "Earth" ].push_back( std::make_shared< AccelerationSettings >(
//                                                     aerodynamic ) );
    accelerationsOfSatellite2[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >(
                                                     2, 2 ) );
    accelerationMap[ "Satellite1" ] = accelerationsOfSatellite1;
    accelerationMap[ "Satellite2" ] = accelerationsOfSatellite1;

    bodiesToPropagate.push_back( "Satellite1" );
    bodiesToPropagate.push_back( "Satellite2" );

    centralBodies.push_back( "Earth" );
    centralBodies.push_back( "Earth" );

    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodies, accelerationMap, bodiesToPropagate, centralBodies );

    double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );

    // Set Keplerian elements for Asterix.
    Eigen::Vector6d satellite1InitialState;
    satellite1InitialState( semiMajorAxisIndex ) = 6800.0E3;
    satellite1InitialState( eccentricityIndex ) = 0.1;
    satellite1InitialState( inclinationIndex ) = convertDegreesToRadians( 85.3 );
    satellite1InitialState( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 235.7 );
    satellite1InitialState( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 23.4 );
    satellite1InitialState( trueAnomalyIndex ) = convertDegreesToRadians( 139.87 );

    Eigen::Vector6d satellite2InitialState = satellite1InitialState;
    satellite2InitialState( eccentricityIndex ) = 0.0;
    satellite2InitialState( inclinationIndex ) = convertDegreesToRadians( 86.3 );
    satellite2InitialState( trueAnomalyIndex ) = convertDegreesToRadians( 141.87 );

    Eigen::Vector6d satellite3InitialState = satellite1InitialState;
    satellite3InitialState( eccentricityIndex ) = 0.0;
    satellite3InitialState( inclinationIndex ) = convertDegreesToRadians( 86.3 );
    satellite3InitialState( trueAnomalyIndex ) = convertDegreesToRadians( 137.87 );
    bodies.at( "Satellite3" )->setEphemeris( createBodyEphemeris( std::make_shared< KeplerEphemerisSettings >(
                                                 satellite3InitialState, 0.0,earthGravitationalParameter, "Earth", "J2000" ),
                                             "Satellite3" ) );

    Eigen::Vector6d satellite4InitialState = satellite1InitialState;
    satellite4InitialState( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 24.4 );
    bodies.at( "Satellite4" )->setEphemeris( createBodyEphemeris( std::make_shared< KeplerEphemerisSettings >(
                                                 satellite4InitialState, 0.0,earthGravitationalParameter, "Earth", "J2000" ),
                                             "Satellite4" ) );

    Eigen::Vector6d satellite5InitialState = satellite1InitialState;
    satellite5InitialState( semiMajorAxisIndex ) = 6810.0E3;
    bodies.at( "Satellite5" )->setEphemeris( createBodyEphemeris( std::make_shared< KeplerEphemerisSettings >(
                                                 satellite5InitialState, 0.0, earthGravitationalParameter, "Earth", "J2000" ),
                                             "Satellite5" ) );

    createGroundStation( bodies.at( "Earth" ), "Station", ( Eigen::Vector3d( ) << 0.0, 1.5, 2.1 ).finished( ),
                         coordinate_conversions::geodetic_position );

    // Convert Asterix state from Keplerian elements to Cartesian elements.
    Eigen::VectorXd systemInitialState = Eigen::VectorXd::Zero( 12 );
    systemInitialState.segment( 0, 6 ) = convertKeplerianToCartesianElements(
                satellite1InitialState,
                earthGravitationalParameter );
    systemInitialState.segment( 6, 6 ) = convertKeplerianToCartesianElements(
                satellite2InitialState,
                earthGravitationalParameter );

    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4, simulationStartEpoch, fixedStepSize );


    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;

    dependentVariables.push_back(
                std::make_shared< MinimumConstellationDistanceDependentVariableSaveSettings >(
                    "Satellite1", std::vector< std::string >( { "Satellite2", "Satellite3", "Satellite4", "Satellite5" } ) ) );
    dependentVariables.push_back(
                std::make_shared< MinimumConstellationDistanceDependentVariableSaveSettings >(
                    "Satellite2", std::vector< std::string >( { "Satellite1", "Satellite3", "Satellite4", "Satellite5" } ) ) );
    dependentVariables.push_back(
                std::make_shared< MinimumConstellationDistanceDependentVariableSaveSettings >(
                    "Satellite3", std::vector< std::string >( { "Satellite1", "Satellite2", "Satellite4", "Satellite5" } ) ) );
    dependentVariables.push_back(
                std::make_shared< MinimumConstellationStationDistanceDependentVariableSaveSettings >(
                    "Earth", "Station", std::vector< std::string >( { "Satellite1", "Satellite2", "Satellite3", "Satellite4", "Satellite5" } ), 0.0 ) );

    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, simulationEndEpoch, cowell, dependentVariables );


    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodies, integratorSettings, propagatorSettings, true, false, true );
    std::map< double, Eigen::VectorXd > dependentVariableResults = dynamicsSimulator.getDependentVariableHistory( );
    std::pair< int, double > testPair;
    std::tuple< int, double, double > testTuple;
    for( auto it : dependentVariableResults )
    {
        testPair = getClosestSatelliteDistance(
                bodies, "Satellite1", std::vector< std::string >( { "Satellite2", "Satellite3", "Satellite4", "Satellite5" } ), it.first );
        BOOST_CHECK_CLOSE_FRACTION( it.second( 0 ), testPair.second, 1.0E-12 );
        BOOST_CHECK_EQUAL( it.second( 1 ), testPair.first );

        testPair = getClosestSatelliteDistance(
                bodies, "Satellite2", std::vector< std::string >( { "Satellite1", "Satellite3", "Satellite4", "Satellite5" } ), it.first );
        BOOST_CHECK_CLOSE_FRACTION( it.second( 2 ), testPair.second, 1.0E-12 );
        BOOST_CHECK_EQUAL( it.second( 3 ), testPair.first );

        testPair = getClosestSatelliteDistance(
                bodies, "Satellite3", std::vector< std::string >( { "Satellite1", "Satellite2", "Satellite4", "Satellite5" } ), it.first );
        BOOST_CHECK_CLOSE_FRACTION( it.second( 4 ), testPair.second, 1.0E-12 );
        BOOST_CHECK_EQUAL( it.second( 5 ), testPair.first );

        testTuple = getClosestStationSatelliteDistance(
                bodies, bodies.at( "Earth" )->getGroundStation( "Station" ),
                     std::vector< std::string >( { "Satellite1", "Satellite2", "Satellite3", "Satellite4", "Satellite5" } ),
                    it.first );
        if( std::get< 0 >( testTuple ) == -1 )
        {
            BOOST_CHECK_EQUAL( it.second( 7 ), -1 );
            BOOST_CHECK_EQUAL( ( it.second( 6 ) != it.second( 6 ) ), true );
            BOOST_CHECK_EQUAL( ( it.second( 8 ) != it.second( 8 ) ), true );
        }
        else
        {
            BOOST_CHECK_EQUAL( it.second( 7 ), std::get< 0 >( testTuple ) );
            BOOST_CHECK_CLOSE_FRACTION( it.second( 6 ), std::get< 1 >( testTuple ), 1.0E-12 );
            BOOST_CHECK_CLOSE_FRACTION( it.second( 8 ), std::get< 2 >( testTuple ), 1.0E-12 );
        }

    }

}

BOOST_AUTO_TEST_SUITE_END( )

}

}
