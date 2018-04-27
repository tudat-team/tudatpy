/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>
#include <boost/bind.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include "Tudat/Astrodynamics/Aerodynamics/UnitTests/testApolloCapsuleCoefficients.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/geodeticCoordinateConversions.h"
#include "Tudat/Basics/testMacros.h"
#include "Tudat/SimulationSetup/PropagationSetup/dynamicsSimulator.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/SimulationSetup/PropagationSetup/createNumericalSimulator.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.h"
#include "Tudat/InputOutput/basicInputOutput.h"
#include <limits>
#include <string>

#include <Eigen/Core>

namespace tudat
{

namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_dependent_variable_output )

//! Propagate entry of Apollo capsule, and save a list of dependent variables during entry. The saved dependent variables
//! are compared against theoretical/manual values in this test.
BOOST_AUTO_TEST_CASE( testDependentVariableOutput )
{
    using namespace tudat;
    using namespace ephemerides;
    using namespace interpolators;
    using namespace numerical_integrators;
    using namespace spice_interface;
    using namespace simulation_setup;
    using namespace basic_astrodynamics;
    using namespace orbital_element_conversions;
    using namespace propagators;
    using namespace aerodynamics;
    using namespace basic_mathematics;
    using namespace input_output;

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

    for( unsigned int testCase = 0; testCase < 4; testCase++ )
    {
        // Define simulation body settings.
        std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
                getDefaultBodySettings( { "Earth", "Moon" }, simulationStartEpoch - 10.0 * fixedStepSize,
                                        simulationEndEpoch + 10.0 * fixedStepSize );
        bodySettings[ "Earth" ]->gravityFieldSettings =
                boost::make_shared< simulation_setup::GravityFieldSettings >( central_spice );

        if( testCase >= 2 )
        {
            bodySettings[ "Earth" ]->atmosphereSettings =
                    boost::make_shared< simulation_setup::AtmosphereSettings >( nrlmsise00 );
        }

        bool isOblateSpheroidUsed = 0;
        double oblateSpheroidEquatorialRadius = 6378.0E3;
        double oblateSpheroidFlattening = 1.0 / 300.0;

        if( testCase% 2 == 0 )
        {
            isOblateSpheroidUsed = 1;
            bodySettings[ "Earth" ]->shapeModelSettings =
                    boost::make_shared< simulation_setup::OblateSphericalBodyShapeSettings >(
                        oblateSpheroidEquatorialRadius, oblateSpheroidFlattening );
        }

        // Create Earth object
        simulation_setup::NamedBodyMap bodyMap = simulation_setup::createBodies( bodySettings );

        // Create vehicle objects.
        bodyMap[ "Apollo" ] = boost::make_shared< simulation_setup::Body >( );

        // Create vehicle aerodynamic coefficients
        bodyMap[ "Apollo" ]->setAerodynamicCoefficientInterface(
                    unit_tests::getApolloCoefficientInterface( ) );
        bodyMap[ "Apollo" ]->setConstantBodyMass( 5.0E3 );
        bodyMap[ "Apollo" ]->setEphemeris(
                    boost::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                        boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector6d  > >( ),
                        "Earth" ) );
        boost::shared_ptr< system_models::VehicleSystems > vehicleSystems =
                boost::make_shared< system_models::VehicleSystems >( );

        double noseRadius = 0.7;
        double wallEmissivity = 0.7;
        vehicleSystems->setNoseRadius( noseRadius );
        vehicleSystems->setWallEmissivity( wallEmissivity );

        bodyMap[ "Apollo" ]->setVehicleSystems( vehicleSystems );


        // Finalize body creation.
        setGlobalFrameBodyEphemerides( bodyMap, "Earth", "ECLIPJ2000" );

        // Define propagator settings variables.
        SelectedAccelerationMap accelerationMap;
        std::vector< std::string > bodiesToPropagate;
        std::vector< std::string > centralBodies;

        // Define acceleration model settings.
        std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfApollo;
        accelerationsOfApollo[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
        accelerationsOfApollo[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( aerodynamic ) );
        accelerationsOfApollo[ "Moon" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
        accelerationMap[ "Apollo" ] = accelerationsOfApollo;

        bodiesToPropagate.push_back( "Apollo" );
        centralBodies.push_back( "Earth" );

        // Set initial state
        Eigen::Vector6d systemInitialState = apolloInitialState;

        // Define list of dependent variables to save.
        std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
        dependentVariables.push_back(
                    boost::make_shared< SingleDependentVariableSaveSettings >( mach_number_dependent_variable, "Apollo" ) );
        dependentVariables.push_back(
                    boost::make_shared< SingleDependentVariableSaveSettings >( altitude_dependent_variable,
                                                                               "Apollo", "Earth" ) );
        dependentVariables.push_back(
                    boost::make_shared< SingleDependentVariableSaveSettings >( relative_distance_dependent_variable,
                                                                               "Apollo", "Earth" ) );
        dependentVariables.push_back(
                    boost::make_shared< SingleDependentVariableSaveSettings >( relative_speed_dependent_variable,
                                                                               "Apollo", "Earth" ) );
        dependentVariables.push_back(
                    boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                        central_gravity, "Apollo", "Earth", 1 ) );

        dependentVariables.push_back(
                    boost::make_shared< SingleDependentVariableSaveSettings >( total_aerodynamic_g_load_variable,
                                                                               "Apollo", "Earth" ) );
        dependentVariables.push_back(
                    boost::make_shared< SingleDependentVariableSaveSettings >( stagnation_point_heat_flux_dependent_variable,
                                                                               "Apollo", "Earth" ) );

        dependentVariables.push_back(
                    boost::make_shared< SingleDependentVariableSaveSettings >( local_temperature_dependent_variable,
                                                                               "Apollo", "Earth" ) );
        dependentVariables.push_back(
                    boost::make_shared< SingleDependentVariableSaveSettings >( geodetic_latitude_dependent_variable,
                                                                               "Apollo", "Earth" ) );

        dependentVariables.push_back(
                    boost::make_shared< SingleDependentVariableSaveSettings >( local_density_dependent_variable,
                                                                               "Apollo", "Earth" ) );

        dependentVariables.push_back(
                    boost::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                        "Apollo", reference_frames::latitude_angle ) );
        dependentVariables.push_back(
                    boost::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                        "Apollo", reference_frames::longitude_angle ) );

        dependentVariables.push_back(
                    boost::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                        "Apollo", reference_frames::angle_of_attack ) );
        dependentVariables.push_back(
                    boost::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                        "Apollo", reference_frames::angle_of_sideslip ) );
        dependentVariables.push_back(
                    boost::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                        "Apollo", reference_frames::bank_angle ) );




        dependentVariables.push_back(
                    boost::make_shared< SingleDependentVariableSaveSettings >( relative_position_dependent_variable,
                                                                               "Apollo", "Earth" ) );
        dependentVariables.push_back(
                    boost::make_shared< SingleDependentVariableSaveSettings >( relative_velocity_dependent_variable,
                                                                               "Apollo", "Earth" ) );
        dependentVariables.push_back(
                    boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                        central_gravity, "Apollo", "Earth", 0 ) );
        dependentVariables.push_back(
                    boost::make_shared< SingleDependentVariableSaveSettings >(
                        total_acceleration_dependent_variable, "Apollo" ) );
        dependentVariables.push_back(
                    boost::make_shared< SingleDependentVariableSaveSettings >(
                        aerodynamic_moment_coefficients_dependent_variable, "Apollo" ) );
        dependentVariables.push_back(
                    boost::make_shared< SingleDependentVariableSaveSettings >(
                        aerodynamic_force_coefficients_dependent_variable, "Apollo" ) );
        dependentVariables.push_back(
                    boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                        aerodynamic, "Apollo", "Earth", 0 ) );
        dependentVariables.push_back(
                    boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                        central_gravity, "Apollo", "Moon", 0 ) );
        dependentVariables.push_back(
                    boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                        third_body_central_gravity, "Apollo", "Moon", 0 ) );
        dependentVariables.push_back(
                    boost::make_shared< SingleDependentVariableSaveSettings >(
                        keplerian_state_dependent_variable,  "Apollo", "Earth" ) );
        dependentVariables.push_back(
                    boost::make_shared< SingleDependentVariableSaveSettings >(
                        modified_equinocial_state_dependent_variable,  "Apollo", "Earth" ) );


        // Create acceleration models and propagation settings.
        basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

        setTrimmedConditions( bodyMap.at( "Apollo" ) );

        boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                boost::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState,
                  boost::make_shared< propagators::PropagationTimeTerminationSettings >( 3200.0 ), cowell,
                  boost::make_shared< DependentVariableSaveSettings >( dependentVariables ) );
        boost::shared_ptr< IntegratorSettings< > > integratorSettings =
                boost::make_shared< IntegratorSettings< > >
                ( rungeKutta4, simulationStartEpoch, fixedStepSize );

        // Create simulation object and propagate dynamics.
        SingleArcDynamicsSimulator< > dynamicsSimulator(
                    bodyMap, integratorSettings, propagatorSettings, true, false, false );

        // Retrieve numerical solutions for state and dependent variables
        std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > numericalSolution =
                dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
        std::map< double, Eigen::VectorXd > dependentVariableSolution =
                dynamicsSimulator.getDependentVariableHistory( );

        // Iterate over results for dependent variables, and check against manually retrieved values.
        Eigen::Vector6d currentStateDerivative;
        Eigen::Vector3d manualCentralGravity;
        boost::shared_ptr< ephemerides::RotationalEphemeris > earthRotationModel =
                bodyMap.at( "Earth" )->getRotationalEphemeris( );
        boost::shared_ptr< aerodynamics::AtmosphereModel > earthAtmosphereModel =
                bodyMap.at( "Earth" )->getAtmosphereModel( );
        boost::shared_ptr< aerodynamics::FlightConditions > apolloFlightConditions =
                bodyMap.at( "Apollo" )->getFlightConditions( );
        boost::shared_ptr< aerodynamics::AerodynamicCoefficientInterface > apolloCoefficientInterface =
                bodyMap.at( "Apollo" )->getAerodynamicCoefficientInterface( );

        for( std::map< double, Eigen::VectorXd >::iterator variableIterator = dependentVariableSolution.begin( );
             variableIterator != dependentVariableSolution.end( ); variableIterator++ )
        {
            double machNumber = variableIterator->second( 0 );
            double altitude = variableIterator->second( 1 );
            double relativeDistance = variableIterator->second( 2 );
            double relativeSpeed= variableIterator->second( 3 );
            double gravitationalAccelerationNorm = variableIterator->second( 4 );
            double gLoad = variableIterator->second( 5 );
            double stagnationPointHeatFlux = variableIterator->second( 6 );
            double freestreamTemperature = variableIterator->second( 7 );
            double geodeticLatitude = variableIterator->second( 8 );
            double freestreamDensity = variableIterator->second( 9 );
            double latitude = variableIterator->second( 10 );
            double longitude = variableIterator->second( 11 );
            double angleOfAttack = variableIterator->second( 12 );
            double sideslipAngle = variableIterator->second( 13 );
            double bankAngle = variableIterator->second( 14 );

            Eigen::Vector3d relativePosition = variableIterator->second.segment( 15, 3 );
            Eigen::Vector3d computedBodyFixedPosition =
                    earthRotationModel->getRotationToTargetFrame( variableIterator->first ) * relativePosition;
            Eigen::Vector3d computedSphericalBodyFixedPosition =
                    coordinate_conversions::convertCartesianToSpherical( computedBodyFixedPosition );

            Eigen::Vector3d relativeVelocity = variableIterator->second.segment( 18, 3 );
            Eigen::Vector3d gravitationalAcceleration = variableIterator->second.segment( 21, 3 );
            Eigen::Vector3d totalAcceleration = variableIterator->second.segment( 24, 3 );
            Eigen::Vector3d momentCoefficients = variableIterator->second.segment( 27, 3 );
            Eigen::Vector3d forceCoefficients = variableIterator->second.segment( 30, 3 );
            Eigen::Vector3d aerodynamicAcceleration = variableIterator->second.segment( 33, 3 );
            Eigen::Vector3d moonAcceleration1 = variableIterator->second.segment( 36, 3 );
            Eigen::Vector3d moonAcceleration2 = variableIterator->second.segment( 39, 3 );

            Eigen::Vector6d keplerElements =  variableIterator->second.segment( 42, 6 );
            Eigen::Vector6d modifiedEquinoctialElements =  variableIterator->second.segment( 48, 6 );

            currentStateDerivative = dynamicsSimulator.getDynamicsStateDerivative( )->computeStateDerivative(
                        variableIterator->first, numericalSolution.at( variableIterator->first ) );

            // Manually compute central gravity.
            manualCentralGravity =
                    -bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( ) *
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

            // Check total acceleration (tolerance is not epsilon due to numerical root finding for trim)
            for( unsigned int i = 0; i < 3; i++ )
            {
                BOOST_CHECK_SMALL(
                            std::fabs( currentStateDerivative( 3 + i ) - totalAcceleration( i ) ), 1.0E-13 );
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



        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )


}

}
