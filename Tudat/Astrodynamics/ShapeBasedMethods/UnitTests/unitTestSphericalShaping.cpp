/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Wakker, K. F. (2007), Lecture Notes Astrodynamics II (Chapter 18), TU Delft course AE4-874,
 *          Delft University of technology, Delft, The Netherlands.
 *
 */

#define BOOST_TEST_MAIN

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <Eigen/Dense>
#include <math.h>
#include <iostream>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"
#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"
#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositions.h"
#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/celestialBodyConstants.h"
#include "Tudat/Astrodynamics/ShapeBasedMethods/baseFunctionsSphericalShaping.h"
#include "Tudat/Astrodynamics/ShapeBasedMethods/compositeFunctionSphericalShaping.h"
#include "Tudat/Astrodynamics/ShapeBasedMethods/sphericalShaping.h"

namespace tudat
{
namespace unit_tests
{

//! Test spherical shaping implementation.
BOOST_AUTO_TEST_SUITE( test_spherical_shaping )

//! Test.
BOOST_AUTO_TEST_CASE( test_spherical_shaping1 )
{

    double numberOfRevolutions = 1.0;
    double julianDate = 8174.5 * physical_constants::JULIAN_DAY;
    double  timeOfFlight = 580.0;

    // Ephemeris departure body.
    ephemerides::EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions>(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );

    // Ephemeris arrival body.
    ephemerides::EphemerisPointer pointerToArrivalBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions >(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mars );

    std::cout << "cartesian state departure body = Earth: " << pointerToDepartureBodyEphemeris->getCartesianState( julianDate  ) << "\n\n";
    std::cout << "cartesian state arrival body = Mars: "
            << pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight * physical_constants::JULIAN_DAY ) << "\n\n";

    Eigen::VectorXd radialFunctionCoefficients = ( Eigen::Vector7d() << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 ).finished();
    Eigen::Vector4d elevationFunctionCoefficients = ( Eigen::Vector4d() << 1.0, 1.0, 1.0, 1.0 ).finished();

    Eigen::Vector6d normalisedInitialState;
    normalisedInitialState.segment( 0, 3 ) = pointerToDepartureBodyEphemeris->getCartesianState( julianDate ).segment( 0, 3 )
            / physical_constants::ASTRONOMICAL_UNIT;
    normalisedInitialState.segment( 3, 3 ) = pointerToDepartureBodyEphemeris->getCartesianState( julianDate ).segment( 3, 3 )
             * physical_constants::JULIAN_YEAR / physical_constants::ASTRONOMICAL_UNIT;

    Eigen::Vector6d normalisedFinalState;
    normalisedFinalState.segment( 0, 3 ) = pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight * tudat::physical_constants::JULIAN_DAY ).segment( 0, 3 )
            / physical_constants::ASTRONOMICAL_UNIT;
    normalisedFinalState.segment( 3, 3 ) = pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight * tudat::physical_constants::JULIAN_DAY ).segment( 3, 3 )
             * physical_constants::JULIAN_YEAR / physical_constants::ASTRONOMICAL_UNIT;

    shape_based_methods::SphericalShaping sphericalShaping = shape_based_methods::SphericalShaping(
                normalisedInitialState,
                normalisedFinalState,
                timeOfFlight * tudat::physical_constants::JULIAN_DAY  / physical_constants::JULIAN_YEAR, 0.000703 /*- 0.00000000000007*/,
                radialFunctionCoefficients, elevationFunctionCoefficients,
                celestial_body_constants::SUN_GRAVITATIONAL_PARAMETER  * std::pow( physical_constants::JULIAN_YEAR, 2.0 )
                / std::pow( physical_constants::ASTRONOMICAL_UNIT, 3.0 ),
                root_finders::bisection_root_finder, 1.0e-6, 1.0e-1, TUDAT_NAN, 30, 1.0e-6 );

    std::cout << "test initial spherical position: " << sphericalShaping.computeCurrentSphericalState( sphericalShaping.getInitialAzimuthalAngle() /*2.25776*/ ).transpose() << "\n\n";
    std::cout << "test final spherical position: " << sphericalShaping.computeCurrentSphericalState( sphericalShaping.getFinalAzimuthalAngle() /*12.8632*/ ).transpose() << "\n\n";

    std::cout << "test initial cartesian position: " << sphericalShaping.computeCurrentCartesianState( sphericalShaping.getInitialAzimuthalAngle() /*2.25776*/ ).transpose() << "\n\n";
    std::cout << "test final cartesian position: " << sphericalShaping.computeCurrentCartesianState( sphericalShaping.getFinalAzimuthalAngle() /*12.8632*/ ).transpose() << "\n\n";

    std::cout << "difference final cartesian state: " << ( sphericalShaping.computeCurrentCartesianState( sphericalShaping.getFinalAzimuthalAngle() )
                 - pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight * tudat::physical_constants::JULIAN_DAY ) ).transpose() << "\n\n";

    std::map< double, Eigen::Vector6d > mapPosition;
    std::map< double, Eigen::Vector3d > mapAcceleration;

    // Compute step size.
    double stepSize = ( sphericalShaping.getFinalAzimuthalAngle() - sphericalShaping.getInitialAzimuthalAngle() ) / 5000.0;

    // Initialise peak acceleration.
    double peakThrustAcceleration = 0.0;

    // Check that the trajectory is feasible, ie curved toward the central body.
    for ( int i = 0 ; i <= 5000 ; i++ )
    {
        double currentThetaAngle = sphericalShaping.getInitialAzimuthalAngle() + i * stepSize;

        Eigen::Vector6d currentUnnormalisedState;
        currentUnnormalisedState.segment( 0, 3 ) = sphericalShaping.computeCurrentCartesianState( currentThetaAngle ).segment( 0, 3 )
                * physical_constants::ASTRONOMICAL_UNIT;
        currentUnnormalisedState.segment( 3, 3 ) = sphericalShaping.computeCurrentCartesianState( currentThetaAngle ).segment( 3, 3 )
                * physical_constants::ASTRONOMICAL_UNIT / physical_constants::JULIAN_YEAR;

        mapPosition[ currentThetaAngle ] = currentUnnormalisedState;
        mapAcceleration[ currentThetaAngle ] = sphericalShaping.computeSphericalControlAccelerationVector( currentThetaAngle )
                * physical_constants::ASTRONOMICAL_UNIT / std::pow( physical_constants::JULIAN_YEAR, 2.0 );

        if ( sphericalShaping.computeSphericalControlAccelerationVector( currentThetaAngle ).norm() >  peakThrustAcceleration )
        {
            peakThrustAcceleration = sphericalShaping.computeSphericalControlAccelerationVector( currentThetaAngle ).norm();
        }
    }

    tudat::input_output::writeDataMapToTextFile( mapPosition,
                                          "mapPosition.dat",
                                           "C:/Users/chamb/Documents/Master_2/SOCIS/SphericalShapingTest/",
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    tudat::input_output::writeDataMapToTextFile( mapAcceleration,
                                          "mapAcceleration.dat",
                                           "C:/Users/chamb/Documents/Master_2/SOCIS/SphericalShapingTest/",
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    std::cout << "deltaV: " << sphericalShaping.computeDeltav() * physical_constants::ASTRONOMICAL_UNIT / physical_constants::JULIAN_YEAR << "\n\n";
    std::cout << "peak acceleration: " << peakThrustAcceleration * physical_constants::ASTRONOMICAL_UNIT / std::pow( physical_constants::JULIAN_YEAR, 2.0 ) << "\n\n";
    std::cout << "free coefficient: " << sphericalShaping.getRadialCompositionFunctionCoefficients().transpose() << "\n\n";



    std::map< double, Eigen::VectorXd > fullPropagationResults;
    std::map< double, Eigen::VectorXd > shapingMethodResults;
    std::map< double, Eigen::VectorXd > dependentVariablesHistory;

    spice_interface::loadStandardSpiceKernels( );

    // Create central, departure and arrival bodies.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Mars" );
    bodiesToCreate.push_back( "Jupiter" );


    std::map< std::string, std::shared_ptr< simulation_setup::BodySettings > > bodySettings =
            simulation_setup::getDefaultBodySettings( bodiesToCreate );

    std::string frameOrigin = "SSB";
    std::string frameOrientation = "ECLIPJ2000";


    // Define central body ephemeris settings.
    bodySettings[ "Sun" ]->ephemerisSettings = std::make_shared< simulation_setup::ConstantEphemerisSettings >(
                ( Eigen::Vector6d( ) << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( ), frameOrigin, frameOrientation );

    bodySettings[ "Sun" ]->ephemerisSettings->resetFrameOrientation( frameOrientation );
    bodySettings[ "Sun" ]->rotationModelSettings->resetOriginalFrame( frameOrientation );


    // Create body map.
    simulation_setup::NamedBodyMap bodyMap = createBodies( bodySettings );

    bodyMap[ "Vehicle" ] = std::make_shared< simulation_setup::Body >( );
    bodyMap.at( "Vehicle" )->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                                                         std::shared_ptr< interpolators::OneDimensionalInterpolator
                                                         < double, Eigen::Vector6d > >( ), frameOrigin, frameOrientation ) );


    setGlobalFrameBodyEphemerides( bodyMap, frameOrigin, frameOrientation );



    std::vector< std::string > bodiesToPropagate;
    bodiesToPropagate.push_back( "Vehicle" );
    std::vector< std::string > centralBodies;
    centralBodies.push_back( "Sun" );

    // Acceleration from the central body.
    std::map< std::string, std::vector< std::shared_ptr< simulation_setup::AccelerationSettings > > > bodyToPropagateAccelerations;
    bodyToPropagateAccelerations[ "Sun" ].push_back( std::make_shared< simulation_setup::AccelerationSettings >(
                                                                basic_astrodynamics::central_gravity ) );
//    bodyToPropagateAccelerations[ "Mars" ].push_back( std::make_shared< simulation_setup::AccelerationSettings >(
//                                                          basic_astrodynamics::central_gravity ) );
//    bodyToPropagateAccelerations[ "Earth" ].push_back( std::make_shared< simulation_setup::AccelerationSettings >(
//                                                          basic_astrodynamics::central_gravity ) );
//    bodyToPropagateAccelerations[ "Jupiter" ].push_back( std::make_shared< simulation_setup::AccelerationSettings >(
//                                                          basic_astrodynamics::central_gravity ) );

    simulation_setup::SelectedAccelerationMap accelerationMap;
    accelerationMap[ "Vehicle" ] = bodyToPropagateAccelerations;

    // Create the acceleration map.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    // Define integrator settings
    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
            std::make_shared< numerical_integrators::IntegratorSettings< double > > ( numerical_integrators::rungeKutta4, 0.0,
                                                                                      timeOfFlight * physical_constants::JULIAN_DAY / ( 1000.0 ) /*/ 400.0*/ );

    // Define mass function of the vehicle.
    std::function< double( const double ) > newMassFunction = [ = ]( const double currentTime )
    {
        return 200.0 - 50.0 / ( timeOfFlight * physical_constants::JULIAN_DAY ) * currentTime ;
    };
    bodyMap[ "Vehicle" ]->setBodyMassFunction( newMassFunction );


    // Define specific impulse function.
    std::function< double( const double ) > specificImpulseFunction = [ = ]( const double currentTime )
    {
        return 200.0;
    };

    // Define list of dependent variables to save.
    std::vector< std::shared_ptr< propagators::SingleDependentVariableSaveSettings > > dependentVariablesList;
    dependentVariablesList.push_back( std::make_shared< propagators::SingleAccelerationDependentVariableSaveSettings >(
                        basic_astrodynamics::thrust_acceleration, "Vehicle", "Vehicle", 0 ) );

    // Create object with list of dependent variables
    std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< propagators::DependentVariableSaveSettings >( dependentVariablesList );

    // Create termination conditions settings.
    std::pair< std::shared_ptr< propagators::PropagationTerminationSettings >, std::shared_ptr< propagators::PropagationTerminationSettings > > terminationConditions;

    terminationConditions.first = std::make_shared< propagators::PropagationTimeTerminationSettings >( 0.0, true );
    terminationConditions.second = std::make_shared< propagators::PropagationTimeTerminationSettings >( timeOfFlight * physical_constants::JULIAN_DAY, true );

    // Compute shaped trajectory and propagated trajectory.
    sphericalShaping.computeShapingTrajectoryAndFullPropagation( bodyMap, accelerationModelMap, "Sun", "Vehicle", specificImpulseFunction,
                                                                      integratorSettings, terminationConditions,
                                                                      fullPropagationResults, shapingMethodResults, dependentVariablesHistory,
                                                                      propagators::cowell, dependentVariablesToSave );

    tudat::input_output::writeDataMapToTextFile( fullPropagationResults,
                                          "fullPropagationResults.dat",
                                           "C:/Users/chamb/Documents/Master_2/SOCIS/SphericalShapingTest/",
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    tudat::input_output::writeDataMapToTextFile( shapingMethodResults,
                                          "shapingMethodResults.dat",
                                           "C:/Users/chamb/Documents/Master_2/SOCIS/SphericalShapingTest/",
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    tudat::input_output::writeDataMapToTextFile( dependentVariablesHistory,
                                          "dependentVariables.dat",
                                           "C:/Users/chamb/Documents/Master_2/SOCIS/SphericalShapingTest/",
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    std::cout << "test conversion azimuth angle -> time: " << sphericalShaping.computeCurrentTimeFromAzimuthAngle( sphericalShaping.getFinalAzimuthalAngle() )
                 * physical_constants::JULIAN_YEAR << "\n\n";


}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
