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

#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <Eigen/Dense>
#include <cmath>
#include <iostream>

#include <tudat/basics/testMacros.h>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/math/basic/coordinateConversions.h"
#include "tudat/astro/ephemerides/approximatePlanetPositions.h"
#include "tudat/simulation/simulation.h"
#include "tudat/simulation/propagation_setup/propagationLowThrustProblem.h"
#include "tudat/astro/basic_astro/celestialBodyConstants.h"
#include "tudat/astro/low_thrust/shape_based/sphericalShapingLeg.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::simulation_setup;
using namespace tudat::ephemerides;
using namespace tudat::shape_based_methods;
using namespace tudat::root_finders;
using namespace tudat::propagators;

//! Test spherical shaping implementation.
BOOST_AUTO_TEST_SUITE( test_spherical_shaping )

//! Test.
BOOST_AUTO_TEST_CASE( test_spherical_shaping_earth_mars_transfer )
{
    spice_interface::loadStandardSpiceKernels( );

    int numberOfRevolutions = 1;
    double julianDate = 8174.5 * physical_constants::JULIAN_DAY;
    double timeOfFlight = 580.0 * physical_constants::JULIAN_DAY;

    // Ephemeris departure body.
    EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ApproximateJplEphemeris>( "Earth" );
    EphemerisPointer pointerToArrivalBodyEphemeris = std::make_shared< ApproximateJplEphemeris >( "Mars" );
    Eigen::Vector6d initialState = pointerToDepartureBodyEphemeris->getCartesianState( julianDate );
    Eigen::Vector6d finalState = pointerToArrivalBodyEphemeris->getCartesianState(
                julianDate + timeOfFlight );

    // Define root finder settings (used to update the updated value of the free coefficient, so that it
    // matches the required time of flight).
    std::shared_ptr< RootFinderSettings > rootFinderSettings =
            tudat::root_finders::bisectionRootFinderSettings( 1.0E-6, TUDAT_NAN, TUDAT_NAN, 30 );

    SphericalShapingLeg *sphericalShapingLegPointer;
    for ( unsigned int creationType = 0; creationType < 2; creationType++ )
    {
        if ( creationType == 0)
        {
            sphericalShapingLegPointer = new SphericalShapingLeg(
                    pointerToDepartureBodyEphemeris, pointerToArrivalBodyEphemeris,
                    spice_interface::getBodyGravitationalParameter("Sun"),
                    rootFinderSettings, 1.0e-6, 1.0e-1);
        }
        else if ( creationType == 1)
        {
            std::function< Eigen::Vector3d( ) > departureVelocityFunction = [=]( ){ return initialState.segment( 3, 3 ); };
            std::function< Eigen::Vector3d( ) > arrivalVelocityFunction = [=]( ){ return finalState.segment( 3, 3 ); };

            sphericalShapingLegPointer = new SphericalShapingLeg(
                    pointerToDepartureBodyEphemeris, pointerToArrivalBodyEphemeris,
                    spice_interface::getBodyGravitationalParameter("Sun"),
                    departureVelocityFunction, arrivalVelocityFunction,
                    rootFinderSettings, 1.0e-6, 1.0e-1);
        }

        sphericalShapingLegPointer->updateLegParameters((
                Eigen::Vector3d( )<< julianDate, julianDate + timeOfFlight, numberOfRevolutions ).finished( ) );

        // Initialise peak acceleration
        double peakThrustAcceleration = 0.0;

        std::map< double, Eigen::Vector3d > thrustAccelerationsAlongTrajectory;
        sphericalShapingLegPointer->getThrustAccelerationsAlongTrajectory(thrustAccelerationsAlongTrajectory, 5000 );

        // Check initial and final time on output list
        BOOST_CHECK_CLOSE_FRACTION( thrustAccelerationsAlongTrajectory.begin( )->first, julianDate, 1.0E-14 );
        BOOST_CHECK_CLOSE_FRACTION( thrustAccelerationsAlongTrajectory.rbegin( )->first, julianDate + timeOfFlight, 1.0E-14 );

        // Loop over computed acceleration values and check peak acceleration
        for( auto it : thrustAccelerationsAlongTrajectory )
        {
            Eigen::Vector3d currentCartesianThrustAcceleration = it.second;
            if ( currentCartesianThrustAcceleration.norm() > peakThrustAcceleration )
            {
                peakThrustAcceleration = currentCartesianThrustAcceleration.norm();
            }
        }


        // Check results consistency w.r.t. Roegiers, T., Application of the Spherical Shaping Method to a Low-Thrust
        // Multiple Asteroid Rendezvous Mission, TU Delft (MSc thesis), 2014
        double expectedDeltaV = 5700.0;
        double expectedPeakAcceleration = 2.4e-4;

        // DeltaV provided with a precision of 5 m/s
        BOOST_CHECK_SMALL(std::fabs(sphericalShapingLegPointer->getLegDeltaV() - expectedDeltaV ), 5.0 );
        // Peak acceleration provided with a precision 2.0e-6 m/s^2
        BOOST_CHECK_SMALL(std::fabs(peakThrustAcceleration - expectedPeakAcceleration ), 2.0e-6 );


        // Retrieve state history
        std::map< double, Eigen::Vector6d > statesAlongTrajectory;
        sphericalShapingLegPointer->getStatesAlongTrajectory(statesAlongTrajectory, 10 );

        // Check initial and final time on output list
        BOOST_CHECK_CLOSE_FRACTION( statesAlongTrajectory.begin( )->first, julianDate, 1.0E-14 );
        BOOST_CHECK_CLOSE_FRACTION( statesAlongTrajectory.rbegin( )->first, julianDate + timeOfFlight, 1.0E-14 );

        // Check initial and final state on output list
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( initialState, statesAlongTrajectory.begin( )->second, 1.0E-5 );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( finalState, statesAlongTrajectory.rbegin( )->second, 1.0E-5 );
    }

}


//! Test.
BOOST_AUTO_TEST_CASE( test_spherical_shaping_earth_1989ML_transfer )
{
    spice_interface::loadStandardSpiceKernels( );

    int numberOfRevolutions = 1;
    double julianDate = 7799.5 * physical_constants::JULIAN_DAY;
    double timeOfFlight = 600.0 * physical_constants::JULIAN_DAY;

    // Ephemeris departure body.
    EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ApproximateJplEphemeris >( "Earth" );
    Eigen::Vector6d initialState = pointerToDepartureBodyEphemeris->getCartesianState( julianDate );

    // Final state derived from ML1989 ephemeris (from Spice).
    Eigen::Vector6d finalState = (
                Eigen::Vector6d() <<
                1.197701029846094E+00 * physical_constants::ASTRONOMICAL_UNIT,
                1.653518856610793E-01 * physical_constants::ASTRONOMICAL_UNIT,
                - 9.230177854743750E-02 * physical_constants::ASTRONOMICAL_UNIT,
                - 4.891080912584867E-05 * physical_constants::ASTRONOMICAL_UNIT / physical_constants::JULIAN_DAY,
                1.588950249593135E-02 * physical_constants::ASTRONOMICAL_UNIT / physical_constants::JULIAN_DAY,
                - 2.980245580772588E-04 * physical_constants::ASTRONOMICAL_UNIT / physical_constants::JULIAN_DAY ).finished();
    EphemerisPointer pointerToArrivalBodyEphemeris = std::make_shared< ConstantEphemeris >(finalState );

    // Define root finder settings (used to update the updated value of the free coefficient, so that it matches the required time of flight).
    std::shared_ptr< RootFinderSettings > rootFinderSettings =
            tudat::root_finders::bisectionRootFinderSettings( 1.0E-6, TUDAT_NAN, TUDAT_NAN, 30 );

    SphericalShapingLeg *sphericalShapingLegPointer;
    for ( unsigned int creationType = 0; creationType < 2; creationType++ )
    {
        if (creationType == 0) {
            sphericalShapingLegPointer = new SphericalShapingLeg(
                    pointerToDepartureBodyEphemeris, pointerToArrivalBodyEphemeris,
                    spice_interface::getBodyGravitationalParameter("Sun"),
                    rootFinderSettings, -1.0e-2, 1.0e-2);
        }
        else if (creationType == 1)
        {
            std::function< Eigen::Vector3d () > departureVelocityFunction = [=] () {return initialState.segment(3, 3);};
            std::function< Eigen::Vector3d () > arrivalVelocityFunction = [=] () {return finalState.segment(3, 3);};

            sphericalShapingLegPointer = new SphericalShapingLeg(
                    pointerToDepartureBodyEphemeris, pointerToArrivalBodyEphemeris,
                    spice_interface::getBodyGravitationalParameter("Sun"),
                    departureVelocityFunction, arrivalVelocityFunction,
                    rootFinderSettings, -1.0e-2, 1.0e-2);
        }

        sphericalShapingLegPointer->updateLegParameters(
                ( Eigen::Vector3d( ) << julianDate, julianDate + timeOfFlight, numberOfRevolutions ).finished( ));

        // Initialise peak acceleration
        double peakThrustAcceleration = 0.0;

        std::map< double, Eigen::Vector3d > thrustAccelerationsAlongTrajectory;
        sphericalShapingLegPointer->getThrustAccelerationsAlongTrajectory(thrustAccelerationsAlongTrajectory, 5000);

        // Check initial and final time on output list
        BOOST_CHECK_CLOSE_FRACTION(thrustAccelerationsAlongTrajectory.begin( )->first, julianDate, 1.0E-14);
        BOOST_CHECK_CLOSE_FRACTION(thrustAccelerationsAlongTrajectory.rbegin( )->first, julianDate + timeOfFlight,
                                   1.0E-14);

        // Loop over computed acceleration values and check peak acceleration
        for (auto it: thrustAccelerationsAlongTrajectory) {
            Eigen::Vector3d currentCartesianThrustAcceleration = it.second;
            if (currentCartesianThrustAcceleration.norm( ) > peakThrustAcceleration) {
                peakThrustAcceleration = currentCartesianThrustAcceleration.norm( );
            }
        }


        // Check results consistency w.r.t. Roegiers, T., Application of the Spherical Shaping Method to a Low-Thrust
        // Multiple Asteroid Rendezvous Mission, TU Delft (MSc thesis), 2014
        // The expected differences are a bit larger than for Earth-Mars transfer due to the higher uncertainty in 1989ML's ephemeris.
        double expectedDeltaV = 4530.0;
        double expectedPeakAcceleration = 1.8e-4;

        // DeltaV provided with a precision of 0.1 km/s
        BOOST_CHECK_SMALL(std::fabs(sphericalShapingLegPointer->getLegDeltaV( ) - expectedDeltaV), 100.0);
        // Peak acceleration provided with a precision 1.0e-5 m/s^2
        BOOST_CHECK_SMALL(std::fabs(peakThrustAcceleration - expectedPeakAcceleration), 1e-5);

        // Retrieve state history
        std::map< double, Eigen::Vector6d > statesAlongTrajectory;
        sphericalShapingLegPointer->getStatesAlongTrajectory(statesAlongTrajectory, 10);

        // Check initial and final time on output list
        BOOST_CHECK_CLOSE_FRACTION(statesAlongTrajectory.begin( )->first, julianDate, 1.0E-14);
        BOOST_CHECK_CLOSE_FRACTION(statesAlongTrajectory.rbegin( )->first, julianDate + timeOfFlight, 1.0E-14);

        // Check initial and final state on output list
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(initialState, statesAlongTrajectory.begin( )->second, 1.0E-5);
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(finalState, statesAlongTrajectory.rbegin( )->second, 1.0E-5);
    }
}


//! Test spherical shaping method with various number of revolutions.
BOOST_AUTO_TEST_CASE( test_spherical_shaping_earth_mars_transfer_multi_revolutions )
{
    spice_interface::loadStandardSpiceKernels( );

    double JD = physical_constants::JULIAN_DAY;
    double julianDate = 8174.5 * JD;

    std::vector< int > numberOfRevolutionsVector = { 0, 1, 2 };
    std::vector< double > timeOfFlightVector = { 300.0*JD, 580.0*JD, 750.0*JD };

    // Ephemeris departure body.
    EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ApproximateJplEphemeris>( "Earth" );
    EphemerisPointer pointerToArrivalBodyEphemeris = std::make_shared< ApproximateJplEphemeris >( "Mars" );

    // Define root finder settings (used to update the updated value of the free coefficient, so that it matches the required time of flight).
    std::shared_ptr< RootFinderSettings > rootFinderSettings =
            tudat::root_finders::bisectionRootFinderSettings( 1.0E-6, TUDAT_NAN, TUDAT_NAN, 30 );

    // Bounds for the free parameter.
    std::vector< double > freeParameterLowerBoundVector = { -1.0, 1.0e-6, -1.0e-2 };
    std::vector< double > freeParameterUpperBoundVector = { 5.0e-1, 1.0e-1, 1.0e-2 };

    // Define initial state.
    Eigen::Vector6d initialState = pointerToDepartureBodyEphemeris->getCartesianState( julianDate );

    for ( unsigned int currentTestCase = 0 ; currentTestCase < numberOfRevolutionsVector.size( ) ; currentTestCase++ )
    {
        // Define final state.
        Eigen::Vector6d finalState = pointerToArrivalBodyEphemeris->getCartesianState(
                    julianDate + timeOfFlightVector.at( currentTestCase ) );

        // Compute shaped trajectory.
        SphericalShapingLeg sphericalShapingLeg = SphericalShapingLeg(
                pointerToDepartureBodyEphemeris, pointerToArrivalBodyEphemeris,
                spice_interface::getBodyGravitationalParameter("Sun"), rootFinderSettings,
                freeParameterLowerBoundVector.at(currentTestCase), freeParameterUpperBoundVector.at(currentTestCase) );
        sphericalShapingLeg.updateLegParameters( (
                Eigen::Vector3d( )<< julianDate, julianDate + timeOfFlightVector.at( currentTestCase ),
                numberOfRevolutionsVector.at(currentTestCase) ).finished( ) );

        // Check consistency of final azimuth angle value with required number of revolutions.
        double initialAzimuthAngle = sphericalShapingLeg.getInitialValueAzimuth();
        double finalAzimuthAngle = sphericalShapingLeg.getFinalValueAzimuth();

        double expectedInitialAzimuthAngle = coordinate_conversions::convertCartesianToSphericalState( initialState )[ 1 ];
        if ( expectedInitialAzimuthAngle < 0.0 )
        {
            expectedInitialAzimuthAngle += 2.0 * mathematical_constants::PI;
        }
        double expectedFinalAzimuthAngle = coordinate_conversions::convertCartesianToSphericalState( finalState )[ 1 ];
        if ( expectedFinalAzimuthAngle < 0.0 )
        {
            expectedFinalAzimuthAngle += 2.0 * mathematical_constants::PI;
        }

        if ( expectedFinalAzimuthAngle - expectedInitialAzimuthAngle < 0 )
        {
            expectedFinalAzimuthAngle += 2.0 * mathematical_constants::PI * ( numberOfRevolutionsVector.at( currentTestCase )+ 1 );
        }
        else
        {
            expectedFinalAzimuthAngle += 2.0 * mathematical_constants::PI * numberOfRevolutionsVector.at( currentTestCase );
        }

        BOOST_CHECK_SMALL(std::fabs(initialAzimuthAngle - expectedInitialAzimuthAngle ), 1.0e-15 );
        BOOST_CHECK_SMALL(std::fabs(finalAzimuthAngle - expectedFinalAzimuthAngle ), 1.0e-15 );

        // Retrieve state history
        std::map< double, Eigen::Vector6d > statesAlongTrajectory;
        sphericalShapingLeg.getStatesAlongTrajectory( statesAlongTrajectory, 10 );

        // Check initial and final time on output list
        BOOST_CHECK_CLOSE_FRACTION( statesAlongTrajectory.begin( )->first, julianDate, 1.0E-14 );
        BOOST_CHECK_CLOSE_FRACTION( statesAlongTrajectory.rbegin( )->first, julianDate + timeOfFlightVector.at( currentTestCase ), 1.0E-14 );

        // Check initial and final state on output list
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( initialState, statesAlongTrajectory.begin( )->second, 1.0E-12 );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( finalState, statesAlongTrajectory.rbegin( )->second, 1.0E-12 );
    }
}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
