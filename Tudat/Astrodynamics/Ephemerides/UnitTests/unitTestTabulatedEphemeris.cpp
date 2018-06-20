/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#define BOOST_TEST_MAIN

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Basics/testMacros.h"

#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositions.h"
#include "Tudat/Astrodynamics/Ephemerides/tabulatedEphemeris.h"
#include "Tudat/Basics/basicTypedefs.h"
#include "Tudat/Mathematics/Interpolators/cubicSplineInterpolator.h"

namespace tudat
{
namespace unit_tests
{

//! Test the functionality of the tabulated ephemeris
BOOST_AUTO_TEST_SUITE( test_tabulated_ephemeris )

//! Function to create state history map from any ephemeris.
/*!
 *  Function to create state history map from any ephemeris. Uses arbitrary hardcoded input times
 *  and is used for comparison of tabulated and other ephemerides.
 *  \param originalEphemeris Ephemeris from which table is to be generated
 *  \return State history map generated from given ephemeris
 */
std::map< double, Eigen::Vector6d > getStateHistoryMap(
        const std::shared_ptr< ephemerides::Ephemeris > originalEphemeris  )
{
    std::map< double, Eigen::Vector6d > stateHistoryMap;

    // Define time limits and step
    double startTime = 0.0;
    double finalTime = 1.0E7;
    double timeStep = 1000.0;

    // Create time history map.
    double currentTime = startTime;
    while( currentTime <= finalTime )
    {
        stateHistoryMap[ currentTime ] = originalEphemeris->getCartesianState(
                    currentTime);
        currentTime += timeStep;
    }

    return stateHistoryMap;
}


//! Test the functionality of the tabulated ephemeris
BOOST_AUTO_TEST_CASE( testTabulatedEphemeris )
{
    // Create ephemeris from which table is to be generated; used as input to tabulated ephemeris
    using namespace ephemerides;
    std::shared_ptr< ApproximatePlanetPositions > marsNominalEphemeris =
            std::make_shared< ApproximatePlanetPositions >(
                ApproximatePlanetPositionsBase::mars );

    // Generate state history map from ephemeris
    std::map< double, Eigen::Vector6d > marsStateHistoryMap = getStateHistoryMap(
                marsNominalEphemeris );

    // Create interpolator from state history map.
    std::shared_ptr< interpolators::OneDimensionalInterpolator
            < double, Eigen::Vector6d > > marsStateInterpolator =
            std::make_shared< interpolators::CubicSplineInterpolator
            < double, Eigen::Vector6d > >( marsStateHistoryMap );

    // Create tabulated epehemeris from interpolator
    std::shared_ptr< TabulatedCartesianEphemeris< > > tabulatedEphemeris =
            std::make_shared< TabulatedCartesianEphemeris< > >(
                marsStateInterpolator, "SSB", "J2000");

    // Compare interpolated and tabulated ephemeris state at dummy time.
    double testTime = 1.9337E5;
    Eigen::Vector6d interpolatorState = marsStateInterpolator->interpolate( testTime );
    Eigen::Vector6d ephemerisState = tabulatedEphemeris->getCartesianState(
                testTime);
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( interpolatorState, ephemerisState, 0.0 );

    // Compare direct and tabulated ephemeris state at dummy time (comparison not equal due to
    // interpolation errors).
    Eigen::Vector6d directState = marsNominalEphemeris->getCartesianState(
                testTime);
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( directState, ephemerisState, 1.0E-10 );
    testTime = 5.836392E6;

    // Compare interpolated and tabulated ephemeris state at second dummy time.
    interpolatorState = marsStateInterpolator->interpolate( testTime );
    ephemerisState = tabulatedEphemeris->getCartesianState(
                testTime);
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( interpolatorState, ephemerisState, 0.0 );

    // Compare direct and tabulated ephemeris state at second dummy time (comparison not equal due
    // to interpolation errors).
    directState = marsNominalEphemeris->getCartesianState(
                testTime);
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( directState, ephemerisState, 1.0E-10 );

    // Create new ephemeris from which table is to be generated; used to reset input to tabulated
    // ephemeris
    std::shared_ptr< ApproximatePlanetPositions > jupiterNominalEphemeris =
            std::make_shared< ApproximatePlanetPositions >(
                ApproximatePlanetPositionsBase::jupiter );

    // Reset tabulated ephemeris data.
    std::map< double, Eigen::Vector6d > jupiterStateHistoryMap = getStateHistoryMap(
                jupiterNominalEphemeris );
    std::shared_ptr< interpolators::OneDimensionalInterpolator
            < double, Eigen::Vector6d > > jupiterStateInterpolator =
            std::make_shared< interpolators::CubicSplineInterpolator
            < double, Eigen::Vector6d > >( jupiterStateHistoryMap );
    tabulatedEphemeris->resetInterpolator( jupiterStateInterpolator );

    // Test tabulated ephemeris with reset data.
    interpolatorState = jupiterStateInterpolator->interpolate( testTime );
    ephemerisState = tabulatedEphemeris->getCartesianState(
                testTime);
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( interpolatorState, ephemerisState, 0.0 );
    directState = jupiterNominalEphemeris->getCartesianState(
                testTime);
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( directState, ephemerisState, 1.0E-10 );


    // Check whether getting of interpolator is correct
    BOOST_CHECK_EQUAL( tabulatedEphemeris->getInterpolator( ), jupiterStateInterpolator );


}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
