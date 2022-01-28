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

#include "tudat/basics/testMacros.h"

#include "tudat/astro/ephemerides/approximatePlanetPositions.h"
#include "tudat/astro/ephemerides/tabulatedEphemeris.h"
#include "tudat/basics/basicTypedefs.h"
#include "tudat/math/interpolators/cubicSplineInterpolator.h"

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
template< typename VectorType = Eigen::Vector6d >
std::map< double, VectorType > getStateHistoryMap(
        const std::shared_ptr< ephemerides::Ephemeris > originalEphemeris  )
{
    std::map< double, VectorType > stateHistoryMap;

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
    std::shared_ptr< ApproximateJplEphemeris > marsNominalEphemeris =
            std::make_shared< ApproximateJplEphemeris >(
                "Mars" );

    // Generate state history map from ephemeris
    std::map< double, Eigen::Vector6d > marsStateHistoryMap = getStateHistoryMap< Eigen::Vector6d >(
                marsNominalEphemeris );
    std::map< double, Eigen::VectorXd > marsDynamicSizeStateHistoryMap = getStateHistoryMap< Eigen::VectorXd >(
                marsNominalEphemeris );

    // Create interpolator from state history map.
    std::shared_ptr< interpolators::OneDimensionalInterpolator
            < double, Eigen::Vector6d > > marsStateInterpolator =
            std::make_shared< interpolators::CubicSplineInterpolator
            < double, Eigen::Vector6d > >( marsStateHistoryMap );
    std::shared_ptr< interpolators::OneDimensionalInterpolator
            < double, Eigen::VectorXd > > marsDynamicSizeStateInterpolator =
            std::make_shared< interpolators::CubicSplineInterpolator
            < double, Eigen::VectorXd > >( marsDynamicSizeStateHistoryMap );

    // Create tabulated epehemeris from interpolator
    std::shared_ptr< TabulatedCartesianEphemeris< > > tabulatedEphemeris =
            std::make_shared< TabulatedCartesianEphemeris< > >(
                marsStateInterpolator, "SSB", "J2000");
    std::shared_ptr< TabulatedCartesianEphemeris< > > tabulatedDynamicSizeEphemeris =
            std::make_shared< TabulatedCartesianEphemeris< > >(
                marsDynamicSizeStateInterpolator, "SSB", "J2000");

    // Compare interpolated and tabulated ephemeris state at dummy time.
    double testTime = 1.9337E5;
    Eigen::Vector6d interpolatorState = marsStateInterpolator->interpolate( testTime );
    Eigen::VectorXd interpolatorDynamicSizeState = marsStateInterpolator->interpolate( testTime );

    Eigen::Vector6d ephemerisState = tabulatedEphemeris->getCartesianState(
                testTime);
    Eigen::Vector6d ephemerisDynamicSizeState = tabulatedDynamicSizeEphemeris->getCartesianState(
                testTime);
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( interpolatorState, ephemerisState, 0.0 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( interpolatorDynamicSizeState, ephemerisDynamicSizeState, 0.0 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( interpolatorDynamicSizeState, interpolatorState, 0.0 );

    // Compare direct and tabulated ephemeris state at dummy time (comparison not equal due to
    // interpolation errors).
    Eigen::Vector6d directState = marsNominalEphemeris->getCartesianState(
                testTime);
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( directState, ephemerisState, 1.0E-10 );
    testTime = 5.836392E6;

    // Compare interpolated and tabulated ephemeris state at second dummy time.
    interpolatorState = marsStateInterpolator->interpolate( testTime );
    interpolatorDynamicSizeState = marsStateInterpolator->interpolate( testTime );
    ephemerisState = tabulatedEphemeris->getCartesianState( testTime);
    ephemerisDynamicSizeState = tabulatedDynamicSizeEphemeris->getCartesianState( testTime );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( interpolatorState, ephemerisState, 0.0 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( interpolatorDynamicSizeState, ephemerisDynamicSizeState, 0.0 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( interpolatorDynamicSizeState, interpolatorState, 0.0 );

    // Compare direct and tabulated ephemeris state at second dummy time (comparison not equal due
    // to interpolation errors).
    directState = marsNominalEphemeris->getCartesianState(
                testTime);
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( directState, ephemerisState, 1.0E-10 );

    // Create new ephemeris from which table is to be generated; used to reset input to tabulated
    // ephemeris
    std::shared_ptr< ApproximateJplEphemeris > jupiterNominalEphemeris =
            std::make_shared< ApproximateJplEphemeris >(
                "Jupiter" );

    // Reset tabulated ephemeris data.
    std::map< double, Eigen::Vector6d > jupiterStateHistoryMap = getStateHistoryMap(
                jupiterNominalEphemeris );
    std::map< double, Eigen::VectorXd > jupiterDynamicSizeStateHistoryMap = getStateHistoryMap< Eigen::VectorXd >(
                jupiterNominalEphemeris );
    std::shared_ptr< interpolators::OneDimensionalInterpolator
            < double, Eigen::Vector6d > > jupiterStateInterpolator =
            std::make_shared< interpolators::CubicSplineInterpolator
            < double, Eigen::Vector6d > >( jupiterStateHistoryMap );
    std::shared_ptr< interpolators::OneDimensionalInterpolator
            < double, Eigen::VectorXd > > jupiterDynamicSizeStateInterpolator =
            std::make_shared< interpolators::CubicSplineInterpolator
            < double, Eigen::VectorXd > >( jupiterDynamicSizeStateHistoryMap );

    tabulatedEphemeris->resetInterpolator( jupiterStateInterpolator );
    tabulatedDynamicSizeEphemeris->resetInterpolator( jupiterDynamicSizeStateInterpolator );

    // Test tabulated ephemeris with reset data.
    interpolatorState = jupiterStateInterpolator->interpolate( testTime );
    interpolatorDynamicSizeState = jupiterDynamicSizeStateInterpolator->interpolate( testTime );
    ephemerisState = tabulatedEphemeris->getCartesianState( testTime);
    ephemerisDynamicSizeState = tabulatedDynamicSizeEphemeris->getCartesianState( testTime );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( interpolatorState, ephemerisState, 0.0 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( interpolatorDynamicSizeState, ephemerisDynamicSizeState, 0.0 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( interpolatorDynamicSizeState, interpolatorState, 0.0 );

    directState = jupiterNominalEphemeris->getCartesianState(
                testTime);
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( directState, ephemerisState, 1.0E-10 );


    // Check whether getting of interpolator is correct
    BOOST_CHECK_EQUAL( tabulatedEphemeris->getInterpolator( ), jupiterStateInterpolator );


}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
