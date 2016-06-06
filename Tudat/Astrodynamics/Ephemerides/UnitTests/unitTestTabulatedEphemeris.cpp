/*    Copyright (c) 2010-2015, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author              Comment
 *      141105    D. Dirkx            File created.
 *
 *    References
 *
 *    Notes
 */

#define BOOST_TEST_MAIN

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Basics/testMacros.h"

#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositions.h"
#include "Tudat/Astrodynamics/Ephemerides/tabulatedEphemeris.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h"
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
std::map< double, basic_mathematics::Vector6d > getStateHistoryMap(
        const boost::shared_ptr< ephemerides::Ephemeris > originalEphemeris  )
{
    std::map< double, basic_mathematics::Vector6d > stateHistoryMap;

    // Define time limits and step
    double startTime = 0.0;
    double finalTime = 1.0E7;
    double timeStep = 1000.0;

    // Create time history map.
    double currentTime = startTime;
    while( currentTime <= finalTime )
    {
        stateHistoryMap[ currentTime ] = originalEphemeris->getCartesianStateFromEphemeris(
                    currentTime, basic_astrodynamics::JULIAN_DAY_ON_J2000 );
        currentTime += timeStep;
    }

    return stateHistoryMap;
}


//! Test the functionality of the tabulated ephemeris
BOOST_AUTO_TEST_CASE( testTabulatedEphemeris )
{
    // Create ephemeris from which table is to be generated; used as input to tabulated ephemeris
    using namespace ephemerides;
    boost::shared_ptr< ApproximatePlanetPositions > marsNominalEphemeris =
            boost::make_shared< ApproximatePlanetPositions >(
                ApproximatePlanetPositionsBase::mars );

    // Generate state history map from ephemeris
    std::map< double, basic_mathematics::Vector6d > marsStateHistoryMap = getStateHistoryMap(
                marsNominalEphemeris );

    // Create interpolator from state history map.
    boost::shared_ptr< interpolators::OneDimensionalInterpolator
            < double, basic_mathematics::Vector6d > > marsStateInterpolator =
            boost::make_shared< interpolators::CubicSplineInterpolator
            < double, basic_mathematics::Vector6d > >( marsStateHistoryMap );

    // Create tabulated epehemeris from interpolator
    boost::shared_ptr< TabulatedCartesianEphemeris< > > tabulatedEphemeris =
            boost::make_shared< TabulatedCartesianEphemeris< > >(
                marsStateInterpolator, "SSB", "J2000", basic_astrodynamics::JULIAN_DAY_ON_J2000 );

    // Compare interpolated and tabulated ephemeris state at dummy time.
    double testTime = 1.9337E5;
    basic_mathematics::Vector6d interpolatorState = marsStateInterpolator->interpolate( testTime );
    basic_mathematics::Vector6d ephemerisState = tabulatedEphemeris->getCartesianStateFromEphemeris(
                testTime, basic_astrodynamics::JULIAN_DAY_ON_J2000 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( interpolatorState, ephemerisState, 0.0 );

    // Compare direct and tabulated ephemeris state at dummy time (comparison not equal due to
    // interpolation errors).
    basic_mathematics::Vector6d directState = marsNominalEphemeris->getCartesianStateFromEphemeris(
                testTime, basic_astrodynamics::JULIAN_DAY_ON_J2000 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( directState, ephemerisState, 1.0E-10 );
    testTime = 5.836392E6;

    // Compare interpolated and tabulated ephemeris state at second dummy time.
    interpolatorState = marsStateInterpolator->interpolate( testTime );
    ephemerisState = tabulatedEphemeris->getCartesianStateFromEphemeris(
                testTime, basic_astrodynamics::JULIAN_DAY_ON_J2000 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( interpolatorState, ephemerisState, 0.0 );

    // Compare direct and tabulated ephemeris state at second dummy time (comparison not equal due
    // to interpolation errors).
    directState = marsNominalEphemeris->getCartesianStateFromEphemeris(
                testTime, basic_astrodynamics::JULIAN_DAY_ON_J2000 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( directState, ephemerisState, 1.0E-10 );

    // Create new ephemeris from which table is to be generated; used to reset input to tabulated
    // ephemeris
    boost::shared_ptr< ApproximatePlanetPositions > jupiterNominalEphemeris =
            boost::make_shared< ApproximatePlanetPositions >(
                ApproximatePlanetPositionsBase::jupiter );

    // Reset tabulated ephemeris data.
    std::map< double, basic_mathematics::Vector6d > jupiterStateHistoryMap = getStateHistoryMap(
                jupiterNominalEphemeris );
    boost::shared_ptr< interpolators::OneDimensionalInterpolator
            < double, basic_mathematics::Vector6d > > jupiterStateInterpolator =
            boost::make_shared< interpolators::CubicSplineInterpolator
            < double, basic_mathematics::Vector6d > >( jupiterStateHistoryMap );
    tabulatedEphemeris->resetInterpolator( jupiterStateInterpolator,
                                           basic_astrodynamics::JULIAN_DAY_ON_J2000 );

    // Test tabulated ephemeris with reset data.
    interpolatorState = jupiterStateInterpolator->interpolate( testTime );
    ephemerisState = tabulatedEphemeris->getCartesianStateFromEphemeris(
                testTime, basic_astrodynamics::JULIAN_DAY_ON_J2000 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( interpolatorState, ephemerisState, 0.0 );
    directState = jupiterNominalEphemeris->getCartesianStateFromEphemeris(
                testTime, basic_astrodynamics::JULIAN_DAY_ON_J2000 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( directState, ephemerisState, 1.0E-10 );

    // Test whether check of reference epoch is succesfull. The following is expected to fail.
    bool isErrorCaught = 0;
    try
    {
        ephemerisState = tabulatedEphemeris->getCartesianStateFromEphemeris(
                    testTime, basic_astrodynamics::JULIAN_DAY_ON_J2000 - 1.0 );
    }
    catch( std::runtime_error )
    {
        // Store the fact that a runtime error occurred, such that the values will be stored.
        isErrorCaught = 1;
    }
    BOOST_CHECK_EQUAL( isErrorCaught, 1 );

    // Check whether getting of interpolator is correct
    BOOST_CHECK_EQUAL( tabulatedEphemeris->getInterpolator( ), jupiterStateInterpolator );


}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
