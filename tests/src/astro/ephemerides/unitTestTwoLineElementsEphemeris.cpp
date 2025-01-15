/*    Copyright (c) 2010-2020, Delft University of Technology
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

#include "tudat/astro/ephemerides/tleEphemeris.h"
#include "tudat/astro/basic_astro/dateTime.h"
#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"

namespace tudat
{
namespace unit_tests
{

//! Test the functionality of the two-line elements ephemeris
BOOST_AUTO_TEST_SUITE( test_two_line_elements_ephemeris )


//! Test the functionality of the two-line elements ephemeris
BOOST_AUTO_TEST_CASE( testTwoLineElementsEphemerisVallado )
{
	using namespace tudat;
	using namespace tudat::ephemerides;

	// Dummy two line element set from Vallado (2013), page 234
	std::string elements =  "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753\n"
							"2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667";

	// Create TLE object from element set
	std::shared_ptr< Tle > tlePtr = std::make_shared< Tle >( elements );

	// Create TLE Ephemeris object; output state in J2000 frame fixed at the centre of the Earth
	std::shared_ptr< TleEphemeris > tleEphemeris = std::make_shared< TleEphemeris >( "Earth", "J2000", tlePtr, false );

	// Propagate TLE for 3 Julian days to be able to compare state to the one given in Vallado
    double utcEpoch =  3.0 * tudat::physical_constants::JULIAN_DAY + tlePtr->getEpoch();
    double tdbEpoch = earth_orientation::createDefaultTimeConverter( )->getCurrentTime( basic_astrodynamics::utc_scale, basic_astrodynamics::tdb_scale, utcEpoch );

	Eigen::Vector6d propagatedState = tleEphemeris->getCartesianState( tdbEpoch );
	Eigen::Vector3d propagatedPosition = propagatedState.head( 3 );
	Eigen::Vector3d propagatedVelocity = propagatedState.tail( 3 );

	// Position and velocity vectors as found in Vallado
	Eigen::Vector3d verificationPosition;
	Eigen::Vector3d verificationVelocity;
	verificationPosition << -9059941.3786, 4659697.2000, 813958.8875;
	verificationVelocity << -2233.348094, -4110.136162, -3157.394074;

	// Check if difference within tolerances
    BOOST_CHECK_SMALL( ( verificationPosition - propagatedPosition ).norm( ), 50.0 );
    BOOST_CHECK_SMALL( ( propagatedVelocity - verificationVelocity ).norm( ), 0.05 );


}

BOOST_AUTO_TEST_CASE( testTwoLineElementsEphemerisGehly )
{
    using namespace tudat;
    using namespace tudat::ephemerides;

    //Retrieve the initial state of Delfi-C3 using Two-Line-Elements (TLEs)
        std::string elements = "1 32789U 07021G   08119.60740078 -.00000054  00000-0  00000+0 0  9999\n"
                               "2 32789 098.0082 179.6267 0015321 307.2977 051.0656 14.81417433    68";


    std::shared_ptr< Tle > delfiTle = std::make_shared< Tle >( elements );

    // Create TLE Ephemeris object; output state in J2000 frame fixed at the centre of the Earth
    std::shared_ptr< TleEphemeris > tleEphemeris = std::make_shared< TleEphemeris >( "Earth", "J2000", delfiTle, false );


    basic_astrodynamics::DateTime currentDate = basic_astrodynamics::DateTime(
        2008, 4, 28, 14, 34, 39.427392 );

    // Create time scale converter object
    auto timeScaleConverter = earth_orientation::createDefaultTimeConverter( );

    // Set the epoch in UTC scale (for instance from the above example using DateTime)
    double tdbTime = timeScaleConverter->getCurrentTime< double >(
        basic_astrodynamics::utc_scale, basic_astrodynamics::tdb_scale,
        currentDate.epoch< double >( ) );

    Eigen::VectorXd currentState  = tleEphemeris->getCartesianState( tdbTime );
    Eigen::VectorXd referenceState =
        ( Eigen::Vector6d( ) << -6.99682490e+06,  3.09944881e+04, -1.90322410e+05,
        -1.99894020e+02,  1.05296158e+03,  7.47516550e+03 ).finished( );

    // Check if difference within tolerances
    BOOST_CHECK_SMALL( ( currentState - referenceState ).segment( 0, 3 ).norm( ), 50.0 );
    BOOST_CHECK_SMALL( ( currentState - referenceState ).segment( 3, 3 ).norm( ), 0.05 );

}

BOOST_AUTO_TEST_CASE( testTwoLineElementsEphemerisLangbroke )
{
    using namespace tudat;
    using namespace tudat::ephemerides;

    std::vector< std::string > elementsList1;
    elementsList1.push_back( "1  6073U 72023E   25002.25894564  .00281450  96007-5  26991-3 0  9999");
    elementsList1.push_back( "1  6073U 72023E   23002.27791083  .00058038  63497-5  17711-3 0  9992");
    elementsList1.push_back( "1  6073U 72023E   24002.18650274  .00066887  72730-5  14914-3 0  9992");
    elementsList1.push_back( "1  6073U 72023E   22002.17600020  .00039319  66505-5  14175-3 0  9995");
    elementsList1.push_back( "1 06073U 72023E   01002.11505704  .00021795  12501-5  21396-3 0  9998");
    elementsList1.push_back( "1 06073U 72023E   00002.21547433 +.00023639 +15357-5 +24010-3 0  9990");

    std::vector< std::string > elementsList2;
    elementsList2.push_back( "2  6073  51.9788 162.9677 0455091 357.5947   2.2848 15.25405973731352" );
    elementsList2.push_back( "2  6073  52.0272 335.7777 1067120 299.8132  50.0834 13.74925322626299" );
    elementsList2.push_back( "2  6073  52.0129 331.6919 0834584 268.0601  82.4676 14.30697539677464" );
    elementsList2.push_back( "2  6073  52.0268 243.1202 1221888  42.8126 326.1665 13.38718442576770" );
    elementsList2.push_back( "2 06073 052.1205 048.6650 2751755 317.6011 024.1509 10.01527228685102" );
    elementsList2.push_back( "2 06073 052.0987 075.0950 2837480 139.3635 245.7845 09.83846471648763" );


    std::vector< Eigen::Vector6d > j2000States;
    j2000States.push_back( 1000.0 * ( Eigen::Vector6d( ) << -006257.80116,  +001955.41773,  +000015.12941,  -01.45696,  -04.68443,  +06.28855 ).finished( ) );
    j2000States.push_back( 1000.0 * ( Eigen::Vector6d( ) << +006288.31206,  -002867.99337,  -000013.82935,  +02.62122,  +04.08080,  +06.14099 ).finished( ) );
    j2000States.push_back( 1000.0 * ( Eigen::Vector6d( ) << +006267.43452,  -003419.32884,  -000014.42114,  +02.75547,  +03.73557,  +05.87787 ).finished( ) );
    j2000States.push_back( 1000.0 * ( Eigen::Vector6d( ) << -003091.90412,  -006026.12550,  +000006.62649,  +04.67822,  -01.70085,  +06.30628 ).finished( ) );
    j2000States.push_back( 1000.0 * ( Eigen::Vector6d( ) << +004615.97301,  +005243.96807,  -000001.45652,  -02.97643,  +04.31717,  +06.54252  ).finished( ) );
    j2000States.push_back( 1000.0 * ( Eigen::Vector6d( ) << +002772.31830,  +010415.43162,  +000000.35614,  -03.52556,  -00.37951,  +04.25257  ).finished( ) );

    std::vector< basic_astrodynamics::DateTime > evaluationTimes;
    evaluationTimes.push_back( basic_astrodynamics::DateTime( 2025, 01, 02, 06, 12, 52.90 ) );
    evaluationTimes.push_back( basic_astrodynamics::DateTime( 2023, 01, 02, 06, 40, 11.50 ) );
    evaluationTimes.push_back( basic_astrodynamics::DateTime( 2024, 01, 02, 04, 28, 33.84 ) );
    evaluationTimes.push_back( basic_astrodynamics::DateTime( 2022, 01, 02, 04, 13, 26.42 ) );
    evaluationTimes.push_back( basic_astrodynamics::DateTime( 2001, 01, 02, 02, 45, 40.93) );
    evaluationTimes.push_back( basic_astrodynamics::DateTime( 2000, 01, 02, 05, 10, 16.98) );

//    evaluationTimes.push_back( basic_astrodynamics::DateTime( 2025, 01, 02, 06, 12, 52.90 ) );

    // Create time scale converter object
    auto timeScaleConverter = earth_orientation::createDefaultTimeConverter( );


    for( int i = 0; i < evaluationTimes.size( ); i++ )
    {
        std::cout<<i<<std::endl;
        std::shared_ptr<Tle> currentTle = std::make_shared<Tle>( elementsList1.at( i ), elementsList2.at( i ) );

        // Create TLE Ephemeris object; output state in J2000 frame fixed at the centre of the Earth
        std::shared_ptr< TleEphemeris > tleEphemeris = std::make_shared<TleEphemeris>( "Earth", "J2000", currentTle, false );

        basic_astrodynamics::DateTime currentDate = evaluationTimes.at( i );

        // Set the epoch in UTC scale (for instance from the above example using DateTime)
        double tdbTime = timeScaleConverter->getCurrentTime<double>(
            basic_astrodynamics::utc_scale, basic_astrodynamics::tdb_scale,
            currentDate.epoch<double>( ) );

        Eigen::VectorXd currentState = tleEphemeris->getCartesianState( tdbTime );
        Eigen::VectorXd referenceState = j2000States.at( i );

        // Check if difference within tolerances
        BOOST_CHECK_SMALL(( currentState - referenceState ).segment( 0, 3 ).norm( ), 1.0 );
        BOOST_CHECK_SMALL(( currentState - referenceState ).segment( 3, 3 ).norm( ), 0.2 );
    }

}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
