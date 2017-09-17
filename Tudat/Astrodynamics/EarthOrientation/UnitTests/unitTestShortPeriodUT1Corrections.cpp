#define BOOST_TEST_MAIN #include <iostream> #include <boost/test/unit_test.hpp>

#include "Tudat/Basics/testMacros.h"

#include "Tudat/Astrodynamics/EarthOrientation/shortPeriodEarthOrientationCorrectionCalculator.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"

namespace tudat
{

namespace unit_tests
{

using namespace tudat::basic_astrodynamics;

BOOST_AUTO_TEST_SUITE( test_short_period_polar_motion )

//! Test short-periodic ut1-utc variations by comparing to test output of UTLIBR.f file
BOOST_AUTO_TEST_CASE( testShortPeriodLibrationalPolarMotion)
{
    using namespace tudat::earth_orientation;

    ShortPeriodEarthOrientationCorrectionCalculator < double > ut1CorrectionCalculator(
                1.0E-6, 0.0,
    { tudat::input_output::getEarthOrientationDataFilesPath( ) + "utcLibrationAmplitudes.txt" },
    { tudat::input_output::getEarthOrientationDataFilesPath( ) + "utcLibrationDoodsonMultipliers.txt" } );


    double testMjd = 44239.1 ;
    double testJulianDay = testMjd + JULIAN_DAY_AT_0_MJD;
    double testUtc = convertJulianDayToSecondsSinceEpoch( testJulianDay,
                                                          basic_astrodynamics::JULIAN_DAY_ON_J2000 );

    double testEphemerisTime = sofa_interface::convertUTCtoTT( testUtc );

    double ut1Correction = ut1CorrectionCalculator.getCorrections( testEphemerisTime );

    BOOST_CHECK_SMALL( ut1Correction - 2.441143834386761746E-6, 1.0E-8 );

    testMjd = 55227.4 ;
    testJulianDay = testMjd + JULIAN_DAY_AT_0_MJD;
    testUtc = convertJulianDayToSecondsSinceEpoch( testJulianDay,
                                                   basic_astrodynamics::JULIAN_DAY_ON_J2000 );

    testEphemerisTime = sofa_interface::convertUTCtoTT( testUtc );

    ut1Correction = ut1CorrectionCalculator.getCorrections( testEphemerisTime );

    std::cout<<ut1Correction + 2.655705844335680244E-6<<std::endl;
    BOOST_CHECK_SMALL( ut1Correction + 2.655705844335680244E-6, 1.0E-8 );

}

//! Test short-periodic ut1-utc variations by comparing to test output of ORTHO_EOP.f file
BOOST_AUTO_TEST_CASE( testShortPeriodOceanTidePolarMotion )
{
    using namespace tudat::earth_orientation;

    ShortPeriodEarthOrientationCorrectionCalculator < double > ut1CorrectionCalculator =
            ShortPeriodEarthOrientationCorrectionCalculator < double >(
                1.0E-6, 0.0,
    { tudat::input_output::getEarthOrientationDataFilesPath( ) + "utcOceanTidesAmplitudes.txt", },
    { tudat::input_output::getEarthOrientationDataFilesPath( ) + "utcOceanTidesDoodsonMultipliers.txt" } );

    double testMjd = 47100.0;
    double testJulianDay = testMjd + JULIAN_DAY_AT_0_MJD;
    double testUtc = convertJulianDayToSecondsSinceEpoch( testJulianDay,
                                                          basic_astrodynamics::JULIAN_DAY_ON_J2000 );

    double testEphemerisTime = sofa_interface::convertUTCtoTT( testUtc );

    double ut1Correction = ut1CorrectionCalculator.getCorrections( testEphemerisTime );

    std::cout<<ut1Correction + 23.39092370609808214E-6<<std::endl;

    BOOST_CHECK_SMALL( ut1Correction + 23.39092370609808214E-6, 5.0E-8 );

}

BOOST_AUTO_TEST_SUITE_END( )

}

}




