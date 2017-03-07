#define BOOST_TEST_MAIN

#include <iostream>

#include <boost/test/unit_test.hpp>

#include "Tudat/Basics/testMacros.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/InputOutput/basicInputOutput.h"

#include "Tudat/Astrodynamics/EarthOrientation/shortPeriodEarthOrientationCorrectionCalculator.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/External/SofaInterface/sofaTimeConversions.h"

namespace tudat
{

namespace unit_tests
{

using namespace tudat::basic_astrodynamics;

BOOST_AUTO_TEST_SUITE( test_short_period_polar_motion )


BOOST_AUTO_TEST_CASE( testShortPeriodLibrationalPolarMotion)
{
    using namespace tudat::earth_orientation;

    ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > polarMotionCalculator(
                unit_conversions::convertArcSecondsToRadians< double >( 1.0E-6 ), 0.0,
    { tudat::input_output::getDataFilesRootPath( ) + "EarthOrientation/polarMotionLibrationAmplitudesQuasiDiurnalOnly.txt" },
    { tudat::input_output::getDataFilesRootPath( ) + "EarthOrientation/polarMotionLibrationDoodsonMultipliersQuasiDiurnalOnly.txt" } );

    double testMjd = 54335.0;

    double testJulianDay = testMjd + JULIAN_DAY_AT_0_MJD;

    double testUtc = convertJulianDayToSecondsSinceEpoch( testJulianDay,
                                                          basic_astrodynamics::JULIAN_DAY_ON_J2000 );

    double testEphemerisTime = sofa_interface::convertUTCtoTT( testUtc );
    Eigen::Vector2d polarMotionCorrections = polarMotionCalculator.getCorrections( testEphemerisTime );

    double microAsToRadians =  mathematical_constants::PI / ( 180.0 * 1.0E6 * 3600.0 );

    // Difference between IERS code and this code occurs since IERS example does not distinguish between UTC and TDB for calculating slow arguments.
    // The difference that is introduced is < nas for this case.
    BOOST_CHECK_SMALL( polarMotionCorrections( 0 )/microAsToRadians - 24.83144238273364834, 5.0E-3 );
    BOOST_CHECK_SMALL( polarMotionCorrections( 1 )/microAsToRadians + 14.09240692041837661, 5.0E-3 );

}

BOOST_AUTO_TEST_CASE( testShortPeriodOceanTidePolarMotion)
{
    using namespace tudat::earth_orientation;

    ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > polarMotionCalculator(
                unit_conversions::convertArcSecondsToRadians< double >( 1.0E-6 ), 0.0,
    { tudat::input_output::getDataFilesRootPath( ) + "EarthOrientation/polarMotionOceanTidesAmplitudes.txt", },
    { tudat::input_output::getDataFilesRootPath( ) + "EarthOrientation/polarMotionOceanTidesDoodsonMultipliers.txt" } );

    double testMjd = 47100.0;

    double testJulianDay = testMjd + JULIAN_DAY_AT_0_MJD;

    double testUtc = convertJulianDayToSecondsSinceEpoch( testJulianDay,
                                                          basic_astrodynamics::JULIAN_DAY_ON_J2000 );

    double testEphemerisTime = sofa_interface::convertUTCtoTT( testUtc );

    double microAsToRadians =  mathematical_constants::PI / ( 180.0 * 1.0E6 * 3600.0 );

    Eigen::Vector2d polarMotionCorrections = polarMotionCalculator.getCorrections( testEphemerisTime ) / microAsToRadians;

    // Difference between result and IERS code is due to different implementation (ortho-weights vs. explicit arguments), difference is still at << for Earth stations).
    BOOST_CHECK_SMALL( polarMotionCorrections.x( ) + 162.8386373279636530, 1.0 );
    BOOST_CHECK_SMALL( polarMotionCorrections.y( ) - 117.7907525842668974, 1.0 );


}

BOOST_AUTO_TEST_SUITE_END( )

}

}



