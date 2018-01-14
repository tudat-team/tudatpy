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
#include <boost/make_shared.hpp>

#include "Tudat/Basics/testMacros.h"

#include "Tudat/Astrodynamics/Ephemerides/itrsToGcrsRotationModel.h"
#include "Tudat/Astrodynamics/EarthOrientation/UnitTests/sofaEarthOrientationCookbookExamples.h"
#include "Tudat/External/SpiceInterface/spiceInterface.h"

namespace tudat
{
namespace unit_tests
{

using namespace ephemerides;
using namespace earth_orientation;
using namespace basic_astrodynamics;

BOOST_AUTO_TEST_SUITE( test_itrs_to_gcrs_rotation )

//! Test ITRS <-> GCRS rotation by compariong against Spice
BOOST_AUTO_TEST_CASE( test_ItrsToGcrsRotationAgainstSpice )
{

    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "naif0012.tls" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "earth_latest_high_prec.bpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "earth_fixed.tf" );

    // Create rotation model
    boost::shared_ptr< GcrsToItrsRotationModel > earthRotationModel =
            boost::make_shared< GcrsToItrsRotationModel >(
                earth_orientation::createStandardEarthOrientationCalculator( ) );

    // Compare spice vs. Tudat for list of evaluation times
    std::vector< double > testTimes;
    testTimes.push_back( 1.0E8 );
    testTimes.push_back( 1.0E7 );
    testTimes.push_back( 0.0 );
    for( unsigned test = 0; test < testTimes.size( ); test++ )
    {
        Eigen::Matrix3d sofaRotation = earthRotationModel->getRotationToBaseFrame( testTimes.at( test ) ).toRotationMatrix( );
        Eigen::Matrix3d sofaRotationDerivative = earthRotationModel->getDerivativeOfRotationToBaseFrame( testTimes.at( test ) );
        Eigen::Matrix3d spiceRotation = spice_interface::computeRotationQuaternionBetweenFrames(
                    "ITRF93", "J2000", testTimes.at( test ) ).toRotationMatrix( );
        Eigen::Matrix3d spiceRotationDerivative = spice_interface::computeRotationMatrixDerivativeBetweenFrames(
                    "ITRF93", "J2000", testTimes.at( test ) );

        // Check whether Spice and Tudat give same result. Note that Spice model is not accurate up to IERS standards. Comparison
        // is done at 10 cm position difference on Earth surface (per component).
        double tolerance = 0.1 / 6378.0E3;
        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL( sofaRotation( i, j ) - spiceRotation( i, j ), tolerance );
                BOOST_CHECK_SMALL( sofaRotationDerivative( i, j ) - spiceRotationDerivative( i, j ), 5.0E-12 );

            }
        }

        std::cout << sofaRotationDerivative << std::endl << std::endl
                  << sofaRotationDerivative - spiceRotationDerivative << std::endl << std::endl << std::endl << std::endl;
    }
}

//! Test ITRS <-> GCRS rotation by compariong against Sofa
BOOST_AUTO_TEST_CASE( test_ItrsToGcrsRotationAgainstSofaCookbook )
{

    // Get UTC time for evaluation
    int year = 2007;
    int month = 4;
    int day = 5;
    int hour = 12;
    int minutes = 0;
    double seconds = 0.0;
    double sofaCookbookTime = convertCalendarDateToJulianDaysSinceEpoch(
                year, month, day, hour, minutes, seconds, JULIAN_DAY_ON_J2000 ) *
            physical_constants::JULIAN_DAY;

    // Create Earth rotation model
    boost::shared_ptr< GcrsToItrsRotationModel > earthRotationModelFromUtc =
            boost::make_shared< GcrsToItrsRotationModel >(
                earth_orientation::createStandardEarthOrientationCalculator( ),
                utc_scale );

    // Test Tudat vs. Sofa implementations, with default Sofa EOP corrections (as defined in cookbook
    {
        Eigen::Matrix3d sofaCookbookResult = getSofaEarthOrientationExamples( 3 ).transpose( );
        Eigen::Matrix3d tudatResult = earthRotationModelFromUtc->getRotationToBaseFrame( sofaCookbookTime ).toRotationMatrix( );

        // Check sofa against Tudat result, small difference due to slightly different values of EOP corrections.
        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL( sofaCookbookResult( i, j ) - tudatResult( i, j ), 2.0E-9 );
            }
        }
    }

    // Test Tudat vs. Sofa implementations with identical EOP corrections.
    {
        // Set current time in UTC and TT
        double interpolationUtc = sofaCookbookTime;
        double interpolationTt = earthRotationModelFromUtc->getAnglesCalculator( )->getTerrestrialTimeScaleConverter( )->
                getCurrentTime( utc_scale, tt_scale, interpolationUtc );

        // Get EOP corrections
        double Xcorrection = earthRotationModelFromUtc->getAnglesCalculator( )->getPrecessionNutationCalculator( )->
                getDailyCorrectionInterpolator( )->interpolate( interpolationUtc ).x( );
        double Ycorrection = earthRotationModelFromUtc->getAnglesCalculator( )->getPrecessionNutationCalculator( )->
                getDailyCorrectionInterpolator( )->interpolate( interpolationUtc ).y( );
        double xPolarMotion = earthRotationModelFromUtc->getAnglesCalculator( )->getPolarMotionCalculator( )->
                getPositionOfCipInItrs( interpolationTt, interpolationUtc ).x( );
        double yPolarMotion = earthRotationModelFromUtc->getAnglesCalculator( )->getPolarMotionCalculator( )->
                getPositionOfCipInItrs( interpolationTt, interpolationUtc ).y( );
        double ut1Correction = earthRotationModelFromUtc->getAnglesCalculator( )->getTerrestrialTimeScaleConverter( )->
                getUt1Correction( utc_scale, Time( sofaCookbookTime ) );

        // Compute Sofa rotation matrix
        Eigen::Matrix3d  sofaCookbookResult = getSofaEarthOrientationExamples(
                    3, unit_conversions::convertRadiansToArcSeconds( Xcorrection ) * 1000.0,
                    unit_conversions::convertRadiansToArcSeconds( Ycorrection ) * 1000.0,
                    unit_conversions::convertRadiansToArcSeconds( xPolarMotion ),
                    unit_conversions::convertRadiansToArcSeconds( yPolarMotion ), ut1Correction ).transpose( );

        // Compute Tudat rotation matrix
        Eigen::Matrix3d tudatResult = earthRotationModelFromUtc->getRotationToBaseFrame( sofaCookbookTime ).toRotationMatrix( );

        // Check sofa against Tudat result, small difference due to rounding errors, in particular in Earth rotation angle
        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                if( i < 2 && j < 2 )
                {
                    BOOST_CHECK_SMALL( sofaCookbookResult( i, j ) - tudatResult( i, j ), 1.0E-11 );
                }
                else
                {
                    BOOST_CHECK_SMALL( sofaCookbookResult( i, j ) - tudatResult( i, j ), 1.0E-14 );
                }
            }
        }

        // Compute Tudat rotation matrix with high-precision time input
        long double sofaCookbookExtendedTime = convertCalendarDateToJulianDaysSinceEpoch< long double >(
                    year, month, day, hour, minutes, seconds, JULIAN_DAY_ON_J2000 ) *
                physical_constants::JULIAN_DAY_LONG;
        Eigen::Matrix3d tudatResultPrecise = earthRotationModelFromUtc->getRotationToBaseFrameFromExtendedTime(
                    Time( sofaCookbookExtendedTime ) ).toRotationMatrix( );

        // Check sofa against Tudat result
        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL( sofaCookbookResult( i, j ) - tudatResultPrecise( i, j ), 1.0E-15 );
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}

}
