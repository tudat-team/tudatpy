#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>
#include <iostream>

#include "tudat/io/basicInputOutput.h"
#include "tudat/astro/ground_stations/oceanTideEarthDeformation.h"
#include "tudat/astro/basic_astro/dateTime.h"
#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"

namespace tudat
{
namespace unit_tests
{

using namespace basic_astrodynamics;
using namespace earth_orientation;

BOOST_AUTO_TEST_SUITE( test_ocean_tide_site_displacement )

BOOST_AUTO_TEST_CASE( testOceanTideSiteDisplacement )
{

    std::string testFile = paths::getEarthDeformationDataFilesPath( ) + "/oceanTidetest.blq";

    OceanTideEarthDeformation oceanTideDeformationModel = OceanTideEarthDeformation( testFile );

    double timeStep = 3600.0;
    double testTime = timeFromDecomposedDateTime<double>( 2009, 06, 25, 01, 10, 45 );

    Eigen::Matrix<double, 24, 3> onsalaResults;
    Eigen::Matrix<double, 24, 3> icelandResults;

    for ( int i = 0; i < 24; i++ )
    {
        onsalaResults.block( i, 0, 1, 3 ) =
            oceanTideDeformationModel.calculateDisplacementInEnuFrame( testTime, "ONSALA" ).transpose( );
        icelandResults.block( i, 0, 1, 3 ) =
            oceanTideDeformationModel.calculateDisplacementInEnuFrame( testTime, "REYKJAVIK" ).transpose( );
        testTime += timeStep;
    }

    Eigen::Matrix<double, 24, 3> onsalaRawTestMatrix;
    onsalaRawTestMatrix << 0.003094, -0.001538, -0.000895,
        0.001812, -0.000950, -0.000193,
        0.000218, -0.000248, 0.000421,
        -0.001104, 0.000404, 0.000741,
        -0.001668, 0.000863, 0.000646,
        -0.001209, 0.001042, 0.000137,
        0.000235, 0.000926, -0.000667,
        0.002337, 0.000580, -0.001555,
        0.004554, 0.000125, -0.002278,
        0.006271, -0.000291, -0.002615,
        0.006955, -0.000537, -0.002430,
        0.006299, -0.000526, -0.001706,
        0.004305, -0.000244, -0.000559,
        0.001294, 0.000245, 0.000793,
        -0.002163, 0.000819, 0.002075,
        -0.005375, 0.001326, 0.003024,
        -0.007695, 0.001622, 0.003448,
        -0.008669, 0.001610, 0.003272,
        -0.008143, 0.001262, 0.002557,
        -0.006290, 0.000633, 0.001477,
        -0.003566, -0.000155, 0.000282,
        -0.000593, -0.000941, -0.000766,
        0.001992, -0.001561, -0.001457,
        0.003689, -0.001889, -0.001680;


    Eigen::Matrix<double, 24, 3> icelandRawTestMatrix;
    icelandRawTestMatrix << -0.005940, -0.001245, -0.000278,
        0.013516, -0.001086, 0.003212,
        0.029599, -0.000353, 0.005483,
        0.038468, 0.000699, 0.005997,
        0.038098, 0.001721, 0.004690,
        0.028780, 0.002363, 0.001974,
        0.013016, 0.002371, -0.001369,
        -0.005124, 0.001653, -0.004390,
        -0.021047, 0.000310, -0.006225,
        -0.030799, -0.001383, -0.006313,
        -0.032056, -0.003048, -0.004549,
        -0.024698, -0.004288, -0.001314,
        -0.010814, -0.004794, 0.002623,
        0.005849, -0.004416, 0.006291,
        0.020857, -0.003208, 0.008766,
        0.030226, -0.001413, 0.009402,
        0.031437, 0.000594, 0.007996,
        0.024079, 0.002389, 0.004844,
        0.009945, 0.003606, 0.000663,
        -0.007426, 0.004022, -0.003581,
        -0.023652, 0.003601, -0.006911,
        -0.034618, 0.002505, -0.008585,
        -0.037515, 0.001044, -0.008270,
        -0.031544, -0.000402, -0.006125;

    Eigen::Matrix<double, 24, 3> onsalaTestMatrix = onsalaRawTestMatrix;
    onsalaTestMatrix.block(0, 0, 24, 1 ) = -onsalaRawTestMatrix.block(0, 2, 24, 1 );
    onsalaTestMatrix.block(0, 1, 24, 1 ) = -onsalaRawTestMatrix.block(0, 1, 24, 1 );
    onsalaTestMatrix.block(0, 2, 24, 1 ) = onsalaRawTestMatrix.block(0, 0, 24, 1 );

    std::cout<<onsalaTestMatrix<<std::endl<<std::endl<<onsalaResults<<std::endl;

    Eigen::Matrix<double, 24, 3> icelandTestMatrix = icelandRawTestMatrix;
//
    icelandTestMatrix.block(0, 0, 24, 1 ) = -icelandRawTestMatrix.block(0, 2, 24, 1 );
    icelandTestMatrix.block(0, 1, 24, 1 ) = -icelandRawTestMatrix.block(0, 1, 24, 1 );
    icelandTestMatrix.block(0, 2, 24, 1 ) = icelandRawTestMatrix.block(0, 0, 24, 1 );

    BOOST_CHECK_SMALL(( onsalaResults - onsalaTestMatrix ).block( 0, 0, 24, 1 ).cwiseAbs( ).maxCoeff( ), 0.25E-3 );
    BOOST_CHECK_SMALL(( onsalaResults - onsalaTestMatrix ).block( 0, 1, 24, 1 ).cwiseAbs( ).maxCoeff( ), 0.25E-3 );
    BOOST_CHECK_SMALL(( onsalaResults - onsalaTestMatrix ).block( 0, 2, 24, 1 ).cwiseAbs( ).maxCoeff( ), 0.75E-3 );

    BOOST_CHECK_SMALL(( icelandResults - icelandTestMatrix ).block( 0, 0, 24, 1 ).cwiseAbs( ).maxCoeff( ), 1.0E-3 );
    BOOST_CHECK_SMALL(( icelandResults - icelandTestMatrix ).block( 0, 1, 24, 1 ).cwiseAbs( ).maxCoeff( ), 1.0E-3 );
    BOOST_CHECK_SMALL(( icelandResults - icelandTestMatrix ).block( 0, 2, 24, 1 ).cwiseAbs( ).maxCoeff( ), 4.0E-3 );
}

BOOST_AUTO_TEST_SUITE_END( )

}

}
