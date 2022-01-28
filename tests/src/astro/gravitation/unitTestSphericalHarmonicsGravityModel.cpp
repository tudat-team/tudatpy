/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Mathworks. gravitysphericalharmonic, Implement spherical harmonic representation of
 *        planetary gravity, help documentation of Aerospace Toolbox of MATLAB R2012a, 2012.
 *
 *    Notes
 *      In future, more tests should be added here to test the completeness of the functions
 *      implemented. In particular, tests should be added to ascertain the maximum degree and order
 *      to which the functions are still able to produce accelerations. Further, the tests are
 *      currently all based off of data generated with MATLAB (Mathworks, 2012). Ideally, at least
 *      one other source of benchmark data should be included to thoroughly test the code and
 *      minimize the risks of bugs being present. The runtime errors thrown are also not tested;
 *      this behaviour needs to be tested rigorously too.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <cmath>
#include <limits>

#include <boost/lambda/lambda.hpp>
#include <boost/make_shared.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "tudat/basics/testMacros.h"

#include "tudat/astro/gravitation/sphericalHarmonicsGravityModel.h"
#include "tudat/math/basic/sphericalHarmonics.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_SphericalHarmonicsGravity )

// Check single harmonics term of degree = 2 and order = 0.
BOOST_AUTO_TEST_CASE( test_SphericalHarmonicsGravitationalAcceleration_Demo1 )
{
    // Define gravitational parameter of Earth [m^3 s^-2]. The value is obtained from the Earth
    // Gravitational Model 2008 as described by Mathworks [2012].
    const double gravitationalParameter = 3.986004418e14;

    // Define radius of Earth [m]. The value is obtained from the Earth Gravitational Model 2008 as
    // described by Mathworks [2012].
    const double planetaryRadius = 6378137.0;

    // Define degree and order.
    const int degree = 2;
    const int order = 0;

    // Define geodesy-normalized coefficients for degree 2 and order 0. The values are obtained
    // from the Earth Gravitational Model 2008 as described by Mathworks [2012].
    const double cosineCoefficient = -4.841651437908150e-4;
    const double sineCoefficient = 0.0;

    // Define arbitrary Cartesian position [m].
    const Eigen::Vector3d position( 7.0e6, 8.0e6, 9.0e6 );

    // Compute acceleration [m s^-2].
    const Eigen::Vector3d acceleration
            = gravitation::computeSingleGeodesyNormalizedGravitationalAcceleration(
                position,
                gravitationalParameter,
                planetaryRadius,
                degree,
                order,
                cosineCoefficient,
                sineCoefficient,
                std::make_shared< basic_mathematics::SphericalHarmonicsCache >( 3, 1 ) );

    // Define expected acceleration according to the MATLAB function 'gravitysphericalharmonic'
    // described by Mathworks [2012] [m s^-2].
    const Eigen::Vector3d expectedAcceleration(
                3.824456141317033e-4, 4.370807018648038e-4, -4.124819656816540e-4 );

    // Check if expected result matches computed result.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedAcceleration, acceleration, 1.0e-14 );
}

// Check single harmonics term of degree = 2 and order = 1.
BOOST_AUTO_TEST_CASE( test_SphericalHarmonicsGravitationalAcceleration_Demo2 )
{
    // Define gravitational parameter of Earth [m^3 s^-2]. The value is obtained from the Earth
    // Gravitational Model 2008 as described by Mathworks [2012].
    const double gravitationalParameter = 3.986004418e14;

    // Define radius of Earth [m]. The value is obtained from the Earth Gravitational Model 2008 as
    // described by Mathworks [2012].
    const double planetaryRadius = 6378137.0;

    // Define degree and order.
    const int degree = 2;
    const int order = 1;

    // Define geodesy-normalized coefficients for degree 2 and order 1. The values are obtained
    // from the Earth Gravitational Model 2008 as described by Mathworks [2012].
    const double cosineCoefficient = -2.066155090741760e-10;
    const double sineCoefficient = 1.384413891379790e-9;

    // Define arbitrary Cartesian position [m].
    const Eigen::Vector3d position( 7.0e6, 8.0e6, 9.0e6 );

    // Compute acceleration for 2,1 term [m s^-2].
    const Eigen::Vector3d acceleration
            = gravitation::computeSingleGeodesyNormalizedGravitationalAcceleration(
                position,
                gravitationalParameter,
                planetaryRadius,
                degree,
                order,
                cosineCoefficient,
                sineCoefficient,
                std::make_shared< basic_mathematics::SphericalHarmonicsCache >( 3, 2 ) );

    // Define expected acceleration according to the MATLAB function 'gravitysphericalharmonic'
    // described by Mathworks [2012] [m s^-2].
    const Eigen::Vector3d expectedAcceleration(
                -2.095860391422327e-9, -6.479563983539470e-10, -1.254667924094711e-9 );

    // Check if expected result matches computed result.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedAcceleration, acceleration, 2.0e-15 );
}

// Check single harmonics term of degree = 2 and order = 2.
BOOST_AUTO_TEST_CASE( test_SphericalHarmonicsGravitationalAcceleration_Demo3 )
{
    // Define gravitational parameter of Earth [m^3 s^-2]. The value is obtained from the Earth
    // Gravitational Model 2008 as described by Mathworks [2012].
    const double gravitationalParameter = 3.986004418e14;

    // Define radius of Earth [m]. The value is obtained from the Earth Gravitational Model 2008 as
    // described by Mathworks [2012].
    const double planetaryRadius = 6378137.0;

    // Define degree and order.
    const int degree = 2;
    const int order = 2;

    // Define geodesy-normalized coefficients for degree 2 and order 2. The values are obtained
    // from the Earth Gravitational Model 2008 as described by Mathworks [2012].
    const double cosineCoefficient = 2.439383573283130e-6;
    const double sineCoefficient = -1.400273703859340e-6;

    // Define arbitrary Cartesian position [m].
    const Eigen::Vector3d position( 7.0e6, 8.0e6, 9.0e6 );

    // Compute acceleration for 2,2 term [m s^-2].
    const Eigen::Vector3d acceleration
            = gravitation::computeSingleGeodesyNormalizedGravitationalAcceleration(
                position,
                gravitationalParameter,
                planetaryRadius,
                degree,
                order,
                cosineCoefficient,
                sineCoefficient,
                std::make_shared< basic_mathematics::SphericalHarmonicsCache >( 3, 3 ) );

    // Define expected acceleration according to the MATLAB function 'gravitysphericalharmonic'
    // described by Mathworks [2012] [m s^-2].
    const Eigen::Vector3d expectedAcceleration(
                2.793956087356544e-6, -1.123346383523296e-6, 2.687522426265113e-6 );

    // Check if expected result matches computed result.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedAcceleration, acceleration, 1.0e-15 );
}

// Check the sum of all harmonics terms up to degree = 5 and order = 5.
BOOST_AUTO_TEST_CASE( test_SphericalHarmonicsGravitationalAcceleration_Demo4 )
{
    // Define gravitational parameter of Earth [m^3 s^-2]. The value is obtained from the Earth
    // Gravitational Model 2008 as described by Mathworks [2012].
    const double gravitationalParameter = 3.986004418e14;

    // Define radius of Earth [m]. The value is obtained from the Earth Gravitational Model 2008 as
    // described by Mathworks [2012].
    const double planetaryRadius = 6378137.0;

    // Define geodesy-normalized coefficients up to degree 5 and order 5. The values are obtained
    // from the Earth Gravitational Model 2008 as described by Mathworks [2012].
    const Eigen::MatrixXd cosineCoefficients =
            ( Eigen::MatrixXd( 6, 6 ) <<
              1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              -4.841651437908150e-4, -2.066155090741760e-10, 2.439383573283130e-6, 0.0, 0.0, 0.0,
              9.571612070934730e-7, 2.030462010478640e-6, 9.047878948095281e-7,
              7.213217571215680e-7, 0.0, 0.0, 5.399658666389910e-7, -5.361573893888670e-7,
              3.505016239626490e-7, 9.908567666723210e-7, -1.885196330230330e-7, 0.0,
              6.867029137366810e-8, -6.292119230425290e-8, 6.520780431761640e-7,
              -4.518471523288430e-7, -2.953287611756290e-7, 1.748117954960020e-7
              ).finished( );

    const Eigen::MatrixXd sineCoefficients =
            ( Eigen::MatrixXd( 6, 6 ) <<
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 1.384413891379790e-9, -1.400273703859340e-6, 0.0, 0.0, 0.0,
              0.0, 2.482004158568720e-7, -6.190054751776180e-7, 1.414349261929410e-6, 0.0, 0.0,
              0.0, -4.735673465180860e-7, 6.624800262758290e-7, -2.009567235674520e-7,
              3.088038821491940e-7, 0.0, 0.0, -9.436980733957690e-8, -3.233531925405220e-7,
              -2.149554083060460e-7, 4.980705501023510e-8, -6.693799351801650e-7
              ).finished( );

    // Define arbitrary Cartesian position [m].
    const Eigen::Vector3d position( 7.0e6, 8.0e6, 9.0e6 );

    // Compute resultant acceleration [m s^-2].
    std::map< std::pair< int, int >, Eigen::Vector3d > dummyMap;
    const Eigen::Vector3d acceleration
            = gravitation::computeGeodesyNormalizedGravitationalAccelerationSum(
                position,
                gravitationalParameter,
                planetaryRadius,
                cosineCoefficients,
                sineCoefficients,
                std::make_shared< basic_mathematics::SphericalHarmonicsCache >( 6, 6 ),
                dummyMap );

    // Define expected acceleration according to the MATLAB function 'gravitysphericalharmonic'
    // described by Mathworks [2012] [m s^-2].
    const Eigen::Vector3d expectedAcceleration(
                -1.032215878106932, -1.179683946769393, -1.328040277155269 );

    // Check if expected result matches computed result.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedAcceleration, acceleration, 1.0e-15 );
}

// Check the sum of all harmonics terms up to degree = 5 and order = 5 using the wrapper class.
BOOST_AUTO_TEST_CASE( test_SphericalHarmonicsGravitationalAccelerationWrapperClass )
{
    // Short-cuts.
    using namespace gravitation;

    // Define gravitational parameter of Earth [m^3 s^-2]. The value is obtained from the Earth
    // Gravitational Model 2008 as described by Mathworks [2012].
    const double gravitationalParameter = 3.986004418e14;

    // Define radius of Earth [m]. The value is obtained from the Earth Gravitational Model 2008 as
    // described by Mathworks [2012].
    const double planetaryRadius = 6378137.0;

    // Define geodesy-normalized coefficients up to degree 5 and order 5. The values are obtained
    // from the Earth Gravitational Model 2008 as described by Mathworks [2012].
    const Eigen::MatrixXd cosineCoefficients =
            ( Eigen::MatrixXd( 6, 6 ) <<
              1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              -4.841651437908150e-4, -2.066155090741760e-10, 2.439383573283130e-6, 0.0, 0.0, 0.0,
              9.571612070934730e-7, 2.030462010478640e-6, 9.047878948095281e-7,
              7.213217571215680e-7, 0.0, 0.0, 5.399658666389910e-7, -5.361573893888670e-7,
              3.505016239626490e-7, 9.908567666723210e-7, -1.885196330230330e-7, 0.0,
              6.867029137366810e-8, -6.292119230425290e-8, 6.520780431761640e-7,
              -4.518471523288430e-7, -2.953287611756290e-7, 1.748117954960020e-7
              ).finished( );

    const Eigen::MatrixXd sineCoefficients =
            ( Eigen::MatrixXd( 6, 6 ) <<
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 1.384413891379790e-9, -1.400273703859340e-6, 0.0, 0.0, 0.0,
              0.0, 2.482004158568720e-7, -6.190054751776180e-7, 1.414349261929410e-6, 0.0, 0.0,
              0.0, -4.735673465180860e-7, 6.624800262758290e-7, -2.009567235674520e-7,
              3.088038821491940e-7, 0.0, 0.0, -9.436980733957690e-8, -3.233531925405220e-7,
              -2.149554083060460e-7, 4.980705501023510e-8, -6.693799351801650e-7
              ).finished( );

    // Define arbitrary Cartesian position [m].
    const Eigen::Vector3d position( 7.0e6, 8.0e6, 9.0e6 );

    // Declare spherical harmonics gravitational acceleration class object.
    SphericalHarmonicsGravitationalAccelerationModelPointer earthGravity
            = std::make_shared< SphericalHarmonicsGravitationalAccelerationModel >(
                [ & ]( Eigen::Vector3d& input ){ input = position; }, gravitationalParameter, planetaryRadius,
                cosineCoefficients, sineCoefficients );
    earthGravity->updateMembers( );

    // Compute resultant acceleration [m s^-2].
    const Eigen::Vector3d acceleration = earthGravity->getAcceleration( );

    // Define expected acceleration according to the MATLAB function 'gravitysphericalharmonic'
    // described by Mathworks [2012] [m s^-2].
    const Eigen::Vector3d expectedAcceleration(
                -1.032215878106932, -1.179683946769393, -1.328040277155269 );

    // Check if expected result matches computed result.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedAcceleration, acceleration, 1.0e-15 );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
