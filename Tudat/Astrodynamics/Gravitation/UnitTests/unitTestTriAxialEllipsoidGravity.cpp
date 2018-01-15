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

#include <boost/test/unit_test.hpp>

#include "Tudat/Basics/testMacros.h"
#include "Tudat/Astrodynamics/Gravitation/triAxialEllipsoidGravity.h"

namespace tudat
{

namespace unit_tests
{

using namespace tudat::gravitation;


//! Function to get theoretical values of zonal gravity field coefficient of ellipsoid of revolution
double getZonalTermForEllipsoidOfRevolution(
        const double axisA, const double axisC, const double referenceRadius, const int degree )
{
    double zonalTerm = 0.0;
    if( degree % 2 == 0 )
    {
        zonalTerm = 3.0 / std::pow( referenceRadius, degree ) *
                std::pow( axisA * axisA - axisC * axisC, degree / 2 ) *
                std::pow( -1.0, degree / 2 ) /
                ( ( static_cast< double >( degree ) + 1.0 ) *
                  ( static_cast< double >( degree ) + 3.0 ) );
    }

    return zonalTerm;
}

BOOST_AUTO_TEST_SUITE( test_tri_axial_ellipsoid_gravity )

//! Test gravity field coefficients of homogeneous triaxial ellipsoid
BOOST_AUTO_TEST_CASE( testTriAxialEllipsoidGravity )
{
    double axisA = 26.8E3;
    double axisB = 22.4E3;
    double axisC = 18.4E3;

    // Compute coefficients
    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > sphericalHarmonicCoefficients =
            createTriAxialEllipsoidSphericalHarmonicCoefficients( axisA, axisB, axisC, 4, 4 );
    double referenceRadius = calculateTriAxialEllipsoidReferenceRadius( axisA, axisB, axisC );

    // Compute expected low-degree coefficients
    double expectedC20 = 1.0 / ( 5.0 * referenceRadius * referenceRadius ) *
            ( axisC * axisC - ( axisA * axisA + axisB * axisB ) / 2.0 );
    double expectedC22 = 1.0 / ( 20.0 * referenceRadius * referenceRadius ) *
            ( axisA * axisA - axisB * axisB );
    double expectedC40 = 15.0 / 7.0 * ( expectedC20 * expectedC20 +
                                        2.0 * expectedC22 * expectedC22 );
    double expectedC42 = 5.0 / 7.0 * ( expectedC20 * expectedC22 );
    double expectedC44 = 5.0 / 28.0 * ( expectedC22 * expectedC22 );

    BOOST_CHECK_CLOSE_FRACTION(
                expectedC20, sphericalHarmonicCoefficients.first( 2, 0 ),
                2.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION(
                expectedC22, sphericalHarmonicCoefficients.first( 2, 2 ),
                2.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION(
                expectedC40, sphericalHarmonicCoefficients.first( 4, 0 ),
                2.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION(
                expectedC42, sphericalHarmonicCoefficients.first( 4, 2 ),
                2.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION(
                expectedC44, sphericalHarmonicCoefficients.first( 4, 4 ),
                2.0 * std::numeric_limits< double >::epsilon( ) );

    // Check whether C-coefficients with odd degree or order, and all S-coefficients, are zero.
    for( unsigned int i = 0; i < 5; i++ )
    {
        for( unsigned int j = 0; j < 5; j++ )
        {
            if( ( i % 2 != 0 ) || ( j % 2 ) != 0 )
            {
                BOOST_CHECK_EQUAL( sphericalHarmonicCoefficients.first( i, j ), 0.0 );
            }
            BOOST_CHECK_EQUAL( sphericalHarmonicCoefficients.second( i, j ), 0.0 );

        }
    }

    // Compute gravity field coefficients for ellipsoid of revolution
    axisB = axisA;
    referenceRadius = calculateTriAxialEllipsoidReferenceRadius( axisA, axisB, axisC );
    sphericalHarmonicCoefficients = createTriAxialEllipsoidSphericalHarmonicCoefficients( axisA, axisB, axisC, 20, 1 );

    // Check computation against theoretical values.
    for( int i = 2; i < 21; i += 2 )
    {
        BOOST_CHECK_CLOSE_FRACTION(
                    sphericalHarmonicCoefficients.first( i, 0 ),
                    getZonalTermForEllipsoidOfRevolution( axisA, axisC, referenceRadius, i ),
                    2.0 * std::numeric_limits< double >::epsilon( ) );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}

}

