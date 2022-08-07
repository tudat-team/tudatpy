/*    Copyright (c) 2010-2022, Delft University of Technology
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

#include <limits>

#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "tudat/astro/basic_astro/polyhedronFuntions.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_polyhedron_functions )

//! Test useful functions to process polyhedron shape.
BOOST_AUTO_TEST_CASE( testPolyhedronUtilities )
{
    // Define cuboid polyhedron
    Eigen::MatrixXd verticesCoordinates(8,3);
    Eigen::MatrixXi verticesDefiningEachFacet(12,3);
    verticesCoordinates <<
        0.0, 0.0, 0.0,
        20.0, 0.0, 0.0,
        0.0, 10.0, 0.0,
        20.0, 10.0, 0.0,
        0.0, 0.0, 10.0,
        20.0, 0.0, 10.0,
        0.0, 10.0, 10.0,
        20.0, 10.0, 10.0;
    verticesDefiningEachFacet <<
        2, 1, 0,
        1, 2, 3,
        4, 2, 0,
        2, 4, 6,
        1, 4, 0,
        4, 1, 5,
        6, 5, 7,
        5, 6, 4,
        3, 6, 7,
        6, 3, 2,
        5, 3, 7,
        3, 5, 1;

    // Computation of volume
    const double expectedVolume = 20.0 * 10.0 * 10.0;
    const double computedVolume = polyhedron_utilities::computeVolume(verticesCoordinates, verticesDefiningEachFacet);
    BOOST_CHECK_CLOSE_FRACTION( expectedVolume, computedVolume, std::numeric_limits< double >::epsilon( ) );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
