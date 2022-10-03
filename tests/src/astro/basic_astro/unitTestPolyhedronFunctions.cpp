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

#include "tudat/basics/testMacros.h"

#include "tudat/astro/basic_astro/polyhedronFuntions.h"
#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_polyhedron_functions )

//! Test useful functions to process polyhedron shape.
BOOST_AUTO_TEST_CASE( testPolyhedronUtilities )
{
    // Define cuboid polyhedron
    const double w = 10.0; // width
    const double h = 10.0; // height
    const double l = 20.0; // length

    Eigen::MatrixXd verticesCoordinates(8,3);
    Eigen::MatrixXi verticesDefiningEachFacet(12,3);
    verticesCoordinates <<
        0.0, 0.0, 0.0,
        l, 0.0, 0.0,
        0.0, w, 0.0,
        l, w, 0.0,
        0.0, 0.0, h,
        l, 0.0, h,
        0.0, w, h,
        l, w, h;
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

    // Define tolerance
    const double tolerance = 1e-15;

    // Computation of surface area
    double expectedArea = 2 * l * h + 2 * l * w + 2 * w * h;
    double computedArea = basic_astrodynamics::computePolyhedronSurfaceArea( verticesCoordinates,
                                                                             verticesDefiningEachFacet );
    BOOST_CHECK_CLOSE_FRACTION( expectedArea, computedArea, tolerance );

    // Computation of volume
    double expectedVolume = l * w * h;
    double computedVolume = basic_astrodynamics::computePolyhedronVolume( verticesCoordinates,
                                                                          verticesDefiningEachFacet );
    BOOST_CHECK_CLOSE_FRACTION( expectedVolume, computedVolume, tolerance );

    // Mean radius
    double expectedRadius = 7.815926417967719; // Computed by hand
    double computedRadius = basic_astrodynamics::computePolyhedronMeanRadius( verticesCoordinates,
                                                                              verticesDefiningEachFacet );
    BOOST_CHECK_CLOSE_FRACTION( expectedRadius, computedRadius, tolerance );

    // Computation of centroid
    Eigen::Vector3d expectedCentroid = (Eigen::Vector3d() << l / 2.0, w / 2.0, h / 2.0).finished();
    Eigen::Vector3d computedCentroid = basic_astrodynamics::computePolyhedronCentroidPosition( verticesCoordinates,
                                                                                               verticesDefiningEachFacet );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedCentroid, computedCentroid, tolerance );

    // Correction of centroid
    Eigen::Vector3d desiredCentroid = (Eigen::Vector3d() << 3.0, 4.0, -5.0).finished();
    Eigen::MatrixXd correctedVerticesCoordinates = basic_astrodynamics::modifyPolyhedronCentroidPosition(
            verticesCoordinates, verticesDefiningEachFacet, desiredCentroid );
    computedCentroid = basic_astrodynamics::computePolyhedronCentroidPosition( correctedVerticesCoordinates,
                                                                               verticesDefiningEachFacet );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( desiredCentroid, computedCentroid, tolerance );

    // Computation of inertia tensor
    double density = 2.0;
    double mass = density * expectedVolume;
    // Inertia tensor for a cuboid wrt principal axes of inertia
    // http://www2.ece.ohio-state.edu/~zhang/RoboticsClass/docs/LN11_RigidBodyDynamics.pdf
    Eigen::Matrix3d inertiaTensorWrtPrincipalAxes = (Eigen::Matrix3d() <<
            mass / 12.0 * ( std::pow(w, 2) + std::pow(h, 2) ), 0.0, 0.0,
            0.0, mass / 12.0 * ( std::pow(l, 2) + std::pow(h, 2) ), 0.0,
            0.0, 0.0, mass / 12.0 * ( std::pow(l, 2) + std::pow(w, 2) )).finished();
    // Inertia tensor for the cuboid, using parallel axis theorem (Dobrovolskis, 1996)
    Eigen::Matrix3d expectedInertiaTensor = inertiaTensorWrtPrincipalAxes + mass * ( Eigen::Matrix3d() <<
            std::pow(expectedCentroid(1), 2) + std::pow(expectedCentroid(2), 2),
            - expectedCentroid(0) * expectedCentroid(1),
            - expectedCentroid(0) * expectedCentroid(2),
            - expectedCentroid(0) * expectedCentroid(1),
            std::pow(expectedCentroid(0), 2) + std::pow(expectedCentroid(2), 2),
            - expectedCentroid(1) * expectedCentroid(2),
            - expectedCentroid(0) * expectedCentroid(2),
            - expectedCentroid(1) * expectedCentroid(2),
            std::pow(expectedCentroid(0), 2) + std::pow(expectedCentroid(1), 2) ).finished();

    Eigen::Matrix3d computedInertiaTensor = basic_astrodynamics::computePolyhedronInertiaTensor(
            verticesCoordinates, verticesDefiningEachFacet, density );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedInertiaTensor, computedInertiaTensor, tolerance );

}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
