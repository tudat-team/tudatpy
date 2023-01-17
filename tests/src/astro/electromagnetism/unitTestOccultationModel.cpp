/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <iostream>

#include <boost/test/unit_test.hpp>

#include "tudat/astro/electromagnetism/occultationModel.h"


namespace tudat
{
namespace unit_tests
{

using namespace tudat::electromagnetism;

BOOST_AUTO_TEST_SUITE(test_luminosity_model)

// No need to test shadow functions since they are already tested in mission geometry tests

// Test visibility with simple geometric cases
BOOST_AUTO_TEST_CASE( testEvaluateVisibilityWithOccultation )
{
    // Target behind occulting body
    BOOST_CHECK(!evaluateVisibilityWithOccultation(
            Eigen::Vector3d(0, 0, 0),
            Eigen::Vector3d(2, 0, 0),
            1,
            Eigen::Vector3d(3, 1, 0)
        ));

    // Target behind occulting body
    BOOST_CHECK(!evaluateVisibilityWithOccultation(
            Eigen::Vector3d(0, 0, 0),
            Eigen::Vector3d(2, 0, 0),
            1,
            Eigen::Vector3d(3, 1, 1)
        ));

    // Target behind occulting body but to the side
    BOOST_CHECK(evaluateVisibilityWithOccultation(
            Eigen::Vector3d(0, 0, 0),
            Eigen::Vector3d(2, 0, 0),
            1,
            Eigen::Vector3d(3, 2, 0)
        ));

    // Target behind occulting body but to the side
    BOOST_CHECK(evaluateVisibilityWithOccultation(
            Eigen::Vector3d(0, 0, 0),
            Eigen::Vector3d(2, 0, 0),
            1,
            Eigen::Vector3d(3, 0, 2)
        ));

    // Target in front of occulting body
    BOOST_CHECK(evaluateVisibilityWithOccultation(
            Eigen::Vector3d(0, 0, 0),
            Eigen::Vector3d(2, 0, 0),
            1,
            Eigen::Vector3d(0.5, 0, 0)
        ));

    // Target inside occulting body but in front of its center
    // (unrealistic but should be visible by convention)
    BOOST_CHECK(evaluateVisibilityWithOccultation(
            Eigen::Vector3d(0, 0, 0),
            Eigen::Vector3d(2, 0, 0),
            1,
            Eigen::Vector3d(1.5, 0, 0)
        ));

    // Target inside occulting body but behind its center
    // (unrealistic but should not be visible by convention)
    BOOST_CHECK(!evaluateVisibilityWithOccultation(
            Eigen::Vector3d(0, 0, 0),
            Eigen::Vector3d(2, 0, 0),
            1,
            Eigen::Vector3d(2.5, 0, 0)
        ));
}

// Test inverse cases of previous test to check commutativity
// (if the target can see the source, the source can see the target)
BOOST_AUTO_TEST_CASE( testEvaluateVisibilityWithOccultation_Commutativity )
{
    // Target behind occulting body
    BOOST_CHECK(!evaluateVisibilityWithOccultation(
            Eigen::Vector3d(3, 1, 0),
            Eigen::Vector3d(2, 0, 0),
            1,
            Eigen::Vector3d(0, 0, 0)
        ));

    // Target behind occulting body
    BOOST_CHECK(!evaluateVisibilityWithOccultation(
            Eigen::Vector3d(3, 1, 1),
            Eigen::Vector3d(2, 0, 0),
            1,
            Eigen::Vector3d(0, 0, 0)
        ));

    // Target behind occulting body but to the side
    BOOST_CHECK(evaluateVisibilityWithOccultation(
            Eigen::Vector3d(3, 2, 0),
            Eigen::Vector3d(2, 0, 0),
            1,
            Eigen::Vector3d(0, 0, 0)
        ));

    // Target behind occulting body but to the side
    BOOST_CHECK(evaluateVisibilityWithOccultation(
            Eigen::Vector3d(3, 0, 2),
            Eigen::Vector3d(2, 0, 0),
            1,
            Eigen::Vector3d(0, 0, 0)
        ));

    // Target in front of occulting body
    BOOST_CHECK(evaluateVisibilityWithOccultation(
            Eigen::Vector3d(0.5, 0, 0),
            Eigen::Vector3d(2, 0, 0),
            1,
            Eigen::Vector3d(0, 0, 0)
        ));

    // Target inside occulting body but in front of its center
    // (unrealistic but should be visible by convention)
    BOOST_CHECK(evaluateVisibilityWithOccultation(
            Eigen::Vector3d(1.5, 0, 0),
            Eigen::Vector3d(2, 0, 0),
            1,
            Eigen::Vector3d(0, 0, 0)
        ));

    // Target inside occulting body but behind its center
    // (unrealistic but should not be visible by convention)
    BOOST_CHECK(!evaluateVisibilityWithOccultation(
            Eigen::Vector3d(2.5, 0, 0),
            Eigen::Vector3d(2, 0, 0),
            1,
            Eigen::Vector3d(0, 0, 0)
        ));
}

BOOST_AUTO_TEST_SUITE_END()

} // namespace unit_tests
} // namespace tudat
