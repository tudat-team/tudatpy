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

#include <memory>

#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "tudat/basics/testMacros.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/astro/electromagnetism/surfacePropertyDistribution.h"


namespace tudat
{
namespace unit_tests
{

using namespace tudat::electromagnetism;

using mathematical_constants::PI;

BOOST_AUTO_TEST_SUITE(test_surface_property_distribution)

//! Test constant surface property distribution
BOOST_AUTO_TEST_CASE( testConstantSurfacePropertyDistribution )
{
    const auto expectedValue = 42.0;

    ConstantSurfacePropertyDistribution distributionModel(expectedValue);
    distributionModel.updateMembers(TUDAT_NAN);

    const auto actualValue = distributionModel.getValue(PI / 4, PI / 2);

    BOOST_CHECK_CLOSE_FRACTION(actualValue, expectedValue, 1.0e-15);
}

//! Test spherical harmonics surface property distribution
BOOST_AUTO_TEST_CASE( testSphericalHarmonicsSurfacePropertyDistribution )
{

    // Order 0, degree 0 (identical to constant distribution)
    {
        const auto expectedValue = 42.0;
        const auto expectedMaximumDegree = 0;
        const auto expectedMaximumOrder = 0;

        const Eigen::MatrixXd cosineCoefficients = expectedValue * Eigen::MatrixXd::Identity(1, 1);
        const Eigen::MatrixXd sineCoefficients = Eigen::MatrixXd::Identity(1, 1);

        SphericalHarmonicsSurfacePropertyDistribution distributionModel(
                cosineCoefficients, sineCoefficients);
        distributionModel.updateMembers(TUDAT_NAN);

        const auto actualMaximumDegree = distributionModel.getMaximumDegree();
        const auto actualMaximumOrder = distributionModel.getMaximumOrder();

        BOOST_CHECK(actualMaximumDegree == expectedMaximumDegree);
        BOOST_CHECK(actualMaximumOrder == expectedMaximumOrder);

        auto actualValue = distributionModel.getValue(0, PI / 4);
        BOOST_CHECK_CLOSE_FRACTION(actualValue, expectedValue, 1.0e-15);

        actualValue = distributionModel.getValue(PI / 2, PI / 2);
        BOOST_CHECK_CLOSE_FRACTION(actualValue, expectedValue, 1.0e-15);
    }

    // Order 1, degree 1
    {
        const auto expectedMaximumDegree = 1;
        const auto expectedMaximumOrder = 1;

        Eigen::MatrixXd cosineCoefficients(2, 2);
        cosineCoefficients << 0.1, 0.0,
                              0.2, 0.3;

        Eigen::MatrixXd sineCoefficients(2, 2);
        sineCoefficients << 0.0, 0.0,
                            0.0, 0.4;

        SphericalHarmonicsSurfacePropertyDistribution distributionModel(
                cosineCoefficients, sineCoefficients);
        distributionModel.updateMembers(TUDAT_NAN);

        const auto actualMaximumDegree = distributionModel.getMaximumDegree();
        const auto actualMaximumOrder = distributionModel.getMaximumOrder();

        BOOST_CHECK(actualMaximumDegree == expectedMaximumDegree);
        BOOST_CHECK(actualMaximumOrder == expectedMaximumOrder);

        // Expected values calculated manually
        auto expectedValue = 0.1 + (0.3 + 0.4) / std::sqrt(2);
        auto actualValue = distributionModel.getValue(0, PI / 4);
        BOOST_CHECK_CLOSE_FRACTION(actualValue, expectedValue, 1.0e-15);

        expectedValue = 0.1 + 0.2;
        actualValue = distributionModel.getValue(PI / 2, PI / 2);
        BOOST_CHECK_CLOSE_FRACTION(actualValue, expectedValue, 1.0e-15);
    }
}


BOOST_AUTO_TEST_SUITE_END()

}
}
