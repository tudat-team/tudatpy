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
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/unitConversions.h"
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

//! Test if second-degree zonal surface property distribution is zonal
BOOST_AUTO_TEST_CASE( testSecondDegreeZonalPeriodicSurfacePropertyDistribution_Zonality )
{
    spice_interface::loadStandardSpiceKernels();

    const double referenceEpoch =
            spice_interface::convertDateStringToEphemerisTime("1981 DEC 22");
    SecondDegreeZonalPeriodicSurfacePropertyDistribution distributionModel(
            0.34, 0, 0.1, 0, 0.29, referenceEpoch, physical_constants::JULIAN_YEAR_IN_DAYS);
    distributionModel.updateMembers(spice_interface::convertDateStringToEphemerisTime("2005 AUG 19 13:46:17"));

    for (double latitude : {-PI/2, 0., 1.403, PI/2})
    {
        std::vector<double> values;
        for (int longitude = -180; longitude <= 180; longitude += 40)
        {
            values.push_back(distributionModel.getValue(
                    latitude, unit_conversions::convertDegreesToRadians(longitude)));
        }

        // Check that all calculated values are identical
        for (auto& v : values)
        {
            BOOST_CHECK_CLOSE_FRACTION(v, values.front(), 1e-15);
        }
    }
}

//! Test second-degree zonal surface property distribution using albedo from Knocke (1988) by comparison with Python implementation
// https://github.com/DominikStiller/tudelft-hpb-project/blob/ed7637750bd13ac584f37bcd3008b9731044313f/analysis/earth_models.ipynb
BOOST_AUTO_TEST_CASE( testSecondDegreeZonalPeriodicSurfacePropertyDistribution_Albedo )
{
    spice_interface::loadStandardSpiceKernels();

    const double referenceEpoch =
            spice_interface::convertDateStringToEphemerisTime("1981 DEC 22");
    SecondDegreeZonalPeriodicSurfacePropertyDistribution distributionModel(
            0.34, 0, 0.1, 0, 0.29, referenceEpoch, physical_constants::JULIAN_YEAR_IN_DAYS);

    distributionModel.updateMembers(spice_interface::convertDateStringToEphemerisTime("2005 AUG 19 13:46:17"));
    double actualValue = distributionModel.getValue(unit_conversions::convertDegreesToRadians(29.73));
    BOOST_CHECK_CLOSE(actualValue, 0.2752338314886392, 1e-13);

    distributionModel.updateMembers(spice_interface::convertDateStringToEphemerisTime("2012 DEC 21 17:26:17"));
    actualValue = distributionModel.getValue(unit_conversions::convertDegreesToRadians(81.43));
    BOOST_CHECK_CLOSE(actualValue, 0.7192237282075249, 1e-13);

    distributionModel.updateMembers(spice_interface::convertDateStringToEphemerisTime("2022 APR 28 22:49:57"));
    actualValue = distributionModel.getValue(unit_conversions::convertDegreesToRadians(0));
    BOOST_CHECK_CLOSE(actualValue, 0.19500000000000003, 1e-13);

    distributionModel.updateMembers(spice_interface::convertDateStringToEphemerisTime("1977 JAN 08 01:46:13"));
    actualValue = distributionModel.getValue(unit_conversions::convertDegreesToRadians(90.));
    BOOST_CHECK_CLOSE(actualValue, 0.725592271381903, 1e-13);
}

//! Test second-degree zonal surface property distribution using emissivity from Knocke (1988) by comparison with Python implementation
// https://github.com/DominikStiller/tudelft-hpb-project/blob/ed7637750bd13ac584f37bcd3008b9731044313f/analysis/earth_models.ipynb
BOOST_AUTO_TEST_CASE( testSecondDegreeZonalPeriodicSurfacePropertyDistribution_Emissivity )
{
    spice_interface::loadStandardSpiceKernels();

    const double referenceEpoch =
            spice_interface::convertDateStringToEphemerisTime("1981 DEC 22");
    SecondDegreeZonalPeriodicSurfacePropertyDistribution distributionModel(
            0.68, 0, -0.07, 0, -0.18, referenceEpoch, physical_constants::JULIAN_YEAR_IN_DAYS);

    distributionModel.updateMembers(spice_interface::convertDateStringToEphemerisTime("2005 AUG 19 13:46:17"));
    double actualValue = distributionModel.getValue(unit_conversions::convertDegreesToRadians(29.73));
    BOOST_CHECK_CLOSE(actualValue, 0.7223209069278839, 1e-13);

    distributionModel.updateMembers(spice_interface::convertDateStringToEphemerisTime("2012 DEC 21 17:26:17"));
    actualValue = distributionModel.getValue(unit_conversions::convertDegreesToRadians(81.43));
    BOOST_CHECK_CLOSE(actualValue, 0.4367772746818569, 1e-13);

    distributionModel.updateMembers(spice_interface::convertDateStringToEphemerisTime("2022 APR 28 22:49:57"));
    actualValue = distributionModel.getValue(unit_conversions::convertDegreesToRadians(0));
    BOOST_CHECK_CLOSE(actualValue, 0.77, 1e-13);

    distributionModel.updateMembers(spice_interface::convertDateStringToEphemerisTime("1977 JAN 08 01:46:13"));
    actualValue = distributionModel.getValue(unit_conversions::convertDegreesToRadians(90.));
    BOOST_CHECK_CLOSE(actualValue, 0.43308541003266804, 1e-13);
}


BOOST_AUTO_TEST_SUITE_END()

}
}
