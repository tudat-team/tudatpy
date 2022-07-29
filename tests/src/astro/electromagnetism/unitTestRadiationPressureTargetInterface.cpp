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
 *      Ganeff, M.I. Solar radiation pressure benchmark data script,
 *          solarRadiationPressureBenchmarkData.m, available at http://tudat.tudelft.nl, 2012.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <iostream>

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/astro/electromagnetism/radiationPressureTargetInterface.h"


namespace tudat
{
namespace unit_tests
{

using namespace tudat::electromagnetism;

BOOST_AUTO_TEST_SUITE(test_radiation_pressure_target_interface)

//! Check force due to solar constant based on data from (Montenbruck, 2000)
BOOST_AUTO_TEST_CASE( testCannonballRadiationPressureTargetInterfaceAtAU )
{
    // From (Montenbruck, 2000) Eq. (3.69)
    const auto expectedForce = Eigen::Vector3d(4.56e-6, 0, 0);

    CannonballRadiationPressureTargetInterface targetInterface(1, 1);

    const auto actualForce = targetInterface.evaluateRadiationPressureForce(
            1367,
            Eigen::Vector3d::UnitX());

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(actualForce, expectedForce, 0.01);
}

//! Check force with benchmark data is obtained from MATLAB script (Ganeff, 2012)
BOOST_AUTO_TEST_CASE( testCannonballRadiationPressureTargetInterfaceGaneff )
{
    const auto expectedForce = Eigen::Vector3d( -0.975383093968723e-6, -0.975383093968723e-6, 0.0 );

    CannonballRadiationPressureTargetInterface targetInterface(0.5, 1.21);

    const auto sourceIrradiance = 683.5;  // at (1 AU, 1 AU, 0) from Sun
    const auto sourceToTargetDirection = Eigen::Vector3d(-1, -1, 0).normalized();
    const auto actualForce = targetInterface.evaluateRadiationPressureForce(
            sourceIrradiance,
            sourceToTargetDirection);

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(actualForce, expectedForce, 1e-3);
}

//! Check force when radiation pressure coefficient is zero
BOOST_AUTO_TEST_CASE( testCannonballRadiationPressureTargetInterfaceZeroCoeff )
{
    const auto expectedForce = Eigen::Vector3d::Zero();

    CannonballRadiationPressureTargetInterface targetInterface(1, 0);

    const auto sourceIrradiance = 1367;
    const auto sourceToTargetDirection = Eigen::Vector3d(-1, -1, 0).normalized();
    const auto actualForce = targetInterface.evaluateRadiationPressureForce(
            sourceIrradiance,
            sourceToTargetDirection);

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(actualForce, expectedForce, 1e-3);
}

BOOST_AUTO_TEST_SUITE_END()

} // namespace unit_tests
} // namespace tudat
