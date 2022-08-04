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
 *      Noomen, R. AE2230-I Flight and Orbital Mechanics Lecture Notes, Ch. Perturbations (2),
 *             Delft University of Technology, 2022.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <memory>

#include <boost/test/unit_test.hpp>
#include <Eigen/Core>

#include "tudat/math/basic/rotationRepresentations.h"
#include "tudat/basics/testMacros.h"
#include "tudat/astro/electromagnetism/radiationPressureAcceleration.h"


namespace tudat
{
namespace unit_tests
{

using namespace tudat::electromagnetism;

BOOST_AUTO_TEST_SUITE(test_radiation_pressure_acceleration)

//! Test acceleration with all unity values
BOOST_AUTO_TEST_CASE( testRadiationPressureAccelerationUnity )
{
    const auto expectedAcceleration = Eigen::Vector3d::UnitX();

    auto luminosityModel = std::make_shared<IrradianceBasedLuminosityModel>(
            [] () { return physical_constants::SPEED_OF_LIGHT; }, 1);
    auto sourceModel = std::make_shared<IsotropicPointRadiationSourceModel>(
            [] () { return Eigen::Vector3d::Zero(); }, luminosityModel);
    auto targetModel = std::make_shared<CannonballRadiationPressureTargetModel>(
            1, 1);
    RadiationPressureAcceleration accelerationModel(
            sourceModel,
            targetModel,
            []() { return Eigen::Vector3d::UnitX(); },
            []() { return 1; },
            []() { return Eigen::Quaterniond::Identity(); });

    accelerationModel.updateMembers(0);
    const auto actualAcceleration = accelerationModel.getAcceleration();

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(actualAcceleration, expectedAcceleration, 1e-15);
}

//! Test acceleration for idealized GOCE spacecraft (Noomen 2022)
BOOST_AUTO_TEST_CASE( testRadiationPressureAccelerationGOCE )
{
    const auto expectedAcceleration = Eigen::Vector3d(1, 1, 0).normalized() * 5.2e-9;

    auto luminosityModel = std::make_shared<IrradianceBasedLuminosityModel>(
            [] () { return 1371; }, physical_constants::ASTRONOMICAL_UNIT);
    auto sourceModel = std::make_shared<IsotropicPointRadiationSourceModel>(
            [] () { return Eigen::Vector3d::Zero(); }, luminosityModel);
    auto targetModel = std::make_shared<CannonballRadiationPressureTargetModel>(
            1, 1.2);
    RadiationPressureAcceleration accelerationModel(
            sourceModel,
            targetModel,
            []() { return Eigen::Vector3d(1, 1, 0).normalized() * physical_constants::ASTRONOMICAL_UNIT; },
            []() { return 1050; },
            []() { return Eigen::Quaterniond::Identity(); });

    accelerationModel.updateMembers(0);
    const auto actualAcceleration = accelerationModel.getAcceleration();

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(actualAcceleration, expectedAcceleration, 1e-2);
}

//! Test that cannonball acceleration is invariant under target rotation
BOOST_AUTO_TEST_CASE( testRadiationPressureAccelerationCannonballRotationInvariance )
{
    using mathematical_constants::PI;

    std::vector<Eigen::Vector3d> actualAccelerations;

    // Iterate over arbitrary values for Euler angles
    for (double x : {0., 0.984, 2.579, 2*PI})
    {
        for (double y : {0., 0.743, 1.903, PI})
        {
            for (double z : {0., 0.646, 5.634, 2*PI})
            {
                auto rotation = basic_mathematics::getQuaternionFrom313EulerAngles(Eigen::Vector3d(x, y, z));

                auto luminosityModel = std::make_shared<ConstantLuminosityModel>(1);
                auto sourceModel = std::make_shared<IsotropicPointRadiationSourceModel>(
                        [] () { return Eigen::Vector3d::Zero(); }, luminosityModel);
                auto targetModel = std::make_shared<CannonballRadiationPressureTargetModel>(
                        1, 1);
                RadiationPressureAcceleration accelerationModel(
                        sourceModel,
                        targetModel,
                        []() { return Eigen::Vector3d(1, 1, 0).normalized(); },
                        []() { return 1; },
                        [=]() { return rotation; });

                accelerationModel.updateMembers(0);
                actualAccelerations.push_back(accelerationModel.getAcceleration());
            }
        }
    }

    // Check that all calculated irradiances are accelerations
    BOOST_CHECK(std::all_of(
            actualAccelerations.begin(),
            actualAccelerations.end(),
            [&] (Eigen::Vector3d const &e) { return e.isApprox(actualAccelerations.front()); }));
}

BOOST_AUTO_TEST_SUITE_END()

} // namespace unit_tests
} // namespace tudat
