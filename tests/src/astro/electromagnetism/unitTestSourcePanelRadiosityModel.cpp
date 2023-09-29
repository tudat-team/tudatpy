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
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/astro/electromagnetism/sourcePanelRadiosityModel.h"
#include "tudat/astro/electromagnetism/surfacePropertyDistribution.h"


namespace tudat
{
namespace unit_tests
{

using mathematical_constants::PI;

using namespace tudat::electromagnetism;

BOOST_AUTO_TEST_SUITE(test_source_panel_radiosity_model)

//! Test basic cases for constant panel radiosity model
BOOST_AUTO_TEST_CASE( testConstantPanelRadiosityModel )
{
    ConstantSourcePanelRadiosityModel radiosityModel(42);
    radiosityModel.updateMembers(TUDAT_NAN, TUDAT_NAN, TUDAT_NAN);

    {
        // Target in front of panel, 1 away orthogonal to panel
        const auto expectedEmittedIrradiance = 42 / mathematical_constants::PI;
        const auto actualEmittedIrradiance = radiosityModel.evaluateIrradianceAtPosition(
                1,
                Eigen::Vector3d::UnitX(),
                Eigen::Vector3d::UnitX());
        BOOST_CHECK_CLOSE(expectedEmittedIrradiance, actualEmittedIrradiance, 1e-15);
    }

    {
        // Target at 45° angle with panel, 2 away orthogonal to panel
        const auto expectedEmittedIrradiance = 42 / (mathematical_constants::PI * sqrt(2) * 4);
        const auto actualEmittedIrradiance = radiosityModel.evaluateIrradianceAtPosition(
                1,
                Eigen::Vector3d::UnitX(),
                2 * Eigen::Vector3d(1, 0, 1).normalized());
        BOOST_CHECK_CLOSE(expectedEmittedIrradiance, actualEmittedIrradiance, 1e-10);
    }

    {
        // Target in-plane with panel
        const auto expectedEmittedIrradiance = 0;
        const auto actualEmittedIrradiance = radiosityModel.evaluateIrradianceAtPosition(
                1,
                Eigen::Vector3d::UnitX(),
                3 * Eigen::Vector3d::UnitY());
        BOOST_CHECK_CLOSE(expectedEmittedIrradiance, actualEmittedIrradiance, 1e-15);
    }

    {
        // Target behind panel
        const auto expectedEmittedIrradiance = 0;
        const auto actualEmittedIrradiance = radiosityModel.evaluateIrradianceAtPosition(
                1,
                Eigen::Vector3d::UnitX(),
                -Eigen::Vector3d::UnitX());
        BOOST_CHECK_CLOSE(expectedEmittedIrradiance, actualEmittedIrradiance, 1e-15);
    }
}

//! Test basic cases for albedo panel radiosity model
BOOST_AUTO_TEST_CASE( testAlbedoPanelRadiosityModel )
{
    AlbedoSourcePanelRadiosityModel radiosityModel("", std::make_shared<ConstantSurfacePropertyDistribution>(1));
    radiosityModel.updateMembers(TUDAT_NAN, TUDAT_NAN, TUDAT_NAN);

    {
        // Target in front of panel, 1 away orthogonal to panel
        const auto expectedReflectedIrradiance = 1 / mathematical_constants::PI;
        radiosityModel.updateOriginalSourceProperties(1.0, 1.0, -Eigen::Vector3d::UnitX());
        const auto actualReflectedIrradiance = radiosityModel.evaluateIrradianceAtPosition(
                1,
                Eigen::Vector3d::UnitX(),
                Eigen::Vector3d::UnitX() );
        BOOST_CHECK_CLOSE(actualReflectedIrradiance, expectedReflectedIrradiance, 1e-15);
    }

    {
        // Target in front of panel, 1 away at 45° angle
        const auto expectedReflectedIrradiance = 1 / (mathematical_constants::PI * sqrt(2));
        radiosityModel.updateOriginalSourceProperties(1.0, 1.0, -Eigen::Vector3d::UnitX());
        const auto actualReflectedIrradiance = radiosityModel.evaluateIrradianceAtPosition(
                1,
                Eigen::Vector3d::UnitX(),
                Eigen::Vector3d(1, 1, 0).normalized());
        BOOST_CHECK_CLOSE(actualReflectedIrradiance, expectedReflectedIrradiance, 1e-10);
    }

    {
        // Target behind panel
        const auto expectedReflectedIrradiance = 0;
        radiosityModel.updateOriginalSourceProperties(1.0, 1.0, -Eigen::Vector3d::UnitX());
        const auto actualReflectedIrradiance = radiosityModel.evaluateIrradianceAtPosition(
                1,
                Eigen::Vector3d::UnitX(),
                -Eigen::Vector3d::UnitX());
        BOOST_CHECK_CLOSE(actualReflectedIrradiance, expectedReflectedIrradiance, 1e-15);
    }

    {
        // Original source behind panel
        const auto expectedReflectedIrradiance = 0;
        radiosityModel.updateOriginalSourceProperties(1.0, 1.0, Eigen::Vector3d::UnitX());
        const auto actualReflectedIrradiance = radiosityModel.evaluateIrradianceAtPosition(
                1,
                Eigen::Vector3d::UnitX(),
                Eigen::Vector3d::UnitX());
        BOOST_CHECK_CLOSE(actualReflectedIrradiance, expectedReflectedIrradiance, 1e-15);
    }
}

//! Test basic cases for delayed thermal panel radiosity model
BOOST_AUTO_TEST_CASE( testDelayedThermalPanelRadiosityModel )
{
    DelayedThermalSourcePanelRadiosityModel radiosityModel("", std::make_shared<ConstantSurfacePropertyDistribution>(1));
    radiosityModel.updateMembers(TUDAT_NAN, TUDAT_NAN, TUDAT_NAN);

    {
        // Target in front of panel, 1 away orthogonal to panel
        const auto expectedEmittedIrradiance = 1 / mathematical_constants::PI;
        radiosityModel.updateOriginalSourceProperties(4.0,4.0, -Eigen::Vector3d::UnitX());
        const auto actualEmittedIrradiance = radiosityModel.evaluateIrradianceAtPosition(
                1,
                Eigen::Vector3d::UnitX(),
                Eigen::Vector3d::UnitX());
        BOOST_CHECK_CLOSE(expectedEmittedIrradiance, actualEmittedIrradiance, 1e-10);
    }

    {
        // Target at 45° angle with panel, 2 away orthogonal to panel
        const auto expectedEmittedIrradiance = 1 / (mathematical_constants::PI * sqrt(2) * 4);
        radiosityModel.updateOriginalSourceProperties(4.0,4.0, -Eigen::Vector3d(1, 1, 0).normalized());
        const auto actualEmittedIrradiance = radiosityModel.evaluateIrradianceAtPosition(
                1,
                Eigen::Vector3d::UnitX(),
                2 * Eigen::Vector3d(1, 0, 1).normalized());
        BOOST_CHECK_CLOSE(expectedEmittedIrradiance, actualEmittedIrradiance, 1e-10);
    }

    {
        // Target behind panel
        const auto expectedEmittedIrradiance = 0;
        radiosityModel.updateOriginalSourceProperties(1.0, 1.0, -Eigen::Vector3d::UnitX());
        const auto actualEmittedIrradiance = radiosityModel.evaluateIrradianceAtPosition(
                1,
                Eigen::Vector3d::UnitX(),
                -Eigen::Vector3d::UnitX());
        BOOST_CHECK_CLOSE(expectedEmittedIrradiance, actualEmittedIrradiance, 1e-10);
    }
}

//! Test basic cases for angle-based thermal panel radiosity model
BOOST_AUTO_TEST_CASE( testAngleBasedThermalPanelRadiosityModel )
{
    // Based on Lemoine (2013)
    AngleBasedThermalSourcePanelRadiosityModel radiosityModel("",
                                                              100, 375,
                                                              std::make_shared<ConstantSurfacePropertyDistribution>(1));
    radiosityModel.updateMembers(TUDAT_NAN, TUDAT_NAN, TUDAT_NAN);

    {
        // Target in front of panel, 1 away orthogonal to panel
        const auto expectedEmittedIrradiance = 1121.3386912573242 / mathematical_constants::PI;
        radiosityModel.updateOriginalSourceProperties(1.0, 1.0, -Eigen::Vector3d::UnitX());
        const auto actualEmittedIrradiance = radiosityModel.evaluateIrradianceAtPosition(
                1,
                Eigen::Vector3d::UnitX(),
                Eigen::Vector3d::UnitX());
        BOOST_CHECK_CLOSE(expectedEmittedIrradiance, actualEmittedIrradiance, 1e-2);
    }

    {
        // Original source at 45° angle with panel, target in front of panel, 1 away orthogonal to panel
        const auto expectedEmittedIrradiance = 792.9061925949023 / mathematical_constants::PI;
        radiosityModel.updateOriginalSourceProperties(1.0, 1.0, -Eigen::Vector3d(1, 1, 0).normalized());
        const auto actualEmittedIrradiance = radiosityModel.evaluateIrradianceAtPosition(
                1,
                Eigen::Vector3d::UnitX(),
                Eigen::Vector3d::UnitX());
        BOOST_CHECK_CLOSE(expectedEmittedIrradiance, actualEmittedIrradiance, 1e-2);
    }

    {
        // Original source at 45° angle with panel, target at 45° angle with panel, 2 away orthogonal to panel
        const auto expectedEmittedIrradiance = 792.9061925949023 / (mathematical_constants::PI * sqrt(2) * 4);
        radiosityModel.updateOriginalSourceProperties(1.0, 1.0, -Eigen::Vector3d(1, 1, 0).normalized());
        const auto actualEmittedIrradiance = radiosityModel.evaluateIrradianceAtPosition(
                1,
                Eigen::Vector3d::UnitX(),
                2 * Eigen::Vector3d(1, 0, 1).normalized());
        BOOST_CHECK_CLOSE(expectedEmittedIrradiance, actualEmittedIrradiance, 1e-2);
    }

    {
        // Target behind panel
        const auto expectedEmittedIrradiance = 0;
        radiosityModel.updateOriginalSourceProperties(1.0, 1.0, -Eigen::Vector3d::UnitX());
        const auto actualEmittedIrradiance = radiosityModel.evaluateIrradianceAtPosition(
                1,
                Eigen::Vector3d::UnitX(),
                -Eigen::Vector3d::UnitX());
        BOOST_CHECK_CLOSE(expectedEmittedIrradiance, actualEmittedIrradiance, 1e-2);
    }

    {
        // Original source behind panel
        const auto expectedReflectedIrradiance = 5.670374419 / mathematical_constants::PI;
        radiosityModel.updateOriginalSourceProperties(1.0, 1.0, Eigen::Vector3d::UnitX());
        const auto actualReflectedIrradiance = radiosityModel.evaluateIrradianceAtPosition(
                1,
                Eigen::Vector3d::UnitX(),
                Eigen::Vector3d::UnitX());
        BOOST_CHECK_CLOSE(actualReflectedIrradiance, expectedReflectedIrradiance, 1e-2);
    }
}

BOOST_AUTO_TEST_SUITE_END()

}
}
