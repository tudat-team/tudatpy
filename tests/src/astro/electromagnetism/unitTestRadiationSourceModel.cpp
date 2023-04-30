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

#include <limits>
#include <memory>
#include <numeric>

#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "tudat/basics/testMacros.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/math/basic/coordinateConversions.h"
#include "tudat/astro/electromagnetism/radiationSourceModel.h"
#include "tudat/astro/electromagnetism/sourcePanelRadiosityModel.h"
#include "tudat/astro/electromagnetism/surfacePropertyDistribution.h"
#include "tudat/astro/basic_astro/sphericalBodyShapeModel.h"


namespace tudat
{
namespace unit_tests
{

using mathematical_constants::PI;

using namespace tudat::electromagnetism;

BOOST_AUTO_TEST_SUITE(test_radiation_source_model)

//*********************************************************************************************
//   Isotropic point radiation source
//*********************************************************************************************

//! Test if solar constant and source position is returned when using solar constant-based Sun radiation model
BOOST_AUTO_TEST_CASE( testIsotropicPointRadiationSourceModel )
{
    const auto expectedIrradiance = 1360.8;

    auto targetPosition = Eigen::Vector3d(physical_constants::ASTRONOMICAL_UNIT, 0, 0);

    auto luminosityModel = std::make_shared<IrradianceBasedLuminosityModel>(
            [=](double) { return expectedIrradiance; }, physical_constants::ASTRONOMICAL_UNIT);
    auto radiationSourceModel = std::make_shared<IsotropicPointRadiationSourceModel>(luminosityModel);
    radiationSourceModel->updateMembers(TUDAT_NAN);

    const auto actualIrradiance = radiationSourceModel->evaluateIrradianceAtPosition(targetPosition);

    BOOST_CHECK_CLOSE_FRACTION(actualIrradiance, expectedIrradiance, 1.0e-15);
}

//! Test if isotropic point source is invariant w.r.t. position at given distance
BOOST_AUTO_TEST_CASE( testIsotropicPointRadiationSourceModelPositionInvariance )
{
    using mathematical_constants::PI;

    std::vector<double> actualIrradiances;

    const auto radius = 3;

    // Iterate over arbitrary values for zenith and azimuth angles
    for (double zenithAngle : {0., 0.743, 1.903, PI})
    {
        for (double azimuthAngle : {0., 0.984, 2.579, 2*PI})
        {
            auto targetPosition = coordinate_conversions::convertSphericalToCartesian(
                    Eigen::Vector3d(radius, zenithAngle, azimuthAngle));
            auto luminosityModel = std::make_shared<ConstantLuminosityModel>(1);
            auto radiationSourceModel = std::make_shared<IsotropicPointRadiationSourceModel>(luminosityModel);
            radiationSourceModel->updateMembers(TUDAT_NAN);

            const auto actualIrradiance = radiationSourceModel->evaluateIrradianceAtPosition(targetPosition);
            actualIrradiances.push_back(actualIrradiance);
        }
    }

    // Check that all calculated irradiances are identical
    for (auto& e : actualIrradiances)
    {
        BOOST_CHECK_CLOSE_FRACTION(e, actualIrradiances.front(), 1e-15);
    }
}

//*********************************************************************************************
//   Paneled radiation source
//*********************************************************************************************

//! Sanity-check generation of evenly spaced points on sphere
BOOST_AUTO_TEST_CASE( testStaticallyPaneledRadiationSourceModel_Generation )
{
    const auto n = 10;
    const auto radius = 1000;

    const std::vector<std::unique_ptr<SourcePanelRadiosityModel>> radiosityModels{};

    StaticallyPaneledRadiationSourceModel radiationSourceModel(
            "",
            std::make_shared<basic_astrodynamics::SphericalBodyShapeModel>(radius),
            radiosityModels,
            n);
    radiationSourceModel.updateMembers(TUDAT_NAN);
    const auto& panels = radiationSourceModel.getPanels();

    BOOST_CHECK_EQUAL(panels.size(), n);

    // Check radial position and normal of all panels
    for (auto& panel : panels)
    {
        BOOST_CHECK_CLOSE(panel.getRelativeCenter().norm(), radius, 1.0e-10);

        Eigen::Vector3d expectedNormal = panel.getRelativeCenter().normalized();
        TUDAT_CHECK_MATRIX_CLOSE(panel.getSurfaceNormal(), expectedNormal, 1.0e-10);
    }
}

//! Test basic cases for paneled source with albedo panel radiosity model
BOOST_AUTO_TEST_CASE( testStaticallyPaneledRadiationSourceModel_Albedo )
{
    {
        // Target in front of panel, 1 away orthogonal to panel, albedo 1
        const auto expectedIrradiance = 1 / mathematical_constants::PI;
        const auto expectedSourcePosition = Eigen::Vector3d(1, 0, 0);

        std::vector<std::unique_ptr<SourcePanelRadiosityModel>> fullAlbedoRadiosityModel;
        fullAlbedoRadiosityModel.push_back(std::make_unique<AlbedoSourcePanelRadiosityModel>(
                std::make_shared<ConstantSurfacePropertyDistribution>(1)));

        std::vector<PaneledRadiationSourceModel::Panel> panels;
        panels.emplace_back(1,
                            expectedSourcePosition,
                            Eigen::Vector3d(1, 0, 0).normalized(),
                            std::move(fullAlbedoRadiosityModel));

        StaticallyPaneledRadiationSourceModel radiationSourceModel("", std::move(panels));
        radiationSourceModel.updateMembers(TUDAT_NAN);
        const auto irradianceList = radiationSourceModel.evaluateIrradianceAtPosition(
                Eigen::Vector3d(2, 0, 0),
                1,
                -Eigen::Vector3d::UnitX());

        BOOST_CHECK_EQUAL(irradianceList.size(), 1);

        const auto actualIrradiance = std::get<0>(irradianceList.front());
        const auto actualSourcePosition = std::get<1>(irradianceList.front());

        BOOST_CHECK_CLOSE_FRACTION(actualIrradiance, expectedIrradiance, 1e-15);
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(actualSourcePosition, expectedSourcePosition, 1e-15);
    }

    {
        // Target in front of panel, 1 away at 45° angle, area 1, albedo 1
        const auto expectedIrradiance = 1 / (mathematical_constants::PI * sqrt(2));
        const auto expectedSourcePosition = Eigen::Vector3d(-1, 0, 0);

        std::vector<std::unique_ptr<SourcePanelRadiosityModel>> fullAlbedoRadiosityModel;
        fullAlbedoRadiosityModel.push_back(std::make_unique<AlbedoSourcePanelRadiosityModel>(
                std::make_shared<ConstantSurfacePropertyDistribution>(1)));

        std::vector<PaneledRadiationSourceModel::Panel> panels;
        panels.emplace_back(1,
                            expectedSourcePosition,
                            Eigen::Vector3d(1, 0, 0).normalized(),
                            std::move(fullAlbedoRadiosityModel));

        StaticallyPaneledRadiationSourceModel radiationSourceModel("", std::move(panels));
        radiationSourceModel.updateMembers(TUDAT_NAN);
        const auto irradianceList = radiationSourceModel.evaluateIrradianceAtPosition(
                expectedSourcePosition + Eigen::Vector3d(1, 1, 0).normalized(),
                1,
                -Eigen::Vector3d::UnitX());

        BOOST_CHECK_EQUAL(irradianceList.size(), 1);

        const auto actualIrradiance = std::get<0>(irradianceList.front());
        const auto actualSourcePosition = std::get<1>(irradianceList.front());

        BOOST_CHECK_CLOSE_FRACTION(actualIrradiance, expectedIrradiance, 1e-15);
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(actualSourcePosition, expectedSourcePosition, 1e-15);
    }

    {
        // Panels arranged symmetrically around z axis, angled 45° angle upwards, area 1, albedo 0.5
        // Target is above at same distance from all, use panel 1 for expected values
        const auto targetPositionRelativeToPanel = Eigen::Vector3d(-1, 0, 3);
        const auto cosPanelNormalAndTargetPosition =
                targetPositionRelativeToPanel.normalized().dot(Eigen::Vector3d(1, 0, 1).normalized());
        const auto cosPanelNormalAndSourcePosition = 1 / sqrt(2);
        const auto expectedIrradiance = 0.5 * cosPanelNormalAndSourcePosition * cosPanelNormalAndTargetPosition /
                (mathematical_constants::PI * targetPositionRelativeToPanel.squaredNorm());

        const auto expectedSourcePositionPanel1 = Eigen::Vector3d(1, 0, -1);
        const auto expectedSourcePositionPanel2 = Eigen::Vector3d(0, 1, -1);
        const auto expectedSourcePositionPanel3 = Eigen::Vector3d(-1, 0, -1);
        const auto expectedSourcePositionPanel4 = Eigen::Vector3d(0, -1, -1);

        std::vector<std::unique_ptr<SourcePanelRadiosityModel>> halfAlbedoRadiosityModel1;
        halfAlbedoRadiosityModel1.push_back(std::make_unique<AlbedoSourcePanelRadiosityModel>(
                std::make_shared<ConstantSurfacePropertyDistribution>(0.5)));
        std::vector<std::unique_ptr<SourcePanelRadiosityModel>> halfAlbedoRadiosityModel2;
        halfAlbedoRadiosityModel2.push_back(std::make_unique<AlbedoSourcePanelRadiosityModel>(
                std::make_shared<ConstantSurfacePropertyDistribution>(0.5)));
        std::vector<std::unique_ptr<SourcePanelRadiosityModel>> halfAlbedoRadiosityModel3;
        halfAlbedoRadiosityModel3.push_back(std::make_unique<AlbedoSourcePanelRadiosityModel>(
                std::make_shared<ConstantSurfacePropertyDistribution>(0.5)));
        std::vector<std::unique_ptr<SourcePanelRadiosityModel>> halfAlbedoRadiosityModel4;
        halfAlbedoRadiosityModel4.push_back(std::make_unique<AlbedoSourcePanelRadiosityModel>(
                std::make_shared<ConstantSurfacePropertyDistribution>(0.5)));

        std::vector<PaneledRadiationSourceModel::Panel> panels;
        panels.emplace_back(1,
                            expectedSourcePositionPanel1,
                            Eigen::Vector3d(1, 0, 1).normalized(),
                            std::move(halfAlbedoRadiosityModel1));
        panels.emplace_back(1,
                            expectedSourcePositionPanel2,
                            Eigen::Vector3d(0, 1, 1).normalized(),
                            std::move(halfAlbedoRadiosityModel2));
        panels.emplace_back(1,
                            expectedSourcePositionPanel3,
                            Eigen::Vector3d(-1, 0, 1).normalized(),
                            std::move(halfAlbedoRadiosityModel3));
        panels.emplace_back(1,
                            expectedSourcePositionPanel4,
                            Eigen::Vector3d(0, -1, 1).normalized(),
                            std::move(halfAlbedoRadiosityModel4));

        StaticallyPaneledRadiationSourceModel radiationSourceModel("", std::move(panels));
        radiationSourceModel.updateMembers(TUDAT_NAN);
        const auto irradianceList = radiationSourceModel.evaluateIrradianceAtPosition(
                Eigen::Vector3d(0, 0, 2),
                1,
                -Eigen::Vector3d::UnitZ());

        BOOST_CHECK_EQUAL(irradianceList.size(), 4);

        const auto actualIrradiancePanel1 = std::get<0>(irradianceList[0]);
        const auto actualSourcePositionPanel1 = std::get<1>(irradianceList[0]);
        const auto actualIrradiancePanel2 = std::get<0>(irradianceList[1]);
        const auto actualSourcePositionPanel2 = std::get<1>(irradianceList[1]);
        const auto actualIrradiancePanel3 = std::get<0>(irradianceList[2]);
        const auto actualSourcePositionPanel3 = std::get<1>(irradianceList[2]);
        const auto actualIrradiancePanel4 = std::get<0>(irradianceList[3]);
        const auto actualSourcePositionPanel4 = std::get<1>(irradianceList[3]);

        BOOST_CHECK_CLOSE_FRACTION(actualIrradiancePanel1, expectedIrradiance, 1e-15);
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(actualSourcePositionPanel1, expectedSourcePositionPanel1, 1e-15);
        BOOST_CHECK_CLOSE_FRACTION(actualIrradiancePanel2, expectedIrradiance, 1e-15);
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(actualSourcePositionPanel2, expectedSourcePositionPanel2, 1e-15);
        BOOST_CHECK_CLOSE_FRACTION(actualIrradiancePanel3, expectedIrradiance, 1e-15);
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(actualSourcePositionPanel3, expectedSourcePositionPanel3, 1e-15);
        BOOST_CHECK_CLOSE_FRACTION(actualIrradiancePanel4, expectedIrradiance, 1e-15);
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(actualSourcePositionPanel4, expectedSourcePositionPanel4, 1e-15);
    }
}

//! Test polar/azimuth angle to latitude/longitude conversion in constructor
BOOST_AUTO_TEST_CASE( testPaneledRadiationSourceModelPanel )
{
    {
        // North pole
        const auto expectedLatitude = PI / 2;

        std::vector<std::unique_ptr<SourcePanelRadiosityModel>> emptyRadiosityModels {};
        const auto panel = PaneledRadiationSourceModel::Panel(
                1,
                Eigen::Vector3d(0, 0, 100),
                Eigen::Vector3d(0, 0, 1),
                std::move(emptyRadiosityModels)
        );

        const auto actualLatitude = panel.getLatitude();
        BOOST_CHECK_CLOSE_FRACTION(actualLatitude, expectedLatitude, 1e-15);
    }

    {
        // South pole
        const auto expectedLatitude = -PI / 2;

        std::vector<std::unique_ptr<SourcePanelRadiosityModel>> emptyRadiosityModels {};
        const auto panel = PaneledRadiationSourceModel::Panel(
                1,
                Eigen::Vector3d(0, 0, -100),
                Eigen::Vector3d(00, 0, -1),
                std::move(emptyRadiosityModels)
        );

        const auto actualLatitude = panel.getLatitude();
        BOOST_CHECK_CLOSE_FRACTION(actualLatitude, expectedLatitude, 1e-15);
    }

    {
        // Equator at 45° W
        const auto expectedLatitude = 0;
        const auto expectedLongitude = -PI / 4;

        std::vector<std::unique_ptr<SourcePanelRadiosityModel>> emptyRadiosityModels {};
        const auto panel = PaneledRadiationSourceModel::Panel(
                1,
                Eigen::Vector3d(100, -100, 0),
                Eigen::Vector3d(1, -1, 0).normalized(),
                std::move(emptyRadiosityModels)
        );

        const auto actualLatitude = panel.getLatitude();
        const auto actualLongitude = panel.getLongitude();

        BOOST_CHECK_CLOSE_FRACTION(actualLatitude, expectedLatitude, 1e-15);
        BOOST_CHECK_CLOSE_FRACTION(actualLongitude, expectedLongitude, 1e-15);
    }
}

//! Test generation of evenly spaced points on sphere with Saff's algorithm by comparison with Python implementation
// https://github.com/DominikStiller/tudelft-hpb-project/blob/7d27bd188fba6033e5a52c1e95c710d45ea6c09b/analysis/paneling.ipynb
BOOST_AUTO_TEST_CASE( testGenerateEvenlySpacedPoints_Spiraling_Values )
{
    auto n = 10;
    const auto pairOfAngleVectors = generateEvenlySpacedPoints_Spiraling(n);
    const auto polarAngles = std::get<0>(pairOfAngleVectors);
    const auto azimuthAngles = std::get<1>(pairOfAngleVectors);

    BOOST_CHECK_EQUAL(polarAngles.size(), n);
    BOOST_CHECK_EQUAL(azimuthAngles.size(), n);

    // Generated with Python script with visually verified results
    const std::vector<double> actualPolarAngles{
        3.141592653589793,
        2.4619188346815495,
        2.1598272970111707,
        1.9106332362490186,
        1.6821373411358607,
        1.4594553124539327,
        1.2309594173407747,
        0.9817653565786227,
        0.6796738189082441,
        0.0
    };
    const std::vector<double> actualAzimuthAngles{
        0.0,
        1.8112150617748297,
        3.180364954435027,
        4.387841662284913,
        5.53335464780712,
        0.39568232614974086,
        1.6031590339996273,
        2.9723089266598244,
        4.783523988434654,
        0.0
    };

    for (int i = 0; i < n; ++i)
    {
        BOOST_CHECK_CLOSE(polarAngles[i], actualPolarAngles[i], 1e-15);
        BOOST_CHECK_CLOSE(azimuthAngles[i], actualAzimuthAngles[i], 1e-15);
    }
}

//! Test generation of evenly spaced points with Saff's algorithm on sphere by checking bounds
BOOST_AUTO_TEST_CASE( testGenerateEvenlySpacedPoints_Spiraling_Validity )
{
    auto n = 100;
    const auto pairOfAngleVectors = generateEvenlySpacedPoints_Spiraling(n);
    const auto polarAngles = std::get<0>(pairOfAngleVectors);
    const auto azimuthAngles = std::get<1>(pairOfAngleVectors);

    BOOST_CHECK_EQUAL(polarAngles.size(), n);
    BOOST_CHECK_EQUAL(azimuthAngles.size(), n);

    // Check if polar and azimuth angles are valid spherical coordinates
    for (int i = 0; i < n; ++i)
    {
        BOOST_CHECK_GE(polarAngles[i], 0);
        BOOST_CHECK_LE(polarAngles[i], mathematical_constants::PI);

        BOOST_CHECK_GE(azimuthAngles[i], 0);
        BOOST_CHECK_LE(azimuthAngles[i], 2*mathematical_constants::PI);
    }

    // Check if first and last point are poles
    BOOST_CHECK_CLOSE(polarAngles.front(), mathematical_constants::PI, 1e-15);
    BOOST_CHECK_CLOSE(polarAngles.back(), 0, 1e-15);
}

//! Test generation of evenly spaced points on sphere with Wetterer's algorithm by comparison with Python implementation
// https://github.com/DominikStiller/tudelft-hpb-project/blob/7d27bd188fba6033e5a52c1e95c710d45ea6c09b/analysis/paneling.ipynb
BOOST_AUTO_TEST_CASE( testGenerateEvenlySpacedPoints_Staggered_Values )
{
    auto n = 10;
    const auto pairOfAngleVectors = generateEvenlySpacedPoints_Staggered(n);
    const auto actualPolarAngles = std::get<0>(pairOfAngleVectors);
    const auto actualAzimuthAngles = std::get<1>(pairOfAngleVectors);

    BOOST_CHECK_EQUAL(actualPolarAngles.size(), n);
    BOOST_CHECK_EQUAL(actualAzimuthAngles.size(), n);

    // Generated with Python script with visually verified results
    const std::vector<double> expectedPolarAngles{
        2.6905658417935308,
        2.34619382340565,
        2.0943951023931953,
        1.8754889808102941,
        1.6709637479564563,
        1.4706289056333366,
        1.266103672779499,
        1.0471975511965974,
        0.7953988301841433,
        0.4510268117962619
    };
    const std::vector<double> expectedAzimuthAngles{
        0.0,
        2.399963229728653,
        4.799926459457306,
        0.9167043820063725,
        3.3166676117350256,
        5.716630841463679,
        1.833408764012745,
        4.2333719937413985,
        0.35014991629046577,
        2.750113146019119
    };

    for (int i = 0; i < n; ++i)
    {
        BOOST_CHECK_CLOSE(actualPolarAngles[i], expectedPolarAngles[i], 1e-15);
        BOOST_CHECK_CLOSE(actualAzimuthAngles[i], expectedAzimuthAngles[i], 1e-15);
    }
}

//! Test generation of evenly spaced points with Wetterer's algorithm on sphere by checking bounds
BOOST_AUTO_TEST_CASE( testGenerateEvenlySpacedPoints_Staggered_Validity )
{
    auto n = 100;
    const auto pairOfAngleVectors = generateEvenlySpacedPoints_Staggered(n);
    const auto polarAngles = std::get<0>(pairOfAngleVectors);
    const auto azimuthAngles = std::get<1>(pairOfAngleVectors);

    BOOST_CHECK_EQUAL(polarAngles.size(), n);
    BOOST_CHECK_EQUAL(azimuthAngles.size(), n);

    // Check if polar and azimuth angles are valid spherical coordinates
    for (int i = 0; i < n; ++i)
    {
        BOOST_CHECK_GE(polarAngles[i], 0);
        BOOST_CHECK_LE(polarAngles[i], mathematical_constants::PI);

        BOOST_CHECK_GE(azimuthAngles[i], 0);
        BOOST_CHECK_LE(azimuthAngles[i], 2*mathematical_constants::PI);
    }
}

//! Test generation of spherical cap panels using Knocke's algorithm with target at infinity
BOOST_AUTO_TEST_CASE( generatePaneledSphericalCap_FarAway )
{
    const auto expectedNumberOfPanels = 4;

    // Target above North Pole at infinity can see whole Northern Hemisphere
    const Eigen::Vector3d targetPosition(0, 0, std::numeric_limits<double>::max());

    const auto panels = generatePaneledSphericalCap(targetPosition, {3}, 1);
    const auto panelCenters = std::get<0>(panels);
    const auto polarAngles = std::get<1>(panels);
    const auto azimuthAngles = std::get<2>(panels);
    const auto areas = std::get<3>(panels);

    BOOST_CHECK_EQUAL(panelCenters.size(), expectedNumberOfPanels);
    BOOST_CHECK_EQUAL(polarAngles.size(), expectedNumberOfPanels);
    BOOST_CHECK_EQUAL(azimuthAngles.size(), expectedNumberOfPanels);
    BOOST_CHECK_EQUAL(areas.size(), expectedNumberOfPanels);

    // All areas should sum to 2π (one hemisphere)
    const auto actualTotalArea = std::accumulate(areas.begin(), areas.end(), 0.0);
    const auto expectedTotalArea = 2 * PI;
    BOOST_CHECK_CLOSE(actualTotalArea, expectedTotalArea, 1e-15);

    // Check central cap
    BOOST_CHECK_CLOSE(polarAngles[0], 0, 1e-15);
    BOOST_CHECK_CLOSE(azimuthAngles[0], 0, 1e-15);

    // Check ring panels
    // Their center should be at 22.5° latitude (ring stretches from equator to 45°)
    const auto expectedPolarAngle = 3. / 8 * PI;
    BOOST_CHECK_CLOSE(polarAngles[1], expectedPolarAngle, 1e-15);
    BOOST_CHECK_CLOSE(polarAngles[2], expectedPolarAngle, 1e-15);
    BOOST_CHECK_CLOSE(polarAngles[3], expectedPolarAngle, 1e-15);
    // They should be evenly spaced
    const auto expectedAzimuthAngle = 2. / 3 * PI;
    BOOST_CHECK_CLOSE(azimuthAngles[3] - azimuthAngles[2], expectedAzimuthAngle, 1e-13);
    BOOST_CHECK_CLOSE(azimuthAngles[2] - azimuthAngles[1], expectedAzimuthAngle, 1e-13);
    // They should have equal area
    BOOST_CHECK_CLOSE(areas[1], areas[2], 1e-15);
    BOOST_CHECK_CLOSE(areas[2], areas[3], 1e-15);
}

//! Test area of spherical cap panels using Knocke's algorithm
BOOST_AUTO_TEST_CASE( generatePaneledSphericalCap_Area )
{
    const auto radius = 42;

    {
        // Case 1: closer and 2 rings with larger panels
        const Eigen::Vector3d targetPosition(53, 89, 60);
        const auto panels = generatePaneledSphericalCap(targetPosition, {7, 15}, radius);
        const auto areas = std::get<3>(panels);

        // Sum of panel areas should equal spherical cap area
        const auto expectedTotalArea = 2 * PI * radius * radius * (1 - radius / targetPosition.norm());
        const auto actualTotalArea = std::accumulate(areas.begin(), areas.end(), 0.0);
        BOOST_CHECK_CLOSE(actualTotalArea, expectedTotalArea, 1e-13);

        // Central panel area should equal central spherical cap area
        // Central cap angle is 1/3 of total spherical cap angle (there are 2 rings)
        const auto expectedCentralCapArea = 2 * PI * radius * radius * (1 - radius / targetPosition.norm() / 3.);
        const auto actualCentralCapArea = areas[0];
        BOOST_CHECK_CLOSE(actualTotalArea, expectedTotalArea, 1e-13);
    }

    {
        // Case 2: further and 3 rings with smaller panels
        const Eigen::Vector3d targetPosition(100, -129, 98);
        const auto panels = generatePaneledSphericalCap(targetPosition, {30, 60, 80}, radius);
        const auto areas = std::get<3>(panels);

        // Sum of panel areas should equal spherical cap area
        const auto expectedTotalArea = 2 * PI * radius * radius * (1 - radius / targetPosition.norm());
        const auto actualTotalArea = std::accumulate(areas.begin(), areas.end(), 0.0);
        BOOST_CHECK_CLOSE(actualTotalArea, expectedTotalArea, 1e-12);

        // Central panel area should equal central spherical cap area
        // Central cap angle is 1/4 of total spherical cap angle (there are 3 rings)
        const auto expectedCentralCapArea = 2 * PI * radius * radius * (1 - radius / targetPosition.norm() / 4.);
        const auto actualCentralCapArea = areas[0];
        BOOST_CHECK_CLOSE(actualTotalArea, expectedTotalArea, 1e-12);
    }
}

//! Test generation of spherical cap panels using Knocke's algorithm at realistic target position by comparison with
//! Python implementation
// https://github.com/DominikStiller/tudelft-hpb-project/blob/7bf0d9d8f2f0195395edbf86213c2370d78b45d3/analysis/paneling.ipynb
BOOST_AUTO_TEST_CASE( generatePaneledSphericalCap_Realistic )
{
    const auto expectedNumberOfPanels = 18;

    // Lunar radius
    const auto radius = 1736e3;
    // Moon orbiter at 50 km altitude at random position
    const Eigen::Vector3d targetPosition(389737.1519614824, 1558948.6078459297, -779474.3039229648);

    const auto panels = generatePaneledSphericalCap(targetPosition, {6, 11}, radius);
    const auto actualPanelCenters = std::get<0>(panels);
    const auto actualPolarAngles = std::get<1>(panels);
    const auto actualAzimuthAngles = std::get<2>(panels);
    const auto actualAreas = std::get<3>(panels);

    // Generated with Python script with visually verified results
    const std::vector<double> expectedPolarAngles{
            2.0224297678744287,
            2.047877556334641,
            2.136150520476182,
            2.105861502237243,
            1.9905687343731608,
            1.9081869714047708,
            1.9356833964512896,
            2.0609813958629832,
            2.1611485204045766,
            2.2163640745411137,
            2.2037113459078386,
            2.1287065837117263,
            2.0209606333816557,
            1.9161553515699112,
            1.8445432907351633,
            1.825732091599134,
            1.8648422421398858,
            1.951226570528701
    };
    const std::vector<double> expectedAzimuthAngles{
            1.3258176636680326,
            1.1962484139825238,
            1.2872726659199039,
            1.421566020343163,
            1.4518429848343788,
            1.3603079148388213,
            1.2376672636524402,
            1.1081744372204867,
            1.1631836508258644,
            1.2809632745094919,
            1.4180375866156136,
            1.5164475025054232,
            1.5454931543149055,
            1.5068002183873146,
            1.418141127572344,
            1.3051176681093595,
            1.197514780036697,
            1.1241208988411142
    };
    const std::vector<double> expectedAreas{
            59147445289.420876,
            29512138145.00732,
            29512138145.00732,
            29512138145.00732,
            29512138145.00732,
            29512138145.00732,
            29512138145.00732,
            26717454531.532776,
            26717454531.532776,
            26717454531.532776,
            26717454531.532776,
            26717454531.532776,
            26717454531.532776,
            26717454531.532776,
            26717454531.532776,
            26717454531.532776,
            26717454531.532776,
            26717454531.532776
    };

    for (int i = 0; i < expectedNumberOfPanels; ++i)
    {
        BOOST_CHECK_CLOSE(actualPanelCenters[i].norm(), radius, 1e-13);
        BOOST_CHECK_CLOSE(actualPolarAngles[i], expectedPolarAngles[i], 1e-13);
        BOOST_CHECK_CLOSE(actualAzimuthAngles[i], expectedAzimuthAngles[i], 1e-13);
        BOOST_CHECK_CLOSE(actualAreas[i], expectedAreas[i], 1e-13);
    }
}

BOOST_AUTO_TEST_SUITE_END()

}
}
