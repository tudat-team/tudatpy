#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <memory>

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/math/basic/coordinateConversions.h"
#include "tudat/astro/electromagnetism/radiationSourceInterface.h"


namespace tudat
{
namespace unit_tests
{

using namespace tudat::electromagnetism;

BOOST_AUTO_TEST_SUITE(test_radiation_source_interface)

//! Test if solar constant and source position is returned when using solar constant-based Sun radiation model
BOOST_AUTO_TEST_CASE( testIsotropicPointRadiationSourceInterface )
{
    const auto expectedIrradiance = 1360.8;

    auto sourcePosition = Eigen::Vector3d(0, 0, 0);
    auto targetPosition = Eigen::Vector3d(physical_constants::ASTRONOMICAL_UNIT, 0, 0);

    auto radiantPowerModel = std::make_shared<IrradianceBasedRadiantPowerModel>(
            [=]() { return expectedIrradiance; }, physical_constants::ASTRONOMICAL_UNIT);
    auto radiationSourceInterface = std::make_shared<IsotropicPointRadiationSourceInterface>(
            [=]() { return sourcePosition; }, radiantPowerModel);

    const auto ret = radiationSourceInterface->evaluateIrradianceAtPosition(targetPosition).front();
    const auto actualIrradiance = std::get<0>(ret);
    const auto actualSourcePosition = std::get<1>(ret);

    BOOST_CHECK_CLOSE_FRACTION(actualIrradiance, expectedIrradiance, 1.0e-15);
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(sourcePosition, actualSourcePosition, 1.0e-15);
}

//! Test if isotropic point source is invariant w.r.t. position at given distance
BOOST_AUTO_TEST_CASE( testIsotropicPointRadiationSourceInterfacePositionInvariance )
{
    using mathematical_constants::PI;

    std::vector<double> actualIrradiances;

    const auto radius = 1;

    // Iterate over arbitrary values for zenith and azimuth angles
    for (double zenithAngle : {0., 0.743, 1.903, PI})
    {
        for (double azimuthAngle : {0., 0.984, 2.579, 2*PI})
        {
            auto targetPosition = coordinate_conversions::convertSphericalToCartesian(
                    Eigen::Vector3d(radius, zenithAngle, azimuthAngle));
            auto radiantPowerModel = std::make_shared<ConstantRadiantPowerModel>(1);
            auto radiationSourceInterface = std::make_shared<IsotropicPointRadiationSourceInterface>(
                    [=]() { return Eigen::Vector3d(0, 0, 0); }, radiantPowerModel);

            const auto ret = radiationSourceInterface->evaluateIrradianceAtPosition(targetPosition).front();
            actualIrradiances.push_back(std::get<0>(ret));
        }
    }

    // Check that all calculated irradiances are identical
    BOOST_CHECK(std::all_of(
            actualIrradiances.begin(),
            actualIrradiances.end(),
            [&] (double const &e) { return e == actualIrradiances.front(); }));
}

BOOST_AUTO_TEST_SUITE_END()

}
}
