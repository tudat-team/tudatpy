#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <iostream>

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/electromagnetism/radiationSourceInterface.h"


namespace tudat
{
namespace unit_test
{

using namespace tudat::electromagnetism;

BOOST_AUTO_TEST_SUITE(test_radiation_source_interface)

//! Test if solar constant and source position is returned when using solar constant-based Sun radiation model
BOOST_AUTO_TEST_CASE( testIsotropicPointRadiationSourceInterface )
{
    auto expectedIrradiance = 1360.8;

    auto sourcePosition = Eigen::Vector3d(0, 0, 0);
    auto targetPosition = Eigen::Vector3d(physical_constants::ASTRONOMICAL_UNIT, 0, 0);

    auto radiantPowerModel = std::make_shared<IrradianceBasedRadiantPowerModel>(
            [=]() { return expectedIrradiance; }, physical_constants::ASTRONOMICAL_UNIT);
    auto radiationSourceInterface = std::make_shared<IsotropicPointRadiationSourceInterface>(
            [=]() { return sourcePosition; }, radiantPowerModel);

    auto ret = radiationSourceInterface->evaluateIrradianceAtPosition(targetPosition).front();
    auto actualIrradiance = std::get<0>(ret);
    auto actualSourcePosition = std::get<1>(ret);

    BOOST_CHECK_CLOSE_FRACTION(actualIrradiance, expectedIrradiance, 1.0e-15);
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(sourcePosition, actualSourcePosition, 1.0e-15);
}

BOOST_AUTO_TEST_SUITE_END()

}
}
