#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <iostream>

#include <boost/test/unit_test.hpp>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/electromagnetism/radiantPowerModel.h"


namespace tudat
{
namespace unit_tests
{

using namespace tudat::electromagnetism;

BOOST_AUTO_TEST_SUITE(test_radiant_power_model)

BOOST_AUTO_TEST_CASE( testConstantRadiantPowerModel )
{
    const auto expectedRadiantPower = 42;

    ConstantRadiantPowerModel model(expectedRadiantPower);
    const auto actualRadiantPower = model.getRadiantPower();

    BOOST_CHECK_EQUAL(actualRadiantPower, expectedRadiantPower);
}

BOOST_AUTO_TEST_CASE( testIrradianceBasedRadiantPowerModel )
{
    const auto expectedRadiantPower = 3.8e26;

    IrradianceBasedRadiantPowerModel model([]() { return 1360.8; },
                                           physical_constants::ASTRONOMICAL_UNIT);
    const auto actualRadiantPower = model.getRadiantPower();

    BOOST_CHECK_CLOSE(actualRadiantPower, expectedRadiantPower, 0.1e26);
}

BOOST_AUTO_TEST_SUITE_END()

} // namespace unit_tests
} // namespace tudat
