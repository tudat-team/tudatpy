#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <iostream>

#include <boost/test/unit_test.hpp>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/electromagnetism/radiantPowerModel.h"


namespace tudat
{
namespace unit_test
{

using namespace tudat::electromagnetism;

BOOST_AUTO_TEST_SUITE(test_radiant_power_model)

BOOST_AUTO_TEST_CASE( testConstantRadiantPowerModel )
{
    ConstantRadiantPowerModel model(42);

    BOOST_CHECK_EQUAL(model.getRadiantPower(), 42);
}

BOOST_AUTO_TEST_CASE( testIrradianceBasedRadiantPowerModel )
{
    IrradianceBasedRadiantPowerModel model([]() { return 1360.8; },
                                           physical_constants::ASTRONOMICAL_UNIT);

    BOOST_CHECK_CLOSE(model.getRadiantPower(), 3.8e26, 0.1e26);
}

BOOST_AUTO_TEST_SUITE_END()

}
}
