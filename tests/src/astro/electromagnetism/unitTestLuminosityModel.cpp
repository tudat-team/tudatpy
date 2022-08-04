#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <iostream>

#include <boost/test/unit_test.hpp>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/electromagnetism/luminosityModel.h"


namespace tudat
{
namespace unit_tests
{

using namespace tudat::electromagnetism;

BOOST_AUTO_TEST_SUITE(test_radiant_power_model)

BOOST_AUTO_TEST_CASE( testConstantLuminosityModel )
{
    const auto expectedLuminosity = 42;

    ConstantLuminosityModel model(expectedLuminosity);
    const auto actualLuminosity = model.getLuminosity();

    BOOST_CHECK_EQUAL(actualLuminosity, expectedLuminosity);
}

BOOST_AUTO_TEST_CASE( testIrradianceBasedLuminosityModel )
{
    const auto expectedLuminosity = 3.8e26;

    IrradianceBasedLuminosityModel model([]() { return 1360.8; },
                                           physical_constants::ASTRONOMICAL_UNIT);
    const auto actualLuminosity = model.getLuminosity();

    BOOST_CHECK_CLOSE(actualLuminosity, expectedLuminosity, 0.1e26);
}

BOOST_AUTO_TEST_SUITE_END()

} // namespace unit_tests
} // namespace tudat
