#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <iostream>

#include <boost/test/unit_test.hpp>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/celestialBodyConstants.h"
#include "tudat/astro/electromagnetism/luminosityModel.h"


namespace tudat
{
namespace unit_tests
{

using namespace tudat::electromagnetism;

BOOST_AUTO_TEST_SUITE(test_luminosity_model)

BOOST_AUTO_TEST_CASE( testConstantLuminosityModel )
{
    const auto expectedLuminosity = 42;

    ConstantLuminosityModel luminosityModel(expectedLuminosity);
    const auto actualLuminosity = luminosityModel.getLuminosity();

    BOOST_CHECK_EQUAL(actualLuminosity, expectedLuminosity);
}

BOOST_AUTO_TEST_CASE( testIrradianceBasedLuminosityModel )
{
    const auto expectedLuminosity = celestial_body_constants::SUN_LUMINOSITY;

    IrradianceBasedLuminosityModel luminosityModel(1360.8, physical_constants::ASTRONOMICAL_UNIT);
    luminosityModel.updateMembers();
    const auto actualLuminosity = luminosityModel.getLuminosity();

    BOOST_CHECK_CLOSE(actualLuminosity, expectedLuminosity, 0.1);
}

BOOST_AUTO_TEST_SUITE_END()

} // namespace unit_tests
} // namespace tudat
