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

    IrradianceBasedLuminosityModel luminosityModel([](double) { return 1360.8; }, physical_constants::ASTRONOMICAL_UNIT);
    luminosityModel.updateMembers(TUDAT_NAN);
    const auto actualLuminosity = luminosityModel.getLuminosity();

    BOOST_CHECK_CLOSE(actualLuminosity, expectedLuminosity, 0.1);
}

BOOST_AUTO_TEST_SUITE_END()

} // namespace unit_tests
} // namespace tudat
