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
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <limits>
#include "fstream"
#include "iostream"

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/astro/aerodynamics/customAerodynamicCoefficientInterface.h"
#include "tudat/astro/aerodynamics/aerodynamicAcceleration.h"
#include "tudat/astro/reference_frames/aerodynamicAngleCalculator.h"
#include "tudat/simulation/propagation_setup/dynamicsSimulator.h"
#include "tudat/interface/spice/spiceEphemeris.h"
#include "tudat/interface/spice/spiceRotationalEphemeris.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/estimation_setup/createNumericalSimulator.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"

using namespace tudat::aerodynamics;


namespace tudat {
    namespace unit_tests {
        using namespace tudat;
        using namespace aerodynamics;
        using namespace simulation_setup;
        using namespace numerical_integrators;
        using namespace simulation_setup;
        using namespace basic_astrodynamics;
        using namespace propagators;
        using namespace basic_mathematics;
        using namespace basic_astrodynamics;

        BOOST_AUTO_TEST_SUITE(test_mars_dtm_atmosphere)

            BOOST_AUTO_TEST_CASE(testMarsDtmAtmosphere)
            {
                // Define tolerance for equality
                double tolerance = 1.0E-15;
                std::shared_ptr<AtmosphereSettings> marsDtmAtmosphereSettings;
                marsDtmAtmosphereSettings = std::make_shared<MarsDtmAtmosphereSettings>();
                std::shared_ptr<aerodynamics::AtmosphereModel> marsAtmosphereModel = createAtmosphereModel(
                        marsDtmAtmosphereSettings, "Mars");
                std::shared_ptr<MarsDtmAtmosphereModel> atmosphereModel =
                        std::dynamic_pointer_cast<MarsDtmAtmosphereModel>(
                                createAtmosphereModel(marsDtmAtmosphereSettings, "Mars"));
                int alt_km = 400E3;
                int time = 86400;
                double latitude = unit_conversions::convertDegreesToRadians(15.0);
                double longitude = unit_conversions::convertDegreesToRadians(10.0);
                double rho = atmosphereModel->getDensity(alt_km, longitude, latitude, time);

                // reference density for this specific inputs
                double rho_ref = 2.57862799E-14;
                double rho_diff = std::abs(rho - rho_ref);
                // Check density
                BOOST_CHECK_SMALL(rho - rho_ref, tolerance);
            }

        BOOST_AUTO_TEST_SUITE_END()

    } // namespace unit_tests
}//namespace tudat
//{
//namespace unit_tests
//{
//
//BOOST_AUTO_TEST_SUITE( test_mars_dtm_atmosphere )
//
//BOOST_AUTO_TEST_CASE( testMarsDtmAtmosphere )
//{
//
//}
//
//BOOST_AUTO_TEST_SUITE_END( )
//
//} // namespace unit_tests
//} // namespace tudat
