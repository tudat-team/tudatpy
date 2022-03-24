#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <limits>

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/ephemerides/approximatePlanetPositions.h"
#include "tudat/astro/ephemerides/tabulatedEphemeris.h"
#include "tudat/basics/testMacros.h"

#include "tudat/interface/spice/spiceEphemeris.h"
#include "tudat/interface/spice/spiceRotationalEphemeris.h"

#include "tudat/io/basicInputOutput.h"
#include "tudat/math/interpolators/linearInterpolator.h"
#include "tudat/simulation/estimation_setup/createNumericalSimulator.h"
#include "tudat/simulation/environment_setup/createBodies.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/astro/aerodynamics/testApolloCapsuleCoefficients.h"
#include "tudat/astro/ephemerides/constantRotationalEphemeris.h"
#include "tudat/astro/electromagnetism/radiationPressureInterface.h"

namespace tudat
{
namespace unit_tests
{

using namespace simulation_setup;
using namespace basic_astrodynamics;
using namespace input_output;
using namespace reference_frames;

BOOST_AUTO_TEST_SUITE( test_polyhedron )

BOOST_AUTO_TEST_CASE( test_polyhedron_set_up )
{
    const double gravitationalConstant = 5;
    const double density = 42;
    const std::string associatedReferenceFrame = "Frame";

    Eigen::MatrixXd verticesCoordinates(8,3);
    Eigen::MatrixXi verticesDefiningEachFacet(12,3);

    verticesCoordinates << -5.000000000000000000e-01, -5.000000000000000000e-01, -5.000000000000000000e-01,
    5.000000000000000000e-01, -5.000000000000000000e-01, -5.000000000000000000e-01,
    -5.000000000000000000e-01, 5.000000000000000000e-01, -5.000000000000000000e-01,
    5.000000000000000000e-01, 5.000000000000000000e-01, -5.000000000000000000e-01,
    -5.000000000000000000e-01, -5.000000000000000000e-01, 5.000000000000000000e-01,
    5.000000000000000000e-01, -5.000000000000000000e-01, 5.000000000000000000e-01,
    -5.000000000000000000e-01, 5.000000000000000000e-01, 5.000000000000000000e-01,
    5.000000000000000000e-01, 5.000000000000000000e-01, 5.000000000000000000e-01;

    verticesDefiningEachFacet << 2, 1, 0,
        1, 2, 3,
        4, 2, 0,
        2, 4, 6,
        1, 4, 0,
        4, 1, 5,
        6, 5, 7,
        5, 6, 4,
        3, 6, 7,
        6, 3, 2,
        5, 3, 7,
        3, 5, 1;

    PolyhedronGravityFieldSettings(gravitationalConstant, density, verticesCoordinates, verticesDefiningEachFacet,
                                   associatedReferenceFrame);
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
