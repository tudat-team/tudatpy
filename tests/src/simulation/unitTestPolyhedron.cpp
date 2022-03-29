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
#include "tudat/astro/gravitation/polyhedronGravityField.h"

namespace tudat
{
namespace unit_tests
{

using namespace simulation_setup;
using namespace basic_astrodynamics;
using namespace input_output;
using namespace reference_frames;
using namespace gravitation;

BOOST_AUTO_TEST_SUITE( test_polyhedron )

BOOST_AUTO_TEST_CASE( test_polyhedron_set_up )
{
    const double gravitationalConstant = 6.67259e-11;
    const double density = 2670;
    const double volume = 2000;
    const double gravitationalParameter = gravitationalConstant * density * volume;
    const std::string associatedReferenceFrame = "Frame";

    Eigen::MatrixXd verticesCoordinates(8,3);
    Eigen::MatrixXi verticesDefiningEachFacet(12,3);

    verticesCoordinates <<
        0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00,
        20.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00,
        0.000000000000000000e+00, 10.000000000000000000e+00, 0.000000000000000000e+00,
        20.000000000000000000e+00, 10.000000000000000000e+00, 0.000000000000000000e+00,
        0.000000000000000000e+00, 0.000000000000000000e+00, 10.000000000000000000e+00,
        20.000000000000000000e+00, 0.000000000000000000e+00, 10.000000000000000000e+00,
        0.000000000000000000e+00, 10.000000000000000000e+00, 10.000000000000000000e+00,
        20.000000000000000000e+00, 10.000000000000000000e+00, 10.000000000000000000e+00;

    verticesDefiningEachFacet <<
        2, 1, 0,
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

//    PolyhedronGravityFieldSettings gravitySettings = simulation_setup::PolyhedronGravityFieldSettings(
//            gravitationalConstant, density, verticesCoordinates, verticesDefiningEachFacet, associatedReferenceFrame);
    PolyhedronGravityFieldSettings gravitySettings = simulation_setup::PolyhedronGravityFieldSettings(
            gravitationalParameter, verticesCoordinates, verticesDefiningEachFacet, associatedReferenceFrame);

    Eigen::MatrixXi verticesDefiningEachEdge = gravitySettings.getVerticesDefiningEachEdge();
    std::vector< Eigen::Vector3d > facetNormalVectors = gravitySettings.getFacetNormalVectors( );
    std::vector< Eigen::MatrixXd > facetDyads = gravitySettings.getFacetDyads( );
    std::vector< Eigen::MatrixXd > edgeDyads = gravitySettings.getEdgeDyads( );

    Eigen::MatrixXd verticesCoordinatesRelativeToFieldPoint;
    Eigen::Vector3d bodyFixedPosition(0.0,3.0,2.0);

    for (unsigned int positionId : {0,1,2})
    {
        if ( positionId == 0 )
        {
            (bodyFixedPosition << 0.0, 0.0, 0.0).finished();
        }
        else if ( positionId == 1 )
        {
            (bodyFixedPosition << 5.0, 0.0, 0.0).finished();
        }
        else
        {
            (bodyFixedPosition << 0.0, 3.0, 2.0).finished();
        }
        std::cout << std::endl << "###### Point " << positionId << std::endl;

        gravitation::calculatePolyhedronVerticesCoordinatesRelativeToFieldPoint(
                verticesCoordinatesRelativeToFieldPoint, bodyFixedPosition, verticesCoordinates);

        Eigen::VectorXd perFacetFactor;
        gravitation::calculatePolyhedronPerFacetFactor(
                perFacetFactor, verticesCoordinatesRelativeToFieldPoint, verticesDefiningEachFacet);

        Eigen::VectorXd perEdgeFactor;
        gravitation::calculatePolyhedronPerEdgeFactor(
                perEdgeFactor, verticesCoordinatesRelativeToFieldPoint, verticesDefiningEachEdge);

        double laplacian = gravitation::calculatePolyhedronLaplacianOfGravitationalPotential(
                gravitySettings.getGravitationalConstantTimesDensity(), perFacetFactor);
        std::cout << "Laplace equation: " << - laplacian / gravitationalConstant / density / M_PI << " pi" << std::endl;

        double potential = gravitation::calculatePolyhedronGravitationalPotential(
                gravitySettings.getGravitationalConstantTimesDensity(), verticesCoordinatesRelativeToFieldPoint, verticesDefiningEachFacet,
                verticesDefiningEachEdge, facetDyads, edgeDyads, perFacetFactor, perEdgeFactor);
        std::cout << "Potential: " << potential << std::endl;

        Eigen::Vector3d acceleration = gravitation::calculatePolyhedronGradientOfGravitationalPotential(
                gravitySettings.getGravitationalConstantTimesDensity(), verticesCoordinatesRelativeToFieldPoint, verticesDefiningEachFacet,
                verticesDefiningEachEdge, facetDyads, edgeDyads, perFacetFactor, perEdgeFactor);
        std::cout << "Acceleration: " << acceleration.transpose() << std::endl;
    }

}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
