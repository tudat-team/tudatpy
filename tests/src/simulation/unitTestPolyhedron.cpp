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

using namespace orbital_element_conversions;
using namespace numerical_integrators;
using namespace propagators;
using namespace ephemerides;

BOOST_AUTO_TEST_SUITE( test_polyhedron )

BOOST_AUTO_TEST_CASE( test_polyhedron_propagation )
{
    const std::string associatedReferenceFrame = "IAU_Sun";

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

    // Set simulation time settings.
    const double simulationStartEpoch = 0.0;
    const double simulationEndEpoch = tudat::physical_constants::JULIAN_DAY;
    const double environmentTimeBuffer = 0.0;

    // Create bodies in simulation
    std::vector< std::string > bodiesToCreate = { "Sun", "Earth" };
    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );
    BodyListSettings bodySettings = getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - environmentTimeBuffer,
                                    simulationEndEpoch + environmentTimeBuffer,  "Earth", "J2000" );

    bodySettings.get("Sun")->gravityFieldSettings = //centralGravitySettings(1.327e20);
            polyhedronGravitySettings(
            1.327e20,
            verticesCoordinates, verticesDefiningEachFacet, associatedReferenceFrame);

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // Create spacecraft object.
    bodies.createEmptyBody( "Asterix" );
    bodies.at( "Asterix" )->setConstantBodyMass( 0.0 );

    // Define propagation settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfAsterix;
    accelerationsOfAsterix[ "Earth" ] = { pointMassGravityAcceleration( ) };
    accelerationsOfAsterix[ "Sun" ] = { polyhedronAcceleration( ) }; // { polyhedronAcceleration( ) }; //{ pointMassGravityAcceleration( ) };

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationSettingsList;
    accelerationSettingsList[ "Asterix" ] = accelerationsOfAsterix;

    std::vector< std::string > bodiesToPropagate = { "Asterix" };
    std::vector< std::string > centralBodies = { "Earth" };

    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodies, accelerationSettingsList, bodiesToPropagate, centralBodies );

    // Set Keplerian elements for Asterix.
    Eigen::Vector6d asterixInitialStateInKeplerianElements;
    asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7500.0E3;
    asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.1;
    asterixInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
    asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = unit_conversions::convertDegreesToRadians( 235.7 );
    asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = unit_conversions::convertDegreesToRadians( 23.4 );
    asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );

    double earthGravitationalParameter = getBodyGravitationalParameter( bodies, "Earth" );
    const Eigen::Vector6d asterixInitialState = convertKeplerianToCartesianElements(
                asterixInitialStateInKeplerianElements, earthGravitationalParameter );

    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, asterixInitialState, simulationEndEpoch );

    const double fixedStepSize = 10.0;
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > > ( rungeKutta4, 0.0, fixedStepSize );

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator( bodies, integratorSettings, propagatorSettings );
    std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

    dynamicsSimulator.integrateEquationsOfMotion(asterixInitialState);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        PROVIDE OUTPUT TO CONSOLE AND FILES           //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Eigen::VectorXd finalIntegratedState = (--integrationResult.end( ) )->second;
    // Print the position (in km) and the velocity (in km/s) at t = 0.
    std::cout << "The initial position vector of Asterix is [km]:" << std::endl <<
                 asterixInitialState.segment( 0, 3 ).transpose() / 1E3 << std::endl;
    // Print the position (in km) and the velocity (in km/s) at t = 86400.
    std::cout << "After " << simulationEndEpoch <<
                 " seconds, the position vector of Asterix is [km]:" << std::endl <<
                 finalIntegratedState.segment( 0, 3 ).transpose() / 1E3 << std::endl;
}

BOOST_AUTO_TEST_CASE( test_polyhedron_set_up )
{
    std::cout << std::endl << "############################################# Individual functions test" << std::endl;
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

    std::cout << "Center of mass: " << gravitySettings.getCenterOfMassPosition().transpose() << std::endl;
    std::cout << "Volume: " << gravitySettings.getVolume() << std::endl;

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
                gravitySettings.getGravitationalParameter() / gravitySettings.getVolume(), perFacetFactor);
        std::cout << "Laplace equation: " << - laplacian / gravitationalConstant / density / M_PI << " pi" << std::endl;

        double potential = gravitation::calculatePolyhedronGravitationalPotential(
                gravitySettings.getGravitationalParameter() / gravitySettings.getVolume(), verticesCoordinatesRelativeToFieldPoint, verticesDefiningEachFacet,
                verticesDefiningEachEdge, facetDyads, edgeDyads, perFacetFactor, perEdgeFactor);
        std::cout << "Potential: " << potential << std::endl;

        Eigen::Vector3d acceleration = gravitation::calculatePolyhedronGradientOfGravitationalPotential(
                gravitySettings.getGravitationalParameter() / gravitySettings.getVolume(), verticesCoordinatesRelativeToFieldPoint, verticesDefiningEachFacet,
                verticesDefiningEachEdge, facetDyads, edgeDyads, perFacetFactor, perEdgeFactor);
        std::cout << "Acceleration: " << acceleration.transpose() << std::endl;

        if ( positionId == 2 )
        {
            Eigen::Matrix3d hessian = gravitation::calculatePolyhedronHessianOfGravitationalPotential(
                gravitySettings.getGravitationalParameter() / gravitySettings.getVolume(), facetDyads, edgeDyads, perFacetFactor, perEdgeFactor);
            std::cout << "Hessian\n: " << hessian.transpose() << std::endl;
        }
    }

}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
