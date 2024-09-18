/*    Copyright (c) 2010-2019, Delft University of Technology
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

#include <limits>

#include <boost/test/unit_test.hpp>


#include <Eigen/Core>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/celestialBodyConstants.h"
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

BOOST_AUTO_TEST_SUITE( test_acceleration_model_setup )

//! Test set up of point mass gravitational accelerations, both direct and third-body.
BOOST_AUTO_TEST_CASE( test_centralGravityModelSetup )
{
    // Load Spice kernel with gravitational parameters.
    spice_interface::loadStandardSpiceKernels( );

    // Create bodies with gravitational parameters from Spice and JPL approximane positions
    // as ephemerides
    BodyListSettings bodySettings;
    bodySettings.addSettings( "Mars" );
    bodySettings.addSettings( "Jupiter" );
    bodySettings.addSettings( "Sun" );
    bodySettings.at( "Mars" )->ephemerisSettings = std::make_shared< ApproximateJplEphemerisSettings >(
               "Mars", 0 );
    bodySettings.at( "Jupiter" )->ephemerisSettings = std::make_shared< ApproximateJplEphemerisSettings >(
                "Jupiter", 0 );
    bodySettings.at( "Mars" )->gravityFieldSettings =
            std::make_shared< GravityFieldSettings >( central_spice );
    bodySettings.at( "Jupiter" )->gravityFieldSettings =
            std::make_shared< GravityFieldSettings >( central_spice );
    bodySettings.at( "Sun" )->gravityFieldSettings =
            std::make_shared< GravityFieldSettings >( central_spice );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    

    // Defins state of Sun to be all zero.
    std::map< double, Eigen::Vector6d > sunStateHistory;
    sunStateHistory[ -1.0E9 ] = Eigen::Vector6d::Zero( );
    sunStateHistory[ 0.0 ] = Eigen::Vector6d::Zero( );
    sunStateHistory[ -1.0E9 ] = Eigen::Vector6d::Zero( );
    std::shared_ptr< interpolators::LinearInterpolator< double, Eigen::Vector6d > >
            sunStateInterpolaotor = std::make_shared<
            interpolators::LinearInterpolator< double, Eigen::Vector6d > >(
                sunStateHistory );
    bodies.at( "Sun" ) ->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                                           sunStateInterpolaotor ) );

    // Update bodies to current state (normally done by numerical integrator).
    for( auto bodyIterator : bodies.getMap( )  )
    {
        bodyIterator.second->setStateFromEphemeris( 1.0E7 );
    }


    // Define settings for accelerations: point  mass atraction by Jupiter and Sun on Mars
    SelectedAccelerationMap accelerationSettingsMap;
    accelerationSettingsMap[ "Mars" ][ "Sun" ].push_back(
                std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationSettingsMap[ "Mars" ][ "Jupiter" ].push_back(
                std::make_shared< AccelerationSettings >( point_mass_gravity ) );

    // Define origin of integration to be barycenter.
    std::map< std::string, std::string > centralBodies;
    centralBodies[ "Mars" ] = "SSB";

    // Create accelerations
    AccelerationMap accelerationsMap = createAccelerationModelsMap(
                bodies, accelerationSettingsMap, centralBodies );

    // Retrieve created accelerations.
    std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > >
            sunAcceleration = accelerationsMap[ "Mars" ][ "Sun" ][ 0 ];
    std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > >
            jupiterAcceleration = accelerationsMap[ "Mars" ][ "Jupiter" ][ 0 ];

    // Create accelerations manually (point mass inertial).
    std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > >
            manualSunAcceleration =
            std::make_shared< gravitation::CentralGravitationalAccelerationModel< > >(
                std::bind( &Body::getPositionByReference, bodies.at( "Mars" ), std::placeholders::_1 ),
                spice_interface::getBodyGravitationalParameter( "Sun" ),
                std::bind( &Body::getPositionByReference, bodies.at( "Sun" ), std::placeholders::_1 ) );

    std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > >
            manualJupiterAcceleration =
            std::make_shared< gravitation::CentralGravitationalAccelerationModel< > >(
                std::bind( &Body::getPositionByReference, bodies.at( "Mars" ), std::placeholders::_1 ),
                spice_interface::getBodyGravitationalParameter( "Jupiter" ),
                std::bind( &Body::getPositionByReference, bodies.at( "Jupiter" ), std::placeholders::_1 ) );

    // Test equivalence of two acceleration models.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                ( basic_astrodynamics::updateAndGetAcceleration( sunAcceleration ) ),
                ( basic_astrodynamics::updateAndGetAcceleration( manualSunAcceleration ) ),
                std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                ( basic_astrodynamics::updateAndGetAcceleration( jupiterAcceleration ) ),
                ( basic_astrodynamics::updateAndGetAcceleration( manualJupiterAcceleration ) ),
                std::numeric_limits< double >::epsilon( ) );

    // Change central body to Sun, which will result in modified accelerations.
    centralBodies[ "Mars" ] = "Sun";

    // Recreate and retrieve accelerations.
    accelerationsMap = createAccelerationModelsMap(
                bodies, accelerationSettingsMap, centralBodies );
    sunAcceleration = accelerationsMap[ "Mars" ][ "Sun" ][ 0 ];
    jupiterAcceleration = accelerationsMap[ "Mars" ][ "Jupiter" ][ 0 ];

    // Manually create Sun's acceleration on Mars, which now include's Mars'gravitational parameter,
    // since the integration is done w.r.t. the Sun, not the barycenter.
    manualSunAcceleration =
            std::make_shared< gravitation::CentralGravitationalAccelerationModel< > >(
                std::bind( &Body::getPositionByReference, bodies.at( "Mars" ), std::placeholders::_1 ),
                spice_interface::getBodyGravitationalParameter( "Sun" ) +
                spice_interface::getBodyGravitationalParameter( "Mars" ),
                std::bind( &Body::getPositionByReference, bodies.at( "Sun" ), std::placeholders::_1 ) );

    // Manually create Jupiter's acceleration on Mars, which now a third body acceleration,
    // with the Sun the central body.
    manualJupiterAcceleration =
            std::make_shared< gravitation::ThirdBodyAcceleration<
            gravitation::CentralGravitationalAccelerationModel< > > >(
                std::make_shared< gravitation::CentralGravitationalAccelerationModel< > >(
                    std::bind( &Body::getPositionByReference, bodies.at( "Mars" ), std::placeholders::_1 ),
                    spice_interface::getBodyGravitationalParameter( "Jupiter" ),
                    std::bind( &Body::getPositionByReference, bodies.at( "Jupiter" ), std::placeholders::_1 ) ),
                std::make_shared< gravitation::CentralGravitationalAccelerationModel< > >(
                    std::bind( &Body::getPositionByReference, bodies.at( "Sun" ), std::placeholders::_1 ),
                    spice_interface::getBodyGravitationalParameter( "Jupiter" ),
                    std::bind( &Body::getPositionByReference, bodies.at( "Jupiter" ), std::placeholders::_1 ) ), "Jupiter" );

    // Test equivalence of two acceleration models.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                ( basic_astrodynamics::updateAndGetAcceleration( sunAcceleration ) ),
                ( basic_astrodynamics::updateAndGetAcceleration( manualSunAcceleration ) ),
                std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                ( basic_astrodynamics::updateAndGetAcceleration( jupiterAcceleration ) ),
                ( basic_astrodynamics::updateAndGetAcceleration( manualJupiterAcceleration ) ),
                std::numeric_limits< double >::epsilon( ) );
}

//! Test set up of spherical harmonic gravitational accelerations.
BOOST_AUTO_TEST_CASE( test_shGravityModelSetup )
{

    // Create system of bodies
    SystemOfBodies bodies;
    bodies.createEmptyBody( "Earth" );
    bodies.createEmptyBody( "Vehicle" );

    // Set constant state for Earth and Vehicle
    Eigen::Vector6d dummyEarthState =
            ( Eigen::Vector6d ( ) << 1.1E11, 0.5E11, 0.01E11, 0.0, 0.0, 0.0
              ).finished( );
    Eigen::Vector7d unitRotationalState = Eigen::Vector7d::Zero( );
    unitRotationalState.segment( 0, 4 ) = linear_algebra::convertQuaternionToVectorFormat(
                Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ) );
    bodies.at( "Earth" )->setState( dummyEarthState );
    bodies.at( "Vehicle" )->setState(
                ( Eigen::Vector6d ( ) << 7.0e6, 8.0e6, 9.0e6, 0.0, 0.0, 0.0
                  ).finished( ) + dummyEarthState );
    bodies.at( "Earth" )->setCurrentRotationalStateToLocalFrame( unitRotationalState );

    // Define Earth gravity field.
    double gravitationalParameter = 3.986004418e14;
    double planetaryRadius = 6378137.0;
    Eigen::MatrixXd cosineCoefficients =
            ( Eigen::MatrixXd( 6, 6 ) <<
              1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              -4.841651437908150e-4, -2.066155090741760e-10, 2.439383573283130e-6, 0.0, 0.0, 0.0,
              9.571612070934730e-7, 2.030462010478640e-6, 9.047878948095281e-7,
              7.213217571215680e-7, 0.0, 0.0, 5.399658666389910e-7, -5.361573893888670e-7,
              3.505016239626490e-7, 9.908567666723210e-7, -1.885196330230330e-7, 0.0,
              6.867029137366810e-8, -6.292119230425290e-8, 6.520780431761640e-7,
              -4.518471523288430e-7, -2.953287611756290e-7, 1.748117954960020e-7
              ).finished( );
    Eigen::MatrixXd sineCoefficients =
            ( Eigen::MatrixXd( 6, 6 ) <<
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 1.384413891379790e-9, -1.400273703859340e-6, 0.0, 0.0, 0.0,
              0.0, 2.482004158568720e-7, -6.190054751776180e-7, 1.414349261929410e-6, 0.0, 0.0,
              0.0, -4.735673465180860e-7, 6.624800262758290e-7, -2.009567235674520e-7,
              3.088038821491940e-7, 0.0, 0.0, -9.436980733957690e-8, -3.233531925405220e-7,
              -2.149554083060460e-7, 4.980705501023510e-8, -6.693799351801650e-7
              ).finished( );
    bodies.at( "Earth" )->setGravityFieldModel(
                std::make_shared< gravitation::SphericalHarmonicsGravityField >(
                    gravitationalParameter, planetaryRadius, cosineCoefficients,
                    sineCoefficients, "IAU_Earth" ) );
    bodies.at( "Earth" )->setRotationalEphemeris(
                std::make_shared< ephemerides::SpiceRotationalEphemeris >(
                    "ECLIPJ2000", "IAU_Earth" ) );


    // Define settings for acceleration model (spherical harmonic due to Earth up to degree and
    // order 5.
    SelectedAccelerationMap accelerationSettingsMap;
    accelerationSettingsMap[ "Vehicle" ][ "Earth" ].push_back(
                std::make_shared< SphericalHarmonicAccelerationSettings >( 5, 5 ) );

    // Set accelerations to be calculated w.r.t. the Earth.
    std::map< std::string, std::string > centralBodies;
    centralBodies[ "Vehicle" ] = "Earth";

    // Create and retrieve acceleration.
    AccelerationMap accelerationsMap = createAccelerationModelsMap(
                bodies, accelerationSettingsMap, centralBodies );
    std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > >
            directAcceleration = accelerationsMap[ "Vehicle" ][ "Earth" ][ 0 ];

    // Manually create acceleration model.
    std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > >
            manualAcceleration =
            std::make_shared< gravitation::SphericalHarmonicsGravitationalAccelerationModel >(
                std::bind( &Body::getPositionByReference, bodies.at( "Vehicle" ), std::placeholders::_1 ),
                gravitationalParameter,
                planetaryRadius, cosineCoefficients, sineCoefficients,
                std::bind( &Body::getPositionByReference, bodies.at( "Earth" ), std::placeholders::_1 ) );

    // Test equivalence of two acceleration models.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                ( basic_astrodynamics::updateAndGetAcceleration( manualAcceleration ) ),
                ( basic_astrodynamics::updateAndGetAcceleration( directAcceleration ) ),
                std::numeric_limits< double >::epsilon( ) );

    // Set (unrealistically) a gravity field model on the Vehicle, to test its
    // influence on acceleration.
    bodies.at( "Vehicle" )->setGravityFieldModel( std::make_shared< gravitation::GravityFieldModel >(
                                                      0.1 * gravitationalParameter ) );

    // Recreate and retrieve acceleration.
    accelerationsMap = createAccelerationModelsMap(
                bodies, accelerationSettingsMap, centralBodies );
    directAcceleration = accelerationsMap[ "Vehicle" ][ "Earth" ][ 0 ];

    // Manually create acceleration.
    manualAcceleration =
            std::make_shared< gravitation::SphericalHarmonicsGravitationalAccelerationModel >(
                std::bind( &Body::getPositionByReference, bodies.at( "Vehicle" ), std::placeholders::_1 ),
                gravitationalParameter * 1.1,
                planetaryRadius, cosineCoefficients, sineCoefficients,
                std::bind( &Body::getPositionByReference, bodies.at( "Earth" ), std::placeholders::_1 ) );

    // Test equivalence of two acceleration models.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                ( basic_astrodynamics::updateAndGetAcceleration( manualAcceleration ) ),
                ( basic_astrodynamics::updateAndGetAcceleration( directAcceleration ) ),
                std::numeric_limits< double >::epsilon( ) );

}

//! Test set up of polyhedron gravitational accelerations.
BOOST_AUTO_TEST_CASE( test_polyhedronGravityModelSetup )
{

    // Create system of bodies
    SystemOfBodies bodies;
    bodies.createEmptyBody( "Earth" );
    bodies.createEmptyBody( "Vehicle" );

    // Set constant state for Earth and Vehicle
    Eigen::Vector6d dummyEarthState =
            ( Eigen::Vector6d ( ) << 1.1E11, 0.5E11, 0.01E11, 0.0, 0.0, 0.0
              ).finished( );
    Eigen::Vector7d unitRotationalState = Eigen::Vector7d::Zero( );
    unitRotationalState.segment( 0, 4 ) = linear_algebra::convertQuaternionToVectorFormat(
                Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ) );
    bodies.at( "Earth" )->setState( dummyEarthState );
    bodies.at( "Vehicle" )->setState(
                ( Eigen::Vector6d ( ) << 7.0e6, 8.0e6, 9.0e6, 0.0, 0.0, 0.0
                  ).finished( ) + dummyEarthState );
    bodies.at( "Earth" )->setCurrentRotationalStateToLocalFrame( unitRotationalState );

    // Define Earth gravity field.
    const double gravitationalParameter = 3.986004418e14;

    // Define cuboid polyhedron dimensions
    const double w = 6378137.0 / 2; // width
    const double h = 6378137.0 / 2; // height
    const double l = 6378137.0 / 2; // length
    const double volume = w * h * l;
    // Define cuboid
    Eigen::MatrixXd verticesCoordinates(8,3);
    verticesCoordinates <<
        0.0, 0.0, 0.0,
        l, 0.0, 0.0,
        0.0, w, 0.0,
        l, w, 0.0,
        0.0, 0.0, h,
        l, 0.0, h,
        0.0, w, h,
        l, w, h;
    Eigen::MatrixXi verticesDefiningEachFacet(12,3);
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

    bodies.at( "Earth" )->setGravityFieldModel(
            std::make_shared< gravitation::PolyhedronGravityField >(
                    gravitationalParameter, verticesCoordinates,
                    verticesDefiningEachFacet, "IAU_Earth" ) );
    bodies.at( "Earth" )->setRotationalEphemeris(
            std::make_shared< ephemerides::SpiceRotationalEphemeris >( "ECLIPJ2000", "IAU_Earth" ) );


    // Define settings for acceleration model (spherical harmonic due to Earth up to degree and
    // order 5.
    SelectedAccelerationMap accelerationSettingsMap;
    accelerationSettingsMap[ "Vehicle" ][ "Earth" ].push_back( polyhedronAcceleration() );

    // Set accelerations to be calculated w.r.t. the Earth.
    std::map< std::string, std::string > centralBodies;
    centralBodies[ "Vehicle" ] = "Earth";

    // Create and retrieve acceleration.
    AccelerationMap accelerationsMap = createAccelerationModelsMap( bodies, accelerationSettingsMap, centralBodies );
    std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > directAcceleration =
            accelerationsMap[ "Vehicle" ][ "Earth" ][ 0 ];

    // Manually create acceleration model.
    gravitation::PolyhedronGravityField manualGravityField = gravitation::PolyhedronGravityField(
            gravitationalParameter, verticesCoordinates, verticesDefiningEachFacet);

    std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > manualAcceleration =
            std::make_shared< gravitation::PolyhedronGravitationalAccelerationModel >(
                    std::bind( &Body::getPositionByReference, bodies.at( "Vehicle" ), std::placeholders::_1 ),
                    gravitationalParameter, volume, verticesCoordinates,
                    verticesDefiningEachFacet, manualGravityField.getVerticesDefiningEachEdge(),
                    manualGravityField.getFacetDyads(), manualGravityField.getEdgeDyads(),
                    std::bind( &Body::getPositionByReference, bodies.at( "Earth" ), std::placeholders::_1 ));

    // Test equivalence of two acceleration models.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                ( basic_astrodynamics::updateAndGetAcceleration( manualAcceleration ) ),
                ( basic_astrodynamics::updateAndGetAcceleration( directAcceleration ) ),
                std::numeric_limits< double >::epsilon( ) );

    // Set (unrealistically) a gravity field model on the Vehicle, to test its influence on acceleration.
    bodies.at( "Vehicle" )->setGravityFieldModel(
            std::make_shared< gravitation::GravityFieldModel >( 0.1 * gravitationalParameter ) );

    // Recreate and retrieve acceleration.
    accelerationsMap = createAccelerationModelsMap( bodies, accelerationSettingsMap, centralBodies );
    directAcceleration = accelerationsMap[ "Vehicle" ][ "Earth" ][ 0 ];

    // Manually create acceleration.
    manualAcceleration = std::make_shared< gravitation::PolyhedronGravitationalAccelerationModel >(
            std::bind( &Body::getPositionByReference, bodies.at( "Vehicle" ), std::placeholders::_1 ),
            gravitationalParameter * 1.1, volume, verticesCoordinates,
            verticesDefiningEachFacet, manualGravityField.getVerticesDefiningEachEdge(),
            manualGravityField.getFacetDyads(), manualGravityField.getEdgeDyads(),
            std::bind( &Body::getPositionByReference, bodies.at( "Earth" ), std::placeholders::_1 ) );

    // Test equivalence of two acceleration models.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                ( basic_astrodynamics::updateAndGetAcceleration( manualAcceleration ) ),
                ( basic_astrodynamics::updateAndGetAcceleration( directAcceleration ) ),
                std::numeric_limits< double >::epsilon( ) );

}

//! Test radiation pressure acceleration
BOOST_AUTO_TEST_CASE( test_radiationPressureAcceleration )
{
    using namespace tudat::simulation_setup;
    using namespace tudat;

    // Load Spice kernels
    spice_interface::loadStandardSpiceKernels( );

    // Get settings for celestial bodies
    BodyListSettings bodySettings;
    bodySettings.addSettings( getDefaultSingleBodySettings( "Earth", 0.0, 10.0 * 86400.0 ), "Earth" );
    bodySettings.addSettings( getDefaultSingleBodySettings( "Sun", 0.0, 10.0 * 86400.0 ), "Sun" );

    // Get settings for vehicle
    double area = 2.34;
    double coefficient = 1.2;
    double bodyMass = 500.0;
    bodySettings.addSettings( "Vehicle" );
    bodySettings.at("Vehicle")->constantMass = bodyMass;
    bodySettings.at( "Vehicle" )->radiationPressureTargetModelSettings =
            std::make_shared<CannonballRadiationPressureTargetModelSettings>(area, coefficient);
    bodySettings.at( "Vehicle" )->ephemerisSettings =
            std::make_shared< KeplerEphemerisSettings >(
                ( Eigen::Vector6d( ) << 12000.0E3, 0.13, 0.3, 0.0, 0.0, 0.0 ).finished( ),
                0.0, spice_interface::getBodyGravitationalParameter( "Earth" ), "Earth", "ECLIPJ2000" );

    // Create bodies
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    

    // Define settings for accelerations
    SelectedAccelerationMap accelerationSettingsMap;
    accelerationSettingsMap[ "Vehicle" ][ "Sun" ].push_back(radiationPressureAcceleration());

    // Define origin of integration
    std::map< std::string, std::string > centralBodies;
    centralBodies[ "Vehicle" ] = "Earth";

    // Create accelerations
    AccelerationMap accelerationsMap = createAccelerationModelsMap(
                bodies, accelerationSettingsMap, centralBodies );
    auto radiationPressureAcceleration =
            std::dynamic_pointer_cast<electromagnetism::IsotropicPointSourceRadiationPressureAcceleration>(
                    accelerationsMap[ "Vehicle" ][ "Sun" ][ 0 ]);

    // Set (arbitrary) test time.
    double testTime = 5.0 * 86400.0;

    // Update environment to current time.
    bodies.at( "Sun" )->setStateFromEphemeris< double, double >( testTime );
    bodies.at( "Earth" )->setStateFromEphemeris< double, double >( testTime );
    bodies.at( "Vehicle" )->setStateFromEphemeris< double, double >( testTime );
    bodies.at("Sun")->getRadiationSourceModel()->updateMembers(testTime);
    bodies.at( "Vehicle" )->getRadiationPressureTargetModel()->updateMembers(testTime);
    radiationPressureAcceleration->updateMembers(testTime);

    // Get acceleration
    Eigen::Vector3d calculatedAcceleration = radiationPressureAcceleration->getAcceleration();
    double calculatedReceivedIrradiance = radiationPressureAcceleration->getReceivedIrradiance();
    double calculatedReceivedFraction = radiationPressureAcceleration->getSourceToTargetReceivedFraction();

    // Manually calculate acceleration
    Eigen::Vector3d sourceTargetRelativePosition =
            ( bodies.at( "Vehicle" )->getState() -  bodies.at( "Sun" )->getState() ).segment( 0, 3 );
    auto sourceDistance = sourceTargetRelativePosition.norm();
    auto irradiance =
                celestial_body_constants::SUN_LUMINOSITY / (4 * mathematical_constants::PI * sourceDistance * sourceDistance);
    auto radiationPressure = irradiance / physical_constants::SPEED_OF_LIGHT;
    auto expectedForceMagnitude = radiationPressure * area * coefficient;
    auto expectedForceDirection = sourceTargetRelativePosition.normalized();
    Eigen::Vector3d expectedAcceleration = expectedForceDirection * expectedForceMagnitude / bodyMass;

    // Compare results
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                expectedAcceleration, calculatedAcceleration, ( 2.0 * std::numeric_limits< double >::epsilon( ) ) );
    BOOST_CHECK_CLOSE(irradiance, calculatedReceivedIrradiance, 1e-15);
    BOOST_CHECK_CLOSE(1.0, calculatedReceivedFraction, 1e-15);
}

//! Test setup of aerodynamic accelerations (constant coefficients)
BOOST_AUTO_TEST_CASE( test_aerodynamicAccelerationModelSetup )
{
    using namespace tudat::simulation_setup;
    using namespace tudat;

    // Load Spice kernels
    spice_interface::loadStandardSpiceKernels( );

    // Test creation with coefficients positive/negative in body/aerodynamic frame (4 cases).
    for( unsigned int testCase = 0; testCase < 4; testCase++ )
    {
        // Get settings for Earth.
        BodyListSettings bodySettings;
        bodySettings.addSettings( getDefaultSingleBodySettings( "Earth", 0.0, 86400.0 ), "Earth" );
        bodySettings.addSettings( "Vehicle" );

        // Define (arbitrary) aerodynamic coefficient settings.
        Eigen::Vector3d aerodynamicCoefficients =
                ( Eigen::Vector3d( ) << 1.0, 3.0, -2.0 ).finished( );
        double referenceArea = 4.7;

        // Get coefficient settings for current test case
        bool areCoefficientsInAerodynamicFrame = 0;
        bool areCoefficientsInNegativeAxisDirection = 1;
        if( testCase < 2 )
        {
            areCoefficientsInAerodynamicFrame = 1;
        }
        else
        {
            areCoefficientsInAerodynamicFrame = 0;
        }
        if( testCase % 2 == 0 )
        {
            areCoefficientsInNegativeAxisDirection = 1;
        }
        else
        {
            areCoefficientsInNegativeAxisDirection = 0;
        }

        bodySettings.at( "Vehicle" )->aerodynamicCoefficientSettings =
                std::make_shared< ConstantAerodynamicCoefficientSettings >(
                    1.0, referenceArea, Eigen::Vector3d::Zero( ), aerodynamicCoefficients, Eigen::Vector3d::Zero( ),
                    getAerodynamicCoefficientFrame( areCoefficientsInAerodynamicFrame, areCoefficientsInNegativeAxisDirection ),
                    getAerodynamicCoefficientFrame( areCoefficientsInAerodynamicFrame, areCoefficientsInNegativeAxisDirection ) );

        // Create body objects.
        SystemOfBodies bodies = createSystemOfBodies( bodySettings );
        
        double angleOfAttack = 1.232;
        double angleOfSideslip = -0.00322;
        double bankAngle = 2.323432;

        std::shared_ptr< ephemerides::AerodynamicAngleRotationalEphemeris > vehicleRotationModel =
                createAerodynamicAngleBasedRotationModel(
                                "Vehicle", "Earth", bodies,
                                "ECLIPJ2000", "VehicleFixed" );
        vehicleRotationModel->setAerodynamicAngleFunction(
                    [=]( const double ){ return ( Eigen::Vector3d( ) << angleOfAttack, angleOfSideslip, bankAngle ).finished( ); } );
        bodies.at( "Vehicle" )->setRotationalEphemeris( vehicleRotationModel );



        // Define settings for accelerations
        SelectedAccelerationMap accelerationSettingsMap;
        accelerationSettingsMap[ "Vehicle" ][ "Earth" ].push_back(
                    std::make_shared< AccelerationSettings >( aerodynamic ) );

        // Define origin of integration
        std::map< std::string, std::string > centralBodies;
        centralBodies[ "Vehicle" ] = "Earth";

        // Create accelerations
        AccelerationMap accelerationsMap = createAccelerationModelsMap(
                    bodies, accelerationSettingsMap, centralBodies );
        std::shared_ptr< AccelerationModel3d > aerodynamicAcceleration = accelerationsMap[ "Vehicle" ][ "Earth" ][ 0 ];

        // Define expected aerodynamic angles (see testAerodynamicAngleCalculator)
        double testHeadingAngle = 1.229357188236127;
        double testFlightPathAngle = -0.024894033070522;
        double testLatitude = -0.385027359562548;
        double testLongitude = -1.849449608688977;


        // Set vehicle body-fixed state (see testAerodynamicAngleCalculator)
        Eigen::Vector6d vehicleBodyFixedState =
                ( Eigen::Vector6d( ) << -1656517.23153109, -5790058.28764025, -2440584.88186829,
                  6526.30784888051, -2661.34558272018, 2377.09572383163 ).finished( );

        double testTime = 0.5E7;

        // Convert vehicle state to inertial frame.
        Eigen::Vector6d vehicleInertialState =
                ephemerides::transformStateToFrameFromRotations(
                    vehicleBodyFixedState,
                    bodies.at( "Earth" )->getRotationalEphemeris( )->getRotationToBaseFrame( testTime ),
                    bodies.at( "Earth" )->getRotationalEphemeris( )->getDerivativeOfRotationToBaseFrame( testTime ) );

        // Set states in environment.
        vehicleRotationModel->setIsBodyInPropagation( true );
        bodies.at( "Earth" )->setState( Eigen::Vector6d::Zero( ) );
        bodies.at( "Earth" )->setCurrentRotationalStateToLocalFrameFromEphemeris( testTime );
        bodies.at( "Vehicle" )->setState( vehicleInertialState );

        // Set vehicle mass
        double bodyMass = 500.0;
        bodies.at( "Vehicle" )->setBodyMassFunction( [ & ]( const double ){ return bodyMass; } );
        bodies.at( "Vehicle" )->updateMass( testTime );

        // Update flight conditions.
        std::shared_ptr< aerodynamics::FlightConditions > vehicleFlightConditions =
                bodies.at( "Vehicle" )->getFlightConditions( );
        vehicleFlightConditions->updateConditions( testTime );

        // Check whether flight conditions object has been correctly automatically created
        // (see testAerodynamicAngleCalculator)
        {
            BOOST_CHECK_SMALL(
                        std::fabs( vehicleFlightConditions->getAerodynamicAngleCalculator( )->
                                   getAerodynamicAngle( latitude_angle ) -
                                   testLatitude), 10.0 * std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_SMALL(
                        std::fabs( vehicleFlightConditions->getAerodynamicAngleCalculator( )->
                                   getAerodynamicAngle( longitude_angle ) -
                                   testLongitude), 10.0 * std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_SMALL(
                        std::fabs( vehicleFlightConditions->getAerodynamicAngleCalculator( )->
                                   getAerodynamicAngle( heading_angle ) -
                                   testHeadingAngle), 10.0 * std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_SMALL(
                        std::fabs( vehicleFlightConditions->getAerodynamicAngleCalculator( )->
                                   getAerodynamicAngle( flight_path_angle ) -
                                   testFlightPathAngle), 10.0 * std::numeric_limits< double >::epsilon( ) );

            BOOST_CHECK_SMALL(
                        std::fabs( vehicleFlightConditions->getAerodynamicAngleCalculator( )->
                                   getAerodynamicAngle( angle_of_attack ) -
                                   angleOfAttack ), 10.0 * std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_SMALL(
                        std::fabs( vehicleFlightConditions->getAerodynamicAngleCalculator( )->
                                   getAerodynamicAngle( angle_of_sideslip ) -
                                   angleOfSideslip), 10.0 * std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_SMALL(
                        std::fabs( vehicleFlightConditions->getAerodynamicAngleCalculator( )->
                                   getAerodynamicAngle( bank_angle ) -
                                   bankAngle), 10.0 * std::numeric_limits< double >::epsilon( ) );

            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        vehicleFlightConditions->getCurrentBodyCenteredBodyFixedState( ), vehicleBodyFixedState,
                        ( 2.0 * std::numeric_limits< double >::epsilon( ) ) );
        }

        // Define current frame for aerodynamic coefficients.
        AerodynamicsReferenceFrames coefficientFrame = aerodynamic_frame;
        if( testCase < 2 )
        {
            coefficientFrame = aerodynamic_frame;
        }
        else
        {
            coefficientFrame = body_frame;
        }

        // Get rotation from coefficient to propagation frame.
        Eigen::Quaterniond rotationToPropagationFrame =
                bodies.at( "Earth" )->getCurrentRotationToGlobalFrame( ) *
                vehicleFlightConditions->getAerodynamicAngleCalculator( )->getRotationQuaternionBetweenFrames(
                    coefficientFrame, corotating_frame );

        // Calculate aerodynamic force manually.
        double dynamicPressure = 0.5 * bodies.at( "Earth" )->getAtmosphereModel( )->getDensity(
                    vehicleFlightConditions->getCurrentAltitude( ), testLongitude, testLatitude, testTime ) *
                std::pow( vehicleBodyFixedState.segment( 3, 3 ).norm( ), 2.0 );
        Eigen::Vector3d expectedAerodynamicForce = dynamicPressure * referenceArea  *
                ( rotationToPropagationFrame * aerodynamicCoefficients );
        if( testCase % 2 == 0 )
        {
            expectedAerodynamicForce *= -1.0;
        }

        // Get automatic aerodynamic force.
        Eigen::Vector3d computedAerodynamicForce =
                bodyMass * updateAndGetAcceleration< Eigen::Vector3d >( aerodynamicAcceleration );

        // Compare results
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    expectedAerodynamicForce, computedAerodynamicForce, ( 10.0 *  std::numeric_limits< double >::epsilon( ) ) );
    }
}

//! Test setup of aerodynamic accelerations (non-constant coefficients)
BOOST_AUTO_TEST_CASE( test_aerodynamicAccelerationModelSetupWithCoefficientIndependentVariables )
{
    using namespace tudat::simulation_setup;
    using namespace tudat;

    // Load Spice kernels
    spice_interface::loadStandardSpiceKernels( );
    // Get settings for Earth.
    BodyListSettings bodySettings;
    bodySettings.addSettings( getDefaultSingleBodySettings( "Earth", 0.0, 86400.0 ), "Earth" );
    bodySettings.addSettings( "Vehicle" );

    // Create body objects.
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    

    // Create vehicle aerodynamic coefficients
    bodies.at( "Vehicle" )->setAerodynamicCoefficientInterface(
                getApolloCoefficientInterface( ) );

    // Define settings for accelerations: point
    SelectedAccelerationMap accelerationSettingsMap;
    accelerationSettingsMap[ "Vehicle" ][ "Earth" ].push_back(
                std::make_shared< AccelerationSettings >( aerodynamic ) );

    // Define origin of integration
    std::map< std::string, std::string > centralBodies;
    centralBodies[ "Vehicle" ] = "Earth";

    // Create accelerations
    AccelerationMap accelerationsMap = createAccelerationModelsMap(
                bodies, accelerationSettingsMap, centralBodies );
    std::shared_ptr< AccelerationModel3d > aerodynamicAcceleration = accelerationsMap[ "Vehicle" ][ "Earth" ][ 0 ];


    // Retrieve flight conditions and orientation angles
    std::shared_ptr< aerodynamics::FlightConditions > vehicleFlightConditions =
            bodies.at( "Vehicle" )->getFlightConditions( );
    std::shared_ptr< aerodynamics::AerodynamicCoefficientInterface > coefficientInterface =
            bodies.at( "Vehicle" )->getAerodynamicCoefficientInterface( );

    // Define orientation angles.
    double angleOfAttack = 1.232;
    double angleOfSideslip = -0.00322;
    double bankAngle = 2.323432;


    std::shared_ptr< ephemerides::AerodynamicAngleRotationalEphemeris > vehicleRotationModel =
            createAerodynamicAngleBasedRotationModel(
                            "Vehicle", "Earth", bodies,
                            "ECLIPJ2000", "VehicleFixed" );
    vehicleRotationModel->setAerodynamicAngleFunction(
                [=]( const double ){ return ( Eigen::Vector3d( ) << angleOfAttack, angleOfSideslip, bankAngle ).finished( ); } );
    bodies.at( "Vehicle" )->setRotationalEphemeris( vehicleRotationModel );


    // Update environment to current time.
    double testTime = 0.5E7;
    bodies.at( "Earth" )->setCurrentRotationalStateToLocalFrameFromEphemeris( testTime );
    bodies.at( "Earth" )->setState( Eigen::Vector6d::Zero( ) );
    double bodyMass = 500.0;

    // Set vehicle mass
    bodies.at( "Vehicle" )->setBodyMassFunction( [ & ]( const double ){ return bodyMass; } );
    bodies.at( "Vehicle" )->updateMass( testTime );

    // Test aerodynamic coefficients for various cases of independent variables.
    for( unsigned int i = 0; i < 4; i++ )
    {
        // Define body-fixed vehicle state.
        Eigen::Vector6d vehicleBodyFixedState =
                ( Eigen::Vector6d( ) << -1656517.23153109, -5790058.28764025, -2440584.88186829,
                  6526.30784888051, -2661.34558272018, 2377.09572383163 ).finished( );
        if( i > 0 )
        {
            vehicleBodyFixedState.segment( 3, 3 ) += ( Eigen::Vector3d( ) << -3234.2, 2456.2, 33.245 ).finished( );
        }

        // Define vehicle inertial state.
        Eigen::Vector6d vehicleInertialState =
                ephemerides::transformStateToFrameFromRotations(
                    vehicleBodyFixedState,
                    bodies.at( "Earth" )->getRotationalEphemeris( )->getRotationToBaseFrame( testTime ),
                    bodies.at( "Earth" )->getRotationalEphemeris( )->getDerivativeOfRotationToBaseFrame( testTime ) );
        bodies.at( "Vehicle" )->setState( vehicleInertialState );

        // Define orientation angles.
        if( i > 1 )
        {
            angleOfAttack -= 0.6465;
        }
        if( i > 2  )
        {
            angleOfSideslip += 0.00123;
        }

        vehicleRotationModel->setIsBodyInPropagation( true );
        vehicleRotationModel->resetCurrentTime( );
        vehicleRotationModel->setAerodynamicAngleFunction(
                    [=]( const double ){ return ( Eigen::Vector3d( ) << angleOfAttack, angleOfSideslip, bankAngle ).finished( ); } );

        // Update flight conditions
        vehicleFlightConditions->resetCurrentTime( );
        vehicleFlightConditions->updateConditions( testTime );

        // Calculate Mach number
        double velocity = vehicleBodyFixedState.segment( 3, 3 ).norm( );
        double speedOfSound = bodies.at( "Earth" )->getAtmosphereModel( )->getSpeedOfSound(
                    vehicleFlightConditions->getCurrentAltitude( ), 0.0, 0.0, 0.0 );
        double machNumber = velocity / speedOfSound;

        // Get manual and automatic coefficients and compare.
        Eigen::Vector3d automaticCoefficients = coefficientInterface->getCurrentForceCoefficients( );
        coefficientInterface->updateFullCurrentCoefficients(
        { machNumber, angleOfAttack, angleOfSideslip } );
        Eigen::Vector3d manualCoefficients = coefficientInterface->getCurrentForceCoefficients( );


        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    automaticCoefficients, manualCoefficients, ( 5.0 *  std::numeric_limits< double >::epsilon( ) ) );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat


