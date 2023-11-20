#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/basics/testMacros.h"
#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/interface/spice/spiceInterface.h"

#include "tudat/astro/ephemerides/keplerEphemeris.h"
#include "tudat/astro/ephemerides/simpleRotationalEphemeris.h"
#include "tudat/astro/ephemerides/constantRotationalEphemeris.h"
#include "tudat/simulation/environment_setup/createBodies.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/propagation_setup/createAccelerationModels.h"

namespace tudat
{

using namespace tudat::basic_astrodynamics;
using namespace tudat::simulation_setup;
using namespace tudat::ephemerides;
using namespace tudat::basic_astrodynamics;
using namespace tudat::electromagnetism;

namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_PanelledRadiationPressure )

//! Test panelled radiation acceleration model with various simple boxes-and-wings models.
BOOST_AUTO_TEST_CASE( testSimpleGeometryPanelledRadiationPressure )
{

    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Create bodies needed in simulation
    double initialEphemerisTime = -86400.0;
    double finalEphemerisTime = 1.1 * 365.25 * 86400.0;
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Sun" );
    SystemOfBodies bodies = createSystemOfBodies(
                getDefaultBodySettings( bodyNames,initialEphemerisTime, finalEphemerisTime ) );

    // Create vehicle
    double vehicleMass = 2500.0;
    bodies.createEmptyBody( "Vehicle" );
    bodies.at( "Vehicle" )->setConstantBodyMass( vehicleMass );

    // Put vehicle on circular orbit around Sun
    Eigen::Vector6d initialStateInKeplerianElements = Eigen::Vector6d::Zero( );
    initialStateInKeplerianElements[ 0 ] = physical_constants::ASTRONOMICAL_UNIT;
    bodies.at( "Vehicle" )->setEphemeris(
                std::make_shared< KeplerEphemeris >(
                    initialStateInKeplerianElements, 0.0, spice_interface::getBodyGravitationalParameter( "Sun" ),
                    "Sun", "ECLIPJ2000", 1 ) );
    bodies.processBodyFrameDefinitions( );

    // Define constant rotational ephemeris
    Eigen::Vector7d rotationalStateVehicle;
    rotationalStateVehicle.segment( 0, 4 ) = linear_algebra::convertQuaternionToVectorFormat( Eigen::Quaterniond( Eigen::Matrix3d::Identity() ));
    rotationalStateVehicle.segment( 4, 3 ) = Eigen::Vector3d::Zero();
    bodies.at( "Vehicle" )->setRotationalEphemeris(
                std::make_shared< ConstantRotationalEphemeris >(
                    rotationalStateVehicle, "ECLIPJ2000", "VehicleFixed" ) );



    for( unsigned int test = 0; test <= 2; test++ )
    {
        // Create radiation pressure properties
        std::vector< double > areas;
        areas.push_back( 2.532 );
        areas.push_back( 3.254 );
        areas.push_back( 8.654 );
        areas.push_back( 1.346 );
        areas.push_back( 2.454 );
        areas.push_back( 5.345 );

        std::vector< double > emissivities;
        if( test == 0 || test == 2 )
        {
            emissivities.push_back( 0.0 );
            emissivities.push_back( 1.0 );
            emissivities.push_back( 0.0 );
            emissivities.push_back( 1.0 );
        }
        else if( test == 1 )
        {
            emissivities.push_back( 0.2 );
            emissivities.push_back( 0.3 );
            emissivities.push_back( 0.4 );
            emissivities.push_back( 0.5 );
        }
        emissivities.push_back( 1.0 );
        emissivities.push_back( 0.0 );


        std::vector< double > diffuseReflectionCoefficients;
        if( test == 0 || test == 2 )
        {
            diffuseReflectionCoefficients.push_back( 0.0 );
            diffuseReflectionCoefficients.push_back( 0.0 );
            diffuseReflectionCoefficients.push_back( 0.0 );
            diffuseReflectionCoefficients.push_back( 0.0 );
            diffuseReflectionCoefficients.push_back( 0.0 );
            diffuseReflectionCoefficients.push_back( 0.0 );
        }
        else
        {
            diffuseReflectionCoefficients.push_back( 0.6 );
            diffuseReflectionCoefficients.push_back( 0.5 );
            diffuseReflectionCoefficients.push_back( 0.4 );
            diffuseReflectionCoefficients.push_back( 0.3 );
            diffuseReflectionCoefficients.push_back( 0.2 );
            diffuseReflectionCoefficients.push_back( 0.1 );
        }

        if( test == 2 )
        {
            // Define simple rotational ephemeris
            bodies.at( "Vehicle" )->setRotationalEphemeris(
                        std::make_shared< tudat::ephemerides::SimpleRotationalEphemeris >(
                            0.2, 0.4, -0.2, 1.0E-5, 0.0, "ECLIPJ2000", "VehicleFixed" ) );
        }

        std::vector< Eigen::Vector3d > panelSurfaceNormals;
        panelSurfaceNormals.push_back( -Eigen::Vector3d::UnitX( ) );
        panelSurfaceNormals.push_back( -Eigen::Vector3d::UnitY( ) );
        panelSurfaceNormals.push_back( Eigen::Vector3d::UnitX( ) );
        panelSurfaceNormals.push_back( Eigen::Vector3d::UnitY( ) );
        panelSurfaceNormals.push_back( Eigen::Vector3d::UnitZ( ) );
        panelSurfaceNormals.push_back( -Eigen::Vector3d::UnitZ( ) );

        std::shared_ptr< PanelledRadiationPressureInterfaceSettings > radiationPressureInterfaceSettings =
                std::make_shared< PanelledRadiationPressureInterfaceSettings >(
                    "Sun", emissivities, areas, diffuseReflectionCoefficients, panelSurfaceNormals );
        std::shared_ptr< PanelledRadiationPressureInterface > radiationPressureInterface =
                std::dynamic_pointer_cast< PanelledRadiationPressureInterface >(
                    createRadiationPressureInterface( radiationPressureInterfaceSettings, "Vehicle", bodies ) );
        bodies.at( "Vehicle" )->setRadiationPressureInterface( "Sun", radiationPressureInterface );


        // Define accelerations
        SelectedAccelerationMap accelerationMap;
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOnVehicle;
        accelerationsOnVehicle[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                       panelled_radiation_pressure_acceleration ) );
        accelerationMap[ "Vehicle" ] = accelerationsOnVehicle;
        std::map< std::string, std::string > centralBodyMap;
        centralBodyMap[ "Vehicle" ] = "Sun";
        AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodies, accelerationMap, centralBodyMap );
        std::shared_ptr< AccelerationModel< Eigen::Vector3d > > accelerationModel =
                accelerationModelMap.at( "Vehicle" ).at( "Sun" ).at( 0 );

        // Compute spacecraft orbital period, and compute test times
        double orbitalPeriod =
                2.0 * mathematical_constants::PI * std::sqrt(
                    std::pow( physical_constants::ASTRONOMICAL_UNIT, 3.0 ) /
                    spice_interface::getBodyGravitationalParameter( "Sun" ) );
        std::vector< double > testTimes = { 0.0, orbitalPeriod / 4.0, orbitalPeriod / 2.0, 3.0 * orbitalPeriod / 4.0 };

        // Compute panelled radiation pressure for various relative Sun positions
        Eigen::Vector3d calculatedAcceleration, expectedAcceleration;
        Eigen::Vector3d sunCenteredVehiclePosition;

        std::shared_ptr< Ephemeris > vehicleEphemeris = bodies.at( "Vehicle" )->getEphemeris( );
        for( unsigned int i = 0; i < testTimes.size( ); i++ )
        {
            // Update environment and acceleration
            bodies.at( "Sun" )->setStateFromEphemeris( testTimes[ i ] );
            bodies.at( "Vehicle" )->setStateFromEphemeris( testTimes[ i ] );
            bodies.at( "Vehicle" )->setCurrentRotationToLocalFrameFromEphemeris( testTimes[ i ] );
            bodies.at( "Vehicle" )->updateMass( testTimes[ i ] );
            radiationPressureInterface->updateInterface( testTimes[ i ] );
            accelerationModel->updateMembers( testTimes[ i ] );

            // Retrieve acceleration
            calculatedAcceleration = accelerationModel->getAcceleration( );

            // Manually compute acceleration
            double radiationPressure = radiationPressureInterface->getCurrentRadiationPressure( );
            sunCenteredVehiclePosition = vehicleEphemeris->getCartesianState( testTimes[ i ] ).segment( 0, 3 );

            // Case with emissivities equal to either 0 or 1 and no diffuse reflection.
            if( test == 0 )
            {
                if( i == 0 || i == 2 )
                {
                    expectedAcceleration = radiationPressure / vehicleMass  *
                            areas.at( i ) * sunCenteredVehiclePosition.normalized( );
                }
                else if( i == 1 || i == 3 )
                {
                    expectedAcceleration = -radiationPressure / vehicleMass * (
                                2.0 * areas.at( i ) * panelSurfaceNormals.at( i ) );
                }
            }

            // Case with more complex panel characteristics.
            else if ( test == 1 )
            {
                expectedAcceleration = radiationPressure / vehicleMass *
                        ( 1.0 + emissivities.at( i ) + 2.0 * diffuseReflectionCoefficients.at( i ) / 3.0 ) *
                        areas.at( i ) * sunCenteredVehiclePosition.normalized( );
            }

            if( test == 0 || test == 1 )
            {
                // Check computed acceleration
                for( unsigned int j = 0; j < 3; j++ )
                {
                    BOOST_CHECK_SMALL( std::fabs( calculatedAcceleration( j ) - expectedAcceleration( j ) ), 2.0E-23 );
                }
            }


            // Case with non constant rotational ephemeris for the spacecraft
            else
            {
                std::shared_ptr< PanelledRadiationPressureAcceleration > panelledRadiationPressureAcceleration =
                        std::dynamic_pointer_cast< PanelledRadiationPressureAcceleration >( accelerationModel );
                Eigen::Quaterniond currentRotationToInertialFrame =
                        bodies.at( "Vehicle" )->getRotationalEphemeris( )->getRotationToBaseFrame( testTimes.at( i ) );

                expectedAcceleration = Eigen::Vector3d::Zero();

                for( unsigned int j = 0; j < panelSurfaceNormals.size( ); j++ )
                {
                   Eigen::Vector3d currentPanelNormal =
                           panelledRadiationPressureAcceleration->getCurrentPanelSurfaceNormalInPropagationFrame( j );
                   Eigen::Vector3d expectedPanelNormal =
                           currentRotationToInertialFrame * panelSurfaceNormals.at( j );

                   // Check surface normals orientation.
                   for( unsigned int k = 0; k < 3; k++ )
                   {
                       BOOST_CHECK_SMALL( std::fabs( currentPanelNormal( k ) - expectedPanelNormal( k ) ), 1.0E-15 );
                   }


                   // Compute expected acceleration.
                    double cosinusPanelInclination =  - sunCenteredVehiclePosition.normalized( ).dot( currentPanelNormal );

                    if ( cosinusPanelInclination > 0.0 ){
                        expectedAcceleration += - radiationPressure / vehicleMass  * areas.at( j ) * cosinusPanelInclination *
                                ( ( 1.0 - emissivities.at( j ) ) * - sunCenteredVehiclePosition.normalized(  )
                                    + 2.0 * emissivities.at( j ) * cosinusPanelInclination * currentPanelNormal );
                    }
                }

                // Check computed acceleration
                for( unsigned int j = 0; j < 3; j++ )
                {
                    BOOST_CHECK_SMALL( std::fabs( calculatedAcceleration( j ) - expectedAcceleration( j ) ), 2.0E-23 );
                }

            }
        }
    }

}


//! Test panelled radiation acceleration model for a spacecraft in various orbits with respect to the Sun.
BOOST_AUTO_TEST_CASE( testPanelledRadiationPressureMontenbruckModel )
{

    // Box-and-wings model is partially obtained from Montenbruck, O., 2017.
    // Semi-analytical solar radiation pressure modeling for QZS-1 orbit-normal and yaw-steering attitude.
    // Advances in Space Research, 59(8), pp.2088-2100..

    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Create bodies needed in simulation
    double initialEphemerisTime = 0.0;
    double finalEphemerisTime = 1.1 * 365.25 * 86400.0;
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Sun" );
    SystemOfBodies bodies = createSystemOfBodies(
                getDefaultBodySettings( bodyNames,initialEphemerisTime, finalEphemerisTime ) );

    // Create vehicle
    double vehicleMass = 2000.0;
    bodies.createEmptyBody( "Vehicle" );
    bodies.at( "Vehicle" )->setConstantBodyMass( vehicleMass );


    for ( int testCase = 0 ; testCase < 4 ; testCase++){

        // Put vehicle on circular orbit around the Sun with i = 0 deg
        if ( testCase == 0 || testCase == 3 ){
            Eigen::Vector6d initialStateInKeplerianElements = Eigen::Vector6d::Zero( );
            initialStateInKeplerianElements[ 0 ] = physical_constants::ASTRONOMICAL_UNIT;
            bodies.at( "Vehicle" )->setEphemeris( std::make_shared< KeplerEphemeris >( initialStateInKeplerianElements, 0.0,
                                                         spice_interface::getBodyGravitationalParameter( "Sun" ), "Sun", "ECLIPJ2000", 1 ) );
        }

        // Put vehicle in circular orbit around the Sun with i = 90 deg
        else if ( testCase == 1 ){
            Eigen::Vector6d initialStateInKeplerianElements = Eigen::Vector6d::Zero( );
            initialStateInKeplerianElements[ 0 ] = physical_constants::ASTRONOMICAL_UNIT;
            initialStateInKeplerianElements[ orbital_element_conversions::inclinationIndex ] = unit_conversions::convertDegreesToRadians( 90.0 );
            bodies.at( "Vehicle" )->setEphemeris( std::make_shared< KeplerEphemeris >( initialStateInKeplerianElements, 0.0,
                                                         spice_interface::getBodyGravitationalParameter( "Sun" ), "Sun", "ECLIPJ2000", 1 ) );
        }

        // Put vehicle in circular orbit around the Sun with arbitrary chosen inclination
        else if ( testCase == 2 ){
            Eigen::Vector6d initialStateInKeplerianElements = Eigen::Vector6d::Zero( );
            initialStateInKeplerianElements[ 0 ] = physical_constants::ASTRONOMICAL_UNIT;
            initialStateInKeplerianElements[ orbital_element_conversions::inclinationIndex ] = unit_conversions::convertDegreesToRadians( 20.0 );
            bodies.at( "Vehicle" )->setEphemeris( std::make_shared< KeplerEphemeris >( initialStateInKeplerianElements, 0.0,
                                                         spice_interface::getBodyGravitationalParameter( "Sun" ), "Sun", "ECLIPJ2000", 1 ) );
        }


        // Set-up rotational ephemeris for vehicle.
        if ( testCase < 3 ){
            // Define constant rotational model.
            Eigen::Vector7d rotationalStateVehicle;
            rotationalStateVehicle.segment( 0, 4 ) = linear_algebra::convertQuaternionToVectorFormat( Eigen::Quaterniond( Eigen::Matrix3d::Identity() ));
            rotationalStateVehicle.segment( 4, 3 ) = Eigen::Vector3d::Zero();
            bodies.at( "Vehicle" )->setRotationalEphemeris( std::make_shared< ConstantRotationalEphemeris >( rotationalStateVehicle, "ECLIPJ2000",
                                                                                                           "VehicleFixed" ) );
        }
        else if ( testCase == 3 ){
            // Define simple rotational model.
            bodies.at( "Vehicle" )->setRotationalEphemeris( std::make_shared< tudat::ephemerides::SimpleRotationalEphemeris >(
                            0.2, 0.4, -0.2, 1.0E-5, 0.0, "ECLIPJ2000", "VehicleFixed" ) );
        }
        bodies.processBodyFrameDefinitions( );





        std::vector< double > areas;
        std::vector< double > emissivities;
        std::vector< double > diffuseReflectionCoefficients;
        std::vector< Eigen::Vector3d > panelSurfaceNormals;


        // Define panels properties for test cases with constant rotational model.
        if ( testCase < 3 ){

            // Panels areas
            areas.push_back( 2.0 );
            areas.push_back( 4.0 );
            areas.push_back( 6.0 );
            areas.push_back( 9.9 );
            areas.push_back( 2.3 );
            areas.push_back( 9.9 );
            areas.push_back( 2.3 );
            areas.push_back( 4.6 );
            areas.push_back( 2.7 );
            areas.push_back( 5.8 );
            areas.push_back( 2.7 );

            // Emissivities
            emissivities.push_back( 0.0 );
            emissivities.push_back( 0.1 );
            emissivities.push_back( 0.0 );
            emissivities.push_back( 0.0 );
            emissivities.push_back( 0.1 );
            emissivities.push_back( 0.0 );
            emissivities.push_back( 0.1 );
            emissivities.push_back( 0.0 );
            emissivities.push_back( 0.1 );
            emissivities.push_back( 0.0 );
            emissivities.push_back( 0.1 );

            //Diffuse reflection coefficients
            diffuseReflectionCoefficients.push_back( 0.06 );
            diffuseReflectionCoefficients.push_back( 0.46 );
            diffuseReflectionCoefficients.push_back( 0.06 );
            diffuseReflectionCoefficients.push_back( 0.06 );
            diffuseReflectionCoefficients.push_back( 0.46 );
            diffuseReflectionCoefficients.push_back( 0.06 );
            diffuseReflectionCoefficients.push_back( 0.46 );
            diffuseReflectionCoefficients.push_back( 0.06 );
            diffuseReflectionCoefficients.push_back( 0.46 );
            diffuseReflectionCoefficients.push_back( 0.06 );
            diffuseReflectionCoefficients.push_back( 0.46 );

            // Panels surface normals
            panelSurfaceNormals.push_back( Eigen::Vector3d::UnitZ( ) );
            panelSurfaceNormals.push_back( Eigen::Vector3d::UnitZ( ) );
            panelSurfaceNormals.push_back( - Eigen::Vector3d::UnitZ( ) );
            panelSurfaceNormals.push_back( Eigen::Vector3d::UnitX( ) );
            panelSurfaceNormals.push_back( Eigen::Vector3d::UnitX( ) );
            panelSurfaceNormals.push_back( - Eigen::Vector3d::UnitX( ) );
            panelSurfaceNormals.push_back( - Eigen::Vector3d::UnitX( ) );
            panelSurfaceNormals.push_back( Eigen::Vector3d::UnitY( ) );
            panelSurfaceNormals.push_back( Eigen::Vector3d::UnitY( ) );
            panelSurfaceNormals.push_back( - Eigen::Vector3d::UnitY( ) );
            panelSurfaceNormals.push_back( - Eigen::Vector3d::UnitY( ) );
        }


        // Define panels properties for test cases with constant rotational model (simpler boxes-and-wings model).
        else if ( testCase == 3 ){

            // Panels areas
            areas.push_back( 9.9 );
            areas.push_back( 2.3 );
            areas.push_back( 9.9 );
            areas.push_back( 2.3 );

            // Emissivities
            emissivities.push_back( 0.0 );
            emissivities.push_back( 0.1 );
            emissivities.push_back( 0.0 );
            emissivities.push_back( 0.1 );

            //Diffuse reflection coefficients
            diffuseReflectionCoefficients.push_back( 0.06 );
            diffuseReflectionCoefficients.push_back( 0.46 );
            diffuseReflectionCoefficients.push_back( 0.06 );
            diffuseReflectionCoefficients.push_back( 0.46 );

            // Panels surface normals
            panelSurfaceNormals.push_back( Eigen::Vector3d::UnitX( ) );
            panelSurfaceNormals.push_back( Eigen::Vector3d::UnitX( ) );
            panelSurfaceNormals.push_back( - Eigen::Vector3d::UnitX( ) );
            panelSurfaceNormals.push_back( - Eigen::Vector3d::UnitX( ) );
        }



        // Create panelled radiation pressure interface.
        std::shared_ptr< PanelledRadiationPressureInterfaceSettings > radiationPressureInterfaceSettings =
                std::make_shared< PanelledRadiationPressureInterfaceSettings >(
                    "Sun", emissivities, areas, diffuseReflectionCoefficients, panelSurfaceNormals );

        std::shared_ptr< PanelledRadiationPressureInterface > radiationPressureInterface =
                std::dynamic_pointer_cast< PanelledRadiationPressureInterface >(
                    createRadiationPressureInterface( radiationPressureInterfaceSettings, "Vehicle", bodies ) );

        bodies.at( "Vehicle" )->setRadiationPressureInterface( "Sun", radiationPressureInterface );



        // Define accelerations.
        SelectedAccelerationMap accelerationMap;
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOnVehicle;
        accelerationsOnVehicle[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                       panelled_radiation_pressure_acceleration ) );

        accelerationMap[ "Vehicle" ] = accelerationsOnVehicle;

        std::map< std::string, std::string > centralBodyMap;
        centralBodyMap[ "Vehicle" ] = "Sun";
        AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodies, accelerationMap, centralBodyMap );
        std::shared_ptr< AccelerationModel< Eigen::Vector3d > > accelerationModel =
                accelerationModelMap.at( "Vehicle" ).at( "Sun" ).at( 0 );


        // Compute spacecraft orbital period, and compute test times
        double orbitalPeriod = 2.0 * mathematical_constants::PI * std::sqrt( std::pow( physical_constants::ASTRONOMICAL_UNIT, 3.0 ) /
                                                                             spice_interface::getBodyGravitationalParameter( "Sun" ) );
        std::vector< double > testTimes = { 0.0, orbitalPeriod / 4.0, orbitalPeriod / 2.0, 3.0 * orbitalPeriod / 4.0 };


        // Compute panelled radiation pressure for various relative Sun positions.
        Eigen::Vector3d calculatedAcceleration, expectedAcceleration;

        Eigen::Vector3d sunCenteredVehiclePosition;
        std::shared_ptr< Ephemeris > vehicleEphemeris = bodies.at( "Vehicle" )->getEphemeris( );


        for( unsigned int i = 0; i < testTimes.size( ) ; i++ )
        {
            // Update environment and acceleration
            bodies.at( "Sun" )->setStateFromEphemeris( testTimes[ i ] );
            bodies.at( "Vehicle" )->setStateFromEphemeris( testTimes[ i ] );
            bodies.at( "Vehicle" )->setCurrentRotationToLocalFrameFromEphemeris( testTimes[ i ] );
            bodies.at( "Vehicle" )->updateMass( testTimes[ i ] );
            radiationPressureInterface->updateInterface( testTimes[ i ] );
            accelerationModel->updateMembers( testTimes[ i ] );

            // Retrieve acceleration.
            calculatedAcceleration = accelerationModel->getAcceleration( );

            double radiationPressure = radiationPressureInterface->getCurrentRadiationPressure( );
            Eigen::Vector3d expectedVehicleToSunNormalisedVector = ( bodies.at( "Sun" )->getState( ) - bodies.at( "Vehicle" )->getState( ) ).segment( 0, 3 )
                    .normalized();


            // Calculated accelerations.

            // Calculated acceleration if Sun-Vehicle vector aligned with +X-axis (acceleration generated by panels whose normals are
            // along +X-axis only)
            double cosinusPanelInclinationPositiveXaxis = expectedVehicleToSunNormalisedVector.dot( Eigen::Vector3d::UnitX( ) );
            Eigen::Vector3d accelerationPositiveXaxisOrientedPanels = - cosinusPanelInclinationPositiveXaxis
                    * radiationPressure / bodies.at( "Vehicle" )->getBodyMass()
                    * ( areas[5] * ( ( 1.0 - emissivities[5] ) * expectedVehicleToSunNormalisedVector
                    + 2.0 / 3.0 * diffuseReflectionCoefficients[5] * Eigen::Vector3d::UnitX( ) )
                    + areas[6] * ( ( 1.0 - emissivities[6] ) * expectedVehicleToSunNormalisedVector
                    + ( 2.0 / 3.0 * diffuseReflectionCoefficients[6] + 2.0 * cosinusPanelInclinationPositiveXaxis * emissivities[6] )
                    * Eigen::Vector3d::UnitX( ) )  );


            // Calculated acceleration generated by the panels whose normals are along -X-axis.
            double cosinusPanelInclinationNegativeXaxis = expectedVehicleToSunNormalisedVector.dot( - Eigen::Vector3d::UnitX() );
            Eigen::Vector3d accelerationNegativeXaxisOrientedPanels = - cosinusPanelInclinationNegativeXaxis
                    * radiationPressure / bodies.at( "Vehicle" )->getBodyMass()
                    * ( areas[3] * ( ( 1.0 - emissivities[3] ) * expectedVehicleToSunNormalisedVector
                    + 2.0 / 3.0 * diffuseReflectionCoefficients[3] * - Eigen::Vector3d::UnitX( ) )
                    + areas[4] * ( ( 1.0 - emissivities[4] ) * expectedVehicleToSunNormalisedVector
                    + ( 2.0 / 3.0 * diffuseReflectionCoefficients[4] + 2.0 * emissivities[4] * cosinusPanelInclinationNegativeXaxis )
                    * - Eigen::Vector3d::UnitX( ) ) );

            // Calculated acceleration generated by the panels whose normals are along +Y-axis.
            double cosinusPanelInclinationPositiveYaxis = expectedVehicleToSunNormalisedVector.dot( Eigen::Vector3d::UnitY() );
            Eigen::Vector3d accelerationPositiveYaxisOrientedPanels = - cosinusPanelInclinationPositiveYaxis
                    * radiationPressure / bodies.at( "Vehicle" )->getBodyMass()
                    * ( areas[7] * ( ( 1.0 - emissivities[7] ) * expectedVehicleToSunNormalisedVector
                    + 2.0 / 3.0 * diffuseReflectionCoefficients[7] * Eigen::Vector3d::UnitY( ) )
                    + areas[8] * ( ( 1.0 - emissivities[8] ) * expectedVehicleToSunNormalisedVector
                    + ( 2.0 / 3.0 * diffuseReflectionCoefficients[8] + 2.0 * emissivities[8] * cosinusPanelInclinationPositiveYaxis )
                    * Eigen::Vector3d::UnitY( ) ) );

            // Calculated acceleration generated by the panels whose normals are along -Y-axis.
            double cosinusPanelInclinationNegativeYaxis = expectedVehicleToSunNormalisedVector.dot( - Eigen::Vector3d::UnitY() );
            Eigen::Vector3d accelerationNegativeYaxisOrientedPanels = - cosinusPanelInclinationNegativeYaxis
                    * radiationPressure / bodies.at( "Vehicle" )->getBodyMass()
                    * ( areas[9] * ( ( 1.0 - emissivities[9] ) * expectedVehicleToSunNormalisedVector
                    + 2.0 / 3.0 * diffuseReflectionCoefficients[9] * - Eigen::Vector3d::UnitY( ) )
                    + areas[10] * ( ( 1.0 - emissivities[10] ) * expectedVehicleToSunNormalisedVector
                    + ( 2.0 / 3.0 * diffuseReflectionCoefficients[10] + 2.0 * emissivities[10] * cosinusPanelInclinationNegativeYaxis )
                    * - Eigen::Vector3d::UnitY( ) ) );

            // Calculated acceleration generated by the panels whose normals are along +Z-axis.
            double cosinusPanelInclinationPositiveZaxis = expectedVehicleToSunNormalisedVector.dot( Eigen::Vector3d::UnitZ() );
            Eigen::Vector3d accelerationPositiveZaxisOrientedPanels = - cosinusPanelInclinationPositiveZaxis
                    * radiationPressure / bodies.at( "Vehicle" )->getBodyMass()
                    * ( areas[0] * ( ( 1.0 - emissivities[0] ) * expectedVehicleToSunNormalisedVector
                    + 2.0 / 3.0 * diffuseReflectionCoefficients[0] * Eigen::Vector3d::UnitZ( ) )
                    + areas[1] * ( ( 1.0 - emissivities[1] ) * expectedVehicleToSunNormalisedVector
                    + ( 2.0 / 3.0 * diffuseReflectionCoefficients[1] + 2.0 * cosinusPanelInclinationPositiveZaxis * emissivities[1] )
                    * Eigen::Vector3d::UnitZ( ) ) );

            // Calculated acceleration generated by the panels whose normals are along -Z-axis.
            double cosinusPanelInclinationNegativeZaxis = expectedVehicleToSunNormalisedVector.dot( - Eigen::Vector3d::UnitZ() );
            Eigen::Vector3d accelerationNegativeZaxisOrientedPanels = - cosinusPanelInclinationNegativeZaxis
                    * radiationPressure / bodies.at( "Vehicle" )->getBodyMass()
                    * areas[2] * ( ( 1.0 - emissivities[2] ) * expectedVehicleToSunNormalisedVector
                    + 2.0 / 3.0 * diffuseReflectionCoefficients[2] * - Eigen::Vector3d::UnitZ( ) );




            // Equatorial orbit case: vehicle in circular orbit around the Sun with i = 0 deg.
            if ( testCase == 0 ){

                // At t = 0: Sun-Vehicle vector expected along +X-axis
                // The acceleration is expected to be generated by the panels whose normals are along -X-axis only.
                if ( i == 0 ){
                    expectedAcceleration = accelerationNegativeXaxisOrientedPanels;
                }

                // At t = period/4: Sun-Vehicle vector expected along +Y-axis
                // The acceleration is expected to be generated by the panels whose normals are along -Y-axis only.
                else if ( i == 1 ){
                    expectedAcceleration = accelerationNegativeYaxisOrientedPanels;
                }

                // At t = period/2: Sun-Vehicle vector expected along -X-axis
                // The acceleration is expected to be generated by the panels whose normals are along +X-axis only.
                else if ( i == 2 ){
                    expectedAcceleration = accelerationPositiveXaxisOrientedPanels;
                }

                // At t = 3 * period/4: Sun-Vehicle vector expected along -Y-axis
                // The acceleration is expected to be generated by the panels whose normals are along +Y-axis only.
                else if ( i == 3 ){
                    expectedAcceleration = accelerationPositiveYaxisOrientedPanels;
                }
            }


            // Polar orbit case: vehicle in circular orbit around the Sun with i = 90 deg
            if ( testCase == 1 ){

                // At t = 0: Sun-Vehicle vector expected along +X-axis.
                // The acceleration is expected to be generated by the panels whose normals are along -X-axis only.
                if ( i == 0 ){
                    expectedAcceleration = accelerationNegativeXaxisOrientedPanels;
                }

                // At t = period/4: Sun-Vehicle vector expected along +Z-axis.
                // The acceleration is expected to be generated by the panels whose normals are along -Z-axis only.
                else if ( i == 1 ){
                    expectedAcceleration = accelerationNegativeZaxisOrientedPanels;
                }

                // At t = period/2: Sun-Vehicle vector expected along -X-axis.
                // The acceleration is expected to be generated by the panels whose normals are along +X-axis only.
                else if ( i == 2 ){
                    expectedAcceleration = accelerationPositiveXaxisOrientedPanels;
                }

                // At t = 3*period/4: Sun-Vehicle vector expected along -Z-axis.
                // The acceleration is expected to be generated by the panels whose normals are along +Z-axis only.
                else if ( i == 3 ){
                    expectedAcceleration = accelerationPositiveZaxisOrientedPanels;
                }
            }


            // Circular orbit with arbitrary chosen inclination
            if ( testCase == 2 ){

                // At t = 0: Sun-Vehicle vector expected along +X-axis.
                // The acceleration is expected to be generated by the panels whose normals are along -X-axis only.
                if ( i == 0 ){
                    expectedAcceleration = accelerationNegativeXaxisOrientedPanels;
                }

                // At t = period/4: Sun-vehicle vector expected to have components along +Z and +Y axes.
                // The acceleration is expected to be generated by the panels whose normals are along -Y and -Z axes.
                else if ( i == 1 ){
                    expectedAcceleration = accelerationNegativeZaxisOrientedPanels + accelerationNegativeYaxisOrientedPanels;
                }

                // At t = period/2: Sun-Vehicle vector expected along -X-axis.
                // The acceleration is expected to be generated by the panels whose normals are along +X-axis only.
                else if ( i == 2 ){
                    expectedAcceleration = accelerationPositiveXaxisOrientedPanels;
                }

                // At t = period/4: Sun-vehicle vector expected to have components along -Z and -Y axes.
                // The acceleration is expected to be generated by the panels whose normals are along +Y and +Z axes.
                else if ( i == 3 ){
                    expectedAcceleration = accelerationPositiveZaxisOrientedPanels + accelerationPositiveYaxisOrientedPanels;
                }

            }


            // Case with non constant rotational model for the spacecraft.
            if ( testCase == 3 )
            {
               Eigen::Quaterniond currentRotationToInertialFrame = bodies.at( "Vehicle" )->getRotationalEphemeris( )
                        ->getRotationToBaseFrame( testTimes.at( i ) );
               Eigen::Vector3d expectedPanelNormalPositiveXaxis = currentRotationToInertialFrame * Eigen::Vector3d::UnitX();
               Eigen::Vector3d expectedPanelNormalNegativeXaxis = currentRotationToInertialFrame * - Eigen::Vector3d::UnitX();


               double cosinusPanelInclinationPositiveXaxis = expectedVehicleToSunNormalisedVector.dot( expectedPanelNormalPositiveXaxis );

               // Determine the normal of the panels that are actually generating a radiation pressure acceleration,
               // depending on their current inertial orientation.
               Eigen::Vector3d expectedPanelNormal;
               if ( cosinusPanelInclinationPositiveXaxis >= 0.0 ){
                   expectedPanelNormal = expectedPanelNormalPositiveXaxis;
               }
               else{
                   expectedPanelNormal = expectedPanelNormalNegativeXaxis;
               }

               // Calculate expected acceleration (same characteristics for the panels oriented along the -X and +X axes).
               expectedAcceleration = - radiationPressure / bodies.at( "Vehicle" )->getBodyMass()
                       * ( areas[0] * fabs( cosinusPanelInclinationPositiveXaxis ) * ( ( 1.0 - emissivities[0] ) * expectedVehicleToSunNormalisedVector
                       + 2.0 / 3.0 * diffuseReflectionCoefficients[0] * expectedPanelNormal )
                       + areas[1] * fabs( cosinusPanelInclinationPositiveXaxis ) * ( ( 1.0 - emissivities[1] ) * expectedVehicleToSunNormalisedVector
                       + ( 2.0 / 3.0 * diffuseReflectionCoefficients[1] + 2.0 * fabs(cosinusPanelInclinationPositiveXaxis) * emissivities[1] )
                       * expectedPanelNormal ) );
            }


            for( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL( std::fabs( calculatedAcceleration( j ) - expectedAcceleration( j ) ), 3.0E-23 );
            }

        }
    }

}


//! Test panelled radiation acceleration model with time-varying solar panels orientation.
BOOST_AUTO_TEST_CASE( testPanelledRadiationPressureTimeVaryingPanelOrientation )
{

    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Create bodies needed in simulation
    double initialEphemerisTime = 0.0;
    double finalEphemerisTime = 1.1 * 365.25 * 86400.0;
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Sun" );
    SystemOfBodies bodies = createSystemOfBodies(
                getDefaultBodySettings( bodyNames,initialEphemerisTime, finalEphemerisTime ) );

    // Create vehicle
    double vehicleMass = 2000.0;
    bodies.createEmptyBody( "Vehicle" );
    bodies.at( "Vehicle" )->setConstantBodyMass( vehicleMass );


    // Compute spacecraft orbital period, and compute test times
    double orbitalPeriod = 2.0 * mathematical_constants::PI * std::sqrt( std::pow( physical_constants::ASTRONOMICAL_UNIT, 3.0 ) /
                                                                         spice_interface::getBodyGravitationalParameter( "Sun" ) );
    std::vector< double > testTimes = { 0.0, orbitalPeriod / 4.0, orbitalPeriod / 2.0, 3.0 * orbitalPeriod / 4.0 };

    // Put vehicle on circular orbit around Sun.
    Eigen::Vector6d initialStateInKeplerianElements = Eigen::Vector6d::Zero( );
    initialStateInKeplerianElements[ 0 ] = physical_constants::ASTRONOMICAL_UNIT;
    bodies.at( "Vehicle" )->setEphemeris( std::make_shared< KeplerEphemeris >( initialStateInKeplerianElements, 0.0,
                                                 spice_interface::getBodyGravitationalParameter( "Sun" ), "Sun", "ECLIPJ2000", 1 ) );


    // Create radiation pressure properties.
    std::vector< double > areas;
    areas.push_back( 2.0 );
    areas.push_back( 4.0 );

    std::vector< double > emissivities;
    emissivities.push_back( 0.0 );
    emissivities.push_back( 0.1 );

    std::vector< double > diffuseReflectionCoefficients;
    diffuseReflectionCoefficients.push_back( 0.06 );
    diffuseReflectionCoefficients.push_back( 0.46 );

    // Define (constant) panel surface orientations.
    std::vector< Eigen::Vector3d > panelSurfaceNormals;
    panelSurfaceNormals.push_back( Eigen::Vector3d::UnitX( ) );
    panelSurfaceNormals.push_back( - Eigen::Vector3d::UnitX( ) );



    // Define rotational model parameters.
    std::vector< double > rightAscensionPole;
    rightAscensionPole.push_back( 0.0 );
    rightAscensionPole.push_back( 0.2 );

    std::vector< double > declinationPole;
    declinationPole.push_back( mathematical_constants::PI / 2.0 );
    declinationPole.push_back( 0.4 );

    std::vector< double > primeMeridianLongitude;
    primeMeridianLongitude.push_back( - mathematical_constants::PI / 2.0 );
    primeMeridianLongitude.push_back( - 0.2 );

    std::vector< double > rotationalRate;
    rotationalRate.push_back( 1.0E-5 );
    rotationalRate.push_back( 1.0E-5 );

    std::vector< double > numberSecondsSinceEpoch;
    numberSecondsSinceEpoch.push_back( 0.0 );
    numberSecondsSinceEpoch.push_back( 0.0 );


    // Case 0: vehicle-fixed axes aligned with inertial ones at time t = 0.
    // Case 1: arbitrary chosen rotational model for the vehicle.

    for ( unsigned int testCase = 0 ; testCase < 2 ; testCase++){


        /// First calculation with simple rotational ephemeris and constant panel orientation.

        // Define simple rotational ephemeris.
        bodies.at( "Vehicle" )->setRotationalEphemeris( std::make_shared< tudat::ephemerides::SimpleRotationalEphemeris >(
                        rightAscensionPole[ testCase ], declinationPole[ testCase ],  primeMeridianLongitude[ testCase ],
                        rotationalRate[ testCase ], numberSecondsSinceEpoch[ testCase ], "ECLIPJ2000", "VehicleFixed" ) );




        // Create panelled radiation pressure interface.
        std::shared_ptr< PanelledRadiationPressureInterfaceSettings > radiationPressureInterfaceSettings =
                std::make_shared< PanelledRadiationPressureInterfaceSettings >(
                    "Sun", emissivities, areas, diffuseReflectionCoefficients, panelSurfaceNormals );

        std::shared_ptr< PanelledRadiationPressureInterface > radiationPressureInterface =
                std::dynamic_pointer_cast< PanelledRadiationPressureInterface >(
                    createRadiationPressureInterface( radiationPressureInterfaceSettings, "Vehicle", bodies ) );

        bodies.at( "Vehicle" )->setRadiationPressureInterface( "Sun", radiationPressureInterface );


        // Define accelerations.
        SelectedAccelerationMap accelerationMap;
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOnVehicle;
        accelerationsOnVehicle[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                       panelled_radiation_pressure_acceleration ) );

        accelerationMap[ "Vehicle" ] = accelerationsOnVehicle;

        std::map< std::string, std::string > centralBodyMap;
        centralBodyMap[ "Vehicle" ] = "Sun";
        AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodies, accelerationMap, centralBodyMap );
        std::shared_ptr< AccelerationModel< Eigen::Vector3d > > accelerationModel =
                accelerationModelMap.at( "Vehicle" ).at( "Sun" ).at( 0 );


        // Compute radiation pressure acceleration for different Sun positions.
        std::vector< Eigen::Vector3d > calculatedAcceleration;

        Eigen::Vector3d sunCenteredVehiclePosition;
        std::shared_ptr< Ephemeris > vehicleEphemeris = bodies.at( "Vehicle" )->getEphemeris( );


        for( unsigned int i = 0; i < testTimes.size( ) ; i++ )
        {
            // Update environment and acceleration
            bodies.at( "Sun" )->setStateFromEphemeris( testTimes[ i ] );
            bodies.at( "Vehicle" )->setStateFromEphemeris( testTimes[ i ] );
            bodies.at( "Vehicle" )->setCurrentRotationToLocalFrameFromEphemeris( testTimes[ i ] );
            bodies.at( "Vehicle" )->updateMass( testTimes[ i ] );
            radiationPressureInterface->updateInterface( testTimes[ i ] );
            accelerationModel->updateMembers( testTimes[ i ] );

            // Retrieve acceleration.
            calculatedAcceleration.push_back( accelerationModel->getAcceleration( ) );

        }



        /// Second calculation with constant rotational ephemeris and time-varying panel orientation

        // Define constant rotational ephemeris
        Eigen::Vector7d rotationalStateVehicle;
        rotationalStateVehicle.segment( 0, 4 ) = linear_algebra::convertQuaternionToVectorFormat( Eigen::Quaterniond( Eigen::Matrix3d::Identity() ));
        rotationalStateVehicle.segment( 4, 3 ) = Eigen::Vector3d::Zero();
        bodies.at( "Vehicle" )->setRotationalEphemeris( std::make_shared< ConstantRotationalEphemeris >( rotationalStateVehicle, "ECLIPJ2000",
                                                                                                       "VehicleFixed" ) );

        // Define time-varying panel orientation.
        std::vector< std::function< Eigen::Vector3d ( const double ) > > timeVaryingPanelSurfaceNormals;
        for ( unsigned int j = 0 ; j < panelSurfaceNormals.size() ; j++ )
        {
            timeVaryingPanelSurfaceNormals.push_back( [ = ]( const double currentTime ){

                Eigen::Vector3d currentPanelSurfaceOrientation = (
                            reference_frames::getInertialToPlanetocentricFrameTransformationQuaternion(
                                basic_mathematics::computeModulo( ( currentTime - numberSecondsSinceEpoch[ testCase ] )
                                                                  * rotationalRate[ testCase ], 2.0 * mathematical_constants::PI ) )
                          * reference_frames::getInertialToPlanetocentricFrameTransformationQuaternion( declinationPole[ testCase ],
                                                                                                        rightAscensionPole[ testCase ],
                                                                                                        primeMeridianLongitude[ testCase ] ) )
                          .toRotationMatrix().inverse() * panelSurfaceNormals[ j ];

                return currentPanelSurfaceOrientation; } );
        }

        // Create panelled radiation pressure interface.
        std::shared_ptr< PanelledRadiationPressureInterfaceSettings > radiationPressureInterfaceSettingsWithTimeVaryingPanelSurfaceNormal =
                std::make_shared< PanelledRadiationPressureInterfaceSettings >(
                    "Sun", emissivities, areas, diffuseReflectionCoefficients, timeVaryingPanelSurfaceNormals );

        std::shared_ptr< PanelledRadiationPressureInterface > radiationPressureInterfaceTimeVaryingSurfaceNormal =
                std::dynamic_pointer_cast< PanelledRadiationPressureInterface >(
                    createRadiationPressureInterface( radiationPressureInterfaceSettingsWithTimeVaryingPanelSurfaceNormal, "Vehicle", bodies ) );

        bodies.at( "Vehicle" )->setRadiationPressureInterface( "Sun", radiationPressureInterfaceTimeVaryingSurfaceNormal );

        accelerationModelMap = createAccelerationModelsMap( bodies, accelerationMap, centralBodyMap );

        std::shared_ptr< AccelerationModel< Eigen::Vector3d > > accelerationModelTimeVaryingPanelSurfaceNormal =
                accelerationModelMap.at( "Vehicle" ).at( "Sun" ).at( 0 );


        // Compute radiation pressure acceleration for different Sun positions.
        std::vector< Eigen::Vector3d > calculatedAccelerationTimeVaryingPanelOrientation;

        for( unsigned int i = 0; i < testTimes.size( ) ; i++ )
        {
            // Update environment and acceleration
            bodies.at( "Sun" )->setStateFromEphemeris( testTimes[ i ] );
            bodies.at( "Vehicle" )->setStateFromEphemeris( testTimes[ i ] );
            bodies.at( "Vehicle" )->setCurrentRotationToLocalFrameFromEphemeris( testTimes[ i ] );
            bodies.at( "Vehicle" )->updateMass( testTimes[ i ] );
            radiationPressureInterfaceTimeVaryingSurfaceNormal->updateInterface( testTimes[ i ] );
            accelerationModelTimeVaryingPanelSurfaceNormal->updateMembers( testTimes[ i ] );

            // Retrieve acceleration.
            calculatedAccelerationTimeVaryingPanelOrientation.push_back( accelerationModelTimeVaryingPanelSurfaceNormal->getAcceleration( ) );
        }


        for( unsigned int j = 0; j < testTimes.size() ; j++ )
        {
            for ( unsigned int i = 0 ; i < 3 ; i++ ){
                BOOST_CHECK_SMALL( std::fabs(
                                       calculatedAcceleration[ j ][ i ] - calculatedAccelerationTimeVaryingPanelOrientation[ j ][ i ] ), 3.0E-23 );
            }
        }
    }
}




BOOST_AUTO_TEST_SUITE_END( )


}

}
