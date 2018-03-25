#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Basics/testMacros.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/External/SpiceInterface/spiceInterface.h"

#include "Tudat/Astrodynamics/ElectroMagnetism/panelledRadiationPressure.h"
#include "Tudat/Astrodynamics/Ephemerides/keplerEphemeris.h"
#include "Tudat/Astrodynamics/Ephemerides/constantRotationalEphemeris.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createBodies.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.h"
#include "Tudat/SimulationSetup/PropagationSetup/createAccelerationModels.h"

namespace tudat
{

using namespace tudat::basic_astrodynamics;
using namespace tudat::simulation_setup;
using namespace tudat::ephemerides;
using namespace tudat::basic_astrodynamics;
using namespace tudat::electro_magnetism;

namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_PanelledRadiationPressure )

BOOST_AUTO_TEST_CASE( testSimpleGeometryRadiationPressure )
{
    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Sun" );

    // Specify initial time
    double initialEphemerisTime = -86400.0;
    double finalEphemerisTime = 1.1 * 365.25 * 86400.0;

    // Create bodies needed in simulation
    NamedBodyMap bodyMap = createBodies(
                getDefaultBodySettings( bodyNames,initialEphemerisTime, finalEphemerisTime ) );

    double vehicleMass = 2500.0;
    bodyMap[ "Vehicle" ] = boost::make_shared< Body >( );
    bodyMap[ "Vehicle" ]->setConstantBodyMass( vehicleMass );

    Eigen::Vector6d initialStateInKeplerianElements = Eigen::Vector6d::Zero( );
    initialStateInKeplerianElements[ 0 ] = physical_constants::ASTRONOMICAL_UNIT;
    bodyMap[ "Vehicle" ]->setEphemeris(
                boost::make_shared< KeplerEphemeris >(
                    initialStateInKeplerianElements, 0.0, spice_interface::getBodyGravitationalParameter( "Sun" ), "Sun", "ECLIPJ2000", 1 ) );
    bodyMap[ "Vehicle" ]->setRotationalEphemeris( boost::make_shared< ConstantRotationalEphemeris >(
                                                      Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ), "ECLIPJ2000", "VehicleFixed" ) );

    std::vector< double > areas;
    areas.push_back( 2.532 );
    areas.push_back( 3.254 );
    areas.push_back( 8.654 );
    areas.push_back( 1.346 );

    std::vector< double > emissivities;
    emissivities.push_back( 0.0 );
    emissivities.push_back( 1.0 );
    emissivities.push_back( 0.0 );
    emissivities.push_back( 1.0 );

    std::vector< double > diffuseReflectionCoefficients;
    diffuseReflectionCoefficients.push_back( 0.0 );
    diffuseReflectionCoefficients.push_back( 0.0 );
    diffuseReflectionCoefficients.push_back( 0.0 );
    diffuseReflectionCoefficients.push_back( 0.0 );


    std::vector< Eigen::Vector3d > panelSurfaceNormals;
    panelSurfaceNormals.push_back( -Eigen::Vector3d::UnitX( ) );
    panelSurfaceNormals.push_back( -Eigen::Vector3d::UnitY( ) );
    panelSurfaceNormals.push_back( Eigen::Vector3d::UnitX( ) );
    panelSurfaceNormals.push_back( Eigen::Vector3d::UnitY( ) );

    boost::shared_ptr< PanelledRadiationPressureInterfaceSettings > radiationPressureInterfaceSettings =
            boost::make_shared< PanelledRadiationPressureInterfaceSettings >(
                "Sun", areas, emissivities, diffuseReflectionCoefficients, panelSurfaceNormals );
    boost::shared_ptr< PanelledRadiationPressureInterface > radiationPressureInterface =
            boost::dynamic_pointer_cast< PanelledRadiationPressureInterface >(
                createRadiationPressureInterface( radiationPressureInterfaceSettings, "Vehicle", bodyMap ) );

    bodyMap[ "Vehicle" ]->setRadiationPressureInterface( "Sun", radiationPressureInterface );

    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    SelectedAccelerationMap accelerationMap;

    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOnVehicle;
    accelerationsOnVehicle[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >(
                                                   panelled_radiation_pressure_acceleration ) );
    accelerationMap[ "Vehicle" ] = accelerationsOnVehicle;

    std::map< std::string, std::string > centralBodyMap;
    centralBodyMap[ "Vehicle" ] = "Sun";

    AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodyMap, accelerationMap, centralBodyMap );
    boost::shared_ptr< AccelerationModel< Eigen::Vector3d > > accelerationModel =
            accelerationModelMap.at( "Vehicle" ).at( "Sun" ).at( 0 );

    boost::shared_ptr< Ephemeris > vehicleEphemeris = bodyMap[ "Vehicle" ]->getEphemeris( );

    double orbitalPeriod = 2.0 * mathematical_constants::PI * std::sqrt( std::pow( physical_constants::ASTRONOMICAL_UNIT, 3.0 ) /
                                                              spice_interface::getBodyGravitationalParameter( "Sun" ) );

    std::vector< double > testTimes = boost::assign::list_of( 0.0 )( orbitalPeriod / 4.0 )( orbitalPeriod / 2.0 )( 3.0 * orbitalPeriod / 4.0 );

    Eigen::Vector3d calculatedAcceleration, expectedAcceleration;
    Eigen::Vector3d sunCenteredVehiclePosition;

    for( unsigned int i = 0; i < testTimes.size( ); i++ )
    {
        bodyMap[ "Sun" ]->setStateFromEphemeris( testTimes[ i ] );
        bodyMap[ "Vehicle" ]->setStateFromEphemeris( testTimes[ i ] );
        bodyMap[ "Vehicle" ]->setCurrentRotationToLocalFrameFromEphemeris( testTimes[ i ] );

        radiationPressureInterface->updateInterface( 0.0 );
        accelerationModel->updateMembers( 0.0 );

        calculatedAcceleration = accelerationModel->getAcceleration( );

        double radiationPressure = radiationPressureInterface->getCurrentRadiationPressure( );
        sunCenteredVehiclePosition = vehicleEphemeris->getCartesianState( testTimes[ i ] ).segment( 0, 3 );

        if( i == 0 || i == 2 )
        {
            expectedAcceleration = radiationPressure / vehicleMass * ( 1.0 - radiationPressureInterface->getEmmisivity( i ) ) *
                    radiationPressureInterface->getArea( i ) * sunCenteredVehiclePosition.normalized( );
        }
        else if( i == 1 || i == 3 )
        {
            expectedAcceleration = -radiationPressure / vehicleMass * (
                        2.0 * radiationPressureInterface->getArea( i ) * radiationPressureInterface->getCurrentSurfaceNormal( i ) );
        }

        std::cout<<calculatedAcceleration.transpose( )<<" "<<expectedAcceleration.transpose( )<<std::endl<<std::endl;

        for( unsigned int j = 0; j < 3; j++ )
        {
            BOOST_CHECK_SMALL( std::fabs( calculatedAcceleration( j ) - expectedAcceleration( j ) ), 2.0E-23 );
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )


}

}
