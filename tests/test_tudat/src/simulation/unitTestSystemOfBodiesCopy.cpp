/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <filesystem>
#include <thread>

#include <boost/test/included/unit_test.hpp>

#include "tudat/astro/aerodynamics/comaModel.h"
#include "tudat/astro/aerodynamics/comaWindModel.h"
#include "tudat/astro/aerodynamics/exponentialAtmosphere.h"
#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/astro/electromagnetism/radiationPressureTargetModel.h"
#include "tudat/astro/electromagnetism/radiationSourceModel.h"
#include "tudat/astro/ephemerides/aeordynamicAngleRotationalEphemeris.h"
#include "tudat/astro/ephemerides/approximatePlanetPositions.h"
#include "tudat/astro/ephemerides/directionBasedRotationalEphemeris.h"
#include "tudat/astro/ephemerides/keplerEphemeris.h"
#include "tudat/astro/ephemerides/multiArcEphemeris.h"
#include "tudat/astro/ephemerides/simpleRotationalEphemeris.h"
#include "tudat/astro/ephemerides/synchronousRotationalEphemeris.h"
#include "tudat/astro/ephemerides/tabulatedRotationalEphemeris.h"
#include "tudat/astro/gravitation/timeDependentSphericalHarmonicsGravityField.h"
#include "tudat/math/integrators/rungeKuttaCoefficients.h"
#include "tudat/math/interpolators/lagrangeInterpolator.h"
#include "tudat/math/interpolators/linearInterpolator.h"
#include "tudat/simulation/environment_setup/copySystemOfBodies.h"
#include "tudat/simulation/environment_setup/createAtmosphereModel.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/createBodyDeformationModel.h"
#include "tudat/simulation/environment_setup/createFlightConditions.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/simulation/environment_setup/createMetric.h"
#include "tudat/simulation/environment_setup/createSystemModel.h"
#include "tudat/simulation/environment_setup/environmentModelClosure.h"
#include "tudat/simulation/environment_setup/thrustSettings.h"
#include "tudat/simulation/propagation_setup/createAccelerationModels.h"
#include "tudat/simulation/propagation_setup/multiArcDynamicsSimulator.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::simulation_setup;
using namespace tudat::ephemerides;
using namespace tudat::gravitation;
using namespace tudat::aerodynamics;
using namespace tudat::electromagnetism;

BOOST_AUTO_TEST_SUITE( test_system_of_bodies_copy )

const double initialTime = 1.0E7;
const double finalTime = initialTime + 86400.0;

// Create tabulated states, from a Kepler orbit, covering the test interval
std::map< double, Eigen::Vector6d > getKeplerOrbitStates( const Eigen::Vector6d& keplerElements, const double gravitationalParameter )
{
    KeplerEphemeris keplerEphemeris( keplerElements, initialTime, gravitationalParameter );
    std::map< double, Eigen::Vector6d > stateHistory;
    for( double time = initialTime - 3600.0; time <= finalTime + 3600.0; time += 120.0 )
    {
        stateHistory[ time ] = keplerEphemeris.getCartesianState( time );
    }
    return stateHistory;
}

// Create tabulated rotational states, from a uniform rotation, covering the test interval
std::map< double, Eigen::Vector7d > getUniformRotationStates( const double rotationRate )
{
    SimpleRotationalEphemeris uniformRotation( Eigen::Quaterniond::Identity( ), rotationRate, initialTime );
    std::map< double, Eigen::Vector7d > rotationalStateHistory;
    for( double time = initialTime - 3600.0; time <= finalTime + 3600.0; time += 300.0 )
    {
        rotationalStateHistory[ time ] = uniformRotation.getRotationStateVector( time );
    }
    return rotationalStateHistory;
}

// Create bodies, with as many different (types of) environment models as possible, none of which use SPICE
BodyListSettings getTestBodySettings( )
{
    BodyListSettings bodySettings( "SSB", "J2000" );
    const double sunGravitationalParameter = 1.32712440018E20;
    const double earthGravitationalParameter = 3.986004418E14;
    const double moonGravitationalParameter = 4.9028E12;

    bodySettings.addSettings( "Sun" );
    bodySettings.at( "Sun" )->ephemerisSettings = constantEphemerisSettings( Eigen::Vector6d::Zero( ), "SSB", "J2000" );
    bodySettings.at( "Sun" )->gravityFieldSettings = centralGravitySettings( sunGravitationalParameter );
    bodySettings.at( "Sun" )->rotationModelSettings =
            simpleRotationModelSettings( "J2000", "IAU_Sun", Eigen::Quaterniond::Identity( ), initialTime, 2.8E-6 );
    bodySettings.at( "Sun" )->radiationSourceModelSettings =
            isotropicPointRadiationSourceModelSettings( constantLuminosityModelSettings( 3.828E26 ) );
    bodySettings.at( "Sun" )->shapeModelSettings = sphericalBodyShapeSettings( 6.957E8 );

    // Earth: Kepler ephemeris, GCRS-ITRS rotation, spherical harmonics with direct and interpolated tides, tabulated atmosphere
    Eigen::MatrixXd earthCosineCoefficients = Eigen::MatrixXd::Zero( 5, 5 );
    Eigen::MatrixXd earthSineCoefficients = Eigen::MatrixXd::Zero( 5, 5 );
    earthCosineCoefficients( 0, 0 ) = 1.0;
    earthCosineCoefficients( 2, 0 ) = -4.84E-4;
    earthCosineCoefficients( 2, 2 ) = 2.4E-6;
    earthSineCoefficients( 2, 2 ) = -1.4E-6;
    earthCosineCoefficients( 3, 1 ) = 2.0E-6;
    earthSineCoefficients( 4, 3 ) = -2.0E-7;
    bodySettings.addSettings( "Earth" );
    bodySettings.at( "Earth" )->ephemerisSettings =
            keplerEphemerisSettings( ( Eigen::Vector6d( ) << 1.496E11, 0.0167, 0.4, 1.0, 2.0, 0.5 ).finished( ),
                                     initialTime,
                                     sunGravitationalParameter,
                                     "Sun",
                                     "J2000" );
    bodySettings.at( "Earth" )->rotationModelSettings = gcrsToItrsRotationModelSettings( basic_astrodynamics::iau_2006, "J2000" );
    bodySettings.at( "Earth" )->gravityFieldSettings = sphericalHarmonicsGravitySettings(
            earthGravitationalParameter, 6378137.0, earthCosineCoefficients, earthSineCoefficients, "ITRS" );
    std::map< int, std::vector< std::complex< double > > > earthLoveNumbers;
    earthLoveNumbers[ 2 ] = { 0.30, 0.29, 0.28 };
    bodySettings.at( "Earth" )->gravityFieldVariationSettings.push_back( std::make_shared< BasicSolidBodyGravityFieldVariationSettings >(
            std::vector< std::string >{ "Moon", "Sun" }, earthLoveNumbers ) );
    bodySettings.at( "Earth" )->gravityFieldVariationSettings.push_back( std::make_shared< BasicSolidBodyGravityFieldVariationSettings >(
            std::vector< std::string >{ "Sun" },
            earthLoveNumbers,
            std::make_shared< ModelInterpolationSettings >( initialTime - 1800.0, finalTime + 1800.0, 600.0 ) ) );
    // Tabulated atmosphere
    bodySettings.at( "Earth" )->atmosphereSettings =
            tabulatedAtmosphereSettings( paths::getAtmosphereTablesPath( ) + "/USSA1976Until100kmPer100mUntil1000kmPer1000m.dat" );
    bodySettings.at( "Earth" )->shapeModelSettings = oblateSphericalBodyShapeSettings( 6378137.0, 1.0 / 298.257 );
    bodySettings.at( "Earth" )->radiationSourceModelSettings = extendedRadiationSourceModelSettings(
            { albedoPanelRadiosityModelSettings( 0.3, "Sun" ), delayedThermalPanelRadiosityModelSettings( 0.95, "Sun" ) },
            { 6, 12 },
            { "Moon" } );

    // Moon: tabulated ephemeris, synchronous rotation, spherical harmonics with polynomial and periodic variations
    Eigen::MatrixXd moonCosineCoefficients = Eigen::MatrixXd::Zero( 4, 4 );
    Eigen::MatrixXd moonSineCoefficients = Eigen::MatrixXd::Zero( 4, 4 );
    moonCosineCoefficients( 0, 0 ) = 1.0;
    moonCosineCoefficients( 2, 0 ) = -9.1E-5;
    moonCosineCoefficients( 2, 2 ) = 3.5E-5;
    moonSineCoefficients( 3, 1 ) = 5.0E-6;
    bodySettings.addSettings( "Moon" );
    bodySettings.at( "Moon" )->ephemerisSettings = tabulatedEphemerisSettings(
            getKeplerOrbitStates( ( Eigen::Vector6d( ) << 3.844E8, 0.0549, 0.09, 0.3, 0.2, 1.0 ).finished( ), earthGravitationalParameter ),
            "Earth",
            "J2000" );
    bodySettings.at( "Moon" )->rotationModelSettings = synchronousRotationModelSettings( "Earth", "J2000", "IAU_Moon" );
    bodySettings.at( "Moon" )->gravityFieldSettings = sphericalHarmonicsGravitySettings(
            moonGravitationalParameter, 1737.4E3, moonCosineCoefficients, moonSineCoefficients, "IAU_Moon" );
    bodySettings.at( "Moon" )->gravityFieldVariationSettings.push_back( polynomialGravityFieldVariationsSettingsSinglePower(
            Eigen::MatrixXd::Constant( 2, 3, 1.0E-12 ), Eigen::MatrixXd::Constant( 2, 3, -2.0E-12 ), 1, initialTime ) );
    bodySettings.at( "Moon" )->gravityFieldVariationSettings.push_back(
            periodicGravityFieldVariationsSettingsSingleFrequency( Eigen::MatrixXd::Constant( 2, 3, 1.0E-9 ),
                                                                   Eigen::MatrixXd::Constant( 2, 3, 2.0E-9 ),
                                                                   Eigen::MatrixXd::Constant( 2, 3, -1.0E-9 ),
                                                                   Eigen::MatrixXd::Constant( 2, 3, 3.0E-9 ),
                                                                   2.0 * mathematical_constants::PI / 86400.0,
                                                                   initialTime ) );
    bodySettings.at( "Moon" )->shapeModelSettings = sphericalBodyShapeSettings( 1737.4E3 );

    // Mars: multi-arc ephemeris, IAU rotation, ring gravity field, scaled exponential atmosphere
    bodySettings.addSettings( "Mars" );
    std::map< double, std::shared_ptr< EphemerisSettings > > marsArcEphemerides;
    marsArcEphemerides[ initialTime - 3600.0 ] = tabulatedEphemerisSettings(
            getKeplerOrbitStates( ( Eigen::Vector6d( ) << 2.28E11, 0.093, 0.03, 0.4, 0.8, 2.0 ).finished( ), sunGravitationalParameter ),
            "Sun",
            "J2000" );
    marsArcEphemerides[ initialTime + 43200.0 ] = tabulatedEphemerisSettings(
            getKeplerOrbitStates( ( Eigen::Vector6d( ) << 2.28E11, 0.093, 0.03, 0.4, 0.8, 2.001 ).finished( ), sunGravitationalParameter ),
            "Sun",
            "J2000" );
    bodySettings.at( "Mars" )->ephemerisSettings = multiArcEphemerisSettings( marsArcEphemerides, "Sun", "J2000" );
    bodySettings.at( "Mars" )->rotationModelSettings = iauRotationModelSettings( "J2000",
                                                                                 "IAU_Mars",
                                                                                 3.4,
                                                                                 ( Eigen::Vector2d( ) << 5.5, 0.9 ).finished( ),
                                                                                 7.088E-5,
                                                                                 ( Eigen::Vector2d( ) << 1.0E-12, -2.0E-12 ).finished( ),
                                                                                 { { 1.0E-6, { 1.0E-5, 0.1 } } },
                                                                                 {} );
    bodySettings.at( "Mars" )->gravityFieldSettings = ringGravitySettings( 4.2828E13, 3.3895E6, "IAU_Mars" );
    bodySettings.at( "Mars" )->atmosphereSettings = scaledAtmosphereSettings(
            exponentialAtmosphereSettings( 11.1E3, 0.02, 210.0 ), []( const double time ) { return 1.0 + 1.0E-8 * time; }, false );
    bodySettings.at( "Mars" )->shapeModelSettings = sphericalBodyShapeSettings( 3389.5E3 );

    // Asteroid: scaled Kepler ephemeris, tabulated rotation, polyhedron gravity field
    Eigen::MatrixXd verticesCoordinates( 8, 3 );
    verticesCoordinates << 0.0, 0.0, 0.0, 2.0E3, 0.0, 0.0, 0.0, 1.0E3, 0.0, 2.0E3, 1.0E3, 0.0, 0.0, 0.0, 1.0E3, 2.0E3, 0.0, 1.0E3, 0.0,
            1.0E3, 1.0E3, 2.0E3, 1.0E3, 1.0E3;
    Eigen::MatrixXi verticesDefiningEachFacet( 12, 3 );
    verticesDefiningEachFacet << 2, 1, 0, 1, 2, 3, 4, 2, 0, 2, 4, 6, 1, 4, 0, 4, 1, 5, 6, 5, 7, 5, 6, 4, 3, 6, 7, 6, 3, 2, 5, 3, 7, 3, 5, 1;
    bodySettings.addSettings( "Asteroid" );
    bodySettings.at( "Asteroid" )->ephemerisSettings =
            scaledEphemerisSettings( keplerEphemerisSettings( ( Eigen::Vector6d( ) << 4.0E11, 0.2, 0.1, 0.4, 1.0, 3.0 ).finished( ),
                                                              initialTime,
                                                              sunGravitationalParameter,
                                                              "Sun",
                                                              "J2000" ),
                                     1.0001,
                                     false );
    bodySettings.at( "Asteroid" )->rotationModelSettings =
            tabulatedRotationModelSettings( getUniformRotationStates( 3.0E-4 ), "J2000", "AsteroidFixed" );
    bodySettings.at( "Asteroid" )->gravityFieldSettings =
            polyhedronGravitySettingsFromMu( 1.0E3, verticesCoordinates, verticesDefiningEachFacet, "AsteroidFixed" );

    // Vehicle: tabulated ephemeris w.r.t. Earth, aerodynamic angle based rotation, tabulated aerodynamic coefficients
    bodySettings.addSettings( "Vehicle" );
    bodySettings.at( "Vehicle" )->ephemerisSettings = tabulatedEphemerisSettings(
            getKeplerOrbitStates( ( Eigen::Vector6d( ) << 6378137.0 + 300.0E3, 0.001, 1.0, 0.3, 0.2, 0.1 ).finished( ),
                                  earthGravitationalParameter ),
            "Earth",
            "J2000" );
    bodySettings.at( "Vehicle" )->rotationModelSettings = aerodynamicAngleRotationSettings( "Earth", "J2000", "VehicleFixed" );
    bodySettings.at( "Vehicle" )->rigidBodyPropertiesSettings = constantRigidBodyPropertiesSettings(
            500.0, ( Eigen::Vector3d( ) << 0.1, 0.0, 0.0 ).finished( ), Eigen::Matrix3d::Identity( ) * 100.0 );
    bodySettings.at( "Vehicle" )->aerodynamicCoefficientSettings =
            oneDimensionalTabulatedAerodynamicCoefficientSettings( { 1.0, 10.0, 20.0, 30.0, 40.0 },
                                                                   { Eigen::Vector3d( 1.0, 0.1, 0.2 ),
                                                                     Eigen::Vector3d( 1.2, 0.1, 0.3 ),
                                                                     Eigen::Vector3d( 1.4, 0.2, 0.3 ),
                                                                     Eigen::Vector3d( 1.5, 0.2, 0.4 ),
                                                                     Eigen::Vector3d( 1.6, 0.3, 0.4 ) },
                                                                   { Eigen::Vector3d( 0.0, 0.01, 0.0 ),
                                                                     Eigen::Vector3d( 0.0, 0.02, 0.0 ),
                                                                     Eigen::Vector3d( 0.0, 0.02, 0.01 ),
                                                                     Eigen::Vector3d( 0.0, 0.03, 0.01 ),
                                                                     Eigen::Vector3d( 0.0, 0.03, 0.02 ) },
                                                                   1.5,
                                                                   2.0,
                                                                   aerodynamics::mach_number_dependent );
    bodySettings.at( "Vehicle" )->radiationPressureTargetModelSettings =
            cannonballRadiationPressureTargetModelSettings( 2.0, 1.3, { "Earth" } );

    // Vehicle2: tabulated ephemeris w.r.t. Moon, orbital-state based rotation, scaled constant aerodynamic coefficients
    bodySettings.addSettings( "Vehicle2" );
    bodySettings.at( "Vehicle2" )->ephemerisSettings = tabulatedEphemerisSettings(
            getKeplerOrbitStates( ( Eigen::Vector6d( ) << 1837.4E3, 0.01, 1.4, 0.3, 0.2, 0.1 ).finished( ), moonGravitationalParameter ),
            "Moon",
            "J2000" );
    bodySettings.at( "Vehicle2" )->rotationModelSettings =
            orbitalStateBasedRotationSettings( "Moon", true, false, "J2000", "Vehicle2Fixed" );
    bodySettings.at( "Vehicle2" )->rigidBodyPropertiesSettings = constantRigidBodyPropertiesSettings( 200.0 );
    bodySettings.at( "Vehicle2" )->aerodynamicCoefficientSettings = scaledAerodynamicCoefficientSettings(
            constantAerodynamicCoefficientSettings( 3.0, Eigen::Vector3d( 2.0, 0.0, 0.1 ) ), 1.1, 1.0, false );
    bodySettings.at( "Vehicle2" )->radiationPressureTargetModelSettings = cannonballRadiationPressureTargetModelSettings( 3.0, 1.1 );

    // Vehicle3: panelled body (with rotating solar array) with panelled radiation pressure target model, and engine
    bodySettings.addSettings( "Vehicle3" );
    bodySettings.at( "Vehicle3" )->ephemerisSettings = tabulatedEphemerisSettings(
            getKeplerOrbitStates( ( Eigen::Vector6d( ) << 6378137.0 + 800.0E3, 0.002, 1.2, 0.4, 0.2, 0.3 ).finished( ),
                                  earthGravitationalParameter ),
            "Earth",
            "J2000" );
    bodySettings.at( "Vehicle3" )->rotationModelSettings =
            simpleRotationModelSettings( "J2000", "Vehicle3Fixed", Eigen::Quaterniond::Identity( ), initialTime, 1.0E-3 );
    bodySettings.at( "Vehicle3" )->rigidBodyPropertiesSettings = constantRigidBodyPropertiesSettings( 1000.0 );
    std::shared_ptr< BodyPanelReflectionLawSettings > reflectionLawSettings = specularDiffuseBodyPanelReflectionLawSettings( 0.3, 0.2 );
    bodySettings.at( "Vehicle3" )->bodyExteriorPanelSettings_ = fullPanelledBodySettings(
            { bodyPanelSettings( frameFixedPanelGeometry( Eigen::Vector3d::UnitX( ), 2.0 ), reflectionLawSettings ),
              bodyPanelSettings( frameFixedPanelGeometry( -Eigen::Vector3d::UnitZ( ), 3.0 ), reflectionLawSettings ),
              bodyPanelSettings( frameFixedPanelGeometry( Eigen::Vector3d::UnitY( ), 10.0, "SolarArray" ), reflectionLawSettings ),
              bodyPanelSettings( bodyTrackingPanelGeometry( "Sun", true, 5.0 ), reflectionLawSettings ) },
            { { "SolarArray", simpleRotationModelSettings( "", "SolarArray", Eigen::Quaterniond::Identity( ), initialTime, 2.0E-3 ) } } );
    bodySettings.at( "Vehicle3" )->radiationPressureTargetModelSettings = paneledRadiationPressureTargetModelSettings( );

    // Moment contribution of aerodynamic force (which requires flight conditions to be closed)
    bodySettings.at( "Vehicle" )->aerodynamicCoefficientSettings->setAddForceContributionToMoments( true );

    return bodySettings;
}

SystemOfBodies createTestBodies( const BodyListSettings& bodySettings = getTestBodySettings( ) )
{
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    addEngineModel( "Vehicle3", "MainEngine", constantThrustMagnitudeSettings( 10.0, 300.0 ), bodies );

    // Flight conditions (after which aerodynamic moment contribution is closed)
    addAtmosphericFlightConditions( bodies, "Vehicle", "Earth" );
    addFlightConditions( bodies, "Vehicle2", "Moon" );
    std::dynamic_pointer_cast< AerodynamicAngleRotationalEphemeris >( bodies.at( "Vehicle" )->getRotationalEphemeris( ) )
            ->setIsBodyInPropagation( true );

    return bodies;
}

// Total irradiance of a radiation source at a given position
double getTotalIrradiance( const std::shared_ptr< electromagnetism::RadiationSourceModel > sourceModel, const Eigen::Vector3d& position )
{
    double totalIrradiance = 0.0;
    for( auto irradiance : sourceModel->evaluateIrradianceAtPosition( position ) )
    {
        totalIrradiance += irradiance.first;
    }
    return totalIrradiance;
}

// Create copy of the test bodies (aerodynamic angle based rotation can only be evaluated during propagation, as set for the original)
SystemOfBodies createTestBodiesCopy( const SystemOfBodies& bodies, const bool printWarnings = true )
{
    SystemOfBodies bodiesCopy = createSystemOfBodiesCopy( bodies, printWarnings );
    bodiesCopy.at( "Vehicle" )->getRotationalEphemeris( )->setIsBodyInPropagation( true );
    return bodiesCopy;
}

const std::vector< std::string > bodyNames = { "Sun", "Earth", "Moon", "Mars", "Asteroid", "Vehicle", "Vehicle2", "Vehicle3" };

void appendVector( std::vector< double >& values, const Eigen::VectorXd& vector )
{
    for( int i = 0; i < vector.rows( ); i++ )
    {
        values.push_back( vector( i ) );
    }
}

// Set the current states and rotations of all bodies from their ephemerides
void setBodyStates( const SystemOfBodies& bodies, const double time )
{
    for( const std::string& bodyName : bodyNames )
    {
        // States set manually by a test (at the same time) are overwritten
        bodies.at( bodyName )->recomputeStateOnNextCall( );
        bodies.at( bodyName )->setStateFromEphemeris< double, double >( time );
    }
    for( const std::string& bodyName : bodyNames )
    {
        if( bodyName != "Vehicle" )
        {
            bodies.at( bodyName )->setCurrentRotationalStateToLocalFrameFromEphemeris( time );
        }
    }
}

// Evaluate all environment models of the test bodies (using the current states of the bodies, set by setBodyStates)
std::vector< double > evaluateEnvironment( const SystemOfBodies& bodies, const double time )
{
    std::vector< double > values;
    const Eigen::Vector3d bodyFixedPosition( 7.0E6, 2.0E6, -1.0E6 );

    for( const std::string& bodyName : bodyNames )
    {
        appendVector( values, bodies.at( bodyName )->getEphemeris( )->getCartesianState( time ) );
        appendVector( values, bodies.at( bodyName )->getState( ) );
        if( bodyName != "Vehicle" )
        {
            appendVector( values, bodies.at( bodyName )->getRotationalEphemeris( )->getRotationToBaseFrame( time ).coeffs( ) );
            appendVector( values, bodies.at( bodyName )->getCurrentRotationToGlobalFrame( ).coeffs( ) );
        }

        std::shared_ptr< GravityFieldModel > gravityField = bodies.at( bodyName )->getGravityFieldModel( );
        if( gravityField != nullptr )
        {
            std::shared_ptr< TimeDependentSphericalHarmonicsGravityField > timeDependentField =
                    std::dynamic_pointer_cast< TimeDependentSphericalHarmonicsGravityField >( gravityField );
            if( timeDependentField != nullptr )
            {
                timeDependentField->update( time );
                appendVector( values,
                              Eigen::Map< const Eigen::VectorXd >( timeDependentField->getCosineCoefficients( ).data( ),
                                                                   timeDependentField->getCosineCoefficients( ).size( ) ) );
                appendVector( values,
                              Eigen::Map< const Eigen::VectorXd >( timeDependentField->getSineCoefficients( ).data( ),
                                                                   timeDependentField->getSineCoefficients( ).size( ) ) );
            }
            appendVector( values, gravityField->getGradientOfPotential( bodyFixedPosition ) );
            values.push_back( gravityField->getGravitationalPotential( bodyFixedPosition ) );
        }
        if( bodies.at( bodyName )->getAtmosphereModel( ) != nullptr )
        {
            values.push_back( bodies.at( bodyName )->getAtmosphereModel( )->getDensity( 300.0E3, 0.3, 0.2, time ) );
            values.push_back( bodies.at( bodyName )->getAtmosphereModel( )->getTemperature( 300.0E3, 0.3, 0.2, time ) );
        }
        if( bodies.at( bodyName )->getMassProperties( ) != nullptr )
        {
            bodies.at( bodyName )->getMassProperties( )->update( time );
            values.push_back( bodies.at( bodyName )->getMassProperties( )->getCurrentMass( ) );
        }
        for( auto targetModel : bodies.at( bodyName )->getRadiationPressureTargetModels( ) )
        {
            if( std::dynamic_pointer_cast< CannonballRadiationPressureTargetModel >( targetModel ) != nullptr )
            {
                values.push_back( std::dynamic_pointer_cast< CannonballRadiationPressureTargetModel >( targetModel )->getCoefficient( ) );
            }
            else
            {
                bodies.at( bodyName )->getVehicleSystems( )->updatePartOrientations( time );
                targetModel->updateRadiationPressureForcing( 1361.0, Eigen::Vector3d( 0.6, -0.48, 0.64 ), true, "Sun" );
                appendVector( values, targetModel->getCurrentRadiationPressureForce( "Sun" ) );
            }
        }
        for( auto engineModel : bodies.at( bodyName )->getVehicleSystems( )->getEngineModels( ) )
        {
            engineModel.second->updateEngineModel( time );
            values.push_back( engineModel.second->getCurrentThrust( ) );
            values.push_back( engineModel.second->getCurrentMassRate( ) );
            appendVector( values, engineModel.second->getBodyFixedThrustDirection( ) );
        }
    }

    values.push_back( std::dynamic_pointer_cast< IsotropicPointRadiationSourceModel >( bodies.at( "Sun" )->getRadiationSourceModel( ) )
                              ->getLuminosityModel( )
                              ->getLuminosity( ) );
    bodies.at( "Earth" )->getRadiationSourceModel( )->updateMembers( time );
    for( auto irradiance : bodies.at( "Earth" )->getRadiationSourceModel( )->evaluateIrradianceAtPosition(
                 ( Eigen::Vector3d( ) << 5.0E6, 4.0E6, 3.0E6 ).finished( ) ) )
    {
        values.push_back( irradiance.first );
        appendVector( values, irradiance.second );
    }

    // Flight conditions, aerodynamic coefficients (with moment contribution) and aerodynamic angle based rotation
    std::shared_ptr< AtmosphericFlightConditions > flightConditions =
            std::dynamic_pointer_cast< AtmosphericFlightConditions >( bodies.at( "Vehicle" )->getFlightConditions( ) );
    flightConditions->resetCurrentTime( );
    flightConditions->updateConditions( time );
    values.push_back( flightConditions->getCurrentAltitude( ) );
    values.push_back( flightConditions->getCurrentDensity( ) );
    values.push_back( flightConditions->getCurrentAirspeed( ) );
    values.push_back( flightConditions->getCurrentMachNumber( ) );
    appendVector( values, bodies.at( "Vehicle" )->getRotationalEphemeris( )->getRotationToBaseFrame( time ).coeffs( ) );
    std::shared_ptr< AerodynamicCoefficientInterface > coefficientInterface = bodies.at( "Vehicle" )->getAerodynamicCoefficientInterface( );
    coefficientInterface->updateCurrentCoefficients( { flightConditions->getCurrentMachNumber( ) }, time );
    appendVector( values, coefficientInterface->getCurrentForceCoefficients( ) );
    appendVector( values, coefficientInterface->getCurrentMomentCoefficients( ) );

    std::shared_ptr< FlightConditions > flightConditions2 = bodies.at( "Vehicle2" )->getFlightConditions( );
    flightConditions2->resetCurrentTime( );
    flightConditions2->updateConditions( time );
    values.push_back( flightConditions2->getCurrentAltitude( ) );
    bodies.at( "Vehicle2" )->getAerodynamicCoefficientInterface( )->updateCurrentCoefficients( {}, time );
    appendVector( values, bodies.at( "Vehicle2" )->getAerodynamicCoefficientInterface( )->getCurrentForceCoefficients( ) );

    return values;
}

std::vector< double > evaluateEnvironmentHistory( const SystemOfBodies& bodies, const std::vector< double >& times )
{
    std::vector< double > values;
    for( const double time : times )
    {
        setBodyStates( bodies, time );
        std::vector< double > currentValues = evaluateEnvironment( bodies, time );
        values.insert( values.end( ), currentValues.begin( ), currentValues.end( ) );
    }
    return values;
}

std::vector< double > getTestTimes( const unsigned int numberOfTimes )
{
    std::vector< double > times;
    for( unsigned int i = 0; i < numberOfTimes; i++ )
    {
        // Non-monotonic times, so that interpolator look-up indices change by varying amounts
        times.push_back( initialTime +
                         static_cast< double >( ( i * 7919 ) % numberOfTimes ) / numberOfTimes * ( finalTime - initialTime ) );
    }
    return times;
}

// Check that a model exists in the copy if, and only if, it exists in the original, and that it is a different object
template< typename ModelType >
void checkModelIsCopied( const std::shared_ptr< ModelType > model, const std::shared_ptr< ModelType > modelCopy )
{
    BOOST_CHECK_EQUAL( model == nullptr, modelCopy == nullptr );
    if( model != nullptr )
    {
        BOOST_CHECK( model != modelCopy );
    }
}

// Check that a model exists in the copy if, and only if, it exists in the original, and that it is of the same type
template< typename ModelType >
void checkModelTypeIsEqual( const std::shared_ptr< ModelType > model, const std::shared_ptr< ModelType > modelCopy )
{
    BOOST_CHECK_EQUAL( model == nullptr, modelCopy == nullptr );
    if( model != nullptr && modelCopy != nullptr )
    {
        BOOST_CHECK( typeid( *model ) == typeid( *modelCopy ) );
    }
}

//! Test closure of environment models (links to other models): which models require closure, when it is applied, and that models
//! without closure cannot be used
BOOST_AUTO_TEST_CASE( testClosureOfEnvironmentModels )
{
    SystemOfBodies bodies = createSystemOfBodies( getTestBodySettings( ) );

    // All closures are applied, except the moment contribution of the aerodynamic force, which requires flight conditions
    std::map< std::string, std::map< EnvironmentModelClosureType, bool > > expectedClosureStates;
    expectedClosureStates[ "Earth" ] = { { gravity_field_variations_closure, true }, { radiation_source_closure, true } };
    expectedClosureStates[ "Moon" ] = { { rotation_model_closure, true } };
    expectedClosureStates[ "Vehicle" ] = { { rotation_model_closure, true }, { aerodynamic_coefficients_closure, false } };
    expectedClosureStates[ "Vehicle2" ] = { { rotation_model_closure, true } };
    expectedClosureStates[ "Vehicle3" ] = { { vehicle_systems_closure, true }, { radiation_pressure_target_closure, true } };
    for( auto bodyIterator : bodies.getMap( ) )
    {
        std::map< EnvironmentModelClosureType, bool > expectedStates;
        if( expectedClosureStates.count( bodyIterator.first ) > 0 )
        {
            expectedStates = expectedClosureStates.at( bodyIterator.first );
        }
        BOOST_CHECK( bodyIterator.second->getEnvironmentModelClosureStates( ) == expectedStates );
    }
    BOOST_CHECK( bodies.at( "Vehicle" )->getAerodynamicCoefficientInterface( )->getMomentContributionInterface( ) == nullptr );

    // Adding flight conditions applies their closure, and that of the moment contribution
    addAtmosphericFlightConditions( bodies, "Vehicle", "Earth" );
    BOOST_CHECK_EQUAL( getEnvironmentModelsWithoutClosure( bodies ).size( ), 0 );
    BOOST_CHECK( bodies.at( "Vehicle" )->getEnvironmentModelClosureStates( ).at( flight_conditions_closure ) );
    BOOST_CHECK( bodies.at( "Vehicle" )->getEnvironmentModelClosureStates( ).at( aerodynamic_coefficients_closure ) );
    BOOST_CHECK( bodies.at( "Vehicle" )->getAerodynamicCoefficientInterface( )->getMomentContributionInterface( ) != nullptr );

    // Copy requires and applies the same closures
    SystemOfBodies bodiesCopy = createTestBodiesCopy( bodies, false );
    for( auto bodyIterator : bodies.getMap( ) )
    {
        BOOST_CHECK( bodiesCopy.at( bodyIterator.first )->getEnvironmentModelClosureStates( ) ==
                     bodyIterator.second->getEnvironmentModelClosureStates( ) );
    }

    // Models without closure throw an exception when used
    BOOST_CHECK_THROW( SynchronousRotationalEphemeris( nullptr, "Earth", "J2000", "IAU_Moon" ).getRotationToBaseFrame( initialTime ),
                       std::runtime_error );
    AerodynamicAngleRotationalEphemeris aerodynamicAngleRotation( "Earth", "J2000", "VehicleFixed" );
    aerodynamicAngleRotation.setIsBodyInPropagation( true );
    BOOST_CHECK_THROW( aerodynamicAngleRotation.getRotationToBaseFrame( initialTime ), std::runtime_error );
    BOOST_CHECK_THROW( StateBasedBodyFixedDirectionCalculator( "Moon", true, false, nullptr ).update( initialTime ), std::runtime_error );
    BOOST_CHECK_THROW( FlightConditions( "Earth" ).updateConditions( initialTime ), std::runtime_error );
    BOOST_CHECK_THROW( AtmosphericFlightConditions( "Earth" ).updateConditions( initialTime ), std::runtime_error );
    BOOST_CHECK_THROW(
            PaneledRadiationPressureTargetModel( std::map< std::string, std::vector< std::string > >( ), std::map< std::string, int >( ) )
                    .updateRadiationPressureForcing( 1361.0, Eigen::Vector3d::UnitX( ), true, "Sun" ),
            std::runtime_error );
    BOOST_CHECK_THROW( system_models::VehicleExteriorPanel( std::function< Eigen::Vector3d( ) >( nullptr ), 1.0, "Sun" )
                               .updatePanel( Eigen::Quaterniond::Identity( ) ),
                       std::runtime_error );

    std::map< int, std::vector< std::complex< double > > > loveNumbers;
    loveNumbers[ 2 ] = { 0.3, 0.3, 0.3 };
    std::shared_ptr< GravityFieldVariations > tidalVariation = createGravityFieldVariationsModel(
            std::make_shared< BasicSolidBodyGravityFieldVariationSettings >( std::vector< std::string >{ "Moon" }, loveNumbers ),
            "Earth",
            bodies,
            false );
    Eigen::MatrixXd cosineCoefficients = Eigen::MatrixXd::Zero( 5, 5 ), sineCoefficients = Eigen::MatrixXd::Zero( 5, 5 );
    BOOST_CHECK_THROW( tidalVariation->addSphericalHarmonicsCorrections( initialTime, sineCoefficients, cosineCoefficients ),
                       std::runtime_error );

    std::shared_ptr< electromagnetism::RadiationSourceModel > paneledSource =
            createRadiationSourceModelWithoutClosure( getTestBodySettings( ).at( "Earth" )->radiationSourceModelSettings, "Earth", bodies );
    BOOST_CHECK_THROW( paneledSource->updateMembers( initialTime ), std::runtime_error );

    std::shared_ptr< AerodynamicCoefficientInterface > coefficientInterface = createAerodynamicCoefficientInterfaceWithoutClosure(
            getTestBodySettings( ).at( "Vehicle" )->aerodynamicCoefficientSettings, "Vehicle", bodies );
    BOOST_CHECK_THROW( coefficientInterface->updateFullCurrentCoefficients( { 10.0 }, {}, initialTime ), std::runtime_error );
}

//! Test which environment models are copied, and which are shared, by the copy
BOOST_AUTO_TEST_CASE( testCopiedAndSharedModels )
{
    SystemOfBodies bodies = createTestBodies( );
    SystemOfBodies bodiesCopy = createTestBodiesCopy( bodies );

    BOOST_CHECK_EQUAL( bodiesCopy.getMap( ).size( ), bodies.getMap( ).size( ) );
    BOOST_CHECK_EQUAL( bodiesCopy.getFrameOrigin( ), bodies.getFrameOrigin( ) );
    BOOST_CHECK_EQUAL( bodiesCopy.getFrameOrientation( ), bodies.getFrameOrientation( ) );

    for( const std::string& bodyName : bodyNames )
    {
        std::shared_ptr< Body > body = bodies.at( bodyName );
        std::shared_ptr< Body > bodyCopy = bodiesCopy.at( bodyName );
        BOOST_CHECK( body != bodyCopy );

        // Models without data that changes during propagation are shared
        const bool isEphemerisShared = ( bodyName == "Sun" );
        BOOST_CHECK_EQUAL( body->getEphemeris( ) == bodyCopy->getEphemeris( ), isEphemerisShared );
        const bool isRotationShared = ( bodyName == "Sun" || bodyName == "Vehicle3" );
        BOOST_CHECK_EQUAL( body->getRotationalEphemeris( ) == bodyCopy->getRotationalEphemeris( ), isRotationShared );
        BOOST_CHECK( body->getShapeModel( ) == bodyCopy->getShapeModel( ) );

        const bool isGravityFieldShared = ( bodyName == "Sun" );
        if( body->getGravityFieldModel( ) != nullptr )
        {
            BOOST_CHECK_EQUAL( body->getGravityFieldModel( ) == bodyCopy->getGravityFieldModel( ), isGravityFieldShared );
        }
        checkModelTypeIsEqual( body->getGravityFieldModel( ), bodyCopy->getGravityFieldModel( ) );
        checkModelTypeIsEqual( body->getEphemeris( ), bodyCopy->getEphemeris( ) );
        checkModelTypeIsEqual( body->getRotationalEphemeris( ), bodyCopy->getRotationalEphemeris( ) );

        // Models with data that changes during propagation are copied
        checkModelIsCopied( body->getAtmosphereModel( ), bodyCopy->getAtmosphereModel( ) );
        checkModelIsCopied( body->getAerodynamicCoefficientInterface( ), bodyCopy->getAerodynamicCoefficientInterface( ) );
        checkModelIsCopied( body->getFlightConditions( ), bodyCopy->getFlightConditions( ) );
        checkModelIsCopied( body->getMassProperties( ), bodyCopy->getMassProperties( ) );
        checkModelIsCopied( body->getVehicleSystems( ), bodyCopy->getVehicleSystems( ) );
        checkModelIsCopied( body->getRadiationSourceModel( ), bodyCopy->getRadiationSourceModel( ) );
        checkModelIsCopied( body->getGravityFieldVariationSet( ), bodyCopy->getGravityFieldVariationSet( ) );
        BOOST_CHECK_EQUAL( body->getRadiationPressureTargetModels( ).size( ), bodyCopy->getRadiationPressureTargetModels( ).size( ) );
        for( unsigned int i = 0; i < body->getRadiationPressureTargetModels( ).size( ); i++ )
        {
            BOOST_CHECK( body->getRadiationPressureTargetModels( ).at( i ) != bodyCopy->getRadiationPressureTargetModels( ).at( i ) );
        }
    }

    // Variation objects of time-dependent gravity fields are copied, and are those used by the copied gravity field
    for( const std::string bodyName : { "Earth", "Moon" } )
    {
        std::vector< std::shared_ptr< GravityFieldVariations > > variations =
                bodies.at( bodyName )->getGravityFieldVariationSet( )->getVariationObjects( );
        std::vector< std::shared_ptr< GravityFieldVariations > > variationsCopy =
                bodiesCopy.at( bodyName )->getGravityFieldVariationSet( )->getVariationObjects( );
        BOOST_CHECK_EQUAL( variations.size( ), 2 );
        BOOST_CHECK_EQUAL( variationsCopy.size( ), variations.size( ) );
        for( unsigned int i = 0; i < variations.size( ); i++ )
        {
            BOOST_CHECK( variations.at( i ) != variationsCopy.at( i ) );
            checkModelTypeIsEqual( variations.at( i ), variationsCopy.at( i ) );
        }
    }

    // Links between models are recreated in the copy
    std::shared_ptr< AerodynamicAngleRotationalEphemeris > vehicleRotationCopy =
            std::dynamic_pointer_cast< AerodynamicAngleRotationalEphemeris >( bodiesCopy.at( "Vehicle" )->getRotationalEphemeris( ) );
    BOOST_CHECK( vehicleRotationCopy->getAerodynamicAngleCalculator( ) ==
                 bodiesCopy.at( "Vehicle" )->getFlightConditions( )->getAerodynamicAngleCalculator( ) );
    BOOST_CHECK( vehicleRotationCopy->getAerodynamicAngleCalculator( ) !=
                 bodies.at( "Vehicle" )->getFlightConditions( )->getAerodynamicAngleCalculator( ) );
    BOOST_CHECK( std::dynamic_pointer_cast< AtmosphericFlightConditions >( bodiesCopy.at( "Vehicle" )->getFlightConditions( ) )
                         ->getAerodynamicCoefficientInterface( ) == bodiesCopy.at( "Vehicle" )->getAerodynamicCoefficientInterface( ) );
    BOOST_CHECK( bodiesCopy.at( "Vehicle" )->getAerodynamicCoefficientInterface( )->getMomentContributionInterface( ) !=
                 bodies.at( "Vehicle" )->getAerodynamicCoefficientInterface( )->getMomentContributionInterface( ) );
    BOOST_CHECK( bodiesCopy.at( "Vehicle" )->getAerodynamicCoefficientInterface( )->getMomentContributionInterface( ) != nullptr );
    // Paneled radiation source is recreated, with its original source (Sun) from the copy
    std::shared_ptr< DynamicallyPaneledRadiationSourceModel > paneledSourceCopy =
            std::dynamic_pointer_cast< DynamicallyPaneledRadiationSourceModel >( bodiesCopy.at( "Earth" )->getRadiationSourceModel( ) );
    BOOST_CHECK( paneledSourceCopy != nullptr );
    BOOST_CHECK( paneledSourceCopy != bodies.at( "Earth" )->getRadiationSourceModel( ) );
    BOOST_CHECK_EQUAL(
            paneledSourceCopy->getNumberOfPanels( ),
            std::dynamic_pointer_cast< DynamicallyPaneledRadiationSourceModel >( bodies.at( "Earth" )->getRadiationSourceModel( ) )
                    ->getNumberOfPanels( ) );
    BOOST_CHECK( paneledSourceCopy->getSourcePanelRadiosityModelUpdater( )->getOriginalSourceToSourceOccultingBodyNames( ) ==
                 std::set< std::string >( { "Moon" } ) );

    // Panels, part rotation models, panelled radiation pressure target and engine are copied, and linked to the copied vehicle systems
    std::shared_ptr< system_models::VehicleSystems > vehicleSystems = bodies.at( "Vehicle3" )->getVehicleSystems( );
    std::shared_ptr< system_models::VehicleSystems > vehicleSystemsCopy = bodiesCopy.at( "Vehicle3" )->getVehicleSystems( );
    BOOST_CHECK_EQUAL( vehicleSystemsCopy->getAllPanels( ).size( ), 4 );
    BOOST_CHECK_EQUAL( vehicleSystemsCopy->getAllPanels( ).size( ), vehicleSystems->getAllPanels( ).size( ) );
    for( unsigned int i = 0; i < vehicleSystems->getAllPanels( ).size( ); i++ )
    {
        BOOST_CHECK( vehicleSystemsCopy->getAllPanels( ).at( i ) != vehicleSystems->getAllPanels( ).at( i ) );
    }
    BOOST_CHECK( vehicleSystemsCopy->getVehiclePartOrientation( ).at( "SolarArray" ) ==
                 vehicleSystems->getVehiclePartOrientation( ).at( "SolarArray" ) );
    BOOST_CHECK( vehicleSystemsCopy->getEngineModels( ).at( "MainEngine" ) != vehicleSystems->getEngineModels( ).at( "MainEngine" ) );
    BOOST_CHECK( vehicleSystemsCopy->getEngineModels( ).at( "MainEngine" )->getThrustMagnitudeWrapper( ) !=
                 vehicleSystems->getEngineModels( ).at( "MainEngine" )->getThrustMagnitudeWrapper( ) );
    std::shared_ptr< PaneledRadiationPressureTargetModel > paneledTargetCopy =
            std::dynamic_pointer_cast< PaneledRadiationPressureTargetModel >(
                    bodiesCopy.at( "Vehicle3" )->getRadiationPressureTargetModels( ).at( 0 ) );
    BOOST_CHECK( paneledTargetCopy != nullptr );
    BOOST_CHECK_EQUAL( paneledTargetCopy->getAllPanels( ).size( ), 4 );
    for( unsigned int i = 0; i < paneledTargetCopy->getAllPanels( ).size( ); i++ )
    {
        BOOST_CHECK( paneledTargetCopy->getAllPanels( ).at( i ) == vehicleSystemsCopy->getAllPanels( ).at( i ) );
    }
    BOOST_CHECK( paneledTargetCopy->getSourceToTargetOccultingBodies( ) ==
                 bodies.at( "Vehicle3" )->getRadiationPressureTargetModels( ).at( 0 )->getSourceToTargetOccultingBodies( ) );

    BOOST_CHECK( std::dynamic_pointer_cast< DirectionBasedRotationalEphemeris >( bodiesCopy.at( "Vehicle2" )->getRotationalEphemeris( ) )
                         ->getInertialBodyAxisDirectionCalculator( ) !=
                 std::dynamic_pointer_cast< DirectionBasedRotationalEphemeris >( bodies.at( "Vehicle2" )->getRotationalEphemeris( ) )
                         ->getInertialBodyAxisDirectionCalculator( ) );
}

//! Test whether the copied environment gives results identical to the original (also for a copy of a copy)
BOOST_AUTO_TEST_CASE( testCopiedModelEvaluation )
{
    SystemOfBodies bodies = createTestBodies( );
    SystemOfBodies bodiesCopy = createTestBodiesCopy( bodies );
    SystemOfBodies bodiesCopyOfCopy = createTestBodiesCopy( bodiesCopy, false );

    std::vector< double > times = getTestTimes( 50 );
    std::vector< double > originalValues = evaluateEnvironmentHistory( bodies, times );
    std::vector< double > copyValues = evaluateEnvironmentHistory( bodiesCopy, times );
    std::vector< double > copyOfCopyValues = evaluateEnvironmentHistory( bodiesCopyOfCopy, times );

    BOOST_CHECK( originalValues.size( ) > 100 * times.size( ) );
    BOOST_CHECK_EQUAL( copyValues.size( ), originalValues.size( ) );
    BOOST_CHECK_EQUAL( copyOfCopyValues.size( ), originalValues.size( ) );
    for( unsigned int i = 0; i < originalValues.size( ); i++ )
    {
        BOOST_CHECK( !std::isnan( originalValues.at( i ) ) );
        BOOST_CHECK_EQUAL( copyValues.at( i ), originalValues.at( i ) );
        BOOST_CHECK_EQUAL( copyOfCopyValues.at( i ), originalValues.at( i ) );
    }
}

//! Test whether models in the copy are linked to the copied bodies only, and changes to one system of bodies do not affect the other
BOOST_AUTO_TEST_CASE( testIndependenceOfCopy )
{
    SystemOfBodies bodies = createTestBodies( );
    SystemOfBodies bodiesCopy = createTestBodiesCopy( bodies );

    const double time = initialTime + 20000.0;
    setBodyStates( bodies, time );
    const std::vector< double > originalValues = evaluateEnvironment( bodies, time );

    // Different states in the copy, linked models in the copy must use them, models in the original must not
    setBodyStates( bodiesCopy, time );
    const Eigen::Vector6d vehicleStateChange = ( Eigen::Vector6d( ) << 1.0E3, -2.0E3, 3.0E3, 1.0, -2.0, 3.0 ).finished( );
    bodiesCopy.at( "Vehicle" )->setState( bodiesCopy.at( "Vehicle" )->getState( ) + vehicleStateChange );
    bodiesCopy.at( "Vehicle2" )->setState( bodiesCopy.at( "Vehicle2" )->getState( ) + vehicleStateChange );
    bodiesCopy.at( "Moon" )->setState( bodiesCopy.at( "Moon" )->getState( ) + 1.0E3 * vehicleStateChange );

    std::shared_ptr< FlightConditions > flightConditions = bodies.at( "Vehicle" )->getFlightConditions( );
    std::shared_ptr< FlightConditions > flightConditionsCopy = bodiesCopy.at( "Vehicle" )->getFlightConditions( );
    flightConditions->resetCurrentTime( );
    flightConditionsCopy->resetCurrentTime( );
    flightConditionsCopy->updateConditions( time );
    flightConditions->updateConditions( time );
    BOOST_CHECK_EQUAL( flightConditions->getCurrentAltitude( ), originalValues.at( originalValues.size( ) - 18 ) );
    BOOST_CHECK( std::fabs( flightConditionsCopy->getCurrentAltitude( ) - flightConditions->getCurrentAltitude( ) ) > 100.0 );

    // Vehicle2 rotation model (orbital-state based) and Moon rotation model (synchronous, during propagation)
    bodiesCopy.at( "Vehicle2" )->getRotationalEphemeris( )->resetCurrentTime( );
    BOOST_CHECK( !bodiesCopy.at( "Vehicle2" )
                          ->getRotationalEphemeris( )
                          ->getRotationToBaseFrame( time )
                          .isApprox( bodies.at( "Vehicle2" )->getRotationalEphemeris( )->getRotationToBaseFrame( time ), 1.0E-6 ) );
    std::dynamic_pointer_cast< SynchronousRotationalEphemeris >( bodies.at( "Moon" )->getRotationalEphemeris( ) )
            ->setIsBodyInPropagation( true );
    std::dynamic_pointer_cast< SynchronousRotationalEphemeris >( bodiesCopy.at( "Moon" )->getRotationalEphemeris( ) )
            ->setIsBodyInPropagation( true );
    Eigen::Quaterniond moonRotation = bodies.at( "Moon" )->getRotationalEphemeris( )->getRotationToBaseFrame( time );
    BOOST_CHECK( !moonRotation.isApprox( bodiesCopy.at( "Moon" )->getRotationalEphemeris( )->getRotationToBaseFrame( time ), 1.0E-8 ) );
    std::dynamic_pointer_cast< SynchronousRotationalEphemeris >( bodies.at( "Moon" )->getRotationalEphemeris( ) )
            ->setIsBodyInPropagation( false );

    // Tidal variations of Earth use Moon state and gravitational parameter of the copy
    std::shared_ptr< TimeDependentSphericalHarmonicsGravityField > earthField =
            std::dynamic_pointer_cast< TimeDependentSphericalHarmonicsGravityField >( bodies.at( "Earth" )->getGravityFieldModel( ) );
    std::shared_ptr< TimeDependentSphericalHarmonicsGravityField > earthFieldCopy =
            std::dynamic_pointer_cast< TimeDependentSphericalHarmonicsGravityField >( bodiesCopy.at( "Earth" )->getGravityFieldModel( ) );
    earthField->update( time );
    earthFieldCopy->update( time );
    const Eigen::MatrixXd earthCosineCoefficients = earthField->getCosineCoefficients( );
    BOOST_CHECK( ( earthFieldCopy->getCosineCoefficients( ) - earthCosineCoefficients ).cwiseAbs( ).maxCoeff( ) > 1.0E-15 );
    setBodyStates( bodiesCopy, time );
    earthFieldCopy->update( time );
    BOOST_CHECK( earthFieldCopy->getCosineCoefficients( ) == earthCosineCoefficients );
    bodiesCopy.at( "Moon" )->getGravityFieldModel( )->resetGravitationalParameter( 2.0 *
                                                                                   bodies.at( "Moon" )->getGravitationalParameter( ) );
    earthField->update( time );
    earthFieldCopy->update( time );
    BOOST_CHECK( earthField->getCosineCoefficients( ) == earthCosineCoefficients );
    BOOST_CHECK( ( earthFieldCopy->getCosineCoefficients( ) - earthCosineCoefficients ).cwiseAbs( ).maxCoeff( ) > 1.0E-12 );

    // Moment contribution uses center of mass of copy
    addRigidBodyProperties( bodiesCopy,
                            "Vehicle",
                            constantRigidBodyPropertiesSettings(
                                    500.0, ( Eigen::Vector3d( ) << 1.0, 0.0, 0.0 ).finished( ), Eigen::Matrix3d::Identity( ) * 100.0 ) );
    BOOST_CHECK(
            bodies.at( "Vehicle" )->getAerodynamicCoefficientInterface( )->getMomentContributionInterface( )->centerOfMassPosition_( ) ==
            ( Eigen::Vector3d( ) << 0.1, 0.0, 0.0 ).finished( ) );
    BOOST_CHECK( bodiesCopy.at( "Vehicle" )
                         ->getAerodynamicCoefficientInterface( )
                         ->getMomentContributionInterface( )
                         ->centerOfMassPosition_( ) == ( Eigen::Vector3d( ) << 1.0, 0.0, 0.0 ).finished( ) );

    // Tracking panel and paneled radiation source of the copy use the Sun position of the copy
    std::shared_ptr< system_models::VehicleExteriorPanel > trackingPanel =
            bodies.at( "Vehicle3" )->getVehicleSystems( )->getAllPanels( ).at( 2 );
    std::shared_ptr< system_models::VehicleExteriorPanel > trackingPanelCopy =
            bodiesCopy.at( "Vehicle3" )->getVehicleSystems( )->getAllPanels( ).at( 2 );
    BOOST_CHECK_EQUAL( trackingPanelCopy->getTrackedBody( ), "Sun" );
    const Eigen::Vector3d trackingPanelNormal = trackingPanel->getFrameFixedSurfaceNormal( )( );
    BOOST_CHECK( trackingPanelCopy->getFrameFixedSurfaceNormal( )( ) == trackingPanelNormal );
    bodies.at( "Earth" )->getRadiationSourceModel( )->updateMembers( time + 1.0 );
    bodiesCopy.at( "Earth" )->getRadiationSourceModel( )->updateMembers( time + 1.0 );
    const Eigen::Vector3d irradiancePosition = ( Eigen::Vector3d( ) << 5.0E6, 4.0E6, 3.0E6 ).finished( );
    const double earthIrradiance = getTotalIrradiance( bodies.at( "Earth" )->getRadiationSourceModel( ), irradiancePosition );
    BOOST_CHECK_EQUAL( getTotalIrradiance( bodiesCopy.at( "Earth" )->getRadiationSourceModel( ), irradiancePosition ), earthIrradiance );

    const Eigen::Vector6d sunStateChange = ( Eigen::Vector6d( ) << 3.0E10, -4.0E10, 2.0E10, 0.0, 0.0, 0.0 ).finished( );
    bodiesCopy.at( "Sun" )->setState( bodiesCopy.at( "Sun" )->getState( ) + sunStateChange );
    BOOST_CHECK( ( trackingPanelCopy->getFrameFixedSurfaceNormal( )( ) - trackingPanelNormal ).norm( ) > 1.0E-2 );
    BOOST_CHECK( trackingPanel->getFrameFixedSurfaceNormal( )( ) == trackingPanelNormal );
    bodies.at( "Earth" )->getRadiationSourceModel( )->updateMembers( time + 2.0 );
    bodiesCopy.at( "Earth" )->getRadiationSourceModel( )->updateMembers( time + 2.0 );
    BOOST_CHECK_EQUAL( getTotalIrradiance( bodies.at( "Earth" )->getRadiationSourceModel( ), irradiancePosition ), earthIrradiance );
    BOOST_CHECK( std::fabs( getTotalIrradiance( bodiesCopy.at( "Earth" )->getRadiationSourceModel( ), irradiancePosition ) -
                            earthIrradiance ) > 1.0E-6 * std::fabs( earthIrradiance ) );
    bodiesCopy.at( "Sun" )->setState( bodiesCopy.at( "Sun" )->getState( ) - sunStateChange );

    // Orientation of panels in the copy is updated independently of the original
    bodies.at( "Vehicle3" )->getVehicleSystems( )->updatePartOrientations( time );
    bodiesCopy.at( "Vehicle3" )->getVehicleSystems( )->updatePartOrientations( time + 600.0 );
    std::shared_ptr< RadiationPressureTargetModel > paneledTarget = bodies.at( "Vehicle3" )->getRadiationPressureTargetModels( ).at( 0 );
    std::shared_ptr< RadiationPressureTargetModel > paneledTargetCopy =
            bodiesCopy.at( "Vehicle3" )->getRadiationPressureTargetModels( ).at( 0 );
    paneledTargetCopy->updateRadiationPressureForcing( 1361.0, Eigen::Vector3d( 0.0, 1.0, 0.0 ), true, "Sun" );
    paneledTarget->updateRadiationPressureForcing( 1361.0, Eigen::Vector3d( 0.0, 1.0, 0.0 ), true, "Sun" );
    const Eigen::Vector3d paneledForce = paneledTarget->getCurrentRadiationPressureForce( "Sun" );
    BOOST_CHECK( ( paneledTargetCopy->getCurrentRadiationPressureForce( "Sun" ) - paneledForce ).norm( ) > 1.0E-3 * paneledForce.norm( ) );
    bodiesCopy.at( "Vehicle3" )->getVehicleSystems( )->updatePartOrientations( time );
    paneledTargetCopy->updateRadiationPressureForcing( 1361.0, Eigen::Vector3d( 0.0, 1.0, 0.0 ), true, "Sun" );
    BOOST_CHECK( paneledTargetCopy->getCurrentRadiationPressureForce( "Sun" ) == paneledForce );

    // Changes to model parameters in the copy do not affect the original, and vice versa
    Eigen::MatrixXd moonCosineCoefficients =
            std::dynamic_pointer_cast< TimeDependentSphericalHarmonicsGravityField >( bodies.at( "Moon" )->getGravityFieldModel( ) )
                    ->getNominalCosineCoefficients( );
    Eigen::MatrixXd perturbedCoefficients = moonCosineCoefficients;
    perturbedCoefficients( 2, 0 ) *= 2.0;
    std::dynamic_pointer_cast< TimeDependentSphericalHarmonicsGravityField >( bodiesCopy.at( "Moon" )->getGravityFieldModel( ) )
            ->setNominalCosineCoefficients( perturbedCoefficients );
    BOOST_CHECK( std::dynamic_pointer_cast< TimeDependentSphericalHarmonicsGravityField >( bodies.at( "Moon" )->getGravityFieldModel( ) )
                         ->getNominalCosineCoefficients( )( 2, 0 ) == moonCosineCoefficients( 2, 0 ) );

    std::shared_ptr< ExponentialAtmosphere > marsAtmosphere = std::dynamic_pointer_cast< ExponentialAtmosphere >(
            std::dynamic_pointer_cast< ScaledAtmosphereModel >( bodies.at( "Mars" )->getAtmosphereModel( ) )->getBaseAtmosphere( ) );
    std::shared_ptr< ExponentialAtmosphere > marsAtmosphereCopy = std::dynamic_pointer_cast< ExponentialAtmosphere >(
            std::dynamic_pointer_cast< ScaledAtmosphereModel >( bodiesCopy.at( "Mars" )->getAtmosphereModel( ) )->getBaseAtmosphere( ) );
    BOOST_CHECK( marsAtmosphere != marsAtmosphereCopy );
    marsAtmosphereCopy->resetBaseDensity( 0.04 );
    BOOST_CHECK_EQUAL( marsAtmosphere->getDensity( 0.0, 0.0, 0.0, time ), 0.02 );
    marsAtmosphere->resetScaleHeight( 5.0E3 );
    BOOST_CHECK_EQUAL( marsAtmosphereCopy->getScaleHeight( ), 11.1E3 );

    std::shared_ptr< CannonballRadiationPressureTargetModel > targetModel =
            std::dynamic_pointer_cast< CannonballRadiationPressureTargetModel >(
                    bodies.at( "Vehicle" )->getRadiationPressureTargetModels( ).at( 0 ) );
    std::shared_ptr< CannonballRadiationPressureTargetModel > targetModelCopy =
            std::dynamic_pointer_cast< CannonballRadiationPressureTargetModel >(
                    bodiesCopy.at( "Vehicle" )->getRadiationPressureTargetModels( ).at( 0 ) );
    targetModelCopy->resetCoefficient( 1.8 );
    BOOST_CHECK_EQUAL( targetModel->getCoefficient( ), 1.3 );
    BOOST_CHECK_EQUAL( targetModelCopy->getCoefficient( ), 1.8 );

    bodiesCopy.at( "Mars" )->getGravityFieldModel( )->resetGravitationalParameter( 1.0 );
    BOOST_CHECK_EQUAL( bodies.at( "Mars" )->getGravitationalParameter( ), 4.2828E13 );

    // Original environment gives the same results as before the changes to the copy
    flightConditions->resetCurrentTime( );
    std::dynamic_pointer_cast< AerodynamicAngleRotationalEphemeris >( bodies.at( "Vehicle" )->getRotationalEphemeris( ) )
            ->resetCurrentTime( );
    marsAtmosphere->resetScaleHeight( 11.1E3 );
    setBodyStates( bodies, time );
    std::vector< double > originalValuesAfterChanges = evaluateEnvironment( bodies, time );
    BOOST_CHECK_EQUAL( originalValuesAfterChanges.size( ), originalValues.size( ) );
    for( unsigned int i = 0; i < originalValues.size( ); i++ )
    {
        BOOST_CHECK_EQUAL( originalValuesAfterChanges.at( i ), originalValues.at( i ) );
    }
}

//! Test closure and copy of ground station (with deformation and relativistic motion), tidal deformation, space-time metric,
//! parameterized thrust and approximate ephemeris
BOOST_AUTO_TEST_CASE( testCopyOfStationsDeformationMetricAndOtherModels )
{
    BodyListSettings bodySettings = getTestBodySettings( );
    const Eigen::Vector3d stationPosition( 4.0E6, 3.0E6, 3.5E6 );
    bodySettings.at( "Earth" )->bodyDeformationSettings.push_back(
            basicTidalBodyShapeDeformation( { "Moon", "Sun" }, { { 2, { 0.6, 0.08 } } } ) );
    bodySettings.at( "Earth" )->groundStationSettings.push_back(
            groundStationSettings( "Station",
                                   stationPosition,
                                   coordinate_conversions::cartesian_position,
                                   { bodyDeformationStationMotionSettings( ), bodycentricToBarycentricStationMotionSettings( "Sun" ) } ) );
    SystemOfBodies bodies = createTestBodies( bodySettings );
    createBaseMetric( std::make_shared< SolarSystemSpaceTimeMetricSettings >( std::vector< std::string >{ "Sun", "Earth", "Moon" } ),
                      bodies );
    addEngineModel( "Vehicle",
                    "AltitudeEngine",
                    std::make_shared< ParameterizedThrustMagnitudeSettings >(
                            std::make_shared< interpolators::LinearInterpolator< double, double > >(
                                    std::map< double, double >( { { 0.0, 1.0 }, { 1.0E6, 2.0 } } ) ),
                            std::vector< propulsion::ThrustIndependentVariables >{ propulsion::altitude_dependent_thrust },
                            300.0 ),
                    bodies );

    SystemOfBodies bodiesCopy = createTestBodiesCopy( bodies, false );

    // Closures of the models in the copy are registered and applied
    std::map< EnvironmentModelClosureType, bool > earthClosureStates = bodiesCopy.at( "Earth" )->getEnvironmentModelClosureStates( );
    BOOST_CHECK( earthClosureStates.at( ground_stations_closure ) );
    BOOST_CHECK( earthClosureStates.at( body_deformation_closure ) );
    BOOST_CHECK( bodiesCopy.at( "Vehicle" )->getEnvironmentModelClosureStates( ).at( vehicle_systems_closure ) );
    BOOST_CHECK_EQUAL( getEnvironmentModelsWithoutClosure( bodiesCopy ).size( ), 0 );

    // Models are copied
    std::shared_ptr< ground_stations::GroundStation > station = bodies.at( "Earth" )->getGroundStation( "Station" );
    std::shared_ptr< ground_stations::GroundStation > stationCopy = bodiesCopy.at( "Earth" )->getGroundStation( "Station" );
    BOOST_CHECK( station != stationCopy );
    BOOST_CHECK( station->getNominalStationState( ) != stationCopy->getNominalStationState( ) );
    BOOST_CHECK( station->getPointingAnglesCalculator( ) != stationCopy->getPointingAnglesCalculator( ) );
    BOOST_CHECK( station->getNominalStationState( )->getStationMotionModel( ) !=
                 stationCopy->getNominalStationState( )->getStationMotionModel( ) );
    std::shared_ptr< basic_astrodynamics::BodyDeformationModel > deformation = bodies.at( "Earth" )->getBodyDeformationModels( ).at( 0 );
    std::shared_ptr< basic_astrodynamics::BodyDeformationModel > deformationCopy =
            bodiesCopy.at( "Earth" )->getBodyDeformationModels( ).at( 0 );
    BOOST_CHECK( deformation != deformationCopy );
    std::shared_ptr< SpaceTimeProperties > spaceTimeProperties = bodies.getSpaceTimeProperties( );
    std::shared_ptr< SpaceTimeProperties > spaceTimePropertiesCopy = bodiesCopy.getSpaceTimeProperties( );
    BOOST_CHECK( spaceTimeProperties != spaceTimePropertiesCopy );
    BOOST_CHECK( spaceTimeProperties->getPpnParameterSet( ) != spaceTimePropertiesCopy->getPpnParameterSet( ) );
    BOOST_CHECK( spaceTimePropertiesCopy->getBaseMetric( ) != nullptr );
    BOOST_CHECK( spaceTimeProperties->getBaseMetric( ) != spaceTimePropertiesCopy->getBaseMetric( ) );
    std::shared_ptr< system_models::EngineModel > engine =
            bodies.at( "Vehicle" )->getVehicleSystems( )->getEngineModels( ).at( "AltitudeEngine" );
    std::shared_ptr< system_models::EngineModel > engineCopy =
            bodiesCopy.at( "Vehicle" )->getVehicleSystems( )->getEngineModels( ).at( "AltitudeEngine" );
    BOOST_CHECK( engine->getThrustMagnitudeWrapper( ) != engineCopy->getThrustMagnitudeWrapper( ) );

    // Evaluation of the copied models gives the same results as the original
    const double time = initialTime + 30000.0;
    const Eigen::Vector3d inertialVector( 1.0, -2.0, 0.5 );
    std::vector< std::vector< double > > values;
    for( const SystemOfBodies& currentBodies : { bodies, bodiesCopy } )
    {
        setBodyStates( currentBodies, time );
        std::vector< double > currentValues;
        std::shared_ptr< ground_stations::GroundStation > currentStation = currentBodies.at( "Earth" )->getGroundStation( "Station" );
        appendVector( currentValues, currentStation->getNominalStationState( )->getCartesianStateInTime( time, "SSB" ) );
        currentValues.push_back(
                currentStation->getPointingAnglesCalculator( )->calculateElevationAngleFromInertialVector( inertialVector, time ) );
        appendVector( currentValues,
                      currentBodies.at( "Earth" )->getBodyDeformationModels( ).at( 0 )->calculateDisplacement( time, stationPosition ) );

        currentBodies.at( "Vehicle" )->getFlightConditions( )->resetCurrentTime( );
        currentBodies.at( "Vehicle" )->getFlightConditions( )->updateConditions( time );
        std::shared_ptr< system_models::EngineModel > currentEngine =
                currentBodies.at( "Vehicle" )->getVehicleSystems( )->getEngineModels( ).at( "AltitudeEngine" );
        currentEngine->updateEngineModel( time );
        currentValues.push_back( currentEngine->getCurrentThrust( ) );

        values.push_back( currentValues );
    }
    BOOST_CHECK_EQUAL( values.at( 0 ).size( ), 11 );
    BOOST_CHECK_EQUAL( values.at( 1 ).size( ), values.at( 0 ).size( ) );
    for( unsigned int i = 0; i < values.at( 0 ).size( ); i++ )
    {
        BOOST_CHECK( !std::isnan( values.at( 0 ).at( i ) ) );
        BOOST_CHECK_CLOSE_FRACTION( values.at( 1 ).at( i ), values.at( 0 ).at( i ), 1.0E-14 );
    }
    BOOST_CHECK( values.at( 0 ).at( 10 ) > 1.0 );

    // Copied models are linked to the copied bodies: changes in the copy change the copied models only (if the original were linked
    // to the copied bodies, the original and copied models would give equal results)
    const Eigen::Vector3d displacement = deformation->calculateDisplacement( time + 1.0, stationPosition );
    bodiesCopy.at( "Moon" )->getGravityFieldModel( )->resetGravitationalParameter(
            2.0 * bodiesCopy.at( "Moon" )->getGravityFieldModel( )->getGravitationalParameter( ) );
    BOOST_CHECK( ( deformationCopy->calculateDisplacement( time + 1.0, stationPosition ) - displacement ).norm( ) >
                 1.0E-3 * displacement.norm( ) );
    BOOST_CHECK( ( deformation->calculateDisplacement( time + 2.0, stationPosition ) -
                   deformationCopy->calculateDisplacement( time + 2.0, stationPosition ) )
                         .norm( ) > 1.0E-3 * displacement.norm( ) );

    // (evaluation of tidal deformation above resets the current states of the bodies involved)
    setBodyStates( bodies, time );
    setBodyStates( bodiesCopy, time );
    const double thrust = engine->getCurrentThrust( );
    const Eigen::Vector3d vehicleRadialDirection =
            ( bodiesCopy.at( "Vehicle" )->getPosition( ) - bodiesCopy.at( "Earth" )->getPosition( ) ).normalized( );
    bodiesCopy.at( "Vehicle" )
            ->setState( bodiesCopy.at( "Vehicle" )->getState( ) +
                        ( Eigen::Vector6d( ) << 1.0E4 * vehicleRadialDirection, Eigen::Vector3d::Zero( ) ).finished( ) );
    for( const SystemOfBodies& currentBodies : { bodies, bodiesCopy } )
    {
        currentBodies.at( "Vehicle" )->getFlightConditions( )->resetCurrentTime( );
        currentBodies.at( "Vehicle" )->getFlightConditions( )->updateConditions( time );
        currentBodies.at( "Vehicle" )->getVehicleSystems( )->getEngineModels( ).at( "AltitudeEngine" )->resetCurrentTime( );
        currentBodies.at( "Vehicle" )->getVehicleSystems( )->getEngineModels( ).at( "AltitudeEngine" )->updateEngineModel( time );
    }
    BOOST_CHECK_CLOSE_FRACTION( engineCopy->getCurrentThrust( ) - thrust, 1.0E4 / 1.0E6, 1.0E-2 );
    BOOST_CHECK_EQUAL( engine->getCurrentThrust( ), thrust );

    spaceTimePropertiesCopy->getPpnParameterSet( )->setParameterGamma( 0.9 );
    BOOST_CHECK_EQUAL( spaceTimeProperties->getPpnParameterSet( )->getParameterGamma( ), 1.0 );

    // Approximate JPL ephemeris (defined w.r.t. ecliptic) is recreated in the copy
    BodyListSettings eclipticBodySettings( "SSB", "ECLIPJ2000" );
    eclipticBodySettings.addSettings( "Sun" );
    eclipticBodySettings.at( "Sun" )->ephemerisSettings = constantEphemerisSettings( Eigen::Vector6d::Zero( ), "SSB", "ECLIPJ2000" );
    eclipticBodySettings.addSettings( "Jupiter" );
    eclipticBodySettings.at( "Jupiter" )->ephemerisSettings = approximateJplEphemerisSettings( "Jupiter" );
    SystemOfBodies eclipticBodies = createSystemOfBodies( eclipticBodySettings );
    SystemOfBodies eclipticBodiesCopy = createSystemOfBodiesCopy( eclipticBodies, false );
    BOOST_CHECK( eclipticBodies.at( "Jupiter" )->getEphemeris( ) != eclipticBodiesCopy.at( "Jupiter" )->getEphemeris( ) );
    BOOST_CHECK( std::dynamic_pointer_cast< ApproximateJplEphemeris >( eclipticBodiesCopy.at( "Jupiter" )->getEphemeris( ) ) != nullptr );
    const Eigen::Vector6d jupiterState = eclipticBodies.at( "Jupiter" )->getEphemeris( )->getCartesianState( time );
    const Eigen::Vector6d jupiterStateCopy = eclipticBodiesCopy.at( "Jupiter" )->getEphemeris( )->getCartesianState( time );
    for( int i = 0; i < 6; i++ )
    {
        BOOST_CHECK_EQUAL( jupiterStateCopy( i ), jupiterState( i ) );
    }
}

//! Test closure and copy of coma (density and wind) model
BOOST_AUTO_TEST_CASE( testCopyOfComaModel )
{
    const std::filesystem::path comaDataDirectory = std::filesystem::path( paths::getTudatTestDataPath( ) ) / "coma";
    const std::filesystem::path windDataDirectory = comaDataDirectory / "wind" / "polynomial";
    ComaWindModelFileProcessor windFileProcessor( std::vector< std::string >{ ( windDataDirectory / "input_poly_wind_x.txt" ).string( ) },
                                                  std::vector< std::string >{ ( windDataDirectory / "input_poly_wind_y.txt" ).string( ) },
                                                  std::vector< std::string >{ ( windDataDirectory / "input_poly_wind_z.txt" ).string( ) } );

    BodyListSettings bodySettings( "SSB", "J2000" );
    bodySettings.addSettings( "Sun" );
    bodySettings.at( "Sun" )->ephemerisSettings = constantEphemerisSettings( Eigen::Vector6d::Zero( ), "SSB", "J2000" );
    bodySettings.addSettings( "Comet" );
    bodySettings.at( "Comet" )->ephemerisSettings = keplerEphemerisSettings(
            ( Eigen::Vector6d( ) << 4.5E11, 0.6, 0.1, 0.4, 1.0, 0.5 ).finished( ), initialTime, 1.32712440018E20, "Sun", "J2000" );
    bodySettings.at( "Comet" )->rotationModelSettings =
            simpleRotationModelSettings( "J2000", "CometFixed", Eigen::Quaterniond::Identity( ), initialTime, 1.0E-4 );
    bodySettings.at( "Comet" )->atmosphereSettings =
            comaSettings( ComaPolyDatasetReader::readFromFiles(
                                  { ( comaDataDirectory / "density" / "polynomial" / "input_poly_coef_test_file.txt" ).string( ) } ),
                          0.018 );
    bodySettings.at( "Comet" )->atmosphereSettings->setWindSettings( comaWindModelSettings( windFileProcessor.createPolyCoefDatasets( ) ) );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    SystemOfBodies bodiesCopy = createSystemOfBodiesCopy( bodies, false );

    // Coma and wind models are copied, and closed in the copy
    BOOST_CHECK( bodies.at( "Comet" )->getEnvironmentModelClosureStates( ).at( atmosphere_closure ) );
    BOOST_CHECK( bodiesCopy.at( "Comet" )->getEnvironmentModelClosureStates( ).at( atmosphere_closure ) );
    std::shared_ptr< AtmosphereModel > comaModel = bodies.at( "Comet" )->getAtmosphereModel( );
    std::shared_ptr< AtmosphereModel > comaModelCopy = bodiesCopy.at( "Comet" )->getAtmosphereModel( );
    BOOST_CHECK( std::dynamic_pointer_cast< ComaModel >( comaModelCopy ) != nullptr );
    BOOST_CHECK( comaModelCopy != comaModel );
    BOOST_CHECK( std::dynamic_pointer_cast< ComaWindModel >( comaModelCopy->getWindModel( ) ) != nullptr );
    BOOST_CHECK( comaModelCopy->getWindModel( ) != comaModel->getWindModel( ) );
    BOOST_CHECK( std::dynamic_pointer_cast< ComaModel >( comaModelCopy )->getSurfaceSphericalHarmonicsCalculator( ) !=
                 std::dynamic_pointer_cast< ComaModel >( comaModel )->getSurfaceSphericalHarmonicsCalculator( ) );

    // Copy gives the same density and wind as the original
    const double radius = 6000.0, longitude = 0.3, latitude = 0.2;
    // Times within the period covered by the coma data file (2015/07/21 - 2015/08/21)
    std::vector< double > times = { 490708800.0 + 1000.0, 490708800.0 + 50000.0 };
    for( const SystemOfBodies& currentBodies : { bodies, bodiesCopy } )
    {
        for( const std::string bodyName : { "Sun", "Comet" } )
        {
            currentBodies.at( bodyName )->setStateFromEphemeris< double, double >( times.at( 0 ) );
        }
        currentBodies.at( "Comet" )->setCurrentRotationalStateToLocalFrameFromEphemeris( times.at( 0 ) );
    }
    const double density = comaModel->getDensity( radius, longitude, latitude, times.at( 0 ) );
    const Eigen::Vector3d windVelocity =
            comaModel->getWindModel( )->getCurrentBodyFixedCartesianWindVelocity( radius, longitude, latitude, times.at( 0 ) );
    BOOST_CHECK( density > 0.0 );
    BOOST_CHECK( windVelocity.norm( ) > 0.0 );
    BOOST_CHECK_EQUAL( comaModelCopy->getDensity( radius, longitude, latitude, times.at( 0 ) ), density );
    BOOST_CHECK( comaModelCopy->getWindModel( )->getCurrentBodyFixedCartesianWindVelocity( radius, longitude, latitude, times.at( 0 ) ) ==
                 windVelocity );

    // Copied models use the orientation of the copied comet (if the original were linked to the copied comet, the original and
    // copied models would give equal results)
    for( const SystemOfBodies& currentBodies : { bodies, bodiesCopy } )
    {
        for( const std::string bodyName : { "Sun", "Comet" } )
        {
            currentBodies.at( bodyName )->setStateFromEphemeris< double, double >( times.at( 1 ) );
        }
    }
    bodies.at( "Comet" )->setCurrentRotationalStateToLocalFrameFromEphemeris( times.at( 1 ) );
    bodiesCopy.at( "Comet" )->setCurrentRotationalStateToLocalFrameFromEphemeris( times.at( 1 ) + 20000.0 );
    BOOST_CHECK( comaModelCopy->getDensity( radius, longitude, latitude, times.at( 1 ) ) !=
                 comaModel->getDensity( radius, longitude, latitude, times.at( 1 ) ) );
    BOOST_CHECK( ( comaModelCopy->getWindModel( )->getCurrentBodyFixedCartesianWindVelocity( radius, longitude, latitude, times.at( 1 ) ) -
                   comaModel->getWindModel( )->getCurrentBodyFixedCartesianWindVelocity( radius, longitude, latitude, times.at( 1 ) ) )
                         .norm( ) > 0.0 );
}

//! Test whether the original and its copies can be evaluated concurrently, giving results identical to serial evaluation
BOOST_AUTO_TEST_CASE( testConcurrentEvaluationOfCopies )
{
    // First use of Runge-Kutta coefficients (initialized on first use) from several threads at once
    std::vector< std::thread > coefficientThreads;
    std::vector< numerical_integrators::RungeKuttaCoefficients > threadCoefficients( 4 );
    for( unsigned int i = 0; i < threadCoefficients.size( ); i++ )
    {
        coefficientThreads.push_back( std::thread( [ &threadCoefficients, i ]( ) {
            threadCoefficients.at( i ) =
                    numerical_integrators::RungeKuttaCoefficients::get( numerical_integrators::CoefficientSets::rungeKuttaFehlberg78 );
        } ) );
    }
    for( std::thread& thread : coefficientThreads )
    {
        thread.join( );
    }
    for( unsigned int i = 0; i < threadCoefficients.size( ); i++ )
    {
        BOOST_CHECK( threadCoefficients.at( i ).bCoefficients ==
                     numerical_integrators::RungeKuttaCoefficients::get( numerical_integrators::CoefficientSets::rungeKuttaFehlberg78 )
                             .bCoefficients );
    }

    SystemOfBodies bodies = createTestBodies( );
    const unsigned int numberOfThreads = 4;
    std::vector< SystemOfBodies > bodyCopies;
    for( unsigned int i = 0; i < numberOfThreads; i++ )
    {
        bodyCopies.push_back( createTestBodiesCopy( bodies, false ) );
    }

    std::vector< double > times = getTestTimes( 400 );
    std::vector< double > serialValues = evaluateEnvironmentHistory( bodies, times );

    // Evaluate the original (in this thread) and all copies (in separate threads) at the same time
    for( unsigned int iteration = 0; iteration < 2; iteration++ )
    {
        std::vector< std::vector< double > > threadValues( numberOfThreads );
        std::vector< std::thread > threads;
        for( unsigned int i = 0; i < numberOfThreads; i++ )
        {
            threads.push_back( std::thread( [ &threadValues, &bodyCopies, &times, i ]( ) {
                threadValues.at( i ) = evaluateEnvironmentHistory( bodyCopies.at( i ), times );
            } ) );
        }
        std::vector< double > originalValues = evaluateEnvironmentHistory( bodies, times );
        for( std::thread& thread : threads )
        {
            thread.join( );
        }

        BOOST_CHECK( originalValues == serialValues );
        for( unsigned int i = 0; i < numberOfThreads; i++ )
        {
            BOOST_CHECK_EQUAL( threadValues.at( i ).size( ), serialValues.size( ) );
            BOOST_CHECK( threadValues.at( i ) == serialValues );
        }
    }
}

//! Test whether a single interpolator can be evaluated by several threads at once (look-up scheme shared by threads)
BOOST_AUTO_TEST_CASE( testConcurrentEvaluationOfSharedInterpolator )
{
    std::map< double, double > dataPoints;
    for( int i = 0; i < 10000; i++ )
    {
        dataPoints[ static_cast< double >( i ) ] = std::sin( 0.01 * i ) + 1.0E-3 * i;
    }
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > interpolator =
            std::make_shared< interpolators::LagrangeInterpolator< double, double > >( dataPoints, 8 );

    const unsigned int numberOfThreads = 8;
    std::vector< std::vector< double > > evaluationTimes( numberOfThreads );
    std::vector< std::vector< double > > serialValues( numberOfThreads );
    for( unsigned int i = 0; i < numberOfThreads; i++ )
    {
        for( int j = 0; j < 200000; j++ )
        {
            evaluationTimes.at( i ).push_back(
                    static_cast< double >( ( static_cast< long long >( j ) * ( 2 * i + 1 ) * 104729 ) % 9990000 ) / 1000.0 + 1.0E-4 * i );
            serialValues.at( i ).push_back( interpolator->interpolate( evaluationTimes.at( i ).back( ) ) );
        }
    }

    std::vector< std::vector< double > > threadValues( numberOfThreads );
    std::vector< std::thread > threads;
    for( unsigned int i = 0; i < numberOfThreads; i++ )
    {
        threads.push_back( std::thread( [ &threadValues, &evaluationTimes, &interpolator, i ]( ) {
            for( const double time : evaluationTimes.at( i ) )
            {
                threadValues.at( i ).push_back( interpolator->interpolate( time ) );
            }
        } ) );
    }
    for( std::thread& thread : threads )
    {
        thread.join( );
    }
    for( unsigned int i = 0; i < numberOfThreads; i++ )
    {
        BOOST_CHECK( threadValues.at( i ) == serialValues.at( i ) );
    }
}

// Propagate vehicles in the test environment (with as many types of environment models as possible) in a number of arcs
std::vector< std::map< double, Eigen::VectorXd > > propagateTestVehicleArcs( const unsigned int numberOfThreads )
{
    using namespace propagators;
    using namespace basic_astrodynamics;

    SystemOfBodies bodies = createTestBodies( );
    const std::vector< std::string > bodiesToPropagate = { "Vehicle", "Vehicle3" };
    const std::vector< std::string > centralBodies = { "Earth", "Earth" };

    // Arc initial states from tabulated ephemerides, which are then replaced by (empty) multi-arc ephemerides
    std::vector< double > arcStartTimes;
    std::vector< Eigen::VectorXd > arcInitialStates;
    for( unsigned int i = 0; i < 5; i++ )
    {
        arcStartTimes.push_back( initialTime + 5400.0 * i );
        Eigen::VectorXd arcInitialState = Eigen::VectorXd( 12 );
        arcInitialState << bodies.at( "Vehicle" )->getEphemeris( )->getCartesianState( arcStartTimes.back( ) ),
                bodies.at( "Vehicle3" )->getEphemeris( )->getCartesianState( arcStartTimes.back( ) );
        arcInitialStates.push_back( arcInitialState );
    }
    for( const std::string& bodyName : bodiesToPropagate )
    {
        addEmptyEphemeris< double, double >( bodies.at( bodyName ), "Earth", "J2000", true, true );
    }

    SelectedAccelerationMap accelerationMap;
    accelerationMap[ "Vehicle" ][ "Earth" ] = { sphericalHarmonicAcceleration( 4, 4 ), aerodynamicAcceleration( ) };
    accelerationMap[ "Vehicle" ][ "Sun" ] = { pointMassGravityAcceleration( ), radiationPressureAcceleration( ) };
    accelerationMap[ "Vehicle" ][ "Moon" ] = { sphericalHarmonicAcceleration( 2, 2 ) };
    accelerationMap[ "Vehicle3" ][ "Earth" ] = { sphericalHarmonicAcceleration( 4, 4 ), radiationPressureAcceleration( ) };
    accelerationMap[ "Vehicle3" ][ "Sun" ] = { pointMassGravityAcceleration( ), radiationPressureAcceleration( ) };
    accelerationMap[ "Vehicle3" ][ "Vehicle3" ] = { thrustAcceleration( { "MainEngine" } ) };

    AccelerationMap accelerationModelMap;
    if( numberOfThreads == 1 )
    {
        accelerationModelMap = createAccelerationModelsMap( bodies, accelerationMap, bodiesToPropagate, centralBodies );
    }

    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > arcPropagationSettingsList;
    for( unsigned int i = 0; i < arcStartTimes.size( ); i++ )
    {
        std::shared_ptr< numerical_integrators::IntegratorSettings<> > integratorSettings =
                std::make_shared< numerical_integrators::IntegratorSettings<> >( numerical_integrators::rungeKutta4, 30.0 );
        std::shared_ptr< PropagationTerminationSettings > terminationSettings =
                std::make_shared< PropagationTimeTerminationSettings >( arcStartTimes.at( i ) + 3600.0 );
        if( numberOfThreads == 1 )
        {
            arcPropagationSettingsList.push_back(
                    std::make_shared< TranslationalStatePropagatorSettings< double > >( centralBodies,
                                                                                        accelerationModelMap,
                                                                                        bodiesToPropagate,
                                                                                        arcInitialStates.at( i ),
                                                                                        arcStartTimes.at( i ),
                                                                                        integratorSettings,
                                                                                        terminationSettings ) );
        }
        else
        {
            arcPropagationSettingsList.push_back(
                    std::make_shared< TranslationalStatePropagatorSettings< double > >( centralBodies,
                                                                                        accelerationMap,
                                                                                        bodiesToPropagate,
                                                                                        arcInitialStates.at( i ),
                                                                                        arcStartTimes.at( i ),
                                                                                        integratorSettings,
                                                                                        terminationSettings ) );
        }
    }
    std::shared_ptr< MultiArcPropagatorSettings< double > > multiArcPropagatorSettings =
            std::make_shared< MultiArcPropagatorSettings< double > >( arcPropagationSettingsList );
    multiArcPropagatorSettings->getOutputSettings( )->setIntegratedResult( true );
    multiArcPropagatorSettings->setNumberOfThreads( numberOfThreads );

    MultiArcDynamicsSimulator<> dynamicsSimulator( bodies, multiArcPropagatorSettings, false );
    dynamicsSimulator.integrateEquationsOfMotion( multiArcPropagatorSettings->getInitialStates( ) );
    return dynamicsSimulator.getMultiArcPropagationResults( )->getConcatenatedEquationsOfMotionResults( );
}

//! Test whether multi-threaded multi-arc propagation in the test environment gives the same results as single-threaded propagation
BOOST_AUTO_TEST_CASE( testParallelPropagationInCopiedEnvironment )
{
    std::vector< std::map< double, Eigen::VectorXd > > singleThreadResults = propagateTestVehicleArcs( 1 );
    BOOST_CHECK_EQUAL( singleThreadResults.size( ), 5 );

    for( unsigned int numberOfThreads : { 2, 3 } )
    {
        std::vector< std::map< double, Eigen::VectorXd > > multiThreadResults = propagateTestVehicleArcs( numberOfThreads );
        BOOST_CHECK_EQUAL( multiThreadResults.size( ), singleThreadResults.size( ) );
        for( unsigned int i = 0; i < singleThreadResults.size( ); i++ )
        {
            BOOST_CHECK_EQUAL( multiThreadResults.at( i ).size( ), singleThreadResults.at( i ).size( ) );
            BOOST_CHECK( singleThreadResults.at( i ).size( ) > 100 );
            for( auto stateIterator : singleThreadResults.at( i ) )
            {
                BOOST_CHECK( !stateIterator.second.hasNaN( ) );
                BOOST_CHECK_SMALL( ( multiThreadResults.at( i ).at( stateIterator.first ) - stateIterator.second )
                                           .segment( 0, 3 )
                                           .cwiseAbs( )
                                           .maxCoeff( ),
                                   1.0E-6 );
                BOOST_CHECK_SMALL( ( multiThreadResults.at( i ).at( stateIterator.first ) - stateIterator.second )
                                           .segment( 6, 3 )
                                           .cwiseAbs( )
                                           .maxCoeff( ),
                                   1.0E-6 );
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
