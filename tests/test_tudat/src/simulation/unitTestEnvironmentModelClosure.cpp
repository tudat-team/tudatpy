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

#include <boost/test/included/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/astro/aerodynamics/aerodynamicCoefficientInterface.h"
#include "tudat/astro/basic_astro/astrodynamicsFunctions.h"
#include "tudat/astro/basic_astro/celestialBodyConstants.h"
#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/astro/electromagnetism/occultationModel.h"
#include "tudat/astro/electromagnetism/radiationPressureAcceleration.h"
#include "tudat/astro/ephemerides/aeordynamicAngleRotationalEphemeris.h"
#include "tudat/astro/ephemerides/tleEphemeris.h"
#include "tudat/astro/gravitation/gravityFieldVariations.h"
#include "tudat/astro/ground_stations/bodyDeformationModel.h"
#include "tudat/astro/reference_frames/referenceFrameTransformations.h"
#include "tudat/astro/system_models/vehicleSystems.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/math/integrators/createNumericalIntegrator.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/createFlightConditions.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/propagation_setup/createAccelerationModels.h"
#include "tudat/simulation/propagation_setup/propagationOutputSettings.h"
#include "tudat/simulation/propagation_setup/propagationSettings.h"
#include "tudat/simulation/propagation_setup/propagationTerminationSettings.h"
#include "tudat/simulation/propagation_setup/singleArcDynamicsSimulator.h"

namespace tudat
{

namespace unit_tests
{

using namespace simulation_setup;

BOOST_AUTO_TEST_SUITE( test_environment_model_closure )

//! Aerodynamic force coefficients, as function of angle of attack (second independent variable)
Eigen::Vector3d getForceCoefficients( const std::vector< double >& independentVariables )
{
    return ( Eigen::Vector3d( ) << 1.2, 0.0, 0.3 + 0.5 * independentVariables.at( 1 ) ).finished( );
}

//! Aerodynamic moment coefficients, as function of angle of attack (second independent variable), with zero pitch moment at 0.2 rad
Eigen::Vector3d getMomentCoefficients( const std::vector< double >& independentVariables )
{
    return ( Eigen::Vector3d( ) << 0.0, 0.02 - 0.1 * independentVariables.at( 1 ), 0.0 ).finished( );
}

//! Test that adding flight conditions applies the closures that wait for them, and does not modify models that are closed
BOOST_AUTO_TEST_CASE( testAddFlightConditionsClosure )
{
    spice_interface::loadStandardSpiceKernels( );

    // Create Earth with atmosphere, and two vehicles of which the aerodynamic moment contribution requires flight conditions. One
    // vehicle has an aerodynamic angle based rotation model, the other a trimmed rotation model (which also requires flight
    // conditions)
    BodyListSettings bodySettings = getDefaultBodySettings( { "Earth" }, "Earth", "ECLIPJ2000" );
    bodySettings.at( "Earth" )->atmosphereSettings = exponentialAtmosphereSettings( 7.2E3, 1.225 );

    const std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > independentVariables = {
        aerodynamics::mach_number_dependent, aerodynamics::angle_of_attack_dependent
    };
    const Eigen::Vector3d momentReferencePoint = ( Eigen::Vector3d( ) << 0.1, 0.0, 0.0 ).finished( );
    for( const std::string& vehicleName : std::vector< std::string >( { "Vehicle", "TrimmedVehicle" } ) )
    {
        bodySettings.addSettings( vehicleName );
        bodySettings.at( vehicleName )->constantMass = 1000.0;
        bodySettings.at( vehicleName )->aerodynamicCoefficientSettings =
                customAerodynamicCoefficientSettings( &getForceCoefficients,
                                                      &getMomentCoefficients,
                                                      2.0,
                                                      4.0,
                                                      independentVariables,
                                                      aerodynamics::negative_aerodynamic_frame_coefficients,
                                                      aerodynamics::body_fixed_frame_coefficients,
                                                      momentReferencePoint );
    }
    bodySettings.at( "Vehicle" )->rotationModelSettings = aerodynamicAngleRotationSettings( "Earth", "ECLIPJ2000", "VehicleFixed" );
    bodySettings.at( "TrimmedVehicle" )->rotationModelSettings = pitchTrimRotationSettings( "Earth", "ECLIPJ2000", "TrimmedVehicleFixed" );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // Closure of aerodynamic angle based rotation model is applied on creation, closure of moment contribution and trim conditions
    // wait for flight conditions
    std::shared_ptr< ephemerides::AerodynamicAngleRotationalEphemeris > vehicleRotationModel =
            std::dynamic_pointer_cast< ephemerides::AerodynamicAngleRotationalEphemeris >(
                    bodies.at( "Vehicle" )->getRotationalEphemeris( ) );
    std::shared_ptr< ephemerides::AerodynamicAngleRotationalEphemeris > trimmedVehicleRotationModel =
            std::dynamic_pointer_cast< ephemerides::AerodynamicAngleRotationalEphemeris >(
                    bodies.at( "TrimmedVehicle" )->getRotationalEphemeris( ) );
    std::shared_ptr< aerodynamics::AerodynamicCoefficientInterface > vehicleCoefficientInterface =
            bodies.at( "Vehicle" )->getAerodynamicCoefficientInterface( );
    std::shared_ptr< aerodynamics::AerodynamicCoefficientInterface > trimmedVehicleCoefficientInterface =
            bodies.at( "TrimmedVehicle" )->getAerodynamicCoefficientInterface( );

    BOOST_CHECK( vehicleRotationModel->getAerodynamicAngleCalculator( ) != nullptr );
    BOOST_CHECK( vehicleCoefficientInterface->getMomentContributionInterface( ) == nullptr );
    BOOST_CHECK( trimmedVehicleRotationModel->getAerodynamicAngleFunction( ) == nullptr );
    BOOST_CHECK( trimmedVehicleCoefficientInterface->getMomentContributionInterface( ) == nullptr );

    // Set moment contribution of first vehicle by user, so that its closure is applied before adding flight conditions
    std::shared_ptr< aerodynamics::AerodynamicMomentContributionInterface > userMomentContributionInterface =
            std::make_shared< aerodynamics::AerodynamicMomentContributionInterface >( []( ) { return Eigen::Matrix3d::Identity( ); },
                                                                                      []( ) { return Eigen::Matrix3d::Identity( ); },
                                                                                      []( ) { return Eigen::Vector3d::Zero( ); },
                                                                                      1.0 );
    vehicleCoefficientInterface->setMomentContributionInterface( userMomentContributionInterface );
    std::shared_ptr< reference_frames::AerodynamicAngleCalculator > vehicleAngleCalculator =
            vehicleRotationModel->getAerodynamicAngleCalculator( );

    addFlightConditions( bodies, "Vehicle", "Earth" );
    addFlightConditions( bodies, "TrimmedVehicle", "Earth" );

    // Models of first vehicle of which closure was applied are not modified; flight conditions use angle calculator of rotation model
    BOOST_CHECK( vehicleRotationModel->getAerodynamicAngleCalculator( ) == vehicleAngleCalculator );
    BOOST_CHECK( bodies.at( "Vehicle" )->getFlightConditions( )->getAerodynamicAngleCalculator( ) == vehicleAngleCalculator );
    BOOST_CHECK( vehicleCoefficientInterface->getMomentContributionInterface( ) == userMomentContributionInterface );

    // Closures of second vehicle that waited for flight conditions are applied
    BOOST_CHECK( trimmedVehicleRotationModel->getAerodynamicAngleFunction( ) != nullptr );
    BOOST_CHECK( trimmedVehicleCoefficientInterface->getMomentContributionInterface( ) != nullptr );
}

//! Direction of the body-fixed x-axis in the TNW frame (along the velocity)
Eigen::Vector3d getAlongTrackDirection( const double )
{
    return Eigen::Vector3d::UnitX( );
}

//! Test that a direction-based rotation model created through body settings, with the direction defined in the TNW frame w.r.t.
//! Earth, aligns the body-fixed x-axis with the velocity w.r.t. Earth
BOOST_AUTO_TEST_CASE( testDirectionBasedRotationClosure )
{
    spice_interface::loadStandardSpiceKernels( );

    BodyListSettings bodySettings = getDefaultBodySettings( { "Earth" }, "Earth", "ECLIPJ2000" );
    bodySettings.addSettings( "Vehicle" );
    bodySettings.at( "Vehicle" )->ephemerisSettings = keplerEphemerisSettings(
            ( Eigen::Vector6d( ) << 7.0E6, 0.1, 0.5, 1.0, 2.0, 0.3 ).finished( ), 0.0, 3.986004418E14, "Earth", "ECLIPJ2000" );
    bodySettings.at( "Vehicle" )->rotationModelSettings = std::make_shared< BodyFixedDirectionBasedRotationSettings >(
            &getAlongTrackDirection,
            "ECLIPJ2000",
            "VehicleFixed",
            nullptr,
            std::make_pair( ephemerides::tnw_satellite_based_frame, std::string( "Earth" ) ) );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    for( const double time : { 0.0, 1.0E3, 4.0E3 } )
    {
        bodies.at( "Earth" )->setStateFromEphemeris( time );
        bodies.at( "Vehicle" )->setStateFromEphemeris( time );
        const Eigen::Vector3d bodyFixedXAxis =
                bodies.at( "Vehicle" )->getRotationalEphemeris( )->getRotationToBaseFrame( time ) * Eigen::Vector3d::UnitX( );
        const Eigen::Vector3d unitVelocity =
                ( bodies.at( "Vehicle" )->getState( ) - bodies.at( "Earth" )->getState( ) ).segment( 3, 3 ).normalized( );
        for( int i = 0; i < 3; i++ )
        {
            BOOST_CHECK_SMALL( bodyFixedXAxis( i ) - unitVelocity( i ), 1.0E-12 );
        }
    }
}

//! Constant angle of attack, sideslip angle and bank angle (in radians)
Eigen::Vector3d getConstantAerodynamicAngles( const double )
{
    return ( Eigen::Vector3d( ) << 0.1, -0.02, 0.3 ).finished( );
}

//! Function to get settings of Earth with an exponential atmosphere, and of a vehicle with mass 1000 kg, reference length 2 meters,
//! reference area 4 meters^2 and aerodynamic coefficients given by getForceCoefficients and getMomentCoefficients
BodyListSettings getVehicleInAtmosphereBodySettings( const Eigen::Vector3d& momentReferencePoint = Eigen::Vector3d::Constant( TUDAT_NAN ) )
{
    BodyListSettings bodySettings = getDefaultBodySettings( { "Earth" }, "Earth", "ECLIPJ2000" );
    bodySettings.at( "Earth" )->atmosphereSettings = exponentialAtmosphereSettings( "Earth" );
    bodySettings.addSettings( "Vehicle" );
    bodySettings.at( "Vehicle" )->constantMass = 1000.0;
    bodySettings.at( "Vehicle" )->aerodynamicCoefficientSettings =
            customAerodynamicCoefficientSettings( &getForceCoefficients,
                                                  &getMomentCoefficients,
                                                  2.0,
                                                  4.0,
                                                  { aerodynamics::mach_number_dependent, aerodynamics::angle_of_attack_dependent },
                                                  aerodynamics::negative_aerodynamic_frame_coefficients,
                                                  aerodynamics::body_fixed_frame_coefficients,
                                                  momentReferencePoint );
    return bodySettings;
}

//! Function to propagate the vehicle around Earth for 10 s (RK4, 1 s step) with point mass gravity and aerodynamic acceleration,
//! and to return the history of the given dependent variables
std::map< double, Eigen::VectorXd > propagateVehicleInAtmosphere(
        const SystemOfBodies& bodies,
        const std::vector< std::shared_ptr< propagators::SingleDependentVariableSaveSettings > >& dependentVariables )
{
    // Creating the aerodynamic acceleration adds the flight conditions, which applies the closures of the models of the vehicle that
    // wait for them
    SelectedAccelerationMap accelerationSettings;
    accelerationSettings[ "Vehicle" ][ "Earth" ] = { pointMassGravityAcceleration( ), aerodynamicAcceleration( ) };
    const Eigen::Vector6d initialState =
            ( Eigen::Vector6d( ) << spice_interface::getAverageRadius( "Earth" ) + 120.0E3, 0.0, 0.0, 0.0, 7.8E3, 0.0 ).finished( );
    propagators::SingleArcDynamicsSimulator<> dynamicsSimulator(
            bodies,
            propagators::translationalStatePropagatorSettings< double >(
                    { "Earth" },
                    createAccelerationModelsMap( bodies, accelerationSettings, { "Vehicle" }, { "Earth" } ),
                    { "Vehicle" },
                    initialState,
                    0.0,
                    numerical_integrators::rungeKutta4Settings< double >( 1.0 ),
                    propagators::propagationTimeTerminationSettings( 10.0 ),
                    propagators::cowell,
                    dependentVariables ) );
    return dynamicsSimulator.getDependentVariableHistory( );
}

//! Test that an aerodynamic angle based rotation model created through body settings imposes the aerodynamic angles of its angle
//! function during propagation
BOOST_AUTO_TEST_CASE( testAerodynamicAngleRotationClosure )
{
    spice_interface::loadStandardSpiceKernels( );

    BodyListSettings bodySettings = getVehicleInAtmosphereBodySettings( );
    bodySettings.at( "Vehicle" )->rotationModelSettings =
            aerodynamicAngleRotationSettings( "Earth", "ECLIPJ2000", "VehicleFixed", &getConstantAerodynamicAngles );
    for( auto it : propagateVehicleInAtmosphere( createSystemOfBodies( bodySettings ),
                                                 { propagators::angleOfAttackDependentVariable( "Vehicle", "Earth" ),
                                                   propagators::sideslipAngleDependentVariable( "Vehicle", "Earth" ),
                                                   propagators::bankAngleDependentVariable( "Vehicle", "Earth" ) } ) )
    {
        for( int i = 0; i < 3; i++ )
        {
            BOOST_CHECK_SMALL( it.second( i ) - getConstantAerodynamicAngles( it.first )( i ), 1.0E-14 );
        }
    }
}

//! Test that a pitch trim rotation model created through body settings trims the vehicle during propagation: the pitch moment
//! coefficient of getMomentCoefficients is zero at an angle of attack of 0.2 rad
BOOST_AUTO_TEST_CASE( testPitchTrimRotationClosure )
{
    spice_interface::loadStandardSpiceKernels( );

    BodyListSettings bodySettings = getVehicleInAtmosphereBodySettings( );
    bodySettings.at( "Vehicle" )->rotationModelSettings = pitchTrimRotationSettings( "Earth", "ECLIPJ2000", "VehicleFixed" );
    for( auto it : propagateVehicleInAtmosphere( createSystemOfBodies( bodySettings ),
                                                 { propagators::angleOfAttackDependentVariable( "Vehicle", "Earth" ) } ) )
    {
        // Tolerance above the absolute tolerance (1.0E-13 rad) of the root finder of aerodynamics::TrimOrientationCalculator
        BOOST_CHECK_SMALL( it.second( 0 ) - 0.2, 1.0E-12 );
    }
}

//! Test that the contribution of the aerodynamic force to the moment coefficients, with a center of mass from rigid body properties,
//! is added for coefficients created through body settings
BOOST_AUTO_TEST_CASE( testAerodynamicMomentContributionClosure )
{
    spice_interface::loadStandardSpiceKernels( );

    const Eigen::Vector3d momentReferencePoint( 0.1, 0.0, 0.0 );
    const Eigen::Vector3d centerOfMass( 0.0, 0.05, -0.02 );
    BodyListSettings bodySettings = getVehicleInAtmosphereBodySettings( momentReferencePoint );
    bodySettings.at( "Vehicle" )->rigidBodyPropertiesSettings = constantRigidBodyPropertiesSettings( 1000.0, centerOfMass );
    bodySettings.at( "Vehicle" )->rotationModelSettings =
            aerodynamicAngleRotationSettings( "Earth", "ECLIPJ2000", "VehicleFixed", &getConstantAerodynamicAngles );

    // Coefficients at the imposed angle of attack (Mach number is not used by the coefficient functions). The force coefficients are
    // in the negative aerodynamic frame, the moment coefficients in the body-fixed frame, so that the force coefficients in the
    // body-fixed frame are -R_AB C_F (as in testCombinedAerodynamicForceAndMoment in unitTestAerodynamicMomentAndAerodynamicForce.cpp)
    const Eigen::Vector3d aerodynamicAngles = getConstantAerodynamicAngles( 0.0 );
    const std::vector< double > independentVariables = { 0.0, aerodynamicAngles( 0 ) };
    const Eigen::Vector3d bodyFixedForceCoefficients = -reference_frames::getAirspeedBasedAerodynamicToBodyFrameTransformationMatrix(
                                                               aerodynamicAngles( 0 ), aerodynamicAngles( 1 ) ) *
            getForceCoefficients( independentVariables );
    const Eigen::Vector3d expectedMomentCoefficients = getMomentCoefficients( independentVariables ) +
            ( momentReferencePoint - centerOfMass ).cross( bodyFixedForceCoefficients ) / 2.0;

    for( auto it : propagateVehicleInAtmosphere( createSystemOfBodies( bodySettings ),
                                                 { propagators::aerodynamicMomentCoefficientDependentVariable( "Vehicle" ),
                                                   propagators::aerodynamicForceCoefficientDependentVariable( "Vehicle" ) } ) )
    {
        for( int i = 0; i < 3; i++ )
        {
            BOOST_CHECK_SMALL( it.second( i ) - expectedMomentCoefficients( i ), 1.0E-14 );
            BOOST_CHECK_SMALL( it.second( i + 3 ) - getForceCoefficients( independentVariables )( i ), 1.0E-14 );
        }
    }
}

//! Test that mode-coupled tidal gravity field variations created through body settings, with only diagonal Love numbers, are equal
//! to basic tidal gravity field variations with the same Love numbers
BOOST_AUTO_TEST_CASE( testModeCoupledTidalGravityFieldVariationClosure )
{
    spice_interface::loadStandardSpiceKernels( );

    const double loveNumber = 0.3;
    std::map< std::pair< int, int >, std::map< std::pair< int, int >, double > > modeCoupledLoveNumbers;
    for( int order = 0; order <= 2; order++ )
    {
        modeCoupledLoveNumbers[ std::make_pair( 2, order ) ][ std::make_pair( 2, order ) ] = loveNumber;
    }

    // Basic (first) and mode-coupled (second) tidal variations of the gravity field of Earth due to the Moon
    BodyListSettings bodySettings = getDefaultBodySettings( { "Earth", "Moon" }, "Earth", "ECLIPJ2000" );
    bodySettings.at( "Earth" )->gravityFieldVariationSettings = {
        fixedSingleDegreeLoveNumberGravityFieldVariationSettings( "Moon", loveNumber, 2 ),
        modeCoupledSolidBodyGravityFieldVariationSettings( { "Moon" }, modeCoupledLoveNumbers )
    };
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    std::vector< std::shared_ptr< gravitation::GravityFieldVariations > > variations =
            bodies.at( "Earth" )->getGravityFieldVariationSet( )->getVariationObjects( );

    for( const double time : { 0.0, 1.0E5, 1.0E7 } )
    {
        bodies.at( "Earth" )->setStateFromEphemeris( time );
        bodies.at( "Earth" )->setCurrentRotationalStateToLocalFrameFromEphemeris( time );
        bodies.at( "Moon" )->setStateFromEphemeris( time );

        const std::pair< Eigen::MatrixXd, Eigen::MatrixXd > basicCorrections =
                variations.at( 0 )->calculateSphericalHarmonicsCorrections( time );
        const std::pair< Eigen::MatrixXd, Eigen::MatrixXd > modeCoupledCorrections =
                variations.at( 1 )->calculateSphericalHarmonicsCorrections( time );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( modeCoupledCorrections.first, basicCorrections.first, 1.0E-12 );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( modeCoupledCorrections.second, basicCorrections.second, 1.0E-12 );
    }
}

//! Test that panelled aerodynamic coefficients created through body settings, with a Newton gas-surface interaction model, are
//! equal to the analytical Newton result
BOOST_AUTO_TEST_CASE( testPanelledAerodynamicCoefficientsClosure )
{
    // Single panel with an area of 0.5 meters^2 and a surface normal along the body-fixed x-axis, and a reference area of
    // 1.0 meters^2 (as in test_panelled_coefficients in unitTestAerodynamicMomentAndAerodynamicForce.cpp)
    BodyListSettings bodySettings;
    bodySettings.addSettings( "Vehicle" );
    bodySettings.at( "Vehicle" )->rotationModelSettings =
            constantRotationModelSettings( "ECLIPJ2000", "VehicleFixed", Eigen::Matrix3d::Identity( ) );
    bodySettings.at( "Vehicle" )->bodyExteriorPanelSettings_ = fullPanelledBodySettings( { bodyPanelSettings(
            std::make_shared< FrameFixedBodyPanelGeometrySettings >(
                    Eigen::Vector3d::Zero( ), Eigen::Vector3d::UnitY( ), Eigen::Vector3d::UnitZ( ), Eigen::Vector3d::UnitX( ), 0.5 ),
            nullptr ) } );
    bodySettings.at( "Vehicle" )->aerodynamicCoefficientSettings = panelledAerodynamicCoefficientSettings( aerodynamics::newton, 1.0 );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // Set body-fixed orientation of panels (done by VehicleSystems::updatePartOrientations during propagation)
    bodies.at( "Vehicle" )->getVehicleSystems( )->updatePartOrientations( 0.0 );
    std::shared_ptr< aerodynamics::AerodynamicCoefficientInterface > coefficientInterface =
            bodies.at( "Vehicle" )->getAerodynamicCoefficientInterface( );
    for( const double angleOfAttack : { 0.0, 0.3 } )
    {
        for( const double sideslipAngle : { 0.0, 0.2 } )
        {
            coefficientInterface->updateFullCurrentCoefficients( { angleOfAttack, sideslipAngle } );

            // Newton: C_p = 2 cos^2( delta ), with cos( delta ) the inner product of the panel normal and the reversed incoming flow
            // direction, and C_F = -C_p n A_panel / A_ref
            const double cosineDelta = std::cos( angleOfAttack ) * std::cos( sideslipAngle );
            const Eigen::Vector3d expectedForceCoefficients = -2.0 * cosineDelta * cosineDelta * Eigen::Vector3d::UnitX( ) * 0.5 / 1.0;
            for( int i = 0; i < 3; i++ )
            {
                BOOST_CHECK_SMALL( coefficientInterface->getCurrentForceCoefficients( )( i ) - expectedForceCoefficients( i ),
                                   std::numeric_limits< double >::epsilon( ) );
            }
        }
    }
}

//! Test that a paneled radiation source of Earth created through body settings gives the albedo and infrared accelerations of
//! LAGEOS of Knocke (1988), with the set-up and checks of
//! testRadiationPressureAcceleration_DynamicallyPanelunitestradedSource_CannonballTarget_LAGEOS in
//! unitTestRadiationPressureAcceleration.cpp, of which the source is created directly
BOOST_AUTO_TEST_CASE( testPaneledRadiationSourceClosure )
{
    spice_interface::loadStandardSpiceKernels( );

    // Earth source with Knocke albedo and emissivity distributions (same settings as the LAGEOS test), and cannonball target
    BodyListSettings bodySettings = getDefaultBodySettings( { "Earth", "Sun" }, "Earth", "J2000" );
    bodySettings.at( "Earth" )->radiationSourceModelSettings = getKnockeEarthRadiationPressureSettings( );
    bodySettings.addSettings( "Vehicle" );
    bodySettings.at( "Vehicle" )->radiationPressureTargetModelSettings = cannonballRadiationPressureTargetModelSettings( 0.28483, 1.13 );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    std::shared_ptr< electromagnetism::PaneledRadiationSourceModel > sourceModel =
            std::dynamic_pointer_cast< electromagnetism::PaneledRadiationSourceModel >( bodies.at( "Earth" )->getRadiationSourceModel( ) );

    ephemerides::TleEphemeris ephemeris(
            "Earth",
            "J2000",
            std::make_shared< ephemerides::Tle >( "1 08820U 76039  A 77047.52561960  .00000002 +00000-0 +00000-0 0  9994\n"
                                                  "2 08820 109.8332 127.3884 0044194 201.3006 158.6132 06.38663945018402" ) );
    const double startTime = spice_interface::convertDateStringToEphemerisTime( "1977-02-22" );
    const double orbitalPeriod = basic_astrodynamics::computeKeplerOrbitalPeriod(
            orbital_element_conversions::convertCartesianToKeplerianElements(
                    ephemeris.getCartesianState( startTime ),
                    celestial_body_constants::EARTH_GRAVITATIONAL_PARAMETER )( orbital_element_conversions::semiMajorAxisIndex ),
            celestial_body_constants::EARTH_GRAVITATIONAL_PARAMETER,
            0 );

    for( int i = 0; i <= 30; ++i )
    {
        const double t = startTime + i / 10.0 * orbitalPeriod;

        bodies.at( "Earth" )->setStateFromEphemeris( t );
        bodies.at( "Earth" )->setCurrentRotationToLocalFrameFromEphemeris( t );
        bodies.at( "Sun" )->setStateFromEphemeris( t );

        const Eigen::Quaterniond earthToGlobalRotation = bodies.at( "Earth" )->getCurrentRotationToGlobalFrame( );
        const Eigen::Vector3d position = ephemeris.getCartesianPosition( t );
        const Eigen::Vector3d velocity = ephemeris.getCartesianVelocity( t );

        electromagnetism::PaneledSourceRadiationPressureAcceleration accelerationModel(
                sourceModel,
                []( ) { return Eigen::Vector3d::Zero( ); },
                [ = ]( ) { return earthToGlobalRotation; },
                bodies.at( "Vehicle" )->getRadiationPressureTargetModels( ).at( 0 ),
                [ = ]( ) { return position; },
                []( ) { return Eigen::Quaterniond::Identity( ); },
                []( ) { return 406.9; },
                std::make_shared< electromagnetism::NoOccultingBodyOccultationModel >( ) );
        sourceModel->updateMembers( t );
        accelerationModel.updateMembers( t );
        const Eigen::Vector3d acceleration = accelerationModel.getAcceleration( );

        const Eigen::Vector3d radialUnit = position.normalized( );
        const Eigen::Vector3d alongTrackUnit = ( velocity - radialUnit * velocity.dot( radialUnit ) ).normalized( );
        const Eigen::Vector3d crossTrackUnit = radialUnit.cross( alongTrackUnit );

        // Compare with Knocke (1988), Fig. 4 (ranges copied from unitTestRadiationPressureAcceleration.cpp, lines 1470-1480)
        BOOST_CHECK_GT( acceleration.dot( radialUnit ), 1.5e-10 );
        BOOST_CHECK_LT( acceleration.dot( radialUnit ), 3.9e-10 );

        BOOST_CHECK_GT( acceleration.dot( alongTrackUnit ), -4.0e-11 );
        BOOST_CHECK_LT( acceleration.dot( alongTrackUnit ), 4.0e-11 );

        BOOST_CHECK_GT( acceleration.dot( crossTrackUnit ), -3.0e-11 );
        BOOST_CHECK_LT( acceleration.dot( crossTrackUnit ), 1.0e-11 );
    }
}

//! Test that a coma atmosphere created through body settings gives the reference densities of the coma model
BOOST_AUTO_TEST_CASE( testComaAtmosphereClosure )
{
    // Stokes dataset (radii from 4 to 10 km, solar longitudes from 0 to 360 degrees) and molecular weight (in kg/mol) as in
    // test_coma_model_density_validation_from_python in unitTestComaModel.cpp
    std::vector< double > radii;
    for( int i = 0; i < 100; ++i )
    {
        radii.push_back( 4000.0 + i * ( 10000.0 - 4000.0 ) / 99.0 );
    }
    std::vector< double > solarLongitudes;
    for( int i = 0; i < 37; ++i )
    {
        solarLongitudes.push_back( i * 360.0 / 36.0 );
    }
    const ComaModelFileProcessor processor(
            std::vector< std::string >( { paths::getTudatTestDataPath( ) + "coma/density/polynomial/input_poly_coef_test_file.txt" } ) );

    BodyListSettings bodySettings;
    bodySettings.addSettings( "Sun" );
    bodySettings.addSettings( "Comet" );
    bodySettings.at( "Comet" )->rotationModelSettings =
            constantRotationModelSettings( "ECLIPJ2000", "CometFixed", Eigen::Matrix3d::Identity( ) );
    bodySettings.at( "Comet" )->atmosphereSettings =
            comaSettings( processor.createSHDataset( radii, solarLongitudes, 10, 10 ), 0.018, 10, 10 );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // Time [s], radial distance [m], latitude [rad], longitude [rad], solar longitude [rad] and density [kg/m^3], copied from
    // tests/test_tudat/data/coma/reference_values.txt, lines 2 and 700
    const std::vector< std::vector< double > > referencePoints = {
        { 490708800.000000, 10000.000000, 0.1548961939, -2.3127640117, 5.9867226630, 4.429418284212357e-10 },
        { 490729740.000000, 11285.333314, 0.3405541931, 1.4972023270, 2.9119497845, 6.768902722269597e-11 }
    };
    for( const std::vector< double >& point : referencePoints )
    {
        // Comet at the origin, Sun at 1.0E11 meters in the comet equatorial plane at the solar longitude of the reference point (as
        // in test_coma_model_density_validation_from_python)
        bodies.at( "Sun" )->setState(
                ( Eigen::Vector6d( ) << 1.0E11 * std::cos( point.at( 4 ) ), 1.0E11 * std::sin( point.at( 4 ) ), 0.0, 0.0, 0.0, 0.0 )
                        .finished( ) );
        bodies.at( "Comet" )->setState( Eigen::Vector6d::Zero( ) );
        bodies.at( "Comet" )->setCurrentRotationalStateToLocalFrameFromEphemeris( point.at( 0 ) );

        BOOST_CHECK_CLOSE_FRACTION(
                bodies.at( "Comet" )->getAtmosphereModel( )->getDensity( point.at( 1 ), point.at( 3 ), point.at( 2 ), point.at( 0 ) ),
                point.at( 5 ),
                1.0E-8 );
    }
}

//! Test that a basic tidal deformation model created through body settings gives the degree-2 displacement of IERS Conventions
//! 2010, eq. 7.5
BOOST_AUTO_TEST_CASE( testBasicTidalBodyDeformationClosure )
{
    spice_interface::loadStandardSpiceKernels( );

    const double loveNumber = 0.6078;
    const double shidaNumber = 0.0847;
    const double referenceRadius = 6378137.0;
    BodyListSettings bodySettings = getDefaultBodySettings( { "Earth", "Moon" }, "Earth", "ECLIPJ2000" );
    bodySettings.at( "Earth" )->bodyDeformationSettings.push_back(
            degreeTwoBasicTidalBodyShapeDeformation( { "Moon" }, loveNumber, shidaNumber, referenceRadius ) );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    const Eigen::Vector3d stationPosition( 3.9E6, 3.0E5, 5.0E6 );
    const Eigen::Vector3d stationUnitVector = stationPosition.normalized( );
    const double gravitationalParameterRatio = bodies.at( "Moon" )->getGravityFieldModel( )->getGravitationalParameter( ) /
            bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
    for( const double time : { 0.0, 1.0E5, 1.0E7 } )
    {
        // Position of the Moon w.r.t. Earth, in the body-fixed frame of Earth
        const Eigen::Vector3d inertialMoonPosition = bodies.at( "Moon" )->getPositionInBaseFrameFromEphemeris( time ) -
                bodies.at( "Earth" )->getPositionInBaseFrameFromEphemeris( time );
        const Eigen::Vector3d moonPosition =
                bodies.at( "Earth" )->getRotationalEphemeris( )->getRotationToTargetFrame( time ) * inertialMoonPosition;
        const Eigen::Vector3d moonUnitVector = moonPosition.normalized( );
        const double innerProduct = moonUnitVector.dot( stationUnitVector );

        // IERS Conventions 2010, eq. 7.5
        const Eigen::Vector3d expectedDisplacement = gravitationalParameterRatio * std::pow( referenceRadius, 4 ) /
                std::pow( moonPosition.norm( ), 3 ) *
                ( loveNumber * stationUnitVector * ( 1.5 * innerProduct * innerProduct - 0.5 ) +
                  3.0 * shidaNumber * innerProduct * ( moonUnitVector - innerProduct * stationUnitVector ) );
        const Eigen::Vector3d displacement =
                bodies.at( "Earth" )->getBodyDeformationModels( ).at( 0 )->calculateDisplacement( time, stationPosition );
        for( int i = 0; i < 3; i++ )
        {
            BOOST_CHECK_SMALL( displacement( i ) - expectedDisplacement( i ), 1.0E-12 * expectedDisplacement.norm( ) );
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
