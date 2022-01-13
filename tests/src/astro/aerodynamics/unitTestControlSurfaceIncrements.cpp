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

#include <boost/array.hpp>
#include <boost/make_shared.hpp>
#include <memory>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/astro/aerodynamics/hypersonicLocalInclinationAnalysis.h"
#include "tudat/astro/aerodynamics/customAerodynamicCoefficientInterface.h"
#include "tudat/basics/basicTypedefs.h"
#include "tudat/math/geometric/capsule.h"
#include "tudat/math/geometric/sphereSegment.h"
#include "tudat/simulation/simulation.h"

#include "tudat/astro/aerodynamics/testApolloCapsuleCoefficients.h"

namespace tudat
{

namespace unit_tests
{

//! Dummy class used to test the use of control surface deflections in numerical propagation.
//! In this class, a single control surface named "TestSurface" is set to a time-dependent deflection at each time step,
//! as is the angle of attack. The update is performed by the AerodynamicAngleCalculator, as linked by this class'
//! constructor
class DummyGuidanceSystem
{
public:
    DummyGuidanceSystem(
            const std::function< void( const std::string&, const double ) > controlSurfaceFunction,
            const std::shared_ptr< reference_frames::AerodynamicAngleCalculator > angleCalculator ):
        controlSurfaceFunction_( controlSurfaceFunction ), angleCalculator_( angleCalculator ),
        currentAngleOfAttack_( 0.0 ), currentSurfaceDeflection_( 0.0 )
    {
        angleCalculator_->setOrientationAngleFunctions(
                    std::bind( &DummyGuidanceSystem::getCurrentAngleOfAttack, this ),
                    std::function< double( ) >( ),
                    std::function< double( ) >( ),
                    std::bind( &DummyGuidanceSystem::updateGuidance, this, std::placeholders::_1 ) );

        controlSurfaceFunction_( "TestSurface", 0.2 );
    }

    ~DummyGuidanceSystem( ){ }

    void updateGuidance( const double currentTime )
    {
        currentAngleOfAttack_ = 0.3 * ( 1.0 - currentTime / 1000.0 );
        currentSurfaceDeflection_ = -0.02 + 0.04 * currentTime / 1000.0;
        controlSurfaceFunction_( "TestSurface", currentSurfaceDeflection_ );
    }

    double getCurrentAngleOfAttack( )
    {
        return currentAngleOfAttack_;
    }

    double getCurrentSurfaceDeflection( )
    {
        return currentSurfaceDeflection_;
    }

private:
    std::function< void( const std::string&, const double ) > controlSurfaceFunction_;

    std::shared_ptr< reference_frames::AerodynamicAngleCalculator > angleCalculator_;

    double currentAngleOfAttack_;

    double currentSurfaceDeflection_;

};

using Eigen::Vector6d;
using mathematical_constants::PI;

BOOST_AUTO_TEST_SUITE( test_control_surface_increments )

//! Function to return dummy control increments as a function of 2 independnt variables
Eigen::Vector6d dummyControlIncrements(
        const std::vector< double > independentVariables )
{
    Eigen::Vector6d randomControlIncrements =
            ( Eigen::Vector6d( ) << 1.0, -3.5, 2.1, 0.4, -0.75, 1.3 ).finished( );
    for( unsigned int i = 0; i < 6; i++ )
    {
        randomControlIncrements( i ) *= (
                    0.01 * independentVariables.at( 0 ) + static_cast< double >( i ) * 0.005 * independentVariables.at( 1 ) );
    }

    BOOST_CHECK_EQUAL( independentVariables.size( ), 2 );

    return randomControlIncrements;
}

//! Test update and retrieval of control surface aerodynamic coefficient increments, outside of the numerical propagation
BOOST_AUTO_TEST_CASE( testControlSurfaceIncrementInterface )
{
    // Create aerodynamic coefficient interface without control increments.
    std::shared_ptr< AerodynamicCoefficientInterface > coefficientInterfaceWithoutIncrements =
            getApolloCoefficientInterface( );

    // Create aerodynamic coefficient interface with control increments.
    std::shared_ptr< AerodynamicCoefficientInterface > coefficientInterfaceWithIncrements =
            getApolloCoefficientInterface( );
    std::shared_ptr< ControlSurfaceIncrementAerodynamicInterface > controlSurfaceInterface =
            std::make_shared< CustomControlSurfaceIncrementAerodynamicInterface >(
                &dummyControlIncrements,
                std::vector< AerodynamicCoefficientsIndependentVariables >{ angle_of_attack_dependent, control_surface_deflection_dependent } );
    std::map< std::string, std::shared_ptr< ControlSurfaceIncrementAerodynamicInterface >  > controlSurfaceList;
    controlSurfaceList[ "TestSurface" ] = controlSurfaceInterface;
    coefficientInterfaceWithIncrements->setControlSurfaceIncrements( controlSurfaceList );

    // Define values of independent variables of body aerodynamics
    std::vector< double > independentVariables;
    independentVariables.push_back( 10.0 );
    independentVariables.push_back( 0.1 );
    independentVariables.push_back( -0.01 );

    // Define values of independent variables of control surface aerodynamics
    std::map< std::string, std::vector< double > > controlSurfaceIndependentVariables;
    controlSurfaceIndependentVariables[ "TestSurface" ].push_back( 0.1 );
    controlSurfaceIndependentVariables[ "TestSurface" ].push_back( 0.0 );

    // Declare test variables.
    Eigen::Vector3d forceWithIncrement, forceWithoutIncrement;
    Eigen::Vector3d momentWithIncrement, momentWithoutIncrement;

    Eigen::Vector6d manualControlIncrements;

    // Test coefficient interfaces for range of independent variables.
    for( double angleOfAttack = -0.4; angleOfAttack < 0.4; angleOfAttack += 0.02 )
    {
        for( double deflectionAngle = -0.05; deflectionAngle < 0.05; deflectionAngle += 0.001 )
        {
            // Set indepdnent variables.
            controlSurfaceIndependentVariables[ "TestSurface" ][ 0 ] = angleOfAttack;
            controlSurfaceIndependentVariables[ "TestSurface" ][ 1 ] = deflectionAngle;
            independentVariables[ 1 ] = angleOfAttack;

            // Update coefficients
            coefficientInterfaceWithoutIncrements->updateFullCurrentCoefficients(
                        independentVariables );
            coefficientInterfaceWithIncrements->updateFullCurrentCoefficients(
                        independentVariables, controlSurfaceIndependentVariables );

            // Retrieve coefficients.
            forceWithIncrement = coefficientInterfaceWithIncrements->getCurrentForceCoefficients( );
            forceWithoutIncrement = coefficientInterfaceWithoutIncrements->getCurrentForceCoefficients( );

            momentWithIncrement = coefficientInterfaceWithIncrements->getCurrentMomentCoefficients( );
            momentWithoutIncrement = coefficientInterfaceWithoutIncrements->getCurrentMomentCoefficients( );

            // Test coefficients
            manualControlIncrements = dummyControlIncrements( controlSurfaceIndependentVariables[ "TestSurface" ] );

            for( unsigned int i = 0; i < 3; i++ )
            {
                BOOST_CHECK_SMALL( std::fabs( forceWithIncrement( i ) -
                                              forceWithoutIncrement( i ) - manualControlIncrements( i ) ), 1.0E-14 );
                BOOST_CHECK_SMALL( std::fabs( momentWithIncrement( i ) -
                                              momentWithoutIncrement( i ) - manualControlIncrements( i + 3 ) ), 1.0E-14 );
            }
        }
    }
}

//! Test use of control surface deflections in a full numerical propagation, with a dummy (e.g. non-physical) model
//! for aerodynamic and control surface guidance. Test case uses Apollo capsule entry and coefficients.
BOOST_AUTO_TEST_CASE( testControlSurfaceIncrementInterfaceInPropagation )
{
    using namespace tudat;
    using namespace ephemerides;
    using namespace interpolators;
    using namespace numerical_integrators;
    using namespace spice_interface;
    using namespace simulation_setup;
    using namespace basic_astrodynamics;
    using namespace orbital_element_conversions;
    using namespace propagators;
    using namespace aerodynamics;
    using namespace basic_mathematics;
    using namespace input_output;

    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );


    // Set simulation start epoch.
    const double simulationStartEpoch = 0.0;

    // Set simulation end epoch.
    const double simulationEndEpoch = 3300.0;

    // Set numerical integration fixed step size.
    const double fixedStepSize = 1.0;


    // Set initial Keplerian elements for vehicle.
    Vector6d apolloInitialStateInKeplerianElements;
    apolloInitialStateInKeplerianElements( semiMajorAxisIndex ) = spice_interface::getAverageRadius( "Earth" ) + 120.0E3;
    apolloInitialStateInKeplerianElements( eccentricityIndex ) = 0.005;
    apolloInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
    apolloInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
            = unit_conversions::convertDegreesToRadians( 235.7 );
    apolloInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
            = unit_conversions::convertDegreesToRadians( 23.4 );
    apolloInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );

    // Convert apollo state from Keplerian elements to Cartesian elements.
    const Vector6d apolloInitialState = convertKeplerianToCartesianElements(
                apolloInitialStateInKeplerianElements,
                getBodyGravitationalParameter( "Earth" ) );

    // Define simulation body settings.
    BodyListSettings bodySettings =
            getDefaultBodySettings( { "Earth", "Moon" }, simulationStartEpoch - 10.0 * fixedStepSize,
                                    simulationEndEpoch + 10.0 * fixedStepSize );
    bodySettings.at( "Earth" )->gravityFieldSettings =
            std::make_shared< simulation_setup::GravityFieldSettings >( central_spice );

    // Create Earth object
    simulation_setup::SystemOfBodies bodies = simulation_setup::createSystemOfBodies( bodySettings );

    // Create vehicle objects.
    bodies.createEmptyBody( "Apollo" );

    // Create vehicle aerodynamic coefficients
    bodies.at( "Apollo" )->setAerodynamicCoefficientInterface(
                unit_tests::getApolloCoefficientInterface( ) );
    std::shared_ptr< ControlSurfaceIncrementAerodynamicInterface > controlSurfaceInterface =
            std::make_shared< CustomControlSurfaceIncrementAerodynamicInterface >(
                &dummyControlIncrements,
                std::vector< AerodynamicCoefficientsIndependentVariables >{ angle_of_attack_dependent, control_surface_deflection_dependent } );
    std::map< std::string, std::shared_ptr< ControlSurfaceIncrementAerodynamicInterface >  > controlSurfaceList;
    controlSurfaceList[ "TestSurface" ] = controlSurfaceInterface;
    bodies.at( "Apollo" )->getAerodynamicCoefficientInterface( )->setControlSurfaceIncrements( controlSurfaceList );


    bodies.at( "Apollo" )->setConstantBodyMass( 5.0E3 );
    std::shared_ptr< system_models::VehicleSystems > apolloSystems = std::make_shared< system_models::VehicleSystems >( );
    bodies.at( "Apollo" )->setVehicleSystems( apolloSystems );


    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Define acceleration model settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfApollo;
    accelerationsOfApollo[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationsOfApollo[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( aerodynamic ) );
    accelerationsOfApollo[ "Moon" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationMap[ "Apollo" ] = accelerationsOfApollo;

    bodiesToPropagate.push_back( "Apollo" );
    centralBodies.push_back( "Earth" );

    // Set initial state
    Eigen::Vector6d systemInitialState = apolloInitialState;

    // Define list of dependent variables to save.
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
    dependentVariables.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    mach_number_dependent_variable, "Apollo" ) );
    dependentVariables.push_back(
                std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    "Apollo", reference_frames::angle_of_attack ) );
    dependentVariables.push_back(
                std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    "Apollo", reference_frames::angle_of_sideslip ) );
    dependentVariables.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    control_surface_deflection_dependent_variable, "Apollo", "TestSurface" ) );
    dependentVariables.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    aerodynamic_moment_coefficients_dependent_variable, "Apollo" ) );
    dependentVariables.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    aerodynamic_force_coefficients_dependent_variable, "Apollo" ) );


    // Create acceleration models
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodies, accelerationMap, bodiesToPropagate, centralBodies );

    // Set update function for body orientation and control surface deflections
    std::shared_ptr< DummyGuidanceSystem > dummyGuidanceSystem = std::make_shared< DummyGuidanceSystem >(
                std::bind( &system_models::VehicleSystems::setCurrentControlSurfaceDeflection, apolloSystems, std::placeholders::_1, std::placeholders::_2 ),
            bodies.at( "Apollo" )->getFlightConditions( )->getAerodynamicAngleCalculator( ) );

    // Create propagation and integrtion settings.
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState,
              std::make_shared< propagators::PropagationTimeTerminationSettings >( 1000.0 ), cowell,
              std::make_shared< DependentVariableSaveSettings >( dependentVariables ) );
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4, simulationStartEpoch, fixedStepSize );

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodies, integratorSettings, propagatorSettings, true, false, false );

    // Retrieve numerical solutions for state and dependent variables
    std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > numericalSolution =
            dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > dependentVariableSolution =
            dynamicsSimulator.getDependentVariableHistory( );


    // Declare test variables.
    double currentAngleOfAttack, currentSideslipAngle, currentMachNumber, currentSurfaceDeflection, currentTime;
    Eigen::Vector3d currentForceCoefficients, currentMomentCoefficients;
    Eigen::Vector3d expectedForceCoefficients, expectedMomentCoefficients;

    std::vector< double > currentAerodynamicsIndependentVariables;
    currentAerodynamicsIndependentVariables.resize( 3 );

    std::map< std::string, std::vector< double > > currentAerodynamicsControlIndependentVariables;
    currentAerodynamicsControlIndependentVariables[ "TestSurface" ].resize( 2 );

    // Iterate over saved variables and compare to expected values
    std::shared_ptr< AerodynamicCoefficientInterface > coefficientInterface =
            bodies.at( "Apollo" )->getAerodynamicCoefficientInterface( );
    for( std::map< double, Eigen::VectorXd >::iterator variableIterator = dependentVariableSolution.begin( );
         variableIterator != dependentVariableSolution.end( ); variableIterator++ )
    {
        // Retrieve dependent variables
        currentTime = variableIterator->first;
        currentMachNumber = variableIterator->second( 0 );
        currentAngleOfAttack = variableIterator->second( 1 );
        currentSideslipAngle = variableIterator->second( 2 );
        currentSurfaceDeflection = variableIterator->second( 3 );
        currentMomentCoefficients = variableIterator->second.segment( 4, 3 );
        currentForceCoefficients = variableIterator->second.segment( 7, 3  );

        // Test angles of attack and sideslip, and control surface deflection, against expectec values
        BOOST_CHECK_SMALL( std::fabs( currentAngleOfAttack - 0.3 * ( 1.0 - currentTime / 1000.0 ) ), 1.0E-14 );
        BOOST_CHECK_SMALL( std::fabs( currentSideslipAngle ), 1.0E-14 );
        BOOST_CHECK_SMALL( std::fabs( currentSurfaceDeflection - ( -0.02 + 0.04 * currentTime / 1000.0 ) ), 1.0E-14 );


        // Set current aerodynamic coefficient independent variables and retrieve coefficients.c.
        currentAerodynamicsIndependentVariables[ 0 ] = currentMachNumber;
        currentAerodynamicsIndependentVariables[ 1 ] = currentAngleOfAttack;
        currentAerodynamicsIndependentVariables[ 2 ] = currentSideslipAngle;

        currentAerodynamicsControlIndependentVariables[ "TestSurface" ][ 0 ] = currentAngleOfAttack;
        currentAerodynamicsControlIndependentVariables[ "TestSurface" ][ 1 ] = currentSurfaceDeflection;

        coefficientInterface->updateFullCurrentCoefficients(
                    currentAerodynamicsIndependentVariables, currentAerodynamicsControlIndependentVariables );

        expectedForceCoefficients = coefficientInterface->getCurrentForceCoefficients( );
        expectedMomentCoefficients = coefficientInterface->getCurrentMomentCoefficients( );

        // Test expected against actual aerodynamic coefficients.
        for( unsigned int i = 0; i < 3; i++ )
        {
            BOOST_CHECK_SMALL( std::fabs( expectedForceCoefficients( i ) - currentForceCoefficients( i ) ), 1.0E-14 );
            BOOST_CHECK_SMALL( std::fabs( expectedMomentCoefficients( i ) - currentMomentCoefficients( i ) ), 1.0E-14 );
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )


} // namespace unit_tests

} // namespace tudat

