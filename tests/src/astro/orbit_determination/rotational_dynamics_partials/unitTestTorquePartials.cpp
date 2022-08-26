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
#include <string>
#include "tudat/basics/testMacros.h"
#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/astro/basic_astro/unitConversions.h"

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>
#include <boost/lambda/lambda.hpp>

#include "tudat/astro/basic_astro/sphericalStateConversions.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/astro/ephemerides/keplerEphemeris.h"
#include "tudat/astro/relativity/metric.h"
#include "tudat/astro/orbit_determination/acceleration_partials/numericalAccelerationPartial.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/gravitationalParameter.h"
#include "tudat/simulation/estimation_setup/createTorquePartials.h"
#include "tudat/simulation/environment_setup/createBodies.h"
#include "tudat/simulation/propagation_setup/createTorqueModel.h"
#include "tudat/simulation/estimation_setup/createEstimatableParameters.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"

namespace tudat
{

namespace unit_tests
{

using namespace tudat::relativity;
using namespace tudat::gravitation;
using namespace tudat::aerodynamics;
using namespace tudat::ephemerides;
using namespace tudat::simulation_setup;
using namespace tudat::orbital_element_conversions;
using namespace tudat::unit_conversions;
using namespace tudat::orbit_determination;
using namespace tudat::acceleration_partials;
using namespace tudat::spice_interface;
using namespace tudat::orbit_determination;
using namespace tudat::estimatable_parameters;
using namespace tudat::electromagnetism;
using namespace tudat::basic_astrodynamics;

BOOST_AUTO_TEST_SUITE( test_torque_partials )

//! Test if partial derivatives of degree 2 torque are correctly implemented
BOOST_AUTO_TEST_CASE( testSecondDegreeGravitationalTorquePartials )
{
    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Create bodies
    SystemOfBodies bodies = SystemOfBodies( "Mars", "ECLIPJ2000" );
    bodies.createEmptyBody( "Mars", false );
    bodies.at( "Mars" )->setEphemeris( std::make_shared< ephemerides::ConstantEphemeris >(
                                         [ = ]( ){ return Eigen::Vector6d::Zero( ); } ) );
    bodies.at( "Mars" )->setGravityFieldModel(
                std::make_shared< gravitation::GravityFieldModel >(
                    spice_interface::getBodyGravitationalParameter( "Mars" ) ) );
    bodies.createEmptyBody( "Phobos" );
    std::shared_ptr< Body > phobos = bodies.at( "Phobos" );
    std::shared_ptr< Body > mars = bodies.at( "Mars" );

    // Set inertia tensor
    Eigen::Matrix3d phobosInertiaTensor = Eigen::Matrix3d::Zero( );
    phobosInertiaTensor( 0, 0 ) = 0.3615;
    phobosInertiaTensor( 1, 1 ) = 0.4265;
    phobosInertiaTensor( 2, 2 ) = 0.5024;

    phobosInertiaTensor *= ( 11.27E3 * 11.27E3 * 1.0659E16 );
    bodies.at( "Phobos" )->setBodyInertiaTensor(
                phobosInertiaTensor, ( 0.3615 + 0.4265 + 0.5024 ) / 3.0 );

    // Create gravity field
    double phobosGravitationalParameter = 1.0659E16 * physical_constants::GRAVITATIONAL_CONSTANT;
    double phobosReferenceRadius = 11.27E3;

    Eigen::MatrixXd phobosCosineGravityFieldCoefficients = Eigen::MatrixXd::Zero( 6, 6 ),
            phobosSineGravityFieldCoefficients = Eigen::MatrixXd::Zero( 6, 6 );
    double phobosScaledMeanMomentOfInertia;
    gravitation::getDegreeTwoSphericalHarmonicCoefficients(
                phobosInertiaTensor, phobosGravitationalParameter, phobosReferenceRadius, true,
                phobosCosineGravityFieldCoefficients, phobosSineGravityFieldCoefficients, phobosScaledMeanMomentOfInertia );
    bodies.at( "Phobos" )->setGravityFieldModel(
                std::make_shared< gravitation::SphericalHarmonicsGravityField >(
                    phobosGravitationalParameter, phobosReferenceRadius, phobosCosineGravityFieldCoefficients,
                    phobosSineGravityFieldCoefficients, "Phobos_Fixed",
                    std::bind( &Body::setBodyInertiaTensorFromGravityFieldAndExistingMeanMoment, bodies.at( "Phobos" ), true ) ) );

    Eigen::Quaterniond noRotationQuaternion = Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) );
    Eigen::Matrix< double, 7, 1 > unitRotationState = Eigen::Matrix< double, 7, 1 >::Zero( );
    unitRotationState( 0 ) = noRotationQuaternion.w( );
    unitRotationState( 1 ) = noRotationQuaternion.x( );
    unitRotationState( 2 ) = noRotationQuaternion.y( );
    unitRotationState( 3 ) = noRotationQuaternion.z( );

    // Set rotation and ephemeris models
    std::map< double, Eigen::Matrix< double, 7, 1 > > dummyRotationMap;
    dummyRotationMap[ -1.0E100 ] = unitRotationState;
    dummyRotationMap[ 1.0E100 ] = unitRotationState;
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< double, 7, 1 > > > dummyInterpolator =
            std::make_shared< interpolators::LinearInterpolator< double, Eigen::Matrix< double, 7, 1 > > >( dummyRotationMap );
    bodies.at( "Phobos" )->setRotationalEphemeris( std::make_shared< TabulatedRotationalEphemeris< double, double > >(
                                                     dummyInterpolator, "ECLIPJ2000", "Phobos_Fixed" ) );


    Eigen::Vector6d phobosKeplerElements = Eigen::Vector6d::Zero( );
    double phobosSemiMajorAxis = 9376.0E3;
    phobosKeplerElements( 0 ) = phobosSemiMajorAxis;
    bodies.at( "Phobos" )->setEphemeris( std::make_shared< ephemerides::KeplerEphemeris >(
                                           phobosKeplerElements, 0.0, spice_interface::getBodyGravitationalParameter( "Mars" ),
                                           "Mars", "ECLIPJ2000" ) );


    // Update Phobos and Mars to current state
    double testTime = 1000.0;
    phobos->setStateFromEphemeris( testTime );
    mars->setStateFromEphemeris( testTime );

    // Set new Phobos rotational state
    Eigen::Vector7d phobosRotationalState = Eigen::Vector7d::Zero( );
    phobosRotationalState.segment( 0, 4 ) = tudat::linear_algebra::convertQuaternionToVectorFormat(
                Eigen::Quaterniond( Eigen::AngleAxisd( 0.4343, Eigen::Vector3d::UnitZ( ) ) *
                                    Eigen::AngleAxisd( 2.4354, Eigen::Vector3d::UnitX( ) ) *
                                    Eigen::AngleAxisd( 1.2434, Eigen::Vector3d::UnitY( ) ) ) );
    phobos->setCurrentRotationalStateToLocalFrame( phobosRotationalState );

    // Create torque due to mars on phobos.
    std::shared_ptr< SecondDegreeGravitationalTorqueModel > gravitationalTorque =
            createSecondDegreeGravitationalTorqueModel( bodies.at( "Phobos" ), bodies.at( "Mars" ), "Phobos", "Mars" );

    // Create torque partial.
    std::shared_ptr< TorquePartial > torquePartial =
            createAnalyticalTorquePartial( gravitationalTorque, std::make_pair( "Phobos", phobos ),
                                           std::make_pair( "Mars", mars ) );

    // Create parameter objects
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Mars", gravitational_parameter) );
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Phobos", mean_moment_of_inertia ) );
    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                  1, 0, 3, 3, "Phobos", spherical_harmonics_cosine_coefficient_block ) );
    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                  1, 1, 3, 3, "Phobos", spherical_harmonics_sine_coefficient_block ) );
    std::shared_ptr< EstimatableParameterSet< double > > parameterSet =
            createParametersToEstimate( parameterNames, bodies );

    std::shared_ptr< EstimatableParameter< double > > marsGravitationalParameterParameter =
            parameterSet->getEstimatedDoubleParameters( ).at( 0 );
    std::shared_ptr< EstimatableParameter< double > > phobosMeanMomentOfInertia =
            parameterSet->getEstimatedDoubleParameters( ).at( 1 );
    std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > phobosCosineCoefficientsParameter =
            parameterSet->getEstimatedVectorParameters( ).at( 0 );
    std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > phobosSineCoefficientsParameter =
            parameterSet->getEstimatedVectorParameters( ).at( 1 );

    // Calculate analytical partials.
    torquePartial->update( testTime );

    Eigen::MatrixXd partialWrtPhobosOrientation = Eigen::MatrixXd::Zero( 3, 4 );
    torquePartial->wrtOrientationOfAcceleratedBody( partialWrtPhobosOrientation.block( 0, 0, 3, 4 ) );

    Eigen::MatrixXd partialWrtPhobosRotationalVelocity = Eigen::Matrix3d::Zero( );
    torquePartial->wrtRotationalVelocityOfAcceleratedBody( partialWrtPhobosRotationalVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );

    Eigen::MatrixXd partialWrtMarsOrientation = Eigen::MatrixXd::Zero( 3, 4  );
    torquePartial->wrtOrientationOfAcceleratingBody( partialWrtMarsOrientation.block( 0, 0, 3, 4 ) );

    Eigen::MatrixXd partialWrtMarsRotationalVelocity = Eigen::Matrix3d::Zero( );
    torquePartial->wrtRotationalVelocityOfAcceleratingBody( partialWrtMarsRotationalVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );

    Eigen::MatrixXd partialWrtPhobosState = Eigen::MatrixXd::Zero( 3, 6 );
    torquePartial->wrtNonRotationalStateOfAdditionalBody(
                partialWrtPhobosState.block( 0, 0, 3, 6 ), std::make_pair( "Phobos", "" ), propagators::translational_state );
    Eigen::MatrixXd partialWrtMarsState = Eigen::MatrixXd::Zero( 3, 6 );
    torquePartial->wrtNonRotationalStateOfAdditionalBody(
                partialWrtMarsState.block( 0, 0, 3, 6 ), std::make_pair( "Mars", "" ), propagators::translational_state );

    Eigen::Vector3d partialWrtMarsGravitationalParameter = torquePartial->wrtParameter(
                marsGravitationalParameterParameter );
    Eigen::Vector3d partialWrtMeanMomentOfInertia = torquePartial->wrtParameter(
                phobosMeanMomentOfInertia );

    Eigen::MatrixXd partialWrtPhobosCosineCoefficients = torquePartial->wrtParameter(
                phobosCosineCoefficientsParameter );
    Eigen::MatrixXd partialWrtPhobosSineCoefficients = torquePartial->wrtParameter(
                phobosSineCoefficientsParameter );


    // Declare numerical partials.
    Eigen::Matrix3d testPartialWrtPhobosRotationalVelocity = Eigen::Matrix3d::Zero( );
    Eigen::Matrix< double, 3, 4 > testPartialWrtMarsOrientation = Eigen::Matrix< double, 3, 4 >::Zero( );
    Eigen::Matrix3d testPartialWrtMarsRotationalVelocity = Eigen::Matrix3d::Zero( );

    Eigen::Matrix3d testPartialWrtMarsPosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtMarsVelocity = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtPhobosPosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtPhobosVelocity = Eigen::Matrix3d::Zero( );

    // Define perturbations in orientation for numerical partial/
    Eigen::Vector4d orientationPerturbation;
    orientationPerturbation << 1.0E-9, 1.0E-9, 1.0E-9, 1.0E-9;
    Eigen::Vector3d rotationalVelocityPerturbation;
    rotationalVelocityPerturbation <<  1.0E-6, 1.0E-6, 1.0E-6;

    // Create state access/modification functions for bodies.
    std::function< void( Eigen::Vector7d ) > phobosRotationalStateSetFunction =
            std::bind( &Body::setCurrentRotationalStateToLocalFrame, phobos, std::placeholders::_1 );
    std::function< void( Eigen::Vector7d ) > marsRotationalStateSetFunction =
            std::bind( &Body::setCurrentRotationalStateToLocalFrame, mars, std::placeholders::_1 );

    // Calculate numerical partials wrt rotational state.
    std::vector< Eigen::Vector4d > appliedQuaternionPerturbation;
    Eigen::MatrixXd torqueDeviations = calculateTorqueDeviationDueToOrientationChange(
                phobosRotationalStateSetFunction, gravitationalTorque, phobos->getRotationalStateVector( ), orientationPerturbation,
                appliedQuaternionPerturbation );

    testPartialWrtPhobosRotationalVelocity =  calculateTorqueWrtRotationalStatePartials(
                phobosRotationalStateSetFunction, gravitationalTorque, phobos->getRotationalStateVector( ), rotationalVelocityPerturbation, 4, 3 );
    testPartialWrtMarsOrientation =  calculateTorqueWrtRotationalStatePartials(
                marsRotationalStateSetFunction, gravitationalTorque, mars->getRotationalStateVector( ), orientationPerturbation, 0, 4 );
    testPartialWrtMarsRotationalVelocity =  calculateTorqueWrtRotationalStatePartials(
                marsRotationalStateSetFunction, gravitationalTorque, mars->getRotationalStateVector( ), rotationalVelocityPerturbation, 4, 3 );

    std::function< void( Eigen::Vector6d ) > phobosStateSetFunction =
            std::bind( &Body::setState,  phobos, std::placeholders::_1 );
    std::function< void( Eigen::Vector6d ) > marsStateSetFunction =
            std::bind( &Body::setState, mars, std::placeholders::_1 );

    Eigen::Vector3d positionPerturbation;
    positionPerturbation << 1.0, 1.0, 100.0;
    Eigen::Vector3d velocityPerturbation;
    velocityPerturbation << 1.0E-3, 1.0E-3, 1.0E-3;

    // Calculate numerical partials wrt translational state.
    testPartialWrtMarsPosition = calculateTorqueWrtTranslationalStatePartials(
                marsStateSetFunction, gravitationalTorque, mars->getState( ), positionPerturbation, 0 );
    testPartialWrtMarsVelocity = calculateTorqueWrtTranslationalStatePartials(
                marsStateSetFunction, gravitationalTorque, mars->getState( ), velocityPerturbation, 3 );
    testPartialWrtPhobosPosition = calculateTorqueWrtTranslationalStatePartials(
                phobosStateSetFunction, gravitationalTorque, phobos->getState( ), positionPerturbation, 0 );
    testPartialWrtPhobosVelocity = calculateTorqueWrtTranslationalStatePartials(
                phobosStateSetFunction, gravitationalTorque, phobos->getState( ), velocityPerturbation, 3 );

    // Calculate numerical partials wrt parameters
    Eigen::Vector3d testPartialWrtMarsGravitationalParameter = calculateTorqueWrtParameterPartials(
                marsGravitationalParameterParameter, gravitationalTorque, 1.0E12 );
    Eigen::Vector3d testPartialWrtMeanMomentOfInertia = calculateTorqueWrtParameterPartials(
                phobosMeanMomentOfInertia, gravitationalTorque, 1.0E-1 );
    std::function< void( ) > updateFunction = &emptyFunction;
            //boost::bind( &Body::setBodyInertiaTensorFromGravityFieldAndExistingMeanMoment, bodies.at( "Phobos" ), true );
    Eigen::MatrixXd testPartialWrtPhobosCosineCoefficients = calculateTorqueWrtParameterPartials(
                phobosCosineCoefficientsParameter, gravitationalTorque,
                Eigen::VectorXd::Constant( phobosCosineCoefficientsParameter->getParameterSize( ), 1.0E-6 ), updateFunction );
    Eigen::MatrixXd testPartialWrtPhobosSineCoefficients = calculateTorqueWrtParameterPartials(
                phobosSineCoefficientsParameter, gravitationalTorque,
                Eigen::VectorXd::Constant( phobosSineCoefficientsParameter->getParameterSize( ), 1.0E-6 ), updateFunction );

    // Compare numerical and analytical results.
    for( int index = 1; index < 4; index++ )
    {
        Eigen::Vector3d numericalChangeInTorque = torqueDeviations.block( 0, index - 1, 3, 1 );
        Eigen::Vector3d analyticalChangeInTorque =
                partialWrtPhobosOrientation.block( 0, 0, 3, 1 ) * appliedQuaternionPerturbation[ index ]( 0 ) +
                partialWrtPhobosOrientation.block( 0, index, 3, 1 ) * appliedQuaternionPerturbation[ index ]( index );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( numericalChangeInTorque,
                                           analyticalChangeInTorque, 1.0E-6 );
    }

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPhobosRotationalVelocity,
                                       partialWrtPhobosRotationalVelocity, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtMarsOrientation,
                                       partialWrtMarsOrientation, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtMarsRotationalVelocity,
                                       partialWrtMarsRotationalVelocity, std::numeric_limits< double >::epsilon( ) );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(  partialWrtPhobosState.block( 0, 0, 3, 3 ),
                                        testPartialWrtPhobosPosition, 1.0E-8 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(  partialWrtMarsState.block( 0, 0, 3, 3 ),
                                        testPartialWrtMarsPosition, 1.0E-8 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(  partialWrtPhobosState.block( 0, 3, 3, 3 ),
                                        testPartialWrtPhobosVelocity, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(  partialWrtMarsState.block( 0, 3, 3, 3 ),
                                        testPartialWrtMarsVelocity, std::numeric_limits< double >::epsilon( ) );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtMarsGravitationalParameter,
                                       partialWrtMarsGravitationalParameter, 1.0E-6 );

    // Check derivative of z-component w.r.t. C20 separately: value is slightly non-zero due to rounding error
    BOOST_CHECK_SMALL( std::fabs( testPartialWrtPhobosCosineCoefficients( 2, 2 ) ), 1.0E6 );
    testPartialWrtPhobosCosineCoefficients( 2, 2 ) = 0.0;

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPhobosCosineCoefficients,
                                       partialWrtPhobosCosineCoefficients, 1.0E-9 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPhobosSineCoefficients,
                                       partialWrtPhobosSineCoefficients, 1.0E-9 );

    for( int i = 0; i < 3; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( partialWrtMeanMomentOfInertia( i, 0 ) - testPartialWrtMeanMomentOfInertia( i, 0 ) ), 1.0E2 );
    }
}

BOOST_AUTO_TEST_CASE( testInertialTorquePartials )
{
    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Create bodies
    SystemOfBodies bodies = SystemOfBodies( "Mars", "ECLIPJ2000" );
    bodies.createEmptyBody( "Mars", false );
    bodies.at( "Mars" )->setEphemeris( std::make_shared< ephemerides::ConstantEphemeris >(
                                         [ = ]( ){ return Eigen::Vector6d::Zero( ); } ) );
    bodies.at( "Mars" )->setGravityFieldModel(
                std::make_shared< gravitation::GravityFieldModel >(
                    spice_interface::getBodyGravitationalParameter( "Mars" ) ) );
    bodies.createEmptyBody( "Phobos" );
    std::shared_ptr< Body > phobos = bodies.at( "Phobos" );
    std::shared_ptr< Body > mars = bodies.at( "Mars" );

    // Set inertia tensor
    Eigen::Matrix3d phobosInertiaTensor = Eigen::Matrix3d::Zero( );
    phobosInertiaTensor( 0, 0 ) = 0.3615;
    phobosInertiaTensor( 1, 1 ) = 0.4265;
    phobosInertiaTensor( 2, 2 ) = 0.5024;

    phobosInertiaTensor *= ( 11.27E3 * 11.27E3 * 1.0659E16 );
    bodies.at( "Phobos" )->setBodyInertiaTensor(
                phobosInertiaTensor, ( 0.3615 + 0.4265 + 0.5024 ) / 3.0 );

    // Create gravity field
    double phobosGravitationalParameter = 1.0659E16 * physical_constants::GRAVITATIONAL_CONSTANT;
    double phobosReferenceRadius = 11.27E3;

    Eigen::MatrixXd phobosCosineGravityFieldCoefficients = Eigen::MatrixXd::Zero( 6, 6 ),
            phobosSineGravityFieldCoefficients = Eigen::MatrixXd::Zero( 6, 6 );
    double phobosScaledMeanMomentOfInertia;
    gravitation::getDegreeTwoSphericalHarmonicCoefficients(
                phobosInertiaTensor, phobosGravitationalParameter, phobosReferenceRadius, true,
                phobosCosineGravityFieldCoefficients, phobosSineGravityFieldCoefficients, phobosScaledMeanMomentOfInertia );

    bodies.at( "Phobos" )->setGravityFieldModel(
                std::make_shared< gravitation::SphericalHarmonicsGravityField >(
                    phobosGravitationalParameter, phobosReferenceRadius, phobosCosineGravityFieldCoefficients,
                    phobosSineGravityFieldCoefficients, "Phobos_Fixed" ) );

    Eigen::Quaterniond noRotationQuaternion = Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) );
    Eigen::Matrix< double, 7, 1 > unitRotationState = Eigen::Matrix< double, 7, 1 >::Zero( );
    unitRotationState( 0 ) = noRotationQuaternion.w( );
    unitRotationState( 1 ) = noRotationQuaternion.x( );
    unitRotationState( 2 ) = noRotationQuaternion.y( );
    unitRotationState( 3 ) = noRotationQuaternion.z( );
    unitRotationState( 4 ) = 1.0E-5;
    unitRotationState( 5 ) = 2.0E-5;
    unitRotationState( 6 ) = -3.2E-5;

    // Set rotation and ephemeris models
    std::map< double, Eigen::Matrix< double, 7, 1 > > dummyRotationMap;
    dummyRotationMap[ -1.0E100 ] = unitRotationState;
    dummyRotationMap[ 1.0E100 ] = unitRotationState;
    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< double, 7, 1 > > > dummyInterpolator =
            std::make_shared< interpolators::LinearInterpolator< double, Eigen::Matrix< double, 7, 1 > > >( dummyRotationMap );
    bodies.at( "Phobos" )->setRotationalEphemeris( std::make_shared< TabulatedRotationalEphemeris< double, double > >(
                                                     dummyInterpolator, "ECLIPJ2000", "Phobos_Fixed" ) );

    Eigen::Vector6d phobosKeplerElements = Eigen::Vector6d::Zero( );
    double phobosSemiMajorAxis = 9376.0E3;
    phobosKeplerElements( 0 ) = phobosSemiMajorAxis;
    bodies.at( "Phobos" )->setEphemeris( std::make_shared< ephemerides::KeplerEphemeris >(
                                           phobosKeplerElements, 0.0, spice_interface::getBodyGravitationalParameter( "Mars" ),
                                           "Mars", "ECLIPJ2000" ) );

    // Update Phobos and Mars to current state
    double testTime = 1000.0;
    phobos->setStateFromEphemeris( testTime );
    mars->setStateFromEphemeris( testTime );

    // Set new Phobos rotational state
    Eigen::Vector7d phobosRotationalState = Eigen::Vector7d::Zero( );
    phobosRotationalState.segment( 0, 4 ) = tudat::linear_algebra::convertQuaternionToVectorFormat(
                Eigen::Quaterniond( Eigen::AngleAxisd( 0.4343, Eigen::Vector3d::UnitZ( ) ) *
                                    Eigen::AngleAxisd( 2.4354, Eigen::Vector3d::UnitX( ) ) *
                                    Eigen::AngleAxisd( 1.2434, Eigen::Vector3d::UnitY( ) ) ) );
    phobosRotationalState( 4 ) = 1.0E-5;
    phobosRotationalState( 5 ) = 2.0E-5;
    phobosRotationalState( 6 ) = -3.2E-5;
    phobos->setCurrentRotationalStateToLocalFrame( phobosRotationalState );

    // Create torque due to mars on phobos.
    std::shared_ptr< InertialTorqueModel > inertialTorqueModel =
            createInertialTorqueModel( bodies.at( "Phobos" ), "Phobos" );
    inertialTorqueModel->updateMembers( 0.0 );

    SingleBodyTorqueModelMap torqueList;
    torqueList[ "Phobos" ].push_back( inertialTorqueModel );

    // Create shared_ptr partial.
    std::shared_ptr< TorquePartial > torquePartial =
            createAnalyticalTorquePartial( inertialTorqueModel, std::make_pair( "Phobos", phobos ),
                                           std::make_pair( "Phobos", phobos ) );


    // Create parameter objects.
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Phobos", gravitational_parameter) );
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Phobos", mean_moment_of_inertia ) );

    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                  1, 0, 3, 3, "Phobos", spherical_harmonics_cosine_coefficient_block ) );
    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                  1, 1, 3, 3, "Phobos", spherical_harmonics_sine_coefficient_block ) );

    std::shared_ptr< EstimatableParameterSet< double > > parameterSet =
            createParametersToEstimate( parameterNames, bodies );

    std::shared_ptr< EstimatableParameter< double > > phobosGravitationalParameterParameter =
            parameterSet->getEstimatedDoubleParameters( ).at( 0 );
    std::shared_ptr< EstimatableParameter< double > > meanMomentOfInertiaParameter =
            parameterSet->getEstimatedDoubleParameters( ).at( 1 );
    std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > phobosCosineCoefficientsParameter =
            parameterSet->getEstimatedVectorParameters( ).at( 0 );
    std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > phobosSineCoefficientsParameter =
            parameterSet->getEstimatedVectorParameters( ).at( 1 );

    // Calculate analytical partials.
    torquePartial->update( testTime );

    Eigen::MatrixXd partialWrtPhobosOrientation = Eigen::MatrixXd::Zero( 3, 4 );
    torquePartial->wrtOrientationOfAcceleratedBody( partialWrtPhobosOrientation.block( 0, 0, 3, 4 ) );
    Eigen::MatrixXd partialWrtPhobosRotationalVelocity = Eigen::Matrix3d::Zero( );
    torquePartial->wrtRotationalVelocityOfAcceleratedBody( partialWrtPhobosRotationalVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );
    Eigen::MatrixXd partialWrtMarsOrientation = Eigen::MatrixXd::Zero( 3, 4  );
    torquePartial->wrtOrientationOfAcceleratingBody( partialWrtMarsOrientation.block( 0, 0, 3, 4 ) );
    Eigen::MatrixXd partialWrtMarsRotationalVelocity = Eigen::Matrix3d::Zero( );
    torquePartial->wrtRotationalVelocityOfAcceleratingBody( partialWrtMarsRotationalVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );

    Eigen::MatrixXd partialWrtPhobosState = Eigen::MatrixXd::Zero( 3, 6 );
    torquePartial->wrtNonRotationalStateOfAdditionalBody(
                partialWrtPhobosState.block( 0, 0, 3, 6 ), std::make_pair( "Phobos", "" ), propagators::translational_state );
    Eigen::MatrixXd partialWrtMarsState = Eigen::MatrixXd::Zero( 3, 6 );
    torquePartial->wrtNonRotationalStateOfAdditionalBody(
                partialWrtMarsState.block( 0, 0, 3, 6 ), std::make_pair( "Mars", "" ), propagators::translational_state );

    Eigen::Vector3d partialWrtPhobosGravitationalParameter = torquePartial->wrtParameter(
                phobosGravitationalParameterParameter );
    Eigen::Vector3d partialWrtMeanMomentOfInertia = torquePartial->wrtParameter(
                meanMomentOfInertiaParameter );

    Eigen::MatrixXd partialWrtPhobosCosineCoefficients = torquePartial->wrtParameter(
                phobosCosineCoefficientsParameter );
    Eigen::MatrixXd partialWrtPhobosSineCoefficients = torquePartial->wrtParameter(
                phobosSineCoefficientsParameter );


    // Declare numerical partials.
    Eigen::Matrix3d testPartialWrtPhobosRotationalVelocity = Eigen::Matrix3d::Zero( );
    Eigen::Matrix< double, 3, 4 > testPartialWrtMarsOrientation = Eigen::Matrix< double, 3, 4 >::Zero( );
    Eigen::Matrix3d testPartialWrtMarsRotationalVelocity = Eigen::Matrix3d::Zero( );

    Eigen::Matrix3d testPartialWrtMarsPosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtMarsVelocity = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtPhobosPosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtPhobosVelocity = Eigen::Matrix3d::Zero( );

    // Declare perturbations in orientation for numerical partial
    Eigen::Vector4d orientationPerturbation;
    orientationPerturbation << 1.0E-9, 1.0E-9, 1.0E-9, 1.0E-9;
    Eigen::Vector3d rotationalVelocityPerturbation;
    rotationalVelocityPerturbation <<  1.0E-6, 1.0E-6, 1.0E-6;

    // Create state access/modification functions for bodies.
    std::function< void( Eigen::Vector7d ) > phobosRotationalStateSetFunction =
            std::bind( &Body::setCurrentRotationalStateToLocalFrame, phobos, std::placeholders::_1 );
    std::function< void( Eigen::Vector7d ) > marsRotationalStateSetFunction =
            std::bind( &Body::setCurrentRotationalStateToLocalFrame, mars, std::placeholders::_1 );

    // Calculate numerical partials.
    std::vector< Eigen::Vector4d > appliedQuaternionPerturbation;
    Eigen::MatrixXd torqueDeviations = calculateTorqueDeviationDueToOrientationChange(
                phobosRotationalStateSetFunction, inertialTorqueModel, phobos->getRotationalStateVector( ), orientationPerturbation,
                appliedQuaternionPerturbation );

    testPartialWrtPhobosRotationalVelocity =  calculateTorqueWrtRotationalStatePartials(
                phobosRotationalStateSetFunction, inertialTorqueModel, phobos->getRotationalStateVector( ), rotationalVelocityPerturbation, 4, 3 );
    testPartialWrtMarsOrientation =  calculateTorqueWrtRotationalStatePartials(
                marsRotationalStateSetFunction, inertialTorqueModel, mars->getRotationalStateVector( ), orientationPerturbation, 0, 4 );
    testPartialWrtMarsRotationalVelocity =  calculateTorqueWrtRotationalStatePartials(
                marsRotationalStateSetFunction, inertialTorqueModel, mars->getRotationalStateVector( ), rotationalVelocityPerturbation, 4, 3 );

    std::function< void( Eigen::Vector6d ) > phobosStateSetFunction =
            std::bind( &Body::setState,  phobos, std::placeholders::_1 );
    std::function< void( Eigen::Vector6d ) > marsStateSetFunction =
            std::bind( &Body::setState, mars, std::placeholders::_1 );

    Eigen::Vector3d positionPerturbation;
    positionPerturbation << 1.0, 1.0, 100.0;
    Eigen::Vector3d velocityPerturbation;
    velocityPerturbation << 1.0E-3, 1.0E-3, 1.0E-3;

    // Calculate numerical partials.
    testPartialWrtMarsPosition = calculateTorqueWrtTranslationalStatePartials(
                marsStateSetFunction, inertialTorqueModel, mars->getState( ), positionPerturbation, 0 );
    testPartialWrtMarsVelocity = calculateTorqueWrtTranslationalStatePartials(
                marsStateSetFunction, inertialTorqueModel, mars->getState( ), velocityPerturbation, 3 );
    testPartialWrtPhobosPosition = calculateTorqueWrtTranslationalStatePartials(
                phobosStateSetFunction, inertialTorqueModel, phobos->getState( ), positionPerturbation, 0 );
    testPartialWrtPhobosVelocity = calculateTorqueWrtTranslationalStatePartials(
                phobosStateSetFunction, inertialTorqueModel, phobos->getState( ), velocityPerturbation, 3 );


    std::function< void( ) > updateFunction = //&emptyFunction;
            std::bind( &Body::setBodyInertiaTensorFromGravityFieldAndExistingMeanMoment, bodies.at( "Phobos" ), true );
    Eigen::Vector3d testPartialWrtPhobosGravitationalParameter = calculateTorqueWrtParameterPartials(
                phobosGravitationalParameterParameter, inertialTorqueModel, 1.0E8, updateFunction );
    Eigen::MatrixXd testPartialWrtMeanMomentOfInertia = calculateTorqueWrtParameterPartials(
                meanMomentOfInertiaParameter, inertialTorqueModel, 1.0E-1, updateFunction );
    Eigen::MatrixXd testPartialWrtPhobosCosineCoefficients = calculateTorqueWrtParameterPartials(
                phobosCosineCoefficientsParameter, inertialTorqueModel,
                Eigen::VectorXd::Constant( phobosCosineCoefficientsParameter->getParameterSize( ), 1.0E-6 ), updateFunction );
    Eigen::MatrixXd testPartialWrtPhobosSineCoefficients = calculateTorqueWrtParameterPartials(
                phobosSineCoefficientsParameter, inertialTorqueModel,
                Eigen::VectorXd::Constant( phobosSineCoefficientsParameter->getParameterSize( ), 1.0E-6 ), updateFunction );

    //    // Compare numerical and analytical results.
    for( int index = 1; index < 4; index++ )
    {
        Eigen::Vector3d numericalChangeInTorque = torqueDeviations.block( 0, index - 1, 3, 1 );
        Eigen::Vector3d analyticalChangeInTorque =
                partialWrtPhobosOrientation.block( 0, 0, 3, 1 ) * appliedQuaternionPerturbation[ index ]( 0 ) +
                partialWrtPhobosOrientation.block( 0, index, 3, 1 ) * appliedQuaternionPerturbation[ index ]( index );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( numericalChangeInTorque,
                                           analyticalChangeInTorque, std::numeric_limits< double >::epsilon( ) );
    }

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPhobosRotationalVelocity,
                                       partialWrtPhobosRotationalVelocity, 1.0E-10 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtMarsOrientation,
                                       partialWrtMarsOrientation, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtMarsRotationalVelocity,
                                       partialWrtMarsRotationalVelocity, std::numeric_limits< double >::epsilon( ) );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(  partialWrtPhobosState.block( 0, 0, 3, 3 ),
                                        testPartialWrtPhobosPosition, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(  partialWrtMarsState.block( 0, 0, 3, 3 ),
                                        testPartialWrtMarsPosition, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(  partialWrtPhobosState.block( 0, 3, 3, 3 ),
                                        testPartialWrtPhobosVelocity, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(  partialWrtMarsState.block( 0, 3, 3, 3 ),
                                        testPartialWrtMarsVelocity, std::numeric_limits< double >::epsilon( ) );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPhobosGravitationalParameter,
                                       partialWrtPhobosGravitationalParameter, 1.0E-6 );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtPhobosCosineCoefficients.block( 0, 0, 3, 2 ),
                                       Eigen::MatrixXd::Zero( 3, 2 ), std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtPhobosCosineCoefficients.block( 0, 5, 3, 4 ),
                                       Eigen::MatrixXd::Zero( 3, 4 ), std::numeric_limits< double >::epsilon( ) );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtPhobosSineCoefficients.block( 0, 0, 3, 1 ),
                                       Eigen::MatrixXd::Zero( 3, 1 ), std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtPhobosSineCoefficients.block( 0, 3, 3, 3 ),
                                       Eigen::MatrixXd::Zero( 3, 3 ), std::numeric_limits< double >::epsilon( ) );

    // Check derivative of z-component w.r.t. C20 separately: value is slightly non-zero due to rounding error
    BOOST_CHECK_SMALL( std::fabs( testPartialWrtPhobosCosineCoefficients( 2, 2 ) ), 1.0E5 );
    testPartialWrtPhobosCosineCoefficients( 2, 2 ) = 0.0;

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPhobosCosineCoefficients,
                                       partialWrtPhobosCosineCoefficients, 1.0E-9 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPhobosSineCoefficients,
                                       partialWrtPhobosSineCoefficients, 1.0E-9 );

    for( int i = 0; i < 3; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( partialWrtMeanMomentOfInertia( i, 0 ) - testPartialWrtMeanMomentOfInertia( i, 0 ) ), 1.0E1 );
    }
}

class EffectiveTorqueModel: public basic_astrodynamics::TorqueModel
{
public:
    EffectiveTorqueModel(
            const std::function< Eigen::Matrix3d( ) > inertiaTensorFunction,
            const SingleBodyTorqueModelMap& torqueList ):
        inertiaTensorFunction_( inertiaTensorFunction ), torqueList_( torqueList ){ }


    Eigen::Vector3d getTorque( )
    {
        return effectiveTorque_;
    }

    void updateMembers( const double currentTime )
    {
        Eigen::Vector3d totalTorque = Eigen::Vector3d::Zero( );
        for( auto it = torqueList_.begin( ); it != torqueList_.end( ); it++ )
        {
            for( unsigned int i = 0; i < it->second.size( ); i++ )
            {
                totalTorque += it->second.at( i )->getTorque( );
            }
        }

        effectiveTorque_ = inertiaTensorFunction_( ).inverse( ) * totalTorque;
    }

private:

    std::function< Eigen::Matrix3d( ) > inertiaTensorFunction_;

    SingleBodyTorqueModelMap torqueList_;

    Eigen::Vector3d effectiveTorque_;
};

BOOST_AUTO_TEST_CASE( testConstantTorquePartials )
{
    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    SystemOfBodies bodies = SystemOfBodies( "Mars", "ECLIPJ2000" );
    bodies.createEmptyBody( "Mars", false );
    bodies.at( "Mars" )->setEphemeris( std::make_shared< ephemerides::ConstantEphemeris >(
                                         [ = ]( ){ return Eigen::Vector6d::Zero( ); } ) );
    bodies.at( "Mars" )->setGravityFieldModel(
                std::make_shared< gravitation::GravityFieldModel >(
                    spice_interface::getBodyGravitationalParameter( "Mars" ) ) );
    bodies.createEmptyBody( "Phobos" );

    Eigen::Matrix3d phobosInertiaTensor = Eigen::Matrix3d::Zero( );
    phobosInertiaTensor( 0, 0 ) = 0.3615;
    phobosInertiaTensor( 1, 1 ) = 0.4265;
    phobosInertiaTensor( 2, 2 ) = 0.5024;


    phobosInertiaTensor *= ( 11.27E3 * 11.27E3 * 1.0659E16 );
    bodies.at( "Phobos" )->setBodyInertiaTensor(
                phobosInertiaTensor, ( 0.3615 + 0.4265 + 0.5024 ) / 3.0 );

    double phobosGravitationalParameter = 1.0659E16 * physical_constants::GRAVITATIONAL_CONSTANT;
    double phobosReferenceRadius = 11.27E3;

    Eigen::MatrixXd phobosCosineGravityFieldCoefficients = Eigen::MatrixXd::Zero( 6, 6 ),
            phobosSineGravityFieldCoefficients = Eigen::MatrixXd::Zero( 6, 6 );
    double phobosScaledMeanMomentOfInertia;
    gravitation::getDegreeTwoSphericalHarmonicCoefficients(
                phobosInertiaTensor, phobosGravitationalParameter, phobosReferenceRadius, true,
                phobosCosineGravityFieldCoefficients, phobosSineGravityFieldCoefficients, phobosScaledMeanMomentOfInertia );

    bodies.at( "Phobos" )->setGravityFieldModel(
                std::make_shared< gravitation::SphericalHarmonicsGravityField >(
                    phobosGravitationalParameter, phobosReferenceRadius, phobosCosineGravityFieldCoefficients,
                    phobosSineGravityFieldCoefficients, "Phobos_Fixed",
                    std::bind( &Body::setBodyInertiaTensorFromGravityFieldAndExistingMeanMoment, bodies.at( "Phobos" ), true ) ) );

    Eigen::Quaterniond noRotationQuaternion = Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) );
    Eigen::Matrix< double, 7, 1 > unitRotationState = Eigen::Matrix< double, 7, 1 >::Zero( );
    unitRotationState( 0 ) = noRotationQuaternion.w( );
    unitRotationState( 1 ) = noRotationQuaternion.x( );
    unitRotationState( 2 ) = noRotationQuaternion.y( );
    unitRotationState( 3 ) = noRotationQuaternion.z( );
    unitRotationState( 4 ) = 1.0E-5;
    unitRotationState( 5 ) = 2.0E-5;
    unitRotationState( 6 ) = -3.2E-5;

    std::map< double, Eigen::Matrix< double, 7, 1 > > dummyRotationMap;
    dummyRotationMap[ -1.0E100 ] = unitRotationState;
    dummyRotationMap[ 1.0E100 ] = unitRotationState;

    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< double, 7, 1 > > > dummyInterpolator =
            std::make_shared< interpolators::LinearInterpolator< double, Eigen::Matrix< double, 7, 1 > > >( dummyRotationMap );
    bodies.at( "Phobos" )->setRotationalEphemeris( std::make_shared< TabulatedRotationalEphemeris< double, double > >(
                                                     dummyInterpolator, "ECLIPJ2000", "Phobos_Fixed" ) );



    Eigen::Vector6d phobosKeplerElements = Eigen::Vector6d::Zero( );
    double phobosSemiMajorAxis = 9376.0E3;
    phobosKeplerElements( 0 ) = phobosSemiMajorAxis;

    bodies.at( "Phobos" )->setEphemeris( std::make_shared< ephemerides::KeplerEphemeris >(
                                           phobosKeplerElements, 0.0, spice_interface::getBodyGravitationalParameter( "Mars" ),
                                           "Mars", "ECLIPJ2000" ) );


    // Create empty bodies, phobos and mars.
    std::shared_ptr< Body > phobos = bodies.at( "Phobos" );
    std::shared_ptr< Body > mars = bodies.at( "Mars" );

    double testTime = 1000.0;
    phobos->setStateFromEphemeris( testTime );

    Eigen::Vector7d phobosRotationalState = Eigen::Vector7d::Zero( );
    phobosRotationalState.segment( 0, 4 ) = tudat::linear_algebra::convertQuaternionToVectorFormat(
                Eigen::Quaterniond( Eigen::AngleAxisd( 0.4343, Eigen::Vector3d::UnitZ( ) ) *
                                    Eigen::AngleAxisd( 2.4354, Eigen::Vector3d::UnitX( ) ) *
                                    Eigen::AngleAxisd( 1.2434, Eigen::Vector3d::UnitY( ) ) ) );
    phobosRotationalState( 4 ) = 1.0E-5;
    phobosRotationalState( 5 ) = 2.0E-5;
    phobosRotationalState( 6 ) = -3.2E-5;

    phobos->setCurrentRotationalStateToLocalFrame( phobosRotationalState );

    mars->setStateFromEphemeris( testTime );
    //    mars->setCurrentRotationalStateToLocalFrameFromEphemeris( 0.0 );


    // Create acceleration due to mars on phobos.
    std::shared_ptr< SecondDegreeGravitationalTorqueModel > gravitationalTorque =
            createSecondDegreeGravitationalTorqueModel( bodies.at( "Phobos" ), bodies.at( "Mars" ), "Phobos", "Mars" );
    std::shared_ptr< InertialTorqueModel > inertialTorqueModel =
            createInertialTorqueModel( bodies.at( "Phobos" ), "Phobos" );
    gravitationalTorque->updateMembers( 0.0 );
    inertialTorqueModel->updateMembers( 0.0 );

    SingleBodyTorqueModelMap torqueList;
    torqueList[ "Mars" ].push_back( gravitationalTorque );
    torqueList[ "Phobos" ].push_back( inertialTorqueModel );

    std::shared_ptr< EffectiveTorqueModel > effectiveTorqueModel =
            std::make_shared< EffectiveTorqueModel >(
                std::bind( &Body::getBodyInertiaTensor, bodies.at( "Phobos" ) ), torqueList );
    effectiveTorqueModel->updateMembers( 0.0 );

    // Create central gravity partial.
    std::shared_ptr< TorquePartial > torquePartial =
            createConstantTorqueRotationalDynamicsPartial(
                std::make_pair( "Phobos", bodies.at( "Phobos" ) ), torqueList );

    // Create gravitational parameter object.
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Phobos", gravitational_parameter) );
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Phobos", mean_moment_of_inertia ) );

    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                  1, 0, 3, 3, "Phobos", spherical_harmonics_cosine_coefficient_block ) );
    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                  1, 1, 3, 3, "Phobos", spherical_harmonics_sine_coefficient_block ) );

    std::shared_ptr< EstimatableParameterSet< double > > parameterSet =
            createParametersToEstimate( parameterNames, bodies );

    std::shared_ptr< EstimatableParameter< double > > phobosGravitationalParameterParameter =
            parameterSet->getEstimatedDoubleParameters( ).at( 0 );
    std::shared_ptr< EstimatableParameter< double > > meanMomentOfInertiaParameter =
            parameterSet->getEstimatedDoubleParameters( ).at( 1 );
    std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > phobosCosineCoefficientsParameter =
            parameterSet->getEstimatedVectorParameters( ).at( 0 );
    std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > phobosSineCoefficientsParameter =
            parameterSet->getEstimatedVectorParameters( ).at( 1 );

    // Calculate analytical partials.
    torquePartial->update( testTime );

    Eigen::MatrixXd partialWrtPhobosOrientation = Eigen::MatrixXd::Zero( 3, 4 );
    torquePartial->wrtOrientationOfAcceleratedBody( partialWrtPhobosOrientation.block( 0, 0, 3, 4 ) );
    Eigen::MatrixXd partialWrtPhobosRotationalVelocity = Eigen::Matrix3d::Zero( );
    torquePartial->wrtRotationalVelocityOfAcceleratedBody( partialWrtPhobosRotationalVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );
    Eigen::MatrixXd partialWrtMarsOrientation = Eigen::MatrixXd::Zero( 3, 4  );
    torquePartial->wrtOrientationOfAcceleratingBody( partialWrtMarsOrientation.block( 0, 0, 3, 4 ) );
    Eigen::MatrixXd partialWrtMarsRotationalVelocity = Eigen::Matrix3d::Zero( );
    torquePartial->wrtRotationalVelocityOfAcceleratingBody( partialWrtMarsRotationalVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );

    Eigen::MatrixXd partialWrtPhobosState = Eigen::MatrixXd::Zero( 3, 6 );
    torquePartial->wrtNonRotationalStateOfAdditionalBody(
                partialWrtPhobosState.block( 0, 0, 3, 6 ), std::make_pair( "Phobos", "" ), propagators::translational_state );
    Eigen::MatrixXd partialWrtMarsState = Eigen::MatrixXd::Zero( 3, 6 );
    torquePartial->wrtNonRotationalStateOfAdditionalBody(
                partialWrtMarsState.block( 0, 0, 3, 6 ), std::make_pair( "Mars", "" ), propagators::translational_state );

    Eigen::Vector3d partialWrtPhobosGravitationalParameter = torquePartial->wrtParameter(
                phobosGravitationalParameterParameter );
    Eigen::Vector3d partialWrtMeanMomentOfInertia = torquePartial->wrtParameter(
                meanMomentOfInertiaParameter );

    Eigen::MatrixXd partialWrtPhobosCosineCoefficients = torquePartial->wrtParameter(
                phobosCosineCoefficientsParameter );
    Eigen::MatrixXd partialWrtPhobosSineCoefficients = torquePartial->wrtParameter(
                phobosSineCoefficientsParameter );

    // Declare numerical partials.
    Eigen::Matrix3d testPartialWrtPhobosRotationalVelocity = Eigen::Matrix3d::Zero( );
    Eigen::Matrix< double, 3, 4 > testPartialWrtMarsOrientation = Eigen::Matrix< double, 3, 4 >::Zero( );
    Eigen::Matrix3d testPartialWrtMarsRotationalVelocity = Eigen::Matrix3d::Zero( );

    Eigen::Matrix3d testPartialWrtMarsPosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtMarsVelocity = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtPhobosPosition = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d testPartialWrtPhobosVelocity = Eigen::Matrix3d::Zero( );

    // Declare perturbations in orientation for numerical partial/
    Eigen::Vector4d orientationPerturbation;
    orientationPerturbation << 1.0E-9, 1.0E-9, 1.0E-9, 1.0E-9;
    Eigen::Vector3d rotationalVelocityPerturbation;
    rotationalVelocityPerturbation <<  1.0E-6, 1.0E-6, 1.0E-6;

    // Create state access/modification functions for bodies.
    std::function< void( Eigen::Vector7d ) > phobosRotationalStateSetFunction =
            std::bind( &Body::setCurrentRotationalStateToLocalFrame, phobos, std::placeholders::_1 );
    std::function< void( Eigen::Vector7d ) > marsRotationalStateSetFunction =
            std::bind( &Body::setCurrentRotationalStateToLocalFrame, mars, std::placeholders::_1 );

    //    // Calculate numerical partials.
    std::vector< Eigen::Vector4d > appliedQuaternionPerturbation;
    Eigen::MatrixXd torqueDeviations = calculateTorqueDeviationDueToOrientationChange(
                phobosRotationalStateSetFunction, effectiveTorqueModel, phobos->getRotationalStateVector( ), orientationPerturbation,
                appliedQuaternionPerturbation );

    testPartialWrtPhobosRotationalVelocity =  calculateTorqueWrtRotationalStatePartials(
                phobosRotationalStateSetFunction, effectiveTorqueModel, phobos->getRotationalStateVector( ), rotationalVelocityPerturbation, 4, 3 );
    testPartialWrtMarsOrientation =  calculateTorqueWrtRotationalStatePartials(
                marsRotationalStateSetFunction, effectiveTorqueModel, mars->getRotationalStateVector( ), orientationPerturbation, 0, 4 );
    testPartialWrtMarsRotationalVelocity =  calculateTorqueWrtRotationalStatePartials(
                marsRotationalStateSetFunction, effectiveTorqueModel, mars->getRotationalStateVector( ), rotationalVelocityPerturbation, 4, 3 );



    std::function< void( Eigen::Vector6d ) > phobosStateSetFunction =
            std::bind( &Body::setState,  phobos, std::placeholders::_1 );
    std::function< void( Eigen::Vector6d ) > marsStateSetFunction =
            std::bind( &Body::setState, mars, std::placeholders::_1 );

    Eigen::Vector3d positionPerturbation;
    positionPerturbation << 1.0, 1.0, 100.0;
    Eigen::Vector3d velocityPerturbation;
    velocityPerturbation << 1.0E-3, 1.0E-3, 1.0E-3;

    // Calculate numerical partials.
    testPartialWrtMarsPosition = calculateTorqueWrtTranslationalStatePartials(
                marsStateSetFunction, effectiveTorqueModel, mars->getState( ), positionPerturbation, 0 );
    testPartialWrtMarsVelocity = calculateTorqueWrtTranslationalStatePartials(
                marsStateSetFunction, effectiveTorqueModel, mars->getState( ), velocityPerturbation, 3 );
    testPartialWrtPhobosPosition = calculateTorqueWrtTranslationalStatePartials(
                phobosStateSetFunction, effectiveTorqueModel, phobos->getState( ), positionPerturbation, 0 );
    testPartialWrtPhobosVelocity = calculateTorqueWrtTranslationalStatePartials(
                phobosStateSetFunction, effectiveTorqueModel, phobos->getState( ), velocityPerturbation, 3 );


    std::function< void( ) > updateFunction = &emptyFunction;
           // std::bind( &Body::setBodyInertiaTensorFromGravityFieldAndExistingMeanMoment, bodies.at( "Phobos" ), true );
    Eigen::Vector3d testPartialWrtPhobosGravitationalParameter = calculateTorqueWrtParameterPartials(
                phobosGravitationalParameterParameter, effectiveTorqueModel, 1.0E2, updateFunction );
    Eigen::MatrixXd testPartialWrtMeanMomentOfInertia = calculateTorqueWrtParameterPartials(
                meanMomentOfInertiaParameter, effectiveTorqueModel, 1.0E-4, updateFunction );
    Eigen::MatrixXd testPartialWrtPhobosCosineCoefficients = calculateTorqueWrtParameterPartials(
                phobosCosineCoefficientsParameter, effectiveTorqueModel,
                Eigen::VectorXd::Constant( phobosCosineCoefficientsParameter->getParameterSize( ), 1.0E-6 ), updateFunction );
    Eigen::MatrixXd testPartialWrtPhobosSineCoefficients = calculateTorqueWrtParameterPartials(
                phobosSineCoefficientsParameter, effectiveTorqueModel,
                Eigen::VectorXd::Constant( phobosSineCoefficientsParameter->getParameterSize( ), 1.0E-6 ), updateFunction );

    //    // Compare numerical and analytical results.
    for( int index = 1; index < 4; index++ )
    {
        Eigen::Vector3d numericalChangeInTorque = torqueDeviations.block( 0, index - 1, 3, 1 );
        Eigen::Vector3d analyticalChangeInTorque =
                partialWrtPhobosOrientation.block( 0, 0, 3, 1 ) * appliedQuaternionPerturbation[ index ]( 0 ) +
                partialWrtPhobosOrientation.block( 0, index, 3, 1 ) * appliedQuaternionPerturbation[ index ]( index );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( numericalChangeInTorque,
                                           analyticalChangeInTorque, std::numeric_limits< double >::epsilon( ) );
    }

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPhobosRotationalVelocity,
                                       partialWrtPhobosRotationalVelocity, 1.0E-10 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtMarsOrientation,
                                       partialWrtMarsOrientation, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtMarsRotationalVelocity,
                                       partialWrtMarsRotationalVelocity, std::numeric_limits< double >::epsilon( ) );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(  partialWrtPhobosState.block( 0, 0, 3, 3 ),
                                        testPartialWrtPhobosPosition, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(  partialWrtMarsState.block( 0, 0, 3, 3 ),
                                        testPartialWrtMarsPosition, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(  partialWrtPhobosState.block( 0, 3, 3, 3 ),
                                        testPartialWrtPhobosVelocity, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(  partialWrtMarsState.block( 0, 3, 3, 3 ),
                                        testPartialWrtMarsVelocity, std::numeric_limits< double >::epsilon( ) );

    Eigen::Vector3d sclaedPartialWrtPhobosGravitationalParameter = phobosInertiaTensor * testPartialWrtPhobosGravitationalParameter;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( sclaedPartialWrtPhobosGravitationalParameter,
                                       partialWrtPhobosGravitationalParameter, 1.0E-6 );

    Eigen::Vector3d scaledPartialWrtMeanMomentOfInertia = phobosInertiaTensor * testPartialWrtMeanMomentOfInertia;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( scaledPartialWrtMeanMomentOfInertia,
                                       partialWrtMeanMomentOfInertia, 1.0E-6 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtPhobosCosineCoefficients.block( 0, 0, 3, 2 ),
                                       Eigen::MatrixXd::Zero( 3, 2 ), std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtPhobosCosineCoefficients.block( 0, 5, 3, 4 ),
                                       Eigen::MatrixXd::Zero( 3, 4 ), std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtPhobosSineCoefficients.block( 0, 0, 3, 1 ),
                                       Eigen::MatrixXd::Zero( 3, 1 ), std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtPhobosSineCoefficients.block( 0, 3, 3, 3 ),
                                       Eigen::MatrixXd::Zero( 3, 3 ), std::numeric_limits< double >::epsilon( ) );

    // Check derivative of z-component w.r.t. C20 separately: value is slightly non-zero due to rounding error
    Eigen::MatrixXd scaledPartialWrtPhobosCosineCoefficients = phobosInertiaTensor * testPartialWrtPhobosCosineCoefficients;
    BOOST_CHECK_SMALL( std::fabs( scaledPartialWrtPhobosCosineCoefficients( 2, 4 ) ), 1.0E5 );
    scaledPartialWrtPhobosCosineCoefficients( 2, 4 ) = 0.0;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( scaledPartialWrtPhobosCosineCoefficients, partialWrtPhobosCosineCoefficients, 1.0E-9 );

    Eigen::MatrixXd scaledPartialWrtPhobosSineCoefficients = phobosInertiaTensor * testPartialWrtPhobosSineCoefficients;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( scaledPartialWrtPhobosSineCoefficients, partialWrtPhobosSineCoefficients, 1.0E-9 );

}

BOOST_AUTO_TEST_SUITE_END( )


} // namespace unit_tests

} // namespace tudat

