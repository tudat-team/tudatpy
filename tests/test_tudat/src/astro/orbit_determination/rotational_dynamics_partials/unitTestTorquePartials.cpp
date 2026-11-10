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

#include <limits>
#include <string>
#include "tudat/basics/testMacros.h"
#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/astro/basic_astro/unitConversions.h"

#include <boost/test/included/unit_test.hpp>

#include "tudat/astro/basic_astro/sphericalStateConversions.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/astro/ephemerides/keplerEphemeris.h"
#include "tudat/astro/relativity/metric.h"
#include "tudat/astro/orbit_determination/acceleration_partials/numericalAccelerationPartial.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/gravitationalParameter.h"
#include "tudat/simulation/estimation_setup/createTorquePartials.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/propagation_setup/createTorqueModel.h"
#include "tudat/simulation/estimation_setup/createEstimatableParametersFactory.h"
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

    // Second-degree and spherical harmonic torque, with the rotational state given by a quaternion and by its negative
    for( int test = 0; test < 4; test++ )
    {
        const int torqueModelTest = test % 2;
        const double quaternionSign = ( test < 2 ) ? 1.0 : -1.0;

        // Create bodies
        SystemOfBodies bodies = SystemOfBodies( "Mars", "ECLIPJ2000" );
        bodies.createEmptyBody( "Mars", false );
        bodies.at( "Mars" )->setEphemeris(
                std::make_shared< ephemerides::ConstantEphemeris >( [ = ]( ) { return Eigen::Vector6d::Zero( ); } ) );
        bodies.at( "Mars" )->setGravityFieldModel(
                std::make_shared< gravitation::GravityFieldModel >( spice_interface::getBodyGravitationalParameter( "Mars" ) ) );
        bodies.createEmptyBody( "Phobos" );
        std::shared_ptr< Body > phobos = bodies.at( "Phobos" );
        std::shared_ptr< Body > mars = bodies.at( "Mars" );

        // Set inertia tensor
        Eigen::Matrix3d phobosInertiaTensor = Eigen::Matrix3d::Zero( );
        phobosInertiaTensor( 0, 0 ) = 0.3615;
        phobosInertiaTensor( 1, 1 ) = 0.4265;
        phobosInertiaTensor( 2, 2 ) = 0.5024;
        double scaledMeanMomentOfInertia =
                ( phobosInertiaTensor( 0, 0 ) + phobosInertiaTensor( 1, 1 ) + phobosInertiaTensor( 2, 2 ) ) / 3.0;

        // Create gravity field
        double phobosMass = 1.0659E16;
        double phobosGravitationalParameter = phobosMass * physical_constants::GRAVITATIONAL_CONSTANT;
        double phobosReferenceRadius = 11.27E3;

        phobosInertiaTensor *= ( phobosMass * phobosReferenceRadius * phobosReferenceRadius );

        Eigen::MatrixXd phobosCosineGravityFieldCoefficients = Eigen::MatrixXd::Zero( 6, 6 ),
                        phobosSineGravityFieldCoefficients = Eigen::MatrixXd::Zero( 6, 6 );
        gravitation::getDegreeTwoSphericalHarmonicCoefficients( phobosInertiaTensor,
                                                                phobosGravitationalParameter,
                                                                phobosReferenceRadius,
                                                                true,
                                                                phobosCosineGravityFieldCoefficients,
                                                                phobosSineGravityFieldCoefficients,
                                                                scaledMeanMomentOfInertia );
        const std::shared_ptr< gravitation::SphericalHarmonicsGravityField > phobosGravityField =
                std::make_shared< gravitation::SphericalHarmonicsGravityField >( phobosGravitationalParameter,
                                                                                 phobosReferenceRadius,
                                                                                 phobosCosineGravityFieldCoefficients,
                                                                                 phobosSineGravityFieldCoefficients,
                                                                                 "Phobos_Fixed" );
        bodies.at( "Phobos" )->setGravityFieldModel( phobosGravityField );
        bodies.at( "Phobos" )
                ->setMassProperties( std::make_shared< simulation_setup::FromGravityFieldRigidBodyProperties >(
                        phobosGravityField, scaledMeanMomentOfInertia ) );
        double testTime = 1000.0;
        bodies.at( "Phobos" )->getMassProperties( )->update( testTime );
        std::cout << bodies.at( "Phobos" )->getBodyInertiaTensor( ) << std::endl;

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
        bodies.at( "Phobos" )
                ->setRotationalEphemeris( std::make_shared< TabulatedRotationalEphemeris< double, double > >(
                        dummyInterpolator, "ECLIPJ2000", "Phobos_Fixed" ) );

        Eigen::Vector6d phobosKeplerElements = Eigen::Vector6d::Zero( );
        double phobosSemiMajorAxis = 9376.0E3;
        phobosKeplerElements( 0 ) = phobosSemiMajorAxis;
        bodies.at( "Phobos" )
                ->setEphemeris( std::make_shared< ephemerides::KeplerEphemeris >(
                        phobosKeplerElements, 0.0, spice_interface::getBodyGravitationalParameter( "Mars" ), "Mars", "ECLIPJ2000" ) );

        // Update Phobos and Mars to current state
        phobos->setStateFromEphemeris( testTime );
        phobos->getMassProperties( )->update( testTime );
        mars->setStateFromEphemeris( testTime );

        // Set new Phobos rotational state
        Eigen::Vector7d phobosRotationalState = Eigen::Vector7d::Zero( );
        phobosRotationalState.segment( 0, 4 ) = tudat::linear_algebra::convertQuaternionToVectorFormat( Eigen::Quaterniond(
                Eigen::AngleAxisd( 0.4343, Eigen::Vector3d::UnitZ( ) ) * Eigen::AngleAxisd( 2.4354, Eigen::Vector3d::UnitX( ) ) *
                Eigen::AngleAxisd( 1.2434, Eigen::Vector3d::UnitY( ) ) ) );
        phobosRotationalState.segment( 0, 4 ) *= quaternionSign;
        phobos->setCurrentRotationalStateToLocalFrame( phobosRotationalState );

        // Create torque due to mars on phobos.
        std::shared_ptr< TorqueModel > gravitationalTorque;
        if( torqueModelTest == 0 )
        {
            gravitationalTorque =
                    createSecondDegreeGravitationalTorqueModel( bodies.at( "Phobos" ), bodies.at( "Mars" ), "Phobos", "Mars" );
        }
        else if( torqueModelTest == 1 )
        {
            gravitationalTorque =
                    createSphericalHarmonicGravitationalTorqueModel( bodies.at( "Phobos" ),
                                                                     bodies.at( "Mars" ),
                                                                     std::make_shared< SphericalHarmonicTorqueSettings >( 2, 2 ),
                                                                     "Phobos",
                                                                     "Mars" );
        }

        // Create parameter objects
        std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Mars", gravitational_parameter ) );
        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Phobos", mean_moment_of_inertia ) );
        parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                1, 0, 3, 3, "Phobos", spherical_harmonics_cosine_coefficient_block ) );
        parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                1, 1, 3, 3, "Phobos", spherical_harmonics_sine_coefficient_block ) );
        std::shared_ptr< EstimatableParameterSet< double > > parameterSet = createParametersToEstimate( parameterNames, bodies );

        std::shared_ptr< EstimatableParameter< double > > marsGravitationalParameterParameter =
                parameterSet->getEstimatedDoubleParameters( ).at( 0 );
        std::shared_ptr< EstimatableParameter< double > > phobosMeanMomentOfInertia = parameterSet->getEstimatedDoubleParameters( ).at( 1 );
        std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > phobosCosineCoefficientsParameter =
                parameterSet->getEstimatedVectorParameters( ).at( 0 );
        std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > phobosSineCoefficientsParameter =
                parameterSet->getEstimatedVectorParameters( ).at( 1 );

        // Create torque partial.
        std::shared_ptr< TorquePartial > torquePartial = createAnalyticalTorquePartial( gravitationalTorque,
                                                                                        std::make_pair( "Phobos", phobos ),
                                                                                        std::make_pair( "Mars", mars ),
                                                                                        basic_astrodynamics::SingleBodyTorqueModelMap( ),
                                                                                        bodies,
                                                                                        parameterSet );

        // Calculate analytical partials.
        torquePartial->update( testTime );

        Eigen::MatrixXd partialWrtPhobosOrientation = Eigen::MatrixXd::Zero( 3, 4 );
        torquePartial->wrtOrientationOfAcceleratedBody( partialWrtPhobosOrientation.block( 0, 0, 3, 4 ) );

        Eigen::MatrixXd partialWrtPhobosRotationalVelocity = Eigen::Matrix3d::Zero( );
        torquePartial->wrtRotationalVelocityOfAcceleratedBody( partialWrtPhobosRotationalVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );

        Eigen::MatrixXd partialWrtMarsOrientation = Eigen::MatrixXd::Zero( 3, 4 );
        torquePartial->wrtOrientationOfAcceleratingBody( partialWrtMarsOrientation.block( 0, 0, 3, 4 ) );

        Eigen::MatrixXd partialWrtMarsRotationalVelocity = Eigen::Matrix3d::Zero( );
        torquePartial->wrtRotationalVelocityOfAcceleratingBody( partialWrtMarsRotationalVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );

        Eigen::MatrixXd partialWrtPhobosState = Eigen::MatrixXd::Zero( 3, 6 );
        torquePartial->wrtNonRotationalStateOfAdditionalBody(
                partialWrtPhobosState.block( 0, 0, 3, 6 ), std::make_pair( "Phobos", "" ), propagators::translational_state );
        Eigen::MatrixXd partialWrtMarsState = Eigen::MatrixXd::Zero( 3, 6 );
        torquePartial->wrtNonRotationalStateOfAdditionalBody(
                partialWrtMarsState.block( 0, 0, 3, 6 ), std::make_pair( "Mars", "" ), propagators::translational_state );

        Eigen::Vector3d partialWrtMarsGravitationalParameter = torquePartial->wrtParameter( marsGravitationalParameterParameter );
        Eigen::Vector3d partialWrtMeanMomentOfInertia = torquePartial->wrtParameter( phobosMeanMomentOfInertia );

        Eigen::MatrixXd partialWrtPhobosCosineCoefficients = torquePartial->wrtParameter( phobosCosineCoefficientsParameter );
        Eigen::MatrixXd partialWrtPhobosSineCoefficients = torquePartial->wrtParameter( phobosSineCoefficientsParameter );

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
        rotationalVelocityPerturbation << 1.0E-6, 1.0E-6, 1.0E-6;

        // Create state access/modification functions for bodies.
        std::function< void( Eigen::Vector7d ) > phobosRotationalStateSetFunction =
                std::bind( &Body::setCurrentRotationalStateToLocalFrame, phobos, std::placeholders::_1 );
        std::function< void( Eigen::Vector7d ) > marsRotationalStateSetFunction =
                std::bind( &Body::setCurrentRotationalStateToLocalFrame, mars, std::placeholders::_1 );

        // Calculate numerical partials wrt rotational state.
        std::vector< Eigen::Vector4d > appliedQuaternionPerturbation;
        Eigen::MatrixXd torqueDeviations = calculateTorqueDeviationDueToOrientationChange( phobosRotationalStateSetFunction,
                                                                                           gravitationalTorque,
                                                                                           phobos->getRotationalStateVector( ),
                                                                                           orientationPerturbation,
                                                                                           appliedQuaternionPerturbation );

        testPartialWrtPhobosRotationalVelocity = calculateTorqueWrtRotationalStatePartials( phobosRotationalStateSetFunction,
                                                                                            gravitationalTorque,
                                                                                            phobos->getRotationalStateVector( ),
                                                                                            rotationalVelocityPerturbation,
                                                                                            4,
                                                                                            3 );
        testPartialWrtMarsOrientation = calculateTorqueWrtRotationalStatePartials(
                marsRotationalStateSetFunction, gravitationalTorque, mars->getRotationalStateVector( ), orientationPerturbation, 0, 4 );
        testPartialWrtMarsRotationalVelocity = calculateTorqueWrtRotationalStatePartials( marsRotationalStateSetFunction,
                                                                                          gravitationalTorque,
                                                                                          mars->getRotationalStateVector( ),
                                                                                          rotationalVelocityPerturbation,
                                                                                          4,
                                                                                          3 );

        std::function< void( Eigen::Vector6d ) > phobosStateSetFunction = std::bind( &Body::setState, phobos, std::placeholders::_1 );
        std::function< void( Eigen::Vector6d ) > marsStateSetFunction = std::bind( &Body::setState, mars, std::placeholders::_1 );

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
        Eigen::Vector3d testPartialWrtMarsGravitationalParameter =
                calculateTorqueWrtParameterPartials( marsGravitationalParameterParameter, gravitationalTorque, 1.0E12, [ = ]( ) {
                    mars->getMassProperties( )->update( testTime );
                } );
        Eigen::Vector3d testPartialWrtMeanMomentOfInertia =
                calculateTorqueWrtParameterPartials( phobosMeanMomentOfInertia, gravitationalTorque, 1.0E-1 );
        std::function< void( ) > updateFunction =
                std::bind( &RigidBodyProperties::update, bodies.at( "Phobos" )->getMassProperties( ), testTime );
        // std::bind( &Body::setBodyInertiaTensorFromGravityFieldAndExistingMeanMoment, bodies.at( "Phobos" ), true );
        Eigen::MatrixXd testPartialWrtPhobosCosineCoefficients = calculateTorqueWrtParameterPartials(
                phobosCosineCoefficientsParameter,
                gravitationalTorque,
                Eigen::VectorXd::Constant( phobosCosineCoefficientsParameter->getParameterSize( ), 1.0E-6 ),
                updateFunction );
        Eigen::MatrixXd testPartialWrtPhobosSineCoefficients = calculateTorqueWrtParameterPartials(
                phobosSineCoefficientsParameter,
                gravitationalTorque,
                Eigen::VectorXd::Constant( phobosSineCoefficientsParameter->getParameterSize( ), 1.0E-6 ),
                updateFunction );

        // Compare numerical and analytical results.
        for( int index = 1; index < 4; index++ )
        {
            Eigen::Vector3d numericalChangeInTorque = torqueDeviations.block( 0, index - 1, 3, 1 );
            Eigen::Vector3d analyticalChangeInTorque =
                    partialWrtPhobosOrientation.block( 0, 0, 3, 1 ) * appliedQuaternionPerturbation[ index ]( 0 ) +
                    partialWrtPhobosOrientation.block( 0, index, 3, 1 ) * appliedQuaternionPerturbation[ index ]( index );

            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( numericalChangeInTorque, analyticalChangeInTorque, 1.0E-6 );
        }

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                testPartialWrtPhobosRotationalVelocity, partialWrtPhobosRotationalVelocity, std::numeric_limits< double >::epsilon( ) );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                testPartialWrtMarsOrientation, partialWrtMarsOrientation, std::numeric_limits< double >::epsilon( ) );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                testPartialWrtMarsRotationalVelocity, partialWrtMarsRotationalVelocity, std::numeric_limits< double >::epsilon( ) );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtPhobosState.block( 0, 0, 3, 3 ), testPartialWrtPhobosPosition, 1.0E-7 );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtMarsState.block( 0, 0, 3, 3 ), testPartialWrtMarsPosition, 1.0E-7 );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                partialWrtPhobosState.block( 0, 3, 3, 3 ), testPartialWrtPhobosVelocity, std::numeric_limits< double >::epsilon( ) );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                partialWrtMarsState.block( 0, 3, 3, 3 ), testPartialWrtMarsVelocity, std::numeric_limits< double >::epsilon( ) );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtMarsGravitationalParameter, partialWrtMarsGravitationalParameter, 1.0E-6 );

        // Check partial w.r.t. Phobos gravitational parameter, to which its gravity-derived inertia tensor is proportional
        std::shared_ptr< EstimatableParameter< double > > phobosGravitationalParameterParameter =
                std::make_shared< GravitationalParameter >( phobos->getGravityFieldModel( ), "Phobos" );
        const Eigen::Vector3d partialWrtPhobosGravitationalParameter = torquePartial->wrtParameter( phobosGravitationalParameterParameter );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                calculateTorqueWrtParameterPartials( phobosGravitationalParameterParameter, gravitationalTorque, 1.0E2, updateFunction ),
                partialWrtPhobosGravitationalParameter,
                1.0E-6 );

        // Check derivative of z-component w.r.t. C20 separately: value is slightly non-zero due to rounding error
        BOOST_CHECK_SMALL( std::fabs( testPartialWrtPhobosCosineCoefficients( 2, 2 ) ), 2.0E6 );
        testPartialWrtPhobosCosineCoefficients( 2, 2 ) = 0.0;

        if( torqueModelTest == 1 )
        {
            BOOST_CHECK_SMALL( std::fabs( partialWrtPhobosCosineCoefficients( 2, 2 ) ), 2.0E6 );
            partialWrtPhobosCosineCoefficients( 2, 2 ) = 0.0;

            BOOST_CHECK_SMALL( std::fabs( testPartialWrtPhobosCosineCoefficients( 0, 1 ) ), 2.0E6 );
            testPartialWrtPhobosCosineCoefficients( 0, 1 ) = 0.0;

            BOOST_CHECK_SMALL( std::fabs( partialWrtPhobosCosineCoefficients( 0, 1 ) ), 2.0E6 );
            partialWrtPhobosCosineCoefficients( 0, 1 ) = 0.0;

            BOOST_CHECK_SMALL( std::fabs( testPartialWrtPhobosSineCoefficients( 1, 0 ) ), 2.0E6 );
            testPartialWrtPhobosSineCoefficients( 1, 0 ) = 0.0;

            BOOST_CHECK_SMALL( std::fabs( partialWrtPhobosSineCoefficients( 1, 0 ) ), 2.0E6 );
            partialWrtPhobosSineCoefficients( 1, 0 ) = 0.0;

            BOOST_CHECK_SMALL( std::fabs( testPartialWrtPhobosCosineCoefficients( 2, 0 ) ), 2.0E6 );
            testPartialWrtPhobosCosineCoefficients( 2, 0 ) = 0.0;

            BOOST_CHECK_SMALL( std::fabs( partialWrtPhobosCosineCoefficients( 2, 0 ) ), 2.0E6 );
            partialWrtPhobosCosineCoefficients( 2, 0 ) = 0.0;
        }

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPhobosCosineCoefficients, partialWrtPhobosCosineCoefficients, 1.0E-9 );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPhobosSineCoefficients, partialWrtPhobosSineCoefficients, 1.0E-9 );

        for( int i = 0; i < 3; i++ )
        {
            BOOST_CHECK_SMALL( std::fabs( partialWrtMeanMomentOfInertia( i, 0 ) - testPartialWrtMeanMomentOfInertia( i, 0 ) ), 1.0E2 );
        }
    }
}

BOOST_AUTO_TEST_CASE( testInertialTorquePartials )
{
    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Create bodies
    SystemOfBodies bodies = SystemOfBodies( "Mars", "ECLIPJ2000" );
    bodies.createEmptyBody( "Mars", false );
    bodies.at( "Mars" )->setEphemeris(
            std::make_shared< ephemerides::ConstantEphemeris >( [ = ]( ) { return Eigen::Vector6d::Zero( ); } ) );
    bodies.at( "Mars" )->setGravityFieldModel(
            std::make_shared< gravitation::GravityFieldModel >( spice_interface::getBodyGravitationalParameter( "Mars" ) ) );
    bodies.createEmptyBody( "Phobos" );
    std::shared_ptr< Body > phobos = bodies.at( "Phobos" );
    std::shared_ptr< Body > mars = bodies.at( "Mars" );

    // Set inertia tensor
    Eigen::Matrix3d phobosInertiaTensor = Eigen::Matrix3d::Zero( );
    phobosInertiaTensor( 0, 0 ) = 0.3615;
    phobosInertiaTensor( 1, 1 ) = 0.4265;
    phobosInertiaTensor( 2, 2 ) = 0.5024;
    double scaledMeanMomentOfInertia = ( phobosInertiaTensor( 0, 0 ) + phobosInertiaTensor( 1, 1 ) + phobosInertiaTensor( 2, 2 ) ) / 3.0;

    // Create gravity field
    double phobosMass = 1.0659E16;
    double phobosGravitationalParameter = phobosMass * physical_constants::GRAVITATIONAL_CONSTANT;
    double phobosReferenceRadius = 11.27E3;

    Eigen::MatrixXd phobosCosineGravityFieldCoefficients = Eigen::MatrixXd::Zero( 6, 6 ),
                    phobosSineGravityFieldCoefficients = Eigen::MatrixXd::Zero( 6, 6 );
    gravitation::getDegreeTwoSphericalHarmonicCoefficients( phobosInertiaTensor,
                                                            phobosGravitationalParameter,
                                                            phobosReferenceRadius,
                                                            true,
                                                            phobosCosineGravityFieldCoefficients,
                                                            phobosSineGravityFieldCoefficients,
                                                            scaledMeanMomentOfInertia );

    const std::shared_ptr< gravitation::SphericalHarmonicsGravityField > phobosGravityField =
            std::make_shared< gravitation::SphericalHarmonicsGravityField >( phobosGravitationalParameter,
                                                                             phobosReferenceRadius,
                                                                             phobosCosineGravityFieldCoefficients,
                                                                             phobosSineGravityFieldCoefficients,
                                                                             "Phobos_Fixed" );
    bodies.at( "Phobos" )->setGravityFieldModel( phobosGravityField );
    bodies.at( "Phobos" )
            ->setMassProperties( std::make_shared< simulation_setup::FromGravityFieldRigidBodyProperties >( phobosGravityField,
                                                                                                            scaledMeanMomentOfInertia ) );

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
    bodies.at( "Phobos" )
            ->setRotationalEphemeris(
                    std::make_shared< TabulatedRotationalEphemeris< double, double > >( dummyInterpolator, "ECLIPJ2000", "Phobos_Fixed" ) );

    Eigen::Vector6d phobosKeplerElements = Eigen::Vector6d::Zero( );
    double phobosSemiMajorAxis = 9376.0E3;
    phobosKeplerElements( 0 ) = phobosSemiMajorAxis;
    bodies.at( "Phobos" )
            ->setEphemeris( std::make_shared< ephemerides::KeplerEphemeris >(
                    phobosKeplerElements, 0.0, spice_interface::getBodyGravitationalParameter( "Mars" ), "Mars", "ECLIPJ2000" ) );

    // Update Phobos and Mars to current state
    double testTime = 1000.0;
    phobos->setStateFromEphemeris( testTime );
    mars->setStateFromEphemeris( testTime );

    // Set new Phobos rotational state
    Eigen::Vector7d phobosRotationalState = Eigen::Vector7d::Zero( );
    phobosRotationalState.segment( 0, 4 ) = tudat::linear_algebra::convertQuaternionToVectorFormat( Eigen::Quaterniond(
            Eigen::AngleAxisd( 0.4343, Eigen::Vector3d::UnitZ( ) ) * Eigen::AngleAxisd( 2.4354, Eigen::Vector3d::UnitX( ) ) *
            Eigen::AngleAxisd( 1.2434, Eigen::Vector3d::UnitY( ) ) ) );
    phobosRotationalState( 4 ) = 1.0E-5;
    phobosRotationalState( 5 ) = 2.0E-5;
    phobosRotationalState( 6 ) = -3.2E-5;
    phobos->setCurrentRotationalStateToLocalFrame( phobosRotationalState );

    // Create torque due to mars on phobos.
    std::shared_ptr< InertialTorqueModel > inertialTorqueModel = createInertialTorqueModel( bodies.at( "Phobos" ), "Phobos" );
    inertialTorqueModel->updateMembers( 0.0 );

    SingleBodyTorqueModelMap torqueList;
    torqueList[ "Phobos" ].push_back( inertialTorqueModel );

    // Create shared_ptr partial.
    std::shared_ptr< TorquePartial > torquePartial =
            createAnalyticalTorquePartial( inertialTorqueModel, std::make_pair( "Phobos", phobos ), std::make_pair( "Phobos", phobos ) );

    // Create parameter objects.
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Phobos", gravitational_parameter ) );
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Phobos", mean_moment_of_inertia ) );

    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
            1, 0, 3, 3, "Phobos", spherical_harmonics_cosine_coefficient_block ) );
    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
            1, 1, 3, 3, "Phobos", spherical_harmonics_sine_coefficient_block ) );

    std::shared_ptr< EstimatableParameterSet< double > > parameterSet = createParametersToEstimate( parameterNames, bodies );

    std::shared_ptr< EstimatableParameter< double > > phobosGravitationalParameterParameter =
            parameterSet->getEstimatedDoubleParameters( ).at( 0 );
    std::shared_ptr< EstimatableParameter< double > > meanMomentOfInertiaParameter = parameterSet->getEstimatedDoubleParameters( ).at( 1 );
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
    Eigen::MatrixXd partialWrtMarsOrientation = Eigen::MatrixXd::Zero( 3, 4 );
    torquePartial->wrtOrientationOfAcceleratingBody( partialWrtMarsOrientation.block( 0, 0, 3, 4 ) );
    Eigen::MatrixXd partialWrtMarsRotationalVelocity = Eigen::Matrix3d::Zero( );
    torquePartial->wrtRotationalVelocityOfAcceleratingBody( partialWrtMarsRotationalVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );

    Eigen::MatrixXd partialWrtPhobosState = Eigen::MatrixXd::Zero( 3, 6 );
    torquePartial->wrtNonRotationalStateOfAdditionalBody(
            partialWrtPhobosState.block( 0, 0, 3, 6 ), std::make_pair( "Phobos", "" ), propagators::translational_state );
    Eigen::MatrixXd partialWrtMarsState = Eigen::MatrixXd::Zero( 3, 6 );
    torquePartial->wrtNonRotationalStateOfAdditionalBody(
            partialWrtMarsState.block( 0, 0, 3, 6 ), std::make_pair( "Mars", "" ), propagators::translational_state );

    Eigen::Vector3d partialWrtPhobosGravitationalParameter = torquePartial->wrtParameter( phobosGravitationalParameterParameter );
    Eigen::Vector3d partialWrtMeanMomentOfInertia = torquePartial->wrtParameter( meanMomentOfInertiaParameter );

    Eigen::MatrixXd partialWrtPhobosCosineCoefficients = torquePartial->wrtParameter( phobosCosineCoefficientsParameter );
    Eigen::MatrixXd partialWrtPhobosSineCoefficients = torquePartial->wrtParameter( phobosSineCoefficientsParameter );

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
    rotationalVelocityPerturbation << 1.0E-6, 1.0E-6, 1.0E-6;

    // Create state access/modification functions for bodies.
    std::function< void( Eigen::Vector7d ) > phobosRotationalStateSetFunction =
            std::bind( &Body::setCurrentRotationalStateToLocalFrame, phobos, std::placeholders::_1 );
    std::function< void( Eigen::Vector7d ) > marsRotationalStateSetFunction =
            std::bind( &Body::setCurrentRotationalStateToLocalFrame, mars, std::placeholders::_1 );

    // Calculate numerical partials.
    std::vector< Eigen::Vector4d > appliedQuaternionPerturbation;
    Eigen::MatrixXd torqueDeviations = calculateTorqueDeviationDueToOrientationChange( phobosRotationalStateSetFunction,
                                                                                       inertialTorqueModel,
                                                                                       phobos->getRotationalStateVector( ),
                                                                                       orientationPerturbation,
                                                                                       appliedQuaternionPerturbation );

    testPartialWrtPhobosRotationalVelocity = calculateTorqueWrtRotationalStatePartials( phobosRotationalStateSetFunction,
                                                                                        inertialTorqueModel,
                                                                                        phobos->getRotationalStateVector( ),
                                                                                        rotationalVelocityPerturbation,
                                                                                        4,
                                                                                        3 );
    testPartialWrtMarsOrientation = calculateTorqueWrtRotationalStatePartials(
            marsRotationalStateSetFunction, inertialTorqueModel, mars->getRotationalStateVector( ), orientationPerturbation, 0, 4 );
    testPartialWrtMarsRotationalVelocity = calculateTorqueWrtRotationalStatePartials(
            marsRotationalStateSetFunction, inertialTorqueModel, mars->getRotationalStateVector( ), rotationalVelocityPerturbation, 4, 3 );

    std::function< void( Eigen::Vector6d ) > phobosStateSetFunction = std::bind( &Body::setState, phobos, std::placeholders::_1 );
    std::function< void( Eigen::Vector6d ) > marsStateSetFunction = std::bind( &Body::setState, mars, std::placeholders::_1 );

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

    std::function< void( ) > updateFunction =
            std::bind( &RigidBodyProperties::update, bodies.at( "Phobos" )->getMassProperties( ), testTime );
    //            std::bind( &Body::setBodyInertiaTensorFromGravityFieldAndExistingMeanMoment, bodies.at( "Phobos" ), true );
    Eigen::Vector3d testPartialWrtPhobosGravitationalParameter =
            calculateTorqueWrtParameterPartials( phobosGravitationalParameterParameter, inertialTorqueModel, 1.0E8, updateFunction );
    Eigen::MatrixXd testPartialWrtMeanMomentOfInertia =
            calculateTorqueWrtParameterPartials( meanMomentOfInertiaParameter, inertialTorqueModel, 1.0E-1, updateFunction );
    Eigen::MatrixXd testPartialWrtPhobosCosineCoefficients = calculateTorqueWrtParameterPartials(
            phobosCosineCoefficientsParameter,
            inertialTorqueModel,
            Eigen::VectorXd::Constant( phobosCosineCoefficientsParameter->getParameterSize( ), 1.0E-6 ),
            updateFunction );
    Eigen::MatrixXd testPartialWrtPhobosSineCoefficients =
            calculateTorqueWrtParameterPartials( phobosSineCoefficientsParameter,
                                                 inertialTorqueModel,
                                                 Eigen::VectorXd::Constant( phobosSineCoefficientsParameter->getParameterSize( ), 1.0E-6 ),
                                                 updateFunction );

    //    // Compare numerical and analytical results.
    for( int index = 1; index < 4; index++ )
    {
        Eigen::Vector3d numericalChangeInTorque = torqueDeviations.block( 0, index - 1, 3, 1 );
        Eigen::Vector3d analyticalChangeInTorque =
                partialWrtPhobosOrientation.block( 0, 0, 3, 1 ) * appliedQuaternionPerturbation[ index ]( 0 ) +
                partialWrtPhobosOrientation.block( 0, index, 3, 1 ) * appliedQuaternionPerturbation[ index ]( index );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( numericalChangeInTorque, analyticalChangeInTorque, std::numeric_limits< double >::epsilon( ) );
    }

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPhobosRotationalVelocity, partialWrtPhobosRotationalVelocity, 1.0E-10 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            testPartialWrtMarsOrientation, partialWrtMarsOrientation, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            testPartialWrtMarsRotationalVelocity, partialWrtMarsRotationalVelocity, std::numeric_limits< double >::epsilon( ) );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            partialWrtPhobosState.block( 0, 0, 3, 3 ), testPartialWrtPhobosPosition, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            partialWrtMarsState.block( 0, 0, 3, 3 ), testPartialWrtMarsPosition, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            partialWrtPhobosState.block( 0, 3, 3, 3 ), testPartialWrtPhobosVelocity, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            partialWrtMarsState.block( 0, 3, 3, 3 ), testPartialWrtMarsVelocity, std::numeric_limits< double >::epsilon( ) );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPhobosGravitationalParameter, partialWrtPhobosGravitationalParameter, 1.0E-6 );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtPhobosCosineCoefficients.block( 0, 0, 3, 2 ),
                                       Eigen::MatrixXd::Zero( 3, 2 ),
                                       std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtPhobosCosineCoefficients.block( 0, 5, 3, 4 ),
                                       Eigen::MatrixXd::Zero( 3, 4 ),
                                       std::numeric_limits< double >::epsilon( ) );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtPhobosSineCoefficients.block( 0, 0, 3, 1 ),
                                       Eigen::MatrixXd::Zero( 3, 1 ),
                                       std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtPhobosSineCoefficients.block( 0, 3, 3, 3 ),
                                       Eigen::MatrixXd::Zero( 3, 3 ),
                                       std::numeric_limits< double >::epsilon( ) );

    // Check derivative of z-component w.r.t. C20 separately: value is slightly non-zero due to rounding error
    BOOST_CHECK_SMALL( std::fabs( testPartialWrtPhobosCosineCoefficients( 2, 2 ) ), 1.0E5 );
    testPartialWrtPhobosCosineCoefficients( 2, 2 ) = 0.0;

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPhobosCosineCoefficients, partialWrtPhobosCosineCoefficients, 1.0E-9 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPhobosSineCoefficients, partialWrtPhobosSineCoefficients, 1.0E-9 );

    for( int i = 0; i < 3; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( partialWrtMeanMomentOfInertia( i, 0 ) - testPartialWrtMeanMomentOfInertia( i, 0 ) ), 1.0E1 );
    }

    // Check that the partials w.r.t. gravity field parameters are zero when the inertia tensor is explicitly defined
    phobos->setMassProperties( createRigidBodyProperties(
            constantRigidBodyPropertiesSettings( phobos->getBodyMass( ), Eigen::Vector3d::Zero( ), phobos->getBodyInertiaTensor( ) ),
            "Phobos",
            bodies ) );
    torquePartial =
            createAnalyticalTorquePartial( inertialTorqueModel, std::make_pair( "Phobos", phobos ), std::make_pair( "Phobos", phobos ) );
    torquePartial->update( testTime );
    BOOST_CHECK_EQUAL( torquePartial->wrtParameter( phobosGravitationalParameterParameter ).norm( ), 0.0 );
    BOOST_CHECK_EQUAL( torquePartial->wrtParameter( phobosCosineCoefficientsParameter ).norm( ), 0.0 );
    BOOST_CHECK_EQUAL( torquePartial->wrtParameter( phobosSineCoefficientsParameter ).norm( ), 0.0 );
}

//! Test partial of inertial torque w.r.t. angular velocity, for a body with a time-varying inertia tensor
//! (torque -omega x ( I omega ) - dI/dt omega), by comparing with numerical partials
BOOST_AUTO_TEST_CASE( testInertialTorquePartialsWithInertiaTensorDerivative )
{
    // Create body with mass-dependent rigid body properties, and a prescribed inertia tensor derivative
    SystemOfBodies bodies = SystemOfBodies( "SSB", "ECLIPJ2000" );
    bodies.createEmptyBody( "Vehicle" );
    std::shared_ptr< Body > vehicle = bodies.at( "Vehicle" );

    Eigen::Matrix3d inertiaTensor;
    inertiaTensor << 5.0E3, 1.0E2, -2.0E2, 1.0E2, 6.0E3, 3.0E2, -2.0E2, 3.0E2, 8.0E3;
    Eigen::Matrix3d inertiaTensorDerivative;
    inertiaTensorDerivative << -2.0E1, 3.0, 1.0, 3.0, -1.0E1, -4.0, 1.0, -4.0, -3.0E1;

    std::shared_ptr< MassDependentRigidBodyProperties > rigidBodyProperties =
            std::make_shared< MassDependentRigidBodyProperties >( 1.0E3, nullptr, [ = ]( const double ) { return inertiaTensor; } );
    rigidBodyProperties->setInertiaTensorDerivativeFunction( [ = ]( const double ) { return inertiaTensorDerivative; } );
    vehicle->setMassProperties( rigidBodyProperties );

    // Set rotational state of body
    Eigen::Vector7d rotationalState = Eigen::Vector7d::Zero( );
    rotationalState.segment( 0, 4 ) = tudat::linear_algebra::convertQuaternionToVectorFormat( Eigen::Quaterniond(
            Eigen::AngleAxisd( 0.4343, Eigen::Vector3d::UnitZ( ) ) * Eigen::AngleAxisd( 2.4354, Eigen::Vector3d::UnitX( ) ) ) );
    rotationalState( 4 ) = 1.0E-5;
    rotationalState( 5 ) = 2.0E-5;
    rotationalState( 6 ) = -3.2E-5;
    vehicle->setCurrentRotationalStateToLocalFrame( rotationalState );

    double testTime = 1000.0;
    rigidBodyProperties->update( testTime );

    // Create inertial torque and associated partial
    std::shared_ptr< InertialTorqueModel > inertialTorqueModel = createInertialTorqueModel( vehicle, "Vehicle" );
    BOOST_CHECK( inertialTorqueModel->isInertiaTensorDerivativeIncluded( ) );
    std::shared_ptr< TorquePartial > torquePartial = createAnalyticalTorquePartial(
            inertialTorqueModel, std::make_pair( "Vehicle", vehicle ), std::make_pair( "Vehicle", vehicle ) );

    // Calculate analytical partial w.r.t. angular velocity
    inertialTorqueModel->updateMembers( testTime );
    torquePartial->update( testTime );
    Eigen::MatrixXd partialWrtRotationalVelocity = Eigen::Matrix3d::Zero( );
    torquePartial->wrtRotationalVelocityOfAcceleratedBody( partialWrtRotationalVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );

    // Calculate numerical partial w.r.t. angular velocity
    Eigen::Vector3d rotationalVelocityPerturbation = Eigen::Vector3d::Constant( 1.0E-6 );
    std::function< void( Eigen::Vector7d ) > rotationalStateSetFunction =
            std::bind( &Body::setCurrentRotationalStateToLocalFrame, vehicle, std::placeholders::_1 );
    Eigen::Matrix3d testPartialWrtRotationalVelocity = calculateTorqueWrtRotationalStatePartials(
            rotationalStateSetFunction, inertialTorqueModel, rotationalState, rotationalVelocityPerturbation, 4, 3 );

    // Check that the inertia tensor derivative dominates the partial, so that the comparison tests its contribution
    BOOST_CHECK_GT( inertiaTensorDerivative.norm( ), 10.0 * ( partialWrtRotationalVelocity + inertiaTensorDerivative ).norm( ) );

    // Compare analytical and numerical partials
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtRotationalVelocity, partialWrtRotationalVelocity, 1.0E-8 );
}

class EffectiveTorqueModel : public basic_astrodynamics::TorqueModel
{
public:
    EffectiveTorqueModel( const std::function< Eigen::Matrix3d( ) > inertiaTensorFunction, const SingleBodyTorqueModelMap& torqueList ):
        inertiaTensorFunction_( inertiaTensorFunction ), torqueList_( torqueList )
    {}

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
    bodies.at( "Mars" )->setEphemeris(
            std::make_shared< ephemerides::ConstantEphemeris >( [ = ]( ) { return Eigen::Vector6d::Zero( ); } ) );
    bodies.at( "Mars" )->setGravityFieldModel(
            std::make_shared< gravitation::GravityFieldModel >( spice_interface::getBodyGravitationalParameter( "Mars" ) ) );
    bodies.createEmptyBody( "Phobos" );

    Eigen::Matrix3d phobosInertiaTensor = Eigen::Matrix3d::Zero( );
    phobosInertiaTensor( 0, 0 ) = 0.3615;
    phobosInertiaTensor( 1, 1 ) = 0.4265;
    phobosInertiaTensor( 2, 2 ) = 0.5024;

    phobosInertiaTensor *= ( 11.27E3 * 11.27E3 * 1.0659E16 );

    double phobosGravitationalParameter = 1.0659E16 * physical_constants::GRAVITATIONAL_CONSTANT;
    double phobosReferenceRadius = 11.27E3;

    Eigen::MatrixXd phobosCosineGravityFieldCoefficients = Eigen::MatrixXd::Zero( 6, 6 ),
                    phobosSineGravityFieldCoefficients = Eigen::MatrixXd::Zero( 6, 6 );
    double phobosScaledMeanMomentOfInertia;
    gravitation::getDegreeTwoSphericalHarmonicCoefficients( phobosInertiaTensor,
                                                            phobosGravitationalParameter,
                                                            phobosReferenceRadius,
                                                            true,
                                                            phobosCosineGravityFieldCoefficients,
                                                            phobosSineGravityFieldCoefficients,
                                                            phobosScaledMeanMomentOfInertia );

    const std::shared_ptr< gravitation::SphericalHarmonicsGravityField > phobosGravityField =
            std::make_shared< gravitation::SphericalHarmonicsGravityField >( phobosGravitationalParameter,
                                                                             phobosReferenceRadius,
                                                                             phobosCosineGravityFieldCoefficients,
                                                                             phobosSineGravityFieldCoefficients,
                                                                             "Phobos_Fixed" );
    bodies.at( "Phobos" )->setGravityFieldModel( phobosGravityField );
    bodies.at( "Phobos" )
            ->setMassProperties( std::make_shared< simulation_setup::FromGravityFieldRigidBodyProperties >(
                    phobosGravityField, phobosScaledMeanMomentOfInertia ) );

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
    bodies.at( "Phobos" )
            ->setRotationalEphemeris(
                    std::make_shared< TabulatedRotationalEphemeris< double, double > >( dummyInterpolator, "ECLIPJ2000", "Phobos_Fixed" ) );

    Eigen::Vector6d phobosKeplerElements = Eigen::Vector6d::Zero( );
    double phobosSemiMajorAxis = 9376.0E3;
    phobosKeplerElements( 0 ) = phobosSemiMajorAxis;

    bodies.at( "Phobos" )
            ->setEphemeris( std::make_shared< ephemerides::KeplerEphemeris >(
                    phobosKeplerElements, 0.0, spice_interface::getBodyGravitationalParameter( "Mars" ), "Mars", "ECLIPJ2000" ) );

    // Create empty bodies, phobos and mars.
    std::shared_ptr< Body > phobos = bodies.at( "Phobos" );
    std::shared_ptr< Body > mars = bodies.at( "Mars" );

    double testTime = 1000.0;
    phobos->setStateFromEphemeris( testTime );

    Eigen::Vector7d phobosRotationalState = Eigen::Vector7d::Zero( );
    phobosRotationalState.segment( 0, 4 ) = tudat::linear_algebra::convertQuaternionToVectorFormat( Eigen::Quaterniond(
            Eigen::AngleAxisd( 0.4343, Eigen::Vector3d::UnitZ( ) ) * Eigen::AngleAxisd( 2.4354, Eigen::Vector3d::UnitX( ) ) *
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
    std::shared_ptr< InertialTorqueModel > inertialTorqueModel = createInertialTorqueModel( bodies.at( "Phobos" ), "Phobos" );
    gravitationalTorque->updateMembers( 0.0 );
    inertialTorqueModel->updateMembers( 0.0 );

    SingleBodyTorqueModelMap torqueList;
    torqueList[ "Mars" ].push_back( gravitationalTorque );
    torqueList[ "Phobos" ].push_back( inertialTorqueModel );

    std::shared_ptr< EffectiveTorqueModel > effectiveTorqueModel =
            std::make_shared< EffectiveTorqueModel >( std::bind( &Body::getBodyInertiaTensor, bodies.at( "Phobos" ) ), torqueList );
    effectiveTorqueModel->updateMembers( 0.0 );

    // Create central gravity partial.
    std::shared_ptr< TorquePartial > torquePartial =
            createConstantTorqueRotationalDynamicsPartial( std::make_pair( "Phobos", bodies.at( "Phobos" ) ), torqueList );

    // Create gravitational parameter object.
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Phobos", gravitational_parameter ) );
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Phobos", mean_moment_of_inertia ) );

    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
            1, 0, 3, 3, "Phobos", spherical_harmonics_cosine_coefficient_block ) );
    parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
            1, 1, 3, 3, "Phobos", spherical_harmonics_sine_coefficient_block ) );

    std::shared_ptr< EstimatableParameterSet< double > > parameterSet = createParametersToEstimate( parameterNames, bodies );

    std::shared_ptr< EstimatableParameter< double > > phobosGravitationalParameterParameter =
            parameterSet->getEstimatedDoubleParameters( ).at( 0 );
    std::shared_ptr< EstimatableParameter< double > > meanMomentOfInertiaParameter = parameterSet->getEstimatedDoubleParameters( ).at( 1 );
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
    Eigen::MatrixXd partialWrtMarsOrientation = Eigen::MatrixXd::Zero( 3, 4 );
    torquePartial->wrtOrientationOfAcceleratingBody( partialWrtMarsOrientation.block( 0, 0, 3, 4 ) );
    Eigen::MatrixXd partialWrtMarsRotationalVelocity = Eigen::Matrix3d::Zero( );
    torquePartial->wrtRotationalVelocityOfAcceleratingBody( partialWrtMarsRotationalVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );

    Eigen::MatrixXd partialWrtPhobosState = Eigen::MatrixXd::Zero( 3, 6 );
    torquePartial->wrtNonRotationalStateOfAdditionalBody(
            partialWrtPhobosState.block( 0, 0, 3, 6 ), std::make_pair( "Phobos", "" ), propagators::translational_state );
    Eigen::MatrixXd partialWrtMarsState = Eigen::MatrixXd::Zero( 3, 6 );
    torquePartial->wrtNonRotationalStateOfAdditionalBody(
            partialWrtMarsState.block( 0, 0, 3, 6 ), std::make_pair( "Mars", "" ), propagators::translational_state );

    Eigen::Vector3d partialWrtPhobosGravitationalParameter = torquePartial->wrtParameter( phobosGravitationalParameterParameter );
    Eigen::Vector3d partialWrtMeanMomentOfInertia = torquePartial->wrtParameter( meanMomentOfInertiaParameter );

    Eigen::MatrixXd partialWrtPhobosCosineCoefficients = torquePartial->wrtParameter( phobosCosineCoefficientsParameter );
    Eigen::MatrixXd partialWrtPhobosSineCoefficients = torquePartial->wrtParameter( phobosSineCoefficientsParameter );

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
    rotationalVelocityPerturbation << 1.0E-6, 1.0E-6, 1.0E-6;

    // Create state access/modification functions for bodies.
    std::function< void( Eigen::Vector7d ) > phobosRotationalStateSetFunction =
            std::bind( &Body::setCurrentRotationalStateToLocalFrame, phobos, std::placeholders::_1 );
    std::function< void( Eigen::Vector7d ) > marsRotationalStateSetFunction =
            std::bind( &Body::setCurrentRotationalStateToLocalFrame, mars, std::placeholders::_1 );

    //    // Calculate numerical partials.
    std::vector< Eigen::Vector4d > appliedQuaternionPerturbation;
    Eigen::MatrixXd torqueDeviations = calculateTorqueDeviationDueToOrientationChange( phobosRotationalStateSetFunction,
                                                                                       effectiveTorqueModel,
                                                                                       phobos->getRotationalStateVector( ),
                                                                                       orientationPerturbation,
                                                                                       appliedQuaternionPerturbation );

    testPartialWrtPhobosRotationalVelocity = calculateTorqueWrtRotationalStatePartials( phobosRotationalStateSetFunction,
                                                                                        effectiveTorqueModel,
                                                                                        phobos->getRotationalStateVector( ),
                                                                                        rotationalVelocityPerturbation,
                                                                                        4,
                                                                                        3 );
    testPartialWrtMarsOrientation = calculateTorqueWrtRotationalStatePartials(
            marsRotationalStateSetFunction, effectiveTorqueModel, mars->getRotationalStateVector( ), orientationPerturbation, 0, 4 );
    testPartialWrtMarsRotationalVelocity = calculateTorqueWrtRotationalStatePartials(
            marsRotationalStateSetFunction, effectiveTorqueModel, mars->getRotationalStateVector( ), rotationalVelocityPerturbation, 4, 3 );

    std::function< void( Eigen::Vector6d ) > phobosStateSetFunction = std::bind( &Body::setState, phobos, std::placeholders::_1 );
    std::function< void( Eigen::Vector6d ) > marsStateSetFunction = std::bind( &Body::setState, mars, std::placeholders::_1 );

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

    std::function< void( ) > updateFunction =
            std::bind( &RigidBodyProperties::update, bodies.at( "Phobos" )->getMassProperties( ), testTime );
    // std::bind( &Body::setBodyInertiaTensorFromGravityFieldAndExistingMeanMoment, bodies.at( "Phobos" ), true );
    Eigen::Vector3d testPartialWrtPhobosGravitationalParameter =
            calculateTorqueWrtParameterPartials( phobosGravitationalParameterParameter, effectiveTorqueModel, 1.0E2, updateFunction );
    Eigen::MatrixXd testPartialWrtMeanMomentOfInertia =
            calculateTorqueWrtParameterPartials( meanMomentOfInertiaParameter, effectiveTorqueModel, 1.0E-4, updateFunction );
    Eigen::MatrixXd testPartialWrtPhobosCosineCoefficients = calculateTorqueWrtParameterPartials(
            phobosCosineCoefficientsParameter,
            effectiveTorqueModel,
            Eigen::VectorXd::Constant( phobosCosineCoefficientsParameter->getParameterSize( ), 1.0E-6 ),
            updateFunction );
    Eigen::MatrixXd testPartialWrtPhobosSineCoefficients =
            calculateTorqueWrtParameterPartials( phobosSineCoefficientsParameter,
                                                 effectiveTorqueModel,
                                                 Eigen::VectorXd::Constant( phobosSineCoefficientsParameter->getParameterSize( ), 1.0E-6 ),
                                                 updateFunction );

    //    // Compare numerical and analytical results.
    for( int index = 1; index < 4; index++ )
    {
        Eigen::Vector3d numericalChangeInTorque = torqueDeviations.block( 0, index - 1, 3, 1 );
        Eigen::Vector3d analyticalChangeInTorque =
                partialWrtPhobosOrientation.block( 0, 0, 3, 1 ) * appliedQuaternionPerturbation[ index ]( 0 ) +
                partialWrtPhobosOrientation.block( 0, index, 3, 1 ) * appliedQuaternionPerturbation[ index ]( index );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( numericalChangeInTorque, analyticalChangeInTorque, std::numeric_limits< double >::epsilon( ) );
    }

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPhobosRotationalVelocity, partialWrtPhobosRotationalVelocity, 1.0E-10 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            testPartialWrtMarsOrientation, partialWrtMarsOrientation, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            testPartialWrtMarsRotationalVelocity, partialWrtMarsRotationalVelocity, std::numeric_limits< double >::epsilon( ) );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            partialWrtPhobosState.block( 0, 0, 3, 3 ), testPartialWrtPhobosPosition, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            partialWrtMarsState.block( 0, 0, 3, 3 ), testPartialWrtMarsPosition, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            partialWrtPhobosState.block( 0, 3, 3, 3 ), testPartialWrtPhobosVelocity, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            partialWrtMarsState.block( 0, 3, 3, 3 ), testPartialWrtMarsVelocity, std::numeric_limits< double >::epsilon( ) );

    Eigen::Vector3d sclaedPartialWrtPhobosGravitationalParameter = phobosInertiaTensor * testPartialWrtPhobosGravitationalParameter;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( sclaedPartialWrtPhobosGravitationalParameter, partialWrtPhobosGravitationalParameter, 1.0E-6 );

    Eigen::Vector3d scaledPartialWrtMeanMomentOfInertia = phobosInertiaTensor * testPartialWrtMeanMomentOfInertia;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( scaledPartialWrtMeanMomentOfInertia, partialWrtMeanMomentOfInertia, 1.0E-6 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtPhobosCosineCoefficients.block( 0, 0, 3, 2 ),
                                       Eigen::MatrixXd::Zero( 3, 2 ),
                                       std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtPhobosCosineCoefficients.block( 0, 5, 3, 4 ),
                                       Eigen::MatrixXd::Zero( 3, 4 ),
                                       std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtPhobosSineCoefficients.block( 0, 0, 3, 1 ),
                                       Eigen::MatrixXd::Zero( 3, 1 ),
                                       std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtPhobosSineCoefficients.block( 0, 3, 3, 3 ),
                                       Eigen::MatrixXd::Zero( 3, 3 ),
                                       std::numeric_limits< double >::epsilon( ) );

    // Check derivative of z-component w.r.t. C20 separately: value is slightly non-zero due to rounding error
    Eigen::MatrixXd scaledPartialWrtPhobosCosineCoefficients = phobosInertiaTensor * testPartialWrtPhobosCosineCoefficients;
    BOOST_CHECK_SMALL( std::fabs( scaledPartialWrtPhobosCosineCoefficients( 2, 4 ) ), 1.0E5 );
    scaledPartialWrtPhobosCosineCoefficients( 2, 4 ) = 0.0;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( scaledPartialWrtPhobosCosineCoefficients, partialWrtPhobosCosineCoefficients, 1.0E-9 );

    Eigen::MatrixXd scaledPartialWrtPhobosSineCoefficients = phobosInertiaTensor * testPartialWrtPhobosSineCoefficients;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( scaledPartialWrtPhobosSineCoefficients, partialWrtPhobosSineCoefficients, 1.0E-9 );
}

//! Test that a spherical harmonic torque partial w.r.t. a tidal Love number is recomputed after a time reset, by comparing a partial that
//! was first evaluated at an earlier time with a newly created one.
BOOST_AUTO_TEST_CASE( testSphericalHarmonicTorquePartialWrtLoveNumberAfterTimeReset )
{
    spice_interface::loadStandardSpiceKernels( );

    // Create the Earth and a Moon with a degree-two gravity field and the tide raised by the Earth
    BodyListSettings bodySettings = getDefaultBodySettings( { "Earth", "Moon" }, "SSB", "ECLIPJ2000" );
    Eigen::MatrixXd cosineCoefficients = Eigen::MatrixXd::Zero( 3, 3 );
    cosineCoefficients( 0, 0 ) = 1.0;
    cosineCoefficients( 2, 0 ) = -9.09E-5;
    cosineCoefficients( 2, 2 ) = 3.47E-5;
    bodySettings.at( "Moon" )->gravityFieldSettings = std::make_shared< SphericalHarmonicsGravityFieldSettings >(
            4.9028E12, 1738.0E3, cosineCoefficients, Eigen::MatrixXd::Zero( 3, 3 ), "IAU_Moon" );
    bodySettings.at( "Moon" )->gravityFieldVariationSettings.push_back( std::make_shared< BasicSolidBodyGravityFieldVariationSettings >(
            std::vector< std::string >( { "Earth" } ),
            std::map< int, std::vector< std::complex< double > > >( { { 2, { 0.024, 0.024, 0.024 } } } ) ) );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    std::shared_ptr< Body > earth = bodies.at( "Earth" );
    std::shared_ptr< Body > moon = bodies.at( "Moon" );

    std::shared_ptr< TorqueModel > torqueModel = createSphericalHarmonicGravitationalTorqueModel(
            moon, earth, std::make_shared< SphericalHarmonicTorqueSettings >( 2, 2 ), "Moon", "Earth" );

    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( std::make_shared< FullDegreeTidalLoveNumberEstimatableParameterSettings >( "Moon", 2, "Earth", false ) );
    std::shared_ptr< EstimatableParameterSet< double > > parameterSet = createParametersToEstimate( parameterNames, bodies );
    std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > loveNumber = parameterSet->getEstimatedVectorParameters( ).at( 0 );

    // Create both partials before any evaluation, since creating a partial resets the Legendre cache shared with the torque model
    std::shared_ptr< TorquePartial > reusedPartial = createAnalyticalTorquePartial( torqueModel,
                                                                                    std::make_pair( "Moon", moon ),
                                                                                    std::make_pair( "Earth", earth ),
                                                                                    SingleBodyTorqueModelMap( ),
                                                                                    bodies,
                                                                                    parameterSet );
    std::shared_ptr< TorquePartial > newPartial = createAnalyticalTorquePartial( torqueModel,
                                                                                 std::make_pair( "Moon", moon ),
                                                                                 std::make_pair( "Earth", earth ),
                                                                                 SingleBodyTorqueModelMap( ),
                                                                                 bodies,
                                                                                 parameterSet );

    // Evaluate one partial at two times, one week apart, with a time reset in between
    const double testTime = 1.0E7 + 7.0 * 86400.0;
    Eigen::Vector3d partialOfReusedObject;
    for( const double time : { 1.0E7, testTime } )
    {
        earth->setStateFromEphemeris( time );
        moon->setStateFromEphemeris( time );
        moon->setCurrentRotationalStateToLocalFrameFromEphemeris( time );
        reusedPartial->resetCurrentTime( );
        reusedPartial->update( time );
        partialOfReusedObject = reusedPartial->wrtParameter( loveNumber );
    }

    // Compare with a partial that has only been evaluated at the last time, which cannot contain results from the earlier time
    newPartial->update( testTime );
    const Eigen::Vector3d partialOfNewObject = newPartial->wrtParameter( loveNumber );

    BOOST_CHECK_GT( partialOfNewObject.norm( ), 0.0 );
    BOOST_CHECK_SMALL( ( partialOfReusedObject - partialOfNewObject ).norm( ) / partialOfNewObject.norm( ), 1.0E-12 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
