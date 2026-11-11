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
#include "tudat/astro/gravitation/timeDependentSphericalHarmonicsGravityField.h"
#include "tudat/simulation/environment_setup/createGravityFieldVariations.h"

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

    std::shared_ptr< MassDependentRigidBodyProperties > rigidBodyProperties = std::make_shared< MassDependentRigidBodyProperties >(
            1.0E3, nullptr, [ = ]( const double ) { return inertiaTensor; }, [ = ]( const double ) { return inertiaTensorDerivative; } );
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

//! Function to retrieve analytical partials of the inertial torque, and of the tidal part of the inertia tensor and of its time
//! derivative, w.r.t. the entries of a parameter. For the gravitational parameter of the deformed body, the part of the inertia tensor
//! and its derivative that is proportional to this parameter is added, so that the full partials are returned.
template< typename ParameterType >
void addTidalInertialTorqueParameterPartials( const std::shared_ptr< EstimatableParameter< ParameterType > > parameter,
                                              const std::shared_ptr< InertialTorquePartial > torquePartial,
                                              const std::shared_ptr< Body > deformedBody,
                                              const std::string& deformedBodyName,
                                              std::vector< Eigen::Vector3d >& torquePartials,
                                              std::vector< Eigen::Matrix3d >& inertiaTensorPartials,
                                              std::vector< Eigen::Matrix3d >& inertiaTensorDerivativePartials,
                                              std::vector< bool >& isDeformedBodyGravitationalParameter )
{
    Eigen::MatrixXd torquePartialMatrix = torquePartial->wrtParameter( parameter );

    std::pair< InertialTorquePartial::TidalInertiaTensorPartialFunction, int > tidalPartialFunction =
            torquePartial->getTidalInertiaTensorPartialFunction( parameter );
    BOOST_REQUIRE_EQUAL( tidalPartialFunction.second, parameter->getParameterSize( ) );
    std::vector< Eigen::Matrix3d > currentInertiaTensorPartials, currentInertiaTensorDerivativePartials;
    tidalPartialFunction.first( currentInertiaTensorPartials, currentInertiaTensorDerivativePartials );

    const bool isCurrentParameterDeformedBodyGravitationalParameter = parameter->getParameterName( ).first == gravitational_parameter &&
            parameter->getParameterName( ).second.first == deformedBodyName;
    if( isCurrentParameterDeformedBodyGravitationalParameter )
    {
        const double gravitationalParameter = deformedBody->getGravityFieldModel( )->getGravitationalParameter( );
        currentInertiaTensorPartials.at( 0 ) += deformedBody->getBodyInertiaTensor( ) / gravitationalParameter;
        currentInertiaTensorDerivativePartials.at( 0 ) += deformedBody->getBodyInertiaTensorDerivative( ) / gravitationalParameter;
    }

    for( int i = 0; i < parameter->getParameterSize( ); i++ )
    {
        torquePartials.push_back( torquePartialMatrix.block( 0, i, 3, 1 ) );
        inertiaTensorPartials.push_back( currentInertiaTensorPartials.at( i ) );
        inertiaTensorDerivativePartials.push_back( currentInertiaTensorDerivativePartials.at( i ) );
        isDeformedBodyGravitationalParameter.push_back( isCurrentParameterDeformedBodyGravitationalParameter );
    }
}

//! Function to test partials of the inertial torque on Phobos w.r.t. its rotational state and the given parameters, for an inertia
//! tensor (and its derivative) derived from a gravity field with the given tidal variations (caused by Mars and/or a second body), by
//! comparing with numerical partials. The partials of the inertia tensor and of its derivative w.r.t. the quaternion and the parameters
//! are each compared with numerical partials, and the contributions of both to the torque partials are checked to be significant.
void testTidalInertialTorquePartials( const std::vector< std::shared_ptr< GravityFieldVariationSettings > >& gravityFieldVariationSettings,
                                      const std::vector< std::shared_ptr< EstimatableParameterSettings > >& parameterNames )
{
    // Create Mars, a second body causing tidal deformation (with a gravitational parameter and position such that its tidal effect is
    // significant), and Phobos with gravity-derived inertia tensor and tidal gravity field variations
    const double marsGravitationalParameter = 4.282837E13;
    BodyListSettings bodySettings( "SSB", "ECLIPJ2000" );
    bodySettings.addSettings( "Mars" );
    bodySettings.at( "Mars" )->ephemerisSettings = constantEphemerisSettings( Eigen::Vector6d::Zero( ) );
    bodySettings.at( "Mars" )->gravityFieldSettings = centralGravitySettings( marsGravitationalParameter );

    bodySettings.addSettings( "Perturber" );
    Eigen::Vector6d perturberState = Eigen::Vector6d::Zero( );
    perturberState.segment( 0, 3 ) << 2.0E7, -1.5E7, 1.0E7;
    bodySettings.at( "Perturber" )->ephemerisSettings = constantEphemerisSettings( perturberState );
    bodySettings.at( "Perturber" )->gravityFieldSettings = centralGravitySettings( 0.2 * marsGravitationalParameter );

    bodySettings.addSettings( "Phobos" );
    Eigen::Vector6d phobosKeplerElements = Eigen::Vector6d::Zero( );
    phobosKeplerElements( 0 ) = 9376.0E3;
    phobosKeplerElements( 1 ) = 0.015;
    phobosKeplerElements( 2 ) = 0.02;
    bodySettings.at( "Phobos" )->ephemerisSettings =
            keplerEphemerisSettings( phobosKeplerElements, 0.0, marsGravitationalParameter, "Mars", "ECLIPJ2000" );
    bodySettings.at( "Phobos" )->rotationModelSettings =
            constantRotationModelSettings( "ECLIPJ2000", "Phobos_Fixed", Eigen::Quaterniond::Identity( ) );
    Eigen::MatrixXd cosineCoefficients = Eigen::MatrixXd::Zero( 3, 3 );
    cosineCoefficients( 0, 0 ) = 1.0;
    cosineCoefficients( 2, 0 ) = -0.05;
    cosineCoefficients( 2, 2 ) = 0.02;
    bodySettings.at( "Phobos" )->gravityFieldSettings =
            sphericalHarmonicsGravitySettings( 7.11E5, 11.27E3, cosineCoefficients, Eigen::MatrixXd::Zero( 3, 3 ), "Phobos_Fixed" );
    bodySettings.at( "Phobos" )->gravityFieldVariationSettings = gravityFieldVariationSettings;
    bodySettings.at( "Phobos" )->rigidBodyPropertiesSettings = fromGravityFieldRigidBodyPropertiesSettings( 0.4 );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    std::shared_ptr< Body > phobos = bodies.at( "Phobos" );

    // Set states of bodies, and update gravity field (and inertia tensor derivative) of Phobos
    const double testTime = 1000.0;
    bodies.at( "Mars" )->setStateFromEphemeris( testTime );
    bodies.at( "Perturber" )->setStateFromEphemeris( testTime );
    phobos->setStateFromEphemeris( testTime );
    Eigen::Vector7d rotationalState = Eigen::Vector7d::Zero( );
    rotationalState.segment( 0, 4 ) = tudat::linear_algebra::convertQuaternionToVectorFormat( Eigen::Quaterniond(
            Eigen::AngleAxisd( 0.4343, Eigen::Vector3d::UnitZ( ) ) * Eigen::AngleAxisd( 2.4354, Eigen::Vector3d::UnitX( ) ) ) );
    rotationalState.segment( 4, 3 ) << 1.0E-4, -2.0E-4, 3.0E-4;
    phobos->setCurrentRotationalStateToLocalFrame( rotationalState );

    std::shared_ptr< gravitation::TimeDependentSphericalHarmonicsGravityField > phobosGravityField =
            std::dynamic_pointer_cast< gravitation::TimeDependentSphericalHarmonicsGravityField >( phobos->getGravityFieldModel( ) );
    std::function< void( ) > updateFunction = [ = ]( ) {
        phobosGravityField->update( testTime );
        phobos->getMassProperties( )->update( testTime );
    };
    updateFunction( );

    // Create inertial torque, parameters and torque partial
    std::shared_ptr< InertialTorqueModel > inertialTorqueModel = createInertialTorqueModel( phobos, "Phobos" );
    inertialTorqueModel->updateMembers( testTime );
    const Eigen::Matrix3d inertiaTensor = phobos->getBodyInertiaTensor( );
    const Eigen::Matrix3d inertiaTensorDerivative = phobos->getBodyInertiaTensorDerivative( );
    const Eigen::Vector3d angularVelocity = rotationalState.segment( 4, 3 );
    const Eigen::Matrix3d angularVelocityCrossProductMatrix = tudat::linear_algebra::getCrossProductMatrix( angularVelocity );

    std::shared_ptr< EstimatableParameterSet< double > > parameterSet = createParametersToEstimate( parameterNames, bodies );
    std::shared_ptr< InertialTorquePartial > torquePartial =
            std::dynamic_pointer_cast< InertialTorquePartial >( createAnalyticalTorquePartial( inertialTorqueModel,
                                                                                               std::make_pair( "Phobos", phobos ),
                                                                                               std::make_pair( "Phobos", phobos ),
                                                                                               SingleBodyTorqueModelMap( ),
                                                                                               bodies,
                                                                                               parameterSet ) );
    BOOST_REQUIRE( torquePartial != nullptr );

    // Calculate analytical partials w.r.t. rotational state, and partials of tidal part of inertia tensor w.r.t. quaternion
    torquePartial->update( testTime );
    Eigen::MatrixXd partialWrtOrientation = Eigen::MatrixXd::Zero( 3, 4 );
    torquePartial->wrtOrientationOfAcceleratedBody( partialWrtOrientation.block( 0, 0, 3, 4 ) );
    Eigen::MatrixXd partialWrtRotationalVelocity = Eigen::Matrix3d::Zero( );
    torquePartial->wrtRotationalVelocityOfAcceleratedBody( partialWrtRotationalVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );
    const std::vector< Eigen::Matrix3d > inertiaTensorPartialsWrtQuaternion =
            torquePartial->getCurrentTidalInertiaTensorPartialsWrtQuaternion( );
    BOOST_REQUIRE_EQUAL( static_cast< int >( inertiaTensorPartialsWrtQuaternion.size( ) ), 4 );

    // Calculate analytical partials w.r.t. parameters, per parameter entry (double parameters precede vector parameters)
    std::vector< Eigen::Vector3d > torquePartialsWrtParameters;
    std::vector< Eigen::Matrix3d > inertiaTensorPartialsWrtParameters, inertiaTensorDerivativePartialsWrtParameters;
    std::vector< bool > isDeformedBodyGravitationalParameter;
    for( std::shared_ptr< EstimatableParameter< double > > parameter : parameterSet->getEstimatedDoubleParameters( ) )
    {
        addTidalInertialTorqueParameterPartials( parameter,
                                                 torquePartial,
                                                 phobos,
                                                 "Phobos",
                                                 torquePartialsWrtParameters,
                                                 inertiaTensorPartialsWrtParameters,
                                                 inertiaTensorDerivativePartialsWrtParameters,
                                                 isDeformedBodyGravitationalParameter );
    }
    for( std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > parameter : parameterSet->getEstimatedVectorParameters( ) )
    {
        addTidalInertialTorqueParameterPartials( parameter,
                                                 torquePartial,
                                                 phobos,
                                                 "Phobos",
                                                 torquePartialsWrtParameters,
                                                 inertiaTensorPartialsWrtParameters,
                                                 inertiaTensorDerivativePartialsWrtParameters,
                                                 isDeformedBodyGravitationalParameter );
    }

    // Calculate numerical partials w.r.t. rotational state (updating the tidal gravity field variations and inertia tensor at each
    // evaluation)
    std::function< void( Eigen::Vector7d ) > rotationalStateSetFunction =
            std::bind( &Body::setCurrentRotationalStateToLocalFrame, phobos, std::placeholders::_1 );
    Eigen::Matrix3d testPartialWrtRotationalVelocity = calculateTorqueWrtRotationalStatePartials( rotationalStateSetFunction,
                                                                                                  inertialTorqueModel,
                                                                                                  rotationalState,
                                                                                                  Eigen::Vector3d::Constant( 1.0E-8 ),
                                                                                                  4,
                                                                                                  3,
                                                                                                  updateFunction,
                                                                                                  testTime );

    std::vector< Eigen::Vector4d > appliedQuaternionPerturbation;
    Eigen::MatrixXd torqueDeviations = calculateTorqueDeviationDueToOrientationChange( rotationalStateSetFunction,
                                                                                       inertialTorqueModel,
                                                                                       rotationalState,
                                                                                       Eigen::Vector4d::Constant( 1.0E-7 ),
                                                                                       appliedQuaternionPerturbation,
                                                                                       updateFunction,
                                                                                       testTime );

    // Calculate changes in inertia tensor due to the same quaternion perturbations
    std::vector< Eigen::Matrix3d > inertiaTensorDeviations;
    for( int index = 1; index < 4; index++ )
    {
        Eigen::Vector7d perturbedRotationalState = rotationalState;
        perturbedRotationalState.segment( 0, 4 ) += appliedQuaternionPerturbation[ index ];
        phobos->setCurrentRotationalStateToLocalFrame( perturbedRotationalState );
        updateFunction( );
        inertiaTensorDeviations.push_back( phobos->getBodyInertiaTensor( ) - inertiaTensor );
    }
    phobos->setCurrentRotationalStateToLocalFrame( rotationalState );
    updateFunction( );

    // Compare partial w.r.t. angular velocity with numerical partial, and check that the tidal contribution (through the
    // dependency of the inertia tensor derivative on angular velocity) is significant
    const Eigen::Matrix3d nonTidalPartialWrtRotationalVelocity = -angularVelocityCrossProductMatrix * inertiaTensor +
            tudat::linear_algebra::getCrossProductMatrix( inertiaTensor * angularVelocity ) - inertiaTensorDerivative;
    BOOST_CHECK_SMALL(
            ( partialWrtRotationalVelocity - testPartialWrtRotationalVelocity ).norm( ) / testPartialWrtRotationalVelocity.norm( ),
            1.0E-6 );
    BOOST_CHECK_GT( ( partialWrtRotationalVelocity - nonTidalPartialWrtRotationalVelocity ).norm( ),
                    1.0E-2 * inertiaTensorDerivative.norm( ) );

    // Compare partials of torque and inertia tensor w.r.t. orientation with numerical partials, and check that the contributions of the
    // inertia tensor partial and of the inertia tensor derivative partial (remainder) to the torque partial are both significant
    for( int index = 1; index < 4; index++ )
    {
        const Eigen::Vector4d& quaternionPerturbation = appliedQuaternionPerturbation[ index ];
        Eigen::Vector3d analyticalChangeInTorque = partialWrtOrientation.block( 0, 0, 3, 1 ) * quaternionPerturbation( 0 ) +
                partialWrtOrientation.block( 0, index, 3, 1 ) * quaternionPerturbation( index );
        Eigen::Matrix3d analyticalChangeInInertiaTensor = inertiaTensorPartialsWrtQuaternion.at( 0 ) * quaternionPerturbation( 0 ) +
                inertiaTensorPartialsWrtQuaternion.at( index ) * quaternionPerturbation( index );
        Eigen::Vector3d inertiaTensorContribution = -angularVelocityCrossProductMatrix * analyticalChangeInInertiaTensor * angularVelocity;

        BOOST_CHECK_SMALL(
                ( analyticalChangeInTorque - torqueDeviations.block( 0, index - 1, 3, 1 ) ).norm( ) / analyticalChangeInTorque.norm( ),
                1.0E-5 );
        BOOST_CHECK_SMALL( ( analyticalChangeInInertiaTensor - inertiaTensorDeviations.at( index - 1 ) ).norm( ) /
                                   analyticalChangeInInertiaTensor.norm( ),
                           1.0E-5 );
        BOOST_CHECK_GT( inertiaTensorContribution.norm( ), 1.0E-2 * analyticalChangeInTorque.norm( ) );
        BOOST_CHECK_GT( ( analyticalChangeInTorque - inertiaTensorContribution ).norm( ), 1.0E-2 * analyticalChangeInTorque.norm( ) );
    }

    // Compare partials of torque, inertia tensor and inertia tensor derivative w.r.t. each parameter entry with numerical partials, and
    // check that the contributions of the inertia tensor and of its derivative to the torque partial are both significant. The inertia
    // tensor derivative is due to the tidal variations only, which do not depend on the gravitational parameter of the deformed body (as
    // the tidal corrections are inversely proportional to it), so that its partial w.r.t. this parameter is zero.
    const Eigen::VectorXd nominalParameters = parameterSet->getFullParameterValues< double >( );
    BOOST_REQUIRE_EQUAL( static_cast< int >( torquePartialsWrtParameters.size( ) ), nominalParameters.rows( ) );
    for( int i = 0; i < nominalParameters.rows( ); i++ )
    {
        // Torque and inertia tensor are linear in the parameters, so that the perturbation only needs to be small w.r.t. the nominal value
        // to limit round-off errors
        const double perturbation = std::max( 1.0E-5 * std::fabs( nominalParameters( i ) ), 1.0E-3 );
        std::vector< Eigen::Vector3d > perturbedTorques;
        std::vector< Eigen::Matrix3d > perturbedInertiaTensors, perturbedInertiaTensorDerivatives;
        for( double perturbationSign : { 1.0, -1.0 } )
        {
            Eigen::VectorXd perturbedParameters = nominalParameters;
            perturbedParameters( i ) += perturbationSign * perturbation;
            parameterSet->resetParameterValues( perturbedParameters );
            updateFunction( );
            inertialTorqueModel->resetCurrentTime( );
            perturbedTorques.push_back( basic_astrodynamics::updateAndGetTorque( inertialTorqueModel, testTime ) );
            perturbedInertiaTensors.push_back( phobos->getBodyInertiaTensor( ) );
            perturbedInertiaTensorDerivatives.push_back( phobos->getBodyInertiaTensorDerivative( ) );
        }
        parameterSet->resetParameterValues( nominalParameters );
        updateFunction( );

        Eigen::Vector3d testTorquePartial = ( perturbedTorques.at( 0 ) - perturbedTorques.at( 1 ) ) / ( 2.0 * perturbation );
        Eigen::Matrix3d testInertiaTensorPartial =
                ( perturbedInertiaTensors.at( 0 ) - perturbedInertiaTensors.at( 1 ) ) / ( 2.0 * perturbation );
        Eigen::Matrix3d testInertiaTensorDerivativePartial =
                ( perturbedInertiaTensorDerivatives.at( 0 ) - perturbedInertiaTensorDerivatives.at( 1 ) ) / ( 2.0 * perturbation );

        BOOST_CHECK_SMALL( ( torquePartialsWrtParameters.at( i ) - testTorquePartial ).norm( ) / testTorquePartial.norm( ), 1.0E-6 );
        BOOST_CHECK_SMALL(
                ( inertiaTensorPartialsWrtParameters.at( i ) - testInertiaTensorPartial ).norm( ) / testInertiaTensorPartial.norm( ),
                1.0E-6 );
        BOOST_CHECK_GT( ( angularVelocityCrossProductMatrix * inertiaTensorPartialsWrtParameters.at( i ) * angularVelocity ).norm( ),
                        1.0E-2 * testTorquePartial.norm( ) );
        if( isDeformedBodyGravitationalParameter.at( i ) )
        {
            const double inertiaTensorDerivativePartialScale = inertiaTensorDerivative.norm( ) / nominalParameters( i );
            BOOST_CHECK_SMALL( inertiaTensorDerivativePartialsWrtParameters.at( i ).norm( ),
                               1.0E-10 * inertiaTensorDerivativePartialScale );
            BOOST_CHECK_SMALL( testInertiaTensorDerivativePartial.norm( ), 1.0E-6 * inertiaTensorDerivativePartialScale );
        }
        else
        {
            BOOST_CHECK_SMALL( ( inertiaTensorDerivativePartialsWrtParameters.at( i ) - testInertiaTensorDerivativePartial ).norm( ) /
                                       testInertiaTensorDerivativePartial.norm( ),
                               1.0E-6 );
            BOOST_CHECK_GT( ( inertiaTensorDerivativePartialsWrtParameters.at( i ) * angularVelocity ).norm( ),
                            1.0E-2 * testTorquePartial.norm( ) );
        }
    }
}

//! Test partials of inertial torque w.r.t. rotational state, tidal Love numbers and gravitational parameters, for a body with an
//! inertia tensor (and its derivative) varying due to tidal gravity field variations, by comparing with numerical partials
BOOST_AUTO_TEST_CASE( testInertialTorquePartialsWithTidalInertiaTensorDerivative )
{
    // Degree two tide raised by Mars and a second body, with equal Love numbers for all orders
    std::map< int, std::vector< std::complex< double > > > loveNumbers;
    loveNumbers[ 2 ] = std::vector< std::complex< double > >( 3, std::complex< double >( 0.1, 0.0 ) );
    std::shared_ptr< GravityFieldVariationSettings > basicTideSettings = std::make_shared< BasicSolidBodyGravityFieldVariationSettings >(
            std::vector< std::string >( { "Mars", "Perturber" } ), loveNumbers );

    // Full-degree (real) Love number, and gravitational parameters of the deformed body and of both bodies causing deformation
    testTidalInertialTorquePartials( { basicTideSettings },
                                     { std::make_shared< FullDegreeTidalLoveNumberEstimatableParameterSettings >(
                                               "Phobos", 2, std::vector< std::string >( { "Mars", "Perturber" } ), false ),
                                       std::make_shared< EstimatableParameterSettings >( "Mars", gravitational_parameter ),
                                       std::make_shared< EstimatableParameterSettings >( "Perturber", gravitational_parameter ),
                                       std::make_shared< EstimatableParameterSettings >( "Phobos", gravitational_parameter ) } );

    // Complex Love numbers for orders one and two (real and imaginary part per order), for a tide with order-dependent complex Love numbers
    std::map< int, std::vector< std::complex< double > > > complexLoveNumbers;
    complexLoveNumbers[ 2 ] = { std::complex< double >( 0.1, 0.0 ),
                                std::complex< double >( 0.12, 0.01 ),
                                std::complex< double >( 0.09, -0.02 ) };
    testTidalInertialTorquePartials(
            { std::make_shared< BasicSolidBodyGravityFieldVariationSettings >( std::vector< std::string >( { "Mars", "Perturber" } ),
                                                                               complexLoveNumbers ) },
            { std::make_shared< SingleDegreeVariableTidalLoveNumberEstimatableParameterSettings >(
                    "Phobos", 2, std::vector< int >( { 1, 2 } ), std::vector< std::string >( { "Mars", "Perturber" } ), true ) } );

    // Mode-coupled Love numbers due to Mars, including degree three forcing producing a degree two response (Love numbers chosen such
    // that this contribution is comparable to that of the degree two forcing), combined with the degree two tide. The gravitational
    // parameter partials include the contributions of both tides.
    std::map< std::pair< int, int >, std::map< std::pair< int, int >, double > > modeCoupledLoveNumbers;
    modeCoupledLoveNumbers[ std::make_pair( 2, 0 ) ][ std::make_pair( 2, 0 ) ] = 0.05;
    modeCoupledLoveNumbers[ std::make_pair( 3, 0 ) ][ std::make_pair( 2, 0 ) ] = 500.0;
    modeCoupledLoveNumbers[ std::make_pair( 3, 1 ) ][ std::make_pair( 2, 1 ) ] = 500.0;
    std::map< std::pair< int, int >, std::vector< std::pair< int, int > > > modeCoupledLoveNumberIndices;
    for( auto forcingIterator : modeCoupledLoveNumbers )
    {
        for( auto responseIterator : forcingIterator.second )
        {
            modeCoupledLoveNumberIndices[ forcingIterator.first ].push_back( responseIterator.first );
        }
    }
    testTidalInertialTorquePartials( { basicTideSettings,
                                       std::make_shared< ModeCoupledSolidBodyGravityFieldVariationSettings >(
                                               std::vector< std::string >( { "Mars" } ), modeCoupledLoveNumbers ) },
                                     { std::make_shared< ModeCoupledTidalLoveNumberEstimatableParameterSettings >(
                                               "Phobos", modeCoupledLoveNumberIndices, std::vector< std::string >( { "Mars" } ), false ),
                                       std::make_shared< EstimatableParameterSettings >( "Mars", gravitational_parameter ),
                                       std::make_shared< EstimatableParameterSettings >( "Phobos", gravitational_parameter ) } );
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

//! Test partial of inertial torque w.r.t. gravitational parameter, for a body with an inertia tensor (and its time derivative) derived
//! from a time-varying gravity field
BOOST_AUTO_TEST_CASE( testInertialTorqueGravitationalParameterPartialWithInertiaTensorDerivative )
{
    // Create body with gravity-derived rigid body properties, with linearly varying C20 and C22
    Eigen::MatrixXd cosineCoefficients = Eigen::MatrixXd::Zero( 3, 3 );
    cosineCoefficients( 0, 0 ) = 1.0;
    cosineCoefficients( 2, 0 ) = -1.0E-2;
    cosineCoefficients( 2, 2 ) = 2.0E-3;
    BodyListSettings bodySettings;
    bodySettings.addSettings( "Body" );
    bodySettings.at( "Body" )->gravityFieldSettings = std::make_shared< SphericalHarmonicsGravityFieldSettings >(
            1.0E6, 1.0E3, cosineCoefficients, Eigen::MatrixXd::Zero( 3, 3 ), "BodyFixed" );
    bodySettings.at( "Body" )->rigidBodyPropertiesSettings = fromGravityFieldRigidBodyPropertiesSettings( 0.4 );
    bodySettings.at( "Body" )->gravityFieldVariationSettings.push_back( polynomialGravityFieldVariationsSettings(
            { { 1, ( Eigen::MatrixXd( 1, 3 ) << 1.0E-4, 0.0, 2.0E-4 ).finished( ) } }, std::map< int, Eigen::MatrixXd >( ), 0.0 ) );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    std::shared_ptr< Body > body = bodies.at( "Body" );

    Eigen::Vector7d rotationalState = Eigen::Vector7d::Zero( );
    rotationalState( 0 ) = 1.0;
    rotationalState.segment( 4, 3 ) << 1.0E-4, -2.0E-4, 3.0E-4;
    body->setCurrentRotationalStateToLocalFrame( rotationalState );
    const double testTime = 1.0E3;
    body->getMassProperties( )->update( testTime );

    // Check that the inertia tensor derivative dominates the inertial torque, so that the comparison tests its contribution
    const Eigen::Vector3d angularVelocity = rotationalState.segment( 4, 3 );
    BOOST_CHECK_GT( ( body->getBodyInertiaTensorDerivative( ) * angularVelocity ).norm( ),
                    10.0 * angularVelocity.cross( body->getBodyInertiaTensor( ) * angularVelocity ).norm( ) );

    // Compare analytical and numerical partials w.r.t. gravitational parameter
    std::shared_ptr< InertialTorqueModel > inertialTorqueModel = createInertialTorqueModel( body, "Body" );
    std::shared_ptr< TorquePartial > torquePartial =
            createAnalyticalTorquePartial( inertialTorqueModel, std::make_pair( "Body", body ), std::make_pair( "Body", body ) );
    std::shared_ptr< EstimatableParameter< double > > gravitationalParameter =
            std::make_shared< GravitationalParameter >( body->getGravityFieldModel( ), "Body" );
    inertialTorqueModel->updateMembers( testTime );
    torquePartial->update( testTime );
    const Eigen::Vector3d partialWrtGravitationalParameter = torquePartial->wrtParameter( gravitationalParameter );
    const Eigen::Vector3d testPartialWrtGravitationalParameter =
            calculateTorqueWrtParameterPartials( gravitationalParameter,
                                                 inertialTorqueModel,
                                                 1.0E3,
                                                 std::bind( &RigidBodyProperties::update, body->getMassProperties( ), testTime ),
                                                 testTime );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtGravitationalParameter, partialWrtGravitationalParameter, 1.0E-8 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
