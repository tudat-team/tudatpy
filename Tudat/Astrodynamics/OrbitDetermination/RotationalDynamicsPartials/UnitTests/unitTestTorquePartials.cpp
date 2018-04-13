/*    Copyright (c) 2010-2018, Delft University of Technology
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
#include "Tudat/Basics/testMacros.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>
#include <boost/lambda/lambda.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/sphericalStateConversions.h"
#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/Astrodynamics/Ephemerides/keplerEphemeris.h"
#include "Tudat/Astrodynamics/Relativity/metric.h"
#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/numericalAccelerationPartial.h"
#include "Tudat/SimulationSetup/EstimationSetup/createTorquePartials.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createBodies.h"
#include "Tudat/SimulationSetup/PropagationSetup/createTorqueModel.h"
#include "Tudat/SimulationSetup/EstimationSetup/createEstimatableParameters.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.h"

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
using namespace tudat::electro_magnetism;
using namespace tudat::basic_astrodynamics;

BOOST_AUTO_TEST_SUITE( test_torque_partials )

BOOST_AUTO_TEST_CASE( testSecondDegreeGravitationalTorquePartials )
{
    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    NamedBodyMap bodyMap;
    bodyMap[ "Mars" ] = boost::make_shared< Body >( );
    bodyMap[ "Mars" ]->setEphemeris( boost::make_shared< ephemerides::ConstantEphemeris >(
                                         boost::lambda::constant( Eigen::Vector6d::Zero( ) ) ) );
    bodyMap[ "Mars" ]->setGravityFieldModel(
                boost::make_shared< gravitation::GravityFieldModel >(
                    spice_interface::getBodyGravitationalParameter( "Mars" ) ) );
    bodyMap[ "Phobos" ] = boost::make_shared< Body >( );

    std::cout<<"Test A"<<std::endl;
    Eigen::Matrix3d phobosInertiaTensor = Eigen::Matrix3d::Zero( );
    phobosInertiaTensor( 0, 0 ) = 0.3615;
    phobosInertiaTensor( 1, 1 ) = 0.4265;
    phobosInertiaTensor( 2, 2 ) = 0.5024;


    phobosInertiaTensor *= (11.27E3 * 11.27E3 * 1.0659E16 );
    bodyMap[ "Phobos" ]->setBodyInertiaTensor( phobosInertiaTensor );

    double phobosGravitationalParameter = 1.0659E16 * physical_constants::GRAVITATIONAL_CONSTANT;
    double phobosReferenceRadius = 11.27E3;
    std::pair< Eigen::Matrix3d, Eigen::Matrix3d > phobosGravityFieldCoefficients =
            gravitation::getDegreeTwoSphericalHarmonicCoefficients(
                phobosInertiaTensor, phobosGravitationalParameter, phobosReferenceRadius, true );
    std::cout<<"Test B"<<std::endl;

    bodyMap[ "Phobos" ]->setGravityFieldModel(
                boost::make_shared< gravitation::SphericalHarmonicsGravityField >(
                    phobosGravitationalParameter, phobosReferenceRadius, phobosGravityFieldCoefficients.first,
                    phobosGravityFieldCoefficients.second, "Phobos_Fixed" ) );
    std::cout<<"Test C"<<std::endl;

    Eigen::Quaterniond noRotationQuaternion = Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) );
    Eigen::Matrix< double, 7, 1 > unitRotationState = Eigen::Matrix< double, 7, 1 >::Zero( );
    unitRotationState( 0 ) = noRotationQuaternion.w( );
    unitRotationState( 1 ) = noRotationQuaternion.x( );
    unitRotationState( 2 ) = noRotationQuaternion.y( );
    unitRotationState( 3 ) = noRotationQuaternion.z( );

    std::map< double, Eigen::Matrix< double, 7, 1 > > dummyRotationMap;
    dummyRotationMap[ -1.0E100 ] = unitRotationState;
    dummyRotationMap[ 1.0E100 ] = unitRotationState;

    boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< double, 7, 1 > > > dummyInterpolator =
            boost::make_shared< interpolators::LinearInterpolator< double, Eigen::Matrix< double, 7, 1 > > >( dummyRotationMap );
    bodyMap[ "Phobos" ]->setRotationalEphemeris( boost::make_shared< TabulatedRotationalEphemeris< double, double > >(
                                                     dummyInterpolator, "ECLIPJ2000", "Phobos_Fixed" ) );

    std::cout<<"Test D"<<std::endl;


    Eigen::Vector6d phobosKeplerElements = Eigen::Vector6d::Zero( );
    double phobosSemiMajorAxis = 9376.0E3;
    phobosKeplerElements( 0 ) = phobosSemiMajorAxis;

    bodyMap[ "Phobos" ]->setEphemeris( boost::make_shared< ephemerides::KeplerEphemeris >(
                                           phobosKeplerElements, 0.0, spice_interface::getBodyGravitationalParameter( "Mars" ),
                                           "Mars", "ECLIPJ2000" ) );

    std::cout<<"Test E"<<std::endl;

    // Create empty bodies, phobos and mars.
    boost::shared_ptr< Body > phobos = bodyMap.at( "Phobos" );
    boost::shared_ptr< Body > mars = bodyMap.at( "Mars" );


    // Create acceleration due to mars on phobos.
    boost::shared_ptr< SecondDegreeGravitationalTorqueModel > gravitationalTorque =
            createSecondDegreeGravitationalTorqueModel( bodyMap.at( "Phobos" ), bodyMap.at( "Mars" ), "Phobos", "Mars" );


    // Create central gravity partial.
    boost::shared_ptr< TorquePartial > torquePartial =
            createAnalyticalTorquePartial( gravitationalTorque, std::make_pair( "Phobos", phobos ),
                                           std::make_pair( "Mars", mars ) );

    std::cout<<"Test F"<<std::endl;

    // Create gravitational parameter object.
    //    boost::shared_ptr< EstimatableParameter< double > > marsGravitationalParameterParameter = boost::make_shared<
    //            GravitationalParameter >( marsGravityFieldModel, "Mars" );
    //    boost::shared_ptr< EstimatableParameter< double > > phobosGravitationalParameterParameter = boost::make_shared<
    //            GravitationalParameter >( phobosGravityFieldModel, "Phobos" );

    // Calculate analytical partials.
    torquePartial->update( 0.0 );
    Eigen::MatrixXd partialWrtPhobosOrientation = Eigen::MatrixXd::Zero( 3, 4 );
    torquePartial->wrtOrientationOfAcceleratedBody( partialWrtPhobosOrientation.block( 0, 0, 3, 4 ) );
    Eigen::MatrixXd partialWrtPhobosRotationalVelocity = Eigen::Matrix3d::Zero( );
    torquePartial->wrtRotationalVelocityOfAcceleratedBody( partialWrtPhobosRotationalVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );
    Eigen::MatrixXd partialWrtMarsOrientation = Eigen::MatrixXd::Zero( 3, 4  );
    torquePartial->wrtOrientationOfAcceleratingBody( partialWrtMarsOrientation.block( 0, 0, 4, 3 ) );
    Eigen::MatrixXd partialWrtMarsRotationalVelocity = Eigen::Matrix3d::Zero( );
    torquePartial->wrtRotationalVelocityOfAcceleratingBody( partialWrtMarsRotationalVelocity.block( 0, 0, 3, 3 ), 1, 0, 0 );
    //    Eigen::Vector3d partialWrtMarsGravitationalParameter = centralGravitationPartial->wrtParameter(
    //                marsGravitationalParameterParameter );
    //    Eigen::Vector3d partialWrtPhobosGravitationalParameter = centralGravitationPartial->wrtParameter(
    //                phobosGravitationalParameterParameter );
    std::cout<<"Test G"<<std::endl;

    // Declare numerical partials.
    Eigen::Matrix< double, 3, 4 > testPartialWrtPhobosOrientation = Eigen::Matrix< double, 3, 4 >::Zero( );
    Eigen::Matrix3d testPartialWrtPhobosRotationalVelocity = Eigen::Matrix3d::Zero( );
    Eigen::Matrix< double, 3, 4 > testPartialWrtMarsOrientation = Eigen::Matrix< double, 3, 4 >::Zero( );
    Eigen::Matrix3d testPartialWrtMarsRotationalVelocity = Eigen::Matrix3d::Zero( );

    // Declare perturbations in orientation for numerical partial/
    Eigen::Vector3d orientationPerturbation;
    orientationPerturbation << 1.0E-6, 1.0E-6, 1.0E-6, 1.0E-6;
    Eigen::Vector3d rotationalVelocityPerturbation;
    rotationalVelocityPerturbation <<  1.0E-6, 1.0E-6, 1.0E-6;

    // Create state access/modification functions for bodies.
    boost::function< void( Eigen::Vector7d ) > phobosRotationalStateSetFunction =
            boost::bind( &Body::setCurrentRotationalStateToLocalFrame, phobos, _1 );
    boost::function< void( Eigen::Vector7d ) > marsRotationalStateSetFunction =
            boost::bind( &Body::setCurrentRotationalStateToLocalFrame, mars, _1 );

    // Calculate numerical partials.
    testPartialWrtPhobosOrientation =  calculateTorqueWrtRotationalStatePartials(
                phobosRotationalStateSetFunction, gravitationalTorque, phobos->getRotationalStateVector( ), orientationPerturbation, 0, 4 );
    testPartialWrtPhobosRotationalVelocity =  calculateTorqueWrtRotationalStatePartials(
                phobosRotationalStateSetFunction, gravitationalTorque, phobos->getRotationalStateVector( ), rotationalVelocityPerturbation, 4, 3 );
    testPartialWrtMarsOrientation =  calculateTorqueWrtRotationalStatePartials(
                marsRotationalStateSetFunction, gravitationalTorque, mars->getRotationalStateVector( ), orientationPerturbation, 0, 4 );
    testPartialWrtMarsRotationalVelocity =  calculateTorqueWrtRotationalStatePartials(
                marsRotationalStateSetFunction, gravitationalTorque, mars->getRotationalStateVector( ), rotationalVelocityPerturbation, 4, 3 );
    //    Eigen::Vector3d testPartialWrtMarsGravitationalParameter = calculateAccelerationWrtParameterPartials(
    //                marsGravitationalParameterParameter, gravitationalTorque, 1.0E12 );
    //    Eigen::Vector3d testPartialWrtPhobosGravitationalParameter = calculateAccelerationWrtParameterPartials(
    //                phobosGravitationalParameterParameter, gravitationalTorque, 1.0E12 );

    // Compare numerical and analytical results.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPhobosOrientation,
                                       partialWrtPhobosOrientation, 1.0E-8 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtPhobosRotationalVelocity,
                                       partialWrtPhobosRotationalVelocity, std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtMarsOrientation,
                                       partialWrtMarsOrientation, 1.0E-8 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtMarsRotationalVelocity,
                                       partialWrtMarsRotationalVelocity, std::numeric_limits< double >::epsilon( ) );
    //    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( testPartialWrtMarsGravitationalParameter,
    //                                       partialWrtMarsGravitationalParameter, 1.0E-6 );
    //    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtPhobosGravitationalParameter,
    //                                       testPartialWrtPhobosGravitationalParameter, 1.0E-6 );
    //    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtPhobosGravitationalParameter,
    //                                       partialWrtMarsGravitationalParameter, std::numeric_limits< double >::epsilon(  ) );
}

BOOST_AUTO_TEST_SUITE_END( )


} // namespace unit_tests

} // namespace tudat




