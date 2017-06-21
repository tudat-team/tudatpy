#define BOOST_TEST_MAIN

#include <limits>
#include <string>
#include <thread>


#include <boost/make_shared.hpp>
#include <boost/format.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Geometry>

#include "Tudat/Basics/testMacros.h"

#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"

#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/Mathematics/NumericalIntegrators/rungeKuttaCoefficients.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/InputOutput/basicInputOutput.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/torqueModelTypes.h"
#include "Tudat/External/SpiceInterface/spiceRotationalEphemeris.h"
#include "Tudat/Astrodynamics/Ephemerides/constantEphemeris.h"
#include "Tudat/Astrodynamics/Ephemerides/keplerEphemeris.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/Astrodynamics/Gravitation/centralGravityModel.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createBodies.h"
#include "Tudat/SimulationSetup/PropagationSetup/createTorqueModel.h"
#include "Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationSettings.h"
#include "Tudat/SimulationSetup/PropagationSetup/dynamicsSimulator.h"

namespace tudat
{
namespace unit_tests
{

//Using declarations.
using namespace tudat::ephemerides;
using namespace tudat::interpolators;
using namespace tudat::numerical_integrators;
using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;
using namespace tudat::basic_astrodynamics;
using namespace tudat::orbital_element_conversions;
using namespace tudat::propagators;

NamedBodyMap getTestBodyMap( const double phobosSemiMajorAxis,
                             const bool useSymmetricEquator = 0 )
{
    NamedBodyMap bodyMap;
    bodyMap[ "Mars" ] = boost::make_shared< Body >( );
    bodyMap[ "Mars" ]->setEphemeris( boost::make_shared< ephemerides::ConstantEphemeris >(
                                         boost::lambda::constant( Eigen::Vector6d::Zero( ) ) ) );
    //    bodyMap[ "Mars" ]->setRotationalEphemeris( boost::make_shared< ephemerides::ConstantRotationalEphemeris >(
    //                                                   boost::lambda::constant( Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ) ),
    //                                                   "ECLIPJ2000", "Mars_Fixed" ) );
    //    boost::dynamic_pointer_cast< simulation_setup::CelestialBody >( bodyMap[ "Mars" ] )->setGravityFieldModel(
    //                boost::make_shared< gravitation::GravityFieldModel >( spice_interface::getBodyGravitationalParameter( "Mars" ) ) );
    bodyMap[ "Phobos" ] = boost::make_shared< Body >( );
    Eigen::Vector6d phobosInitialStateInKeplerianElements;
    phobosInitialStateInKeplerianElements[ semiMajorAxisIndex ] = phobosSemiMajorAxis;
    phobosInitialStateInKeplerianElements[ eccentricityIndex ] = 0.015;

    bodyMap[ "Phobos" ]->setEphemeris( boost::make_shared< ephemerides::KeplerEphemeris >(
                                           phobosInitialStateInKeplerianElements, 1.0E7, getBodyGravitationalParameter( "Mars" ) ) );
    Eigen::Matrix3d phobosInertiaTensor = Eigen::Matrix3d::Zero( );
    phobosInertiaTensor( 0, 0 ) = 0.3615;
    phobosInertiaTensor( 1, 1 ) = 0.4265;
    phobosInertiaTensor( 2, 2 ) = 0.5024;

    if( useSymmetricEquator )
    {
        phobosInertiaTensor( 0, 0 ) = phobosInertiaTensor( 1, 1 );
    }

    phobosInertiaTensor *= (11.27E3 * 11.27E3 * 1.0659E16 );
    bodyMap[ "Phobos" ]->setBodyInertiaTensor( phobosInertiaTensor );

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


    return bodyMap;
}

Eigen::Vector3d calculateEulerAngleRatesFromEulerAngles(
        const Eigen::Vector3d& eulerAngles,
        const double elapsedTime,
        const double eulerFrequency,
        const double rotationFrequency,
        const double perturbationAmplitude,
        const double perturbationPhase )
{
    Eigen::Vector3d eulerAngleRates;
    eulerAngleRates( 0 ) = perturbationAmplitude * std::sin(
                eulerAngles( 2 ) + eulerFrequency * elapsedTime + perturbationPhase ) / std::sin( eulerAngles( 1 ) );
    eulerAngleRates( 1 ) = perturbationAmplitude * std::cos(
                eulerAngles( 2 ) + eulerFrequency * elapsedTime + perturbationPhase );
    eulerAngleRates( 2 ) = rotationFrequency - perturbationAmplitude * std::sin(
                eulerAngles( 2 ) + eulerFrequency * elapsedTime + perturbationPhase ) / std::tan( eulerAngles( 1 ) );
    return eulerAngleRates;
}

Eigen::Vector3d calculateEulerAngleRatesByNumericalDifference(
        const boost::shared_ptr< ephemerides::RotationalEphemeris > rotationalModel,\
        const double evaluationTime,
        const double timeStep )
{
    Eigen::Matrix3d upPerturbedRotationMatrix = rotationalModel->getRotationToBaseFrame( evaluationTime + timeStep ).toRotationMatrix( );
    Eigen::Vector3d upPerturbedEulerAngles = upPerturbedRotationMatrix.eulerAngles( 2, 0, 2 );

    Eigen::Matrix3d downPerturbedRotationMatrix = rotationalModel->getRotationToBaseFrame( evaluationTime - timeStep ).toRotationMatrix( );
    Eigen::Vector3d downPerturbedEulerAngles = downPerturbedRotationMatrix.eulerAngles( 2, 0, 2 );

    Eigen::Vector3d angleRates = ( upPerturbedEulerAngles - downPerturbedEulerAngles ) / ( 2.0 * timeStep );

    // Check domain flip
    if( std::fabs( upPerturbedEulerAngles( 2 ) - downPerturbedEulerAngles( 2 ) ) > 6.0 )
    {
        if( upPerturbedEulerAngles( 2 ) > downPerturbedEulerAngles( 2 ) )
        {
            angleRates( 2 ) = ( upPerturbedEulerAngles( 2 )  - 2.0 * mathematical_constants::PI - downPerturbedEulerAngles( 2 ) ) / ( 2.0 * timeStep );
        }
        else
        {
            angleRates( 2 ) = ( upPerturbedEulerAngles( 2 ) - ( downPerturbedEulerAngles( 2 ) - 2.0 * mathematical_constants::PI ) ) / ( 2.0 * timeStep );

        }
    }
    return angleRates;
}

BOOST_AUTO_TEST_SUITE( test_rotational_dynamics_propagation )

BOOST_AUTO_TEST_CASE( testSimpleRotationalDynamicsPropagation )
{
    //Load spice kernels.
    std::string kernelsPath = input_output::getSpiceKernelPath( );
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de-403-masses.tpc");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "naif0009.tls");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "pck00009.tpc");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de421.bsp");

    for( unsigned axisCase = 0; axisCase < 3; axisCase++ )
    {
        std::cout<<"Axis case: "<<axisCase<<std::endl;
        double phobosSemiMajorAxis = 9376.0E3;

        NamedBodyMap bodyMap = getTestBodyMap( 9376.0E3 );

        double initialEphemerisTime = 1.0E7;
        double finalEphemerisTime = initialEphemerisTime + 10.0 * 86400.0;

        // Set accelerations between bodies that are to be taken into account.
        SelectedTorqueMap torqueMap;
        std::vector< std::string > bodiesToIntegrate;
        bodiesToIntegrate.push_back( "Phobos" );

        double meanMotion = std::sqrt( getBodyGravitationalParameter( "Mars" ) /
                                       std::pow( phobosSemiMajorAxis, 3.0 ) );

        Eigen::Quaterniond initialRotation = Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) );
                reference_frames::getRotatingPlanetocentricToLocalVerticalFrameTransformationQuaternion(
                    0.2, 0.0 );
        Eigen::Matrix3d initialRotationMatrixToBaseFrame = initialRotation.toRotationMatrix( );
        Eigen::Matrix3d initialRotationMatrixToTargetFrame = initialRotationMatrixToBaseFrame.transpose( );

        Eigen::VectorXd systemInitialState = Eigen::VectorXd::Zero( 7 );
        systemInitialState.segment( 0, 4 ) = linear_algebra::convertQuaternionToVectorFormat(
                    initialRotation );
        systemInitialState( 4 + axisCase ) = meanMotion;

        basic_astrodynamics::TorqueModelMap torqueModelMap = createTorqueModelsMap(
                    bodyMap, torqueMap );

        boost::shared_ptr< numerical_integrators::IntegratorSettings< > > integratorSettings =
                boost::make_shared< RungeKuttaVariableStepSizeSettings< > >
                ( rungeKuttaVariableStepSize,
                  initialEphemerisTime, 30.0,
                  RungeKuttaCoefficients::rungeKuttaFehlberg78,
                  2.0, 30.0, 1.0E-13, 1.0E-13 );

        boost::shared_ptr< RotationalStatePropagatorSettings< double > > propagatorSettings =
                boost::make_shared< RotationalStatePropagatorSettings< double > >
                ( torqueModelMap, bodiesToIntegrate, systemInitialState, boost::make_shared< PropagationTimeTerminationSettings >(
                      finalEphemerisTime ) );

        SingleArcDynamicsSimulator< double > dynamicsSimulator(
                    bodyMap, integratorSettings, propagatorSettings, true, false, true );

        boost::shared_ptr< RotationalEphemeris > phobosEphemeris = bodyMap.at( "Phobos" )->getRotationalEphemeris( );

        double startTime = initialEphemerisTime;
        double endTime = finalEphemerisTime - 3600.0;
        double currentTime = startTime;

        boost::shared_ptr< RotationalEphemeris > phobosRotationalEphemeris = bodyMap[ "Phobos" ]->getRotationalEphemeris( );

        Eigen::Vector3d currentRotationalVelocityInTargetFrame, currentRotationalVelocityInBaseFrame;
        Eigen::Vector3d expectedRotationalVelocityVectorInBaseFrame,
                expectedRotationalVelocityVectorInTargetFrame = Eigen::Vector3d::Zero( );
        expectedRotationalVelocityVectorInTargetFrame( axisCase ) = meanMotion;

        Eigen::Matrix3d currentRotationMatrixToTargetFrame, currentRotationMatrixToBaseFrame;
        Eigen::Matrix3d expectedRotationMatrixToTargetFrame, expectedRotationMatrixToBaseFrame;

        Eigen::Matrix3d currentRotationMatrixDerivativeToTargetFrame, currentRotationMatrixDerivativeToBaseFrame, currentIndirectRotationMatrixDerivative;
        Eigen::Matrix3d expectedRotationMatrixDerivativeToTargetFrame, expectedRotationMatrixDerivativeToBaseFrame;

        Eigen::Matrix3d expectedRotationToTargetFrameFromInitialRotation, expectedRotationToBaseFrameFromInitialRotation;

        double timeStep = 600.0;//( endTime - startTime ) / 20.0;
        while ( currentTime < endTime )
        {
            double currentAngle = meanMotion * ( currentTime - initialEphemerisTime );

            Eigen::Vector3d rotationAxis;
            if( axisCase == 0 )
            {
                rotationAxis = Eigen::Vector3d::UnitX( );
            }
            else if( axisCase == 1 )
            {
                rotationAxis = Eigen::Vector3d::UnitY( );
            }
            else if( axisCase == 2 )
            {
                rotationAxis = Eigen::Vector3d::UnitZ( );
            }

            Eigen::Matrix3d baseRotationToTargetFrame =
                    Eigen::AngleAxisd( -1.0 * currentAngle, rotationAxis ).toRotationMatrix( );

            currentRotationMatrixToTargetFrame = phobosRotationalEphemeris->getRotationToTargetFrame( currentTime );
            expectedRotationToTargetFrameFromInitialRotation =
                    baseRotationToTargetFrame * initialRotationMatrixToTargetFrame;
            for( unsigned int i = 0; i < 3; i++ )
            {
                for( unsigned int j = 0; j < 3; j++ )
                {
                    BOOST_CHECK_SMALL(
                                std::fabs( currentRotationMatrixToTargetFrame( i, j ) -
                                           expectedRotationToTargetFrameFromInitialRotation( i, j ) ), 1.0E-10 );
                }
            }

            currentRotationMatrixToBaseFrame = phobosRotationalEphemeris->getRotationToBaseFrame( currentTime );
            expectedRotationToBaseFrameFromInitialRotation =
                    initialRotationMatrixToBaseFrame * baseRotationToTargetFrame.transpose( );
            for( unsigned int i = 0; i < 3; i++ )
            {
                for( unsigned int j = 0; j < 3; j++ )
                {
                    BOOST_CHECK_SMALL(
                                std::fabs( currentRotationMatrixToBaseFrame( i, j ) -
                                           expectedRotationToBaseFrameFromInitialRotation( i, j ) ), 1.0E-10 );
                }
            }


            currentRotationMatrixDerivativeToTargetFrame = phobosRotationalEphemeris->getDerivativeOfRotationToTargetFrame( currentTime );
            Eigen::Matrix3d premultiplierMatrix;

            if( axisCase == 0 )
            {
                premultiplierMatrix = reference_frames::X_AXIS_ROTATION_MATRIX_DERIVATIVE_PREMULTIPLIER;
            }
            else if( axisCase == 1 )
            {
                premultiplierMatrix = reference_frames::Y_AXIS_ROTATION_MATRIX_DERIVATIVE_PREMULTIPLIER ;
            }
            else if( axisCase == 2 )
            {
                premultiplierMatrix = reference_frames::Z_AXIS_ROTATION_MATRIX_DERIVATIVE_PREMULTIPLIER ;
            }


            expectedRotationMatrixDerivativeToTargetFrame =
                     meanMotion * premultiplierMatrix * baseRotationToTargetFrame *
                    initialRotationMatrixToTargetFrame;

            for( unsigned int i = 0; i < 3; i++ )
            {
                for( unsigned int j = 0; j < 3; j++ )
                {
                    BOOST_CHECK_SMALL(
                                std::fabs( currentRotationMatrixDerivativeToTargetFrame( i, j ) -
                                           expectedRotationMatrixDerivativeToTargetFrame( i, j ) ), meanMotion * 1.0E-10 );
                }
            }

            currentRotationMatrixDerivativeToBaseFrame = phobosRotationalEphemeris->getDerivativeOfRotationToBaseFrame( currentTime );
            expectedRotationMatrixDerivativeToBaseFrame = expectedRotationMatrixDerivativeToTargetFrame.transpose( );

            for( unsigned int i = 0; i < 3; i++ )
            {
                for( unsigned int j = 0; j < 3; j++ )
                {
                    BOOST_CHECK_SMALL(
                                std::fabs( currentRotationMatrixDerivativeToBaseFrame( i, j ) -
                                           expectedRotationMatrixDerivativeToBaseFrame( i, j ) ), meanMotion * 1.0E-10 );
                }
            }

            currentRotationalVelocityInTargetFrame = phobosRotationalEphemeris->getRotationalVelocityVectorInTargetFrame( currentTime );
            BOOST_CHECK_SMALL( std::fabs( currentRotationalVelocityInTargetFrame( 0 ) - expectedRotationalVelocityVectorInTargetFrame( 0 ) ),
                               meanMotion * 1.0E-15 );
            BOOST_CHECK_SMALL( std::fabs( currentRotationalVelocityInTargetFrame( 1 ) - expectedRotationalVelocityVectorInTargetFrame( 1 ) ),
                               meanMotion * 1.0E-15 );
            BOOST_CHECK_SMALL( std::fabs( currentRotationalVelocityInTargetFrame( 2 ) - expectedRotationalVelocityVectorInTargetFrame( 2 ) ),
                               meanMotion * 1.0E-15 );


            currentRotationalVelocityInBaseFrame = phobosRotationalEphemeris->getRotationalVelocityVectorInBaseFrame( currentTime );
            expectedRotationalVelocityVectorInBaseFrame =
                    currentRotationMatrixToBaseFrame * expectedRotationalVelocityVectorInTargetFrame;
            BOOST_CHECK_SMALL( std::fabs( currentRotationalVelocityInBaseFrame( 0 ) - expectedRotationalVelocityVectorInBaseFrame( 0 ) ),
                               meanMotion * 1.0E-15 );
            BOOST_CHECK_SMALL( std::fabs( currentRotationalVelocityInBaseFrame( 1 ) - expectedRotationalVelocityVectorInBaseFrame( 1 ) ),
                               meanMotion * 1.0E-15 );
            BOOST_CHECK_SMALL( std::fabs( currentRotationalVelocityInBaseFrame( 2 ) - expectedRotationalVelocityVectorInBaseFrame( 2 ) ),
                               meanMotion * 1.0E-15 );

            currentTime += timeStep;
        }

        Eigen::Matrix3d numericalRotationMatrixDerivativeToBaseFrame, numericalRotationMatrixDerivativeToTargetFrame;
        Eigen::Matrix3d upperturbedMatrix, downperturbedMatrix;

        double timePerturbation = 0.1;
        currentTime = startTime + timeStep;
        while ( currentTime < endTime )
        {
            currentRotationMatrixDerivativeToBaseFrame = phobosRotationalEphemeris->getDerivativeOfRotationToBaseFrame( currentTime );
            currentRotationMatrixDerivativeToTargetFrame = phobosRotationalEphemeris->getDerivativeOfRotationToTargetFrame( currentTime );

            upperturbedMatrix = phobosRotationalEphemeris->getRotationToBaseFrame( currentTime + timePerturbation ).toRotationMatrix( );
            downperturbedMatrix = phobosRotationalEphemeris->getRotationToBaseFrame( currentTime - timePerturbation ).toRotationMatrix( );
            numericalRotationMatrixDerivativeToBaseFrame = ( upperturbedMatrix - downperturbedMatrix ) / ( 2.0 * timePerturbation );

            for( unsigned int i = 0; i < 3; i++ )
            {
                for( unsigned int j = 0; j < 3; j++ )
                {
                    BOOST_CHECK_SMALL(
                                std::fabs( numericalRotationMatrixDerivativeToBaseFrame( i, j ) -
                                           currentRotationMatrixDerivativeToBaseFrame( i, j ) ), 1.0E-12 );
                }
            }

            upperturbedMatrix = phobosRotationalEphemeris->getRotationToTargetFrame( currentTime + timePerturbation ).toRotationMatrix( );
            downperturbedMatrix = phobosRotationalEphemeris->getRotationToTargetFrame( currentTime - timePerturbation ).toRotationMatrix( );
            numericalRotationMatrixDerivativeToTargetFrame = ( upperturbedMatrix - downperturbedMatrix ) / ( 2.0 * timePerturbation );

            for( unsigned int i = 0; i < 3; i++ )
            {
                for( unsigned int j = 0; j < 3; j++ )
                {
                    BOOST_CHECK_SMALL(
                                std::fabs( numericalRotationMatrixDerivativeToTargetFrame( i, j ) -
                                           currentRotationMatrixDerivativeToTargetFrame( i, j ) ), 1.0E-12 );
                }
            }

            currentTime += timeStep;

        }
    }
}

BOOST_AUTO_TEST_CASE( testSimpleRotationalDynamicsPropagationWithObliquity )
{
    //Load spice kernels.
    std::string kernelsPath = input_output::getSpiceKernelPath( );
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de-403-masses.tpc");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "naif0009.tls");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "pck00009.tpc");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de421.bsp");

    double phobosSemiMajorAxis = 9376.0E3;

    NamedBodyMap bodyMap = getTestBodyMap( 9376.0E3, 1 );

    double initialEphemerisTime = 1.0E7;
    double finalEphemerisTime = initialEphemerisTime + 10.0 * 86400.0;

    // Set accelerations between bodies that are to be taken into account.
    SelectedTorqueMap torqueMap;
    std::vector< std::string > bodiesToIntegrate;
    bodiesToIntegrate.push_back( "Phobos" );

    double meanMotion = std::sqrt( getBodyGravitationalParameter( "Mars" ) /
                                   std::pow( phobosSemiMajorAxis, 3.0 ) );

    std::cout<<"Mean motion: "<<meanMotion<<std::endl;

    Eigen::Quaterniond nominalInitialRotation = Eigen::Quaterniond( 1.0, 0.0, 0.0, 0.0 );
    double initialObliquity = 0.0 * mathematical_constants::PI / 180.0;

    Eigen::VectorXd systemInitialState = Eigen::VectorXd::Zero( 7 );
    systemInitialState.segment( 0, 4 ) = linear_algebra::convertQuaternionToVectorFormat(
                Eigen::AngleAxisd( -initialObliquity, Eigen::Vector3d::UnitX( ) ) * nominalInitialRotation );
    double initialXAngularVelocity = 0.1 * meanMotion;
    systemInitialState( 4 ) = initialXAngularVelocity;
    systemInitialState( 5 ) = 0.0 * meanMotion;
    systemInitialState( 6 ) = meanMotion;

    basic_astrodynamics::TorqueModelMap torqueModelMap = createTorqueModelsMap(
                bodyMap, torqueMap );

    boost::shared_ptr< numerical_integrators::IntegratorSettings< > > integratorSettings =
            boost::make_shared< RungeKuttaVariableStepSizeSettings< > >
            ( rungeKuttaVariableStepSize,
              initialEphemerisTime, 60.0,
              RungeKuttaCoefficients::rungeKuttaFehlberg78,
              30.0, 300.0, 1.0E-14, 1.0E-14 );

    boost::shared_ptr< RotationalStatePropagatorSettings< double > > propagatorSettings =
            boost::make_shared< RotationalStatePropagatorSettings< double > >
            ( torqueModelMap, bodiesToIntegrate, systemInitialState, boost::make_shared< PropagationTimeTerminationSettings >(
                  finalEphemerisTime ) );

    SingleArcDynamicsSimulator< double > dynamicsSimulator(
                bodyMap, integratorSettings, propagatorSettings, true, false, true );

    boost::shared_ptr< RotationalEphemeris > phobosEphemeris = bodyMap.at( "Phobos" )->getRotationalEphemeris( );

    double startTime = initialEphemerisTime + 3600.0;
    double endTime = finalEphemerisTime - 3600.0;
    double currentTime = startTime;// + 3600.0;

    boost::shared_ptr< RotationalEphemeris > phobosRotationalEphemeris = bodyMap[ "Phobos" ]->getRotationalEphemeris( );

    Eigen::Vector3d currentRotationalVelocityInBaseFrame, currentRotationalVelocityInTargetFrame, indirectRotationalVelocityInBaseFrame;
    Eigen::Matrix3d currentRotationMatrixToBaseFrame, currentRotationMatrixToTargetFrame;
    Eigen::Matrix3d currentRotationMatrixDerivativeToBaseFrame, currentRotationMatrixDerivativeToTargetFrame;
    Eigen::Matrix3d numericalRotationMatrixDerivativeToBaseFrame, numericalRotationMatrixDerivativeToTargetFrame;

    Eigen::Vector3d eulerAngles, calculatedEulerAngleRates, expectedEulerAngleRates;
    Eigen::Matrix3d upperturbedMatrix, downperturbedMatrix;

    double eulerFrequency = ( 0.5024 - 0.4265 ) / 0.4265 * meanMotion;
    double eulerPhase = 0.0;
    double timeStep = ( endTime - startTime ) / 20.0;

    double timePerturbation = 0.1;
    while ( currentTime < endTime )
    {
        double currentAngle = meanMotion * ( currentTime - initialEphemerisTime );

        currentRotationalVelocityInTargetFrame = phobosRotationalEphemeris->getRotationalVelocityVectorInTargetFrame( currentTime );

        BOOST_CHECK_SMALL( currentRotationalVelocityInTargetFrame( 0 ) - initialXAngularVelocity * std::cos( eulerFrequency * ( currentTime - initialEphemerisTime ) ), 1.0E-15 );
        BOOST_CHECK_SMALL( currentRotationalVelocityInTargetFrame( 1 ) - initialXAngularVelocity * std::sin( eulerFrequency * ( currentTime - initialEphemerisTime ) ), 1.0E-15 );
        BOOST_CHECK_CLOSE_FRACTION( currentRotationalVelocityInTargetFrame( 2 ), meanMotion, 1.0E-15 );

        currentRotationMatrixDerivativeToBaseFrame = phobosRotationalEphemeris->getDerivativeOfRotationToBaseFrame( currentTime );
        currentRotationMatrixDerivativeToTargetFrame = phobosRotationalEphemeris->getDerivativeOfRotationToTargetFrame( currentTime );

        upperturbedMatrix = phobosRotationalEphemeris->getRotationToBaseFrame( currentTime + timePerturbation ).toRotationMatrix( );
        downperturbedMatrix = phobosRotationalEphemeris->getRotationToBaseFrame( currentTime - timePerturbation ).toRotationMatrix( );
        numericalRotationMatrixDerivativeToBaseFrame = ( upperturbedMatrix - downperturbedMatrix ) / ( 2.0 * timePerturbation );

        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL(
                            std::fabs( numericalRotationMatrixDerivativeToBaseFrame( i, j ) -
                                       currentRotationMatrixDerivativeToBaseFrame( i, j ) ), 1.0E-12 );
            }
        }

        upperturbedMatrix = phobosRotationalEphemeris->getRotationToTargetFrame( currentTime + timePerturbation ).toRotationMatrix( );
        downperturbedMatrix = phobosRotationalEphemeris->getRotationToTargetFrame( currentTime - timePerturbation ).toRotationMatrix( );
        numericalRotationMatrixDerivativeToTargetFrame = ( upperturbedMatrix - downperturbedMatrix ) / ( 2.0 * timePerturbation );

        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL(
                            std::fabs( numericalRotationMatrixDerivativeToTargetFrame( i, j ) -
                                       currentRotationMatrixDerivativeToTargetFrame( i, j ) ), 1.0E-12 );
            }
        }

//        indirectRotationalVelocityInBaseFrame =
//                getRotationalVelocityVectorInBaseFrameFromMatrices(
//                    phobosRotationalEphemeris->getRotationToTargetFrame( currentTime ).toRotationMatrix( ),
//                    phobosRotationalEphemeris->getDerivativeOfRotationToBaseFrame( currentTime ) );

//        currentRotationalVelocityInBaseFrame = phobosRotationalEphemeris->getRotationalVelocityVectorInBaseFrame( currentTime );

        currentTime += timeStep;
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}

}

