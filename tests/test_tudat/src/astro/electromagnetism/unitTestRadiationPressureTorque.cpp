/*    Copyright (c) 2010-2026, Delft University of Technology
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

#include "tudat/math/basic/linearAlgebra.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/propagation_setup/createTorqueModel.h"
#include "tudat/simulation/propagation_setup/singleArcDynamicsSimulator.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::simulation_setup;
using namespace tudat::propagators;
using namespace tudat::numerical_integrators;

BOOST_AUTO_TEST_SUITE( test_radiation_pressure_torque )

//! Propagate the rotation of a cannonball and a panelled target under the radiation pressure torque of an isotropic
//! source, and compare the torque with the force about the center of mass, in the body-fixed frame.
BOOST_AUTO_TEST_CASE( testRadiationPressureTorquePropagation )
{
    const double luminosity = 3.828E26;
    const double distance = physical_constants::ASTRONOMICAL_UNIT;
    const double radiationPressure =
            luminosity / ( 4.0 * mathematical_constants::PI * distance * distance * physical_constants::SPEED_OF_LIGHT );

    const Eigen::Vector3d centerOfMass( 0.0, 0.05, 0.5 );
    const Eigen::Vector3d panelPosition( 0.3, 0.2, -0.1 );
    const double area = 2.0;
    const double radiationPressureCoefficient = 1.5;

    for( const bool usePanelledTarget : { false, true } )
    {
        BodyListSettings bodySettings( "SSB", "ECLIPJ2000" );
        bodySettings.addSettings( "Sun" );
        bodySettings.at( "Sun" )->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >( Eigen::Vector6d::Zero( ) );
        bodySettings.at( "Sun" )->radiationSourceModelSettings =
                isotropicPointRadiationSourceModelSettings( constantLuminosityModelSettings( luminosity ) );

        bodySettings.addSettings( "Satellite" );
        bodySettings.at( "Satellite" )->ephemerisSettings =
                std::make_shared< ConstantEphemerisSettings >( ( Eigen::Vector6d( ) << distance, 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( ) );
        bodySettings.at( "Satellite" )->rigidBodyPropertiesSettings =
                constantRigidBodyPropertiesSettings( 100.0, centerOfMass, Eigen::Vector3d( 10.0, 20.0, 30.0 ).asDiagonal( ) );
        if( usePanelledTarget )
        {
            // Required for the panel geometry; the propagated rotation replaces it
            bodySettings.at( "Satellite" )->rotationModelSettings =
                    constantRotationModelSettings( "ECLIPJ2000", "Satellite_Fixed", Eigen::Quaterniond::Identity( ) );
            // One absorbing panel, facing the source along the body-fixed -x axis
            bodySettings.at( "Satellite" )->bodyExteriorPanelSettings_ = fullPanelledBodySettings(
                    { bodyPanelSettings( frameFixedPanelGeometry( -Eigen::Vector3d::UnitX( ), panelPosition, area ),
                                         specularDiffuseBodyPanelReflectionLawSettings( 0.0, 0.0, false ) ) } );
            bodySettings.at( "Satellite" )->radiationPressureTargetModelSettings = paneledRadiationPressureTargetModelSettings( );
        }
        else
        {
            bodySettings.at( "Satellite" )->radiationPressureTargetModelSettings =
                    cannonballRadiationPressureTargetModelSettings( area, radiationPressureCoefficient );
        }
        SystemOfBodies bodies = createSystemOfBodies( bodySettings );

        SelectedTorqueMap torqueSettings;
        torqueSettings[ "Satellite" ][ "Sun" ].push_back(
                std::make_shared< TorqueSettings >( basic_astrodynamics::radiation_pressure_torque ) );
        basic_astrodynamics::TorqueModelMap torqueModels = createTorqueModelsMap( bodies, torqueSettings, { "Satellite" } );

        // Initial attitude and spin such that the attitude changes by about 0.6 rad
        Eigen::VectorXd initialState = Eigen::VectorXd::Zero( 7 );
        initialState.segment( 0, 4 ) = linear_algebra::convertQuaternionToVectorFormat( Eigen::Quaterniond(
                Eigen::AngleAxisd( 0.3, Eigen::Vector3d::UnitZ( ) ) * Eigen::AngleAxisd( 0.2, Eigen::Vector3d::UnitX( ) ) ) );
        initialState.segment( 4, 3 ) = Eigen::Vector3d( 2.0E-4, -1.0E-4, 1.0E-3 );

        std::shared_ptr< RotationalStatePropagatorSettings< double > > propagatorSettings =
                std::make_shared< RotationalStatePropagatorSettings< double > >(
                        torqueModels,
                        std::vector< std::string >{ "Satellite" },
                        initialState,
                        0.0,
                        std::make_shared< IntegratorSettings<> >( rungeKutta4, 10.0 ),
                        std::make_shared< PropagationTimeTerminationSettings >( 600.0 ),
                        quaternions,
                        std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > >{
                                std::make_shared< SingleTorqueDependentVariableSaveSettings >(
                                        basic_astrodynamics::radiation_pressure_torque, "Satellite", "Sun" ) } );

        SingleArcDynamicsSimulator<> dynamicsSimulator( bodies, propagatorSettings );
        const std::map< double, Eigen::VectorXd > stateHistory = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
        const std::map< double, Eigen::VectorXd > torqueHistory = dynamicsSimulator.getDependentVariableHistory( );
        BOOST_CHECK_EQUAL( torqueHistory.size( ), 61 );

        for( const auto& [ time, torque ] : torqueHistory )
        {
            // Direction of the radiation in the body-fixed frame
            const Eigen::Matrix3d rotationToInertialFrame =
                    linear_algebra::convertVectorQuaternionToMatrixFormat( stateHistory.at( time ).segment( 0, 4 ) );
            const Eigen::Vector3d radiationDirection = rotationToInertialFrame.transpose( ) * Eigen::Vector3d::UnitX( );

            Eigen::Vector3d expectedTorque;
            if( usePanelledTarget )
            {
                const Eigen::Vector3d panelForce = radiationPressure * area * radiationDirection.x( ) * radiationDirection;
                expectedTorque = ( panelPosition - centerOfMass ).cross( panelForce );
            }
            else
            {
                const Eigen::Vector3d force = radiationPressure * radiationPressureCoefficient * area * radiationDirection;
                expectedTorque = -centerOfMass.cross( force );
            }

            for( unsigned int i = 0; i < 3; i++ )
            {
                BOOST_CHECK_SMALL( torque( i ) - expectedTorque( i ), 1.0E-12 * expectedTorque.norm( ) );
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
