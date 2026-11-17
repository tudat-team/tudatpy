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
#include "tudat/simulation/estimation_setup/createEstimatableParameters.h"
#include "tudat/simulation/estimation_setup/singleArcVariationalEquationsSolver.h"
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

//! Propagate the orbit and rotation of a cannonball and a panelled target under the radiation pressure of an isotropic source,
//! with scaled radiation pressure, and compare the torque with the scaled force about the center of mass, in the body-fixed frame.
BOOST_AUTO_TEST_CASE( testRadiationPressureTorquePropagation )
{
    const double luminosity = 3.828E26;
    const double distance = physical_constants::ASTRONOMICAL_UNIT;
    const double sunGravitationalParameter = 1.32712440018E20;

    const Eigen::Vector3d centerOfMass( 0.0, 0.05, 0.5 );
    const Eigen::Vector3d panelPosition( 0.3, 0.2, -0.1 );
    const double area = 2.0;
    const double radiationPressureCoefficient = 1.5;
    const double specularReflectivity = 0.5;
    const double sourceDirectionScaling = 1.3;
    const double perpendicularScaling = 0.7;

    for( const bool usePanelledTarget : { false, true } )
    {
        BodyListSettings bodySettings( "SSB", "ECLIPJ2000" );
        bodySettings.addSettings( "Sun" );
        bodySettings.at( "Sun" )->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >( Eigen::Vector6d::Zero( ) );
        bodySettings.at( "Sun" )->gravityFieldSettings = centralGravitySettings( sunGravitationalParameter );
        bodySettings.at( "Sun" )->radiationSourceModelSettings =
                isotropicPointRadiationSourceModelSettings( constantLuminosityModelSettings( luminosity ) );

        bodySettings.addSettings( "Satellite" );
        bodySettings.at( "Satellite" )->rigidBodyPropertiesSettings =
                constantRigidBodyPropertiesSettings( 100.0, centerOfMass, Eigen::Vector3d( 10.0, 20.0, 30.0 ).asDiagonal( ) );
        if( usePanelledTarget )
        {
            // Required for the panel geometry; the propagated rotation replaces it
            bodySettings.at( "Satellite" )->rotationModelSettings =
                    constantRotationModelSettings( "ECLIPJ2000", "Satellite_Fixed", Eigen::Quaterniond::Identity( ) );

            // One panel, facing the source along the body-fixed -x axis, with a force component perpendicular to the source direction
            bodySettings.at( "Satellite" )->bodyExteriorPanelSettings_ = fullPanelledBodySettings(
                    { bodyPanelSettings( frameFixedPanelGeometry( -Eigen::Vector3d::UnitX( ), panelPosition, area ),
                                         specularDiffuseBodyPanelReflectionLawSettings( specularReflectivity, 0.0, false ) ) } );
            bodySettings.at( "Satellite" )->radiationPressureTargetModelSettings = paneledRadiationPressureTargetModelSettings( );
        }
        else
        {
            bodySettings.at( "Satellite" )->radiationPressureTargetModelSettings =
                    cannonballRadiationPressureTargetModelSettings( area, radiationPressureCoefficient );
        }
        SystemOfBodies bodies = createSystemOfBodies( bodySettings );

        // Translational dynamics, about a circular orbit
        SelectedAccelerationMap accelerationSettings;
        accelerationSettings[ "Satellite" ][ "Sun" ] = {
            std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ),
            std::make_shared< AccelerationSettings >( basic_astrodynamics::radiation_pressure )
        };
        Eigen::VectorXd initialTranslationalState = Eigen::VectorXd::Zero( 6 );
        initialTranslationalState( 0 ) = distance;
        initialTranslationalState( 4 ) = std::sqrt( sunGravitationalParameter / distance );
        std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >(
                        std::vector< std::string >{ "Sun" },
                        createAccelerationModelsMap( bodies, accelerationSettings, { "Satellite" }, { "Sun" } ),
                        std::vector< std::string >{ "Satellite" },
                        initialTranslationalState,
                        600.0 );

        // Rotational dynamics, from an initial attitude and spin such that the attitude changes by about 0.6 rad
        SelectedTorqueMap torqueSettings;
        torqueSettings[ "Satellite" ][ "Sun" ].push_back(
                std::make_shared< TorqueSettings >( basic_astrodynamics::radiation_pressure_torque ) );
        Eigen::VectorXd initialRotationalState = Eigen::VectorXd::Zero( 7 );
        initialRotationalState.segment( 0, 4 ) = linear_algebra::convertQuaternionToVectorFormat( Eigen::Quaterniond(
                Eigen::AngleAxisd( 0.3, Eigen::Vector3d::UnitZ( ) ) * Eigen::AngleAxisd( 0.2, Eigen::Vector3d::UnitX( ) ) ) );
        initialRotationalState.segment( 4, 3 ) = Eigen::Vector3d( 2.0E-4, -1.0E-4, 1.0E-3 );
        std::shared_ptr< RotationalStatePropagatorSettings< double > > rotationalPropagatorSettings =
                std::make_shared< RotationalStatePropagatorSettings< double > >(
                        createTorqueModelsMap( bodies, torqueSettings, { "Satellite" } ),
                        std::vector< std::string >{ "Satellite" },
                        initialRotationalState,
                        std::make_shared< PropagationTimeTerminationSettings >( 600.0 ) );

        std::shared_ptr< MultiTypePropagatorSettings< double > > propagatorSettings =
                std::make_shared< MultiTypePropagatorSettings< double > >(
                        std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > >{ translationalPropagatorSettings,
                                                                                                 rotationalPropagatorSettings },
                        std::make_shared< IntegratorSettings<> >( rungeKutta4, 10.0 ),
                        0.0,
                        std::make_shared< PropagationTimeTerminationSettings >( 600.0 ),
                        std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > >{
                                std::make_shared< SingleTorqueDependentVariableSaveSettings >(
                                        basic_astrodynamics::radiation_pressure_torque, "Satellite", "Sun" ) } );

        // Scale the radiation pressure through the estimatable parameters, which also scale the torque
        std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parameters =
                createParametersToEstimate< double, double >(
                        { estimatable_parameters::radiationPressureTargetDirectionScaling( "Satellite", "Sun" ),
                          estimatable_parameters::radiationPressureTargetPerpendicularDirectionScaling( "Satellite", "Sun" ) },
                        bodies,
                        propagatorSettings );
        parameters->resetParameterValues( Eigen::VectorXd( Eigen::Vector2d( sourceDirectionScaling, perpendicularScaling ) ) );

        SingleArcDynamicsSimulator<> dynamicsSimulator( bodies, propagatorSettings );
        const std::map< double, Eigen::VectorXd > stateHistory = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
        const std::map< double, Eigen::VectorXd > torqueHistory = dynamicsSimulator.getDependentVariableHistory( );
        BOOST_CHECK_EQUAL( torqueHistory.size( ), 61 );

        for( const auto& [ time, torque ] : torqueHistory )
        {
            // Source-target direction and radiation pressure, from the translational state
            const Eigen::Vector3d position = stateHistory.at( time ).segment( 0, 3 );
            const double radiationPressure =
                    luminosity / ( 4.0 * mathematical_constants::PI * position.squaredNorm( ) * physical_constants::SPEED_OF_LIGHT );

            // Source-target direction in the body-fixed frame, from the rotational state
            const Eigen::Matrix3d rotationToInertialFrame =
                    linear_algebra::convertVectorQuaternionToMatrixFormat( stateHistory.at( time ).segment( 6, 4 ) );
            const Eigen::Vector3d radiationDirection = rotationToInertialFrame.transpose( ) * position.normalized( );

            Eigen::Vector3d force;
            Eigen::Vector3d momentArm;
            if( usePanelledTarget )
            {
                const double incidenceCosine = radiationDirection.x( );
                force = radiationPressure * area * incidenceCosine *
                        ( ( 1.0 - specularReflectivity ) * radiationDirection +
                          2.0 * specularReflectivity * incidenceCosine * Eigen::Vector3d::UnitX( ) );
                momentArm = panelPosition - centerOfMass;
            }
            else
            {
                force = radiationPressure * radiationPressureCoefficient * area * radiationDirection;
                momentArm = -centerOfMass;
            }
            const Eigen::Vector3d sourceDirectionForce = force.dot( radiationDirection ) * radiationDirection;
            const Eigen::Vector3d scaledForce =
                    sourceDirectionScaling * sourceDirectionForce + perpendicularScaling * ( force - sourceDirectionForce );
            const Eigen::Vector3d expectedTorque = momentArm.cross( scaledForce );

            for( unsigned int i = 0; i < 3; i++ )
            {
                BOOST_CHECK_SMALL( torque( i ) - expectedTorque( i ), 1.0E-12 * expectedTorque.norm( ) );
            }
        }

        // Sensitivity of the final rotational state to the scaling factors: variational equations, using the torque partials,
        // against central differences of propagations with perturbed scaling factors
        std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > > parameterSettings =
                getInitialStateParameterSettings< double, double >( propagatorSettings, bodies );
        parameterSettings.push_back( estimatable_parameters::radiationPressureTargetDirectionScaling( "Satellite", "Sun" ) );
        parameterSettings.push_back( estimatable_parameters::radiationPressureTargetPerpendicularDirectionScaling( "Satellite", "Sun" ) );
        std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > estimatedParameters =
                createParametersToEstimate< double, double >( parameterSettings, bodies, propagatorSettings );
        const Eigen::VectorXd nominalParameters = estimatedParameters->getFullParameterValues< double >( );

        SingleArcVariationalEquationsSolver< double, double > variationalEquationsSolver( bodies, propagatorSettings, estimatedParameters );
        const Eigen::MatrixXd finalSensitivity =
                variationalEquationsSolver.getSingleArcVariationalPropagationResults( )->getSensitivitySolution( ).rbegin( )->second;

        const double scalingPerturbation = 0.01;
        for( unsigned int j = 0; j < 2; j++ )
        {
            std::vector< Eigen::VectorXd > perturbedFinalStates;
            for( const double sign : { 1.0, -1.0 } )
            {
                Eigen::VectorXd perturbedParameters = nominalParameters;
                perturbedParameters( 13 + j ) += sign * scalingPerturbation;
                estimatedParameters->resetParameterValues( perturbedParameters );
                SingleArcDynamicsSimulator<> perturbedSimulator( bodies, propagatorSettings );
                perturbedFinalStates.push_back( perturbedSimulator.getEquationsOfMotionNumericalSolution( ).rbegin( )->second );
            }
            estimatedParameters->resetParameterValues( nominalParameters );

            const Eigen::VectorXd numericalSensitivity =
                    ( perturbedFinalStates.at( 0 ) - perturbedFinalStates.at( 1 ) ) / ( 2.0 * scalingPerturbation );
            const Eigen::VectorXd rotationalSensitivity = finalSensitivity.block( 6, j, 7, 1 );
            // With an absolute floor, as a cannonball has no force perpendicular to the source direction
            BOOST_CHECK_SMALL( ( rotationalSensitivity - numericalSensitivity.segment( 6, 7 ) ).norm( ),
                               1.0E-6 * numericalSensitivity.segment( 6, 7 ).norm( ) + 1.0E-15 );
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
