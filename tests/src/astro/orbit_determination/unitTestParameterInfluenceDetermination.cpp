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

#include <boost/test/unit_test.hpp>

#include "tudat/simulation/estimation_setup/simulateObservations.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"
#include "tudat/simulation/estimation_setup/determinePostFitParameterInfluence.h"
#include "tudat/simulation/simulation.h"
#include "tudat/astro/aerodynamics/testApolloCapsuleCoefficients.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/math/statistics/randomVariableGenerator.h"

namespace tudat
{

namespace unit_tests
{

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
using namespace estimatable_parameters;

BOOST_AUTO_TEST_SUITE( test_parameter_influence_estimtation )

//! This test checks whether the absorption of the solar J2 in propagation of solar system bodies is done correctly.
BOOST_AUTO_TEST_CASE( test_ParameterPostFitResiduals )
{
    // Set list of bodies for which analysis is to be performed
    std::vector< std::string > targetBodies = { "Mercury", "Venus" };

    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation times
    double simulationStartEpoch = 0.0;
    double simulationEndEpoch = 4.0 * physical_constants::JULIAN_YEAR;

    // Create body objects.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Mercury" );
    bodiesToCreate.push_back( "Venus" );
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Moon" );
    bodiesToCreate.push_back( "Mars" );
    bodiesToCreate.push_back( "Jupiter" );

    std::vector< Eigen::VectorXd > parameterCorrections;

    // Test for full list of taret bodies, and then for the full list at the same time
    for( unsigned int testCase = 0; testCase < targetBodies.size( ) + 1; testCase++ )
    {
        BodyListSettings bodySettings =
                getDefaultBodySettings( bodiesToCreate );

        // Set solar J2 value
        double sunNormalizedJ2 = 2.0E-7 / calculateLegendreGeodesyNormalizationFactor( 2, 0 );
        bodySettings.at( "Sun" )->gravityFieldSettings = std::make_shared< SphericalHarmonicsGravityFieldSettings >(
                    getBodyGravitationalParameter( "Sun" ), 695.7E6,
                    ( Eigen::Matrix3d( ) << 1.0, 0.0, 0.0,
                      0.0, 0.0, 0.0,
                      sunNormalizedJ2 , 0.0, 0.0 ).finished( ),
                    Eigen::Matrix3d::Zero( ), "IAU_Sun" );

//        // Setting approximate ephemeris for Jupiter to prevent having to use large Spice kernel.
//        double jupiterGravitationalParameter = getBodyGravitationalParameter( "Jupiter" );
//        bodySettings.at( "Jupiter" )->ephemerisSettings = std::make_shared< KeplerEphemerisSettings >(
//                    convertCartesianToKeplerianElements(
//                        getBodyCartesianStateAtEpoch(
//                            "Jupiter", "SSB", "ECLIPJ2000", "None", simulationStartEpoch ),
//                        jupiterGravitationalParameter ), simulationStartEpoch, jupiterGravitationalParameter );

        // Update environment settings of target body
        std::vector< std::string > currentTargetBodies;
        if( testCase < 2 )
        {
            currentTargetBodies.push_back( targetBodies.at( testCase ) );
        }
        else
        {
            currentTargetBodies = targetBodies;
        }

        for( unsigned int body = 0; body < currentTargetBodies.size( ); body++ )
        {
            bodySettings.at( currentTargetBodies.at( body ) )->ephemerisSettings = std::make_shared< InterpolatedSpiceEphemerisSettings >(
                        simulationStartEpoch - 86400.0, simulationStartEpoch + 86400.0, 300.0, "SSB", "ECLIPJ2000" );
        }

        // Create bodies
        SystemOfBodies bodies = createSystemOfBodies( bodySettings );

        
        

        // Define acceleration model settings and propagated bodies and origins
        SelectedAccelerationMap accelerationMap;
        std::vector< std::string > bodiesToPropagate;
        std::vector< std::string > centralBodies;
        for( unsigned int body = 0; body < currentTargetBodies.size( ); body++ )
        {
            bodiesToPropagate.push_back( currentTargetBodies.at( body ) );
            if( currentTargetBodies.at( body ) != "Moon" )
            {
                centralBodies.push_back( "SSB" );
            }
            else
            {
                centralBodies.push_back( "Earth" );
            }
        }

        for( unsigned int body = 0; body < currentTargetBodies.size( ); body++ )
        {
        for( unsigned int j = 0; j < bodiesToCreate.size( ); j++ )
        {
            if( bodiesToCreate.at( j ) != currentTargetBodies.at( body ) )
            {
                if( bodiesToCreate.at( j )  != "Sun" )
                {
                    accelerationMap[ currentTargetBodies.at( body ) ][ bodiesToCreate.at( j ) ].push_back(
                                std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );
                }
                else
                {
                    accelerationMap[ currentTargetBodies.at( body ) ][ bodiesToCreate.at( j ) ].push_back(
                                std::make_shared< SphericalHarmonicAccelerationSettings >( 2, 0 ) );
                }
            }
        }
        }

        // Create acceleration models and propagation settings.
        basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodies, accelerationMap, bodiesToPropagate, centralBodies );
        Eigen::VectorXd systemInitialState = getInitialStatesOfBodies(
                    bodiesToPropagate, centralBodies, bodies, simulationStartEpoch );
        std::shared_ptr< PropagatorSettings< double > > propagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, simulationEndEpoch, cowell );

        // Create integrator settings.
        std::shared_ptr< IntegratorSettings< > > integratorSettings =
                std::make_shared< RungeKuttaVariableStepSizeSettings< > >
                ( simulationStartEpoch, 12.0 * 3600.0,
                  CoefficientSets::rungeKuttaFehlberg78,
                  3.0 * 3600.0, 12.0 * 3600.0, 1.0E-12, 1.0E-12 );

        // Create settings for parameter that is to be perturbed
        std::shared_ptr< EstimatableParameterSettings > perturbedParameterSettings;
        perturbedParameterSettings = (
                    std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                        2, 0, 2, 0, "Sun", spherical_harmonics_cosine_coefficient_block ) );

        // Generate fit for observations with J2 to model without J2
        std::pair< std::shared_ptr< PodOutput< double > >, Eigen::VectorXd > estimationOutput =
                determinePostfitParameterInfluence(
                    bodies, integratorSettings, propagatorSettings, perturbedParameterSettings,
                    6.0 * 3600.0, { -sunNormalizedJ2 }, { 0 } );

        // Get pre- and postfit residuals with RMS
        Eigen::VectorXd prefitResiduals = estimationOutput.first->residualHistory_.at( 0 );
        Eigen::VectorXd postfitResiduals = estimationOutput.first->residuals_;
        double prefitRms = linear_algebra::getVectorEntryRootMeanSquare(
                    prefitResiduals );
        double postfitRms = linear_algebra::getVectorEntryRootMeanSquare(
                    postfitResiduals );

        // Test behaviour of pre- and postfit residuals.
        int observationVectorSize = prefitResiduals.rows( );
        Eigen::Vector3d initialPrefitDifference = prefitResiduals.segment( 0, 3 );
        Eigen::Vector3d initialPostfitDifference = postfitResiduals.segment( 0, 3 );

        Eigen::Vector3d finalPrefitDifference = prefitResiduals.segment( observationVectorSize - 3, 3 );
        Eigen::Vector3d finalPostfitDifference = postfitResiduals.segment( observationVectorSize - 3, 3 );

        BOOST_CHECK_EQUAL( initialPrefitDifference.norm( ) < 1.0E-2, true );

        if( testCase == 0 )
        {
            BOOST_CHECK_EQUAL( finalPrefitDifference.norm( ) > 400.0, true );

            BOOST_CHECK_EQUAL( initialPostfitDifference.norm( ) > 40.0, true );
            BOOST_CHECK_EQUAL( initialPostfitDifference.norm( ) < 60.0, true );

            BOOST_CHECK_EQUAL( finalPostfitDifference.norm( ) > 40.0, true );
            BOOST_CHECK_EQUAL( finalPostfitDifference.norm( ) < 60.0, true );

            BOOST_CHECK_EQUAL( prefitRms / postfitRms > 8.0, true );

        }
        else if( testCase == 1 )
        {
            BOOST_CHECK_EQUAL( finalPrefitDifference.norm( ) > 100.0, true );

            BOOST_CHECK_EQUAL( initialPostfitDifference.norm( ) > 0.75, true );
            BOOST_CHECK_EQUAL( initialPostfitDifference.norm( ) < 1.0, true );

            BOOST_CHECK_EQUAL( finalPostfitDifference.norm( ) > 0.75, true );
            BOOST_CHECK_EQUAL( finalPostfitDifference.norm( ) < 1.0, true );

            BOOST_CHECK_EQUAL( prefitRms / postfitRms > 8.0, true );
        }
        else if( testCase == 2 )
        {
            Eigen::Vector6d mercuryPositionAdjustmentDifference =
                    parameterCorrections.at( 0 ).segment( 0, 6 ) - estimationOutput.second.segment( 0, 6 );
            Eigen::Vector6d marsPositionAdjustmentDifference =
                    parameterCorrections.at( 1 ).segment( 0, 6 ) - estimationOutput.second.segment( 6, 6 );
            for( unsigned int index = 0; index < 3; index++ )
            {
                BOOST_CHECK_SMALL( std::fabs( mercuryPositionAdjustmentDifference( index ) ), 2.0E-1 );
                BOOST_CHECK_SMALL( std::fabs( mercuryPositionAdjustmentDifference( index + 3 ) ), 2.0E-7 );
                BOOST_CHECK_SMALL( std::fabs( marsPositionAdjustmentDifference( index ) ), 1.0E-1 );
                BOOST_CHECK_SMALL( std::fabs( marsPositionAdjustmentDifference( index + 3 ) ), 1.0E-7 );
            }
        }

        std::cout << "Parameter difference " << estimationOutput.second.transpose( ) << std::endl << std::endl;

        parameterCorrections.push_back( estimationOutput.second );
    }
}

//! This test checks whether the absorption of the solar J2 in re-entry is done correctly.
BOOST_AUTO_TEST_CASE( test_ParameterPostFitResidualsApollo )
{

    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation start and end epoch.
    const double simulationStartEpoch = 0.0;

    // Define simulation body settings.
    BodyListSettings bodySettings =
            getDefaultBodySettings( { "Earth" }, "SSB", "J2000"  );
    bodySettings.at( "Earth" )->ephemerisSettings = std::make_shared< simulation_setup::ConstantEphemerisSettings >(
                Eigen::Vector6d::Zero( ), "SSB", "J2000" );
    bodySettings.at( "Earth" )->rotationModelSettings->resetOriginalFrame( "J2000" );

    // Create Earth object
    simulation_setup::SystemOfBodies bodies = simulation_setup::createSystemOfBodies( bodySettings );

    // Retrieve Earth J2
    double earthC20 =
            std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >(
                bodies.at( "Earth" )->getGravityFieldModel( ) )->getCosineCoefficients( )( 2, 0 );

    // Create vehicle objects.
    bodies.createEmptyBody( "Apollo" );

    bodies.at( "Apollo" )->setAerodynamicCoefficientInterface(
                unit_tests::getApolloCoefficientInterface( ) );
    bodies.at( "Apollo" )->setConstantBodyMass( 5.0E3 );

    double constantAngleOfAttack = 25.0 * mathematical_constants::PI / 180.0;
    std::shared_ptr< ephemerides::RotationalEphemeris > vehicleRotationModel =
            createRotationModel(
                std::make_shared< AerodynamicAngleRotationSettings >(
                    "Earth", "J2000", "VehicleFixed",
                    [=](const double){ return ( Eigen::Vector3d( ) << constantAngleOfAttack, 0.0, 0.0 ).finished( ); } ),
                "Apollo", bodies );

    bodies.at( "Apollo" )->setRotationalEphemeris( vehicleRotationModel );

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Define acceleration model settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfApollo;
    accelerationsOfApollo[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 2, 0 ) );
    accelerationsOfApollo[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( aerodynamic ) );
    accelerationMap[ "Apollo" ] = accelerationsOfApollo;

    bodiesToPropagate.push_back( "Apollo" );
    centralBodies.push_back( "Earth" );

    // Create acceleration models
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodies, accelerationMap, bodiesToPropagate, centralBodies );


    // Set initial spherical elements for Apollo.
    Eigen::Vector6d apolloSphericalEntryState;
    apolloSphericalEntryState( SphericalOrbitalStateElementIndices::radiusIndex ) =
            spice_interface::getAverageRadius( "Earth" ) + 120.0E3;
    apolloSphericalEntryState( SphericalOrbitalStateElementIndices::latitudeIndex ) =
            25.0 * mathematical_constants::PI / 180.0;
    apolloSphericalEntryState( SphericalOrbitalStateElementIndices::longitudeIndex )  =
            25.0 * mathematical_constants::PI / 180.0;
    apolloSphericalEntryState( SphericalOrbitalStateElementIndices::speedIndex ) = 7.7E3;
    apolloSphericalEntryState( SphericalOrbitalStateElementIndices::flightPathIndex ) =
            -1.25 * mathematical_constants::PI / 180.0;
    apolloSphericalEntryState( SphericalOrbitalStateElementIndices::headingAngleIndex )  =
            25.0 * mathematical_constants::PI / 180.0;

    // Convert apollo state from spherical elements to Cartesian elements.
    Eigen::Vector6d systemInitialState = convertSphericalOrbitalToCartesianState(
                apolloSphericalEntryState );
    std::shared_ptr< ephemerides::RotationalEphemeris > earthRotationalEphemeris =
            bodies.at( "Earth" )->getRotationalEphemeris( );
    systemInitialState = transformStateToInertialOrientation( systemInitialState, simulationStartEpoch, earthRotationalEphemeris );

    // Define termination conditions
    std::shared_ptr< SingleDependentVariableSaveSettings > terminationDependentVariable =
            std::make_shared< SingleDependentVariableSaveSettings >(
                altitude_dependent_variable, "Apollo", "Earth" );
    std::shared_ptr< PropagationTerminationSettings > terminationSettings =
            std::make_shared< PropagationDependentVariableTerminationSettings >(
                terminationDependentVariable, 25.0E3, true );

    // Create propagation settings.
    std::shared_ptr< PropagatorSettings< double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState,
              terminationSettings, cowell );

    // Create integrator settings.
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4, simulationStartEpoch, 1.0 );

    // Create settings for parameter that is to be perturbed
    std::shared_ptr< EstimatableParameterSettings > perturbedParameterSettings;
    perturbedParameterSettings = (
                std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                    2, 0, 2, 0, "Earth", spherical_harmonics_cosine_coefficient_block ) );


    // Generate fit for observations with J2 to model without J2
    std::pair< std::shared_ptr< PodOutput< double > >, Eigen::VectorXd > estimationOutput =
            determinePostfitParameterInfluence(
                bodies, integratorSettings, propagatorSettings, perturbedParameterSettings,
                1.0, { -earthC20 }, { 0 } );

    // Get pre- and postfit residuals with RMS
    Eigen::VectorXd prefitResiduals = estimationOutput.first->residualHistory_.at( 0 );
    Eigen::VectorXd postfitResiduals = estimationOutput.first->residuals_;
    double prefitRms = linear_algebra::getVectorEntryRootMeanSquare(
                prefitResiduals );
    double postfitRms = linear_algebra::getVectorEntryRootMeanSquare(
                postfitResiduals );

    // Test behaviour of pre- and postfit residuals.
    int observationVectorSize = prefitResiduals.rows( );
    Eigen::Vector3d initialPrefitDifference = prefitResiduals.segment( 0, 3 );
    Eigen::Vector3d initialPostfitDifference = postfitResiduals.segment( 0, 3 );

    Eigen::Vector3d finalPrefitDifference = prefitResiduals.segment( observationVectorSize - 3, 3 );
    Eigen::Vector3d finalPostfitDifference = postfitResiduals.segment( observationVectorSize - 3, 3 );


    BOOST_CHECK_EQUAL( initialPrefitDifference.norm( ) < 1.0, true );
    BOOST_CHECK_EQUAL( finalPrefitDifference.norm( ) > 20.0E3, true );

    BOOST_CHECK_EQUAL( initialPostfitDifference.norm( ) > 275.0, true );
    BOOST_CHECK_EQUAL( initialPostfitDifference.norm( ) < 325.0, true );

    BOOST_CHECK_EQUAL( finalPostfitDifference.norm( ) > 275.0, true );
    BOOST_CHECK_EQUAL( finalPostfitDifference.norm( ) < 325.0, true );

    BOOST_CHECK_EQUAL( prefitRms / postfitRms > 40.0, true );

    std::cout << "Parameter difference " << estimationOutput.second.transpose( ) << std::endl;

    BOOST_AUTO_TEST_SUITE_END( )

}

}

}
