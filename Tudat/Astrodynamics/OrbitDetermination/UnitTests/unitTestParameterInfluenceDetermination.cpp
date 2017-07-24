/*    Copyright (c) 2010-2017, Delft University of Technology
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

#include <boost/test/unit_test.hpp>

#include "Tudat/Astrodynamics/ObservationModels/simulateObservations.h"
#include "Tudat/Astrodynamics/OrbitDetermination/orbitDeterminationManager.h"
#include "Tudat/Astrodynamics/OrbitDetermination/determinePostFitParameterInfluence.h"
#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
#include "Tudat/Astrodynamics/Aerodynamics/UnitTests/testApolloCapsuleCoefficients.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createGroundStations.h"
#include "Tudat/Mathematics/Statistics/randomVariableGenerator.h"

namespace tudat
{

namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_parameter_influence_estimtation )

//! This test checks whether the input/output of the estimation (weights, a priori covariance, unscaled covariance) are
//! correctly handed
BOOST_AUTO_TEST_CASE( test_ParameterPostFitResiduals )
{
    using namespace tudat::simulation_setup;
    using namespace tudat::estimatable_parameters;
    using namespace tudat::propagators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::basic_mathematics;
    using namespace tudat::unit_conversions;
    using namespace tudat::ephemerides;
    using namespace tudat::spice_interface;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::vector< std::string > targetBodies = { "Mercury", "Venus" };

    // Load Spice kernels.
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de-403-masses.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "pck00009.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "naif0009.tls" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "planetaryOrbitKernel.tm" );

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

    for( unsigned int testCase = 0; testCase < targetBodies.size( ); testCase++ )
    {
        std::string targetBody = targetBodies.at( testCase );
        std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
                getDefaultBodySettings( bodiesToCreate );

        double sunNormalizedJ2 = 2.0E-7 / calculateLegendreGeodesyNormalizationFactor( 2, 0 );
        bodySettings[ "Sun" ]->gravityFieldSettings = boost::make_shared< SphericalHarmonicsGravityFieldSettings >(
                    getBodyGravitationalParameter( "Sun" ), 695.7E6,
                    ( Eigen::Matrix3d( ) << 1.0, 0.0, 0.0,
                      0.0, 0.0, 0.0,
                      sunNormalizedJ2 , 0.0, 0.0 ).finished( ),
                    Eigen::Matrix3d::Zero( ), "IAU_Sun" );
        bodySettings[ targetBody ]->ephemerisSettings = boost::make_shared< InterpolatedSpiceEphemerisSettings >(
                    simulationStartEpoch - 86400.0, simulationStartEpoch + 86400.0, 300.0, "SSB", "ECLIPJ2000" );

        NamedBodyMap bodyMap = createBodies( bodySettings );

        // Finalize body creation.
        setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );


        SelectedAccelerationMap accelerationMap;
        std::vector< std::string > bodiesToPropagate;
        bodiesToPropagate.push_back( targetBody );
        std::vector< std::string > centralBodies;
        if( targetBody != "Moon" )
        {
            centralBodies.push_back( "SSB" );
        }
        else
        {
            centralBodies.push_back( "Earth" );
        }

        for( unsigned int j = 0; j < bodiesToCreate.size( ); j++ )
        {
            if( bodiesToCreate.at( j ) != targetBody )
            {
                if( bodiesToCreate.at( j )  != "Sun" )
                {
                    accelerationMap[ targetBody ][ bodiesToCreate.at( j ) ].push_back(
                                boost::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
                }
                else
                {
                    accelerationMap[ targetBody ][ bodiesToCreate.at( j ) ].push_back(
                                boost::make_shared< SphericalHarmonicAccelerationSettings >( 2, 0 ) );
                }
            }
        }
        // Create acceleration models and propagation settings.
        basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodyMap, accelerationMap, bodiesToPropagate, centralBodies );


        Eigen::VectorXd systemInitialState = getInitialStatesOfBodies(
                    bodiesToPropagate, centralBodies, bodyMap, simulationStartEpoch );

        boost::shared_ptr< PropagatorSettings< double > > propagatorSettings =
                boost::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, simulationEndEpoch, cowell );

        boost::shared_ptr< IntegratorSettings< > > integratorSettings =
                boost::make_shared< RungeKuttaVariableStepSizeSettings< > >
                ( rungeKuttaVariableStepSize, ( simulationStartEpoch ), 12.0 * 3600.0,
                  RungeKuttaCoefficients::CoefficientSets::rungeKuttaFehlberg78,
                  3.0 * 3600.0, 12.0 * 3600.0, 1.0E-12, 1.0E-12 );

        boost::shared_ptr< EstimatableParameterSettings > perturbedParameterSettings;
        perturbedParameterSettings = (
                    boost::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                        2, 0, 2, 0, "Sun", spherical_harmonics_cosine_coefficient_block ) );

        std::pair< boost::shared_ptr< PodOutput< double > >, Eigen::VectorXd > estimationOutput =
                determinePostfitParameterInfluence(
                    bodyMap, integratorSettings, propagatorSettings, perturbedParameterSettings,
                    6.0 * 3600.0, boost::assign::list_of( -sunNormalizedJ2 ), boost::assign::list_of( 0 ) );
        Eigen::VectorXd prefitResiduals = estimationOutput.first->firstIterationResiduals_;
        Eigen::VectorXd postfitResiduals = estimationOutput.first->residuals_;

        double prefitRms = linear_algebra::getVectorEntryRootMeanSquare(
                    prefitResiduals );
        double postfitRms = linear_algebra::getVectorEntryRootMeanSquare(
                    postfitResiduals );


        int observationVectorSize = prefitResiduals.rows( );
        Eigen::Vector3d initialPrefitDifference = prefitResiduals.segment( 0, 3 );
        Eigen::Vector3d initialPostfitDifference = postfitResiduals.segment( 0, 3 );

        Eigen::Vector3d finalPrefitDifference = prefitResiduals.segment( observationVectorSize - 3, 3 );
        Eigen::Vector3d finalPostfitDifference = postfitResiduals.segment( observationVectorSize - 3, 3 );

        BOOST_CHECK_EQUAL( initialPrefitDifference.norm( ) < 1.0E-2, true );

        if( targetBody == "Mercury" )
        {
            BOOST_CHECK_EQUAL( finalPrefitDifference.norm( ) > 400.0, true );

            BOOST_CHECK_EQUAL( initialPostfitDifference.norm( ) > 50.0, true );
            BOOST_CHECK_EQUAL( initialPostfitDifference.norm( ) < 60.0, true );

            BOOST_CHECK_EQUAL( finalPostfitDifference.norm( ) > 50.0, true );
            BOOST_CHECK_EQUAL( finalPostfitDifference.norm( ) < 60.0, true );

            BOOST_CHECK_EQUAL( prefitRms / postfitRms > 8.0, true );

        }
        else if( targetBody == "Venus" )
        {
            BOOST_CHECK_EQUAL( finalPrefitDifference.norm( ) > 100.0, true );

            BOOST_CHECK_EQUAL( initialPostfitDifference.norm( ) > 0.75, true );
            BOOST_CHECK_EQUAL( initialPostfitDifference.norm( ) < 1.0, true );

            BOOST_CHECK_EQUAL( finalPostfitDifference.norm( ) > 0.75, true );
            BOOST_CHECK_EQUAL( finalPostfitDifference.norm( ) < 1.0, true );

            BOOST_CHECK_EQUAL( prefitRms / postfitRms > 8.0, true );

        }

        std::cout<<"Parameter difference "<<estimationOutput.second.transpose( )<<std::endl<<std::endl;
    }
}

BOOST_AUTO_TEST_CASE( test_ParameterPostFitResidualsApollo )
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
    using namespace estimatable_parameters;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ENVIRONMENT            //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Load Spice kernels.
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "pck00009.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de-403-masses.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de421.bsp" );

    // Set simulation start epoch.
    const double simulationStartEpoch = 0.0;

    // Set simulation end epoch.
    const double simulationEndEpoch = 3100.0;

    // Set numerical integration fixed step size.
    const double fixedStepSize = 1.0;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ENVIRONMENT            //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define simulation body settings.
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( { "Earth" }, simulationStartEpoch - 10.0 * fixedStepSize,
                                    simulationEndEpoch + 10.0 * fixedStepSize );
    bodySettings[ "Earth" ]->ephemerisSettings = boost::make_shared< simulation_setup::ConstantEphemerisSettings >(
                Eigen::Vector6d::Zero( ), "SSB", "J2000" );
    bodySettings[ "Earth" ]->rotationModelSettings->resetOriginalFrame( "J2000" );

    // Create Earth object
    simulation_setup::NamedBodyMap bodyMap = simulation_setup::createBodies( bodySettings );

    double earthC20 =
            boost::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >(
                bodyMap.at( "Earth" )->getGravityFieldModel( ) )->getCosineCoefficients( )( 2, 0 );
    // Create vehicle objects.
    bodyMap[ "Apollo" ] = boost::make_shared< simulation_setup::Body >( );
    bodyMap[ "Apollo" ]->setEphemeris( boost::make_shared< TabulatedCartesianEphemeris< > >(
                                           boost::shared_ptr< interpolators::OneDimensionalInterpolator
                                           < double, Eigen::Vector6d > >( ), "Earth", "J2000" ) );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Create vehicle aerodynamic coefficients
    bodyMap[ "Apollo" ]->setAerodynamicCoefficientInterface(
                unit_tests::getApolloCoefficientInterface( ) );
    bodyMap[ "Apollo" ]->setConstantBodyMass( 5.0E3 );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE ACCELERATIONS            ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Define acceleration model settings.
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfApollo;
    accelerationsOfApollo[ "Earth" ].push_back( boost::make_shared< SphericalHarmonicAccelerationSettings >( 2, 0 ) );
    accelerationsOfApollo[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( aerodynamic ) );
    accelerationMap[  "Apollo" ] = accelerationsOfApollo;

    bodiesToPropagate.push_back( "Apollo" );
    centralBodies.push_back( "Earth" );

    // Create acceleration models
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    // Define constant 30 degree angle of attack
    double constantAngleOfAttack = 25.0 * mathematical_constants::PI / 180.0;
    bodyMap.at( "Apollo" )->getFlightConditions( )->getAerodynamicAngleCalculator( )->setOrientationAngleFunctions(
                boost::lambda::constant( constantAngleOfAttack ) );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set spherical elements for Apollo.
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

    boost::shared_ptr< ephemerides::RotationalEphemeris > earthRotationalEphemeris =
            bodyMap.at( "Earth" )->getRotationalEphemeris( );
    systemInitialState = transformStateToGlobalFrame( systemInitialState, simulationStartEpoch, earthRotationalEphemeris );

    // Define termination conditions
    boost::shared_ptr< SingleDependentVariableSaveSettings > terminationDependentVariable =
            boost::make_shared< SingleDependentVariableSaveSettings >(
                altitude_dependent_variable, "Apollo", "Earth" );
    boost::shared_ptr< PropagationTerminationSettings > terminationSettings =
            boost::make_shared< PropagationDependentVariableTerminationSettings >(
                terminationDependentVariable, 25.0E3, true );

    // Create propagation settings.
    boost::shared_ptr< PropagatorSettings< double > > propagatorSettings =
            boost::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState,
              terminationSettings, cowell );
    boost::shared_ptr< IntegratorSettings< > > integratorSettings =
            boost::make_shared< IntegratorSettings< > >
            ( rungeKutta4, simulationStartEpoch, fixedStepSize );

    boost::shared_ptr< EstimatableParameterSettings > perturbedParameterSettings;
    perturbedParameterSettings = (
                boost::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                    2, 0, 2, 0, "Earth", spherical_harmonics_cosine_coefficient_block ) );


    SingleArcDynamicsSimulator< > simulator( bodyMap, integratorSettings, propagatorSettings );

    std::pair< boost::shared_ptr< PodOutput< double > >, Eigen::VectorXd > estimationOutput =
            determinePostfitParameterInfluence(
                bodyMap, integratorSettings, propagatorSettings, perturbedParameterSettings,
                1.0, boost::assign::list_of( -earthC20 ), boost::assign::list_of( 0 ) );

    Eigen::VectorXd prefitResiduals = estimationOutput.first->firstIterationResiduals_;
    Eigen::VectorXd postfitResiduals = estimationOutput.first->residuals_;

    double prefitRms = linear_algebra::getVectorEntryRootMeanSquare(
                prefitResiduals );
    double postfitRms = linear_algebra::getVectorEntryRootMeanSquare(
                postfitResiduals );


    int observationVectorSize = prefitResiduals.rows( );
    Eigen::Vector3d initialPrefitDifference = prefitResiduals.segment( 0, 3 );
    Eigen::Vector3d initialPostfitDifference = postfitResiduals.segment( 0, 3 );

    Eigen::Vector3d finalPrefitDifference = prefitResiduals.segment( observationVectorSize - 3, 3 );
    Eigen::Vector3d finalPostfitDifference = postfitResiduals.segment( observationVectorSize - 3, 3 );

    std::cout<<"Parameter difference "<<estimationOutput.second.transpose( )<<std::endl;

    BOOST_CHECK_EQUAL( initialPrefitDifference.norm( ) < 1.0, true );

    BOOST_CHECK_EQUAL( finalPrefitDifference.norm( ) > 20.0E3, true );

    BOOST_CHECK_EQUAL( initialPostfitDifference.norm( ) > 275.0, true );
    BOOST_CHECK_EQUAL( initialPostfitDifference.norm( ) < 325.0, true );

    BOOST_CHECK_EQUAL( finalPostfitDifference.norm( ) > 275.0, true );
    BOOST_CHECK_EQUAL( finalPostfitDifference.norm( ) < 325.0, true );

    BOOST_CHECK_EQUAL( prefitRms / postfitRms > 40.0, true );

    BOOST_AUTO_TEST_SUITE_END( )

}

}

}



