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

#include <string>
#include <thread>

#include <limits>

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/simulation/estimation.h"

#include <tudat/astro/orbit_determination/acceleration_partials/numericalAccelerationPartial.h>
#include "tudat/simulation/estimation_setup/fitOrbitToEphemeris.h"


namespace tudat
{
namespace unit_tests
{
BOOST_AUTO_TEST_SUITE( test_estimate_arcwise_exponential_atmosphere_parameters )

// Using declarations.
using namespace tudat::observation_models;
using namespace tudat::orbit_determination;
using namespace tudat::estimatable_parameters;
using namespace tudat::interpolators;
using namespace tudat::numerical_integrators;
using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;
using namespace tudat::orbital_element_conversions;
using namespace tudat::ephemerides;
using namespace tudat::propagators;
using namespace tudat::basic_astrodynamics;
using namespace tudat::coordinate_conversions;
using namespace tudat::ground_stations;
using namespace tudat::observation_models;


//! Unit test to check if atmosphere parameters can be re-estimated after pertubed
BOOST_AUTO_TEST_CASE( test_EstimateArcwiseExponentialAtmosphereParameters )
{
    spice_interface::loadStandardSpiceKernels( );

    std::vector <double> arcStartTimes;
    std::vector <Eigen::Vector6d> arcInitialStates;

    // arc1 (Cassini T022, c/a-5min)
    arcStartTimes.push_back( 2.205723865899772644e+08 );
    arcInitialStates.push_back( Eigen::Vector6d{2.166092477145108860e+06, -7.925218812290439382e+05, 3.553645666049534921e+06, 3.304561684998871897e+03, 2.420278071381293103e+03, -4.244915762044645817e+03} );

    // arc 1 (Cassini T068 c/a -5min)
    arcStartTimes.push_back( 3.275976261739959717e+08 );
    arcInitialStates.push_back( Eigen::Vector6d{-3.091501373934821691e+06, 9.134845736499809427e+05, -2.880028589145255741e+06, 4.617009140001175183e+03, 3.555047142774454642e+03, -4.711669772164187862e+02} );

    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Titan" );
    bodyNames.push_back( "Saturn" );

    BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, "Titan" );

    bodySettings.addSettings( "Cassini" );

    std::map< double, Eigen::Vector6d > emptyMap;
    bodySettings.at( "Cassini" )->ephemerisSettings =
        std::make_shared< TabulatedEphemerisSettings >( emptyMap, "Titan", "ECLIPJ2000" );
    bodySettings.at( "Cassini" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );

    // Create bodies needed in simulation
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    double vehicleMass = 2.0E3;
    bodies.at( "Cassini" )->setConstantBodyMass( vehicleMass );

    Eigen::Vector3d forceCoefficients( 2.5, -0.1, 0.05 );
    std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
            std::make_shared< ConstantAerodynamicCoefficientSettings >(
                    20.0, forceCoefficients, aerodynamics::negative_aerodynamic_frame_coefficients );
    bodies.at( "Cassini" )->setAerodynamicCoefficientInterface(
            createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Cassini", bodies ) );

    bodies.at( "Titan" )->setAtmosphereModel(
        createAtmosphereModel( std::make_shared< ExponentialAtmosphereSettings >( 7.58e04, 175, 3.0553447e-04 ), "Titan" ));

    // Set accelerations between bodies that are to be taken into account.
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfSpacecraft;
    accelerationsOfSpacecraft[ "Titan" ].push_back( std::make_shared< AccelerationSettings >( aerodynamic ) );
    accelerationsOfSpacecraft[ "Titan" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationsOfSpacecraft[ "Saturn" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );

    accelerationMap[ "Cassini" ] = accelerationsOfSpacecraft;

    // Set bodies for which initial state is to be estimated and integrated.
    std::vector< std::string > bodiesToEstimate = { "Cassini" };
    std::vector< std::string > centralBodies = { "Titan" };

    AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodies, accelerationMap, bodiesToEstimate, centralBodies );

    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > > additionalParameterNames;

    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;
    dependentVariablesList.push_back( singleAccelerationDependentVariable( basic_astrodynamics::aerodynamic, "Cassini", "Titan" ) );

    std::vector< double > integrationArcStartTimes = arcStartTimes;
    std::vector< double > integrationArcEndTimes;


    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsList;
    for( unsigned int j = 0; j < integrationArcStartTimes.size( ); j++ )
    {
        integrationArcEndTimes.push_back( integrationArcStartTimes.at( j ) + 10*60);
        Eigen::Matrix< double, Eigen::Dynamic, 1 > currentInitialState = arcInitialStates.at( j );

        propagatorSettingsList.push_back( std::make_shared< TranslationalStatePropagatorSettings< double > >(
                centralBodies,
                accelerationModelMap,
                bodiesToEstimate,
                currentInitialState,
                integrationArcStartTimes.at( j ),
                numerical_integrators::rungeKuttaFixedStepSettings( 1.0, numerical_integrators::rungeKuttaFehlberg78 ),
                propagationTimeTerminationSettings( integrationArcEndTimes.at( j ) ),
                cowell,
                dependentVariablesList ) );
    }

    std::shared_ptr< MultiArcPropagatorSettings< double > > propagatorSettings =
            std::make_shared< MultiArcPropagatorSettings< double > >(
                    propagatorSettingsList, false, std::make_shared< MultiArcPropagatorProcessingSettings >( false, true ) );


    std::cout << "************************************* RUNNING TEST  *************************" << std::endl;

    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
            getInitialMultiArcParameterSettings< double >( propagatorSettings, bodies, integrationArcStartTimes );

    additionalParameterNames.push_back( estimatable_parameters::arcwiseExponentialAtmosphereBaseDensity( "Titan" , arcStartTimes ) );
    additionalParameterNames.push_back( estimatable_parameters::arcwiseExponentialAtmosphereScaleHeight( "Titan", arcStartTimes ) );

    parameterNames.insert( parameterNames.end( ), additionalParameterNames.begin( ), additionalParameterNames.end( ) );

    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double, double >( parameterNames, bodies, propagatorSettings );
    printEstimatableParameterEntries( parametersToEstimate );

    std::vector< double > observationTimes;
    double dataPointInterval = 1.;
    for( int k = 0; k < integrationArcStartTimes.size( ); k++ )
    {
        double currentTime = integrationArcStartTimes.at( k ) + 5.0;
        while( currentTime < integrationArcEndTimes.at( k ) - 5.0 )

        {
            observationTimes.push_back( currentTime );
            currentTime += static_cast< double >( dataPointInterval );
        }
    }

    // Create simulation object (but do not propagate dynamics).
    MultiArcDynamicsSimulator<> dynamicsSimulator( bodies, propagatorSettings, true );

    std::pair< std::vector< std::shared_ptr< observation_models::ObservationModelSettings > >,
               std::shared_ptr< observation_models::ObservationCollection< double > > >
            observationCollectionAndModelSettings =
                    simulatePseudoObservations( bodies, bodiesToEstimate, centralBodies, observationTimes );
    std::shared_ptr< observation_models::ObservationCollection< double > > observationCollection =
            observationCollectionAndModelSettings.second;

    std::vector< std::shared_ptr< observation_models::ObservationModelSettings > > observationModelSettingsList =
            observationCollectionAndModelSettings.first;

    Eigen::VectorXd truthParameters = parametersToEstimate->getFullParameterValues< double >( );

    // resetParameterValues
    Eigen::VectorXd perturbation1( 4 );
    perturbation1 << -3e-05, 6e-05, 1e03, -1.5e03;

    Eigen::VectorXd paramVector = parametersToEstimate->getFullParameterValues< double >( );
    paramVector.tail( perturbation1.size( ) ) += perturbation1;
    parametersToEstimate->resetParameterValues( paramVector );


    OrbitDeterminationManager< double > orbitDeterminationManager =
            OrbitDeterminationManager< double >( bodies, parametersToEstimate, observationModelSettingsList, propagatorSettings );

    std::shared_ptr< EstimationInput< double > > estimationInput =
            std::make_shared< EstimationInput< double > >( observationCollection );
    estimationInput->defineEstimationSettings( 0, 1, 0, 1, 1, 1 );

    std::shared_ptr< EstimationOutput<> > estimationOutput = orbitDeterminationManager.estimateParameters( estimationInput );

    // Check if parameters are correctly estimated
    Eigen::VectorXd estimatedParametervalues = estimationOutput->parameterEstimate_;

    // Initial states (2 arcs).
    for( unsigned int l = 0; l < 3; l++ )
    {
        BOOST_CHECK_SMALL( std::fabs( truthParameters( l ) - estimatedParametervalues( l ) ), 0.1 );
        BOOST_CHECK_SMALL( std::fabs( truthParameters( l + 6 ) - estimatedParametervalues( l + 6 ) ), 0.1 );

        BOOST_CHECK_SMALL( std::fabs( truthParameters( l + 3 ) - estimatedParametervalues( l + 3 ) ), 1.0E-6 );
        BOOST_CHECK_SMALL( std::fabs( truthParameters( l + 9 ) - estimatedParametervalues( l + 9 ) ), 1.0E-6 );
    }

    // NOTE: Estimation of atmospheric parameters is may not be particularly stable, because highly correlated...
    // global aerodynamic parameters.
    BOOST_CHECK_SMALL( std::fabs( truthParameters( 12 ) - estimatedParametervalues( 12 ) ) / std::fabs(truthParameters( 12 )), 1.0e-5 );
    BOOST_CHECK_SMALL( std::fabs( truthParameters( 13 ) - estimatedParametervalues( 13 ) ) / std::fabs(truthParameters( 13 )), 1.0e-5 );
    BOOST_CHECK_SMALL( std::fabs( truthParameters( 14 ) - estimatedParametervalues( 14 ) ) / std::fabs(truthParameters( 14 )), 1.0e-5 );
    BOOST_CHECK_SMALL( std::fabs( truthParameters( 15 ) - estimatedParametervalues( 15 ) ) / std::fabs(truthParameters( 15 )), 1.0e-5 );

}


BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
