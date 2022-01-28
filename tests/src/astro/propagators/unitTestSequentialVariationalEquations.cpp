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
#include <boost/make_shared.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/math/basic/linearAlgebra.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/astro/ephemerides/tabulatedEphemeris.h"

#include "tudat/io/basicInputOutput.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/astro/basic_astro/accelerationModel.h"

#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/createBodies.h"
#include "tudat/simulation/estimation_setup/createNumericalSimulator.h"
#include "tudat/simulation/estimation_setup/createEstimatableParameters.h"
#include "tudat/simulation/estimation_setup/variationalEquationsSolver.h"

namespace tudat
{

namespace unit_tests
{

//Using declarations.
using namespace tudat::interpolators;
using namespace tudat::numerical_integrators;
using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;
using namespace tudat::basic_astrodynamics;
using namespace tudat::estimatable_parameters;
using namespace tudat::orbital_element_conversions;
using namespace tudat::ephemerides;
using namespace tudat::propagators;


BOOST_AUTO_TEST_SUITE( test_sequential_variational_equation_integration )


std::pair< std::shared_ptr< CombinedStateTransitionAndSensitivityMatrixInterface >, std::shared_ptr< Ephemeris > >
integrateEquations( const bool performIntegrationsSequentially )
{
    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Moon" );

    // Specify initial time
    double initialEphemerisTime = 1.0E7;
    double finalEphemerisTime = initialEphemerisTime + 14.0 * 86400.0;
    double maximumTimeStep = 600.0;

    double numberOfTimeStepBuffer = 6.0;
    double buffer = numberOfTimeStepBuffer * maximumTimeStep;

    // Create bodies needed in simulation
    BodyListSettings bodySettings =
            getDefaultBodySettings( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    SystemOfBodies bodies =
            createSystemOfBodies( bodySettings );
    std::shared_ptr< Body > lageos = std::make_shared< Body >( );
    bodies.addBody( lageos, "LAGEOS" );

    // Create  body initial state
    Eigen::Vector6d lageosKeplerianElements;
    lageosKeplerianElements[ semiMajorAxisIndex ] = 8000.0E3;
    lageosKeplerianElements[ eccentricityIndex ] = 0.0044;
    lageosKeplerianElements[ inclinationIndex ] = 109.89 * mathematical_constants::PI / 180.0;
    lageosKeplerianElements[ argumentOfPeriapsisIndex ] = 259.35 * mathematical_constants::PI / 180.0;
    lageosKeplerianElements[ longitudeOfAscendingNodeIndex ] = 31.56 * mathematical_constants::PI / 180.0;
    lageosKeplerianElements[ trueAnomalyIndex ] = 1.0;
    Eigen::Vector6d lageosState = convertKeplerianToCartesianElements(
                lageosKeplerianElements, getBodyGravitationalParameter("Earth" ) );

    // Set accelerations between bodies that are to be taken into account.
    SelectedAccelerationMap accelerationMap;

    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfLageos;
    //accelerationsOfLageos[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    //accelerationsOfLageos[ "Earth" ].push_back( std::make_shared< RelativisticCorrectionSettings >( ) );
    //accelerationsOfLageos[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 8, 8 ) );
    accelerationsOfLageos[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationMap[ "LAGEOS" ] = accelerationsOfLageos;

    // Set bodies for which initial state is to be estimated and integrated.
    std::vector< std::string > bodiesToIntegrate;
    bodiesToIntegrate.push_back( "LAGEOS" );
    unsigned int numberOfNumericalBodies = bodiesToIntegrate.size( );

    std::vector< std::string > centralBodies;
    std::map< std::string, std::string > centralBodyMap;

    centralBodies.resize( numberOfNumericalBodies );
    for( unsigned int i = 0; i < numberOfNumericalBodies; i++ )
    {
        centralBodies[ i ] = "Earth";
        centralBodyMap[ bodiesToIntegrate[ i ] ] = centralBodies[ i ];
    }
    AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodies, accelerationMap, centralBodyMap );

    // Define propagator settings.
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToIntegrate, lageosState, finalEphemerisTime );

    // Set parameters that are to be included.
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
            getInitialStateParameterSettings< double >( propagatorSettings, bodies );
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >
                              ( "Earth", gravitational_parameter ) );
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >
                              ( "Moon", gravitational_parameter ) );
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double >( parameterNames, bodies, propagatorSettings );

    // Define integrator settings.
    std::shared_ptr< IntegratorSettings< > > matrixTypeIntegratorSettings =
            std::make_shared< RungeKuttaVariableStepSizeSettings< > >
            ( initialEphemerisTime, 10.0,
              RungeKuttaCoefficients::rungeKuttaFehlberg45, 0.01, 10.0, 1.0E-6, 1.0E-6 );



    // Perform requested propagation
    std::shared_ptr< SingleArcVariationalEquationsSolver< double, double> > variationalEquationSolver;
    if( !performIntegrationsSequentially )
    {
        // Propagate
        variationalEquationSolver = std::make_shared< SingleArcVariationalEquationsSolver< double, double> >(
                    bodies, matrixTypeIntegratorSettings,
                    propagatorSettings, parametersToEstimate );
    }
    else
    {

        // Define integrator settings for vector type.
        std::shared_ptr< IntegratorSettings< > > vectorTypeIntegratorSettings =
                std::make_shared< RungeKuttaVariableStepSizeSettings< > >
                ( initialEphemerisTime, 10.0, RungeKuttaCoefficients::rungeKuttaFehlberg45, 0.01, 10.0, 1.0E-6, 1.0E-6 );

        // Propagate
        variationalEquationSolver = std::make_shared< SingleArcVariationalEquationsSolver< double, double > >(
                    bodies, vectorTypeIntegratorSettings,
                    propagatorSettings, parametersToEstimate, 0,
                    matrixTypeIntegratorSettings );
    }

    return std::make_pair( variationalEquationSolver->getStateTransitionMatrixInterface( ),
                           bodies.at( "LAGEOS" )->getEphemeris( ) );
}

//! Test whether concurrent and sequential propagation of variational equations gives same results.
BOOST_AUTO_TEST_CASE( testSequentialVariationalEquationIntegration )
{
    // Propagate concurrently.
    std::pair< std::shared_ptr< CombinedStateTransitionAndSensitivityMatrixInterface >, std::shared_ptr< Ephemeris > >
            concurrentResult = integrateEquations( 0 );

    // Propagate sequentially.
    std::pair< std::shared_ptr< CombinedStateTransitionAndSensitivityMatrixInterface >, std::shared_ptr< Ephemeris > >
            sequentialResult = integrateEquations( 1 );

    // Test variational equations solution.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                concurrentResult.first->getCombinedStateTransitionAndSensitivityMatrix( 1.0E7 + 14.0 * 80000.0 ),
                sequentialResult.first->getCombinedStateTransitionAndSensitivityMatrix( 1.0E7 + 14.0 * 80000.0 ), 2.0E-6 );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                concurrentResult.first->getCombinedStateTransitionAndSensitivityMatrix( 1.0E7 + 14.0 * 80000.0 ),
                sequentialResult.first->getFullCombinedStateTransitionAndSensitivityMatrix( 1.0E7 + 14.0 * 80000.0 ), 2.0E-6 );

    // Test dynamics solution.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                concurrentResult.second->getCartesianState( 1.0E7 + 14.0 * 80000.0 ),
                sequentialResult.second->getCartesianState( 1.0E7 + 14.0 * 80000.0 ),
                std::numeric_limits< double >::epsilon( ) );
}

BOOST_AUTO_TEST_SUITE_END( )

}

}
