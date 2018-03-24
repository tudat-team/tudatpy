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

#include <string>
#include <thread>
#include <limits>

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>

#include "Tudat/Basics/testMacros.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Astrodynamics/Ephemerides/tabulatedEphemeris.h"

#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"

#include "Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createBodies.h"
#include "Tudat/SimulationSetup/PropagationSetup/createNumericalSimulator.h"
#include "Tudat/SimulationSetup/EstimationSetup/createEstimatableParameters.h"
#include "Tudat/SimulationSetup/PropagationSetup/variationalEquationsSolver.h"

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


std::pair< boost::shared_ptr< CombinedStateTransitionAndSensitivityMatrixInterface >, boost::shared_ptr< Ephemeris > >
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
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    NamedBodyMap bodyMap =
            createBodies( bodySettings );
    boost::shared_ptr< Body > lageos = boost::make_shared< Body >( );
    bodyMap[ "LAGEOS" ] = lageos;

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

    lageos->setEphemeris( boost::make_shared< TabulatedCartesianEphemeris< double, double > >(
                              boost::shared_ptr< interpolators::OneDimensionalInterpolator<
                              double, Eigen::Vector6d > >( ), "Earth" ) );
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    // Set accelerations between bodies that are to be taken into account.
    SelectedAccelerationMap accelerationMap;

    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfLageos;
    //accelerationsOfLageos[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
    //accelerationsOfLageos[ "Earth" ].push_back( boost::make_shared< RelativisticCorrectionSettings >( ) );
    //accelerationsOfLageos[ "Earth" ].push_back( boost::make_shared< SphericalHarmonicAccelerationSettings >( 8, 8 ) );
    accelerationsOfLageos[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
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
                bodyMap, accelerationMap, centralBodyMap );

    // Set parameters that are to be included.
    std::vector< boost::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( boost::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
                                  "LAGEOS", lageosState, "Earth" ) );
    parameterNames.push_back( boost::make_shared< EstimatableParameterSettings >
                              ( "Earth", gravitational_parameter ) );
    parameterNames.push_back( boost::make_shared< EstimatableParameterSettings >
                              ( "Moon", gravitational_parameter ) );
    boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate( parameterNames, bodyMap, accelerationModelMap );

    // Define integrator settings.
    boost::shared_ptr< IntegratorSettings< > > integratorSettings =
            boost::make_shared< RungeKuttaVariableStepSizeSettings< > >
            ( rungeKuttaVariableStepSize, initialEphemerisTime, 10.0,
              RungeKuttaCoefficients::rungeKuttaFehlberg45, 0.01, 10.0, 1.0E-6, 1.0E-6 );

    // Define propagator settings.
    boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            boost::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToIntegrate, lageosState, finalEphemerisTime );

    // Perform requested propagation
    boost::shared_ptr< SingleArcVariationalEquationsSolver< double, double> > variationalEquationSolver;
    if( !performIntegrationsSequentially )
    {
        variationalEquationSolver = boost::make_shared< SingleArcVariationalEquationsSolver< double, double> >(
                    bodyMap, integratorSettings,
                    propagatorSettings, parametersToEstimate );
    }
    else
    {
        variationalEquationSolver = boost::make_shared< SingleArcVariationalEquationsSolver< double, double> >(
                    bodyMap, integratorSettings,
                    propagatorSettings, parametersToEstimate, 0,
                    integratorSettings );
    }

    return std::make_pair( variationalEquationSolver->getStateTransitionMatrixInterface( ),
                           bodyMap[ "LAGEOS" ]->getEphemeris( ) );
}

//! Test whether concurrent and sequential propagation of variational equations gives same results.
BOOST_AUTO_TEST_CASE( testSequentialVariationalEquationIntegration )
{
    // Propagate concurrently.
    std::pair< boost::shared_ptr< CombinedStateTransitionAndSensitivityMatrixInterface >, boost::shared_ptr< Ephemeris > >
            concurrentResult = integrateEquations( 0 );

    // Propagate sequentially.
    std::pair< boost::shared_ptr< CombinedStateTransitionAndSensitivityMatrixInterface >, boost::shared_ptr< Ephemeris > >
            sequentialResult = integrateEquations( 1 );

    // Test variational equations solution.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                concurrentResult.first->getCombinedStateTransitionAndSensitivityMatrix( 1.0E7 + 14.0 * 80000.0 ),
                sequentialResult.first->getCombinedStateTransitionAndSensitivityMatrix( 1.0E7 + 14.0 * 80000.0 ), 2.0E-6 );

    // Test dynamics solution.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                concurrentResult.second->getCartesianState( 1.0E7 + 14.0 * 80000.0 ),
                sequentialResult.second->getCartesianState( 1.0E7 + 14.0 * 80000.0 ),
                std::numeric_limits< double >::epsilon( ) );

}

BOOST_AUTO_TEST_SUITE_END( )

}

}



