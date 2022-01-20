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

#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>

#include "tudat/math/basic/linearAlgebra.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/orbitalElementConversions.h"

#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/math/integrators/rungeKuttaCoefficients.h"
#include "tudat/math/interpolators/lagrangeInterpolator.h"
#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/astro/basic_astro/keplerPropagator.h"
#include "tudat/io/basicInputOutput.h"

#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/astro/propagators/nBodyCowellStateDerivative.h"
#include "tudat/simulation/propagation_setup/dynamicsSimulator.h"
#include "tudat/math/integrators/createNumericalIntegrator.h"
#include "tudat/simulation/environment_setup/createBodies.h"
#include "tudat/simulation/estimation_setup/createNumericalSimulator.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"

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
using namespace tudat;

BOOST_AUTO_TEST_SUITE( test_forwards_backwards_rpopagation )

template< typename TimeType >
std::shared_ptr< IntegratorSettings< TimeType > > getIntegrationSettings(
        const int integratorCase, const int initialTime, const bool propagateForwards )
{
    double initialTimeMultiplier = ( propagateForwards ? 1.0 : -1.0 );
    std::shared_ptr< IntegratorSettings< TimeType > > integratorSettings;
    if( integratorCase == 0 )
    {
        integratorSettings = std::make_shared< IntegratorSettings< TimeType > >
                ( rungeKutta4, initialTime, initialTimeMultiplier * 300.0 );
    }
    else if( integratorCase < 5 )
    {
        RungeKuttaCoefficients::CoefficientSets coefficientSet = RungeKuttaCoefficients::undefinedCoefficientSet;
        if( integratorCase == 1 )
        {
            coefficientSet = RungeKuttaCoefficients::rungeKuttaFehlberg45;
        }
        else if( integratorCase == 2 )
        {
            coefficientSet = RungeKuttaCoefficients::rungeKuttaFehlberg56;

        }
        else if( integratorCase == 3 )
        {
            coefficientSet = RungeKuttaCoefficients::rungeKuttaFehlberg78;

        }
        else if( integratorCase == 4 )
        {
            coefficientSet = RungeKuttaCoefficients::rungeKutta87DormandPrince;

        }
        integratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettings< TimeType > >
                ( initialTime, initialTimeMultiplier * 300.0, coefficientSet, 1.0E-3, 3600.0 );
    }
    return integratorSettings;
}

//! Test to ensure that forward and backward in time integration are properly and consistently performed using different
//! integrators. The state of the Moon is numerically propagated forward, and the back to the original time. The test checks
//! if the state at the midpoint is sufficiently close for the foprward and backward intgrations.
template< typename TimeType = double, typename StateScalarType = double >
Eigen::Matrix< StateScalarType, 6, 1 > propagateForwardBackwards( const int integratorCase )
{
    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Define bodies in simulation.
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Moon" );
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Mars" );
    bodyNames.push_back( "Venus" );

    // Specify initial time
    double initialEphemerisTime = 1.0E7;
    double finalEphemerisTime = initialEphemerisTime + 86400.0;
    double buffer = 3600.0;

    // Create bodies needed in simulation
    BodyListSettings bodySettings =
            getDefaultBodySettings( bodyNames, initialEphemerisTime - 2.0 * buffer, finalEphemerisTime + 2.0 * buffer,
                                    "Earth", "ECLIPJ2000" );
    std::dynamic_pointer_cast< InterpolatedSpiceEphemerisSettings >( bodySettings.at( "Moon" )->ephemerisSettings )->
            resetFrameOrigin( "Earth" );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // Set accelerations between bodies that are to be taken into account.
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfMoon;
    accelerationsOfMoon[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationsOfMoon[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationMap[ "Moon" ] = accelerationsOfMoon;

    // Propagate the moon only
    std::vector< std::string > bodiesToIntegrate;
    bodiesToIntegrate.push_back( "Moon" );
    std::vector< std::string > centralBodies;
    centralBodies.push_back( "Earth" );

    // Define settings for numerical integrator.
    std::shared_ptr< IntegratorSettings< TimeType > > integratorSettings = getIntegrationSettings< TimeType >(
                integratorCase, initialEphemerisTime, true );

    // Propagate forwards
    {


        // Create acceleration models and propagation settings.
        Eigen::Matrix< StateScalarType, 6, 1  > systemInitialState =
                spice_interface::getBodyCartesianStateAtEpoch(
                    bodiesToIntegrate[ 0 ], centralBodies[ 0 ], "ECLIPJ2000", "NONE", initialEphemerisTime ).
                template cast< StateScalarType >( );
        AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodies, accelerationMap, bodiesToIntegrate, centralBodies );
        std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > > propagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< StateScalarType > >
                ( centralBodies, accelerationModelMap, bodiesToIntegrate, systemInitialState, finalEphemerisTime + buffer );

        // Create dynamics simulation object.
        SingleArcDynamicsSimulator< StateScalarType, TimeType > dynamicsSimulator(
                    bodies, integratorSettings, propagatorSettings, true, true, true );
    }

    double testTime = initialEphemerisTime + ( finalEphemerisTime - initialEphemerisTime ) / 2.0;
    Eigen::Vector6d forwardState = bodies.at( "Moon" )->getEphemeris( )->getCartesianState( testTime );

    // Re-define settings for numerical integrator.
    integratorSettings = getIntegrationSettings< TimeType >( integratorCase, finalEphemerisTime, false );

    // Propagate backwards
    {
        // Create acceleration models and propagation settings.
        Eigen::Matrix< StateScalarType, 6, 1  > systemInitialState =
                spice_interface::getBodyCartesianStateAtEpoch(
                    bodiesToIntegrate[ 0 ], centralBodies[ 0 ], "ECLIPJ2000", "NONE", finalEphemerisTime ).
                template cast< StateScalarType >( );
        AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodies, accelerationMap, bodiesToIntegrate, centralBodies );
        std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > > propagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< StateScalarType > >
                ( centralBodies, accelerationModelMap, bodiesToIntegrate, systemInitialState, initialEphemerisTime - buffer );

        // Create dynamics simulation object.
        SingleArcDynamicsSimulator< StateScalarType, TimeType > dynamicsSimulator(
                    bodies, integratorSettings, propagatorSettings, true, true, true );
    }

    Eigen::Vector6d backwardState = bodies.at( "Moon" )->getEphemeris( )->getCartesianState( testTime );

    return forwardState - backwardState;
}

BOOST_AUTO_TEST_CASE( testCowellPropagatorKeplerCompare )
{
    for( unsigned int j = 0; j < 5; j++ )
    {
        Eigen::Vector6d stateDifference = propagateForwardBackwards< double, double >( j );
        for( int i = 0; i < 3; i++ )
        {
            BOOST_CHECK_SMALL( std::fabs( stateDifference( i ) ), 0.1 );
            BOOST_CHECK_SMALL( std::fabs( stateDifference( i + 3 ) ), 1.0E-4 );
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )


}

}
