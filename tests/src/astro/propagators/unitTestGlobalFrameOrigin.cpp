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

BOOST_AUTO_TEST_SUITE( test_global_frame_origin )

//! Test to ensure that a point-mass acceleration on a body produces a Kepler orbit (to within
//! numerical error bounds).
template< typename TimeType, typename StateScalarType >
Eigen::Matrix< StateScalarType, 6, 1 > testGlobalFrameOrigin(
        const std::string& globalFrameOrigin, const std::string& moonEphemerisOrigin  )
{
    std::cout << "Testing with origin: ************************************** " << globalFrameOrigin << std::endl;
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
    double finalEphemerisTime = initialEphemerisTime + 56.0 * 86400.0;
    double maximumTimeStep = 3600.0;
    double buffer = 5.0 * maximumTimeStep;

    // Create bodies needed in simulation
    BodyListSettings bodySettings =
            getDefaultBodySettings( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer,
                                    globalFrameOrigin, "ECLIPJ2000" );

    if( std::is_same< long double, StateScalarType >::value )
    {
        std::dynamic_pointer_cast< InterpolatedSpiceEphemerisSettings >(
                    bodySettings.at( "Moon" )->ephemerisSettings )->setUseLongDoubleStates( 1 );
    }

    // Change ephemeris origins to test full functionality
    std::dynamic_pointer_cast< InterpolatedSpiceEphemerisSettings >( bodySettings.at( "Moon" )->ephemerisSettings )->
            resetFrameOrigin( moonEphemerisOrigin );
    std::dynamic_pointer_cast< InterpolatedSpiceEphemerisSettings >( bodySettings.at( "Earth" )->ephemerisSettings )->
            resetFrameOrigin( "Sun" );
    std::dynamic_pointer_cast< InterpolatedSpiceEphemerisSettings >( bodySettings.at( "Mars" )->ephemerisSettings )->
            resetFrameOrigin( "Sun" );
    std::dynamic_pointer_cast< InterpolatedSpiceEphemerisSettings >( bodySettings.at( "Venus" )->ephemerisSettings )->
            resetFrameOrigin( "SSB" );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // Set accelerations between bodies that are to be taken into account.
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfMoon;
    accelerationsOfMoon[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationsOfMoon[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationsOfMoon[ "Mars" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationsOfMoon[ "Venus" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationMap[ "Moon" ] = accelerationsOfMoon;

    // Propagate the moon only
    std::vector< std::string > bodiesToIntegrate;
    bodiesToIntegrate.push_back( "Moon" );
    std::vector< std::string > centralBodies;
    centralBodies.push_back( "Earth" );

    // Define settings for numerical integrator.
    std::shared_ptr< IntegratorSettings< TimeType > > integratorSettings =
            std::make_shared< IntegratorSettings< TimeType > >
            ( rungeKutta4, initialEphemerisTime, 300.0 );

    // Create acceleration models and propagation settings.
    Eigen::Matrix< StateScalarType, 6, 1  > systemInitialState =
                spice_interface::getBodyCartesianStateAtEpoch(
                  bodiesToIntegrate[ 0 ], centralBodies[ 0 ], "ECLIPJ2000", "NONE", initialEphemerisTime ).
                template cast< StateScalarType >( );
    AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodies, accelerationMap, bodiesToIntegrate, centralBodies );
    std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< StateScalarType > >
            ( centralBodies, accelerationModelMap, bodiesToIntegrate, systemInitialState, finalEphemerisTime );

    // Create dynamics simulation object.
    SingleArcDynamicsSimulator< StateScalarType, TimeType > dynamicsSimulator(
                bodies, integratorSettings, propagatorSettings, true, false, true );

    return dynamicsSimulator.getEquationsOfMotionNumericalSolution( ).rbegin( )->second;
}

BOOST_AUTO_TEST_CASE( testCowellPropagatorKeplerCompare )
{
    std::vector< std::string > origins;

    origins.push_back( "SSB" );
    origins.push_back( "Earth" );
    origins.push_back( "Mars" );
    origins.push_back( "Venus" );
    origins.push_back( "Sun" );

    // Get final state with Earth ephemeris origin and SSB global frame origin
    Eigen::Vector6d benchmarkFinalState = testGlobalFrameOrigin< double, double >( "SSB", "Earth" );

    // Iterate over all combinations of ephemeris/global frame origins
    for( unsigned int i = 0; i < origins.size( ); i++ )
    {
        for( unsigned int j = 0; j < origins.size( ); j++ )
        {
            Eigen::Vector6d currentFinalState = testGlobalFrameOrigin< double, double >(
                        origins.at( i ), origins.at( j ) );
            for( unsigned int k = 0; k < 3 ; k++ )
            {
                BOOST_CHECK_SMALL( std::fabs( benchmarkFinalState( k ) - currentFinalState( k ) ), 1.0E-4 );
                BOOST_CHECK_SMALL( std::fabs( benchmarkFinalState( k + 3 ) - currentFinalState( k + 3 ) ), 1.0E-9 );
            }


        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )


}

}
