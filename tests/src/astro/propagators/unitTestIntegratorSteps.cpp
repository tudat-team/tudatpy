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

BOOST_AUTO_TEST_SUITE( test_cowell_propagator )

BOOST_AUTO_TEST_CASE( testCowellPropagatorKeplerCompare )
{
    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Define bodies in simulation.
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Moon" );

    // Specify initial time
    double initialEphemerisTime = 1.0E7;
    double finalEphemerisTime = initialEphemerisTime + 2.0 * 3600.0;

    // Create bodies needed in simulation
    BodyListSettings bodySettings =
        getDefaultBodySettings(
            bodyNames, "Earth", "ECLIPJ2000" );

    SystemOfBodies bodies = createSystemOfBodies< double, double >( bodySettings );
    bodies.createEmptyBody( "Vehicle" );

    // Set accelerations between bodies that are to be taken into account.
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationsOfVehicle[ "Moon" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;

    // Create acceleration models and propagation settings.
    AccelerationMap accelerationModelMap = createAccelerationModelsMap(
        bodies, accelerationMap,
        std::vector< std::string >( { "Vehicle" } ),
        std::vector< std::string >( { "Earth" } ) );

    // Set Keplerian elements for Vehicle.
    Eigen::Vector6d vehicleInitialStateInKeplerianElements;
    vehicleInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7500.0E3;
    vehicleInitialStateInKeplerianElements( eccentricityIndex ) = 0.1;
    vehicleInitialStateInKeplerianElements( inclinationIndex ) = 0.0;
    vehicleInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = 0.0;
    vehicleInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = 0.0;
    vehicleInitialStateInKeplerianElements( trueAnomalyIndex ) = 0.0;

    // Convert Vehicle state from Keplerian elements to Cartesian elements.
    double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );


    for( unsigned int integratorCase = 0; integratorCase < 8; integratorCase++ )
    {
        for ( unsigned int orbitCase = 0; orbitCase < 3; orbitCase++ )
        {
            // Define settings for numerical integrator.
            std::shared_ptr< IntegratorSettings< double > > integratorSettings;
            if( integratorCase == 0 )
            {
                integratorSettings = rungeKuttaFixedStepSettings(
                    60.0, CoefficientSets::rungeKutta4Classic );
            }
            else if( integratorCase == 1 )
            {
                integratorSettings = rungeKuttaVariableStepSettingsScalarTolerances(
                        60.0, rungeKuttaFehlberg78, 1.0E-3, 3600.0, 1.0E-6, 1.0E-6 );
            }
            else if( integratorCase == 2 )
            {
                integratorSettings = bulirschStoerFixedStepIntegratorSettings(
                    600.0, bulirsch_stoer_sequence, 6 );
            }
            else if( integratorCase == 3 )
            {
                integratorSettings = bulirschStoerIntegratorSettings< double >(
                    600.0, bulirsch_stoer_sequence, 8, 1.0E-3, 7200.0, 1.0E-6, 1.0E-6 );
            }
            else if( integratorCase == 4 )
            {
                integratorSettings = adamsBashforthMoultonSettingsFixedStep(
                    15.0, 1.0E-6, 1.0E-6 );
            }
            else if( integratorCase == 5 )
            {
                 integratorSettings = adamsBashforthMoultonSettings(
                    30.0, 1.0E-3, 1200.0, 1.0E-6, 1.0E-6 );
            }
            else if( integratorCase == 6 )
            {
                integratorSettings = adamsBashforthMoultonSettingsFixedStepFixedOrder(
                    10.0, 6 );
            }
            else if( integratorCase == 7 )
            {
                integratorSettings = adamsBashforthMoultonSettingsFixedOrder(
                    30.0, 1.0E-3, 1200.0, 1.0E-6, 1.0E-6 );
            }

            vehicleInitialStateInKeplerianElements( 1 ) = static_cast< double >( orbitCase ) * 0.45;
            Eigen::VectorXd systemInitialState = convertKeplerianToCartesianElements(
                vehicleInitialStateInKeplerianElements,
                earthGravitationalParameter );

            std::shared_ptr<TranslationalStatePropagatorSettings<double, double> > propagatorSettings =
                translationalStatePropagatorSettings(
                    std::vector< std::string >( { "Earth" } ), accelerationModelMap, std::vector< std::string >( { "Vehicle" } ),
                        systemInitialState, initialEphemerisTime, integratorSettings,
                    propagationTimeTerminationSettings( finalEphemerisTime ) );

            // Create dynamics simulation object.
            auto dynamicsSimulator = std::dynamic_pointer_cast< SingleArcDynamicsSimulator< > >(
                createDynamicsSimulator<double, double>( bodies, propagatorSettings ) );
            std::map< double, Eigen::VectorXd > stateHistory =
                dynamicsSimulator->getSingleArcPropagationResults( )->getEquationsOfMotionNumericalSolution( );

            auto firstIterator = stateHistory.begin( );
            auto secondIterator = stateHistory.begin( );
            secondIterator++;

            double firstStepSize = secondIterator->first - firstIterator->first;
            bool allStepSizesEqual = true;
            while( secondIterator != stateHistory.end( ) )
            {
                double currentStepSize = secondIterator->first - firstIterator->first;
                if( std::fabs( currentStepSize - firstStepSize ) / firstStepSize > 1.0E-14 )
                {
                    allStepSizesEqual = false;
                }
                firstIterator++;
                secondIterator++;
            }

            if( integratorCase % 2 == 0 )
            {
                BOOST_CHECK_EQUAL( allStepSizesEqual, true );
                if( integratorCase == 0 )
                {
                    BOOST_CHECK_CLOSE_FRACTION( firstStepSize, 60.0, 10.0 * std::numeric_limits< double >::epsilon( ) );
                }
                else if( integratorCase == 2 )
                {
                    BOOST_CHECK_CLOSE_FRACTION( firstStepSize, 600.0, 10.0 * std::numeric_limits< double >::epsilon( ) );
                }
                else if( integratorCase == 4 )
                {
                    BOOST_CHECK_CLOSE_FRACTION( firstStepSize, 15.0, 10.0 * std::numeric_limits< double >::epsilon( ) );
                }
                else if( integratorCase == 6 )
                {
                    BOOST_CHECK_CLOSE_FRACTION( firstStepSize, 10.0, 10.0 * std::numeric_limits< double >::epsilon( ) );
                }
            }
            else
            {
                BOOST_CHECK_EQUAL( allStepSizesEqual, false );
            }
        }

    }
}

BOOST_AUTO_TEST_SUITE_END( )


}

}
