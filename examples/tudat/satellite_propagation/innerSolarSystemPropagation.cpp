/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <tudat/simulation/simulation.h>

#include "tudat/io/applicationOutput.h"

//! Simulate the dynamics of the main bodies in the inner solar system
int main( )
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    using namespace tudat;
    using namespace tudat::ephemerides;
    using namespace tudat::interpolators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::spice_interface;
    using namespace tudat::simulation_setup;
    using namespace tudat::basic_astrodynamics;
    using namespace tudat::basic_mathematics;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat:: propagators;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ENVIRONMENT            //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    //Load spice kernels.
    std::string kernelsPath = input_output::getSpiceKernelPath( );
    spice_interface::loadStandardSpiceKernels( );

    // Define bodies in simulation.
    unsigned int totalNumberOfBodies = 6;
    std::vector< std::string > bodyNames;
    bodyNames.resize( totalNumberOfBodies );
    bodyNames[ 0 ] = "Moon";
    bodyNames[ 1 ] = "Earth";
    bodyNames[ 2 ] = "Mars";
    bodyNames[ 3 ] = "Venus";
    bodyNames[ 4 ] = "Mercury";
    bodyNames[ 5 ] = "Sun";

    // Create bodies needed in simulation
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodyNames );
    NamedBodyMap bodyMap = createBodies( bodySettings );
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    // Run simulation for 2 different central body settings (barycentric and hierarchical)
    for( int centralBodySettings = 0; centralBodySettings < 2; centralBodySettings++ )
    {
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             CREATE ACCELERATIONS            ///////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Set accelerations between bodies that are to be taken into account (mutual point mass gravity between all bodies).
        SelectedAccelerationMap accelerationMap;
        for( unsigned int i = 0; i < bodyNames.size( ); i++ )
        {
            std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > currentAccelerations;
            for( unsigned int j = 0; j < bodyNames.size( ); j++ )
            {
                // Create central gravity acceleration between each 2 bodies.
                if( i != j )
                {
                    currentAccelerations[ bodyNames.at( j ) ].push_back(
                                std::make_shared< AccelerationSettings >( central_gravity ) );
                }
            }
            accelerationMap[ bodyNames.at( i ) ] = currentAccelerations;
        }

        // Define list of bodies to propagate
        std::vector< std::string > bodiesToPropagate = bodyNames;
        unsigned int numberOfNumericalBodies = bodiesToPropagate.size( );

        // Define central bodies to use in propagation.
        std::vector< std::string > centralBodies;
        centralBodies.resize( numberOfNumericalBodies );

        // Set central body as Solar System Barycente for each body
        if( centralBodySettings == 0 )
        {
            for( unsigned int i = 0; i < numberOfNumericalBodies; i++ )
            {
                centralBodies[ i ] = "SSB";
            }
        }
        else if( centralBodySettings == 1 )
        {
            for( unsigned int i = 0; i < numberOfNumericalBodies; i++ )
            {
                // Set Earth as central body for Moon
                if( i == 0 )
                {
                    centralBodies[ i ] = "Earth";
                }
                // Set barycenter as central 'body' for Sun
                else if( i == 5 )
                {
                    centralBodies[ i ] = "SSB";
                }
                // Set Sun as central body for all planets
                else
                {
                    centralBodies[ i ] = "Sun";
                }
            }
        }

        // Create acceleration models and propagation settings.
        AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             CREATE PROPAGATION SETTINGS            ///////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Specify initial time.
        double initialEphemerisTime = 1.0E7;
        double finalEphemerisTime = 1.0E7 + 5.0 * physical_constants::JULIAN_YEAR;

        // Get initial state vector as input to integration.
        Eigen::VectorXd systemInitialState = getInitialStatesOfBodies(
                    bodiesToPropagate, centralBodies, bodyMap, initialEphemerisTime );

        // Define propagator settings.
        std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, finalEphemerisTime );

        // Define numerical integrator settings.
        std::shared_ptr< IntegratorSettings< > > integratorSettings =
                std::make_shared< IntegratorSettings< > >
                ( rungeKutta4, initialEphemerisTime, 3600.0 );

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             PROPAGATE ORBITS            ///////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Create simulation object and propagate dynamics.
        SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap, integratorSettings, propagatorSettings );

        std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

        // Retrieve numerically integrated state for each body.
        std::vector< std::map< double, Eigen::VectorXd > > allBodiesPropagationHistory;
        allBodiesPropagationHistory.resize( numberOfNumericalBodies );
        for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = integrationResult.begin( );
             stateIterator != integrationResult.end( ); stateIterator++ )
        {
            for( unsigned int i = 0; i < numberOfNumericalBodies; i++ )
            {
                allBodiesPropagationHistory[ i ][ stateIterator->first ] = stateIterator->second.segment( i * 6, 6 );
            }
        }

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////        PROVIDE OUTPUT TO FILES           ////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        std::string outputSubFolder = "InnerSolarSystemPropagationExample/";

        for( unsigned int i = 0; i < numberOfNumericalBodies; i++ )
        {
            // Write propagation history to file.
            input_output::writeDataMapToTextFile(
                        allBodiesPropagationHistory[ i ],
                        "innerSolarSystemPropagationHistory" + bodyNames.at( i ) +
                        boost::lexical_cast< std::string >( centralBodySettings ) + ".dat",
                        tudat_applications::getOutputPath( ) + outputSubFolder,
                        "",
                        std::numeric_limits< double >::digits10,
                        std::numeric_limits< double >::digits10,
                        "," );
        }
    }

    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}
