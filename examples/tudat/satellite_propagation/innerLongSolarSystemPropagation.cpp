/*    Copyright (c) 2010-2018, Delft University of Technology
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
    time_t tstart, tend;
    tstart = time(0);
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
    spice_interface::loadStandardSpiceKernels( );

    // Define bodies in simulation.
    unsigned int totalNumberOfBodies = 7;
    std::vector< std::string > bodyNames;
    bodyNames.resize( totalNumberOfBodies );
    bodyNames[ 0 ] = "Moon";
    bodyNames[ 1 ] = "Earth";
    bodyNames[ 2 ] = "Mars";
    bodyNames[ 3 ] = "Venus";
    bodyNames[ 4 ] = "Mercury";
    bodyNames[ 5 ] = "Sun";
    bodyNames[ 6 ] = "Jupiter";

    // Create bodies needed in simulation
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodyNames );
    NamedBodyMap bodyMap = createBodies( bodySettings );
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    // Run simulation for 2 different central body settings (barycentric and hierarchical)
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
                                std::make_shared< AccelerationSettings >( central_gravity ) );\
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

        for( unsigned int i = 0; i < numberOfNumericalBodies; i++ )
        {
            // Set Earth as central body for Moon
            if( i == 0 )
            {
                centralBodies[ i ] = "Earth";
            }
            else
            {
                centralBodies[ i ] = "SSB";
            }
        }


        // Create acceleration models and propagation settings.
        AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodyMap, accelerationMap, bodiesToPropagate, centralBodies );



        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             CREATE PROPAGATION SETTINGS            ///////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Specify initial time
        double initialEphemerisTime = 0.0;
        double finalEphemerisTime = 1000.0 * physical_constants::JULIAN_YEAR;

        // Get initial state vector as input to integration.
        Eigen::VectorXd systemInitialState = getInitialStatesOfBodies(
                    bodiesToPropagate, centralBodies, bodyMap, initialEphemerisTime );


        std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, finalEphemerisTime,
                  cowell, std::shared_ptr< DependentVariableSaveSettings >( ) );

        // Define numerical integrator settings.
        std::shared_ptr< IntegratorSettings< > > integratorSettings =
                std::make_shared< BulirschStoerIntegratorSettings< double > >(
                    initialEphemerisTime, 3600.0, bulirsch_stoer_sequence, 6,
                    std::numeric_limits< double >::epsilon( ), std::numeric_limits< double >::infinity( ),
                    1.0E-10, 1.0E-8, 10 );

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             PROPAGATE ORBITS            ///////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Create simulation object and propagate dynamics.
        SingleArcDynamicsSimulator< > dynamicsSimulator(
                    bodyMap, integratorSettings, propagatorSettings, true, false, false );

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
                        "innerLongSolarSystemPropagationHistory" + bodyNames.at( i ) + ".dat",
                        tudat_applications::getOutputPath( ) + outputSubFolder,
                        "",
                        std::numeric_limits< double >::digits10,
                        std::numeric_limits< double >::digits10,
                        "," );
        }
    }

    tend = time(0);
    std::cout << "It took "<< difftime(tend, tstart) <<" second(s)."<< std::endl;

    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}

