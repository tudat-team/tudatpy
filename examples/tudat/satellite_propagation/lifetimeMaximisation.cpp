/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <tudat/interface/json/jsonInterface.h>

#include <tudat/io/applicationOutput.h>


//! Execute propagation of orbits of Apollo during entry using the JSON Interface.
int main( )
{
    const std::string cppFilePath( __FILE__ );
    const std::string cppFolder = cppFilePath.substr( 0, cppFilePath.find_last_of("/\\") + 1 );
    tudat::json_interface::JsonSimulationManager< > jsonSimulationManager( cppFolder + "lifetimeMaximisation.json" );

    const std::string outputDirectory = tudat_applications::getOutputPath( ) + "LifetimeMaximisation/";
    const unsigned int numberOfCases = 365;
    for ( unsigned int i = 0; i < numberOfCases; ++i )
    {
        // Notify on propagation start
        std::cout << "Running propagation " << i + 1 << " of " << numberOfCases << "..." << std::endl;

        // Define the initial and final epochs
        const double initialEpoch = i * tudat::physical_constants::JULIAN_DAY;
        jsonSimulationManager[ "initialEpoch" ] = initialEpoch;
        jsonSimulationManager[ "finalEpoch" ] = initialEpoch + tudat::physical_constants::JULIAN_YEAR;

        // Define the output file
        jsonSimulationManager[ "export" ][ 0 ][ "file" ] = outputDirectory + "day" + std::to_string( i + 1 ) + ".dat";

        // Create settings objects
        jsonSimulationManager.updateSettings( );

        // Propagate
        jsonSimulationManager.runPropagation( );

        // Export results
        jsonSimulationManager.exportResults( );

        // Silence unused key warnings after first propagation
        // In this example, when determining the initial conditions, the property "bodies.satellite.initialState" is
        // converted from Keplerian to Cartesian and assigned to "propagators[0].initialStates". Thus, after the first
        // propagation, both keys are defined, but only one is actually used, resulting in warnings about unused keys.
        if ( i == 0 )
        {
            jsonSimulationManager[ "options" ][ "unusedKey" ] = tudat::json_interface::continueSilently;
        }
    }

    return EXIT_SUCCESS;
}

