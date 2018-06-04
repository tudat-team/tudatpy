/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <getopt.h>

#include "Tudat/JsonInterface/jsonInterface.h"

void printHelp( )
{
    std::cout <<
                 "Usage:\n"
                 "\n"
                 "tudat [options] [path]\n"
                 "\n"
                 "path: absolute or relative path to a JSON input file or directory containing a main.json file. "
                 "If not provided, a main.json file will be looked for in the current directory.\n"
                 "\n"
                 "Options:\n"
                 "-h, --help       Show help\n"
              << std::endl;
    exit( EXIT_FAILURE );
}

//! Execute propagation of orbit of Asterix around the Earth.
int main( int argumentCount, char* arguments[ ] )
{
    int currentOption;
    int optionCount = 0;
    const char* const shortOptions = "h";
    const option longOptions[ ] =
    {
        { "help", no_argument, nullptr, 'h' },
        { nullptr, 0, nullptr, 0 }
    };

    while ( ( currentOption = getopt_long( argumentCount, arguments, shortOptions, longOptions, nullptr ) ) != -1 )
    {
        switch ( currentOption )
        {
        case 'h':
        case '?':
        default:
            printHelp( );
        }
        optionCount++;
    }

    const int nonOptionArgumentCount = argumentCount - optionCount - 1;
    if ( nonOptionArgumentCount > 1 )
    {
        printHelp( );
    }
    const std::string inputPath = nonOptionArgumentCount == 1 ? arguments[ argumentCount - 1 ] : "";

    // FIXME: Get binary path (not working on Mac OS)
    // boost::filesystem::path full_path( boost::filesystem::initial_path< boost::filesystem::path >( ) );
    // full_path = boost::filesystem::system_complete( boost::filesystem::path( arguments[ 0 ] ) );
    // std::cout << full_path << std::endl;

    tudat::json_interface::JsonSimulationManager< > jsonSimulationManager( inputPath );
    jsonSimulationManager.updateSettings( );
    jsonSimulationManager.runPropagation( );
    jsonSimulationManager.exportResults( );

    return EXIT_SUCCESS;
}
