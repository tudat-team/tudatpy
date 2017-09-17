/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Astrodynamics/EarthOrientation/eopReader.h"

namespace tudat
{

EOPReader::EOPReader( std::string eopFile,
             std::string format,
             IAUConventions nutationTheory )
{
    if( format != "C04" )
    {
        std::cerr<<"Warning, only C04 EOP file format currently supported by reader "<<std::endl;
    }
    if( !( nutationTheory == iau_2000_a || nutationTheory == iau_2006 ) )
    {
        std::cerr<<"Warning, only IAU2000 nutation theory format currently supported by reader "<<std::endl;
    }
    readEopFile( eopFile );
}

void EOPReader::readEopFile( std::string fileName )
{
    using namespace tudat::unit_conversions;

    // Open file and create file stream.
    std::fstream stream( fileName.c_str( ), std::ios::in );

    // Check if file opened correctly.
    if ( stream.fail( ) )
    {
        boost::throw_exception( std::runtime_error( boost::str(
                                                        boost::format( "Data file '%s' could not be opened." ) % fileName.c_str( ) ) ) );
    }

    // Initialize boolean that gets set to true once the file header is passed.
    bool isHeaderPassed = 0;


    // Line based parsing
    std::string line;
    std::vector< std::string > vectorOfIndividualStrings;
    while ( !stream.fail( ) && !stream.eof( ) )
    {
        // Get line from stream
        std::getline( stream, line );

        // Trim input string (removes all leading and trailing whitespaces).
        boost::algorithm::trim( line );

        // Split string into multiple strings, each containing one element from a line from the
        // data file.
        boost::algorithm::split( vectorOfIndividualStrings,
                                 line,
                                 boost::algorithm::is_any_of( " " ),
                                 boost::algorithm::token_compress_on );

        // If first line before data is found, check entry names and set isHeaderPassed to true.
        if( !isHeaderPassed && vectorOfIndividualStrings[ 0 ] == "Date" )
        {
            if( vectorOfIndividualStrings[ 6 ] == "dPsi" )
            {
                std::cerr<<"Warning, found dPsi, expected dX as CIP offset in GCRS, wrong nutation format requested"<<std::endl;
            }
            else if( vectorOfIndividualStrings[ 7 ] == "dEps" )
            {
                std::cerr<<"Warning, found dEps, expected dY as CIP offset in GCRS, wrong nutation format requested"<<std::endl;
            }

            isHeaderPassed = 1;
        }

        // Read line of data and set values.
        if( isHeaderPassed && vectorOfIndividualStrings.size( ) == 16 )
        {

            // Read pole positions.
            cipInItrs[ boost::lexical_cast< double >( vectorOfIndividualStrings[ 3 ] ) ].x( ) =
                    convertArcSecondsToRadians< double >( boost::lexical_cast< double >( vectorOfIndividualStrings[ 4 ] ) );
            cipInItrs[ boost::lexical_cast< double >( vectorOfIndividualStrings[ 3 ] ) ].y( ) =
                    convertArcSecondsToRadians< double >( boost::lexical_cast< double >( vectorOfIndividualStrings[ 5 ] ) );

            // Read UTC-UT1 and LOD corrections
            ut1MinusUtc[ boost::lexical_cast< double >( vectorOfIndividualStrings[ 3 ] ) ] =
                    boost::lexical_cast< double >( vectorOfIndividualStrings[ 6 ] );
            lengthOfDayOffset[ boost::lexical_cast< double >( vectorOfIndividualStrings[ 3 ] ) ] =
                    boost::lexical_cast< double >( vectorOfIndividualStrings[ 7 ] );

            // Read precession-nutation corrections.
            cipInGcrsCorrection[ boost::lexical_cast< double >( vectorOfIndividualStrings[ 3 ] ) ].x( ) =
                    convertArcSecondsToRadians< double >( boost::lexical_cast< double >( vectorOfIndividualStrings[ 8 ] ) );
            cipInGcrsCorrection[ boost::lexical_cast< double >( vectorOfIndividualStrings[ 3 ] ) ].y( ) =
                    convertArcSecondsToRadians< double >( boost::lexical_cast< double >( vectorOfIndividualStrings[ 9 ] ) );

        }
    }
}

}
