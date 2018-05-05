/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_MULTI_DIMENSIONAL_ARRAY_WRITER_H
#define TUDAT_MULTI_DIMENSIONAL_ARRAY_WRITER_H

#include <map>

#include <boost/multi_array.hpp>

#include <Eigen/Core>

#include "Tudat/Basics/basicTypedefs.h"

namespace tudat
{

namespace input_output
{

template< unsigned int NumberOfIndependentVariables, unsigned int NumberOfCoefficients >
void writeMultiArrayAndIndependentVariablesToFiles( const std::map< int, std::string >& fileNamesMap,
                                                    const std::vector< std::vector< double > >& independentVariables,
                                                    const boost::multi_array< Eigen::Matrix< double, NumberOfCoefficients,
                                                    1 >, NumberOfIndependentVariables >& dependentVariables )
{
    // Check consistency of inputs
    // Check that number of coefficients matches number of files
    if ( fileNamesMap.size( ) != NumberOfCoefficients )
    {
        throw std::runtime_error( "Error: inconsistency in inputs. Number of coefficients does not match "
                                  "number of input files." );
    }
    // Check that number of independent variables matches vector size
    if ( independentVariables.size( ) != NumberOfIndependentVariables )
    {
        throw std::runtime_error( "Error: inconsistency in inputs. Number of independent variables does not match "
                                  "size of independent variables vector." );
    }

    switch ( NumberOfIndependentVariables )
    {
    case 1:
    {
        // Open file
        FILE * fileIdentifier = std::fopen( fileNamesMap.at( 0 ).c_str( ), "w" );

        // Loop over each coefficient
        for ( unsigned int j = 0; j < independentVariables.at( 0 ).size( ); j++ ) // loop over rows
        {
            for ( unsigned int i = 0; i < NumberOfCoefficients + 1; i++ )
            {
                if ( i == 0 )
                {
                    fprintf( fileIdentifier, "%f ", independentVariables.at( 0 ).at( j ) );
                }
                else
                {
                    fprintf( fileIdentifier, "%f ", dependentVariables[ j ][ i - 1 ] );
                }
            }
            fprintf( fileIdentifier, "\n" );
        }
        break;
    }
    default:
    {
        // Loop over each file
        for ( unsigned int i = 0; i < NumberOfCoefficients; i++ )
        {
            // Open file
            FILE * fileIdentifier = std::fopen( fileNamesMap.at( i ).c_str( ), "w" );

            // Print number of independent variables
            fprintf( fileIdentifier, "%d\n\n", NumberOfIndependentVariables );

            // Print independent variables
            for ( unsigned int j = 0; j < NumberOfIndependentVariables; j++ )
            {
                for ( unsigned int k = 0; k < independentVariables.at( j ).size( ); k++ )
                {
                    fprintf( fileIdentifier, "%f ", independentVariables.at( j ).at( k ) );
                }
                fprintf( fileIdentifier, "\n" );
            }
            fprintf( fileIdentifier, "\n" );

            // Print dependent variables
            switch ( NumberOfIndependentVariables )
            {
            case 2:
            {
                for ( unsigned int k = 0; k < independentVariables.at( 0 ).size( ); k++ )
                {
                    for ( unsigned int l = 0; l < independentVariables.at( 1 ).size( ); l++ )
                    {
                        fprintf( fileIdentifier, "%f ", dependentVariables[ k ][ l ][ i ] );
                    }
                    fprintf( fileIdentifier, "\n" );
                }
                break;
            }
            case 3:
            {
                for ( unsigned int j = 0; j < independentVariables.at( 2 ).size( ); j++ )
                {
                    for ( unsigned int k = 0; k < independentVariables.at( 0 ).size( ); k++ )
                    {
                        for ( unsigned int l = 0; l < independentVariables.at( 1 ).size( ); l++ )
                        {
                            fprintf( fileIdentifier, "%f ", dependentVariables[ k ][ l ][ j ][ i ] );
                        }
                        fprintf( fileIdentifier, "\n" );
                    }
                    fprintf( fileIdentifier, "\n" );
                }
                break;
            }
            default:
                throw std::runtime_error( "Error: number of independent variables not currently supported "
                                          "by the writer of multi-dimensional arrays." );
            }

            // Close file
            std::fclose( fileIdentifier );
        }
        break;
    }
    }
}

} // namespace input_output

} // namespace tudat

#endif // TUDAT_MULTI_DIMENSIONAL_ARRAY_WRITER_H
