/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/InputOutput/tabulatedAtmosphereReader.h"

namespace tudat
{

namespace input_output
{

//! Function to compare if two lists of aerodynamic coefficient independent variables are equal
bool compareIndependentVariables( const std::vector< std::vector< double > >& list1,
                                  const std::vector< std::vector< double > >& list2 )
{
    // Check if sizes of outer vectors are equal
    bool isEqual = 1;
    if( list1.size( ) != list2.size( ) )
    {
        isEqual = 0;
    }

    // Check if sizes of inner vectors are equal
    if( isEqual )
    {
        for( unsigned int i = 0; i < list1.size( ); i++ )
        {
            if( list1.at( i ).size( ) != list2.at( i ).size( ) )
            {
                isEqual = 0;
            }
            else
            {
                for( unsigned int j = 0; j < list1.at( i ).size( ); j++ )
                {
                    if( list1.at( i ).at( j ) != list2.at( i ).at( j ) )
                    {
                        isEqual = 0;
                    }
                }
            }
        }
    }

    return isEqual;
}

}

}
