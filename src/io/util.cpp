//
// Created by ggarrett on 19-05-20.
//

#include "tudat/io/util.h"

namespace tudat
{
namespace input_output{

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