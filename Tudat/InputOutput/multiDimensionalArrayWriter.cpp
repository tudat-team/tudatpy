/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <iostream>
#include <string>

#include "Tudat/InputOutput/multiDimensionalArrayWriter.h"

namespace tudat
{

namespace input_output
{

//! Function to print to console which aerodynamic coefficients are being saved.
void informUserOnSavedCoefficient( const int index )
{
    // Set order of coefficients
    std::vector< std::string > coefficientNames;
    coefficientNames.push_back( "Drag" );
    coefficientNames.push_back( "Side" );
    coefficientNames.push_back( "Lift" );
    coefficientNames.push_back( "X-Moment" );
    coefficientNames.push_back( "Y-Moment" );
    coefficientNames.push_back( "Z-Moment" );

    // Inform user on which variable is being saved
    std::cout << "Saving " + coefficientNames.at( index ) << " aerodynamic coefficient." << std::endl;
}

} // namespace input_output

} // namespace tudat
