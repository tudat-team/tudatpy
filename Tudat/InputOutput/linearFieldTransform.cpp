/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "Tudat/InputOutput/linearFieldTransform.h"

namespace tudat
{
namespace input_output
{

//! Transform input string.
std::shared_ptr< std::string > LinearFieldTransform::transform( const std::string& input )
{
    // Transform string to double.
    const double number = std::stod( input );

    // Perform transformation.
    const double result = slope * number + intercept;

    // Return pointer to transformed value, in string format.
    return std::shared_ptr< std::string >( new std::string( std::to_string( result ) ) );
}

} // namespace input_output
} // namespace tudat
