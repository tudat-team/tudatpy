/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SPICE_EXCEPTIONS_H
#define TUDAT_SPICE_EXCEPTIONS_H

#include <iostream>
#include <string>

namespace tudat
{

namespace exceptions
{

/// @brief Generic exception class for SPICE.
///
/// This class is a generic base class for all SPICE exceptions.
/// It inherits from std::runtime_error.
class SpiceError : public std::runtime_error
{
private:
public:
    /// @brief Default constructor for SpiceError.
    /// @param errorMessage Error message to be displayed.
    SpiceError( const std::string& shortMessage, const std::string& longMessage, const std::string& explanation ):
        std::runtime_error( shortMessage + "\n" + longMessage + "\n" + explanation )
    { }
    ~SpiceError( ) { }
};

}  // namespace exceptions

}  // namespace tudat

#endif  // TUDAT_SPICE_EXCEPTIONS_H