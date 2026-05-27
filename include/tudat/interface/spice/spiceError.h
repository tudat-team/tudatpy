/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SPICE_ERROR_H
#define TUDAT_SPICE_ERROR_H

#include <stdexcept>
#include <string>

namespace tudat
{

namespace exceptions
{

class SpiceError : public std::runtime_error
{
public:
    SpiceError( const std::string& shortMessage,
                const std::string& explanation,
                const std::string& longMessage,
                const std::string& traceback ):
        std::runtime_error( shortMessage + "\n" + explanation + "\n" + longMessage + "\n" + "Traceback: " + traceback )
    {}

    ~SpiceError( ) {}
};

void throwSpiceException( const std::string& shortMessage,
                          const std::string& explanation,
                          const std::string& longMessage,
                          const std::string& traceback );

}  // namespace exceptions

}  // namespace tudat

#endif  // TUDAT_SPICE_ERROR_H
