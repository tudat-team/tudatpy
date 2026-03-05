/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/interface/spice/spiceExceptionTypes.h"

namespace tudat
{

namespace exceptions
{

void throwSpiceException( const std::string& shortMessage,
                          const std::string& explanation,
                          const std::string& longMessage,
                          const std::string& traceback )
{
#define TUDAT_SPICE_EXCEPTION( short_message, exception_type ) \
    if( shortMessage == short_message )                         \
    {                                                           \
        throw exception_type( shortMessage, explanation, longMessage, traceback ); \
    }

#include "tudat/interface/spice/spiceExceptionList.def"

#undef TUDAT_SPICE_EXCEPTION

    throw SpiceError( shortMessage, explanation, longMessage, traceback );
}

}  // namespace exceptions

}  // namespace tudat
