/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SPICE_EXCEPTION_TYPES_H
#define TUDAT_SPICE_EXCEPTION_TYPES_H

#include "tudat/interface/spice/spiceError.h"

namespace tudat
{

namespace exceptions
{

#define TUDAT_SPICE_EXCEPTION( short_message, exception_type ) \
    class exception_type : public SpiceError                      \
    {                                                             \
    public:                                                       \
        using SpiceError::SpiceError;                            \
    };

#include "tudat/interface/spice/spiceExceptionList.def"

#undef TUDAT_SPICE_EXCEPTION

}  // namespace exceptions

}  // namespace tudat

#endif  // TUDAT_SPICE_EXCEPTION_TYPES_H
