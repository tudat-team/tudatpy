/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDATPY_SCALAR_TYPES_H
#define TUDATPY_SCALAR_TYPES_H

#include "tudat/config.hpp"
#include "tudat/basics/timeType.h"
#include "quadPrecisionTypeCasters.h"

using tudat::Time;

#ifndef STATE_SCALAR_TYPE
#define STATE_SCALAR_TYPE double
#endif

#ifndef TUDATPY_STATE_SCALAR_BINDING_NAMESPACE
#define TUDATPY_STATE_SCALAR_BINDING_NAMESPACE double_precision
#endif

#ifndef TUDATPY_STATE_SCALAR_SUPPORTS_SERIALIZATION
#define TUDATPY_STATE_SCALAR_SUPPORTS_SERIALIZATION 1
#endif

// Cereal cannot serialize Boost.Multiprecision's cpp_bin_float backend: its
// intrusive Boost.Serialization implementation passes boost::serialization
// wrappers to the cereal archive. Keep serialization on the double bindings,
// while leaving the unsupported methods out of quad-dependent Python classes.
#if TUDATPY_STATE_SCALAR_SUPPORTS_SERIALIZATION
#define TUDATPY_DEF_STATE_SCALAR_PICKLE( ... ) TUDATPY_DEF_PICKLE( __VA_ARGS__ )
#define TUDATPY_DEF_STATE_SCALAR_PICKLE_POLYMORPHIC( ... ) TUDATPY_DEF_PICKLE_POLYMORPHIC( __VA_ARGS__ )
#define TUDATPY_DEF_STATE_SCALAR_PICKLE_POLYMORPHIC_DERIVED( ... ) TUDATPY_DEF_PICKLE_POLYMORPHIC_DERIVED( __VA_ARGS__ )
#define TUDATPY_DEF_STATE_SCALAR_BINARY_IO( ... ) TUDATPY_DEF_BINARY_IO( __VA_ARGS__ )
#define TUDATPY_DEF_STATE_SCALAR_BINARY_IO_POLYMORPHIC( ... ) TUDATPY_DEF_BINARY_IO_POLYMORPHIC( __VA_ARGS__ )
#define TUDATPY_DEF_STATE_SCALAR_FILE_IO_POLYMORPHIC( ... ) TUDATPY_DEF_FILE_IO_POLYMORPHIC( __VA_ARGS__ )
#else
#define TUDATPY_DEF_STATE_SCALAR_PICKLE( ... )
#define TUDATPY_DEF_STATE_SCALAR_PICKLE_POLYMORPHIC( ... )
#define TUDATPY_DEF_STATE_SCALAR_PICKLE_POLYMORPHIC_DERIVED( ... )
#define TUDATPY_DEF_STATE_SCALAR_BINARY_IO( ... )
#define TUDATPY_DEF_STATE_SCALAR_BINARY_IO_POLYMORPHIC( ... )
#define TUDATPY_DEF_STATE_SCALAR_FILE_IO_POLYMORPHIC( ... )
#endif

#define TIME_TYPE Time
#define INTERPOLATOR_TIME_TYPE Time

#define TUDATPY_TYPE_ID_double 1
#define TUDATPY_TYPE_ID_long_double 2
#define TUDATPY_TYPE_ID_Time 3

#define TUDATPY_TYPE_ID_INDIR( x ) TUDATPY_TYPE_ID_##x
#define TUDATPY_TYPE_ID( x ) TUDATPY_TYPE_ID_INDIR( x )

#endif  // TUDATPY_SCALAR_TYPES_H
