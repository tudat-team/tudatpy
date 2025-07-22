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

#include "tudat/basics/timeType.h"

using tudat::Time;

#define STATE_SCALAR_TYPE double  // long double
#define TIME_TYPE Time
#define INTERPOLATOR_TIME_TYPE Time

#define TUDATPY_TYPE_ID_double 1
#define TUDATPY_TYPE_ID_long_double 2
#define TUDATPY_TYPE_ID_Time 3

#define TUDATPY_TYPE_ID_INDIR( x ) TUDATPY_TYPE_ID_##x
#define TUDATPY_TYPE_ID( x ) TUDATPY_TYPE_ID_INDIR( x )

#endif  // TUDATPY_SCALAR_TYPES_H
