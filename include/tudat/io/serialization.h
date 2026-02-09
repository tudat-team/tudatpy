/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SERIALIZATION_H
#define TUDAT_SERIALIZATION_H

/**
 * @file serialization.h
 * @brief Master header for all Tudat Boost serialization support.
 *
 * Include this single header to get access to all serialization infrastructure
 * and polymorphic type registrations. Individual sub-headers can also be
 * included directly if only a subset is needed.
 *
 * Sub-headers:
 *   - serialization/base.h          Core infrastructure (Eigen, helpers, archives)
 *   - serialization/observations.h  Observation classes + polymorphic exports
 */

#include "tudat/io/serialization/base.h"
#include "tudat/io/serialization/observations.h"

#endif  // TUDAT_SERIALIZATION_H
