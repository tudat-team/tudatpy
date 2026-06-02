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
 * @brief Master header for Tudat cereal serialization infrastructure.
 *
 * Provides core serialization helpers (Eigen, archives, serialize/deserialize
 * utilities) backed by cereal. Polymorphic type registrations are handled
 * via CEREAL_REGISTER_TYPE in the individual class headers, so no separate
 * "implement" headers are needed.
 *
 * Sub-headers:
 *   - serialization/base.h   Core infrastructure (Eigen, helpers, archives)
 *   - serialization/class_versions.h   Central manifest of cereal class versions
 */

#include "tudat/io/serialization/base.h"
#include "tudat/io/serialization/class_versions.h"

#endif  // TUDAT_SERIALIZATION_H
