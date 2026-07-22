/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SERIALIZATION_REGISTRATIONS_H
#define TUDAT_SERIALIZATION_REGISTRATIONS_H

/**
 * @file registrations.h
 * @brief Master header for all centralized cereal polymorphic registrations.
 *
 * Includes all enabled domain-specific registration-retention headers. Include
 * this convenience header in a translation unit that serializes polymorphic
 * types spanning several Tudat domains.
 *
 * Individual domain headers (registrations_*.h) are also available for
 * finer-grained inclusion if only a subset of types is needed.
 *
 * @note All CEREAL_REGISTER_TYPE and CEREAL_REGISTER_POLYMORPHIC_RELATION
 *       calls are instantiated once in the serialization implementation
 *       library. These domain headers only force the matching registration
 *       object to be retained and initialized.
 */

#include <tudat/config.hpp>

#include "tudat/io/serialization/registrations_acceleration.h"
#include "tudat/io/serialization/registrations_environment.h"
#if TUDAT_BUILD_WITH_ESTIMATION_TOOLS
#include "tudat/io/serialization/registrations_estimation.h"
#endif
#include "tudat/io/serialization/registrations_math.h"
#include "tudat/io/serialization/registrations_propagation.h"

#endif  // TUDAT_SERIALIZATION_REGISTRATIONS_H
