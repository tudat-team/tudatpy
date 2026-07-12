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
 * Includes all domain-specific registration headers.  Include this header
 * (once) in any translation unit that performs serialization of any
 * Tudat polymorphic type.
 *
 * Individual domain headers (registrations_*.h) are also available for
 * finer-grained inclusion if only a subset of types is needed.
 *
 * @note All CEREAL_REGISTER_TYPE and CEREAL_REGISTER_POLYMORPHIC_RELATION
 *       calls are centralized here.  Do not place them in individual class
 *       headers.
 */

#include "tudat/io/serialization/registrations_acceleration.h"
#include "tudat/io/serialization/registrations_environment.h"
#include "tudat/io/serialization/registrations_estimation.h"
#include "tudat/io/serialization/registrations_math.h"
#include "tudat/io/serialization/registrations_propagation.h"

#endif  // TUDAT_SERIALIZATION_REGISTRATIONS_H