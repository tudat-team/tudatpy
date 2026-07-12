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
 * utilities) backed by cereal.
 *
 * Sub-headers:
 *   - serialization/base.h                Core infrastructure (Eigen, helpers, archives)
 *   - serialization/pybind_helpers.h      Pybind11 pickle / file-IO helpers
 *   - serialization/registrations.h       Centralized CEREAL_REGISTER_TYPE / ... calls
 *
 * All CEREAL_REGISTER_TYPE and CEREAL_REGISTER_POLYMORPHIC_RELATION calls
 * have been centralized in the registrations_*.h headers under
 * tudat/io/serialization/.  Include "registrations.h" (or the specific domain
 * header) in exactly one translation unit per registration group, or include
 * the master header below.
 */

#include "tudat/io/serialization/base.h"
#include "tudat/io/serialization/registrations.h"

#endif  // TUDAT_SERIALIZATION_H
