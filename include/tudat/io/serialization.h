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
 * Provides archive/file helpers and retains all Tudat polymorphic registration
 * domains. Prefer a narrower sub-header in implementation code.
 *
 * Sub-headers:
 *   - serialization/core.h                Archive-independent Cereal support
 *   - serialization/eigen.h               Eigen matrix adapter
 *   - serialization/archives.h            Supported archive implementations
 *   - serialization/file_io.h             File/string helpers and definition macros
 *   - serialization/pybind_helpers.h      Pybind11 pickle / file-I/O helpers
 *   - serialization/registrations.h       Registration-retention declarations
 *
 * CEREAL_REGISTER_TYPE and CEREAL_REGISTER_POLYMORPHIC_RELATION definitions are
 * instantiated once in the serialization library. The public registrations_*.h
 * headers only force the matching registration objects to be retained.
 */

#if TUDAT_BUILD_WITH_SERIALIZATION
#include "tudat/io/serialization/base.h"
#include "tudat/io/serialization/registrations.h"
#endif

#endif  // TUDAT_SERIALIZATION_H
