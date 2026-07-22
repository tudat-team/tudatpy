/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved.
 */

#ifndef TUDAT_SERIALIZATION_BASE_H
#define TUDAT_SERIALIZATION_BASE_H

/**
 * @file serialization/base.h
 * @brief Backward-compatible umbrella for Tudat serialization file I/O.
 *
 * Internal code should include core.h, eigen.h, archives.h, file_io_declarations.h, or
 * file_io.h directly according to the dependency layer it needs.
 */

#if TUDAT_BUILD_WITH_SERIALIZATION
#include "tudat/io/serialization/file_io.h"
#else
// Open-template model headers use these definition macros inside class declarations. Keeping
// no-op definitions here removes those methods, and avoids parsing archive/file-I/O machinery,
// when serialization is disabled.
#define TUDAT_DEFINE_FILE_IO( ... )
#define TUDAT_DEFINE_BINARY_IO( ... )
#define TUDAT_DEFINE_BINARY_IO_POLYMORPHIC( ... )
#define TUDAT_DEFINE_FILE_IO_POLYMORPHIC( ... )
#define TUDAT_IMPLEMENT_FILE_IO( ... )
#define TUDAT_IMPLEMENT_BINARY_IO( ... )
#define TUDAT_IMPLEMENT_BINARY_IO_POLYMORPHIC( ... )
#define TUDAT_IMPLEMENT_FILE_IO_POLYMORPHIC( ... )
#endif

#endif  // TUDAT_SERIALIZATION_BASE_H
