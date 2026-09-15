/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved.
 */

#ifndef TUDAT_SERIALIZATION_REGISTRATIONS_MATH_H
#define TUDAT_SERIALIZATION_REGISTRATIONS_MATH_H

/**
 * @file registrations_math.h
 * @brief Retains and initializes math polymorphic registrations.
 */

#if TUDAT_BUILD_WITH_SERIALIZATION
#include <cereal/types/polymorphic.hpp>
CEREAL_FORCE_DYNAMIC_INIT( tudat_serialization_math )
#endif

#endif  // TUDAT_SERIALIZATION_REGISTRATIONS_MATH_H
