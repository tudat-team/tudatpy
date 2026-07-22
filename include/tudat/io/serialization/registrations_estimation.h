/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved.
 */

#ifndef TUDAT_SERIALIZATION_REGISTRATIONS_ESTIMATION_H
#define TUDAT_SERIALIZATION_REGISTRATIONS_ESTIMATION_H

/**
 * @file registrations_estimation.h
 * @brief Retains and initializes estimation polymorphic registrations.
 */

#if TUDAT_BUILD_WITH_SERIALIZATION
#include <cereal/types/polymorphic.hpp>
CEREAL_FORCE_DYNAMIC_INIT( tudat_serialization_estimation )
#endif

#endif  // TUDAT_SERIALIZATION_REGISTRATIONS_ESTIMATION_H
