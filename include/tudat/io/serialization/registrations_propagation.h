/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved.
 */

#ifndef TUDAT_SERIALIZATION_REGISTRATIONS_PROPAGATION_H
#define TUDAT_SERIALIZATION_REGISTRATIONS_PROPAGATION_H

/**
 * @file registrations_propagation.h
 * @brief Retains and initializes propagation polymorphic registrations.
 */

#include <tudat/config.hpp>

#if TUDAT_BUILD_WITH_SERIALIZATION
#include <cereal/types/polymorphic.hpp>
CEREAL_FORCE_DYNAMIC_INIT( tudat_serialization_propagation )
#endif

#endif  // TUDAT_SERIALIZATION_REGISTRATIONS_PROPAGATION_H
