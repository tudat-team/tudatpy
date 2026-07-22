/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved.
 */

#ifndef TUDAT_SERIALIZATION_REGISTRATIONS_ACCELERATION_H
#define TUDAT_SERIALIZATION_REGISTRATIONS_ACCELERATION_H

/**
 * @file registrations_acceleration.h
 * @brief Retains and initializes acceleration/torque polymorphic registrations.
 */

#include <cereal/types/polymorphic.hpp>

CEREAL_FORCE_DYNAMIC_INIT( tudat_serialization_acceleration )

#endif  // TUDAT_SERIALIZATION_REGISTRATIONS_ACCELERATION_H
