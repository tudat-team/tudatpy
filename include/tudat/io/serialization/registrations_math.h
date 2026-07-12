/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SERIALIZATION_REGISTRATIONS_MATH_H
#define TUDAT_SERIALIZATION_REGISTRATIONS_MATH_H

/**
 * @file registrations_math.h
 * @brief Centralized cereal polymorphic registrations for math-related types.
 *
 * Includes registrations for:
 *   - InterpolatorSettings hierarchy (createInterpolator.h)
 *
 * Include this header (once) in any translation unit that needs to serialize
 * any of these polymorphic types.
 *
 * @note This header is part of the centralized serialization registration
 *       system.  Do not place CEREAL_REGISTER_TYPE or
 *       CEREAL_REGISTER_POLYMORPHIC_RELATION calls in individual class headers
 *       any more — add them here instead.
 */

#include "tudat/math/interpolators/createInterpolator.h"

// =====================================================================
// InterpolatorSettings hierarchy
// =====================================================================

CEREAL_REGISTER_TYPE( tudat::interpolators::InterpolatorSettings )
CEREAL_REGISTER_TYPE( tudat::interpolators::LagrangeInterpolatorSettings )

CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::interpolators::InterpolatorSettings, tudat::interpolators::LagrangeInterpolatorSettings )

#endif  // TUDAT_SERIALIZATION_REGISTRATIONS_MATH_H