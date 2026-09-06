/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <cereal/types/polymorphic.hpp>

#include "tudat/io/serialization/archives.h"

/**
 * @file registrationsMath.cpp
 * @brief Centralized cereal polymorphic registrations for math-related types.
 *
 * Includes registrations for:
 *   - InterpolatorSettings hierarchy (createInterpolator.h)
 *
 * Consumers should include registrations_math.h to retain this translation unit.
 */

#include "tudat/math/interpolators/createInterpolator.h"

// =====================================================================
// InterpolatorSettings hierarchy
// =====================================================================

CEREAL_REGISTER_TYPE( tudat::interpolators::InterpolatorSettings )
CEREAL_REGISTER_TYPE( tudat::interpolators::LagrangeInterpolatorSettings )

CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::interpolators::InterpolatorSettings, tudat::interpolators::LagrangeInterpolatorSettings )

CEREAL_REGISTER_DYNAMIC_INIT( tudat_serialization_math )
