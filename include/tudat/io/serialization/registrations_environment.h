/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SERIALIZATION_REGISTRATIONS_ENVIRONMENT_H
#define TUDAT_SERIALIZATION_REGISTRATIONS_ENVIRONMENT_H

/**
 * @file registrations_environment.h
 * @brief Centralized cereal polymorphic registrations for environment-related types.
 *
 * Includes registrations for:
 *   - GravityFieldVariationSettings hierarchy (createGravityFieldVariations.h)
 *
 * Include this header (once) in any translation unit that needs to serialize
 * any of these polymorphic types.
 *
 * @note This header is part of the centralized serialization registration
 *       system.  Do not place CEREAL_REGISTER_TYPE or
 *       CEREAL_REGISTER_POLYMORPHIC_RELATION calls in individual class headers
 *       any more — add them here instead.
 */

#include "tudat/simulation/environment_setup/createGravityFieldVariations.h"

// =====================================================================
// GravityFieldVariationSettings hierarchy
// =====================================================================

CEREAL_REGISTER_TYPE( tudat::simulation_setup::GravityFieldVariationSettings )
CEREAL_REGISTER_TYPE( tudat::simulation_setup::BasicSolidBodyGravityFieldVariationSettings )
CEREAL_REGISTER_TYPE( tudat::simulation_setup::ModeCoupledSolidBodyGravityFieldVariationSettings )
CEREAL_REGISTER_TYPE( tudat::simulation_setup::TabulatedGravityFieldVariationSettings )
CEREAL_REGISTER_TYPE( tudat::simulation_setup::PeriodicGravityFieldVariationsSettings )
CEREAL_REGISTER_TYPE( tudat::simulation_setup::PolynomialGravityFieldVariationsSettings )

CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::simulation_setup::GravityFieldVariationSettings,
                                      tudat::simulation_setup::BasicSolidBodyGravityFieldVariationSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::simulation_setup::GravityFieldVariationSettings,
                                      tudat::simulation_setup::ModeCoupledSolidBodyGravityFieldVariationSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::simulation_setup::GravityFieldVariationSettings,
                                      tudat::simulation_setup::TabulatedGravityFieldVariationSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::simulation_setup::GravityFieldVariationSettings,
                                      tudat::simulation_setup::PeriodicGravityFieldVariationsSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::simulation_setup::GravityFieldVariationSettings,
                                      tudat::simulation_setup::PolynomialGravityFieldVariationsSettings )

#endif  // TUDAT_SERIALIZATION_REGISTRATIONS_ENVIRONMENT_H