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
 * @file registrationsEnvironment.cpp
 * @brief Centralized cereal polymorphic registrations for environment-related types.
 *
 * Includes registrations for:
 *   - GravityFieldVariationSettings hierarchy (createGravityFieldVariations.h)
 *
 * Consumers should include registrations_environment.h to retain this translation unit.
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

CEREAL_REGISTER_DYNAMIC_INIT( tudat_serialization_environment )
