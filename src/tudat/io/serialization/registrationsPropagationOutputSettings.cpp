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

#include <cereal/archives/binary.hpp>
#include <cereal/archives/json.hpp>

/**
 * @file registrationsPropagationOutputSettings.cpp
 * @brief Centralized cereal polymorphic registrations for propagation output settings.
 *
 * Includes registrations for:
 *   - VariableSettings / SingleDependentVariableSaveSettings hierarchy (propagationOutputSettings.h)
 *
 * Consumers should include registrations_propagation.h to retain this translation unit.
 */

#include "tudat/simulation/propagation_setup/propagationOutputSettings.h"

// =====================================================================
// VariableSettings / DependentVariableSaveSettings hierarchy  (propagationOutputSettings.h)
// =====================================================================

CEREAL_REGISTER_TYPE( tudat::propagators::VariableSettings )
CEREAL_REGISTER_TYPE( tudat::propagators::SingleDependentVariableSaveSettings )
CEREAL_REGISTER_TYPE( tudat::propagators::SingleAccelerationDependentVariableSaveSettings )
CEREAL_REGISTER_TYPE( tudat::propagators::SphericalHarmonicAccelerationTermsDependentVariableSaveSettings )
CEREAL_REGISTER_TYPE( tudat::propagators::SingleTorqueDependentVariableSaveSettings )
CEREAL_REGISTER_TYPE( tudat::propagators::IntermediateAerodynamicRotationVariableSaveSettings )
CEREAL_REGISTER_TYPE( tudat::propagators::BodyAerodynamicAngleVariableSaveSettings )
CEREAL_REGISTER_TYPE( tudat::propagators::LocalWindVelocityDependentVariableSaveSettings )
CEREAL_REGISTER_TYPE( tudat::propagators::ControlSurfaceCoefficientDependentVariableSettings )
CEREAL_REGISTER_TYPE( tudat::propagators::SingleVariationSphericalHarmonicAccelerationSaveSettings )
CEREAL_REGISTER_TYPE( tudat::propagators::SingleVariationSingleTermSphericalHarmonicAccelerationSaveSettings )
CEREAL_REGISTER_TYPE( tudat::propagators::AccelerationPartialWrtStateSaveSettings )
CEREAL_REGISTER_TYPE( tudat::propagators::TotalAccelerationPartialWrtStateSaveSettings )
CEREAL_REGISTER_TYPE( tudat::propagators::MinimumConstellationDistanceDependentVariableSaveSettings )
CEREAL_REGISTER_TYPE( tudat::propagators::MinimumConstellationStationDistanceDependentVariableSaveSettings )
CEREAL_REGISTER_TYPE( tudat::propagators::CustomDependentVariableSaveSettings )
CEREAL_REGISTER_TYPE( tudat::propagators::TotalGravityFieldVariationSettings )
CEREAL_REGISTER_TYPE( tudat::propagators::IlluminatedPanelFractionDependentVariableSaveSettings )
CEREAL_REGISTER_TYPE( tudat::propagators::CrossSectionDependentVariableSaveSettings )

// VariableSettings -> derived
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::VariableSettings, tudat::propagators::SingleDependentVariableSaveSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::VariableSettings,
                                      tudat::propagators::SingleAccelerationDependentVariableSaveSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::VariableSettings,
                                      tudat::propagators::SphericalHarmonicAccelerationTermsDependentVariableSaveSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::VariableSettings, tudat::propagators::SingleTorqueDependentVariableSaveSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::VariableSettings,
                                      tudat::propagators::IntermediateAerodynamicRotationVariableSaveSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::VariableSettings, tudat::propagators::BodyAerodynamicAngleVariableSaveSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::VariableSettings,
                                      tudat::propagators::LocalWindVelocityDependentVariableSaveSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::VariableSettings,
                                      tudat::propagators::ControlSurfaceCoefficientDependentVariableSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::VariableSettings,
                                      tudat::propagators::SingleVariationSphericalHarmonicAccelerationSaveSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::VariableSettings,
                                      tudat::propagators::SingleVariationSingleTermSphericalHarmonicAccelerationSaveSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::VariableSettings, tudat::propagators::AccelerationPartialWrtStateSaveSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::VariableSettings,
                                      tudat::propagators::TotalAccelerationPartialWrtStateSaveSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::VariableSettings,
                                      tudat::propagators::MinimumConstellationDistanceDependentVariableSaveSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::VariableSettings,
                                      tudat::propagators::MinimumConstellationStationDistanceDependentVariableSaveSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::VariableSettings, tudat::propagators::CustomDependentVariableSaveSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::VariableSettings, tudat::propagators::TotalGravityFieldVariationSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::VariableSettings,
                                      tudat::propagators::IlluminatedPanelFractionDependentVariableSaveSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::VariableSettings, tudat::propagators::CrossSectionDependentVariableSaveSettings )

// SingleDependentVariableSaveSettings -> further derived
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::SingleDependentVariableSaveSettings,
                                      tudat::propagators::SingleAccelerationDependentVariableSaveSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::SingleDependentVariableSaveSettings,
                                      tudat::propagators::SphericalHarmonicAccelerationTermsDependentVariableSaveSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::SingleDependentVariableSaveSettings,
                                      tudat::propagators::SingleTorqueDependentVariableSaveSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::SingleDependentVariableSaveSettings,
                                      tudat::propagators::IntermediateAerodynamicRotationVariableSaveSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::SingleDependentVariableSaveSettings,
                                      tudat::propagators::BodyAerodynamicAngleVariableSaveSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::SingleDependentVariableSaveSettings,
                                      tudat::propagators::LocalWindVelocityDependentVariableSaveSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::SingleDependentVariableSaveSettings,
                                      tudat::propagators::ControlSurfaceCoefficientDependentVariableSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::SingleDependentVariableSaveSettings,
                                      tudat::propagators::SingleVariationSphericalHarmonicAccelerationSaveSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::SingleDependentVariableSaveSettings,
                                      tudat::propagators::SingleVariationSingleTermSphericalHarmonicAccelerationSaveSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::SingleDependentVariableSaveSettings,
                                      tudat::propagators::AccelerationPartialWrtStateSaveSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::SingleDependentVariableSaveSettings,
                                      tudat::propagators::TotalAccelerationPartialWrtStateSaveSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::SingleDependentVariableSaveSettings,
                                      tudat::propagators::MinimumConstellationDistanceDependentVariableSaveSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::SingleDependentVariableSaveSettings,
                                      tudat::propagators::MinimumConstellationStationDistanceDependentVariableSaveSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::SingleDependentVariableSaveSettings,
                                      tudat::propagators::CustomDependentVariableSaveSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::SingleDependentVariableSaveSettings,
                                      tudat::propagators::TotalGravityFieldVariationSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::SingleDependentVariableSaveSettings,
                                      tudat::propagators::IlluminatedPanelFractionDependentVariableSaveSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::SingleDependentVariableSaveSettings,
                                      tudat::propagators::CrossSectionDependentVariableSaveSettings )

CEREAL_REGISTER_DYNAMIC_INIT( tudat_serialization_propagation_output_settings )
