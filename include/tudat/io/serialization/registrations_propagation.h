/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SERIALIZATION_REGISTRATIONS_PROPAGATION_H
#define TUDAT_SERIALIZATION_REGISTRATIONS_PROPAGATION_H

/**
 * @file registrations_propagation.h
 * @brief Centralized cereal polymorphic registrations for propagation-related types.
 *
 * Includes registrations for:
 *   - SimulationResults hierarchy (propagationResults.h)
 *   - VariableSettings / SingleDependentVariableSaveSettings hierarchy (propagationOutputSettings.h)
 *   - PropagationTerminationSettings hierarchy (propagationTerminationSettings.h)
 *   - PropagationTerminationDetails hierarchy (propagationTermination.h)
 *   - PropagatorProcessingSettings hierarchy (propagationProcessingSettings.h)
 *   - DependentVariablesInterface hierarchy (dependentVariablesInterface.h)
 *
 * Include this header (once) in any translation unit that needs to serialize
 * any of these polymorphic types.
 *
 * @note This header is part of the centralized serialization registration
 *       system.  Do not place CEREAL_REGISTER_TYPE or
 *       CEREAL_REGISTER_POLYMORPHIC_RELATION calls in individual class headers
 *       any more — add them here instead.
 */

#include "tudat/simulation/propagation_setup/propagationResults.h"
#include "tudat/simulation/propagation_setup/propagationOutputSettings.h"
#include "tudat/simulation/propagation_setup/propagationTerminationSettings.h"
#include "tudat/simulation/propagation_setup/propagationTermination.h"
#include "tudat/simulation/propagation_setup/propagationProcessingSettings.h"
#include "tudat/simulation/propagation_setup/dependentVariablesInterface.h"

// =====================================================================
// SimulationResults hierarchy  (propagationResults.h)
// =====================================================================

CEREAL_REGISTER_TYPE( tudat::propagators::SingleArcDynamicsResults )
CEREAL_REGISTER_TYPE( tudat::propagators::SingleArcVariationalResults )
CEREAL_REGISTER_TYPE( tudat::propagators::MultiArcDynamicsResults )
CEREAL_REGISTER_TYPE( tudat::propagators::MultiArcVariationalResults )
CEREAL_REGISTER_TYPE( tudat::propagators::HybridArcDynamicsResults )
CEREAL_REGISTER_TYPE( tudat::propagators::HybridArcVariationalResults )

// <double, Time> variants
CEREAL_REGISTER_TYPE( tudat::propagators::SingleArcDynamicsResultsDT )
CEREAL_REGISTER_TYPE( tudat::propagators::SingleArcVariationalResultsDT )
CEREAL_REGISTER_TYPE( tudat::propagators::MultiArcDynamicsResultsDT )
CEREAL_REGISTER_TYPE( tudat::propagators::MultiArcVariationalResultsDT )
CEREAL_REGISTER_TYPE( tudat::propagators::HybridArcDynamicsResultsDT )
CEREAL_REGISTER_TYPE( tudat::propagators::HybridArcVariationalResultsDT )

// Polymorphic relationships (using typedefs to avoid macro comma issues)
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::SimulationResultsDD, tudat::propagators::SingleArcDynamicsResults )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::SimulationResultsDD, tudat::propagators::SingleArcVariationalResults )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::SimulationResultsDD, tudat::propagators::MultiArcDynamicsResults )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::SimulationResultsDD, tudat::propagators::MultiArcVariationalResults )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::SimulationResultsDD, tudat::propagators::HybridArcDynamicsResults )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::SimulationResultsDD, tudat::propagators::HybridArcVariationalResults )

CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::SimulationResultsDT, tudat::propagators::SingleArcDynamicsResultsDT )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::SimulationResultsDT, tudat::propagators::SingleArcVariationalResultsDT )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::SimulationResultsDT, tudat::propagators::MultiArcDynamicsResultsDT )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::SimulationResultsDT, tudat::propagators::MultiArcVariationalResultsDT )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::SimulationResultsDT, tudat::propagators::HybridArcDynamicsResultsDT )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::SimulationResultsDT, tudat::propagators::HybridArcVariationalResultsDT )

// =====================================================================
// PropagatorProcessingSettings hierarchy
// =====================================================================

CEREAL_REGISTER_TYPE( tudat::propagators::SingleArcPropagatorProcessingSettings )
CEREAL_REGISTER_TYPE( tudat::propagators::MultiArcPropagatorProcessingSettings )
CEREAL_REGISTER_TYPE( tudat::propagators::HybridArcPropagatorProcessingSettings )

CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::PropagatorProcessingSettings,
                                      tudat::propagators::SingleArcPropagatorProcessingSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::PropagatorProcessingSettings,
                                      tudat::propagators::MultiArcPropagatorProcessingSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::PropagatorProcessingSettings,
                                      tudat::propagators::HybridArcPropagatorProcessingSettings )

// =====================================================================
// DependentVariablesInterface hierarchy
// =====================================================================

namespace tudat
{
namespace propagators
{
using DependentVariablesInterfaceDouble = DependentVariablesInterface< double >;
using SingleArcDependentVariablesInterfaceDouble = SingleArcDependentVariablesInterface< double >;
using MultiArcDependentVariablesInterfaceDouble = MultiArcDependentVariablesInterface< double >;
using HybridArcDependentVariablesInterfaceDouble = HybridArcDependentVariablesInterface< double >;
using DependentVariablesInterfaceTime = DependentVariablesInterface< tudat::Time >;
using SingleArcDependentVariablesInterfaceTime = SingleArcDependentVariablesInterface< tudat::Time >;
using MultiArcDependentVariablesInterfaceTime = MultiArcDependentVariablesInterface< tudat::Time >;
using HybridArcDependentVariablesInterfaceTime = HybridArcDependentVariablesInterface< tudat::Time >;
}  // namespace propagators
}  // namespace tudat

CEREAL_REGISTER_TYPE( tudat::propagators::SingleArcDependentVariablesInterfaceDouble )
CEREAL_REGISTER_TYPE( tudat::propagators::MultiArcDependentVariablesInterfaceDouble )
CEREAL_REGISTER_TYPE( tudat::propagators::HybridArcDependentVariablesInterfaceDouble )
CEREAL_REGISTER_TYPE( tudat::propagators::SingleArcDependentVariablesInterfaceTime )
CEREAL_REGISTER_TYPE( tudat::propagators::MultiArcDependentVariablesInterfaceTime )
CEREAL_REGISTER_TYPE( tudat::propagators::HybridArcDependentVariablesInterfaceTime )

CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::DependentVariablesInterfaceDouble,
                                      tudat::propagators::SingleArcDependentVariablesInterfaceDouble )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::DependentVariablesInterfaceDouble,
                                      tudat::propagators::MultiArcDependentVariablesInterfaceDouble )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::DependentVariablesInterfaceDouble,
                                      tudat::propagators::HybridArcDependentVariablesInterfaceDouble )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::DependentVariablesInterfaceTime,
                                      tudat::propagators::SingleArcDependentVariablesInterfaceTime )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::DependentVariablesInterfaceTime,
                                      tudat::propagators::MultiArcDependentVariablesInterfaceTime )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::DependentVariablesInterfaceTime,
                                      tudat::propagators::HybridArcDependentVariablesInterfaceTime )

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

// =====================================================================
// PropagationTerminationSettings hierarchy  (propagationTerminationSettings.h)
// =====================================================================

CEREAL_REGISTER_TYPE( tudat::propagators::PropagationTerminationSettings )
CEREAL_REGISTER_TYPE( tudat::propagators::PropagationTimeTerminationSettings )
CEREAL_REGISTER_TYPE( tudat::propagators::PropagationCPUTimeTerminationSettings )
CEREAL_REGISTER_TYPE( tudat::propagators::PropagationDependentVariableTerminationSettings )
CEREAL_REGISTER_TYPE( tudat::propagators::PropagationCustomTerminationSettings )
CEREAL_REGISTER_TYPE( tudat::propagators::PropagationHybridTerminationSettings )
CEREAL_REGISTER_TYPE( tudat::propagators::NonSequentialPropagationTerminationSettings )

CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::PropagationTerminationSettings,
                                      tudat::propagators::PropagationTimeTerminationSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::PropagationTerminationSettings,
                                      tudat::propagators::PropagationCPUTimeTerminationSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::PropagationTerminationSettings,
                                      tudat::propagators::PropagationDependentVariableTerminationSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::PropagationTerminationSettings,
                                      tudat::propagators::PropagationCustomTerminationSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::PropagationTerminationSettings,
                                      tudat::propagators::PropagationHybridTerminationSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::PropagationTerminationSettings,
                                      tudat::propagators::NonSequentialPropagationTerminationSettings )

// =====================================================================
// PropagationTerminationDetails hierarchy  (propagationTermination.h)
// =====================================================================

CEREAL_REGISTER_TYPE( tudat::propagators::PropagationTerminationDetails )
CEREAL_REGISTER_TYPE( tudat::propagators::PropagationTerminationDetailsFromHybridCondition )

CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::propagators::PropagationTerminationDetails,
                                      tudat::propagators::PropagationTerminationDetailsFromHybridCondition )

#endif  // TUDAT_SERIALIZATION_REGISTRATIONS_PROPAGATION_H
