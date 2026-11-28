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
 * @file registrationsPropagation.cpp
 * @brief Centralized cereal polymorphic registrations for propagation-related types.
 *
 * Includes registrations for:
 *   - SimulationResults hierarchy (propagationResults.h)
 *   - PropagatorProcessingSettings hierarchy (propagationProcessingSettings.h)
 *   - DependentVariablesInterface hierarchy (dependentVariablesInterface.h)
 *
 * Consumers should include registrations_propagation.h to retain this translation unit.
 */

#include "tudat/simulation/propagation_setup/propagationResults.h"
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

CEREAL_REGISTER_DYNAMIC_INIT( tudat_serialization_propagation )
