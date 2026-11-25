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
 * @file registrationsPropagationTermination.cpp
 * @brief Centralized cereal polymorphic registrations for propagation termination types.
 *
 * Includes registrations for:
 *   - PropagationTerminationSettings hierarchy (propagationTerminationSettings.h)
 *   - PropagationTerminationDetails hierarchy (propagationTermination.h)
 *
 * Consumers should include registrations_propagation.h to retain this translation unit.
 */

#include "tudat/simulation/propagation_setup/propagationTerminationSettings.h"
#include "tudat/simulation/propagation_setup/propagationTermination.h"

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

CEREAL_REGISTER_DYNAMIC_INIT( tudat_serialization_propagation_termination )
