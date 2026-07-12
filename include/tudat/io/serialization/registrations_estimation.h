/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SERIALIZATION_REGISTRATIONS_ESTIMATION_H
#define TUDAT_SERIALIZATION_REGISTRATIONS_ESTIMATION_H

/**
 * @file registrations_estimation.h
 * @brief Centralized cereal polymorphic registrations for estimation-related types.
 *
 * Includes registrations for:
 *   - ObservationDependentVariableSettings hierarchy (observationOutputSettings.h)
 *   - CovarianceAnalysisOutput / EstimationOutput hierarchy (podInputOutputTypes.h)
 *
 * Include this header (once) in any translation unit that needs to serialize
 * any of these polymorphic types.
 *
 * @note This header is part of the centralized serialization registration
 *       system.  Do not place CEREAL_REGISTER_TYPE or
 *       CEREAL_REGISTER_POLYMORPHIC_RELATION calls in individual class headers
 *       any more — add them here instead.
 */

#include "tudat/simulation/estimation_setup/observationOutputSettings.h"
#include "tudat/astro/orbit_determination/podInputOutputTypes.h"

// =====================================================================
// ObservationDependentVariableSettings hierarchy  (observationOutputSettings.h)
// =====================================================================

CEREAL_REGISTER_TYPE( tudat::simulation_setup::ObservationDependentVariableSettings )
CEREAL_REGISTER_TYPE( tudat::simulation_setup::StationAngleObservationDependentVariableSettings )
CEREAL_REGISTER_TYPE( tudat::simulation_setup::InterlinkObservationDependentVariableSettings )
CEREAL_REGISTER_TYPE( tudat::simulation_setup::AncillaryObservationDependentVariableSettings )
CEREAL_REGISTER_TYPE( tudat::simulation_setup::LightTimeCorrectionComponentsDependentVariableSettings )

CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::simulation_setup::ObservationDependentVariableSettings,
                                      tudat::simulation_setup::StationAngleObservationDependentVariableSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::simulation_setup::ObservationDependentVariableSettings,
                                      tudat::simulation_setup::InterlinkObservationDependentVariableSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::simulation_setup::ObservationDependentVariableSettings,
                                      tudat::simulation_setup::AncillaryObservationDependentVariableSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::simulation_setup::ObservationDependentVariableSettings,
                                      tudat::simulation_setup::LightTimeCorrectionComponentsDependentVariableSettings )

// =====================================================================
// CovarianceAnalysisOutput / EstimationOutput hierarchy  (podInputOutputTypes.h)
// =====================================================================

CEREAL_REGISTER_TYPE( tudat::simulation_setup::CovarianceAnalysisOutputDD )
CEREAL_REGISTER_TYPE( tudat::simulation_setup::EstimationOutputDD )
CEREAL_REGISTER_TYPE( tudat::simulation_setup::CovarianceAnalysisOutputDT )
CEREAL_REGISTER_TYPE( tudat::simulation_setup::EstimationOutputDT )

CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::simulation_setup::CovarianceAnalysisOutputDD, tudat::simulation_setup::EstimationOutputDD )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::simulation_setup::CovarianceAnalysisOutputDT, tudat::simulation_setup::EstimationOutputDT )

#endif  // TUDAT_SERIALIZATION_REGISTRATIONS_ESTIMATION_H