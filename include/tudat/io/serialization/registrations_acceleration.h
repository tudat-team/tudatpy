/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SERIALIZATION_REGISTRATIONS_ACCELERATION_H
#define TUDAT_SERIALIZATION_REGISTRATIONS_ACCELERATION_H

/**
 * @file registrations_acceleration.h
 * @brief Centralized cereal polymorphic registrations for acceleration and torque settings.
 *
 * Include this header (once) in any translation unit that serializes
 * AccelerationSettings, TorqueSettings, or any of their derived types.
 *
 * @note This header is part of the centralized serialization registration
 *       system.  Do not place CEREAL_REGISTER_TYPE or
 *       CEREAL_REGISTER_POLYMORPHIC_RELATION calls in individual class headers
 *       any more — add them here instead.
 */

#include "tudat/simulation/propagation_setup/accelerationSettings.h"
#include "tudat/simulation/propagation_setup/torqueSettings.h"

// =====================================================================
// AccelerationSettings hierarchy
// =====================================================================

CEREAL_REGISTER_TYPE( tudat::simulation_setup::AccelerationSettings )
CEREAL_REGISTER_TYPE( tudat::simulation_setup::RadiationPressureAccelerationSettings )
CEREAL_REGISTER_TYPE( tudat::simulation_setup::SphericalHarmonicAccelerationSettings )
CEREAL_REGISTER_TYPE( tudat::simulation_setup::MutualSphericalHarmonicAccelerationSettings )
CEREAL_REGISTER_TYPE( tudat::simulation_setup::FullTwoBodySphericalHarmonicAccelerationSettings )
CEREAL_REGISTER_TYPE( tudat::simulation_setup::RelativisticAccelerationCorrectionSettings )
CEREAL_REGISTER_TYPE( tudat::simulation_setup::EmpiricalAccelerationSettings )
CEREAL_REGISTER_TYPE( tudat::simulation_setup::YarkovskyAccelerationSettings )
CEREAL_REGISTER_TYPE( tudat::simulation_setup::ThrustAccelerationSettings )
CEREAL_REGISTER_TYPE( tudat::simulation_setup::CustomAccelerationSettings )
CEREAL_REGISTER_TYPE( tudat::simulation_setup::RTGAccelerationSettings )
CEREAL_REGISTER_TYPE( tudat::simulation_setup::DirectTidalDissipationAccelerationSettings )
CEREAL_REGISTER_TYPE( tudat::simulation_setup::MomentumWheelDesaturationAccelerationSettings )

CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::simulation_setup::AccelerationSettings,
                                      tudat::simulation_setup::RadiationPressureAccelerationSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::simulation_setup::AccelerationSettings,
                                      tudat::simulation_setup::SphericalHarmonicAccelerationSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::simulation_setup::AccelerationSettings,
                                      tudat::simulation_setup::MutualSphericalHarmonicAccelerationSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::simulation_setup::AccelerationSettings,
                                      tudat::simulation_setup::FullTwoBodySphericalHarmonicAccelerationSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::simulation_setup::AccelerationSettings,
                                      tudat::simulation_setup::RelativisticAccelerationCorrectionSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::simulation_setup::AccelerationSettings,
                                      tudat::simulation_setup::EmpiricalAccelerationSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::simulation_setup::AccelerationSettings,
                                      tudat::simulation_setup::YarkovskyAccelerationSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::simulation_setup::AccelerationSettings, tudat::simulation_setup::ThrustAccelerationSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::simulation_setup::AccelerationSettings, tudat::simulation_setup::CustomAccelerationSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::simulation_setup::AccelerationSettings, tudat::simulation_setup::RTGAccelerationSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::simulation_setup::AccelerationSettings,
                                      tudat::simulation_setup::DirectTidalDissipationAccelerationSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::simulation_setup::AccelerationSettings,
                                      tudat::simulation_setup::MomentumWheelDesaturationAccelerationSettings )

// =====================================================================
// TorqueSettings hierarchy
// =====================================================================

CEREAL_REGISTER_TYPE( tudat::simulation_setup::TorqueSettings )
CEREAL_REGISTER_TYPE( tudat::simulation_setup::SphericalHarmonicTorqueSettings )
CEREAL_REGISTER_TYPE( tudat::simulation_setup::FullTwoBodySphericalHarmonicTorqueSettings )
CEREAL_REGISTER_TYPE( tudat::simulation_setup::FourthDegreeFullTwoBodyGravitationalTorqueSettings )
CEREAL_REGISTER_TYPE( tudat::simulation_setup::CustomTorqueSettings )

CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::simulation_setup::TorqueSettings, tudat::simulation_setup::SphericalHarmonicTorqueSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::simulation_setup::TorqueSettings,
                                      tudat::simulation_setup::FullTwoBodySphericalHarmonicTorqueSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::simulation_setup::TorqueSettings,
                                      tudat::simulation_setup::FourthDegreeFullTwoBodyGravitationalTorqueSettings )
CEREAL_REGISTER_POLYMORPHIC_RELATION( tudat::simulation_setup::TorqueSettings, tudat::simulation_setup::CustomTorqueSettings )

#endif  // TUDAT_SERIALIZATION_REGISTRATIONS_ACCELERATION_H
