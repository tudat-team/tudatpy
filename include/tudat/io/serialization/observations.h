/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SERIALIZATION_OBSERVATIONS_H
#define TUDAT_SERIALIZATION_OBSERVATIONS_H

/**
 * @file serialization/observations.h
 * @brief Serialization support for observation-related classes.
 *
 * Registers all polymorphic types in the ObservationDependentVariableSettings
 * hierarchy with Boost.Serialization. This header must be included in exactly
 * one translation unit (the one that performs observation serialization).
 *
 * Usage in pybind11 bindings:
 * @code
 * #include "tudat/io/serialization.h"
 *
 * .def(py::pickle(
 *     [](const ObservationCollection<double, double>& obj) {
 *         return py::bytes(tudat::serialization::serializeToBinaryString(obj));
 *     },
 *     [](py::bytes data) {
 *         return tudat::serialization::deserializeFromBinaryString<
 *             ObservationCollection<double, double>>(data.cast<std::string>());
 *     }
 * ))
 * @endcode
 */

#include "tudat/io/serialization/base.h"
#include "tudat/simulation/estimation_setup/observationCollection.h"
#include "tudat/simulation/estimation_setup/singleObservationSet.h"
#include "tudat/simulation/estimation_setup/observationOutput.h"
#include "tudat/simulation/estimation_setup/observationOutputSettings.h"
#include "tudat/astro/observation_models/observationAncillarySettings.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"

// Register polymorphic types for ObservationDependentVariableSettings hierarchy
// This header must be included in exactly one translation unit
BOOST_CLASS_EXPORT_IMPLEMENT( tudat::simulation_setup::ObservationDependentVariableSettings )
BOOST_CLASS_EXPORT_IMPLEMENT( tudat::simulation_setup::StationAngleObservationDependentVariableSettings )
BOOST_CLASS_EXPORT_IMPLEMENT( tudat::simulation_setup::InterlinkObservationDependentVariableSettings )
BOOST_CLASS_EXPORT_IMPLEMENT( tudat::simulation_setup::AncillaryObservationDependentVariableSettings )

#endif  // TUDAT_SERIALIZATION_OBSERVATIONS_H
