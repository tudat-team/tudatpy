/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SERIALIZATION_PROPAGATION_H
#define TUDAT_SERIALIZATION_PROPAGATION_H

/**
 * @file serialization/propagation.h
 * @brief Serialization support for propagation result classes and estimation output.
 *
 * Registers all polymorphic types in the SimulationResults hierarchy,
 * PropagationTerminationDetails hierarchy, and estimation output types
 * with Boost.Serialization. This header must be included in exactly
 * one translation unit (via the master serialization.h header).
 */

#include "tudat/io/serialization/base.h"
#include "tudat/simulation/propagation_setup/propagationResults.h"
#include "tudat/simulation/propagation_setup/propagationTermination.h"
#include "tudat/astro/orbit_determination/podInputOutputTypes.h"

// --- PropagationTerminationDetails hierarchy ---
BOOST_CLASS_EXPORT_IMPLEMENT( tudat::propagators::PropagationTerminationDetails )
BOOST_CLASS_EXPORT_IMPLEMENT( tudat::propagators::PropagationTerminationDetailsFromHybridCondition )

// --- SimulationResults hierarchy (all concrete <double, double> instantiations) ---
BOOST_CLASS_EXPORT_IMPLEMENT( tudat::propagators::SingleArcDynamicsResults )
BOOST_CLASS_EXPORT_IMPLEMENT( tudat::propagators::SingleArcVariationalResults )
BOOST_CLASS_EXPORT_IMPLEMENT( tudat::propagators::MultiArcDynamicsResults )
BOOST_CLASS_EXPORT_IMPLEMENT( tudat::propagators::MultiArcVariationalResults )
BOOST_CLASS_EXPORT_IMPLEMENT( tudat::propagators::HybridArcDynamicsResults )
BOOST_CLASS_EXPORT_IMPLEMENT( tudat::propagators::HybridArcVariationalResults )

// --- SimulationResults hierarchy (<double, Time> instantiations for Python bindings) ---
BOOST_CLASS_EXPORT_IMPLEMENT( tudat::propagators::SingleArcDynamicsResultsDT )
BOOST_CLASS_EXPORT_IMPLEMENT( tudat::propagators::SingleArcVariationalResultsDT )
BOOST_CLASS_EXPORT_IMPLEMENT( tudat::propagators::MultiArcDynamicsResultsDT )
BOOST_CLASS_EXPORT_IMPLEMENT( tudat::propagators::MultiArcVariationalResultsDT )
BOOST_CLASS_EXPORT_IMPLEMENT( tudat::propagators::HybridArcDynamicsResultsDT )
BOOST_CLASS_EXPORT_IMPLEMENT( tudat::propagators::HybridArcVariationalResultsDT )

// --- Estimation output types ---
BOOST_CLASS_EXPORT_IMPLEMENT( tudat::simulation_setup::CovarianceAnalysisOutputDD )
BOOST_CLASS_EXPORT_IMPLEMENT( tudat::simulation_setup::EstimationOutputDD )
BOOST_CLASS_EXPORT_IMPLEMENT( tudat::simulation_setup::CovarianceAnalysisOutputDT )
BOOST_CLASS_EXPORT_IMPLEMENT( tudat::simulation_setup::EstimationOutputDT )

#endif  // TUDAT_SERIALIZATION_PROPAGATION_H
