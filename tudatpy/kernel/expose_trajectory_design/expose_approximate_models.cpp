/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>

#include "expose_approximate_models.h"


#include "tudat/astro/gravitation/unitConversionsCircularRestrictedThreeBodyProblem.h"
#include "tudat/astro/gravitation/jacobiEnergy.h"
#include "tudat/astro/gravitation/librationPoint.h"
#include "tudat/simulation/propagation_setup/createStateDerivativeModel.h"
#include "tudat/simulation/propagation_setup/propagationCR3BPFullProblem.h"
#include "tudat/simulation/propagation_setup/propagationLambertTargeterFullProblem.h"
#include "tudat/astro/basic_astro/clohessyWiltshirePropagator.h"

#include "tudatpy/docstrings.h"

namespace py = pybind11;
//namespace tms = tudat::mission_segments;
namespace tcrp = tudat::circular_restricted_three_body_problem;
namespace tg = tudat::gravitation;
namespace tba = tudat::basic_astrodynamics;
namespace tp = tudat::propagators;
namespace tss = tudat::simulation_setup;
namespace tpc = tudat::physical_constants;
namespace trf = tudat::root_finders;

namespace tudatpy {
namespace trajectory_design {
namespace approximate_models {

void expose_approximate_models(py::module &m)
{
    m.def("libration_point_position",
          &tcrp::computeLibrationPointPosition,
          py::arg("mass_parameter"),
          py::arg("libration_point_index"),
          py::arg("root_finder") = std::make_shared< trf::NewtonRaphson< > >( 1.0e-14, 1000 ),
          get_docstring("libration_point_position").c_str() );

    m.def("jacobi_energy",
          &tg::computeJacobiEnergy,
          py::arg("mass_parameter"),
          py::arg("normalized_corotating_state"),
          get_docstring("jacobi_energy").c_str() );

    m.def("mass_parameter",
          &tcrp::computeMassParameter,
          py::arg("primary_gravitational_parameter"),
          py::arg("secondary_gravitational_parameter"),
          get_docstring("mass_parameter").c_str() );

    m.def("dimensionless_to_dimensional_time",
          &tcrp::convertDimensionlessTimeToDimensionalTime,
          py::arg("dimensionless_time"),
          py::arg("primary_gravitational_parameter"),
          py::arg("secondary_gravitational_parameter"),
          py::arg("primary_secondary_distance"),
          get_docstring("dimensionless_to_dimensional_time").c_str() );

    m.def("dimensional_to_dimensionless_time",
          &tcrp::convertDimensionalTimeToDimensionlessTime,
          py::arg("dimensional_time"),
          py::arg("primary_gravitational_parameter"),
          py::arg("secondary_gravitational_parameter"),
          py::arg("primary_secondary_distance"),
          get_docstring("dimensional_time_to_dimensionless_time").c_str() );

    m.def("corotating_normalized_to_inertial_cartesian",
          &tcrp::convertCorotatingNormalizedToCartesianCoordinates,
          py::arg("primary_gravitational_parameter"),
          py::arg("secondary_gravitational_parameter"),
          py::arg("primary_secondary_distance"),
          py::arg("normalized_corotating_state"),
          py::arg("normaized_time"),
          get_docstring("corotating_normalized_to_inertial_cartesian").c_str() );

    m.def("inertial_cartesian_to_corotating_normalized",
          &tcrp::convertCartesianToCorotatingNormalizedCoordinates,
          py::arg("primary_gravitational_parameter"),
          py::arg("secondary_gravitational_parameter"),
          py::arg("primary_secondary_distance"),
          py::arg("inertial_cartesian_state"),
          py::arg("unnormaized_time") = false,
          get_docstring("inertial_cartesian_to_corotating_normalized").c_str() );

    m.def("compute_cr3bp_state_derivative",
          &tp::computeCr3bpStateDerivative,
          py::arg("dimensionless_time"),
          py::arg("normalized_corotating_state"),
          py::arg("mass_parameter"),
          get_docstring("compute_cr3bp_state_derivative").c_str() );


    m.def("propagate_cr3bp",
          &tp::performCR3BPIntegration,
          py::arg("integrator_settings"),
          py::arg("mass_parameter"),
          py::arg("initial_normalized_corotating_state"),
          py::arg("final_dimensionless_time"),
          py::arg("propagate_to_exact_final_time") = false,
          get_docstring("propagate_cr3bp").c_str() );

    m.def("propagate_clohessy_wiltshire",
          &tba::propagateClohessyWiltshire,
          py::arg("initial_state"),
          py::arg("propagation_duration"),
          py::arg("central_body_graviational_parameter"),
          py::arg("reference_orbital_radius"),
          get_docstring("propagate_clohessy_wiltshire").c_str() );


    m.def("create_cr3bp_equivalent_body_settings",
          &tp::setupBodySettingsCR3BP,
          py::arg("primary_secondary_distance"),
          py::arg("primary_name"),
          py::arg("secondary_name"),
          py::arg("primary_gravitational_parameter"),
          py::arg("secondary_gravitational_parameter"),
          py::arg("frame_orientation"),
          get_docstring("create_cr3bp_equivalent_body_settings").c_str( ) );



}

} // namespace approximate_models
} // namespace trajectory_design
} // namespace tudatpy
