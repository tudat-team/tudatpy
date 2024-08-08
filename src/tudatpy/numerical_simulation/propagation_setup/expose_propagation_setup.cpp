/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <tudat/simulation/propagation_setup.h>

#include "tudatpy/docstrings.h"

namespace py = pybind11;
namespace tba = tudat::basic_astrodynamics;
namespace tss = tudat::simulation_setup;
namespace tp = tudat::propagators;
namespace tinterp = tudat::interpolators;
namespace te = tudat::ephemerides;
namespace tni = tudat::numerical_integrators;
namespace trf = tudat::reference_frames;
namespace tmrf = tudat::root_finders;

namespace tudatpy {
    namespace numerical_simulation {
        namespace propagation_setup {

            PYBIND11_MODULE(expose_propagation_setup, m) {
                m.def("create_acceleration_models",
                      py::overload_cast<const tss::SystemOfBodies &,
                                        const tss::SelectedAccelerationMap &,
                                        const std::vector<std::string> &,
                                        const std::vector<std::string> &>(
                          &tss::createAccelerationModelsMap),
                      py::arg("body_system"),
                      py::arg("selected_acceleration_per_body"),
                      py::arg("bodies_to_propagate"), py::arg("central_bodies"),
                      get_docstring("create_acceleration_models").c_str());

                m.def("create_torque_models", &tss::createTorqueModelsMap,
                      py::arg("body_system"),
                      py::arg("selected_torque_per_body"),
                      py::arg("bodies_to_propagate"),
                      get_docstring("create_torque_models").c_str());

                m.def("create_mass_rate_models", &tss::createMassRateModelsMap,
                      py::arg("body_system"),
                      py::arg("selected_mass_rates_per_body"),
                      py::arg("acceleration_models") = nullptr,
                      get_docstring("create_mass_rate_models").c_str());
            }
        }  // namespace propagation_setup
    }  // namespace numerical_simulation
}  // namespace tudatpy
