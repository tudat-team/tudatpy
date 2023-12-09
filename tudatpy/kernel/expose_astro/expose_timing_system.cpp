/*    Copyright (c) 2010-2018, Delft University of Technology
  *    All rights reserved
  *
  *    This file is part of the Tudat. Redistribution and use in source and
  *    binary forms, with or without modification, are permitted exclusively
  *    under the terms of the Modified BSD license. You should have received
  *    a copy of the license with this file. If not, please or visit:
  *    http://tudat.tudelft.nl/LICENSE.
  */

#include "tudatpy/docstrings.h"
#include "tudatpy/scalarTypes.h"

#include "expose_timing_system.h"

#include "kernel/expose_numerical_simulation.h"

#include <tudat/astro/system_models/timingSystem.h>
#include <tudat/astro/system_models/vehicleSystems.h>

#include "tudat/basics/timeType.h"

#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>  // to enable conversions from/to C++ standard library types
#include <pybind11/operators.h>


namespace py = pybind11;
namespace tsm = tudat::system_models;
namespace ts = tudat::statistics;

namespace tudatpy {

    namespace astro {
        namespace timing_system {

            void expose_timing_system(py::module &m) {

                m.def("convert_allan_variance_amplitudes_to_phase_noise_amplitudes",
                      &tsm::convertAllanVarianceAmplitudesToPhaseNoiseAmplitudes,
                      py::arg("allan_variance_amplitudes"),
                      py::arg("frequency_domain_cutoff_frequency"),
                      py::arg("is_inverse_square_term_flicker_phase_noise") = 0,
                      get_docstring("convert_allan_variance_amplitudes_to_phase_noise_amplitudes").c_str()
                );

                m.def("generate_clock_noise",
                      &tsm::generateClockNoise,
                      py::arg("allan_variance_amplitudes"),
                      py::arg("start_time"),
                      py::arg("end_time"),
                      py::arg("number_of_time_steps"),
                      py::arg("is_inverse_square_term_flicker_phase_noise") = 0,
                      py::arg("seed") = ts::defaultRandomSeedGenerator->getRandomVariableValue(),
                      get_docstring("generate_clock_noise").c_str()
                );

                m.def("generate_colored_clock_noise",
                      &tsm::generateColoredClockNoise,
                      py::arg("allan_variance_amplitudes"),
                      py::arg("variance_type"),
                      py::arg("start_time"),
                      py::arg("end_time"),
                      py::arg("number_of_time_steps"),
                      py::arg("seed") = ts::defaultRandomSeedGenerator->getRandomVariableValue(),
                      get_docstring("generate_clock_noise").c_str()
                );

                m.def("get_clock_noise_interpolator",
                      &tsm::getClockNoiseInterpolator,
                      py::arg("allan_variance_amplitudes"),
                      py::arg("start_time"),
                      py::arg("end_time"),
                      py::arg("time_step"),
                      py::arg("is_inverse_square_term_flicker_phase_noise") = 0,
                      py::arg("seed") = ts::defaultRandomSeedGenerator->getRandomVariableValue(),
                      get_docstring("get_clock_noise_interpolator").c_str()
                );
                m.def("get_colored_clock_noise_interpolator",
                      &tsm::getColoredClockNoiseInterpolator,
                      py::arg("allan_variance_nodes"),
                      py::arg("variance_type"),
                      py::arg("start_time"),
                      py::arg("end_time"),
                      py::arg("time_step"),
                      py::arg("seed") = ts::defaultRandomSeedGenerator->getRandomVariableValue(),
                      get_docstring("get_clock_noise_interpolator").c_str()
                );

                py::class_<tsm::TimingSystem,
                        std::shared_ptr<tsm::TimingSystem>>
                        (m, "TimingSystem",
                         get_docstring("TimingSystem").c_str())

                        .def( // ctor 1
                                py::init<
                                        const std::vector<tudat::Time>,
                                        const std::vector<double>,
                                        const std::function<std::function<double(const double)>(const double, const double,
                                                                                                const double)>,
                                        const double>(),
                                py::arg("arc_times"),
                                py::arg("all_arcs_polynomial_drift_coefficients") = std::vector<double>(),
                                py::arg("clock_noise_generation_function") = nullptr,
                                py::arg("clock_noise_time_step") = 1.0E-3)

                        .def( // ctor 2
                                py::init<
                                        const std::vector<tudat::Time>,
                                        const std::vector<std::vector<double> >,
                                        const std::function<std::function<double(const double)>(const double, const double,
                                                                                                const double)>,
                                        const double>(),
                                py::arg("arc_times"),
                                py::arg("polynomial_drift_coefficients"),
                                py::arg("clock_noise_generation_function") = nullptr,
                                py::arg("clock_noise_time_step") = 1.0E-3)
                        .def( // ctor 3
                                py::init<
                                        const std::vector<std::vector<double> >,
                                        const std::vector<std::function<double(const double)>>,
                                        const std::vector<tudat::Time> >(),
                                py::arg("polynomial_drift_coefficients"),
                                py::arg("stochastic_clock_noise_functions"),
                                py::arg("arc_times"));




            }
        } // namespace timing_system
    } // namespace astro
}// namespace tudatpy
