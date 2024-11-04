/*    Copyright (c) 2010-2018, Delft University of Technology
  *    All rights reserved
  *
  *    This file is part of the Tudat. Redistribution and use in source and
  *    binary forms, with or without modification, are permitted exclusively
  *    under the terms of the Modified BSD license. You should have received
  *    a copy of the license with this file. If not, please or visit:
  *    http://tudat.tudelft.nl/LICENSE.
  */

#include "expose_statistics.h"

#include "docstrings.h"

#include <tudat/basics/basicTypedefs.h>

#include "tudat/astro/system_models/timingSystem.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

namespace py = pybind11;

namespace ts = tudat::statistics;
namespace tsm = tudat::system_models;

namespace tudatpy {

    void expose_statistics(py::module &m) {

        m.def("calculate_allan_variance_of_dataset", &ts::calculateAllanVarianceOfTimeDataSet,
              py::arg( "timing_errors"),
              py::arg( "time_step_size"),
              get_docstring("calculate_allan_variance_of_dataset").c_str());

        m.def("convert_allan_variance_amplitudes_to_phase_noise_amplitudes",
              &tsm::convertAllanVarianceAmplitudesToPhaseNoiseAmplitudes,
              py::arg("allan_variance_amplitudes"),
              py::arg("frequency_domain_cutoff_frequency"),
              py::arg("is_inverse_square_term_flicker_phase_noise") = 0,
              get_docstring("convert_allan_variance_amplitudes_to_phase_noise_amplitudes").c_str());

#if( TUDAT_BUILD_WITH_FFTW3 )
        m.def("generate_noise_from_allan_deviation",
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
#endif


    };

}// namespace tudatpy
