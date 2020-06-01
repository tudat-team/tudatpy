//
// Created by ggarrett on 24-04-20.
//

#include <pybind11/pybind11.h>
#include "expose_basic_astrodynamics.h"

#include "tudat/astro/ephemerides/constantEphemeris.h"

// Ephemerides.
#include "tudat/math/interpolators/createInterpolator.h"
#include "tudat/math/interpolators/lagrangeInterpolator.h"
#include "tudat/astro/ephemerides/constantEphemeris.h"
#include "tudat/astro/ephemerides/ephemeris.h"

#ifndef EPHEMERIS_H
#include "tudat/astro/ephemerides/ephemeris.h"
#define EPHEMERIS_H
#endif

namespace py = pybind11;
namespace tss = tudat::simulation_setup;
namespace te = tudat::ephemerides;

void expose_ephemeris(py::module &m) {

    py::class_<te::Ephemeris>(m, "Ephemeris")
            .def(py::init<const std::string &, const std::string &>(),
                 py::arg("referenceFrameOrigin") = "",
                 py::arg("referenceFrameOrigin") = "");

    py::enum_<tss::EphemerisType>(te::Ephemeris, "type")
            .value("approximate_planet_positions", tss::approximate_planet_positions)
            .value("direct_spice_ephemeris", tss::direct_spice_ephemeris)
            .value("interpolated_spice", tss::interpolated_spice)
            .value("constant_ephemeris", tss::constant_ephemeris)
            .value("kepler_ephemeris", tss::kepler_ephemeris)
            .value("custom_ephemeris", tss::custom_ephemeris);


};

