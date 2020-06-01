//
// Created by ggarrett on 23-04-20.
//

#ifndef TUDATBUNDLE_EXPOSE_SIMULATION_SETUP_H
#define TUDATBUNDLE_EXPOSE_SIMULATION_SETUP_H

#include "tudat/simulation/simulation.h"
#include <pybind11/pybind11.h>

// Ephemerides.
#include "tudat/math/interpolators/createInterpolator.h"
#include "tudat/math/interpolators/lagrangeInterpolator.h"
#include "tudat/astro/ephemerides/constantEphemeris.h"

#ifndef EPHEMERIS_H
#include "tudat/astro/ephemerides/ephemeris.h"
#define EPHEMERIS_H
#endif


namespace py = pybind11;

namespace tudatpy {
    void expose_simulation_setup(py::module &m);
}

#endif //TUDATBUNDLE_EXPOSE_SIMULATION_SETUP_H

