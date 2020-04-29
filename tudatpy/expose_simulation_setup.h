//
// Created by ggarrett on 23-04-20.
//

#ifndef TUDATBUNDLE_EXPOSE_SIMULATION_SETUP_H
#define TUDATBUNDLE_EXPOSE_SIMULATION_SETUP_H

#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
#include <pybind11/pybind11.h>

// Ephemerides.
#include "Tudat/Mathematics/Interpolators/createInterpolator.h"
#include "Tudat/Mathematics/Interpolators/lagrangeInterpolator.h"
#include "Tudat/Astrodynamics/Ephemerides/constantEphemeris.h"

#ifndef EPHEMERIS_H
#include "Tudat/Astrodynamics/Ephemerides/ephemeris.h"
#define EPHEMERIS_H
#endif


namespace py = pybind11;

namespace tudatpy {
    void expose_simulation_setup(py::module &m);
}

#endif //TUDATBUNDLE_EXPOSE_SIMULATION_SETUP_H

