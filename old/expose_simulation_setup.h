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

// Acceleration model types
#include "tudat/simulation/propagation/createAccelerationModels.h"
#include "tudat/astro/basic_astro/accelerationModelTypes.h"
#include "tudat/astro/basic_astro/accelerationModel.h"

#include "tudat/astro/ephemerides/ephemeris.h"
#include "tudat/astro/basic_astro/stateVectorIndices.h"
#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/astro/gravitation/gravityFieldModel.h"
#include "tudat/simulation/propagation/propagationSettings.h"

namespace py = pybind11;

namespace tudatpy {
    void expose_simulation_setup(py::module &m);
}

#endif //TUDATBUNDLE_EXPOSE_SIMULATION_SETUP_H

