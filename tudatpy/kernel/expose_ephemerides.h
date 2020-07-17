//
// Created by ggarrett on 02-05-20.
//

#ifndef TUDATBUNDLE_EXPOSE_EPHEMERIDES_H
#define TUDATBUNDLE_EXPOSE_EPHEMERIDES_H

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include "tudat/astro/ephemerides/rotationalEphemeris.h"
#include "tudat/astro/ephemerides/approximatePlanetPositions.h"

namespace py = pybind11;

namespace tudatpy {
    void expose_ephemerides(py::module &m);
}

#endif //TUDATBUNDLE_EXPOSE_EPHEMERIDES_H
