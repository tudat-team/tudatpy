//
// Created by ggarrett on 29-04-20.
//

#ifndef TUDATBUNDLE_EXPOSE_CONSTANTS_H
#define TUDATBUNDLE_EXPOSE_CONSTANTS_H

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>

#include "tudat/astro/basic_astro/celestialBodyConstants.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/math/basic/mathematicalConstants.h"

namespace py = pybind11;

namespace tudatpy {
    void expose_constants(py::module &m);
}

#endif //TUDATBUNDLE_EXPOSE_CONSTANTS_H
