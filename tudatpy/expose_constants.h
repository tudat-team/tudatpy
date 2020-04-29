//
// Created by ggarrett on 29-04-20.
//

#ifndef TUDATBUNDLE_EXPOSE_CONSTANTS_H
#define TUDATBUNDLE_EXPOSE_CONSTANTS_H

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>

#include "Tudat/Astrodynamics/BasicAstrodynamics/celestialBodyConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

namespace py = pybind11;

namespace tudatpy {
    void expose_constants(py::module &m);
}

#endif //TUDATBUNDLE_EXPOSE_CONSTANTS_H
