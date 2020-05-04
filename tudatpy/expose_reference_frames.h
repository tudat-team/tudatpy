//
// Created by ggarrett on 02-05-20.
//

#ifndef TUDATBUNDLE_EXPOSE_REFERENCE_FRAMES_H
#define TUDATBUNDLE_EXPOSE_REFERENCE_FRAMES_H


#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include "Tudat/Astrodynamics/ReferenceFrames/aerodynamicAngleCalculator.h"

namespace py = pybind11;

namespace tudatpy {

    void expose_reference_frames(py::module &m);

}


#endif //TUDATBUNDLE_EXPOSE_REFERENCE_FRAMES_H
