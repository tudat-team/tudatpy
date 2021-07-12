//
// Created by Filippo Oggionni on 12/07/21.
//

#ifndef TUDATBUNDLE_EXPOSE_TORQUE_SETUP_H
#define TUDATBUNDLE_EXPOSE_TORQUE_SETUP_H

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "tudat/simulation/propagation_setup.h"

namespace py = pybind11;

namespace tp = tudat::propagators;

namespace tudatpy {

    void expose_torque_setup(py::module &m);


}// namespace tudatpy
#endif //TUDATBUNDLE_EXPOSE_TORQUE_SETUP_H
