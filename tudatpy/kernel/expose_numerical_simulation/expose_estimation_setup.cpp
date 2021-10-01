/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_estimation_setup.h"
#include "expose_estimation_setup/expose_estimated_parameter_setup.h"
#include "expose_estimation_setup/expose_observations_setup.h"


#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>

namespace py = pybind11;

namespace tudatpy {
namespace numerical_simulation {
namespace estimation_setup {


void expose_estimation_setup(py::module &m) {

    auto parameter_setup = m.def_submodule("parameter");
    parameter::expose_estimated_parameter_setup(parameter_setup);

    auto observations_setup = m.def_submodule("observations");
    observation::expose_observations_setup(observations_setup);
}

}
}
}// namespace tudatpy
