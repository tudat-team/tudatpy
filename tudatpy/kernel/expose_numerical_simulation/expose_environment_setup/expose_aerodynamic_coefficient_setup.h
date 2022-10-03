/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#ifndef TUDATPY_EXPOSE_AERODYNAMIC_COEFFICIENT_SETUP_H
#define TUDATPY_EXPOSE_AERODYNAMIC_COEFFICIENT_SETUP_H

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace tudatpy {
namespace numerical_simulation {
namespace environment_setup {
namespace aerodynamic_coefficients {

    void expose_aerodynamic_coefficient_setup(py::module &m);

}// namespace aerodynamic_coefficients
}// namespace environment_setup
}// namespace numerical_simulation
}// namespace tudatpy

#endif //TUDATPY_EXPOSE_AERODYNAMIC_COEFFICIENT_SETUP_H
