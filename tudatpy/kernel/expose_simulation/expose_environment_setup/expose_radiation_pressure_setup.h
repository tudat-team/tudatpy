/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDATBUNDLE_EXPOSE_RADIATION_PRESSURE_SETUP_H
#define TUDATBUNDLE_EXPOSE_RADIATION_PRESSURE_SETUP_H

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace tudatpy {
namespace simulation {
namespace environment_setup {

    void expose_radiation_pressure_setup(py::module &m);

}// namespace environment_setup
}// namespace simulation
}// namespace tudatpy

#endif //TUDATBUNDLE_EXPOSE_RADIATION_PRESSURE_SETUP_H
