/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_{{ submodule_name }}.h"

#include <string>

#include <pybind11/pybind11.h>

namespace {{ project }} {

std::string hello_world(){

    return "Hello World!";

};

void expose_{{ module_name.lower() }}(py::module &m){

    m.def("hello_world", &hello_world);

};

}
