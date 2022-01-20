/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef {{ module_name | upper }}_H
#define {{ module_name | upper }}_H

#include <string>

#include <pybind11/pybind11.h>

namespace {{ project_name }} {

void expose_{{ module_name | lower }}(py::module &m);

}

#endif//{{ module_name | upper }}_H
