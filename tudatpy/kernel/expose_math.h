/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the TudatPy. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDATPY_EXPOSE_MATH_H
#define TUDATPY_EXPOSE_MATH_H

namespace py = pybind11;

namespace tudatpy {

void expose_math(py::module &m);

};// namespace tudatpy

#endif// TUDATPY_EXPOSE_MATH_H
