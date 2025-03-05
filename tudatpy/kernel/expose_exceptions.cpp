/*    Copyright (c) 2010-2020, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "expose_exceptions.h"
#include <tudat/basics/tudatExceptions.h>

namespace py = pybind11;
namespace te = tudat::exceptions;

namespace tudatpy
{
namespace exceptions
{

void expose_exceptions( py::module &m )
{
    py::register_exception< te::TudatError >( m, "TudatError", PyExc_RuntimeError ).doc( ) = R"(
        Base Error thrown by Tudat.
    )";
}

}  // namespace exceptions
}  // namespace tudatpy
