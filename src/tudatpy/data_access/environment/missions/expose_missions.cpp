/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif
#include "expose_missions.h"

#include "grail/expose_grail.h"

#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
// #include <pybind11/native_enum.h>
#include <tudat/io/basicInputOutput.h>
#include <string>
#include <vector>

namespace py = pybind11;
namespace tio = tudat::input_output;

namespace tudatpy
{

namespace data_access
{
namespace environment
{
namespace missions
{

void expose_missions( py::module& m )
{
    auto grail_submodule = m.def_submodule( "grail" );
    grail::expose_grail( grail_submodule );
};

}  // namespace missions
}  // namespace environment
}  // namespace data_access
}  // namespace tudatpy
