/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDATPY_EXPOSE_SP3_H
#define TUDATPY_EXPOSE_SP3_H

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace tudatpy
{
namespace data_input
{
namespace environment_data
{
namespace sp3
{

//! Expose the SP3 file reader and parsed-file container in the data-input SP3 submodule.
void expose_sp3( py::module& m );

}  // namespace sp3
}  // namespace environment_data
}  // namespace data_input
}  // namespace tudatpy

#endif  // TUDATPY_EXPOSE_SP3_H
