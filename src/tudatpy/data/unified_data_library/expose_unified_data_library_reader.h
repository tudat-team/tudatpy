/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDATPY_EXPOSE_UNIFIED_DATA_LIBRARY_READER_H
#define TUDATPY_EXPOSE_UNIFIED_DATA_LIBRARY_READER_H

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace tudatpy
{
namespace data
{
namespace unified_data_library
{

void expose_unified_data_library_reader( py::module& m );

}  // namespace unified_data_library
}  // namespace data
}  // namespace tudatpy

#endif  // TUDATPY_EXPOSE_UNIFIED_DATA_LIBRARY_READER_H