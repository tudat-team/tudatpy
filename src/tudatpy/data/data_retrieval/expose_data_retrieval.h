/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDATPY_EXPOSE_DATA_RETRIEVAL_H
#define TUDATPY_EXPOSE_DATA_RETRIEVAL_H

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace tudatpy
{
namespace data
{
namespace data_retrieval
{

void expose_data_retrieval( py::module& m );

}  // namespace data_retrieval
}  // namespace data
}  // namespace tudatpy

#endif  // TUDATPY_EXPOSE_DATA_RETRIEVAL_H
