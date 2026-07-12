/*    Copyright (c) 2010-2019, Delft University of Technology
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
#include "expose_media_corrections.h"

namespace py = pybind11;

namespace tudatpy
{
namespace data
{
namespace data_retrieval
{
namespace media_corrections
{

void expose_media_corrections( py::module& m ) {}

}  // namespace media_corrections
}  // namespace data_retrieval
}  // namespace data
}  // namespace tudatpy
