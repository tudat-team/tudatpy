/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDATPY_EXPOSE_CLOSE_ENCOUNTERS_H
#define TUDATPY_EXPOSE_CLOSE_ENCOUNTERS_H

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace tudatpy
{
namespace astro
{
namespace close_encounters
{

void expose_close_encounters( py::module& m );

}  // namespace close_encounters
}  // namespace astro
}  // namespace tudatpy

#endif  // TUDATPY_EXPOSE_CLOSE_ENCOUNTERS_H
