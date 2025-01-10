/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDATPY_EXPOSE_GRAVITATION_H
#define TUDATPY_EXPOSE_GRAVITATION_H

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace tudatpy
{

namespace astro
{
namespace gravitation
{

void expose_gravitation( py::module &m );

}  // namespace gravitation
}  // namespace astro
}  // namespace tudatpy

#endif  // TUDATPY_EXPOSE_GRAVITATION_H
