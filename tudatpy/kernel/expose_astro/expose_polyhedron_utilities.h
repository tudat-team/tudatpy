/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDATBUNDLE_EXPOSE_POLYHEDRON_UTILITIES_H
#define TUDATBUNDLE_EXPOSE_POLYHEDRON_UTILITIES_H

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace tudatpy {
namespace astro {
namespace polyhedron_utilities {

void expose_polyhedron_utilities(py::module &m);

} // namespace polyhedron_utilities
} // namespace astro
} // namespace tudatpy


#endif //TUDATBUNDLE_EXPOSE_POLYHEDRON_UTILITIES_H
