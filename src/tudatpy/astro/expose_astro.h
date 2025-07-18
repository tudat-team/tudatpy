/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDATPY_EXPOSE_ASTRO_H
#define TUDATPY_EXPOSE_ASTRO_H

#include <pybind11/pybind11.h>

#include "element_conversion/expose_element_conversion.h"
#include "frame_conversion/expose_frame_conversion.h"
#include "fundamentals/expose_fundamentals.h"
#include "gravitation/expose_gravitation.h"
#include "polyhedron_utilities/expose_polyhedron_utilities.h"
#include "time_representation/expose_time_representation.h"
#include "two_body_dynamics/expose_two_body_dynamics.h"

namespace py = pybind11;

namespace tudatpy
{
namespace astro
{

void expose_astro( py::module &m );

}  // namespace astro
}  // namespace tudatpy

#endif  // TUDATPY_EXPOSE_ASTRO_H
