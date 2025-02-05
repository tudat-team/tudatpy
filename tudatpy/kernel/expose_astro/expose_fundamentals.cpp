/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_fundamentals.h"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <tudat/astro/basic_astro.h>

namespace py = pybind11;
namespace tmg = tudat::mission_geometry;

namespace tudatpy
{
namespace astro
{
namespace fundamentals
{

void expose_fundamentals( py::module &m )
{
    m.def( "compute_shadow_function",
           &tmg::computeShadowFunction,
           py::arg( "occulted_body_position" ),
           py::arg( "occulted_body_radius" ),
           py::arg( "occulting_body_position" ),
           py::arg( "occulting_body_radius" ),
           py::arg( "satellite_position" ),
           R"doc(

Compute the shadow function.

Returns the value of of the shadow function. Returns 0 if the satellite is in umbra, 1 if the
satellite is fully exposed and a value between 0 and 1 if the satellite is in penumbra or antumbra.

The point of view is from the satellite. The occulting body (for example the Earth) is the body
that blocks the light from the occulted body (for example the Sun).

Reference: Section 3.4 from ( Montebruck O, Gill E., 2005) and Fig. 5 from (Zhang et al., 2019).

Parameters
----------
occulted_body_position : numpy.ndarray
    Vector containing Cartesian coordinates of the occulted body.
occulted_body_radius : float
    Mean radius of occulted body.
occulting_body_position : numpy.ndarray
    Vector containing Cartesian coordinates of the occulting body.
occulting_body_radius : float
    Mean radius of occulting body.
satellite_position : numpy.ndarray
    Vector containing Cartesian coordinates of the satellite.
Returns
-------
float
    Shadow function value

    )doc" );
}

}  // namespace fundamentals
}  // namespace astro
}  // namespace tudatpy
