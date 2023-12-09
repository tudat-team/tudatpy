/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_astro.h"

#include "expose_astro/expose_gravitation.h"
#include "expose_astro/expose_element_conversion.h"
#include "expose_astro/expose_frame_conversion.h"
#include "expose_astro/expose_time_conversion.h"
#include "expose_astro/expose_two_body_dynamics.h"
#include "expose_astro/expose_fundamentals.h"
#include "expose_astro/expose_polyhedron_utilities.h"
#include "expose_astro/expose_timing_system.h"

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace tudatpy {

namespace astro {

void expose_astro(py::module &m) {

    auto element_conversion = m.def_submodule("element_conversion");
    element_conversion::expose_element_conversion(element_conversion);

    auto frame_conversion = m.def_submodule("frame_conversion");
    frame_conversion::expose_frame_conversion(frame_conversion);

    auto time_conversion = m.def_submodule("time_conversion");
    time_conversion::expose_time_conversion(time_conversion);

    auto two_body_dynamics = m.def_submodule("two_body_dynamics");
    two_body_dynamics::expose_two_body_dynamics(two_body_dynamics);

    auto gravitation = m.def_submodule("gravitation");
    gravitation::expose_gravitation(gravitation);

    auto fundamentals = m.def_submodule("fundamentals");
    fundamentals::expose_fundamentals(fundamentals);

    auto polyhedron_utilities = m.def_submodule("polyhedron_utilities");
    polyhedron_utilities::expose_polyhedron_utilities( polyhedron_utilities );

    auto timing_system = m.def_submodule("timing_system");
    timing_system::expose_timing_system( timing_system );

}

} // namespace astro

}// namespace tudatpy
