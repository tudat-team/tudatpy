
/*    Copyright (c) 2010-2018, Delft University of Technology
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
#include "expose_astro.h"

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

void expose_astro( py::module& m )
{
    auto element_conversion = m.def_submodule( "element_conversion" );
    element_conversion::expose_element_conversion( element_conversion );

    auto frame_conversion = m.def_submodule( "frame_conversion" );
    frame_conversion::expose_frame_conversion( frame_conversion );

    auto time_representation = m.def_submodule( "time_representation" );
    time_representation::expose_time_representation( time_representation );

    auto two_body_dynamics = m.def_submodule( "two_body_dynamics" );
    two_body_dynamics::expose_two_body_dynamics( two_body_dynamics );

    auto gravitation = m.def_submodule( "gravitation" );
    gravitation::expose_gravitation( gravitation );

    auto fundamentals = m.def_submodule( "fundamentals" );
    fundamentals::expose_fundamentals( fundamentals );

    auto polyhedron_utilities = m.def_submodule( "polyhedron_utilities" );
    polyhedron_utilities::expose_polyhedron_utilities( polyhedron_utilities );
}

}  // namespace astro

}  // namespace tudatpy
