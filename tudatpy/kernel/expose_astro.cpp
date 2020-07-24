/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_astro/expose_aerodynamics.h"
#include "expose_astro/expose_conversion.h"
#include "expose_astro/expose_ephemerides.h"
#include "expose_astro/expose_fundamentals.h"
#include "expose_astro/expose_gravitation.h"
#include "expose_astro/expose_propagators.h"
#include "expose_astro/expose_reference_frames.h"

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace tudatpy {

void expose_astro(py::module &m) {

  auto fundamentals = m.def_submodule("fundamentals");
  expose_fundamentals(fundamentals);

  auto conversion = m.def_submodule("conversion");
  expose_conversion(conversion);

  auto frames = m.def_submodule("frames");
  expose_frames(frames);

  auto aerodynamics = m.def_submodule("aerodynamics");
  expose_aerodynamics(aerodynamics);

  auto ephemerides = m.def_submodule("ephemerides");
  expose_ephemerides(ephemerides);

  auto gravitation = m.def_submodule("gravitation");
  expose_gravitation(gravitation);

  auto propagators = m.def_submodule("propagators");
  expose_propagators(propagators);
};

};// namespace tudatpy
