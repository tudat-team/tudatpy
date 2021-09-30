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
#include "expose_astro/expose_two_body_dynamics.h"

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace tudatpy {
<<<<<<< HEAD
namespace astro {
=======
    namespace astro {
>>>>>>> 430786d34e5d83dfba61cf73c209ba1eca6b921e

        void expose_astro(py::module &m) {

<<<<<<< HEAD
  auto element_conversion = m.def_submodule("element_conversion");
  astro::element_conversion::expose_element_conversion(element_conversion);

  auto frame_conversion = m.def_submodule("frame_conversion");
  astro::frame_conversion::expose_frame_conversion(frame_conversion);
=======
            auto element_conversion = m.def_submodule("element_conversion");
            element_conversion::expose_element_conversion(element_conversion);

            auto frame_conversion = m.def_submodule("frame_conversion");
            frame_conversion::expose_frame_conversion(frame_conversion);
>>>>>>> 430786d34e5d83dfba61cf73c209ba1eca6b921e

            auto two_body_dynamics = m.def_submodule("two_body_dynamics");
            two_body_dynamics::expose_two_body_dynamics(two_body_dynamics);

<<<<<<< HEAD
  auto gravitation = m.def_submodule("gravitation");
  astro::gravitation::expose_gravitation(gravitation);

=======
            auto gravitation = m.def_submodule("gravitation");
            gravitation::expose_gravitation(gravitation);

            auto observations = m.def_submodule("observations");
            observations::expose_observations(observations);
>>>>>>> 430786d34e5d83dfba61cf73c209ba1eca6b921e

        }

    } // namespace astro

}
}// namespace tudatpy
