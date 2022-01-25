/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_trajectory_design.h"

#include "expose_trajectory_design/expose_transfer_trajectory.h"
#include "expose_trajectory_design/expose_shape_based_thrust.h"

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace tudatpy {
namespace trajectory_design{

void expose_trajectory_design(py::module &m) {


  auto shape_based_thrust_submodule = m.def_submodule("shape_based_thrust");
  expose_shape_based_thrust(shape_based_thrust_submodule);

  auto transfer_trajectory = m.def_submodule("transfer_trajectory");
  transfer_trajectory::expose_transfer_trajectory(transfer_trajectory);


}
}

}// namespace tudatpy
