/*    Copyright (c) 2010-2020, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_utils/expose_data.h"

namespace py = pybind11;

namespace tudatpy {
namespace utils {

void expose_utils(py::module &m) {

  auto data = m.def_submodule("data");
  data::expose_data(data);
}

} // namespace utils
} // namespace tudatpy