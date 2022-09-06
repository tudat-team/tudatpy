/*    Copyright (c) 2010-2020, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_data.h"
#include <tudat/utils/data/downloadFile.h>

namespace py = pybind11;

#include <pybind11/stl.h>

using namespace tudat::utils::data;

namespace tudatpy {
namespace utils {
namespace data {

void expose_data(py::module &m) {

    m.def("download_file",
          &download_file,
          py::arg("remote_url"),
          py::arg("cache") = "true",
          py::arg("verbosity ") = 1,
          py::arg("try_unzip ") = true,
          py::arg("prefix ") = py::none());
}

} // namespace data
} // namespace utils
} // namespace tudatpy