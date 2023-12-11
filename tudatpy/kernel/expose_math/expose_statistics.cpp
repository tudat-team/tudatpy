/*    Copyright (c) 2010-2018, Delft University of Technology
  *    All rights reserved
  *
  *    This file is part of the Tudat. Redistribution and use in source and
  *    binary forms, with or without modification, are permitted exclusively
  *    under the terms of the Modified BSD license. You should have received
  *    a copy of the license with this file. If not, please or visit:
  *    http://tudat.tudelft.nl/LICENSE.
  */

#include "expose_statistics.h"

#include "tudatpy/docstrings.h"

#include <tudat/basics/basicTypedefs.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

namespace py = pybind11;

namespace ts = tudat::statistics;

namespace tudatpy {

    void expose_statistics(py::module &m) {

        m.def("calculate_allan_variance_of_dataset", &ts::calculateAllanVarianceOfTimeDataSet,
              py::arg( "timing_errors"),
              py::arg( "time_step_size"),
              get_docstring("calculate_allan_variance_of_dataset").c_str());


    };

}// namespace tudatpy
