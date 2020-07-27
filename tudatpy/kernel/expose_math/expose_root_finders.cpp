/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_root_finders.h"

#include <tudat/math/root_finders.h>

#include <pybind11/pybind11.h>

namespace py = pybind11;
namespace trf = tudat::root_finders;

namespace tudatpy {

void expose_root_finders(py::module &m) {

  /*
   *
   *
   *  root_finders
   *  ├── bisection.h
   *  ├── createRootFinder.h
   *  ├── halleyRootFinder.h
   *  ├── newtonRaphson.h
   *  ├── rootFinder.h
   *  ├── secantRootFinder.h
   *  └── terminationConditions.h
   *
   */

  py::class_<trf::RootFinderCore<double>,
             std::shared_ptr<trf::RootFinderCore<double>>>
      RootFinderCore_(m, "RootFinderCore");

  py::class_<trf::NewtonRaphsonCore<double>,
             std::shared_ptr<trf::NewtonRaphsonCore<double>>,
             trf::RootFinderCore<double>>(m, "NewtonRaphsonCore")
      .def(
          py::init<const double, const unsigned int>(),
          py::arg("x_tol"),
          py::arg("max_iter"));

  py::class_<trf::RootFinderSettings,
             std::shared_ptr<trf::RootFinderSettings>>
      RootFinderSettings_(m, "RootFinderSettings");
}

}// namespace tudatpy
