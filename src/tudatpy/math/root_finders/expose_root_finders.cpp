/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <pybind11/pybind11.h>
#include <tudat/math/root_finders.h>

namespace py = pybind11;
namespace trf = tudat::root_finders;

namespace tudatpy {
    namespace math {

        namespace root_finders {

            PYBIND11_MODULE(expose_root_finders, m) {
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

                py::enum_<trf::MaximumIterationHandling>(
                    m, "MaximumIterationHandling",
"")
                    .value(
                        "accept_result",
                        trf::MaximumIterationHandling::accept_result,
"")
                    .value("accept_result_with_warning",
                           trf::MaximumIterationHandling::
                               accept_result_with_warning,
"")
                    .value("throw_exception",
                           trf::MaximumIterationHandling::throw_exception,
"")
                    .export_values();


                py::class_<trf::RootFinder<double>,
                           std::shared_ptr<trf::RootFinder<double>>>(
                    m, "RootFinderCore");

                py::class_<trf::NewtonRaphson<double>,
                           std::shared_ptr<trf::NewtonRaphson<double>>,
                           trf::RootFinder<double>>(m, "NewtonRaphsonCore")
                    .def(py::init<const double, const unsigned int>(),
                         py::arg("x_tol"), py::arg("max_iter"));

                py::class_<trf::RootFinderSettings,
                           std::shared_ptr<trf::RootFinderSettings>>(
                    m, "RootFinderSettings",
"");

                m.def("bisection", &trf::bisectionRootFinderSettings,
                      py::arg("relative_variable_tolerance") = TUDAT_NAN,
                      py::arg("absolute_variable_tolerance") = TUDAT_NAN,
                      py::arg("root_function_tolerance") = TUDAT_NAN,
                      py::arg("maximum_iteration") = 1000,
                      py::arg("maximum_iteration_handling") =
                          trf::throw_exception,
"");


                m.def("newton_raphson", &trf::newtonRaphsonRootFinderSettings,
                      py::arg("relative_variable_tolerance") = TUDAT_NAN,
                      py::arg("absolute_variable_tolerance") = TUDAT_NAN,
                      py::arg("root_function_tolerance") = TUDAT_NAN,
                      py::arg("maximum_iteration") = 1000,
                      py::arg("maximum_iteration_handling") =
                          trf::throw_exception,
"");

                m.def("halley", &trf::halleyRootFinderSettings,
                      py::arg("relative_variable_tolerance") = TUDAT_NAN,
                      py::arg("absolute_variable_tolerance") = TUDAT_NAN,
                      py::arg("root_function_tolerance") = TUDAT_NAN,
                      py::arg("maximum_iteration") = 1000,
                      py::arg("maximum_iteration_handling") =
                          trf::throw_exception,
"");

                m.def("secant", &trf::secantRootFinderSettings,
                      py::arg("relative_variable_tolerance") = TUDAT_NAN,
                      py::arg("absolute_variable_tolerance") = TUDAT_NAN,
                      py::arg("root_function_tolerance") = TUDAT_NAN,
                      py::arg("maximum_iteration") = 1000,
                      py::arg("maximum_iteration_handling") =
                          trf::throw_exception,
"");
            }

        }  // namespace root_finders
    }  // namespace math

}  // namespace tudatpy
