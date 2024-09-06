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
#include <tudat/math/integrators.h>

namespace tni = tudat::numerical_integrators;
namespace py = pybind11;

// typedef std::function<
//     Eigen::VectorXd(
//         const double,
//         const Eigen::VectorXd &)>
//     StateDerivativeFunction;

typedef std::function<Eigen::Matrix<double, Eigen::Dynamic, 1>(
    const double, const Eigen::Matrix<double, Eigen::Dynamic, 1>&)>
    stateDerivativeFunction;

namespace tudatpy {

    PYBIND11_MODULE(expose_numerical_integrators, m) {
        //  py::class_<tni::NumericalIntegrator<>>(m, "NumericalIntegrator");
        //  //      .def(py::init<>);
        //  //

        //  //
        ////  py::class_<tni::RungeKutta4Integrator < double, Eigen::VectorXd,
        ///Eigen::VectorXd >,//      tni::NumericalIntegrator<double,
        ///Eigen::VectorXd, Eigen::VectorXd, double> /
        ///std::shared_ptr<tni::RungeKutta4Integrator < double, Eigen::VectorXd,
        ///Eigen::VectorXd >> /             >(m, "RungeKutta4Integrator") /
        ///.def(py::init< /           const stateDerivativeFunction &, / const
        ///double, /           const double &>());

        //  // Alias
        //  //  m.def("rk4", m.attr("RungeKutta4Integrator"))
    }

}  // namespace tudatpy
