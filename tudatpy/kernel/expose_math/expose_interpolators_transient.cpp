/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_interpolators.h"

#include <pybind11/pybind11.h>

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <type_traits>
#include <unsupported/Eigen/CXX11/Tensor>

#include <variant>
#include <pybind11/stl.h>

namespace py = pybind11;

namespace ti = tudat::interpolators;

namespace tudatpy {

namespace prototype {

// TODO: Make recommendation to add () operator to OneDimensionalInterpolator which calls OneDimensionalInterpolator.interpolate(x)
// TODO: Add constructor for interpolate which interpolates a Eigen::VectorXd.

//template<typename DependentVariableType, typename IndependentVariableType, typename ReturnDataType>
//std::function<ReturnDataType(std::vector<IndependentVariableType>)> interp1d(std::map<IndependentVariableType, DependentVariableType> &stateMap,
//                                                                             const std::string &kind = "linear",
//                                                                             const int &interpolateOrder = 8) {
//
//  std::shared_ptr<ti::InterpolatorSettings> interpolatorSettings;
//  if (kind == "lagrange") {
//    interpolatorSettings = std::make_shared<ti::LagrangeInterpolatorSettings>(interpolateOrder);
//  } else if (kind == "linear") {
//    interpolatorSettings = std::make_shared<ti::InterpolatorSettings>(ti::linear_interpolator);
//  } else if (kind == "piecewise_constant") {
//    interpolatorSettings = std::make_shared<ti::InterpolatorSettings>(ti::piecewise_constant_interpolator);
//  } else if (kind == "cubic_spline") {
//    interpolatorSettings = std::make_shared<ti::InterpolatorSettings>(ti::cubic_spline_interpolator);
//  } else if (kind == "hermite_spline") {
//    interpolatorSettings = std::make_shared<ti::InterpolatorSettings>(ti::hermite_spline_interpolator);
//  } else {
//    throw std::invalid_argument("Argument for 'kind' is invalid.");
//  }
//  std::shared_ptr<ti::OneDimensionalInterpolator<double, Eigen::VectorXd>> interpolator =
//      ti::createOneDimensionalInterpolator(stateMap, interpolatorSettings);
//
//  if (std::is_same<ReturnDataType, std::map<DependentVariableType, IndependentVariableType>>::value){
//
//  } else if (std::is_same<ReturnDataType, Eigen::MatrixXd>::value) {
//
//  } else if (std::is_same<ReturnDataType, Eigen::Tensor<>>)
//
//  return [=](std::vector<IndependentVariableType> x) {
//    //  return [=](double x) {
//    Eigen::MatrixXd interpolated(stateMap.begin()->second.size(), x.size());
//    for (int i = 0; i < x.size(); ++i) { interpolated.row(i) = interpolator->interpolate(x(i)); }
//    return interpolated;
//  };
//}

//interp1d<DependentVariableType, typename IndependentVariableType, typename ReturnDataType>;
//interp1d<double, Eigen::MatrixXd, typename std::map<double, Eigen::MatrixXd>>;
//interp1d<double, Eigen::VectorXd, typename std::map<double, Eigen::VectorXd>>;

// overload y[m,n]
//std::function<Eigen::MatrixXd(Eigen::VectorXd)> interp1d(const Eigen::VectorXd &x,
//                                                         const Eigen::MatrixXd &y,
//                                                         const std::string &kind = "linear",
//                                                         const int axis = 0,
//                                                         const int &interpolateOrder = 8) {
//  std::map<double, Eigen::VectorXd> stateMap;
//  if (axis == 0) {
//    for (auto row = 0; row < y.rows(); ++row) { stateMap[x(row)] = y.row(row); }
//  } else if (axis == 1) {
//    for (auto col = 0; col < y.cols(); ++col) { stateMap[x(col)] = y.col(col).transpose(); }
//  } else if (axis == -1) {
//    for (auto col = 0; col < y.cols(); ++col) { stateMap[x(col)] = y.col(col).transpose(); }
//  } else {
//    throw std::invalid_argument("Provided 'axis' for interpolation along, is invalid for the provided y argument.");
//  }
//  return interp1d(stateMap, kind, interpolateOrder);
//}

}// namespace prototype

void expose_interpolators(py::module &m) {

  py::class_<
      ti::OneDimensionalInterpolator<double, Eigen::VectorXd>,
      std::shared_ptr<ti::OneDimensionalInterpolator<double, Eigen::VectorXd>>>
      OneDimensionalInterpolator(m, "OneDimensionalInterpolator");

  py::class_<ti::CubicSplineInterpolator<double, Eigen::VectorXd>,
             std::shared_ptr<ti::CubicSplineInterpolator<double, Eigen::VectorXd>>,
             ti::OneDimensionalInterpolator<double, Eigen::VectorXd>>
      CubicSplineInterpolator(m, "CubicSplineInterpolator");

  py::class_<ti::HermiteCubicSplineInterpolator<double, Eigen::VectorXd>,
             std::shared_ptr<ti::HermiteCubicSplineInterpolator<double, Eigen::VectorXd>>,
             ti::OneDimensionalInterpolator<double, Eigen::VectorXd>>
      HermiteCubicSplineInterpolator(m, "HermiteCubicSplineInterpolator");

  py::class_<ti::LagrangeInterpolator<double, Eigen::VectorXd>,
             std::shared_ptr<ti::LagrangeInterpolator<double, Eigen::VectorXd>>,
             ti::OneDimensionalInterpolator<double, Eigen::VectorXd>>
      LagrangeInterpolator(m, "LagrangeInterpolator");

  py::class_<ti::LinearInterpolator<double, Eigen::VectorXd>,
             std::shared_ptr<ti::LinearInterpolator<double, Eigen::VectorXd>>,
             ti::OneDimensionalInterpolator<double, Eigen::VectorXd>>
      LinearInterpolator(m, "LinearInterpolator");

//  m.def("interp1d", py::overload_cast<const Eigen::VectorXd &, const Eigen::MatrixXd &, const std::string &, const int &, const int &>(&prototype::interp1d),
//        py::arg("x"),
//        py::arg("y"),
//        py::arg("kind") = "linear",
//        py::arg("axis") = -1,
//        py::arg("order") = 8);
//
//  m.def("interp1d", py::overload_cast<std::map<double, Eigen::VectorXd> &, const std::string &, const int &>(&prototype::interp1d),
//        py::arg("result"),
//        py::arg("kind") = "linear",
//        py::arg("order") = 8);

  py::enum_<ti::BoundaryInterpolationType>(m, "BoundaryInterpolationType")
      .value("throw_exception_at_boundary", ti::BoundaryInterpolationType::throw_exception_at_boundary)
      .value("use_boundary_value", ti::BoundaryInterpolationType::use_boundary_value)
      .value("use_boundary_value_with_warning", ti::BoundaryInterpolationType::use_boundary_value_with_warning)
      .value("extrapolate_at_boundary", ti::BoundaryInterpolationType::extrapolate_at_boundary)
      .value("extrapolate_at_boundary_with_warning", ti::BoundaryInterpolationType::extrapolate_at_boundary_with_warning)
      .value("use_default_value", ti::BoundaryInterpolationType::use_default_value)
      .value("use_default_value_with_warning", ti::BoundaryInterpolationType::use_default_value_with_warning)
      .export_values();

  py::enum_<ti::AvailableLookupScheme>(m, "AvailableLookupScheme")
      .value("hunting_algorithm", ti::AvailableLookupScheme::huntingAlgorithm)
      .value("binary_search", ti::AvailableLookupScheme::binarySearch)
      .export_values();

  py::enum_<ti::LagrangeInterpolatorBoundaryHandling>(m, "LagrangeInterpolatorBoundaryHandling")
      .value("lagrange_cubic_spline_boundary_interpolation", ti::LagrangeInterpolatorBoundaryHandling::lagrange_no_boundary_interpolation)
      .value("lagrange_no_boundary_interpolation", ti::LagrangeInterpolatorBoundaryHandling::lagrange_no_boundary_interpolation)
      .export_values();

  py::class_<ti::InterpolatorSettings,
             std::shared_ptr<ti::InterpolatorSettings>>
      InterpolatorSettings(m, "InterpolatorSettings");

  py::class_<
      ti::LagrangeInterpolatorSettings,
      std::shared_ptr<ti::LagrangeInterpolatorSettings>,
      ti::InterpolatorSettings>(m,
                                "LagrangeInterpolatorSettings",
                                "Class for providing settings to creating a Lagrange interpolator.")
      .def(py::init<
               const int,
               const bool,
               const ti::AvailableLookupScheme,
               const ti::LagrangeInterpolatorBoundaryHandling,
               const ti::BoundaryInterpolationType>(),
           py::arg("interpolate_order"),
           py::arg("use_long_double_time_step") = 0,
           py::arg("selected_lookup_scheme") = ti::huntingAlgorithm,
           py::arg("lagrange_boundary_handling") = ti::lagrange_cubic_spline_boundary_interpolation,
           py::arg("boundary_handling") = ti::extrapolate_at_boundary);
};

}// namespace tudatpy
