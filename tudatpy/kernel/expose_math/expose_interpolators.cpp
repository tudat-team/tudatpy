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
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

namespace py = pybind11;

namespace ti = tudat::interpolators;

namespace tudat{

namespace interpolators
{

template< typename IndependentVariableType, typename DependentVariableType >
std::shared_ptr< OneDimensionalInterpolator< IndependentVariableType, DependentVariableType > >
createOneDimensionalInterpolatorBasic(
        const std::map< IndependentVariableType, DependentVariableType > dataToInterpolate,
        const std::shared_ptr< InterpolatorSettings > interpolatorSettings )
{
    return createOneDimensionalInterpolator( dataToInterpolate, interpolatorSettings );
}

}

}
namespace tudatpy {

void expose_interpolators(py::module &m) {

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
            std::shared_ptr<ti::InterpolatorSettings>>(m, "InterpolatorSettings");

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



    m.def("linear_interpolation", &ti::linearInterpolation,
          py::arg( "lookup_scheme" ) = ti::huntingAlgorithm,
          py::arg( "boundary_interpolation" ) = ti::extrapolate_at_boundary_with_warning);

    m.def("cubic_spline_interpolation", &ti::cubicSplineInterpolation,
          py::arg( "lookup_scheme" ) = ti::huntingAlgorithm,
          py::arg( "boundary_interpolation" ) = ti::extrapolate_at_boundary_with_warning);

    m.def("piecewise_constant_interpolation", &ti::piecewiseConstantInterpolation,
          py::arg( "lookup_scheme" ) = ti::huntingAlgorithm,
          py::arg( "boundary_interpolation" ) = ti::extrapolate_at_boundary_with_warning);

    m.def("lagrange_interpolation", &ti::lagrangeInterpolation,
          py::arg( "order" ),
          py::arg( "lookup_scheme" ) = ti::huntingAlgorithm,
          py::arg( "boundary_interpolation" ) = ti::extrapolate_at_boundary_with_warning,
          py::arg( "lagrange_boundary_handling" ) = ti::lagrange_cubic_spline_boundary_interpolation );

    m.def("create_one_dimensional_interpolator",
          &ti::createOneDimensionalInterpolatorBasic< double, double >,
          py::arg("data_to_interpolate"),
          py::arg("interpolator_settings") );

    m.def("create_one_dimensional_interpolator",
          &ti::createOneDimensionalInterpolatorBasic< double, Eigen::VectorXd >,
          py::arg("data_to_interpolate"),
          py::arg("interpolator_settings") );

    m.def("create_one_dimensional_interpolator",
          &ti::createOneDimensionalInterpolatorBasic< double, Eigen::MatrixXd >,
          py::arg("data_to_interpolate"),
          py::arg("interpolator_settings") );

    py::class_<
            ti::OneDimensionalInterpolator<double, double>,
            std::shared_ptr<ti::OneDimensionalInterpolator<double, double>>>(m, "OneDimensionalInterpolatorDouble")
            .def("interpolate", py::overload_cast< const double >(
                     &ti::OneDimensionalInterpolator<double, double>::interpolate ),
                 py::arg("independent_variable_value") )
            .def("interpolate", py::overload_cast< const std::vector< double >& >(
                     &ti::OneDimensionalInterpolator<double, double>::interpolate ),
                 py::arg("independent_variable_values") );

    py::class_<
            ti::OneDimensionalInterpolator<double, Eigen::VectorXd>,
            std::shared_ptr<ti::OneDimensionalInterpolator<double, Eigen::VectorXd>>>(m, "OneDimensionalInterpolatorVector")
            .def("interpolate", py::overload_cast< const double >(
                     &ti::OneDimensionalInterpolator<double, Eigen::VectorXd>::interpolate ),
                 py::arg("independent_variable_value") )
            .def("interpolate", py::overload_cast< const std::vector< double >& >(
                     &ti::OneDimensionalInterpolator<double, Eigen::VectorXd>::interpolate ),
                 py::arg("independent_variable_values") );

    py::class_<
            ti::OneDimensionalInterpolator<double, Eigen::MatrixXd>,
            std::shared_ptr<ti::OneDimensionalInterpolator<double, Eigen::MatrixXd>>>(m, "OneDimensionalInterpolatorMatrix")
            .def("interpolate", py::overload_cast< const double >(
                     &ti::OneDimensionalInterpolator<double, Eigen::MatrixXd>::interpolate ),
                 py::arg("independent_variable_value") )
            .def("interpolate", py::overload_cast< const std::vector< double >& >(
                     &ti::OneDimensionalInterpolator<double, Eigen::MatrixXd>::interpolate ),
                 py::arg("independent_variable_values") );

};

}// namespace tudatpy
