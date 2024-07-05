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

#include "tudatpy/docstrings.h"
#include "tudatpy/scalarTypes.h"

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
        const std::shared_ptr< InterpolatorSettings > interpolatorSettings,
        const std::vector< DependentVariableType > firstDerivativesOfDataToIntepolate =
        std::vector< DependentVariableType>( ) )
{
    return createOneDimensionalInterpolator< IndependentVariableType, DependentVariableType >(
                dataToInterpolate, interpolatorSettings,
                std::make_pair( IdentityElement::getAdditionIdentity< DependentVariableType >( ),
                                IdentityElement::getAdditionIdentity< DependentVariableType >( ) ),
                firstDerivativesOfDataToIntepolate );
}

}

}
namespace tudatpy {
namespace math {
namespace interpolators {

void expose_interpolators(py::module &m) {

    py::enum_<ti::BoundaryInterpolationType>(m, "BoundaryInterpolationType", get_docstring("BoundaryInterpolationType").c_str())
            .value("throw_exception_at_boundary", ti::BoundaryInterpolationType::throw_exception_at_boundary,
                   get_docstring("BoundaryInterpolationType.throw_exception_at_boundary").c_str())
            .value("use_boundary_value", ti::BoundaryInterpolationType::use_boundary_value,
                   get_docstring("BoundaryInterpolationType.use_boundary_value").c_str())
            .value("use_boundary_value_with_warning", ti::BoundaryInterpolationType::use_boundary_value_with_warning,
                   get_docstring("BoundaryInterpolationType.use_boundary_value_with_warning").c_str())
            .value("extrapolate_at_boundary", ti::BoundaryInterpolationType::extrapolate_at_boundary,
                   get_docstring("BoundaryInterpolationType.extrapolate_at_boundary").c_str())
            .value("extrapolate_at_boundary_with_warning", ti::BoundaryInterpolationType::extrapolate_at_boundary_with_warning,
                   get_docstring("BoundaryInterpolationType.extrapolate_at_boundary_with_warning").c_str())
//            .value("use_default_value", ti::BoundaryInterpolationType::use_default_value)
//            .value("use_default_value_with_warning", ti::BoundaryInterpolationType::use_default_value_with_warning)
            .export_values();

    py::enum_<ti::AvailableLookupScheme>(m, "AvailableLookupScheme", get_docstring("AvailableLookupScheme").c_str())
            .value("hunting_algorithm", ti::AvailableLookupScheme::huntingAlgorithm,
                   get_docstring("AvailableLookupScheme.hunting_algorithm").c_str() )
            .value("binary_search", ti::AvailableLookupScheme::binarySearch,
                   get_docstring("AvailableLookupScheme.binary_search").c_str() )
            .export_values();

    py::enum_<ti::LagrangeInterpolatorBoundaryHandling>(m, "LagrangeInterpolatorBoundaryHandling", get_docstring("LagrangeInterpolatorBoundaryHandling").c_str())
            .value("lagrange_cubic_spline_boundary_interpolation", ti::LagrangeInterpolatorBoundaryHandling::lagrange_no_boundary_interpolation,
                   get_docstring("LagrangeInterpolatorBoundaryHandling.lagrange_cubic_spline_boundary_interpolation").c_str() )
            .value("lagrange_no_boundary_interpolation", ti::LagrangeInterpolatorBoundaryHandling::lagrange_no_boundary_interpolation,
                   get_docstring("LagrangeInterpolatorBoundaryHandling.lagrange_no_boundary_interpolation").c_str() )
            .export_values();

    py::class_<ti::InterpolatorSettings,
            std::shared_ptr<ti::InterpolatorSettings>>(m, "InterpolatorSettings",
                                                       get_docstring("InterpolatorSettings").c_str());

    py::class_<ti::InterpolatorGenerationSettings<TIME_TYPE>,
        std::shared_ptr<ti::InterpolatorGenerationSettings<TIME_TYPE>>>(m, "InterpolatorGenerationSettings",
                                                   get_docstring("InterpolatorGenerationSettings").c_str());

    py::class_<ti::InterpolatorGenerationSettings<double>,
        std::shared_ptr<ti::InterpolatorGenerationSettings<double>>>(m, "InterpolatorGenerationSettingsFloat",
                                                                        get_docstring("InterpolatorGenerationSettingsFloat").c_str());


    py::class_<
            ti::LagrangeInterpolatorSettings,
            std::shared_ptr<ti::LagrangeInterpolatorSettings>,
            ti::InterpolatorSettings>(m,
                                      "LagrangeInterpolatorSettings",
                                      get_docstring("LagrangeInterpolatorSettings").c_str())
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

      m.def("interpolator_generation_settings", &ti::interpolatorGenerationSettings< TIME_TYPE >,
          py::arg( "interpolator_settings" ),
          py::arg( "initial_time" ),
          py::arg( "final_time" ),
          py::arg( "time_step" ),
          get_docstring("interpolator_generation_settings").c_str());

    m.def("interpolator_generation_settings_float", &ti::interpolatorGenerationSettings< double >,
          py::arg( "interpolator_settings" ),
          py::arg( "initial_time" ),
          py::arg( "final_time" ),
          py::arg( "time_step" ),
          get_docstring("interpolator_generation_settings_float").c_str());

    m.def("linear_interpolation", &ti::linearInterpolation,
          py::arg( "lookup_scheme" ) = ti::huntingAlgorithm,
          py::arg( "boundary_interpolation" ) = ti::extrapolate_at_boundary_with_warning,
          get_docstring("linear_interpolation").c_str());

    m.def("cubic_spline_interpolation", &ti::cubicSplineInterpolation,
          py::arg( "lookup_scheme" ) = ti::huntingAlgorithm,
          py::arg( "boundary_interpolation" ) = ti::extrapolate_at_boundary_with_warning,
          get_docstring("cubic_spline_interpolation").c_str());

    m.def("piecewise_constant_interpolation", &ti::piecewiseConstantInterpolation,
          py::arg( "lookup_scheme" ) = ti::huntingAlgorithm,
          py::arg( "boundary_interpolation" ) = ti::extrapolate_at_boundary_with_warning,
          get_docstring("piecewise_constant_interpolation").c_str());

    m.def("lagrange_interpolation", &ti::lagrangeInterpolation,
          py::arg( "number_of_points" ),
          py::arg( "lookup_scheme" ) = ti::huntingAlgorithm,
          py::arg( "boundary_interpolation" ) = ti::extrapolate_at_boundary_with_warning,
          py::arg( "lagrange_boundary_handling" ) = ti::lagrange_cubic_spline_boundary_interpolation,
          get_docstring("lagrange_interpolation").c_str() );

    m.def("hermite_spline_interpolation", &ti::hermiteInterpolation,
          py::arg( "lookup_scheme" ) = ti::huntingAlgorithm,
          py::arg( "boundary_interpolation" ) = ti::extrapolate_at_boundary_with_warning,
          get_docstring("hermite_spline_interpolation").c_str( ) );

    m.def("hermite_interpolation", &ti::hermiteInterpolation,
          py::arg( "lookup_scheme" ) = ti::huntingAlgorithm,
          py::arg( "boundary_interpolation" ) = ti::extrapolate_at_boundary_with_warning );

    m.def("create_one_dimensional_scalar_interpolator",
          &ti::createOneDimensionalInterpolatorBasic<TIME_TYPE, double >,
          py::arg("data_to_interpolate"),
          py::arg("interpolator_settings"),
          py::arg("data_first_derivatives") = std::vector< double  >( ),
          get_docstring("create_one_dimensional_scalar_interpolator").c_str( ) );

    m.def("create_one_dimensional_vector_interpolator",
          &ti::createOneDimensionalInterpolatorBasic<TIME_TYPE, Eigen::VectorXd >,
          py::arg("data_to_interpolate"),
          py::arg("interpolator_settings"),
          py::arg("data_first_derivatives") = std::vector< Eigen::VectorXd >( ),
          get_docstring("create_one_dimensional_vector_interpolator").c_str( ) );

    m.def("create_one_dimensional_matrix_interpolator",
          &ti::createOneDimensionalInterpolatorBasic<TIME_TYPE, Eigen::MatrixXd >,
          py::arg("data_to_interpolate"),
          py::arg("interpolator_settings"),
          py::arg("data_first_derivatives") = std::vector< Eigen::MatrixXd >( ),
          get_docstring("create_one_dimensional_matrix_interpolator").c_str( ) );

    py::class_<
            ti::OneDimensionalInterpolator<TIME_TYPE, double>,
            std::shared_ptr<ti::OneDimensionalInterpolator<TIME_TYPE, double>>>(m, "OneDimensionalInterpolatorScalar",
                                                                             get_docstring("OneDimensionalInterpolatorScalar").c_str())
            .def("interpolate", py::overload_cast< const TIME_TYPE >(
                     &ti::OneDimensionalInterpolator<TIME_TYPE, double>::interpolate ),
                 py::arg("independent_variable_value"),
                 get_docstring("OneDimensionalInterpolatorScalar.interpolate").c_str()  );

    py::class_<
            ti::OneDimensionalInterpolator<TIME_TYPE, Eigen::VectorXd>,
            std::shared_ptr<ti::OneDimensionalInterpolator<TIME_TYPE, Eigen::VectorXd>>>(m, "OneDimensionalInterpolatorVector",
                                                                                      get_docstring("OneDimensionalInterpolatorVector").c_str())
            .def("interpolate", py::overload_cast< const TIME_TYPE >(
                     &ti::OneDimensionalInterpolator<TIME_TYPE, Eigen::VectorXd>::interpolate ),
                 py::arg("independent_variable_value"),
                 get_docstring("OneDimensionalInterpolatorVector.interpolate").c_str() );

    py::class_<
            ti::OneDimensionalInterpolator<TIME_TYPE, Eigen::MatrixXd>,
            std::shared_ptr<ti::OneDimensionalInterpolator<TIME_TYPE, Eigen::MatrixXd>>>(m, "OneDimensionalInterpolatorMatrix",
                                                                                      get_docstring("OneDimensionalInterpolatorMatrix").c_str())
            .def("interpolate", py::overload_cast< const TIME_TYPE >(
                     &ti::OneDimensionalInterpolator<TIME_TYPE, Eigen::MatrixXd>::interpolate ),
                 py::arg("independent_variable_value"),
                 get_docstring("OneDimensionalInterpolatorMatrix.interpolate").c_str() );


}

}
}
}// namespace tudatpy
