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

    py::class_<ti::InterpolatorSettings,
            std::shared_ptr<ti::InterpolatorSettings>>(m, "InterpolatorSettings");

    py::class_<ti::LagrangeInterpolatorSettings,
            std::shared_ptr<ti::LagrangeInterpolatorSettings>,
            ti::InterpolatorSettings>(m, "LagrangeInterpolatorSettings");



    m.def("linear_interpolation", &ti::linearInterpolation );

    m.def("cubic_spline_interpolation", &ti::cubicSplineInterpolation );

    m.def("piecewise_constant_interpolation", &ti::piecewiseConstantInterpolation );

    m.def("lagrange_interpolation", &ti::lagrangeInterpolation,
          py::arg( "order" ) );


    m.def("create_one_dimensional_interpolator",
          &ti::createOneDimensionalInterpolatorBasic< double, Eigen::VectorXd >,
          py::arg("data_to_interpolate"),
          py::arg("interpolator_settings") );

    m.def("create_one_dimensional_interpolator",
          &ti::createOneDimensionalInterpolatorBasic< double, Eigen::MatrixXd >,
          py::arg("data_to_interpolate"),
          py::arg("interpolator_settings") );



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
