/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <tudat/math/interpolators.h>

#include "tudatpy/docstrings.h"
#include "tudatpy/scalarTypes.h"

namespace py = pybind11;

namespace ti = tudat::interpolators;

namespace tudat {

    namespace interpolators {

        template <typename IndependentVariableType,
                  typename DependentVariableType>
        std::shared_ptr<OneDimensionalInterpolator<IndependentVariableType,
                                                   DependentVariableType>>
        createOneDimensionalInterpolatorBasic(
            const std::map<IndependentVariableType, DependentVariableType>
                dataToInterpolate,
            const std::shared_ptr<InterpolatorSettings> interpolatorSettings,
            const std::vector<DependentVariableType>
                firstDerivativesOfDataToIntepolate =
                    std::vector<DependentVariableType>()) {
            return createOneDimensionalInterpolator<IndependentVariableType,
                                                    DependentVariableType>(
                dataToInterpolate, interpolatorSettings,
                std::make_pair(IdentityElement::getAdditionIdentity<
                                   DependentVariableType>(),
                               IdentityElement::getAdditionIdentity<
                                   DependentVariableType>()),
                firstDerivativesOfDataToIntepolate);
        }

    }  // namespace interpolators

}  // namespace tudat
namespace tudatpy {
    namespace math {
        namespace interpolators {

            PYBIND11_MODULE(expose_interpolators, m) {
                py::enum_<ti::BoundaryInterpolationType>(
                    m, "BoundaryInterpolationType",
R"doc(Enumeration of types of behaviour to be used beyond the edges of the interpolation domain.

	Enumeration of types of behaviour to be used beyond the edges of the interpolation domain. For independent variable
	data in the range :math:`[t_{0}..t_{N}]`, this enum is used to define the behaviour of the interpolator at
	:math:`t<t_{0}` and :math:`t>t_{N}
	

	:member throw_exception_at_boundary:
	:member use_boundary_value:
	:member use_boundary_value_with_warning:
	:member extrapolate_at_boundary:
	:member extrapolate_at_boundary_with_warning:
)doc")
                    .value("throw_exception_at_boundary",
                           ti::BoundaryInterpolationType::
                               throw_exception_at_boundary,
get_docstring("BoundaryInterpolationType.throw_exception_at_boundary").c_str())
                    .value("use_boundary_value",
                           ti::BoundaryInterpolationType::use_boundary_value,
get_docstring("BoundaryInterpolationType.use_boundary_value").c_str())
                    .value("use_boundary_value_with_warning",
                           ti::BoundaryInterpolationType::
                               use_boundary_value_with_warning,
get_docstring("BoundaryInterpolationType.use_boundary_value_with_warning").c_str())
                    .value(
                        "extrapolate_at_boundary",
                        ti::BoundaryInterpolationType::extrapolate_at_boundary,
get_docstring("BoundaryInterpolationType.extrapolate_at_boundary").c_str())
                    .value("extrapolate_at_boundary_with_warning",
                           ti::BoundaryInterpolationType::
                               extrapolate_at_boundary_with_warning,
get_docstring("BoundaryInterpolationType.extrapolate_at_boundary_with_warning").c_str())
                    //            .value("use_default_value",
                    //            ti::BoundaryInterpolationType::use_default_value)
                    //            .value("use_default_value_with_warning",
                    //            ti::BoundaryInterpolationType::use_default_value_with_warning)
                    .export_values();

                py::enum_<ti::AvailableLookupScheme>(
                    m, "AvailableLookupScheme",
R"doc(Enumeration of types of behaviour to be used beyond the edges of the interpolation domain.

	When the interpolation is performed, the interpolator scheme will typically start by finding the nearest neighbor of
	the requested value of the independent variable :math:`t` in the data set :math:`[t_{0}..t_{N}]`.
	The choice of lookup scheme can have a significant influence on computational efficiency for large data sets and/or simple
	interpolation algorithms
	

	:member hunting_algorithm:
	:member binary_search:
)doc")
                    .value(
                        "hunting_algorithm",
                        ti::AvailableLookupScheme::huntingAlgorithm,
get_docstring("AvailableLookupScheme.hunting_algorithm").c_str())
                    .value("binary_search",
                           ti::AvailableLookupScheme::binarySearch,
get_docstring("AvailableLookupScheme.binary_search").c_str())
                    .export_values();

                py::enum_<ti::LagrangeInterpolatorBoundaryHandling>(
                    m, "LagrangeInterpolatorBoundaryHandling",
R"doc(Enumeration of types of behaviour to be used close to the edges of the interpolation domain, for the Lagrange interpolator.

	Enumeration of types of behaviour to be used close to the edges of the interpolation domain, for the Lagrange interpolator.
	As explained for :func:`lagrange_interpolation`, the algorithm for the Lagrange interpolation breaks down at the edges of
	the interpolation domain. This enum provides the available options a user has to deal with this.
	

	:member lagrange_cubic_spline_boundary_interpolation:
	:member lagrange_no_boundary_interpolation:
)doc")
                    .value("lagrange_cubic_spline_boundary_interpolation",
                           ti::LagrangeInterpolatorBoundaryHandling::
                               lagrange_no_boundary_interpolation,
get_docstring("LagrangeInterpolatorBoundaryHandling.lagrange_cubic_spline_boundary_interpolation").c_str())
                    .value("lagrange_no_boundary_interpolation",
                           ti::LagrangeInterpolatorBoundaryHandling::
                               lagrange_no_boundary_interpolation,
get_docstring("LagrangeInterpolatorBoundaryHandling.lagrange_no_boundary_interpolation").c_str())
                    .export_values();

                py::class_<ti::InterpolatorSettings,
                           std::shared_ptr<ti::InterpolatorSettings>>(
                    m, "InterpolatorSettings",
R"doc(Base class to define settings for an interpolator.

)doc");

                py::class_<ti::InterpolatorGenerationSettings<TIME_TYPE>,
                           std::shared_ptr<
                               ti::InterpolatorGenerationSettings<TIME_TYPE>>>(
                    m, "InterpolatorGenerationSettings",
get_docstring("InterpolatorGenerationSettings").c_str());


                py::class_<ti::LagrangeInterpolatorSettings,
                           std::shared_ptr<ti::LagrangeInterpolatorSettings>,
                           ti::InterpolatorSettings>(
                    m, "LagrangeInterpolatorSettings",
R"doc(:class:`InterpolatorSettings`-derived class to define settings for a Lagrange interpolator.

)doc")
                    .def(
                        py::init<const int, const bool,
                                 const ti::AvailableLookupScheme,
                                 const ti::LagrangeInterpolatorBoundaryHandling,
                                 const ti::BoundaryInterpolationType>(),
                        py::arg("interpolate_order"),
                        py::arg("use_long_double_time_step") = 0,
                        py::arg("selected_lookup_scheme") =
                            ti::huntingAlgorithm,
                        py::arg("lagrange_boundary_handling") =
                            ti::lagrange_cubic_spline_boundary_interpolation,
                        py::arg("boundary_handling") =
                            ti::extrapolate_at_boundary);


                m.def("linear_interpolation", &ti::linearInterpolation,
                      py::arg("lookup_scheme") = ti::huntingAlgorithm,
                      py::arg("boundary_interpolation") =
                          ti::extrapolate_at_boundary_with_warning,
R"doc(Function to create settings for linear interpolation.

	Function to create settings for linear interpolation, where the interpolator
	defines a linear curve between each two subsequent intervals of the
	independent variable input data.
	

	:param lookup_scheme:
		Algorithm used to find the nearest neighbor in the independent variable input data when the interpolation scheme is called
	:param boundary_interpolation:
		Interpolator behaviour that is to be used beyond the upper and lower boundaries of the independent variable input data
	:return:
		Linear interpolation settings object
)doc");

                m.def("cubic_spline_interpolation",
                      &ti::cubicSplineInterpolation,
                      py::arg("lookup_scheme") = ti::huntingAlgorithm,
                      py::arg("boundary_interpolation") =
                          ti::extrapolate_at_boundary_with_warning,
R"doc(Function to create settings for cubic spline interpolation.

	Function to create settings for cubic spline interpolation, where the interpolator
	defines a cubic curve polynomial curve between each two subsequent intervals of the
	independent variable input data. The curve has continuous value, first derivative and
	second derivate between subsequent intervals. As boundary condition, the spline has
	a zero second derivative imposed at the upper and lower boundaries of the interpolation
	domain.
	

	:param lookup_scheme:
		Algorithm used to find the nearest neighbor in the independent variable input data when the interpolation scheme is called
	:param boundary_interpolation:
		Interpolator behaviour that is to be used beyond the upper and lower boundaries of the independent variable input data
	:return:
		Cubic spline settings object
)doc");

                m.def(
                    "piecewise_constant_interpolation",
                    &ti::piecewiseConstantInterpolation,
                    py::arg("lookup_scheme") = ti::huntingAlgorithm,
                    py::arg("boundary_interpolation") =
                        ti::extrapolate_at_boundary_with_warning,
R"doc(Function to create settings for piecewise constant interpolation.

	Function to create settings for piecewise constant interpolation. If interpolator
	is to return the value at :math:`t`, and :math:`t_{i}\le t \< t_{i+1}`, the interpolator
	returns :math:`\mathbf{x}_{i}`
	

	:param lookup_scheme:
		Algorithm used to find the nearest neighbor in the independent variable input data when the interpolation scheme is called
	:param boundary_interpolation:
		Interpolator behaviour that is to be used beyond the upper and lower boundaries of the independent variable input data
	:return:
		Piecewise constant interpolation settings object
)doc");

                m.def("lagrange_interpolation", &ti::lagrangeInterpolation,
                      py::arg("number_of_points"),
                      py::arg("lookup_scheme") = ti::huntingAlgorithm,
                      py::arg("boundary_interpolation") =
                          ti::extrapolate_at_boundary_with_warning,
                      py::arg("lagrange_boundary_handling") =
                          ti::lagrange_cubic_spline_boundary_interpolation,
R"doc(Function to create settings for cubic Lagrange interpolation.

	Function to create settings for piecewise cubic Lagrange interpolation.  This is typically the interpolator of highest accuracy that is available.
	The Lagrange interpolator uses :math:`m` consecutive points (input to this function) from the input independent variables :math:`[t_{0}...t_{N}]`
	to create the polynomial of order :math:`m-1` that interpolates these points. From here on, we assume :math:`m` is even.
	The algorithm that is used (see `here <https://mathworld.wolfram.com/LagrangeInterpolatingPolynomial.html>`_ for mathematical details
	on interpolating Lagrange polynomials) works as follows:
	
	* The nearest lower neighbor of the data point :math:`t` at which the state :math:`\mathbf{x}` is to be interpolated is determined, and denoted :math:`t_{i}`.
	* An interpolating Lagrange polynomial is constructed from the consecutive data points :math:`[t_{i-(m/2-1)}...t_{i+m}]`
	* This resulting interpolating polynomial is *only* used in the interval :math:`[t_{i}..t_{i+1}]`, to prevent `Runge's phenomenon <https://en.wikipedia.org/wiki/Runge%27s_phenomenon>`_.
	
	For instance, if :math:`m=8` we use a :math:`7^{th}` order polynomial that interpolates a contiguous set of
	8 data points out of the full data set. Normally, the interpolating polynomial is only used between the
	:math:`4^{th}` and :math:`5^{th}` data point, where it will typically be of good accuracy. Consequently,
	a separate interpolating polynomial (using data over a span of :math:`m` consecutive points) is used for
	each single interval :math:`[t_{i}..t_{i+1}]` (with the exception of the boundaries, see below).
	
	.. warning:: Issues can occur if the data point :math:`t` at which the interpolation is to be
	             performed is close to :math:`t_{0}` or :math:`t_{N}`. In those case, there is not sufficient
	             data to construct the interpolating polynomial *and* to only use this interpolating polynomial
	             between the middle two data points that were used to it. In these cases, the user has a number of
	             options (all defined by an entry of the :class:`LagrangeInterpolatorBoundaryHandling` variable,
	             used as input to this function). In short, interpolation between the first and last :math:`m/2`
	             data points will lead to degraded results, warnings, or termination.
	

	:param number_of_points:
		Number of consecutive data points that are used to construct a single interpolating polynomial.
	:param lookup_scheme:
		Algorithm used to find the nearest neighbor in the independent variable input data when the interpolation scheme is called
	:param boundary_interpolation:
		Interpolator behaviour that is to be used beyond the upper and lower boundaries of the independent variable input data
	:param lagrange_boundary_handling:
		Interpolator behaviour that is to be used at the boundaries of the domain, where the regular algorithm cannot be executed.
	:return:
		Lagrange interpolation settings object
)doc");

                m.def("hermite_spline_interpolation", &ti::hermiteInterpolation,
                      py::arg("lookup_scheme") = ti::huntingAlgorithm,
                      py::arg("boundary_interpolation") =
                          ti::extrapolate_at_boundary_with_warning,
R"doc(Function to create settings for cubic Hermite spline interpolation.

	Function to create settings for piecewise cubic Hermite spline interpolation. To use this
	interpolator, a key-value container of values, and a key-value container of first derivatives,
	must be provided to the function creating an interpolator (:func:`create_one_dimensional_scalar_interpolator`,
	:func:`create_one_dimensional_vector_interpolator`,  :func:`create_one_dimensional_matrix_interpolator`). The resulting
	spline uses the value and first derivatives (four piece of information for each interval) at two subsequent nodes to construct
	a cubic polynomial between each two subsequent nodes. The resulting spline has constant values and first derivatives
	

	:param lookup_scheme:
		Algorithm used to find the nearest neighbor in the independent variable input data when the interpolation scheme is called
	:param boundary_interpolation:
		Interpolator behaviour that is to be used beyond the upper and lower boundaries of the independent variable input data
	:return:
		Hermite spline interpolation settings object
)doc");

                m.def("hermite_interpolation", &ti::hermiteInterpolation,
                      py::arg("lookup_scheme") = ti::huntingAlgorithm,
                      py::arg("boundary_interpolation") =
                          ti::extrapolate_at_boundary_with_warning);

                m.def(
                    "create_one_dimensional_scalar_interpolator",
                    &ti::createOneDimensionalInterpolatorBasic<TIME_TYPE,
                                                               double>,
                    py::arg("data_to_interpolate"),
                    py::arg("interpolator_settings"),
                    py::arg("data_first_derivatives") = std::vector<double>(),
R"doc(Function to create an interpolator for scalar dependent variables.

	Function to create an interpolator for scalar dependent variables, with a single independent
	variable. This function takes the interpolator settings, and the data that is to be interpolated,
	as input to create the object that can perform the actual interpolation
	

	:param data_to_interpolate:
		Key-value container with pairs of independent variables (key) and dependent variables (value) from which the interpolation is to be performed
	:param interpolator_settings:
		Settings that define the type of interpolator that is to be used
	:param data_first_derivatives:
		Key-value container with pairs of independent variables (key) and first derivative dependent variables w.r.t. independent variable (value) from which the interpolation is to be performed. This input is *only* required if the requested interpolation algorithm requires first derivatives as input (such as the Hermite spline interpolator
	:return:
		Interpolator object
)doc");

                m.def(
                    "create_one_dimensional_vector_interpolator",
                    &ti::createOneDimensionalInterpolatorBasic<TIME_TYPE,
                                                               Eigen::VectorXd>,
                    py::arg("data_to_interpolate"),
                    py::arg("interpolator_settings"),
                    py::arg("data_first_derivatives") =
                        std::vector<Eigen::VectorXd>(),
R"doc(Function to create an interpolator for vector dependent variables.

	As :func:`create_one_dimensional_scalar_interpolator`, but with vectors as dependent variables
	

	:param data_to_interpolate:
		Key-value container with pairs of independent variables (key) and dependent variables (value) from which the interpolation is to be performed
	:param interpolator_settings:
		Settings that define the type of interpolator that is to be used
	:param data_first_derivatives:
		Key-value container with pairs of independent variables (key) and first derivative dependent variables w.r.t. independent variable (value) from which the interpolation is to be performed. This input is *only* required if the requested interpolation algorithm requires first derivatives as input (such as the Hermite spline interpolator).
	:return:
		Interpolator object
)doc");

                m.def(
                    "create_one_dimensional_matrix_interpolator",
                    &ti::createOneDimensionalInterpolatorBasic<TIME_TYPE,
                                                               Eigen::MatrixXd>,
                    py::arg("data_to_interpolate"),
                    py::arg("interpolator_settings"),
                    py::arg("data_first_derivatives") =
                        std::vector<Eigen::MatrixXd>(),
R"doc(Function to create an interpolator for matrix dependent variables.

	As :func:`create_one_dimensional_scalar_interpolator`, but with matrices (2-dimensional arrays) as dependent variables
	

	:param data_to_interpolate:
		Key-value container with pairs of independent variables (key) and dependent variables (value) from which the interpolation is to be performed
	:param interpolator_settings:
		Settings that define the type of interpolator that is to be used
	:param data_first_derivatives:
		Key-value container with pairs of independent variables (key) and first derivative dependent variables w.r.t. independent variable (value) from which the interpolation is to be performed. This input is *only* required if the requested interpolation algorithm requires first derivatives as input (such as the Hermite spline interpolator
	:return:
		Interpolator object
)doc");

                py::class_<ti::OneDimensionalInterpolator<TIME_TYPE, double>,
                           std::shared_ptr<ti::OneDimensionalInterpolator<
                               TIME_TYPE, double>>>(
                    m, "OneDimensionalInterpolatorScalar",
R"doc(Object that performs interpolation for scalar independent, and scalar dependent variables.

	Object that performs interpolation for scalar independent, and scalar dependent variables. This object is
	not created manually, but is set up using the :func:`create_one_dimensional_scalar_interpolator` function.
	
)doc")
                    .def("interpolate",
                         py::overload_cast<const TIME_TYPE>(
                             &ti::OneDimensionalInterpolator<
                                 TIME_TYPE, double>::interpolate),
                         py::arg("independent_variable_value"),
R"doc(This function performs the interpolation at the requested independent variable value.

	:param independent_variable_value:
		Value of independent variable at which the interpolation is to bse performed.
		
	:return:
		Interpolated dependent variable value, using implemented algorithm at requested independent variable value
)doc");

                py::class_<
                    ti::OneDimensionalInterpolator<TIME_TYPE, Eigen::VectorXd>,
                    std::shared_ptr<ti::OneDimensionalInterpolator<
                        TIME_TYPE, Eigen::VectorXd>>>(
                    m, "OneDimensionalInterpolatorVector",
R"doc(Object that performs interpolation for vector independent, and vector dependent variables.

	Object that performs interpolation for vector independent, and vector dependent variables. This object is
	not created manually, but is set up using the :func:`create_one_dimensional_vector_interpolator` function.
	
)doc")
                    .def("interpolate",
                         py::overload_cast<const TIME_TYPE>(
                             &ti::OneDimensionalInterpolator<
                                 TIME_TYPE, Eigen::VectorXd>::interpolate),
                         py::arg("independent_variable_value"),
R"doc(This function performs the interpolation at the requested independent variable value.

	:param independent_variable_value:
		Value of independent variable at which the interpolation is to be performed.
		
	:return:
		Interpolated dependent variable value, using implemented algorithm at requested independent variable value
)doc");

                py::class_<
                    ti::OneDimensionalInterpolator<TIME_TYPE, Eigen::MatrixXd>,
                    std::shared_ptr<ti::OneDimensionalInterpolator<
                        TIME_TYPE, Eigen::MatrixXd>>>(
                    m, "OneDimensionalInterpolatorMatrix",
R"doc(Object that performs interpolation for matrix independent, and matrix dependent variables.

	Object that performs interpolation for matrix independent, and matrix dependent variables. This object is
	not created manually, but is set up using the :func:`create_one_dimensional_matrix_interpolator` function.
	
)doc")
                    .def("interpolate",
                         py::overload_cast<const TIME_TYPE>(
                             &ti::OneDimensionalInterpolator<
                                 TIME_TYPE, Eigen::MatrixXd>::interpolate),
                         py::arg("independent_variable_value"),
R"doc(This function performs the interpolation at the requested independent variable value.

	:param independent_variable_value:
		Value of independent variable at which the interpolation is to be performed.
		
	:return:
		Interpolated dependent variable value, using implemented algorithm at requested independent variable value
)doc");
            }

        }  // namespace interpolators
    }  // namespace math
}  // namespace tudatpy
