/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of Tudat.
 */
#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "create_one_dimensional_interpolator_basic.h"
#include "scalarTypes.h"

namespace py = pybind11;
namespace ti = tudat::interpolators;

namespace tudatpy
{
namespace math
{
namespace interpolators
{
namespace TUDATPY_STATE_SCALAR_BINDING_NAMESPACE
{

void expose_interpolators_state_scalar( py::module& m )
{
    py::class_< ti::OneDimensionalInterpolator< double, STATE_SCALAR_TYPE >,
                std::shared_ptr< ti::OneDimensionalInterpolator< double, STATE_SCALAR_TYPE > > >( m,
                                                                                                  "OneDimensionalInterpolatorScalar",
                                                                                                  R"doc(
Object that performs interpolation for scalar dependent variables and float
independent variables.
      )doc" )
            .def( "interpolate",
                  py::overload_cast< const double >( &ti::OneDimensionalInterpolator< double, STATE_SCALAR_TYPE >::interpolate ),
                  py::arg( "independent_variable_value" ),
                  R"doc(
        This function performs the interpolation at the requested independent variable value.

        Parameters
        ----------
        independent_variable_value : astro.time_representation.Time
            Value of independent variable at which the interpolation is to be performed.
        Returns
        -------
        np.array
            Interpolated dependent variable value, using implemented algorithm at requested independent variable value
    )doc" )
            .def_property_readonly( "independent_values",
                                    &ti::OneDimensionalInterpolator< double, STATE_SCALAR_TYPE >::getIndependentValues,
                                    R"doc(

         Returns the independent variable values used by the interpolator.

         Returns the independent variable values used by the interpolator. This is a read-only property.

         Returns
         -------
         list[float]
             Independent variable values used by the interpolator
         )doc" )
            .def_property_readonly( "dependent_values",
                                    &ti::OneDimensionalInterpolator< double, STATE_SCALAR_TYPE >::getDependentValues,
                                    R"doc(

         Returns the dependent variable values used by the interpolator.

         Returns the dependent variable values used by the interpolator. This is a read-only property.

         Returns
         -------
         list[np.ndarray]
             Dependent variable values used by the interpolator
         )doc" );

    py::class_< ti::OneDimensionalInterpolator< Time, STATE_SCALAR_TYPE >,
                std::shared_ptr< ti::OneDimensionalInterpolator< Time, STATE_SCALAR_TYPE > > >(
            m,
            "OneDimensionalInterpolatorScalarTimeObject",
            R"doc(
Same as OneDimensionalInterpolatorScalar, using Time as the independent variable.
      )doc" )
            .def( "interpolate",
                  py::overload_cast< const Time >( &ti::OneDimensionalInterpolator< Time, STATE_SCALAR_TYPE >::interpolate ),
                  py::arg( "independent_variable_value" ),
                  R"doc(
        This function performs the interpolation at the requested independent variable value.

        Parameters
        ----------
        independent_variable_value : astro.time_representation.Time
            Value of independent variable at which the interpolation is to be performed.
        Returns
        -------
        np.array
            Interpolated dependent variable value, using implemented algorithm at requested independent variable value
    )doc" );

    m.def( "create_one_dimensional_scalar_interpolator",
           &ti::createOneDimensionalInterpolatorBasic< double, STATE_SCALAR_TYPE >,
           py::arg( "data_to_interpolate" ),
           py::arg( "interpolator_settings" ),
           py::arg( "data_first_derivatives" ) = std::vector< double >( ),
           R"doc(

 Function to create an interpolator for scalar dependent variables.

 Function to create an interpolator for scalar dependent variables, with a single float independent
 variable. This function takes the interpolator settings, and the data that is to be interpolated,
 as input to create the object that can perform the actual interpolation


 Parameters
 ----------
 data_to_interpolate : dict[float, float]
     Key-value container with pairs of independent variables (key) and dependent variables (value) from which the interpolation is to be performed
 interpolator_settings : InterpolatorSettings
     Settings that define the type of interpolator that is to be used
 data_first_derivatives : list[float] = []
     List of first derivative dependent variables w.r.t. independent variable from which the interpolation is to be performed. Must be of the same size as the number of data points in ``data_to_interpolate``. This input is *only* required if the requested interpolation algorithm requires first derivatives as input (such as the Hermite spline interpolator).
 Returns
 -------
 OneDimensionalInterpolatorScalar
     Interpolator object



     )doc" );

    m.def( "create_one_dimensional_scalar_interpolator_time_object",
           &ti::createOneDimensionalInterpolatorBasic< Time, STATE_SCALAR_TYPE >,
           py::arg( "data_to_interpolate" ),
           py::arg( "interpolator_settings" ),
           py::arg( "data_first_derivatives" ) = std::vector< STATE_SCALAR_TYPE >( ),
           R"doc(

     Creates a one-dimensional scalar interpolator with floating-point independent variables.

     Parameters
     ----------
     data_to_interpolate : dict[float, float]
         Key-value container with pairs of independent variables (key) and dependent variables (value) from which the interpolation is to be performed
     interpolator_settings : InterpolatorSettings
         Settings that define the type of interpolator that is to be used
     data_first_derivatives : list[float] = []
         List of first derivative dependent variables w.r.t. independent variable from which the interpolation is to be performed. Must be of the same size as the number of data points in ``data_to_interpolate``. This input is *only* required if the requested interpolation algorithm requires first derivatives as input (such as the Hermite spline interpolator).
     Returns
     -------
     OneDimensionalInterpolatorScalar
         Interpolator object

    )doc" );
}

}  // namespace TUDATPY_STATE_SCALAR_BINDING_NAMESPACE
}  // namespace interpolators
}  // namespace math
}  // namespace tudatpy
