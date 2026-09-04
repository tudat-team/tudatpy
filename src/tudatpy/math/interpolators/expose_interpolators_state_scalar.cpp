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
                  py::arg( "independent_variable_value" ) )
            .def_property_readonly( "independent_values",
                                    &ti::OneDimensionalInterpolator< double, STATE_SCALAR_TYPE >::getIndependentValues )
            .def_property_readonly( "dependent_values", &ti::OneDimensionalInterpolator< double, STATE_SCALAR_TYPE >::getDependentValues );

    py::class_< ti::OneDimensionalInterpolator< Time, STATE_SCALAR_TYPE >,
                std::shared_ptr< ti::OneDimensionalInterpolator< Time, STATE_SCALAR_TYPE > > >(
            m,
            "OneDimensionalInterpolatorScalarTimeObject",
            R"doc(
Same as OneDimensionalInterpolatorScalar, using Time as the independent variable.
      )doc" )
            .def( "interpolate",
                  py::overload_cast< const Time >( &ti::OneDimensionalInterpolator< Time, STATE_SCALAR_TYPE >::interpolate ),
                  py::arg( "independent_variable_value" ) );

    m.def( "create_one_dimensional_scalar_interpolator",
           &ti::createOneDimensionalInterpolatorBasic< double, STATE_SCALAR_TYPE >,
           py::arg( "data_to_interpolate" ),
           py::arg( "interpolator_settings" ),
           py::arg( "data_first_derivatives" ) = std::vector< double >( ),
           R"doc(
Create an interpolator for scalar dependent variables and float independent variables.
      )doc" );

    m.def( "create_one_dimensional_scalar_interpolator_time_object",
           &ti::createOneDimensionalInterpolatorBasic< Time, STATE_SCALAR_TYPE >,
           py::arg( "data_to_interpolate" ),
           py::arg( "interpolator_settings" ),
           py::arg( "data_first_derivatives" ) = std::vector< STATE_SCALAR_TYPE >( ),
           R"doc(
Create an interpolator for scalar dependent variables and Time independent variables.
      )doc" );
}

}  // namespace TUDATPY_STATE_SCALAR_BINDING_NAMESPACE
}  // namespace interpolators
}  // namespace math
}  // namespace tudatpy
