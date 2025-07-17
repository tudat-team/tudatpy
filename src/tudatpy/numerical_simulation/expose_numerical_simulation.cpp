/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#define PYBIND11_DETAILED_ERROR_MESSAGES
#include "expose_numerical_simulation.h"

#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "scalarTypes.h"
#include "tudat/astro/basic_astro/dateTime.h"
#include "tudat/basics/timeType.h"

namespace py = pybind11;
namespace tp = tudat::propagators;
namespace tss = tudat::simulation_setup;
namespace tni = tudat::numerical_integrators;
namespace tep = tudat::estimatable_parameters;
namespace tom = tudat::observation_models;
namespace tba = tudat::basic_astrodynamics;
namespace pybind11
{
namespace detail
{
template<>
struct type_caster< double > {
public:
    // Register implicit conversion from Time -> float/double
    PYBIND11_TYPE_CASTER( double, _( "float" ) );

    bool load( handle src, bool convert )
    {
        // Handle existing conversions first
        if( PyFloat_Check( src.ptr( ) ) )
        {
            value = PyFloat_AsDouble( src.ptr( ) );
            return !PyErr_Occurred( );
        }

        // Try to convert from Time
        if( convert )
        {
            try
            {
                auto time_type = module_::import( "tudatpy.kernel" ).attr( "Time" );
                if( isinstance( src, time_type ) )
                {
                    value = src.attr( "to_float" )( ).cast< double >( );
                    return true;
                }
            }
            catch( const error_already_set& )
            {
                PyErr_Clear( );
            }
        }
        return false;
    }

    static handle cast( double src, return_value_policy policy, handle parent )
    {
        return PyFloat_FromDouble( src );
    }
};
}  // namespace detail
}  // namespace pybind11
namespace tudatpy
{

namespace numerical_simulation
{

void expose_numerical_simulation( py::module& m )
{
    py::class_< tudat::Time >( m, "Time", R"doc(
        
    Class for defining time with a resolution that is sub-femtosecond for very long periods of time.
    
    Using double or long double precision as a representation of time, the issue of reduced precision will 
    occur over long time periods. For instance, over a period of 10^8 seconds (about 3 years), double and 
    long double representations have resolution of about 10^-8 and 10^-11 s respectively, which is 
    insufficient for various applications. 
    
    This class uses an integer to represent the number of hours since an epoch, and a long double to 
    represent the number of seconds into the present hour. This provides a resolution of < 1 femtosecond, 
    over a range of 2147483647 hours (about 300,000 years), which is more than sufficient for practical 
    applications.
    
    The Time class supports standard arithmetic operations (+, -, *, /) with Time objects and floats, comparison operations, and 
    automatic conversion to/from floating-point types.
        )doc" )
            .def( py::init< const double >( ),
                  py::arg( "seconds_since_j2000" ),
                  R"doc(

     Create a Time object from seconds since J2000.
     
     Parameters
     ----------
     seconds_since_j2000 : float
         Number of seconds since J2000 epoch
     
     Returns
     -------
     Time
         Time object initialized to specified seconds since J2000
     
     Examples
     --------
     >>> from tudatpy.kernel import Time
     >>> t = Time(3600.0)  # 1 hour after J2000
     )doc" )

            // Add docstring for the second constructor
            .def( py::init< const int, const long double >( ),
                  py::arg( "full_periods" ),
                  py::arg( "seconds_into_full_period" ),
                  R"doc(

     Create a Time object from full periods (hours) and seconds into the current period.
     
     Parameters
     ----------
     full_periods : int
         Number of full hours since epoch
     seconds_into_full_period : float
         Number of seconds into current hour. Need not be in range [0, 3600];
         the time representation is normalized automatically.
     
     Returns
     -------
     Time
         Time object initialized to specified time
     
     Examples
     --------
     >>> from tudatpy.kernel import Time
     >>> t = Time(2, 1800.0)  # 2.5 hours after epoch
     )doc" )
            .def( "to_float",
                  &tudat::Time::getSeconds< double >,
                  R"doc(
    Converts the time to a float (double) representing seconds since J2000.
        
    Returns
    -------
    float
        Number of seconds since J2000
    
    Examples
    --------
    In this example, a Time object is converted back to seconds since J2000.
    
    .. code-block:: python
    
        from tudatpy.kernel import Time
        
        # Create a Time object from seconds since J2000
        t = Time(3600.0)  # 1 hour after J2000
        
        # Convert back to seconds
        seconds = t.to_float()
        print(seconds)  # prints 3600.0
    )doc" )
            .def( "__float__", &tudat::Time::getSeconds< double > )
            .def( py::self + py::self )
            .def( py::self + double( ) )
            .def( double( ) + py::self )
            .def( py::self += py::self )
            .def( py::self += double( ) )
            .def( py::self - py::self )
            .def( py::self - double( ) )
            .def( py::self -= py::self )
            .def( py::self -= double( ) )
            .def( double( ) - py::self )
            .def( py::self * double( ) )
            .def( double( ) * py::self )
            .def( py::self *= double( ) )
            .def( py::self / double( ) )
            .def( py::self /= double( ) )
            .def( py::self == py::self )
            .def( double( ) == py::self )
            .def( py::self == double( ) )
            .def( py::self != py::self )
            .def( py::self != double( ) )
            .def( double( ) != py::self )
            .def( py::self < py::self )
            .def( py::self < double( ) )
            .def( double( ) < py::self )
            .def( py::self > py::self )
            .def( py::self > double( ) )
            .def( double( ) > py::self )
            .def( py::self <= py::self )
            .def( py::self <= double( ) )
            .def( double( ) <= py::self )
            .def( py::self >= py::self )
            .def( double( ) >= py::self )
            .def( py::self >= double( ) );

    // Register implicit conversion from float/double -> Time
    py::implicitly_convertible< double, tudat::Time >( );

    auto environment_setup_submodule = m.def_submodule( "environment_setup" );
    environment_setup::expose_environment_setup( environment_setup_submodule );

    auto environment_submodule = m.def_submodule( "environment" );
    environment::expose_environment( environment_submodule );

    auto propagation_setup_submodule = m.def_submodule( "propagation_setup" );
    propagation_setup::expose_propagation_setup( propagation_setup_submodule );

    auto propagation_submodule = m.def_submodule( "propagation" );
    propagation::expose_propagation( propagation_submodule );

    auto estimation_setup_submodule = m.def_submodule( "estimation_setup" );
    estimation_setup::expose_estimation_setup( estimation_setup_submodule );

    auto estimation_submodule = m.def_submodule( "estimation" );

    estimation::expose_estimation_filter_parser( estimation_submodule );
    estimation::expose_estimation( estimation_submodule );
    estimation::expose_estimation_observation_collection( estimation_submodule );
    estimation::expose_estimation_propagated_covariance( estimation_submodule );
    estimation::expose_estimation_single_observation_set( estimation_submodule );

    m.def( "get_integrated_type_and_body_list",
           &tp::getIntegratedTypeAndBodyList< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "propagator_settings" ) );

    m.def( "get_single_integration_size", &tp::getSingleIntegrationSize, py::arg( "state_type" ) );
}

}  // namespace numerical_simulation
}  // namespace tudatpy
