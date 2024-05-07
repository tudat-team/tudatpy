/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_root_finders.h"

#include "docstrings.h"
#include "tudat/math/root_finders.h"

#include <pybind11/pybind11.h>

namespace py = pybind11;
namespace trf = tudat::root_finders;

namespace tudatpy
{
namespace math
{

namespace root_finders
{

void expose_root_finders( py::module &m )
{

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

    py::enum_<trf::MaximumIterationHandling>( m, "MaximumIterationHandling",
                                              get_docstring( "MaximumIterationHandling" ).c_str( ))
        .value( "accept_result", trf::MaximumIterationHandling::accept_result,
                get_docstring( "MaximumIterationHandling.accept_result" ).c_str( ))
        .value( "accept_result_with_warning", trf::MaximumIterationHandling::accept_result_with_warning,
                get_docstring( "MaximumIterationHandling.accept_result_with_warning" ).c_str( ))
        .value( "throw_exception", trf::MaximumIterationHandling::throw_exception,
                get_docstring( "MaximumIterationHandling.throw_exception" ).c_str( )).export_values( );


    py::class_<trf::RootFinder<double>,
        std::shared_ptr<trf::RootFinder<double>>>( m, "RootFinderCore" );

    py::class_<trf::NewtonRaphson<double>,
        std::shared_ptr<trf::NewtonRaphson<double>>,
        trf::RootFinder<double>>( m, "NewtonRaphsonCore" )
        .def(
            py::init<const double, const unsigned int>( ),
            py::arg( "x_tol" ),
            py::arg( "max_iter" ));

    py::class_<trf::RootFinderSettings,
        std::shared_ptr<trf::RootFinderSettings>>( m, "RootFinderSettings",
                                                   get_docstring( "RootFinderSettings" ).c_str( ));

    m.def( "bisection",
           &trf::bisectionRootFinderSettings,
           py::arg( "relative_variable_tolerance" ) = TUDAT_NAN,
           py::arg( "absolute_variable_tolerance" ) = TUDAT_NAN,
           py::arg( "root_function_tolerance" ) = TUDAT_NAN,
           py::arg( "maximum_iteration" ) = 1000,
           py::arg( "maximum_iteration_handling" ) = trf::throw_exception,
           get_docstring( "bisection" ).c_str( ));


    m.def( "newton_raphson",
           &trf::newtonRaphsonRootFinderSettings,
           py::arg( "relative_variable_tolerance" ) = TUDAT_NAN,
           py::arg( "absolute_variable_tolerance" ) = TUDAT_NAN,
           py::arg( "root_function_tolerance" ) = TUDAT_NAN,
           py::arg( "maximum_iteration" ) = 1000,
           py::arg( "maximum_iteration_handling" ) = trf::throw_exception,
           get_docstring( "newton_raphson" ).c_str( ));

    m.def( "halley",
           &trf::halleyRootFinderSettings,
           py::arg( "relative_variable_tolerance" ) = TUDAT_NAN,
           py::arg( "absolute_variable_tolerance" ) = TUDAT_NAN,
           py::arg( "root_function_tolerance" ) = TUDAT_NAN,
           py::arg( "maximum_iteration" ) = 1000,
           py::arg( "maximum_iteration_handling" ) = trf::throw_exception,
           get_docstring( "halley" ).c_str( ));

    m.def( "secant",
           &trf::secantRootFinderSettings,
           py::arg( "relative_variable_tolerance" ) = TUDAT_NAN,
           py::arg( "absolute_variable_tolerance" ) = TUDAT_NAN,
           py::arg( "root_function_tolerance" ) = TUDAT_NAN,
           py::arg( "maximum_iteration" ) = 1000,
           py::arg( "maximum_iteration_handling" ) = trf::throw_exception,
           get_docstring( "secant" ).c_str( ));
}

}
}

}// namespace tudatpy
