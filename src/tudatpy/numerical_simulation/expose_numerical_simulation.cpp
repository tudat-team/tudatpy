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

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "scalarTypes.h"
#include "tudat/astro/basic_astro/dateTime.h"
#include "tudat/basics/timeType.h"
#include "tudat/simulation/environment_setup.h"
#include "tudat/simulation/estimation_setup.h"
#include "tudat/simulation/propagation_setup.h"

namespace py = pybind11;
namespace tp = tudat::propagators;
namespace tss = tudat::simulation_setup;
namespace tni = tudat::numerical_integrators;
namespace tep = tudat::estimatable_parameters;
namespace tom = tudat::observation_models;
namespace tba = tudat::basic_astrodynamics;

namespace tudatpy
{

namespace numerical_simulation
{

PYBIND11_MODULE( expose_numerical_simulation, m )
{
    // Additional imports to solve issues with submodules
    // Expose IntegratorSettings locally to avoid circular imports
    py::class_< tni::IntegratorSettings< TIME_TYPE >,
                std::shared_ptr< tni::IntegratorSettings< TIME_TYPE > > >(
            m,
            "IntegratorSettings",
            py::module_local( ),
            R"doc(Functional base class to define settings for integrators.

             Class to define settings for numerical integrators, for instance for use in
             numerical integration of equations of motion/variational equations. This class can
             be used for simple integrators such as fixed step RK and Euler. Integrators that
             require more settings to define have their own derived class.)doc" )
            .def_readwrite( "initial_time",
                            &tni::IntegratorSettings< TIME_TYPE >::initialTimeDeprecated_ );
    /////////////////////////////////////////////////////////////////////////

    py::class_< tudat::Time >( m, "Time", R"doc(No documentation found.)doc" )
            .def( py::init< const int, const long double >( ),
                  py::arg( "full_periods" ),
                  py::arg( "seconds_into_full_period" ) )
            .def( py::init< const double >( ), py::arg( "seconds_into_full_period" ) )
            .def( "to_float",
                  &tudat::Time::getSeconds< double >,
                  R"doc(No documentation found.)doc" )
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

    m.def( "get_integrated_type_and_body_list",
           &tp::getIntegratedTypeAndBodyList< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "propagator_settings" ) );

    m.def( "get_single_integration_size", &tp::getSingleIntegrationSize, py::arg( "state_type" ) );

    // Estimator, Simulator & Variational
    expose_numerical_simulation_estimator( m );
    expose_numerical_simulation_simulator( m );
    expose_numerical_simulation_variational( m );
};

}  // namespace numerical_simulation
}  // namespace tudatpy
