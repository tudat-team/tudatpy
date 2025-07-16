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
#include "expose_estimation_refactoring.h"

#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "scalarTypes.h"
#include "tudat/astro/basic_astro/dateTime.h"
#include "tudat/basics/timeType.h"

namespace py = pybind11;
// namespace tp = tudat::propagators;
// namespace tss = tudat::simulation_setup;
// namespace tni = tudat::numerical_integrators;
// namespace tep = tudat::estimatable_parameters;
// namespace tom = tudat::observation_models;
// namespace tba = tudat::basic_astrodynamics;

namespace tudatpy
{

namespace estimation_refactoring
{

void expose_estimation_refactoring( py::module &m )
{
    auto observable_models_setup_submodule = m.def_submodule( "observable_models_setup" );
    observable_models_setup::expose_observable_models_setup( observable_models_setup_submodule );

    auto observable_models_submodule = m.def_submodule( "observable_models" );
    observable_models::expose_observable_models( observable_models_submodule );

    auto observations_setup_submodule = m.def_submodule( "observations_setup" );
    observations_setup::expose_observations_setup( observations_setup_submodule );

    auto observations_submodule = m.def_submodule( "observations" );
    observations::expose_observations( observations_submodule );

    auto estimation_analysis_submodule = m.def_submodule( "estimation_analysis" );
    estimation_analysis::expose_estimation_analysis( estimation_analysis_submodule );

    // auto estimation_submodule = m.def_submodule( "estimation" );
    // estimation::expose_estimation_filter_parser( estimation_submodule );
    // estimation::expose_estimation( estimation_submodule );
    // estimation::expose_estimation_observation_collection( estimation_submodule );
    // estimation::expose_estimation_propagated_covariance( estimation_submodule );
    // estimation::expose_estimation_single_observation_set( estimation_submodule );


    // py::class_< tudat::Time >( m, "Time", R"doc(No documentation found.)doc" )
    //         .def( py::init< const int, const long double >( ),
    //               py::arg( "full_periods" ),
    //               py::arg( "seconds_into_full_period" ) )
    //         .def( py::init< const double >( ), py::arg( "seconds_since_j2000" ) )
    //         .def( "to_float",
    //               &tudat::Time::getSeconds< double >,
    //               R"doc(No documentation found.)doc" )
    //         .def( py::self + py::self )
    //         .def( py::self + double( ) )
    //         .def( double( ) + py::self )
    //         .def( py::self += py::self )
    //         .def( py::self += double( ) )
    //         .def( py::self - py::self )
    //         .def( py::self - double( ) )
    //         .def( py::self -= py::self )
    //         .def( py::self -= double( ) )
    //         .def( double( ) - py::self )
    //         .def( py::self * double( ) )
    //         .def( double( ) * py::self )
    //         .def( py::self *= double( ) )
    //         .def( py::self / double( ) )
    //         .def( py::self /= double( ) )
    //         .def( py::self == py::self )
    //         .def( double( ) == py::self )
    //         .def( py::self == double( ) )
    //         .def( py::self != py::self )
    //         .def( py::self != double( ) )
    //         .def( double( ) != py::self )
    //         .def( py::self < py::self )
    //         .def( py::self < double( ) )
    //         .def( double( ) < py::self )
    //         .def( py::self > py::self )
    //         .def( py::self > double( ) )
    //         .def( double( ) > py::self )
    //         .def( py::self <= py::self )
    //         .def( py::self <= double( ) )
    //         .def( double( ) <= py::self )
    //         .def( py::self >= py::self )
    //         .def( double( ) >= py::self )
    //         .def( py::self >= double( ) );

    // m.def( "get_integrated_type_and_body_list",
    //        &tp::getIntegratedTypeAndBodyList< STATE_SCALAR_TYPE, TIME_TYPE >,
    //        py::arg( "propagator_settings" ) );

    // m.def( "get_single_integration_size", &tp::getSingleIntegrationSize, py::arg( "state_type" ) );
};

}  // namespace estimation_refactoring
}  // namespace tudatpy
