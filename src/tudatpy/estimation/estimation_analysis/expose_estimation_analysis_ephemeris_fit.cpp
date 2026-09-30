/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif
#include "expose_estimation_analysis_ephemeris_fit.h"

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "scalarTypes.h"
#include "tudat/simulation/estimation_setup/fitOrbitToEphemeris.h"

namespace py = pybind11;
namespace tss = tudat::simulation_setup;
namespace tep = tudat::estimatable_parameters;
namespace tom = tudat::observation_models;
namespace tp = tudat::propagators;
namespace trf = tudat::reference_frames;

namespace tudatpy
{
namespace estimation
{
namespace estimation_analysis
{

void expose_estimation_analysis_ephemeris_fit( py::module& m )
{
    m.def( "create_best_fit_to_ephemeris",
           &tss::createBestFitToCurrentEphemeris< TIME_TYPE, STATE_SCALAR_TYPE >,
           py::arg( "bodies" ),
           py::arg( "acceleration_models" ),
           py::arg( "observed_bodies" ),
           py::arg( "central_bodies" ),
           py::arg( "integrator_settings" ),
           py::arg( "initial_time" ),
           py::arg( "final_time" ),
           py::arg( "data_point_interval" ),
           py::arg( "additional_parameter_names" ) = std::vector< std::shared_ptr< tep::EstimatableParameterSettings > >( ),
           py::arg( "number_of_iterations" ) = 3,
           py::arg( "reintegrate_variational_equations" ) = true,
           py::arg( "results_print_frequency" ) = 0.0,
           R"doc(

         Fit a numerically propagated orbit to the existing ephemerides of ``observed_bodies``.

         ``bodies``, ``acceleration_models``, ``central_bodies`` and ``integrator_settings`` define the dynamics.
         ``initial_time``, ``final_time`` and ``data_point_interval`` define the fit interval and sampling in seconds
         since J2000. Initial states and any ``additional_parameter_names`` are estimated over ``number_of_iterations``
         iterations. ``reintegrate_variational_equations`` controls derivative updates and ``results_print_frequency``
         controls propagation reporting. Returns EstimationOutput.

      )doc" );
}

}  // namespace estimation_analysis
}  // namespace estimation
}  // namespace tudatpy
