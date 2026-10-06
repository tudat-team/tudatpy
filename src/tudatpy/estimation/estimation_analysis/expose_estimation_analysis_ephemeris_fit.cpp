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
         controls propagation reporting.

         Parameters
         ----------
         bodies : SystemOfBodies
             System of bodies defining the physical environment.
         acceleration_models : dict[str, dict[str, list[AccelerationModel]]]
             Acceleration models used for the propagated ephemeris fit.
         observed_bodies : list[str]
             Names of the bodies whose ephemerides provide the observations.
         central_bodies : list[str]
             Names of the reference bodies, paired with the propagated or observed bodies.
         integrator_settings : IntegratorSettings
             Numerical integration settings for the fit or sampling procedure.
         initial_time : Time
             Initial epoch, in seconds since J2000.
         final_time : Time
             Final epoch, in seconds since J2000.
         data_point_interval : Time
             Interval between ephemeris observations used for the fit, in seconds.
         additional_parameter_names : list[EstimatableParameterSettings], optional
             Settings for additional parameters to estimate alongside the initial states.
         number_of_iterations : int, optional
             Number of estimation iterations.
         reintegrate_variational_equations : bool, optional
             Whether to reintegrate the variational equations at each iteration.
         results_print_frequency : float, optional
             Interval in seconds between printed propagation results.

         Returns
         -------
         EstimationOutput
             Estimation results of the dynamical fit to the ephemeris observations.

      )doc" );
}

}  // namespace estimation_analysis
}  // namespace estimation
}  // namespace tudatpy
