/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#define PYBIND11_DETAILED_ERROR_MESSAGES
#include "expose_estimation_analysis_estimator.h"

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "scalarTypes.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"

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

void expose_estimation_analysis_estimator( py::module &m )
{
    py::class_< tss::OrbitDeterminationManager< STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tss::OrbitDeterminationManager< STATE_SCALAR_TYPE, TIME_TYPE > > >( m,
                                                                                                     "Estimator",
                                                                                                     R"doc(

         Class for consolidating all estimation functionality.

         Class for consolidating all functionality required to perform an estimation.





      )doc" )
            .def( py::init< const tss::SystemOfBodies &,
                            const std::shared_ptr< tep::EstimatableParameterSet< STATE_SCALAR_TYPE > >,
                            const std::vector< std::shared_ptr< tom::ObservationModelSettings > > &,
                            const std::shared_ptr< tp::PropagatorSettings< STATE_SCALAR_TYPE > >,
                            const bool >( ),
                  py::arg( "bodies" ),
                  py::arg( "estimated_parameters" ),
                  py::arg( "observation_settings" ),
                  py::arg( "propagator_settings" ),
                  py::arg( "integrate_on_creation" ) = true,
                  R"doc(

         Class constructor.

         Constructor through which the user can create instances of this class.
         Defines environment, propagation and integrations models, as well as a number of settings related
         to the estimatable parameters and observation settings.

         .. note:: When using default settings, creating an object of
                   this type automatically triggers the propagation


         Parameters
         ----------
         bodies : :class:`~tudatpy.dynamics.environment.SystemOfBodies`
             Object defining the physical environment, with all
             properties of artificial and natural bodies.

         estimated_parameters : :class:`~tudatpy.dynamics.parameters.EstimatableParameterSet`
             Object defining a consolidated set of estimatable parameters,
             linked to the environment and acceleration settings of the simulation.

         observation_settings : :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationModelSettings`
             List of settings objects, each object defining the observation model settings for one
             combination of observable and link geometry that is to be simulated.

         integrator_settings : :class:`~tudatpy.dynamics.propagation_setup.integrator.IntegratorSettings`
             Settings to create the numerical integrator that is to be
             used for the integration of the equations of motion

         propagator_settings : :class:`~tudatpy.dynamics.propagation_setup.propagator.PropagatorSettings`
             Settings to create the propagator that is to be
             used for the propagation of dynamics

         integrate_on_creation : bool, default = True
             Boolean defining whether the propagation should be
             performed immediately (default), or at a later time
             (when calling the :func:`perform_estimation` member function.





     )doc" )
            .def_property_readonly( "observation_simulators",
                                    &tss::OrbitDeterminationManager< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationSimulators,
                                    R"doc(

         **read-only**

         Observation simulators contained in the Estimator object. A single observation simulator hosts
         the functionality for simulating a given observable over the defined link geometry.


         :type: list[ :class:`~tudatpy.estimation.observable_models.observables_simulation.ObservationSimulator` ]
      )doc" )
            .def_property_readonly( "observation_managers",
                                    &tss::OrbitDeterminationManager< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationManagers,
                                    R"doc(

         **read-only**

         Observation managers contained in the Estimator object. A single observation manager can simulate observations and
         calculate observation partials for all link ends involved in the given observable type.


         :type: dict[ :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservableType`, :class:`~tudatpy.estimation.observations.ObservationManager` ]
      )doc" )
            .def_property_readonly(
                    "state_transition_interface",
                    &tss::OrbitDeterminationManager< STATE_SCALAR_TYPE, TIME_TYPE >::getStateTransitionAndSensitivityMatrixInterface,
                    R"doc(

         **read-only**

         State transition and sensitivity matrix interface, in which the numerical solution of the variational equations is stored/updated


         :type: :class:`~tudatpy.dynamics.simulator.CombinedStateTransitionAndSensitivityMatrixInterface`
      )doc" )
            .def( "perform_estimation",
                  &tss::OrbitDeterminationManager< STATE_SCALAR_TYPE, TIME_TYPE >::estimateParameters,
                  py::arg( "estimation_input" ),
                  R"doc(

         Function to trigger the parameter estimation.

         Function to trigger the parameter estimation. Part the process and requirements are similar to those described in the
         :attr:`~tudatpy.estimation.estimation_analysis.Estimator.compute_covariance` function.

         The estimation performs an iterative differential correction of the estimated parameters, where for iteration :math:`i` a correction
         to the parameter vector :math:`\mathbf{p}` is computed according to:

         .. math::

           \Delta\mathbf{p}_{i}&=\mathbf{P}_{i}\left(\mathbf{H}_{i}\mathbf{W}\Delta\mathbf{z}_{i}\right)\\
           \mathbf{p}_{i+1}&=\mathbf{p}_{i}+\Delta\mathbf{p}_{i}

         where :math:`\mathbf{P}` is the covariance (see previous section; where using consider parameters, we have :math:`\mathbf{P}\rightarrow\mathbf{P}^{c}` in the above), and :math:`\Delta\mathbf{z}_{i}` is the observation residual at
         iteration :math:`i`, computed from:

         .. math::

           \Delta\mathbf{z}_{i} = \mathbf{z} - \mathbf{h}(\mathbf{p}_{i})

         with :math:`\mathbf{z}` the vector of all observations provided as input to the data (observed data) and
         :math:`\mathbf{h}(\mathbf{p}_{i})` the vector of all observations, as computed from the current
         estimate of the parameters (computed data). The above procedure is performed iteratively, until convergence has been reached.

         The results of the estimation are stored in an :class:`~tudatpy.estimation.estimation_analysis.EstimationOutput` object
         (which by default contains the results from the iteration where the residual was lowest), with
         two key quantities being:

         * The residual vector of the iteration that had the lowest residual, from the :attr:`~tudatpy.estimation.estimation_analysis.EstimationOutput.final_residuals` attribute of the :class:`~tudatpy.estimation.estimation_analysis.EstimationOutput` class
         * The values of the parameters at the iteration that had the lowest residual, from the :attr:`~tudatpyestimation.estimation_analysis.EstimationOutput.final_parameters` attribute of the :class:`~tudatpy.estimation.estimation_analysis.EstimationOutput` class

         After the estimation is finished, the properties of both the environment (in the ``bodies``) and the estimated parameters
         (in the ``parameters_to_estimate``) are modified as follows:

         * The ephemerides of all propagated/estimated bodies will be set to the propagation results of the last iteration in the estimation. For instance, when estimating the state of body "Delfi-C3", the (tabulated) ephemeris of this body will be set to contain the numerical results of the last iteration of the estimation
         * The values of the parameter values in the ``parameters_to_estimate`` object are those of the last iteration of the estimation. Note that, if the ``apply_final_parameter_correction`` parameter to the :class:`~tudatpy.estimation.estimation_analysis.EstimationInput` is set to ``True``, the parameter correction computed at the end of the last iteration (for which the performance has *not* been computed) has been used to update the parameters vector

         Parameters
         ----------
         estimation_input : :class:`~tudatpy.estimation.estimation_analysis.EstimationInput`
             Object consolidating all relevant settings for the estimation
             This includes foremost the simulated observations, as well as a priori information about the estimatable parameters and convergence criteria for the least squares estimation.

         Returns
         -------
         :class:`~tudatpy.estimation.estimation_analysis.EstimationOutput`
             Object containing all outputs from the estimation process.





     )doc" )
            .def( "compute_covariance",
                  &tss::OrbitDeterminationManager< STATE_SCALAR_TYPE, TIME_TYPE >::computeCovariance,
                  py::arg( "covariance_analysis_input" ),
                  R"doc(

         Function to perform a covariance analysis for the given observations and parameters

         Function to perform a covariance analysis for the given observations and parameters. The observations (which include weights, see TODO)
         are provided through the
         ``covariance_analysis_input`` input and inverse a priori covariance :math:`(\mathbf{P}_{0})^{-1}`.
         Calling this function uses the environment and propagator settings provided to the constructor of this class to simulate
         the dynamics of any relevant bodies for the observations (and associated variational equations for the estimated parameters and states :math:`\mathbf{p}`).
         The observations :math:\mathbf{h}` are then
         computed using the observation models created by the settings provided to the constructor of this ``Estimator`` class, as is the
         associated design matrix :math:`\mathbf{H}(=\frac{\partial\mathbf{h}}{\mathbf{p}})`. This function then produces the covariance :math:`\mathbf{P}`
         (internally using the normalized partials and covariance for stability, see :class:`~tudatpy.estimation.estimation_analysis.CovarianceAnalysisOutput` for details)

         .. math::
            \mathbf{P}=\left(\mathbf{H}^{T}\mathbf{W}\mathbf{H}+(\mathbf{P}_{0})^{-1}\right)^{-1}

         In the presence of consider parameters, an additional term :math:`\Delta\mathbf{P}_{c}` is computed that denotes the contribution
         of the consider covariance to the parameter uncertainties, which is computed from the above as:

         .. math::
           \Delta \mathbf{P}^{c} = \left(\mathbf{P}\mathbf{H}^{T}\mathbf{W}\right)\left(\mathbf{H}_{c}\mathbf{C}\mathbf{H}_{c}^{T}\right)\left(\mathbf{P}\mathbf{H}^{T}\mathbf{W}\right)^{T}
           \mathbf{P}^{c}=\mathbf{P}+ \Delta \mathbf{P}^{c}

         where :math:`\mathbf{H}_{c}` is the design matrix for the consider parameters.

         Note that, although the actual observations are formally not required for a covariance analysis, all additional data (e.g. observation time, type, link ends, etc.)
         are. And, as such, the ``covariance_analysis_input`` does require the full set of observations and associated information, for consistency purposes (e.g., same input as
         ``perform_estimation`` function) .

         To propagate the covariance of the initial states obtained from this function to other epochs, use the :func:`~tudatpy.estimation.estimation_analysis.propagate_covariance` function


         Parameters
         ----------
         covariance_analysis_input : :class:`~tudatpy.estimation.estimation_analysis.CovarianceAnalysisInput`
             Object consolidating all relevant settings for the covariance analysis
             This includes foremost the simulated observations, as well as a priori information about the estimatable parameters

         Returns
         -------
         :class:`~tudatpy.estimation.estimation_analysis.CovarianceAnalysisOutput`
             Object containing all outputs from the estimation process.





     )doc" )
            .def_property_readonly( "variational_solver",
                                    &tss::OrbitDeterminationManager< STATE_SCALAR_TYPE, TIME_TYPE >::getVariationalEquationsSolver,
                                    R"doc(

         **read-only**

         Variational equations solver, which is used to manage and execute the numerical integration of
         equations of motion and variational equations/dynamics in the Estimator object.

         :type: :class:`~tudatpy.dynamics.simulator.SingleArcVariationalSimulator`
      )doc" );
}

}  // namespace estimation_analysis
}  // namespace estimation
}  // namespace tudatpy
