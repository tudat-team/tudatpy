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
#include "expose_parameters.h"

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "scalarTypes.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/initialTranslationalState.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/observationBiasParameter.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/sphericalHarmonicCosineCoefficients.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/sphericalHarmonicSineCoefficients.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/tidalLoveNumber.h"
#include "tudat/simulation/estimation_setup/createEstimatableParametersFactory.h"

namespace py = pybind11;
namespace tep = tudat::estimatable_parameters;

namespace tudatpy
{
namespace dynamics
{
namespace parameters
{

void expose_parameters( py::module& m )
{
    py::class_< tep::EstimatableParameterBase, std::shared_ptr< tep::EstimatableParameterBase > >( m,
                                                                                                   "EstimatableParameter",
                                                                                                   R"doc(

         Base class for all estimatable parameters exposed through an untemplated vector-valued interface.

      )doc" )
            .def_property_readonly( "parameter_identifier",
                                    &tep::EstimatableParameterBase::getParameterName,
                                    R"doc(

         **read-only**

         Full parameter identifier as ``(parameter_type, (associated_body, secondary_identifier))``.

         :type: tuple[ :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterTypes`, tuple[str, str] ]
      )doc" )
            .def_property_readonly( "parameter_description",
                                    &tep::EstimatableParameterBase::getParameterDescription,
                                    R"doc(

         **read-only**

         Human-readable parameter description.

         :type: str
      )doc" )
            .def_property_readonly( "parameter_size",
                                    &tep::EstimatableParameterBase::getParameterSize,
                                    R"doc(

         **read-only**

         Size of this parameter block in the estimation vector.

         :type: int
      )doc" )
            .def_property_readonly( "time_intervals_of_influence",
                                    &tep::EstimatableParameterBase::getTimeIntervalsOfInfluencePerEntry,
                                    R"doc(

         **read-only**

         Time interval ``(start, end)`` in which each entry of the parameter can influence the dynamics, in the order of the
         parameter value. Empty for global parameters, which influence the dynamics at all times.
         Arc-wise dynamical parameters (e.g. from :func:`~tudatpy.dynamics.parameters_setup.arcwise_radiation_pressure_coefficient`,
         :func:`~tudatpy.dynamics.parameters_setup.arcwise_empirical_accelerations` or
         :func:`~tudatpy.dynamics.parameters_setup.arcwise_constant_drag_coefficient`) give the interval of the arc of each entry,
         where the last arc ends at an arbitrarily large placeholder time. Desaturation maneuver delta V parameters give
         the thrust interval of the maneuver of each entry.

         :type: list[tuple[float, float]]
      )doc" )
            .def( "get_parameter_value",
                  &tep::EstimatableParameterBase::getParameterValueBase,
                  R"doc(

         Get the parameter value through the untemplated interface.

         Returns
         -------
         numpy.ndarray[numpy.float64[m, 1]]
             Current value of this parameter block.

      )doc" )
            .def( "set_parameter_value",
                  &tep::EstimatableParameterBase::setParameterValueBase,
                  py::arg( "parameter_value" ),
                  R"doc(

         Set the parameter value through the untemplated interface.

         Parameters
         ----------
         parameter_value : numpy.ndarray[numpy.float64[m, 1]]
             New value for this parameter block.

      )doc" );

    py::class_< tep::SphericalHarmonicsCosineCoefficients,
                std::shared_ptr< tep::SphericalHarmonicsCosineCoefficients >,
                tep::EstimatableParameterBase >( m,
                                                 "SphericalHarmonicsCosineCoefficientsParameter",
                                                 R"doc(

         Estimatable parameter of a block of spherical harmonic cosine coefficients of a body's gravity field.

         Created from the settings of :func:`~tudatpy.dynamics.parameters_setup.spherical_harmonics_c_coefficients`
         or :func:`~tudatpy.dynamics.parameters_setup.spherical_harmonics_c_coefficients_block`.
         Entry ``i`` of the parameter value is the cosine coefficient :math:`\bar{C}_{lm}` (geodesy normalized) at the
         degree and order ``block_indices[i]``. For a gravity field with gravity field variations, these are the nominal coefficients.

      )doc" )
            .def_property_readonly( "block_indices",
                                    &tep::SphericalHarmonicsCosineCoefficients::getBlockIndices,
                                    R"doc(

         **read-only**

         Degree and order :math:`(l, m)` of each entry of the parameter, sorted in the order of the parameter value. For the settings of
         :func:`~tudatpy.dynamics.parameters_setup.spherical_harmonics_c_coefficients_block`, they are sorted in the same order of the 
         block indices given by the user.

         :type: list[tuple[int, int]]
      )doc" );

    py::class_< tep::SphericalHarmonicsSineCoefficients,
                std::shared_ptr< tep::SphericalHarmonicsSineCoefficients >,
                tep::EstimatableParameterBase >( m,
                                                 "SphericalHarmonicsSineCoefficientsParameter",
                                                 R"doc(

         Estimatable parameter of a block of spherical harmonic sine coefficients of a body's gravity field.

         Created from the settings of :func:`~tudatpy.dynamics.parameters_setup.spherical_harmonics_s_coefficients`
         or :func:`~tudatpy.dynamics.parameters_setup.spherical_harmonics_s_coefficients_block`.
         Entry ``i`` of the parameter value is the sine coefficient :math:`\bar{S}_{lm}` (geodesy normalized) at the
         degree and order ``block_indices[i]``. For a gravity field with gravity field variations, these are the nominal coefficients.

      )doc" )
            .def_property_readonly( "block_indices",
                                    &tep::SphericalHarmonicsSineCoefficients::getBlockIndices,
                                    R"doc(

         **read-only**

         Degree and order :math:`(l, m)` of each entry of the parameter, sorted in the order of the parameter value. For the settings of
         :func:`~tudatpy.dynamics.parameters_setup.spherical_harmonics_s_coefficients_block`, they are sorted in the same order of the 
         block indices given by the user.

         :type: list[tuple[int, int]]
      )doc" );

    py::class_< tep::TidalLoveNumber< Eigen::VectorXd >,
                std::shared_ptr< tep::TidalLoveNumber< Eigen::VectorXd > >,
                tep::EstimatableParameterBase >( m,
                                                 "TidalLoveNumberParameter",
                                                 R"doc(

         Base class of the estimatable tidal Love numbers :math:`k_{l}` or :math:`k_{lm}` at a single degree of a body's
         solid body tide.

         Parameters retrieved from a parameter set are instances of one of its derived classes:
         :class:`~tudatpy.dynamics.parameters.FullDegreeTidalLoveNumberParameter` or
         :class:`~tudatpy.dynamics.parameters.SingleDegreeVariableTidalLoveNumberParameter`.

      )doc" )
            .def_property_readonly( "degree",
                                    &tep::TidalLoveNumber< Eigen::VectorXd >::getDegree,
                                    R"doc(

         **read-only**

         Degree :math:`l` of the Love number(s).

         :type: int
      )doc" )
            .def_property_readonly( "orders",
                                    &tep::TidalLoveNumber< Eigen::VectorXd >::getOrders,
                                    R"doc(

         **read-only**

         Orders :math:`m` of the Love number(s). For a :class:`~tudatpy.dynamics.parameters.SingleDegreeVariableTidalLoveNumberParameter`,
         the orders that are estimated separately, in the order of the parameter value. For a
         :class:`~tudatpy.dynamics.parameters.FullDegreeTidalLoveNumberParameter`, all orders :math:`0..l`, to which its single Love
         number applies.

         :type: list[int]
      )doc" )
            .def_property_readonly( "use_complex_love_number",
                                    &tep::TidalLoveNumber< Eigen::VectorXd >::useComplexComponents,
                                    R"doc(

         **read-only**

         Whether the real and imaginary part of each Love number are estimated (True), or only the real part (False).

         :type: bool
      )doc" )
            .def_property_readonly( "deforming_bodies",
                                    &tep::TidalLoveNumber< Eigen::VectorXd >::getDeformingBodies,
                                    R"doc(

         **read-only**

         Bodies raising the tides of which the Love number(s) are estimated.

         :type: list[str]
      )doc" );

    py::class_< tep::FullDegreeTidalLoveNumber,
                std::shared_ptr< tep::FullDegreeTidalLoveNumber >,
                tep::TidalLoveNumber< Eigen::VectorXd > >( m,
                                                           "FullDegreeTidalLoveNumberParameter",
                                                           R"doc(

         Estimatable tidal Love number :math:`k_{l}` at a single degree, equal for all orders.

         Created from the settings of :func:`~tudatpy.dynamics.parameters_setup.order_invariant_k_love_number`.
         The parameter value is :math:`[k_{l}]`, or :math:`[\text{Re}(k_{l}), \text{Im}(k_{l})]` if
         :attr:`~tudatpy.dynamics.parameters.TidalLoveNumberParameter.use_complex_love_number` is True.

      )doc" );

    py::class_< tep::SingleDegreeVariableTidalLoveNumber,
                std::shared_ptr< tep::SingleDegreeVariableTidalLoveNumber >,
                tep::TidalLoveNumber< Eigen::VectorXd > >( m,
                                                           "SingleDegreeVariableTidalLoveNumberParameter",
                                                           R"doc(

         Estimatable tidal Love numbers :math:`k_{lm}` at a single degree, separate for each order.

         Created from the settings of :func:`~tudatpy.dynamics.parameters_setup.order_varying_k_love_number`.
         The parameter value has one entry :math:`k_{lm}` per order :math:`m` in
         :attr:`~tudatpy.dynamics.parameters.TidalLoveNumberParameter.orders`, in that order. If
         :attr:`~tudatpy.dynamics.parameters.TidalLoveNumberParameter.use_complex_love_number` is True, it has two entries per order instead,
         :math:`\text{Re}(k_{lm})` followed by :math:`\text{Im}(k_{lm})`.

      )doc" );

    py::class_< tep::InitialTranslationalStateParameter< STATE_SCALAR_TYPE >,
                std::shared_ptr< tep::InitialTranslationalStateParameter< STATE_SCALAR_TYPE > >,
                tep::EstimatableParameterBase >( m,
                                                 "InitialTranslationalStateParameter",
                                                 R"doc(

         Estimatable parameter of the initial translational state of a body in a single-arc propagation.

         Created from the settings of :func:`~tudatpy.dynamics.parameters_setup.initial_states` for a single-arc propagator.
         The parameter value is the Cartesian state :math:`(x, y, z, v_{x}, v_{y}, v_{z})` w.r.t.
         :attr:`~tudatpy.dynamics.parameters.InitialTranslationalStateParameter.central_body`, at the initial time of the
         propagation.

      )doc" )
            .def_property_readonly( "central_body",
                                    &tep::InitialTranslationalStateParameter< STATE_SCALAR_TYPE >::getCentralBody,
                                    R"doc(

         **read-only**

         Body w.r.t. which the initial state is defined.

         :type: str
      )doc" )
            .def_property_readonly( "frame_orientation",
                                    &tep::InitialTranslationalStateParameter< STATE_SCALAR_TYPE >::getFrameOrientation,
                                    R"doc(

         **read-only**

         Orientation of the frame in which the initial state is defined.

         :type: str
      )doc" );

    py::class_< tep::ArcWiseInitialTranslationalStateParameter< STATE_SCALAR_TYPE >,
                std::shared_ptr< tep::ArcWiseInitialTranslationalStateParameter< STATE_SCALAR_TYPE > >,
                tep::EstimatableParameterBase >( m,
                                                 "ArcWiseInitialTranslationalStateParameter",
                                                 R"doc(

         Estimatable parameter of the initial translational states of a body in the arcs of a multi-arc propagation.

         Created from the settings of :func:`~tudatpy.dynamics.parameters_setup.initial_states` for a multi-arc (or the multi-arc
         part of a hybrid-arc) propagator. The parameter value has six entries per arc, in the order of the arcs: the Cartesian
         state :math:`(x, y, z, v_{x}, v_{y}, v_{z})` of arc ``i`` w.r.t. ``central_bodies[i]``, at
         ``arc_start_times[i]``.

      )doc" )
            .def_property_readonly( "arc_start_times",
                                    &tep::ArcWiseInitialTranslationalStateParameter< STATE_SCALAR_TYPE >::getArcStartTimes,
                                    R"doc(

         **read-only**

         Start times of the arcs, at which the initial states are defined.

         :type: list[float]
      )doc" )
            .def_property_readonly( "central_bodies",
                                    &tep::ArcWiseInitialTranslationalStateParameter< STATE_SCALAR_TYPE >::getCentralBodies,
                                    R"doc(

         **read-only**

         Body w.r.t. which the initial state of each arc is defined.

         :type: list[str]
      )doc" )
            .def_property_readonly( "frame_orientation",
                                    &tep::ArcWiseInitialTranslationalStateParameter< STATE_SCALAR_TYPE >::getFrameOrientation,
                                    R"doc(

         **read-only**

         Orientation of the frame in which the initial states are defined.

         :type: str
      )doc" );

    py::class_< tep::SingleArcObservationBiasParameter,
                std::shared_ptr< tep::SingleArcObservationBiasParameter >,
                tep::EstimatableParameterBase >( m,
                                                 "SingleArcObservationBiasParameter",
                                                 R"doc(

         Estimatable parameter of a constant observation bias of one observable and link.

         Created from the settings of :func:`~tudatpy.dynamics.parameters_setup.absolute_observation_bias` or
         :func:`~tudatpy.dynamics.parameters_setup.relative_observation_bias`. The parameter value has one entry per
         component of the observable.

      )doc" )
            .def_property_readonly( "link_ends",
                                    &tep::SingleArcObservationBiasParameter::getLinkEnds,
                                    R"doc(

         **read-only**

         Link ends of the observations to which the bias applies.

         :type: dict[:class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType`, :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndId`]
      )doc" )
            .def_property_readonly( "observable_type",
                                    &tep::SingleArcObservationBiasParameter::getObservableType,
                                    R"doc(

         **read-only**

         Observable type to which the bias applies.

         :type: :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservableType`
      )doc" );

    py::class_< tep::MultiArcObservationBiasParameter,
                std::shared_ptr< tep::MultiArcObservationBiasParameter >,
                tep::EstimatableParameterBase >( m,
                                                 "MultiArcObservationBiasParameter",
                                                 R"doc(

         Estimatable parameter of arc-wise constant observation biases of one observable and link.

         Created from the settings of :func:`~tudatpy.dynamics.parameters_setup.arcwise_absolute_observation_bias` or
         :func:`~tudatpy.dynamics.parameters_setup.arcwise_relative_observation_bias`. The parameter value has one entry per
         component of the observable for each arc, in the order of the arcs: entries ``i * n`` to ``(i + 1) * n - 1`` are the
         bias of the arc starting at ``arc_start_times[i]``, for an observable of size ``n``.

      )doc" )
            .def_property_readonly( "link_ends",
                                    &tep::MultiArcObservationBiasParameter::getLinkEnds,
                                    R"doc(

         **read-only**

         Link ends of the observations to which the biases apply.

         :type: dict[:class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType`, :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndId`]
      )doc" )
            .def_property_readonly( "observable_type",
                                    &tep::MultiArcObservationBiasParameter::getObservableType,
                                    R"doc(

         **read-only**

         Observable type to which the biases apply.

         :type: :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservableType`
      )doc" )
            .def_property_readonly( "arc_start_times",
                                    &tep::MultiArcObservationBiasParameter::getArcStartTimes,
                                    R"doc(

         **read-only**

         Start times of the bias arcs. An observation belongs to the arc in which the time of the link end given by
         :attr:`~tudatpy.dynamics.parameters.MultiArcObservationBiasParameter.link_end_index` falls.

         :type: list[float]
      )doc" )
            .def_property_readonly( "link_end_index",
                                    &tep::MultiArcObservationBiasParameter::getLinkEndIndex,
                                    R"doc(

         **read-only**

         Index, in the link-end order of the observable, of the link end whose time determines the arc of an observation (the
         ``time_link_end`` of the settings).

         :type: int
      )doc" );

    // ESTIMATABLE PARAMETER SET CLASS

    py::class_< tep::EstimatableParameterSet< STATE_SCALAR_TYPE >, std::shared_ptr< tep::EstimatableParameterSet< STATE_SCALAR_TYPE > > >(
            m,
            "EstimatableParameterSet",
            R"doc(

         Class containing a consolidated set of estimatable parameters.

         Class containing a consolidated set of estimatable parameters, linked to the environment and acceleration settings of the simulation.
         The user typically creates instances of this class via the :func:`~tudatpy.dynamics.parameters_setup.create_parameter_set` function.





      )doc" )
            .def_property_readonly( "parameter_set_size",
                                    &tep::EstimatableParameterSet< STATE_SCALAR_TYPE >::getEstimatedParameterSetSize,
                                    R"doc(

         **read-only**

         Size of the parameter set, i.e. amount of estimatable parameters contained in the set.

         :type: int
      )doc" )
            .def_property_readonly( "initial_states_size",
                                    &tep::EstimatableParameterSet< STATE_SCALAR_TYPE >::getInitialDynamicalStateParameterSize,
                                    R"doc(

         **read-only**

         Amount of initial state parameters contained in the set.

         :type: int
      )doc" )
            .def_property_readonly( "initial_single_arc_states_size",
                                    &tep::EstimatableParameterSet< STATE_SCALAR_TYPE >::getInitialDynamicalSingleArcStateParameterSize,
                                    R"doc(

         **read-only**

         Amount of initial state parameters in the set, which are treated in a single-arc fashion.

         :type: int
      )doc" )
            .def_property_readonly( "initial_multi_arc_states_size",
                                    &tep::EstimatableParameterSet< STATE_SCALAR_TYPE >::getInitialDynamicalMultiArcStateParameterSize,
                                    R"doc(

         **read-only**

         Amount of initial state parameters in the set, which are treated in a multi-arc fashion.

         :type: int
      )doc" )
            .def_property_readonly( "constraints_size",
                                    &tep::EstimatableParameterSet< STATE_SCALAR_TYPE >::getConstraintSize,
                                    R"doc(

         **read-only**

         Total size of linear constraint that is to be applied during estimation.

         :type: int
      )doc" )
            .def_property(
                    "parameter_vector",
                    &tep::EstimatableParameterSet< STATE_SCALAR_TYPE >::getFullParameterValuesWithoutConsiderParameters<
                            STATE_SCALAR_TYPE >,
                    &tep::EstimatableParameterSet< STATE_SCALAR_TYPE >::resetParameterValuesWithoutConsiderParameters< STATE_SCALAR_TYPE >,
                    R"doc(

         Vector containing the parameter values of all parameters in the set (excluding consider parameters).

         :type: numpy.ndarray[m, 1] with the configured state scalar dtype (float64 or longdouble)
      )doc" )
            .def_property(
                    "estimated_and_consider_parameter_vector",
                    &tep::EstimatableParameterSet< STATE_SCALAR_TYPE >::getFullParameterValuesWithConsiderParameters< STATE_SCALAR_TYPE >,
                    &tep::EstimatableParameterSet< STATE_SCALAR_TYPE >::resetParameterValuesWithConsiderParameters< STATE_SCALAR_TYPE >,
                    R"doc(

         Vector containing the parameter values of all parameters in the set (including consider parameters; concatenated after estimated parameters).

         :type: numpy.ndarray[m, 1] with the configured state scalar dtype (float64 or longdouble)
      )doc" )
            .def_property( "reduce_sensitivity_matrix_columns",
                           &tep::EstimatableParameterSet< STATE_SCALAR_TYPE >::getReduceSensitivityMatrixColumns,
                           &tep::EstimatableParameterSet< STATE_SCALAR_TYPE >::setReduceSensitivityMatrixColumns,
                           R"doc(

         Boolean denoting whether sensitivity matrix columns that are zero in a propagation arc are omitted when propagating the
         variational equations (default false). If true, parameters that cannot influence the dynamics (e.g. observation biases), and
         parameter entries that only influence the dynamics outside of a propagation arc (e.g. arc-wise empirical accelerations in other
         arcs) are not propagated, which can substantially reduce the computational cost and memory use for large (multi-arc)
         estimations. The time interval of each propagation arc is only used if it is known before the propagation (time-based
         termination). The full sensitivity matrix is reconstructed when it is retrieved. This setting must be set before the
         variational equations solver (or estimator) is created. The propagated columns are selected at that point: if a propagation
         time interval used for this selection is extended afterwards (by changing the initial time or termination settings), or is
         no longer known before the propagation, an error is thrown when propagating the variational equations.

         :type: bool
      )doc" )
            .def_property_readonly( "consider_parameters",
                                    &tep::EstimatableParameterSet< STATE_SCALAR_TYPE >::getConsiderParameters,
                                    R"doc(
                                   
        **read-only**

        Set of consider parameters that are included in the parameter set.

        :type: :class:`~tudatpy.dynamics.parameters.EstimatableParameterSet`

                                   )doc" )
            .def( "indices_for_parameter_type",
                  &tep::EstimatableParameterSet< STATE_SCALAR_TYPE >::getIndicesForParameterType,
                  py::arg( "parameter_type" ),
                  R"doc(

         Retrieve indices of parameter blocks matching a parameter identifier.

         Matching rules:

         * If either identifier string is non-empty, exact matching is used on the full identifier.
         * If both strings are empty, matching is performed on the enum only.
         * If multiple parameter blocks match, all matches are returned.
         * If no parameter blocks match, an empty list is returned.

         Parameters
         ----------
         parameter_type : tuple[ :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterTypes`, tuple[str, str] ]
             Parameter identifier for which matching indices are retrieved.

         Returns
         -------
         list[tuple[int, int]]
             Matching ``(start_index, size)`` entries in parameter-set order.

     )doc" )
            .def( "parameters_for_parameter_type",
                  &tep::EstimatableParameterSet< STATE_SCALAR_TYPE >::getParametersForParameterIdentifier,
                  py::arg( "parameter_type" ),
                  R"doc(

         Retrieve parameter objects matching a parameter identifier.

         Matching rules are identical to
         :meth:`~tudatpy.dynamics.parameters.EstimatableParameterSet.indices_for_parameter_type`.

         Returned parameter objects are instances of
         :class:`~tudatpy.dynamics.parameters.EstimatableParameter`, or of its derived class for the parameter type, if one is
         exposed (e.g. :class:`~tudatpy.dynamics.parameters.SphericalHarmonicsCosineCoefficientsParameter`).

         Parameters
         ----------
         parameter_type : tuple[ :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterTypes`, tuple[str, str] ]
             Parameter identifier for which matching parameter objects are retrieved.

         Returns
         -------
         List[ :class:`~tudatpy.dynamics.parameters.EstimatableParameter` ]
             Matching parameter objects in parameter-set order, as their derived class where one is exposed.

      )doc" )
            .def( "indices_for_parameter_description",
                  &tep::EstimatableParameterSet< STATE_SCALAR_TYPE >::getIndicesForParameterDescription,
                  py::arg( "parameter_description" ),
                  R"doc(

         Function to retrieve the indices of a given type of parameter.

         Function to retrieve the index of all parameters of a given type from the parameter set.
         This function can be very useful, since the order of parameters within the parameter set does not necessarily correspond to the order in which the elements were added to the set.


         Parameters
         ----------
         parameter_description : str
                Parameter description that is to be searched in full list of parameters.
                Parameter descriptions can be retrieved from the :meth:`tudatpy.dynamics.parameters.EstimatableParameterSet.get_parameter_descriptions` method.

         Returns
         -------
         tuple[int, int]
             Indices of the parameters corresponding to the description. The first element of the tuple is the start index, the second element is the size of the parameter.

     )doc" )
            .def( "get_parameter_descriptions",
                  &tep::EstimatableParameterSet< STATE_SCALAR_TYPE >::getParametersDescriptions,
                  R"doc(

         Function to retrieve the descriptions of all parameters in the set.

         Returns
         -------
         List[ str ]
             List of parameter descriptions formatted as strings.


         )doc" )
            .def( "get_parameter_identifiers",
                  &tep::EstimatableParameterSet< STATE_SCALAR_TYPE >::getParametersIdentifiers,
                  R"doc(

         Function to retrieve the identifiers of all parameters in the set.

         Returns
         -------
         list[tuple[:class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterTypes`, tuple[str, str]]]
             List of parameter identifiers. The first element of the tuple is the parameter type, the second element is a tuple containing the body name and, if applicable, a secondary identifier (e.g., station name), else this is an empty string.


         )doc" );

    m.def( "print_parameter_names",
           &tep::printEstimatableParameterEntries< double >,
           py::arg( "parameter_set" ),
           R"doc(

 Function for printing a list of estimatable parameter names.

 Function that allows you to print a verbose list of all parameters that shall be estimated. Consider parameters are listed separately.


 Parameters
 ----------
 parameter_set : :class:`~tudatpy.dynamics.parameters.EstimatableParameterSet`.
     Instance of :class:`~tudatpy.dynamics.parameters.EstimatableParameterSet` class, consolidating all estimatable parameters and simulation models.

 Returns
 -------
 None
     This function prints the parameter names to standard output.

 Examples
 --------
 .. code-block:: python

     import numpy as np
     from tudatpy.interface import spice
     from tudatpy.dynamics import environment_setup, propagation_setup, parameters_setup, parameters
     from tudatpy
     from tudatpy.astro.time_representation import DateTime

     # Load SPICE kernels
     spice.load_standard_kernels()

     # Set simulation epochs
     simulation_start_epoch = DateTime(2000, 1, 1).epoch()
     simulation_end_epoch = DateTime(2000, 1, 4).epoch()

     # Create bodies
     bodies_to_create = ["Sun", "Earth"]
     global_frame_origin = "Earth"
     global_frame_orientation = "J2000"
     body_settings = environment_setup.get_default_body_settings(
       bodies_to_create, global_frame_origin, global_frame_orientation)

     # Create vehicle
     body_settings.add_empty_settings("Delfi-C3")
     body_settings.get("Delfi-C3").constant_mass = 2.2
     bodies = environment_setup.create_system_of_bodies(body_settings)

     # Define propagation settings
     bodies_to_propagate = ["Delfi-C3"]
     central_bodies = ["Earth"]
     accelerations_settings_delfi_c3 = dict(
       Sun=[propagation_setup.acceleration.point_mass_gravity()],
       Earth=[propagation_setup.acceleration.spherical_harmonic_gravity(8, 8)]
     )
     acceleration_settings = {"Delfi-C3": accelerations_settings_delfi_c3}
     acceleration_models = propagation_setup.create_acceleration_models(
       bodies, acceleration_settings, bodies_to_propagate, central_bodies)
     initial_state = np.zeros(6)  # Use a real initial state if needed

     # Integrator settings
     integrator_settings = propagation_setup.integrator.runge_kutta_fixed_step(60.0,
                                                                                    coefficient_set=propagation_setup.integrator.CoefficientSets.rkdp_87)

     # Create propagator
     termination_condition = propagation_setup.propagator.time_termination(simulation_end_epoch)
     propagator_settings = propagation_setup.propagator.translational(
       central_bodies, acceleration_models, bodies_to_propagate, initial_state,
       simulation_start_epoch, integrator_settings, termination_condition)

     # Define parameters to estimate
     parameter_settings = parameters_setup.initial_states(propagator_settings, bodies)
     parameter_settings.append(parameters_setup.gravitational_parameter("Earth"))
     parameters_to_estimate = parameters_setup.create_parameter_set(parameter_settings, bodies)

     # Print parameter names
     print(parameters.print_parameter_names(parameters_to_estimate))




     )doc" );

    m.def( "set_parameter_correction_constraint",
           &tep::setParameterCorrectionConstraint,
           py::arg( "parameter" ),
           py::arg( "constraint_direction" ),
           R"doc(

 Function to constrain the correction to an estimated parameter to be orthogonal to a given direction.

 Function to constrain the correction :math:`\Delta\mathbf{x}` to an estimated parameter to be orthogonal to a given direction
 :math:`\mathbf{g}`, by imposing :math:`\mathbf{g}^{T}\Delta\mathbf{x}=0` (using Lagrange multipliers) in each iteration of the
 estimation. For instance, when :math:`\mathbf{g}` is the weak direction of the covariance of an initial state from an
 unconstrained solution, this can be used to compute points along the line of variations :cite:p:`milani2005`.
 Any constraint previously set on the parameter is replaced; it can be removed with
 :func:`~tudatpy.dynamics.parameters.remove_parameter_correction_constraint`.
 Estimations and covariance analyses with a constraint raise an error if the parameter set includes consider parameters.
 A constraint set on an initial rotational state is not applied (a warning is given), as this parameter only uses its
 quaternion norm constraint.
 With a constraint, the covariance in the estimation or covariance analysis output has one additional row and column per
 constraint; the parameter covariance is the top-left block with the size of the parameter vector. This block is the one used by
 :func:`~tudatpy.estimation.estimation_analysis.propagate_covariance_from_analysis_objects`.

 Parameters
 ----------
 parameter : :class:`~tudatpy.dynamics.parameters.EstimatableParameter`
     Estimated parameter, as retrieved from the parameter set with
     :meth:`~tudatpy.dynamics.parameters.EstimatableParameterSet.parameters_for_parameter_type`.
 constraint_direction : numpy.ndarray[numpy.float64[n, 1]]
     Finite, non-zero direction (of the same size as the parameter) to which the parameter correction is to be orthogonal.

     )doc" );

    m.def( "set_initial_state_correction_constraint",
           &tep::setInitialStateCorrectionConstraint< STATE_SCALAR_TYPE >,
           py::arg( "parameter_set" ),
           py::arg( "constraint_direction" ),
           R"doc(

 Function to constrain the correction to the initial state in a parameter set to be orthogonal to a given direction.

 Function that applies :func:`~tudatpy.dynamics.parameters.set_parameter_correction_constraint` to the single-arc initial
 translational state in the parameter set. The parameter set must contain exactly one such parameter.

 Parameters
 ----------
 parameter_set : :class:`~tudatpy.dynamics.parameters.EstimatableParameterSet`
     Parameter set containing exactly one single-arc initial translational state.
 constraint_direction : numpy.ndarray[numpy.float64[6, 1]]
     Direction to which the correction to the initial state is to be orthogonal.

     )doc" );

    m.def( "remove_parameter_correction_constraint",
           &tep::removeParameterCorrectionConstraint,
           py::arg( "parameter" ),
           R"doc(

 Function to remove the constraint on the correction to an estimated parameter.

 Function to remove the constraint set with :func:`~tudatpy.dynamics.parameters.set_parameter_correction_constraint`,
 so that subsequent estimations and covariance analyses are unconstrained.

 Parameters
 ----------
 parameter : :class:`~tudatpy.dynamics.parameters.EstimatableParameter`
     Estimated parameter from which the constraint is removed.

     )doc" );

    m.def( "remove_parameter_correction_constraints",
           &tep::removeParameterCorrectionConstraints< STATE_SCALAR_TYPE >,
           py::arg( "parameter_set" ),
           R"doc(

 Function to remove the constraints on the corrections to all parameters in a parameter set.

 Function that applies :func:`~tudatpy.dynamics.parameters.remove_parameter_correction_constraint` to all estimated
 parameters in the parameter set.

 Parameters
 ----------
 parameter_set : :class:`~tudatpy.dynamics.parameters.EstimatableParameterSet`
     Parameter set from which all constraints are removed.

     )doc" );
}

}  // namespace parameters
}  // namespace dynamics
}  // namespace tudatpy
