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
#include "expose_parameters_setup.h"

#include "scalarTypes.h"
#include "tudat/simulation/estimation_setup/createEstimatableParameters.h"

namespace py = pybind11;
namespace tep = tudat::estimatable_parameters;
namespace tp = tudat::propagators;
namespace tss = tudat::simulation_setup;
namespace tba = tudat::basic_astrodynamics;

namespace tudatpy
{
namespace dynamics
{
namespace parameters_setup
{

void expose_parameters_setup( py::module& m )
{
    py::enum_< tep::EstimatebleParametersEnum >( m, "EstimatableParameterTypes", R"doc(

         Enumeration of model parameters that are available for estimation.
         In order to establish a parameter estimation settings for a parameter of a certain type, use the function dedicated to this parameter type.
         Note that not all of the listed types might be accessible via functions in the python interface yet.






      )doc" )
            .value( "arc_wise_initial_body_state_type", tep::EstimatebleParametersEnum::arc_wise_initial_body_state )
            .value( "initial_body_state_type", tep::EstimatebleParametersEnum::initial_body_state )
            .value( "initial_rotational_body_state_type", tep::EstimatebleParametersEnum::initial_rotational_body_state )
            .value( "gravitational_parameter_type", tep::EstimatebleParametersEnum::gravitational_parameter )
            .value( "constant_drag_coefficient_type", tep::EstimatebleParametersEnum::constant_drag_coefficient )
            .value( "radiation_pressure_coefficient_type", tep::EstimatebleParametersEnum::radiation_pressure_coefficient )
            .value( "arc_wise_radiation_pressure_coefficient_type",
                    tep::EstimatebleParametersEnum::arc_wise_radiation_pressure_coefficient )
            .value( "spherical_harmonics_cosine_coefficient_block_type",
                    tep::EstimatebleParametersEnum::spherical_harmonics_cosine_coefficient_block )
            .value( "spherical_harmonics_sine_coefficient_block_type",
                    tep::EstimatebleParametersEnum::spherical_harmonics_sine_coefficient_block )
            .value( "constant_rotation_rate_type", tep::EstimatebleParametersEnum::constant_rotation_rate )
            .value( "rotation_pole_position_type", tep::EstimatebleParametersEnum::rotation_pole_position )
            .value( "constant_additive_observation_bias_type", tep::EstimatebleParametersEnum::constant_additive_observation_bias )
            .value( "arcwise_constant_additive_observation_bias_type",
                    tep::EstimatebleParametersEnum::arcwise_constant_additive_observation_bias )
            .value( "constant_relative_observation_bias_type", tep::EstimatebleParametersEnum::constant_relative_observation_bias )
            .value( "arcwise_constant_relative_observation_bias_type",
                    tep::EstimatebleParametersEnum::arcwise_constant_relative_observation_bias )
            .value( "ppn_parameter_gamma_type", tep::EstimatebleParametersEnum::ppn_parameter_gamma )
            .value( "ppn_parameter_beta_type", tep::EstimatebleParametersEnum::ppn_parameter_beta )
            .value( "ground_station_position_type", tep::EstimatebleParametersEnum::ground_station_position )
            .value( "equivalence_principle_lpi_violation_parameter_type",
                    tep::EstimatebleParametersEnum::equivalence_principle_lpi_violation_parameter )
            .value( "empirical_acceleration_coefficients_type",
                    tep::EstimatebleParametersEnum::empirical_acceleration_coefficients )  // TO
                                                                                           // EXPOSE
            .value( "arc_wise_empirical_acceleration_coefficients_type",
                    tep::EstimatebleParametersEnum::arc_wise_empirical_acceleration_coefficients )  // TO EXPOSE
            .value( "full_degree_tidal_love_number_type",
                    tep::EstimatebleParametersEnum::full_degree_tidal_love_number )  // TO EXPOSE
            .value( "single_degree_variable_tidal_love_number_type",
                    tep::EstimatebleParametersEnum::single_degree_variable_tidal_love_number )
            .value( "direct_dissipation_tidal_time_lag_type", tep::EstimatebleParametersEnum::direct_dissipation_tidal_time_lag )
            .value( "mean_moment_of_inertia_type", tep::EstimatebleParametersEnum::mean_moment_of_inertia )
            .value( "arc_wise_constant_drag_coefficient_type", tep::EstimatebleParametersEnum::arc_wise_constant_drag_coefficient )
            .value( "periodic_spin_variation_type", tep::EstimatebleParametersEnum::periodic_spin_variation )
            .value( "polar_motion_amplitude_type", tep::EstimatebleParametersEnum::polar_motion_amplitude )
            .value( "core_factor_type", tep::EstimatebleParametersEnum::core_factor )
            .value( "free_core_nutation_rate_type", tep::EstimatebleParametersEnum::free_core_nutation_rate )
            .value( "desaturation_delta_v_values_type", tep::EstimatebleParametersEnum::desaturation_delta_v_values )
            .value( "constant_time_drift_observation_bias_type", tep::EstimatebleParametersEnum::constant_time_drift_observation_bias )
            .value( "arc_wise_time_drift_observation_bias_type", tep::EstimatebleParametersEnum::arc_wise_time_drift_observation_bias )
            .value( "global_polynomial_clock_corrections_type", tep::EstimatebleParametersEnum::global_polynomial_clock_corrections )
            .value( "arc_wise_polynomial_clock_corrections_type", tep::EstimatebleParametersEnum::arc_wise_polynomial_clock_corrections )
            .value( "inverse_tidal_quality_factor_type", tep::EstimatebleParametersEnum::inverse_tidal_quality_factor )
            .value( "radiation_pressure_target_perpendicular_direction_scaling_factor_type",
                    tep::EstimatebleParametersEnum::source_perpendicular_direction_radiation_pressure_scaling_factor )
            .value( "radiation_pressure_target_direction_scaling_factor_type",
                    tep::EstimatebleParametersEnum::source_direction_radiation_pressure_scaling_factor )
            .value( "drag_component_scaling_factor_type", tep::EstimatebleParametersEnum::drag_component_scaling_factor )
            .value( "side_component_scaling_factor_type", tep::EstimatebleParametersEnum::side_component_scaling_factor )
            .value( "lift_component_scaling_factor_type", tep::EstimatebleParametersEnum::lift_component_scaling_factor )
            .export_values( );

    py::enum_< tba::EmpiricalAccelerationComponents >(
            m,
            "EmpiricalAccelerationComponents",
            R"doc(Enumeration of the available empirical acceleration components that are available to estimate.
            
            These are used in the :func:`~tudatpy.dynamics.parameters_setup.empirical_accelerations` function to specify which components of the empirical acceleration are to be estimated.
            )doc" )
            .value( "radial_empirical_acceleration_component",
                    tba::EmpiricalAccelerationComponents::radial_empirical_acceleration_component )
            .value( "along_track_empirical_acceleration_component",
                    tba::EmpiricalAccelerationComponents::along_track_empirical_acceleration_component )
            .value( "across_track_empirical_acceleration_component",
                    tba::EmpiricalAccelerationComponents::across_track_empirical_acceleration_component )
            .export_values( );

    py::enum_< tba::EmpiricalAccelerationFunctionalShapes >(
            m,
            "EmpiricalAccelerationFunctionalShapes",
            R"doc(Enumeration of the available empirical acceleration shapes that are available per component
            
            These are used in the :func:`~tudatpy.dynamics.parameters_setup.empirical_accelerations` function to specify the signature of the estimated empirical acceleration component.
            .)doc" )
            .value( "constant_empirical", tba::EmpiricalAccelerationFunctionalShapes::constant_empirical )
            .value( "sine_empirical", tba::EmpiricalAccelerationFunctionalShapes::sine_empirical )
            .value( "cosine_empirical", tba::EmpiricalAccelerationFunctionalShapes::cosine_empirical )
            .export_values( );

    py::class_< tep::EstimatableParameterSettings, std::shared_ptr< tep::EstimatableParameterSettings > >( m,
                                                                                                           "EstimatableParameterSettings",
                                                                                                           R"doc(

         Base class to defining settings of parameter to be estimated.

         Functional (base) class for settings of model parameter to be estimated.
         Settings of simple parameters types are managed via this class, more complex parameter types are handled by specialised derivatives of this class.
         Instances of either base or derived class can be created via dedicated functions.





      )doc" )
            .def_readwrite( "custom_partial_settings", &tep::EstimatableParameterSettings::customPartialSettings_ )
            .def_readwrite( "parameter_identifier",
                            &tep::EstimatableParameterSettings::parameterType_,
                            R"doc(
                            
Type and associated body of the parameter.

The identifier contains the type of the parameter, defined by the :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterTypes` enumeration, the body and (if applicable) the reference point to which the parameter is associated.
The identifier is represented by a tuple of the form ``(parameter_type, (body_name, reference_point_name))``.

:type: tuple[ :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterTypes`, tuple[str, str] ]
                            
                            )doc" );

    // # EstimatableParameterSettings --> EstimatableParameterSet #
    m.def( "create_parameter_set",
           &tss::createParametersToEstimate< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "parameter_settings" ),
           py::arg( "bodies" ),
           py::arg( "propagator_settings" ) = nullptr,
           py::arg( "consider_parameters_names" ) = std::vector< std::shared_ptr< tep::EstimatableParameterSettings > >( ),
           R"doc(

 Function for creating a consolidated parameter from the given estimatable parameter settings.

 Function for creating a consolidated parameter from the given estimatable parameter settings.
 The function checks for consistency between the parameter settings and the models contained in the simulation setup (given by the `bodies` & and `propagator_settings` parameters).


 Parameters
 ----------
 parameter_settings : list( :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` )
     List of objects that define the settings for the parameters that are to be created. Each entry in this list is typically created by a call to a function in the :ref:`parameters_setup` module

 bodies : :class:`~tudatpy.dynamics.environment.SystemOfBodies`
     Object consolidating all bodies and environment models, including ground station models, that constitute the physical environment.

 propagator_settings : :class:`~tudatpy.dynamics.propagation_setup.propagator.PropagatorSettings`
     Object containing the consolidated propagation settings of the simulation.

 Returns
 -------
 :class:`~tudatpy.dynamics.parameters.EstimatableParameterSet`
     Instance of :class:`~tudatpy.dynamics.parameters.EstimatableParameterSet` class, consolidating all estimatable parameters and simulation models.

 Examples
 --------
 .. code-block:: python

     # Create bodies
     bodies = ...
     # Define parameters settings
     parameter_settings = ...
     # Create the parameters that will be estimated
     parameters_to_estimate = dynamics.parameters_setup.create_parameter_set(parameter_settings, bodies)

 This code snippet closely follows what is done in: `Full Estimation Example <https://github.com/tudat-team/tudatpy-examples/blob/master/estimation/full_estimation_example.ipynb>`_.




     )doc" );

    // ###############    Initial States
    // ################################

    m.def( "initial_states",
           &tss::getInitialStateParameterSettings< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "propagator_settings" ),
           py::arg( "bodies" ),
           py::arg( "arc_initial_times" ) = std::vector< double >( ),
           R"doc(

 Function for creating parameter settings for initial state parameters.

 Function for creating a parameter settings object for initial state parameters.
 The function uses the propagator settings to determine which type of initial state parameter (single/multi/hybrid-arc; translational/rotational/... dynamics) is to be estimated,
 e.g. if a single-arc translational state propagator is defined, the function will automatically create the parameters for the associated initial state parameter

 .. note::
 
    This function return lists of parameter settings objects.
    This means that the return of this function cannot simply be added to the parameter settings objects of single parameters in a list creation statement.
    Instead, list concatenation is recommended. Please see the following example:

 .. code-block:: python

    # define single estimatable parameters
    single_parameter_1 = ...
    single_parameter_2 = ...
    ...

    # bad: list creation statement --> will result in nested list, undesired!
    list_of_all_parameters = [dynamics.parameters_setup.initial_states(...), single_parameter_1, single_parameter_2, ...]

    # better: list concatenation --> will result in simple list, desired!
    list_of_all_parameters = dynamics.parameters_setup.initial_states(...) + [single_parameter_1, single_parameter_2, ...]


 Parameters
 ----------
 propagator_settings : :class:`~tudatpy.dynamics.propagation_setup.propagator.PropagatorSettings`
     Object containing the consolidated propagation settings of the simulation in the context of which the given model parameters are to be estimated.

 bodies : :class:`~tudatpy.dynamics.environment.SystemOfBodies`
     Object consolidating all bodies and environment models that constitute the physical environment.

 arc_initial_times : List[ float ] = []
     Initial times of arcs, only required if arc-wise propagation settings are passed via the `propagator_settings` object.

 Returns
 -------
 List[ :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` ]
     List of :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` objects, one per component of each initial state in the simulation.







     )doc" );

    // ###############    Vehicle Model Parameters
    // ################################

    m.def( "constant_drag_coefficient",
           &tep::constantDragCoefficient,
           py::arg( "body" ),
           R"doc(

 Function for creating parameter settings for constant drag coefficients.

 Function for creating parameter settings object for a constant drag coefficient parameter :math:`C_{D}` (see :func:`~tudatpy.dynamics.propagation_setup.acceleration.aerodynamic` ).
 Using the constant drag coefficient as an estimatable parameter requires:

 * A :func:`~tudatpy.dynamics.environment_setup.aerodynamic_coefficients.constant` aerodynamic interface (with :func:`~negative_aerodynamic_frame_coefficients` as input to ``force_coefficients_frame`` ) to be defined for the body specified by the ``body`` parameter
 * The body specified by the ``body`` parameter to undergo :func:`~tudatpy.dynamics.propagation_setup.acceleration.aerodynamic` acceleration

 The estimated parameter modifies the :func:`~tudatpy.dynamics.propagation_setup.acceleration.aerodynamic`

 Parameters
 ----------
 body : str
     Name of the body, with whose drag acceleration model the estimatable parameter is associated.

 Returns
 -------
 :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
     :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` object for the specified body's constant drag coefficient.







     )doc" );

    m.def( "arcwise_constant_drag_coefficient",
           &tep::arcwiseDragCoefficient,
           py::arg( "body" ),
           py::arg( "arc_initial_times" ),
           R"doc(

 Function for creating parameter settings for arc-wise constant drag coefficients.

 Function for creating parameter settings object for arc-wise constant drag coefficients :math:`C_{D}` 
 (arc-wise version of :func:`~tudatpy.dynamics.parameters_setup.constant_drag_coefficient`).
 Using the arc-wise constant drag coefficient as an estimatable parameter requires:

 * A :func:`~tudatpy.dynamics.environment_setup.aerodynamic_coefficients.constant` aerodynamic interface to be defined for the body specified by the ``body`` parameter
 * The body specified by the ``body`` parameter to undergo :func:`~tudatpy.dynamics.propagation_setup.acceleration.aerodynamic` acceleration

 When using this parameter, whenever :math:`C_{D}` is required at a time :math:`t`, the index :math:`i` in the ``arc_initial_times`` ordered list is
 found for which :math:`t_{i}\le t<t_{i+1}` (or, if :math:`t` is larger than the largest value in the list, :math:`i` is set to be last index of the list),
 and the parameter entry representing :math:`C_{D,i}` will be used.

 .. note:: This parameter may be estimated for a single-arc propagation, or a multi-arc propagation. In the latter case, the arcs selected for the estimation of the drag coefficient may, but need not, correspond to the arcs used for a multi-arc propagation.


 Parameters
 ----------
 body : str
     Name of the body, with whose drag acceleration model the estimatable parameter is associated.
 arc_initial_times : List[ float ]
     Ordered list of times at which the arcs over which the drag coefficient is to be estimated will start.

 Returns
 -------
 :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
     Instance of :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` derived :class:`~tudatpy.dynamics.parameters_setup.ArcWiseDragCoefficientEstimatableParameterSettings` class
     for arc-wise treatment of the specified body's constant drag coefficient.


     )doc" );

    m.def( "drag_component_scaling", &tep::dragComponentScaling, py::arg( "body" ),
           R"doc(

 Function for creating parameter settings for aerodynamic drag scaling factor

 Function for creating parameter settings object for a scaling factor :math:`K` (initialized to 1.0) for the aerodynamic force along the drag direction
 (effectively scaling the drag coefficient :math:`C_{D}` (see :func:`~tudatpy.dynamics.propagation_setup.acceleration.aerodynamic` )

 Using the arc-wise constant drag coefficient as an estimatable parameter requires:

 * The body specified by the ``body`` parameter to undergo :func:`~tudatpy.dynamics.propagation_setup.acceleration.aerodynamic` acceleration

 Note that, unlike the :func:`constant_drag_coefficient` parameter, this parameter does not modify the drag coefficient itself, but works
 regardless of the type of aerodynamic coefficients (in any frame, and with any dependencies). Using this parameter, the aerodynamic
 force along the drag directon is scaled (multiplied) by the factor :math:`K` during each function evaluation.

 Parameters
 ----------
 body : str
     Name of the body, with whose aerodynamic acceleration model the estimatable parameter is associated.

 Returns
 -------
 :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
     Instance of :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` class that define the settings. )doc" );

    m.def( "side_component_scaling", &tep::sideComponentScaling, py::arg( "body" ),
           R"doc(

 Function for creating parameter settings for aerodynamic side force scaling factor

 As :func:`~drag_component_scaling`, but scales the force along the :math:`C_{S}` direction rather than the :math:`C_{D}` direction

 Parameters
 ----------
 body : str
     Name of the body, with whose aerodynamic acceleration model the estimatable parameter is associated.

 Returns
 -------
 :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
     Instance of :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` class that define the settings. )doc" );

    m.def( "lift_component_scaling", &tep::liftComponentScaling, py::arg( "body" ),
           R"doc(

 Function for creating parameter settings for aerodynamic lift force scaling factor

 As :func:`~drag_component_scaling`, but scales the force along the :math:`C_{L}` direction rather than the :math:`C_{D}` direction

 Parameters
 ----------
 body : str
     Name of the body, with whose aerodynamic acceleration model the estimatable parameter is associated.

 Returns
 -------
 :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
     Instance of :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` class that define the settings. )doc" );

    m.def( "radiation_pressure_coefficient",
           &tep::radiationPressureCoefficient,
           py::arg( "body" ),
           R"doc(

 Function for creating parameter settings for radiation pressure coefficients.

 Function for creating parameter settings object for a radiation pressure coefficient  :math:`C_{r}`.
 Using the radiation pressure coefficient as an estimatable parameter requires:

 * A :func:`~tudatpy.dynamics.environment_setup.radiation_pressure.cannonball_radiation_target` radiation pressure target model to be defined for the body specified by the ``body`` parameter
 * The body specified by the ``body`` parameter to undergo :func:`~tudatpy.dynamics.propagation_setup.acceleration.radiation_pressure` acceleration (which, if the body has multiple target model defined, has the ``target_type`` input set to :attr:`~tudatpy.dynamics.environment_setup.radiation_pressure.RadiationPressureTargetModelType.cannonball_target`)


 Parameters
 ----------
 body : str
     Name of the body, with whose radiation pressure acceleration model the estimatable parameter is associated.

 Returns
 -------
 :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
     :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` object for the specified body's radiation pressure coefficient.


     )doc" );

    m.def( "arcwise_radiation_pressure_coefficient",
           &tep::arcwiseRadiationPressureCoefficient,
           py::arg( "body" ),
           py::arg( "arc_initial_times" ),
           R"doc(

 Function for creating parameter settings for arc-wise radiation pressure coefficients.

 Function for creating parameter settings object for arc-wise radiation pressure coefficients :math:`C_{r}` (arc-wise version of :func:`~tudatpy.dynamics.parameters_setup.radiation_pressure_coefficient`).
 Using the radiation pressure coefficient as an estimatable parameter requires:

 * A :func:`~tudatpy.dynamics.environment_setup.radiation_pressure.cannonball_radiation_target` radiation pressure target model to be defined for the body specified by the ``body`` parameter
 * The body specified by the ``body`` parameter to undergo :func:`~tudatpy.dynamics.propagation_setup.acceleration.radiation_pressure` acceleration (which, if the body has multiple target model defined, has the ``target_type`` input set to :func:`~tudatpy.dynamics.environment_setup.radiation_pressure.cannonball_target`)

 When using this parameter, whenever :math:`C_{r}` is required at a time :math:`t`, the index math:`i` in the ``arc_initial_times`` ordered list is
 found for which :math:`t_{i}\le t<t_{i+1}` (or, if :math:`t` is larger than the largest value in the list, :math:`i` is set to be last index of the list),
 and the parameter entry representing :math:`C_{r,i}` will be used.

 .. note:: This parameter may be estimated for a single-arc propagation, or a multi-arc propagation. In the latter case, the arcs selected for the estimation of the radiation pressure coefficient may, but need not, correspond to the arcs used for a multi-arc propagation.

 Parameters
 ----------
 body : str
     Name of the body, with whose radiation pressure acceleration model the estimatable parameter is associated.

 arc_initial_times : List[ float ]
     List of times at which the arcs over which the radiation pressure coefficient is to be estimated will start.

 Returns
 -------
 :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
     Instance of :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` derived :class:`~tudatpy.dynamics.parameters_setup.ArcWiseRadiationPressureCoefficientEstimatableParameterSettings` class
     for arc-wise treatment of the specified body's radiation pressure coefficient.
           
     )doc" );

    m.def( "radiation_pressure_target_direction_scaling",
           &tep::radiationPressureTargetDirectionScaling,
           py::arg( "target_body" ),
           py::arg( "source_body" ),
           R"doc(
 Function for creating parameter settings for a radiation pressure acceleration scaling factor in target direction.

 Function for creating parameter settings for scaling the radiation pressure acceleration component in the direction from the body
 undergoing the acceleration to the source model. When using this parameter, the radiation pressure :math:`\mathbf{a}` is decomposed into
 a component :math:`\mathbf{a}_{\parallel}` and :math:`\mathbf{a}_{\perp}, such that :math:`\mathbf{a}=\mathbf{a}_{\parallel}+\mathbf{a}_{\perp}`,
 where the parallel direction is computed as the component parallel with the vector from the center of mass of the source direction to the center of mass of the target direction.
 The radiation pressure model has parameters :math:`c_{\parallel}` and :math:`c_{\perp}` (nominally set to unity) that modify the acceleration as:
           
 .. math::
    \mathbf{a}=c_{\parallel}\mathbf{a}_{\parallel}+c_{\perp}\mathbf{a}_{\perp}

 The present function creates settings for a parameter defining :math:`c_{\parallel}` 

 Using this parameter requires:

 * The body specified by the ``target_body`` parameter to undergo :func:`~tudatpy.dynamics.propagation_setup.acceleration.radiation_pressure` acceleration exerted by ``source_body``
           
 Parameters
 ----------
 target_body : str
     Name of the body on which radiation pressure is exerted
 source_body : str
     Name of the body exerting the radiation pressure

 Returns
 -------
 :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
     Instance of :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` class dedining parallel radiation pressure scaling

     )doc" );

    m.def( "radiation_pressure_target_perpendicular_direction_scaling",
           &tep::radiationPressureTargetPerpendicularDirectionScaling,
           py::arg( "target_body" ),
           py::arg( "source_body" ),
           R"doc(
 Function for creating parameter settings for a radiation pressure acceleration scaling factor perpendicular to target direction.

 Function for creating parameter settings for scaling the radiation pressure acceleration component perpenedicular to the direction from the body
 undergoing the acceleration to the source model. The present function creates settings for a parameter defining :math:`c_{\perp}`,
 see :func:`~radiation_pressure_target_direction_scaling`

 Using this parameter requires:

 * The body specified by the ``target_body`` parameter to undergo :func:`~tudatpy.dynamics.propagation_setup.acceleration.radiation_pressure` acceleration exerted by ``source_body``

 Parameters
 ----------
 target_body : str
     Name of the body on which radiation pressure is exerted
 source_body : str
     Name of the body exerting the radiation pressure

 Returns
 -------
 :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
     Instance of :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` class dedining parallel radiation pressure scaling

     )doc" );

    m.def( "empirical_accelerations",
           &tep::empiricalAccelerationMagnitudes,
           py::arg( "body" ),
           py::arg( "centralBody" ),
           py::arg( "acceleration_components" ),
           R"doc(

 Function for creating parameter settings for empirical acceleration magnitudes.

 Function for creating parameter settings object for empirical acceleration magnitudes.
 Using the empirical acceleration terms as estimatable parameters requires:

 * The body specified by the ``body`` parameter to undergo :func:`~tudatpy.dynamics.propagation_setup.acceleration.empirical` acceleration, which include constant (in RSW frame) terms

 Any subset of the directions and functional shapes  can be estimated. The values in the parameter
 vector are ordered first by functional shape (constant, sine, cosine) and then by component (radial, normal, cross-track)
 For instance, if all nine coefficients are estimated, they will be ordered as :math:`\mathbf{a}_{R,\text{const.}},\mathbf{a}_{R,\text{sine}},\mathbf{a}_{R,\text{cosine}},\mathbf{a}_{S,\text{const.}},\mathbf{a}_{S,\text{sine}},\mathbf{a}_{S,\text{cosine}},\mathbf{a}_{W,\text{const.}},\mathbf{a}_{W,\text{sine}},\mathbf{a}_{W,\text{cosine}}`
 Any non-estimated components will be left to the values at which they were initialized.

 Parameters
 ----------
 body : str
     Name of the body, with whose empirical acceleration model the estimatable parameter is associated.

 centralBody : str
     Name of the central body of the empirical acceleration model (of which the gravitational parameter is extracted to compute the true anomaly, and w.r.t. which the RSW directions are determined). This body is the same as the body considered to be 'exerting' the empirical acceleration

 acceleration_components : dict[ EmpiricalAccelerationComponents, list[ EmpiricalAccelerationFunctionalShapes] ]
     Dictionary of components of the empirical acceleration which are to be estimated. There are two 'degrees of freedom' in these components: the direction of the acceleration (e.g. R, S or W direction) and the temporal signature (constant, sine of true anomaly or cosine of true anomaly). With this input, any subset may be selected. This parameter is a dictionary, with the key denoting the direction of the acceleration, and the value a list of the temporal signatures to estimate for this empirical acceleration direction.

 Returns
 -------
 :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
     Instance of :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` derived :class:`~tudatpy.dynamics.parameters_setup.EmpiricalAccelerationEstimatableParameterSettings` class
     for the specified body's empirical acceleration terms.

     )doc" );

    m.def( "constant_empirical_acceleration_terms",
           &tep::constantEmpiricalAccelerationMagnitudes,
           py::arg( "body" ),
           py::arg( "centralBody" ),
           R"doc(

 Function for creating parameter settings for constant empirical acceleration terms.

 As :func:`~tudatpy.dynamics.parameters_setup.empirical_accelerations`, but only using the constant R, S and W components (no sine or cosine term estimation). This function is added as a function of convenience


 Parameters
 ----------
 body : str
     Name of the body, with whose empirical acceleration model the estimatable parameter is associated.

 centralBody : str
     Name of the central body of the empirical acceleration model (of which the gravitational parameter is extracted to compute the true anomaly, and w.r.t. which the RSW directions are determined). This body is the same as the body considered to be 'exerting' the empirical acceleration

 Returns
 -------
 :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
     Instance of :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` derived :class:`~tudatpy.dynamics.parameters_setup.EmpiricalAccelerationEstimatableParameterSettings` class
     for the specified body's empirical acceleration terms.


     )doc" );

    m.def( "full_empirical_acceleration_terms",
           &tep::empiricalAccelerationMagnitudesFull,
           py::arg( "body" ),
           py::arg( "centralBody" ),
           R"doc(

 Function for creating parameter settings for empirical acceleration magnitudes for all components.

 As :func:`~tudatpy.dynamics.parameters_setup.empirical_accelerations`, but using selecting all nine components. This function is added as a function of convenience

 Parameters
 ----------
 body : str
     Name of the body, with whose empirical acceleration model the estimatable parameter is associated.

 centralBody : str
     Name of the central body of the empirical acceleration model (of which the gravitational parameter is extracted to compute the true anomaly, and w.r.t. which the RSW directions are determined). This body is the same as the body considered to be 'exerting' the empirical acceleration

 Returns
 -------
 :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
     Instance of :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` derived :class:`~tudatpy.dynamics.parameters_setup.EmpiricalAccelerationEstimatableParameterSettings` class
     for the specified body's empirical acceleration terms.

     )doc" );

    m.def( "arcwise_empirical_accelerations",
           &tep::arcWiseEmpiricalAccelerationMagnitudes,
           py::arg( "body" ),
           py::arg( "centralBody" ),
           py::arg( "acceleration_components" ),
           py::arg( "arc_start_times" ),
           R"doc(

 Function for creating parameter settings for arc-wise empirical acceleration magnitudes.

 Function for creating parameter settings object for arc-wise empirical acceleration magnitudes (arc-wise version of :func:`~tudatpy.dynamics.parameters_setup.empirical_accelerations`).
 Using the empirical acceleration terms as estimatable parameters requires:

 * The body specified by the ``body`` parameter to undergo :func:`~tudatpy.dynamics.propagation_setup.acceleration.empirical` acceleration, which include constant (in RSW frame) terms

 When using this parameter, whenever an empirical acceleration is required at a time :math:`t`, the index math:`i` in the ``arc_initial_times`` ordered list is
 found for which :math:`t_{i}\le t<t_{i+1}` (or, if :math:`t` is larger than the largest value in the list, :math:`i` is set to be last index of the list),
 and the parameter values representing empirical acceleration components in arc :math:`i` will be used.

 .. note:: This parameter may be estimated for a single-arc propagation, or a multi-arc propagation. In the latter case, the arcs selected for the estimation of the radiation pressure coefficient may, but need not, correspond to the arcs used for a multi-arc propagation.


 Parameters
 ----------
 body : str
     Name of the body, with whose empirical acceleration model the estimatable parameter is associated.
 centralBody : str
     Name of the central body of the empirical acceleration model (of which the gravitational parameter is extracted to compute the true anomaly, and w.r.t. which the RSW directions are determined). This body is the same as the body considered to be 'exerting' the empirical acceleration
 acceleration_components : Dict[ EmpiricalAccelerationComponents, List[ EmpiricalAccelerationFunctionalShapes] ]
     Dictionary of components of the empirical acceleration which are to be estimated. There are two 'degrees of freedom' in these components: the direction of the acceleration (e.g. R, S or W direction) and the temporal signature (constant, sine of true anomaly or cosine of true anomaly). With this input, any subset may be selected. This parameter is a dictionary, with the key denoting the direction of the acceleration, and the value a list of the temporal signatures to estimate for this empirical acceleration direction.
 arc_initial_times : List[ float ]
     List of times at which the arcs over which the empirical accelerations are to be estimated will start.

 Returns
 -------
 :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
     Instance of :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` derived :class:`~tudatpy.dynamics.parameters_setup.EmpiricalAccelerationEstimatableParameterSettings` class
     for the specified body's arc-wise empirical acceleration terms.



     )doc" );

    m.def( "arcwise_constant_empirical_acceleration_terms",
           &tep::constantArcWiseEmpiricalAccelerationMagnitudes,
           py::arg( "body" ),
           py::arg( "centralBody" ),
           py::arg( "arc_start_times" ),
           R"doc(

 Function for creating parameter settings for arc-wise constant empirical acceleration terms.

 As :func:`~tudatpy.dynamics.parameters_setup.arcwise_empirical_accelerations`, but only using the constant R, S and W components (no sine or cosine term estimation). This function is added as a function of convenience

 Parameters
 ----------
 body : str
     Name of the body, with whose empirical acceleration model the estimatable parameter is associated.

 centralBody : str
     Name of the central body of the empirical acceleration model (of which the gravitational parameter is extracted to compute the true anomaly, and w.r.t. which the RSW directions are determined). This body is the same as the body considered to be 'exerting' the empirical acceleration

 arc_initial_times : List[ float ]
     List of times at which the arcs over which the empirical accelerations are to be estimated will start.

 Returns
 -------
 :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
     Instance of :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` derived :class:`~tudatpy.dynamics.parameters_setup.EmpiricalAccelerationEstimatableParameterSettings` class
     for the specified body's arc-wise constant empirical acceleration terms.


     )doc" );

    m.def( "quasi_impulsive_shots",
           &tep::quasiImpulsiveShots,
           py::arg( "body" ),
           R"doc(

 Function for creating parameter settings for quasi-impulsive shots.

 Function for creating parameter settings object for so-called 'quasi-impulsive shots', such as desaturation maneuvers.
 With this parameter, the total :math:`\Delta \mathbf{V}` vector of a set of such maneuvers can be estimated (see :func:`~tudatpy.dynamics.propagation_setup.acceleration.quasi_impulsive_shots_acceleration` for mathematical details).
 Using the quasi-impulsive shots as an estimatable parameter requires:

 * The body specified by the ``body`` parameter to undergo :func:`~tudatpy.dynamics.propagation_setup.acceleration.quasi_impulsive_shots_acceleration` acceleration

 .. note:: this parameter considers *all* shots/maneuvers used in the above acceleration model, and estimates the value of the 'delta_v_values' input of this acceleration.

 Parameters
 ----------
 body : str
     Name of the body, with which the quasi-impulsive shot estimatable parameter is associated.

 Returns
 -------
 :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
     :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` object for the specified body's quasi-impulsive shots







     )doc" );

    // ###############    Gravity Model Parameters
    // ################################

    m.def( "gravitational_parameter",
           &tep::gravitationalParameter,
           py::arg( "body" ),
           R"doc(

 Function for creating parameter settings for a massive body's gravitational parameter.

 Function for creating parameter settings object for the gravitational parameter of massive bodies.
 Using the gravitational parameter as estimatable parameter requires:

 * The body specified by the ``body`` parameter to be endowed with a gravity field (see :ref:`gravity_field` module for options)
 * Any dynamical or observational model to depend on the gravitational parameter of the body specified by the ``body`` parameter


 Parameters
 ----------
 body : str
     Name of the body, with whose gravitational model the estimatable parameter is associated.

 Returns
 -------
 :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
     :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` object for the specified body's gravitational parameter.







     )doc" );

    m.def( "spherical_harmonics_c_coefficients",
           py::overload_cast< const std::string, const int, const int, const int, const int >( &tep::sphericalHarmonicsCosineBlock ),
           py::arg( "body" ),
           py::arg( "minimum_degree" ),
           py::arg( "minimum_order" ),
           py::arg( "maximum_degree" ),
           py::arg( "maximum_order" ),
           R"doc(

 Function for creating parameter settings for the cosine coefficients of body's spherical harmonics gravitational model.

 Function for creating parameter settings object for the spherical harmonics cosine-coefficients (:math:`\bar{C}_{lm}`) of a body with a spherical harmonic gravity field. Using this function, a 'full' set of spherical harmonic coefficients between an minimum/maximum degree/order are estimated. For instance, for minimum degree/order of 2/0, and maximum degree/order 4/4, all spherical harmonic cosine coefficients of degrees 2, 3 and 4 are estimated. If the maximum degree/order is set to 4/2, only coefficients with an order of 0, 1 and 2 are included. The entries in the parameter are sorted first by degree, and then by order (both in ascending order)
 Using the spherical harmonics cosine coefficients as estimatable parameter requires:

 * A :func:`~tudatpy.dynamics.environment_setup.gravity_field.spherical_harmonic` (or derived) gravity model to be defined for the body specified by the ``body`` parameter
 * Any dynamical or observational model to depend on the estimated cosine coefficients of the body specified by the ``body`` parameter. Typically, this dependency will be a :func:`~tudatpy.dynamics.propagation_setup.acceleration.spherical_harmonic` acceleration


 Parameters
 ----------
 body : str
     Name of the body, with whose gravitational model the estimatable parameters are associated.

 minimum_degree : int
     Minimum degree of c-coefficients to be included.
 minimum_order : int
     Minimum order of c-coefficients to be included.
 maximum_degree : int
     Maximum degree of c-coefficients to be included.
 maximum_order : int
     Maximum order of c-coefficients to be included.
 Returns
 -------
 :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
     Instance of :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` derived :class:`~tudatpy.dynamics.parameters_setup.SphericalHarmonicEstimatableParameterSettings` class
     for the applicable spherical harmonics c-coefficients of the specified body's gravitational model.







     )doc" );

    m.def( "spherical_harmonics_c_coefficients_block",
           py::overload_cast< const std::string, std::vector< std::pair< int, int > > >( &tep::sphericalHarmonicsCosineBlock ),
           py::arg( "body" ),
           py::arg( "block_indices" ),
           R"doc(

 Function for creating parameter settings for the cosine coefficients of body's spherical harmonics gravitational model.

 As :class:`~tudatpy.dynamics.parameters_setup.spherical_harmonics_c_coefficients`, but with a manually defined set of coefficients.


 Parameters
 ----------
 body : str
     Name of the body, with whose gravitational model the estimatable parameters are associated.

 block_indices : List[ Tuple[int, int] ]
     List of block indices. The length of this list can be arbitrary, as long as the pairs are unique.
     For each pair, the first value is the degree and the second the order of the coefficient to be included.

 Returns
 -------
 :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
     Instance of :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` derived :class:`~tudatpy.dynamics.parameters_setup.SphericalHarmonicEstimatableParameterSettings` class
     for the applicable spherical harmonics c-coefficients of the specified body's gravitational model.







     )doc" );

    m.def( "spherical_harmonics_s_coefficients",
           py::overload_cast< const std::string, const int, const int, const int, const int >( &tep::sphericalHarmonicsSineBlock ),
           py::arg( "body" ),
           py::arg( "minimum_degree" ),
           py::arg( "minimum_order" ),
           py::arg( "maximum_degree" ),
           py::arg( "maximum_order" ),
           R"doc(

 Function for creating parameter settings for the sine coefficients of body's spherical harmonics gravitational model.

 Function for creating parameter settings object for the spherical harmonics sine-coefficients (:math:`\bar{S}_{lm}`) of a body with a spherical harmonic gravity field. Using this function, a 'full' set of spherical harmonic coefficients between an minimum/maximum degree/order are estimated. For instance, for minimum degree/order of 2/1 (there is no order 0 sine coefficient), and maximum degree/order 4/4, all spherical harmonic sine coefficients of degrees 2, 3 and 4 are estimated. If the maximum degree/order is set to 4/2, only coefficients with an order of 1 and 2 are included. The entries in the parameter are sorted first by degree, and then by order (both in ascending order)
 Using the spherical harmonics cosine coefficients as estimatable parameter requires:

 * A :func:`~tudatpy.dynamics.environment_setup.gravity_field.spherical_harmonic` (or derived) gravity model to be defined for the body specified by the ``body`` parameter
 * Any dynamical or observational model to depend on the estimated cosine coefficients of the body specified by the ``body`` parameter. Typically, this dependency will be a :func:`~tudatpy.dynamics.propagation_setup.acceleration.spherical_harmonic` acceleration


 Parameters
 ----------
 body : str
     Name of the body, with whose gravitational model the estimatable parameters are associated.

 minimum_degree : int
     Minimum degree of s-coefficients to be included.
 minimum_order : int
     Minimum order of s-coefficients to be included.
 maximum_degree : int
     Maximum degree of s-coefficients to be included.
 maximum_order : int
     Maximum order of s-coefficients to be included.
 Returns
 -------
 :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
     Instance of :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` derived :class:`~tudatpy.dynamics.parameters_setup.SphericalHarmonicEstimatableParameterSettings` class
     for the applicable spherical harmonics s-coefficients of the specified body's gravitational model.







     )doc" );

    m.def( "spherical_harmonics_s_coefficients_block",
           py::overload_cast< const std::string, std::vector< std::pair< int, int > > >( &tep::sphericalHarmonicsSineBlock ),
           py::arg( "body" ),
           py::arg( "block_indices" ),
           R"doc(

 Function for creating parameter settings for the sine coefficients of body's spherical harmonics gravitational model.

 As :class:`~tudatpy.dynamics.parameters_setup.spherical_harmonics_s_coefficients`, but with a manually defined set of coefficients.


 Parameters
 ----------
 body : str
     Name of the body, with whose gravitational model the estimatable parameters are associated.

 block_indices : List[ Tuple[int, int] ]
     List of block indices. The length of this list can be arbitrary, as long as the pairs are unique.
     For each pair, the first value is the degree and the second the order of the coefficient to be included.

 Returns
 -------
 :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
     Instance of :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` derived :class:`~tudatpy.dynamics.parameters_setup.SphericalHarmonicEstimatableParameterSettings` class
     for the applicable spherical harmonics s-coefficients of the specified body's gravitational model.







     )doc" );

    // ###############    Rotation Model Parameters
    // ################################

    m.def( "mean_moment_of_inertia",
           &tep::meanMomentOfInertia,
           py::arg( "body" ),
           R"doc(

 Function for creating parameter settings for a body's mean moment of inertia.

 Function for creating parameter settings object for a body's mean moment of inertia.
 In most cases, the mean moment of inertia will not influence the dynamics/observation directly and sensitivity to this parameter will not be included. The dynamics/observation will be sensitive to this parameter if the rotational dynamics of a relevant body is estimated.
 Using the mean moment of inertia as estimatable parameter requires:

 * The estimation of an initial rotational state of the body specified by the ``body`` parameter


 Parameters
 ----------
 body : str
     Name of the body, with whose body model the estimatable parameter is associated.

 Returns
 -------
 :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
     :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` object for the specified body's mean moment of inertia.







     )doc" );

    m.def( "constant_rotation_rate",
           &tep::constantRotationRate,
           py::arg( "body" ),
           R"doc(

 Function for creating parameter settings for a body's constant rotation rate.

 Function for creating parameter settings object for a body's constant rotation rate parameter.
 Using the constant rotation rate as estimatable parameter requires:

 * A :func:`~tudatpy.dynamics.environment_setup.rotation_model.simple` or :func:`~tudatpy.dynamics.environment_setup.rotation_model.simple_from_spice` rotation model specified by the ``body`` parameter
 * Any dynamical or observational model to depend on the rotation model of the body specified by the ``body`` parameter


 Parameters
 ----------
 body : str
     Name of the body, with whose rotation model the estimatable parameter is associated.

 Returns
 -------
 :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
     :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` object for the specified body's constant spin rate.

     )doc" );

    m.def( "rotation_pole_position",
           &tep::rotationPolePosition,
           py::arg( "body" ),
           R"doc(

 Function for creating parameter settings for a body's rotation pole position.

 Function for creating parameter settings object for a body's rotation pole position, parameterized by the constant pole rotation angles (:math:`\alpha` and :math:`\delta`).
 Using the rotation pole position as estimatable parameter requires:

 * A :func:`~tudatpy.dynamics.environment_setup.rotation_model.simple` or :func:`~tudatpy.dynamics.environment_setup.rotation_model.simple_from_spice` rotation model specified by the ``body`` parameter
 * Any dynamical or observational model to depend on the rotation model of the body specified by the ``body`` parameter


 Parameters
 ----------
 body : str
     Name of the body, with whose rotation model the estimatable parameter is associated.

 Returns
 -------
 :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
     :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` object for the specified body's rotation pole position.







     )doc" );

    m.def( "core_factor",
           &tep::coreFactor,
           py::arg( "body" ),
           R"doc(

 Function for creating parameter settings for a body's core factor.

 Function for creating parameter settings object for a body's core factor.
 Using the core factor as estimatable parameter requires

 * A :func:`~tudatpy.dynamics.environment_setup.rotation_model.mars_high_accuracy` rotation model specified by the ``body`` parameter
 * Any dynamical or observational model to depend on the rotation model of the body specified by the ``body`` parameter


 Parameters
 ----------
 body : str
     Name of the body, with whose rotation model the estimatable parameter is associated.

 Returns
 -------
 :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
     :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` object for the specified body's core factor.







     )doc" );

    m.def( "free_core_nutation_rate",
           &tep::freeCoreNutationRate,
           py::arg( "body" ),
           R"doc(

 Function for creating parameter settings for a body's free core nutation rate.

 Function for creating parameter settings object for a body's free core nutation rate.
 Using the free core nutation rate as estimatable parameter requires

 * A :func:`~tudatpy.dynamics.environment_setup.rotation_model.mars_high_accuracy` rotation model specified by the ``body`` parameter
 * Any dynamical or observational model to depend on the rotation model of the body specified by the ``body`` parameter


 Parameters
 ----------
 body : str
     Name of the body, with whose rotation model the estimatable parameter is associated.

 Returns
 -------
 :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
     :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` object for the specified body's free core nutation rate.







     )doc" );

    m.def( "periodic_spin_variations",
           &tep::periodicSpinVariations,
           py::arg( "body" ),
           R"doc(

 Function for creating parameter settings for a body's periodic spin variations.

 Function for creating parameter settings object for a body's periodic spin variation parameters.
 Using the mean moment of inertia as estimatable parameter requires:

 * A :func:`~tudatpy.dynamics.environment_setup.rotation_model.mars_high_accuracy` rotation model specified by the ``body`` parameter
 * Any dynamical or observational model to depend on the rotation model of the body specified by the ``body`` parameter


 Parameters
 ----------
 body : str
     Name of the body, with whose rotation model the estimatable parameter is associated.

 Returns
 -------
 :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
     :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` object for the specified body's periodic spin variations.







     )doc" );

    m.def( "polar_motion_amplitudes",
           &tep::polarMotionAmplitudes,
           py::arg( "body" ),
           R"doc(

 Function for creating parameter settings for a body's polar motion amplitudes.

 Function for creating parameter settings object for a body's polar motion amplitudes.
 Using the polar motion amplitudes as estimatable parameter requires

 * A :func:`~tudatpy.dynamics.environment_setup.rotation_model.mars_high_accuracy` rotation model specified by the ``body`` parameter
 * Any dynamical or observational model to depend on the rotation model of the body specified by the ``body`` parameter


 Parameters
 ----------
 body : str
     Name of the body, with whose rotation model the estimatable parameter is associated.

 Returns
 -------
 :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
     :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` object for the specified body's polar motion amplitudes.
     )doc" );

        m.def( "iau_rotation_model_pole",
           &tep::iauRotationModelNominalPoleParameterSettings,
           py::arg( "body" ),
           R"doc(

 Function for creating parameter settings for a body's nominal pole position in an IAU rotation model

 Function for creating parameter settings for a body's nominal pole position in an IAU rotation model
 Using this requires:

 * A :func:`~tudatpy.dynamics.environment_setup.rotation_model.iau_rotation_model` rotation model specified by the ``body`` parameter
 * Any dynamical or observational model to depend on the rotation model of the body specified by the ``body`` parameter

 This parameter estimates the :math:`[\alpha_{0},\delta_{0}]` variables of the :func:`~tudatpy.dynamics.environment_setup.rotation_model.iau_rotation_model` rotation model

 Parameters
 ----------
 body : str
     Name of the body, with whose rotation model the estimatable parameter is associated.

 Returns
 -------
 :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`

     :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` object for the specified body's property

     )doc" );


        m.def( "iau_rotation_model_pole_rate",
           &tep::iauRotationModelPoleRateParameterSettings,
           py::arg( "body" ),
           R"doc(

 Function for creating parameter settings for a body's pole precession rate in an IAU rotation model

 Function for creating parameter settings for a body's pole precession rate in an IAU rotation model
 Using this requires:

 * A :func:`~tudatpy.dynamics.environment_setup.rotation_model.iau_rotation_model` rotation model specified by the ``body`` parameter
 * Any dynamical or observational model to depend on the rotation model of the body specified by the ``body`` parameter

 This parameter estimates the :math:`[\dot{\alpha},\dot{\delta}]` variables of the :func:`~tudatpy.dynamics.environment_setup.rotation_model.iau_rotation_model` rotation model

 Parameters
 ----------
 body : str
     Name of the body, with whose rotation model the estimatable parameter is associated.

 Returns
 -------
 :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`

     :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` object for the specified body's property

     )doc" );

        m.def( "iau_rotation_model_longitudinal_librations",
           &tep::iauRotationModelLongitudinalLibrationParameterSettings,
           py::arg( "body" ),
           py::arg( "libration_angular_frequencies" ),
           R"doc(

 Function for creating parameter settings for a body's longitudinal libration amplitudes in an IAU rotation model

 Function for creating parameter settings for a body's longitudinal libration amplitudes in an IAU rotation model
 Using this requires:

 * A :func:`~tudatpy.dynamics.environment_setup.rotation_model.iau_rotation_model` rotation model specified by the ``body`` parameter
 * Any dynamical or observational model to depend on the rotation model of the body specified by the ``body`` parameter

 This parameter estimates a list of :math:`W_{i}` variables of the :func:`~tudatpy.dynamics.environment_setup.rotation_model.iau_rotation_model` rotation model.
 The values of :math:`i` for which :math:`W_{i}` is estimated is defined by the ``libration_angular_frequencies`` input, which defines the
 corresponding :math:`\omega_{W_i}` values for which the librations are to be estimated. Note that the parameters are ordered as in
 the ``libration_angular_frequencies`` vector.

 Parameters
 ----------
 body : str
     Name of the body, with whose rotation model the estimatable parameter is associated.

 Returns
 -------
 :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`

     :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` object for the specified body's property

     )doc" );



    // ###############   Observation Model Parameters
    // ################################

    m.def( "absolute_observation_bias",
           &tep::observationBias,
           py::arg( "link_ends" ),
           py::arg( "observable_type" ),
           R"doc(

 Function for creating parameter settings for an absolute observation bias.

 Function for creating parameter settings object for an observation's absolute bias parameter.
 Using the absolute observation bias as estimatable parameter requires:

 * The observation model (corresponding to the `link_ends` and `observable_type`) to include an absolute bias (:func:`~tudatpy.estimation.observable_models_setup.biases.absolute_bias`)


 Parameters
 ----------
 link_ends : Dict[:class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType`, Tuple[str, str]
     Set of link ends that define the geometry of the biased observations.

 observable_type : ObservableType
     Observable type of the biased observations.
 Returns
 -------
 :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
     Instance of the :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` derived :class:`~tudatpy.dynamics.parameters_setup.ConstantObservationBiasEstimatableParameterSettings`
     for the specified observation's arc-wise absolute bias.







     )doc" );

    m.def( "arcwise_absolute_observation_bias",
           &tep::arcwiseObservationBias,
           py::arg( "link_ends" ),
           py::arg( "observable_type" ),
           py::arg( "arc_start_times" ),
           py::arg( "time_link_end" ),
           R"doc(

 Function for creating parameter settings for arc-wise absolute observation bias.

 Function for creating parameter settings object for the arc-wise treatment of an observation's absolute bias parameter.
 Using the arc-wise absolute observation bias as estimatable parameter requires

 * The observation model (corresponding to the `link_ends` and `observable_type`) to include an arc-wise absolute bias (:func:`~tudatpy.estimation.observable_models_setup.biases.arcwise_absolute_bias`)


 Parameters
 ----------
 link_ends : Dict[:class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType`, Tuple[str, str]
     Set of link ends that define the geometry of the biased observations.

 observable_type : ObservableType
     Observable type of the biased observations.
 arc_start_times : List[ float ]
     List of times at which the arcs over which the bias is to be estimated will start.
 time_link_end : LinkEndType
     The link end type (transmitter, receiver, etc.) at which the arc_start_times is evaluated.
 Returns
 -------
 :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
     Instance of the :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` derived :class:`~tudatpy.dynamics.parameters_setup.ArcWiseConstantObservationBiasEstimatableParameterSettings`
     for the specified observation's arc-wise absolute bias.







     )doc" );

    m.def( "global_polynomial_clock_corrections",
           &tep::globalPolynomialClockCorrections,
           py::arg( "associated_body" ),
           py::arg( "associated_station" ),
           py::arg( "correction_powers" ) );

    m.def( "arc_wise_polynomial_clock_corrections",
           &tep::multiArcPolynomialClockCorrections,
           py::arg( "associated_body" ),
           py::arg( "associated_station" ),
           py::arg( "correction_powers" ),
           py::arg( "arc_indices" ) );

    m.def( "relative_observation_bias",
           &tep::relativeObservationBias,
           py::arg( "link_ends" ),
           py::arg( "observable_type" ),
           R"doc(

 Function for creating parameter settings for an relative observation bias.

 Function for creating parameter settings object for an observation's relative bias parameter.
 Using the relative observation bias as estimatable parameter requires

 * The observation model (corresponding to the `link_ends` and `observable_type`) to include a relative bias (:func:`~tudatpy.estimation.observable_models_setup.biases.relative_bias`)

 Parameters
 ----------
 link_ends : Dict[:class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType`, Tuple[str, str]
     Set of link ends that define the geometry of the biased observations.

 observable_type : ObservableType
     Observable type of the biased observations.
 Returns
 -------
 :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
     Instance of the :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` derived :class:`~tudatpy.dynamics.parameters_setup.ConstantObservationBiasEstimatableParameterSettings`
     for the specified observation's arc-wise relative bias.







     )doc" );

    m.def( "arcwise_relative_observation_bias",
           &tep::arcwiseRelativeObservationBias,
           py::arg( "link_ends" ),
           py::arg( "observable_type" ),
           py::arg( "arc_start_times" ),
           py::arg( "time_link_end" ),
           R"doc(

 Function for creating parameter settings for arc-wise absolute observation bias.

 Function for creating parameter settings object for the arc-wise treatment of an observation's relative bias parameter.
 Using the arc-wise relative observation bias as estimatable parameter requires

 * The observation model (corresponding to the `link_ends` and `observable_type`) to include an arc-wise relative bias (:func:`~tudatpy.estimation.observable_models_setup.biases.arcwise_relative_bias`)

 .. note:: This parameter may be estimated for a single-arc propagation, or a multi-arc propagation. In the latter case, the arcs selected for the estimation of the bias may, but need not, correspond to the arcs used for a multi-arc propagation.


 Parameters
 ----------
 link_ends : Dict[:class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType`, Tuple[str, str]
     Set of link ends that define the geometry of the biased observations.

 observable_type : ObservableType
     Observable type of the biased observations.
 arc_start_times : List[ float ]
     List of times at which the arcs over which the bias is to be estimated will start.
 time_link_end : LinkEndType
     The link end type (transmitter, receiver, etc.) at which the arc_start_times is evaluated.
 Returns
 -------
 :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
     Instance of the :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` derived :class:`~tudatpy.dynamics.parameters_setup.ArcWiseConstantObservationBiasEstimatableParameterSettings`
     for the specified observation's arc-wise relative bias.

     )doc" );

    m.def( "time_drift_observation_bias",
           &tep::timeDriftObservationBias,
           py::arg( "link_ends" ),
           py::arg( "observable_type" ),
           py::arg( "ref_epoch" ),
           py::arg( "time_link_end" ) );

    m.def( "arcwise_time_drift_observation_bias",
           &tep::arcwiseTimeDriftObservationBias,
           py::arg( "link_ends" ),
           py::arg( "observable_type" ),
           py::arg( "arc_start_times" ),
           py::arg( "ref_epochs" ),
           py::arg( "time_link_end" ) );

    m.def( "constant_time_bias",
           &tep::timeObservationBias,
           py::arg( "link_ends" ),
           py::arg( "observable_type" ),
           py::arg( "reference_link_end" ) );

    m.def( "arcwise_time_bias",
           &tep::arcwiseTimeObservationBias,
           py::arg( "link_ends" ),
           py::arg( "observable_type" ),
           py::arg( "arc_start_times" ),
           py::arg( "reference_link_end" ) );

    m.def( "ground_station_position",
           &tep::groundStationPosition,
           py::arg( "body" ),
           py::arg( "ground_station_name" ),
           R"doc(

 Function for creating parameter settings for ground station position bias.

 Function for creating parameter settings object for a ground station's body-fixed Cartesian position.
 Using the ground station position bias as estimatable parameter requires:

 * At least one observation model to rely on the specified ground station


 Parameters
 ----------
 body : str
     Body name identifying the body, with which the ground station is associated.
 ground_station_name : str
     Name which identifies the position-biased ground station.
 Returns
 -------
 :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
     :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` object for the specified ground station's position bias.


     )doc" );

    m.def( "reference_point_position",
           &tep::referencePointPosition,
           py::arg( "body" ),
           py::arg( "reference_point_name" ),
           R"doc(No documentation found.)doc" );

    // ###############  Tidal Model Parameters
    // ################################

    m.def( "direct_tidal_dissipation_time_lag",
           py::overload_cast< const std::string&, const std::string& >( &tep::directTidalDissipationLagTime ),
           py::arg( "body" ),
           py::arg( "deforming_body" ),
           R"doc(No documentation found.)doc" );

    m.def( "direct_tidal_dissipation_time_lag",
           py::overload_cast< const std::string&, const std::vector< std::string >& >( &tep::directTidalDissipationLagTime ),
           py::arg( "body" ),
           py::arg( "deforming_body" ),
           R"doc(No documentation found.)doc" );

    m.def( "inverse_tidal_quality_factor",
           py::overload_cast< const std::string&, const std::string& >( &tep::inverseTidalQualityFactor ),
           py::arg( "body" ),
           py::arg( "deforming_body" ),
           R"doc(No documentation found.)doc" );

    m.def( "inverse_tidal_quality_factor",
           py::overload_cast< const std::string&, const std::vector< std::string >& >( &tep::inverseTidalQualityFactor ),
           py::arg( "body" ),
           py::arg( "deforming_body" ),
           R"doc(No documentation found.)doc" );

    m.def( "order_invariant_k_love_number",
           py::overload_cast< const std::string&, const int, const std::vector< std::string >&, const bool >(
                   &tep::orderInvariantKLoveNumber ),
           py::arg( "deformed_body" ),
           py::arg( "degree" ),
           py::arg( "deforming_bodies" ),
           py::arg( "use_complex_love_number" ) = 0,
           R"doc(

Function for creating parameter settings for a body's :math:`k_{l}` Love number.

Function for creating parameter settings for a body's :math:`k_{l}` Love number. When using this function,
we assume (for the case of degree 2 Love number) that :math:`k_{20}=k_{21}=k_{22}`. The estimation of
the Love number can be limited to a subset of the bodies that raise a tide on the body undergoing tidal deformation.

Using the :math:`k_{l}` Love number as estimatable parameter requires:

* A :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.solid_body_tide` gravity field variation model in the ``deformed_body`` (or one the more complex ones such as :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.solid_body_tide_degree_order_variable_k`). The parameter settings have to match the specifics of the variation model. For instance, if ``use_complex_love_number`` is set to true, the gravity field variation has to have been created using a complex Love number
* Any dynamical model to depend on the gravity field of the body specified by the ``deformed_body`` parameter


Parameters
----------
deformed_body : str
    Name of the body that is undergoing tidal deformation
degree : int
    Degree :math:`l` of the Love number :math:`k_{l}` that is to be estimated
deforming_bodies : list[str]
    List of bodies that raise a tide on ``deformed_body`` for which the single Love number defined by this setting is to be used. If the list is left empty, all tide-raising bodies will be used. By using this parameter, the value of :math:`k_{l}` will be indentical for the tides raised by each body in this list once parameter values are reset, even if they were different upon environment initialization
use_complex_love_number: bool
    Boolean defining whether the estimated Love number is real or imaginary

Returns
-------
:class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
    Object for the specified body's Love number :math:`k_{l}` for the tides raised by the specified bodies


    )doc" );

    m.def( "order_varying_k_love_number",
           py::overload_cast< const std::string&, const int, const std::vector< int >&, const std::vector< std::string >&, const bool >(
                   &tep::orderVaryingKLoveNumber ),
           py::arg( "deformed_body" ),
           py::arg( "degree" ),
           py::arg( "orders" ),
           py::arg( "deforming_bodies" ),
           py::arg( "use_complex_love_number" ),
           R"doc(

Function for creating parameter settings for a body's :math:`k_{lm}` Love numbers.

Function for creating parameter settings for a body's :math:`k_{lm}` Love numbers. When using this function,
we assume (for the case of degree 2 Love number) that :math:`k_{20}\neq k_{21}\neq k_{22}`. The estimation of
the Love numbers can be limited to a subset of the bodies that raise a tide on the body undergoing tidal deformation.

Using the :math:`k_{lm}` Love number as estimatable parameter requires:

* A :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.solid_body_tide_degree_order_variable_k` gravity field variation model in the ``deformed_body`` (or :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.solid_body_tide_degree_order_variable_complex_k` if ``use_complex_love_number`` is set to true)
* Any dynamical model to depend on the gravity field of the body specified by the ``deformed_body`` parameter


Parameters
----------
deformed_body : str
    Name of the body that is undergoing tidal deformation
degree : int
    Degree :math:`l` of the Love numbers :math:`k_{lm}` that are to be estimated
degree : list[int]
    Orders :math:`m` of the Love numbers :math:`k_{lm}` that are to be estimated
deforming_bodies : list[str]
    List of bodies that raise a tide on ``deformed_body`` for which the Love numbers defined by this setting is to be used. If the list is left empty, all tide-raising bodies will be used. By using this parameter, the values of :math:`k_{lm}` will be indentical for the tides raised by each body in this list once parameter values are reset, even if they were different upon environment initialization
use_complex_love_number: bool
    Boolean defining whether the estimated Love number is real or imaginary

Returns
-------
:class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
    Object for the specified body's Love numbers :math:`k_{lm}` for the tides raised by the specified bodies


    )doc" );

    m.def( "mode_coupled_k_love_numbers",
           &tep::modeCoupledTidalLoveNumberEstimatableParameterSettings,
           py::arg( "deformed_body" ),
           py::arg( "love_number_indices" ),
           py::arg( "deforming_bodies" ),
           R"doc(

Function for creating parameter settings for a body's mode-coupled :math:`k_{lm}^{l'm'}` Love numbers.

Function for creating parameter settings for a body's :math:`k_{lm}^{l'm'}` Love numbers (see :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.mode_coupled_solid_body_tide model` details). The estimation of
the Love numbers can be limited to a subset of the bodies that raise a (mode-coupled) tide on the body undergoing tidal deformation.

Using the :math:`k_{lm}` Love number as estimatable parameter requires:

* A :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.mode_coupled_solid_body_tide` gravity field variation model in the ``deformed_body``.
* Any dynamical model to depend on the gravity field of the body specified by the ``deformed_body`` parameter


Parameters
----------
deformed_body : str
    Name of the body that is undergoing tidal deformation
love_number_per_degree : dict[tuple[int, int], list[int,int]]]
    Dictionary of Love number indices for each combination for forcing and response degree and order.
    The first tuple (key) is the forcing degree and order :math:`l,m`, the list of tuples (key) is the list of associated response degrees and orders :math:`l',m'`
    for which the Love numbers are to be estimated (see :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.mode_coupled_solid_body_tide` for mathematical definition))
deforming_bodies : list[str]
    List of bodies that raise a tide on ``deformed_body`` for which the Love numbers defined by this setting is to be used. If the list is left empty, all tide-raising bodies will be used. By using this parameter, the values of :math:`k_{lm}` will be indentical for the tides raised by each body in this list once parameter values are reset, even if they were different upon environment initialization

Returns
-------
:class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
    Object for the specified body's mode-coupled Love numbers :math:`k_{lm}^{l'm'}` for the tides raised by the specified bodies

    )doc" );

    m.def( "polynomial_gravity_field_variation_amplitudes",
           &tep::polynomialGravityFieldVariationParameter,
           py::arg( "body_name" ),
           py::arg( "cosine_indices_per_power" ),
           py::arg( "sine_indices_per_power" ),
           R"doc(

Function for creating parameter settings for a body's polynomial gravity field amplitudes.

Function for creating parameter settings for a body's polynomial gravity field amplitudes :math:`K_{i,\bar{C}_{lm}}` and :math:`K_{i,\bar{S}_{lm}}`,
as defined in :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.polynomial`.

Using this settings as estimatable parameter requires:

* A :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.polynomial` (or :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.single_power_polynomial`) gravity field variation model in the ``body_name``.
* Any dynamical model to depend on the gravity field of the body specified by the ``deformed_body`` parameter

When using this parameter, a subset of all the variation amplitudes defined in the gravity field variation model can be estimated. These are defined in the ``cosine_indices_per_power`` and ``sine_indices_per_power`` inputs

Parameters
----------
body_name : str
    Name of the body that is undergoing gravity field variation
cosine_indices_per_power : dict[int, list[int,int]]
    Dictionary of powers :math:`i` (as keys) with list of combinations of degrees :math:`l` and orders :math:`m` for which to estimate :math:`K_{i,\bar{C}_{lm}}` as values (see :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.polynomial` for mathematical definition)
sine_indices_per_power : dict[int, list[int,int]]
    Dictionary of powers :math:`i` (as keys) with list of combinations of degrees :math:`l` and orders :math:`m` for which to estimate :math:`K_{i,\bar{S}_{lm}}` as values (see :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.polynomial` for mathematical definition)

Returns
-------
:class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
    Object for the specified body's polynomial gravity field variations

    )doc" );

    m.def( "periodic_gravity_field_variation_amplitudes",
           &tep::periodicGravityFieldVariationParameter,
           py::arg( "body_name" ),
           py::arg( "cosine_indices_per_period" ),
           py::arg( "sine_indices_per_period" ),
           R"doc(

Function for creating parameter settings for a body's polynomial gravity field variation amplitudes

Function for creating parameter settings for a body's polynomial gravity field variation amplitudes
:math:`A_{i,\bar{C}_{lm}}`, :math:`B_{i,\bar{C}_{lm}}`, :math:`A_{i,\bar{S}_{lm}}` and :math:`B_{i,\bar{S}_{lm}}`
as defined in :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.single_period_periodic`.

Using this settings as estimatable parameter requires:

* A :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.periodic` (or :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.single_period_periodic`) gravity field variation model in the ``body_name``.
* Any dynamical model to depend on the gravity field of the body specified by the ``deformed_body`` parameter

When using this parameter, a subset of all the variation amplitudes defined in the gravity field variation model can be estimated.
These are defined in the ``cosine_indices_per_period`` and ``sine_indices_per_period`` inputs

Parameters
----------
body_name : str
    Name of the body that is undergoing gravity field variation
cosine_indices_per_period : dict[int, list[int,int]]
    Dictionary of frequency index :math:`i` (as keys; corresponding to frequency :math:`f_{i}`) with list of combinations of degrees :math:`l` and orders :math:`m` for which to estimate :math:`A_{i,\bar{C}_{lm}}` and :math:`B_{i,\bar{C}_{lm}}` as values (see :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.periodic` for mathematical definition)
sine_indices_per_period : dict[int, list[int,int]]
    Dictionary of frequency index :math:`i` (as keys; corresponding to frequency :math:`f_{i}`) with list of combinations of degrees :math:`l` and orders :math:`m` for which to estimate :math:`A_{i,\bar{S}_{lm}}` and :math:`B_{i,\bar{S}_{lm}}` as values (see :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.periodic` for mathematical definition)

Returns
-------
:class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
    Object for the specified body's periodic gravity field variations

    )doc" );

    m.def( "monomial_gravity_field_variation_amplitudes",
           &tep::polynomialSinglePowerGravityFieldVariationParameter,
           py::arg( "body_name" ),
           py::arg( "power" ),
           py::arg( "cosine_indices" ),
           py::arg( "sine_indices" ),
           R"doc(
Function for creating parameter settings for a body's polynomial gravity field amplitudes at a single power.

Identical to :func:`~polynomial`, but for only a single power.

Parameters
----------
body_name : str
    Name of the body that is undergoing gravity field variation
power: int
    Power :math:`i` for which to estimate polynomial gravity field variations
cosine_indices: list[int,int]
    List of combinations of degrees :math:`l` and orders :math:`m` for which to estimate :math:`K_{i,\bar{C}_{lm}}` (see :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.polynomial` for mathematical definition)
sine_indices: list[int,int]
    List of combinations of degrees :math:`l` and orders :math:`m` for which to estimate :math:`K_{i,\bar{S}_{lm}}` (see :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.polynomial` for mathematical definition)

Returns
-------
:class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
    Object for the specified body's polynomial gravity field variations

    )doc" );

    m.def( "monomial_full_block_gravity_field_variation_amplitudes",
           &tep::polynomialSinglePowerFullBlockGravityFieldVariationParameter,
           py::arg( "body_name" ),
           py::arg( "power" ),
           py::arg( "minimum_degree" ),
           py::arg( "minimum_order" ),
           py::arg( "maximum_degree" ),
           py::arg( "maximum_order" ),
           R"doc(

Function for creating parameter settings for a body's polynomial gravity field amplitudes at a single power.

Identical to :func:`~polynomial`, but for only a single power, and a full block of spherical harminic coefficient degrees :math:`l` and orders :math`m`
For each degree :math:`l_{\text{min}}\le l \le l_{\text{max}}`, variations are estimated for all orders :math:`m_{\text{min}}\le m \le \left( \text{min}(m_{\text{max}},l) \right)`

Parameters
----------
body_name : str
    Name of the body that is undergoing gravity field variation
power: int
    Power :math:`i` for which to estimate polynomial gravity field variations
minimum_degree: int
    Minimum degree :math:`l_{\text{min}}` for which :math:`K_{i,\bar{C}_{lm}}` and :math:`K_{i,\bar{S}_{lm}}` are to be estimated (see :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.polynomial` for mathematical definition)
minimum_order: int
    Minimum order :math:`m_{\text{min}}` for which :math:`K_{i,\bar{C}_{lm}}` and :math:`K_{i,\bar{S}_{lm}}` are to be estimated (see :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.polynomial` for mathematical definition)
maximum_degree: int
    Maximum degree :math:`l_{\text{max}}` for which :math:`K_{i,\bar{C}_{lm}}` and :math:`K_{i,\bar{S}_{lm}}` are to be estimated (see :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.polynomial` for mathematical definition)
maximum_order: int
    Maximum degree :math:`m_{\text{max}}` for which :math:`K_{i,\bar{C}_{lm}}` and :math:`K_{i,\bar{S}_{lm}}` are to be estimated (see :func:`~tudatpy.dynamics.environment_setup.gravity_field_variation.polynomial` for mathematical definition)

Returns
-------
:class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
    Object for the specified body's polynomial gravity field variations

    )doc" );

    m.def( "scaled_longitude_libration_amplitude",
           &tep::scaledLongitudeLibrationAmplitude,
           py::arg( "body_name" ),
           R"doc(No documentation found.)doc" );

    m.def( "yarkovsky_parameter",
           &tep::yarkovskyParameter,
           py::arg( "body_name" ),
           py::arg( "central_body_name" ) = "Sun",
           R"doc(

 Function for creating parameter settings for Yarkovsky parameter.

 Function for creating parameter settings for Yarkovsky acceleration parameter :math:`A_{2}` (see :func:`~tudatpy.dynamics.propagation_setup.acceleration.yarkovsky`).

 * The body specified by the ``body`` parameter to undergo :func:`~tudatpy.dynamics.propagation_setup.acceleration.yarkovsky` acceleration with ``central_body_name`` as body exerting the acceleration

 Parameters
 ----------
 body : str
     Name of the body, with whose radiation pressure acceleration model the estimatable parameter is associated.

 Returns
 -------
 :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
     :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` object for the specified body's radiation pressure coefficient.


     )doc" );

    m.def( "area_to_mass_ratio_scaling_parameter",
           &tep::areaToMassScaling,
           py::arg( "body_name" ),
           R"doc(

 Function for creating parameter settings for a scaling factor for a body's area to mass ratio(s)

 Using these parameter settings, a scaling factor for the all acceleration models that contain an area-to-mass scaling factor :math:`A/m`, specifically radiation pressure or aerodynamics.
 Upon initialization, the value of this parameter :math:`p` is equal to 1. When using it, it effectively scales the acceleration formulation such that
 :math:`A/m\rightarrow p(A/m)`. Estimating an area-to-mass ratio is typical in, for instance, orbit estimation of near-Earth space debris.

 However, since the mass of a body is not (necesarilly) a constant, and the reference area of a body is not (necesarilly) related to a physical surface area, so
 we have opted to implement an :math:`A/m` scaling factor as parameter instead. This scaling factor applies to both aerodynamics and radiation pressure, regardless of
 whether their reference areas are identical. It is also by definition comaptible with an object of varying mass.

 Using this settings as estimatable parameter requires:

 * The acceleration models of body ``body_name`` to include one or more accelerations that contains an area to mass ratio :math:`A/m` in its formulation, either :func:`~tudatpy.dynamics.propagation_setup.acceleration.aerodynamic` or :func:`~tudatpy.dynamics.propagation_setup.acceleration.radiation_pressure`. Each such acceleration will be scaled (multiplied) by the value of the parameter during the propagation

 Parameters
 ----------
 body_undergoing_acceleration : str
     Name of the body for which the area-to-mass scaling factor is applied

 Returns
 -------
 :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
     Instance of the :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` class for the specified scaling factor
)doc" );

    m.def( "full_acceleration_scaling_parameter",
           &tep::fullAccelerationScaling,
           py::arg( "body_name_undergoing_acceleration" ),
           py::arg( "body_name_exerting_acceleration" ),
           py::arg( "acceleration_type" ),
           R"doc(

 Function for creating parameter settings for a scaling factor for a single acceleration acting on a body

 Using these parameter settings, a scaling factor :math:`p` is applied to a single acceleration :math:`\mathbf{a}`, increasing or decreasing it according
 to the value of the parameter :math:`p` as :math:`\mathbf{a}\rightarrow p\mathbf{a}`. This parameter is typically used in an estimation
 to absorb a (constant scaling) mismodelling in a single acceleration. The value of :math:`p` is initialized to 1 upon parameter creation

 Using this settings as estimatable parameter requires:

 * The body ``body_undergoing_acceleration`` undergoing an acceleration of exerted by ``body_exerting_acceleration`` of type ``acceleration_type``

 Parameters
 ----------
 body_undergoing_acceleration : str
     Name of the body undergoing the acceleration
 body_exerting_acceleration : str
     Name of the body exerting the acceleration
 acceleration_type : AvailableAcceleration
    Type of exerted acceleration

 Returns
 -------
 :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
     Instance of the :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` class for the specified scaling factor
)doc" );

    m.def( "custom_parameter",
           &tep::customParameterSettings,
           py::arg( "custom_id" ),
           py::arg( "parameter_size" ),
           py::arg( "get_parameter_function" ),
           py::arg( "set_parameter_function" ),
           R"doc(No documentation found.)doc" );

    // ###############  Global (GR) Model Parameters
    // ################################

    m.def( "ppn_parameter_gamma",
           &tep::ppnParameterGamma,
           R"doc(

 Function for creating parameter settings for post-newtonian gamma parameter.

 Function for creating parameter settings object for a global PPN :math:`\gamma` parameter.
 Using the post-newtonian gamma parameter as estimatable parameter requires at least one of the following:

 * An acceleration model depending on this parameter, such as :func:`~tudatpy.dynamics.propagation_setup.acceleration.relativistic_correction`
 * An observation model with a light-time correction depending on this parameter, such as :func:`~tudatpy.estimation.observable_models_setup.light_time_corrections.first_order_relativistic_light_time_correction`

 Returns
 -------
 :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
     :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` object for a global post-newtonian :math:`\gamma` parameter.







     )doc" );

    m.def( "ppn_parameter_beta",
           &tep::ppnParameterBeta,
           R"doc(

 Function for creating parameter settings for post-newtonian beta parameter.

 Function for creating parameter settings object for a global PPN :math:`\beta` parameter.
 Using the post-newtonian gamma parameter as estimatable parameter requires at least one of the following:

 * An acceleration model depending on this parameter, such as :func:`~tudatpy.dynamics.propagation_setup.acceleration.relativistic_correction`
 * An observation model with a light-time correction depending on this parameter (none yet implemented)

 Returns
 -------
 :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings`
     :class:`~tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` object for a global post-newtonian :math:`\beta` parameter.







     )doc" );

    // CUSTOM AND ANALYTICAL ACCELERATION PARTIALS

    py::class_< tep::CustomAccelerationPartialSettings, std::shared_ptr< tep::CustomAccelerationPartialSettings > >(
            m, "CustomAccelerationPartialSettings", R"doc(No documentation found.)doc" );

    m.def( "custom_analytical_partial",
           &tep::analyticalAccelerationPartialSettings,
           py::arg( "analytical_partial_function" ),
           py::arg( "body_undergoing_acceleration" ),
           py::arg( "body_exerting_acceleration" ),
           py::arg( "acceleration_type" ),
           R"doc(No documentation found.)doc" );

    m.def( "custom_numerical_partial",
           &tep::numericalAccelerationPartialSettings,
           py::arg( "parameter_perturbation" ),
           py::arg( "body_undergoing_acceleration" ),
           py::arg( "body_exerting_acceleration" ),
           py::arg( "acceleration_type" ),
           py::arg( "environment_updates" ) = std::map< tp::EnvironmentModelsToUpdate, std::vector< std::string > >( ),
           R"doc(No documentation found.)doc" );
}

}  // namespace parameters_setup
}  // namespace dynamics
}  // namespace tudatpy
