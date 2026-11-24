"""Check that a sparse design matrix is returned to Python as a scipy.sparse matrix equal to the dense one."""

import numpy as np
import scipy.sparse

from tudatpy.dynamics import environment_setup, parameters_setup, propagation_setup
from tudatpy.estimation import estimation_analysis, observations, observations_setup
from tudatpy.estimation.observable_models_setup import links, model_settings
from tudatpy.interface import spice


def test_sparse_design_matrix_in_covariance_output():
    spice.load_standard_kernels()

    # Satellite orbiting a point-mass Earth for one hour
    body_settings = environment_setup.get_default_body_settings(["Earth"], "Earth", "J2000")
    body_settings.add_empty_settings("Delfi")
    bodies = environment_setup.create_system_of_bodies(body_settings)

    acceleration_models = propagation_setup.create_acceleration_models(
        bodies,
        {"Delfi": {"Earth": [propagation_setup.acceleration.point_mass_gravity()]}},
        ["Delfi"],
        ["Earth"],
    )
    propagator_settings = propagation_setup.propagator.translational(
        ["Earth"],
        acceleration_models,
        ["Delfi"],
        np.array([7.0e6, 0.0, 0.0, 0.0, 7.5e3, 0.0]),
        0.0,
        propagation_setup.integrator.runge_kutta_fixed_step(
            60.0, propagation_setup.integrator.CoefficientSets.rk_4
        ),
        propagation_setup.propagator.time_termination(3600.0),
    )

    # Position observations of the satellite every 10 minutes, used to estimate its initial state
    link_definition = links.LinkDefinition(
        {links.observed_body: links.body_origin_link_end_id("Delfi")}
    )
    parameters = parameters_setup.create_parameter_set(
        parameters_setup.initial_states(propagator_settings, bodies), bodies
    )
    estimator = estimation_analysis.Estimator(
        bodies,
        parameters,
        [model_settings.cartesian_position(link_definition)],
        propagator_settings,
    )
    simulated_observations = observations.simulate_observations(
        [
            observations_setup.observations_simulation_settings.tabulated_simulation_settings(
                model_settings.position_observable_type,
                link_definition,
                np.arange(600.0, 3001.0, 600.0),
                reference_link_end_type=links.observed_body,
            )
        ],
        estimator.observation_simulators,
        bodies,
    )

    # Covariance analysis with a dense and with a sparse design matrix
    design_matrices = []
    for use_sparse_design_matrix in [False, True]:
        covariance_input = estimation_analysis.CovarianceAnalysisInput(simulated_observations)
        covariance_input.define_covariance_settings(print_output_to_terminal=False)
        covariance_input.use_sparse_design_matrix = use_sparse_design_matrix
        assert covariance_input.use_sparse_design_matrix == use_sparse_design_matrix
        design_matrices.append(estimator.compute_covariance(covariance_input).design_matrix)

    dense_design_matrix, sparse_design_matrix = design_matrices
    assert isinstance(dense_design_matrix, np.ndarray)
    assert isinstance(sparse_design_matrix, scipy.sparse.csc_matrix)
    np.testing.assert_allclose(
        sparse_design_matrix.toarray(), dense_design_matrix, rtol=1.0e-12, atol=0.0
    )
