import numpy as np
import pytest

from tudatpy.estimation import estimation_analysis


def _linear_problem(group_sizes, true_sigmas, number_of_parameters=6, seed=42):
    """Random linear problem y = A x + noise, with one noise level per group."""
    rng = np.random.default_rng(seed)
    design_matrix = rng.standard_normal((sum(group_sizes), number_of_parameters))
    observations = design_matrix @ np.linspace(1.0, 2.0, number_of_parameters)
    groups = []
    start = 0
    for size, sigma in zip(group_sizes, true_sigmas):
        rows = list(range(start, start + size))
        observations[rows] += sigma * rng.standard_normal(size)
        groups.append(rows)
        start += size
    return design_matrix, observations, groups


def _solve(design_matrix, weights, observations, inverse_apriori_covariance):
    normal_matrix = (
        design_matrix.T @ (weights[:, None] * design_matrix) + inverse_apriori_covariance
    )
    estimate = np.linalg.solve(normal_matrix, design_matrix.T @ (weights * observations))
    return normal_matrix, observations - design_matrix @ estimate


def _iterate(
    design_matrix, observations, groups, weights, inverse_apriori_covariance, use_simplified
):
    """Rescale the group weights until all variance factors are 1."""
    weights = weights.copy()
    for _ in range(200):
        normal_matrix, residuals = _solve(
            design_matrix, weights, observations, inverse_apriori_covariance
        )
        estimate = estimation_analysis.compute_variance_components_from_matrices(
            design_matrix, weights, normal_matrix, residuals, groups
        )
        factors = [
            result.simplified_variance_factor if use_simplified else result.variance_factor
            for result in estimate.group_results
        ]
        for rows, factor in zip(groups, factors):
            weights[rows] /= factor
        if max(abs(factor - 1.0) for factor in factors) < 1.0e-12:
            break
    return weights


def test_redundancy_matches_hat_matrix():
    """Group redundancies equal the sum of the diagonal of I - W A N^-1 A^T over the group rows."""
    design_matrix, observations, groups = _linear_problem([300, 100], [0.5, 20.0])
    weights = np.ones(len(observations))
    inverse_apriori_covariance = np.diag(np.linspace(1.0, 50.0, design_matrix.shape[1]))
    normal_matrix, residuals = _solve(
        design_matrix, weights, observations, inverse_apriori_covariance
    )

    estimate = estimation_analysis.compute_variance_components_from_matrices(
        design_matrix, weights, normal_matrix, residuals, groups
    )

    hat_diagonal = (
        np.einsum("ij,ji->i", design_matrix @ np.linalg.inv(normal_matrix), design_matrix.T)
        * weights
    )
    for rows, result in zip(groups, estimate.group_results):
        assert result.number_of_observations == len(rows)
        assert result.redundancy == pytest.approx(np.sum(1.0 - hat_diagonal[rows]), rel=1.0e-10)
        assert result.weighted_square_sum == pytest.approx(
            np.sum(weights[rows] * residuals[rows] ** 2), rel=1.0e-12
        )
    assert np.allclose(
        estimate.redundancies, [result.redundancy for result in estimate.group_results]
    )


def test_helmert_recovers_true_sigmas():
    """Starting 10x off, Helmert's method and the simplified estimator reach the same, correct weights."""
    true_sigmas = [0.5, 20.0]
    design_matrix, observations, groups = _linear_problem([400, 100], true_sigmas)
    initial_weights = np.empty(len(observations))
    for rows, sigma in zip(groups, [5.0, 2.0]):
        initial_weights[rows] = 1.0 / sigma**2
    inverse_apriori_covariance = 10.0 * np.eye(design_matrix.shape[1])

    helmert_weights = _iterate(
        design_matrix, observations, groups, initial_weights, inverse_apriori_covariance, False
    )
    simplified_weights = _iterate(
        design_matrix, observations, groups, initial_weights, inverse_apriori_covariance, True
    )

    assert np.allclose(helmert_weights, simplified_weights, rtol=1.0e-8)
    for rows, sigma in zip(groups, true_sigmas):
        assert 1.0 / np.sqrt(helmert_weights[rows[0]]) == pytest.approx(sigma, rel=0.15)


def test_block_weights_recover_correlated_noise():
    """Group 0 has noise correlated in pairs, with matching block weights. Its sigma is recovered."""
    scipy_sparse = pytest.importorskip("scipy.sparse")
    rng = np.random.default_rng(7)
    true_sigmas = [0.5, 20.0]
    design_matrix, observations, groups = _linear_problem([400, 100], [0.0, true_sigmas[1]])
    correlation = np.array([[1.0, 0.8], [0.8, 1.0]])
    cholesky = np.linalg.cholesky(correlation)
    for pair in range(200):
        observations[2 * pair : 2 * pair + 2] += true_sigmas[0] * cholesky @ rng.standard_normal(2)

    def weight_matrix(sigma0, sigma1):
        blocks = [np.linalg.inv(correlation) / sigma0**2] * 200 + [np.eye(100) / sigma1**2]
        return scipy_sparse.block_diag(
            [scipy_sparse.coo_array(block) for block in blocks], format="csc"
        )

    sigmas = np.array([5.0, 2.0])
    for _ in range(100):
        weights = weight_matrix(*sigmas)
        normal_matrix = design_matrix.T @ (weights @ design_matrix)
        estimate = np.linalg.solve(normal_matrix, design_matrix.T @ (weights @ observations))
        residuals = observations - design_matrix @ estimate
        result = estimation_analysis.compute_variance_components_from_matrices(
            design_matrix,
            weights.diagonal(),
            normal_matrix,
            residuals,
            groups,
            weights_matrix=weights,
        )
        sigmas = sigmas * np.sqrt(result.variance_factors)
        if np.max(np.abs(result.variance_factors - 1.0)) < 1.0e-12:
            break

    assert sigmas[0] == pytest.approx(true_sigmas[0], rel=0.15)
    assert sigmas[1] == pytest.approx(true_sigmas[1], rel=0.15)


def test_overlapping_groups_are_rejected():
    design_matrix, observations, groups = _linear_problem([50, 50], [1.0, 1.0])
    weights = np.ones(len(observations))
    normal_matrix, residuals = _solve(design_matrix, weights, observations, np.zeros((6, 6)))
    with pytest.raises(RuntimeError):
        estimation_analysis.compute_variance_components_from_matrices(
            design_matrix, weights, normal_matrix, residuals, [groups[0], groups[1] + [0]]
        )


def test_settings_factory_rejects_empty_groups():
    with pytest.raises(RuntimeError):
        estimation_analysis.variance_component_settings([])


def _orbit_estimation_setup():
    """Vehicle around a fixed Earth, observed with range from two stations and with its position."""
    from tudatpy.astro import element_conversion
    from tudatpy.dynamics import environment_setup, parameters_setup, propagation_setup
    from tudatpy.estimation import observable_models_setup, observations_setup
    from tudatpy.estimation import observations as observation_data

    gravitational_parameter = 3.986004418e14
    body_settings = environment_setup.BodyListSettings("SSB", "J2000")
    body_settings.add_empty_settings("Earth")
    earth = body_settings.get("Earth")
    earth.ephemeris_settings = environment_setup.ephemeris.constant(np.zeros(6), "SSB", "J2000")
    earth.gravity_field_settings = environment_setup.gravity_field.central(gravitational_parameter)
    earth.shape_settings = environment_setup.shape.spherical(6.371e6)
    earth.rotation_model_settings = environment_setup.rotation_model.constant_rotation_model(
        "J2000", "IAU_Earth", np.eye(3)
    )
    body_settings.add_empty_settings("Vehicle")
    bodies = environment_setup.create_system_of_bodies(body_settings)
    for i in range(2):
        environment_setup.add_ground_station(
            bodies.get("Earth"),
            environment_setup.ground_station.basic_station(
                f"Station{i}",
                [0.0, 0.3 * i, 1.0 + 2.0 * i],
                element_conversion.PositionElementTypes.geodetic_position_type,
            ),
        )

    acceleration_models = propagation_setup.create_acceleration_models(
        bodies,
        {"Vehicle": {"Earth": [propagation_setup.acceleration.point_mass_gravity()]}},
        ["Vehicle"],
        ["Earth"],
    )
    initial_state = element_conversion.keplerian_to_cartesian(
        np.array([7200.0e3, 0.05, 1.2, 0.5, 0.3, 0.0]), gravitational_parameter
    )
    end_time = 4.0 * 3600.0
    propagator_settings = propagation_setup.propagator.translational(
        ["Earth"],
        acceleration_models,
        ["Vehicle"],
        initial_state,
        0.0,
        propagation_setup.integrator.runge_kutta_fixed_step(
            10.0, propagation_setup.integrator.CoefficientSets.rk_4
        ),
        propagation_setup.propagator.time_termination(end_time),
    )
    parameter_set = parameters_setup.create_parameter_set(
        parameters_setup.initial_states(propagator_settings, bodies), bodies
    )

    links = observable_models_setup.links
    model_settings = observable_models_setup.model_settings
    range_links = []
    for i in range(2):
        range_links.append(
            links.LinkDefinition(
                {
                    links.transmitter: links.body_reference_point_link_end_id(
                        "Earth", f"Station{i}"
                    ),
                    links.receiver: links.body_origin_link_end_id("Vehicle"),
                }
            )
        )
    position_link = links.LinkDefinition(
        {links.observed_body: links.body_origin_link_end_id("Vehicle")}
    )
    observation_models = [model_settings.one_way_range(link) for link in range_links]
    observation_models.append(model_settings.cartesian_position(position_link))

    estimator = estimation_analysis.Estimator(
        bodies, parameter_set, observation_models, propagator_settings
    )

    # Seeded, independent Gaussian noise: 1 m on range, 50 m on each position component.
    rng = np.random.default_rng(1234)
    simulation_settings = observations_setup.observations_simulation_settings
    range_times = np.arange(60.0, end_time - 60.0, 30.0)
    position_times = np.arange(60.0, end_time - 60.0, 600.0)
    settings_list = [
        simulation_settings.tabulated_simulation_settings(
            model_settings.one_way_range_type,
            link,
            range_times,
            links.receiver,
            noise_function=lambda time: rng.normal(0.0, 1.0, 1),
        )
        for link in range_links
    ]
    settings_list.append(
        simulation_settings.tabulated_simulation_settings(
            model_settings.position_observable_type,
            position_link,
            position_times,
            links.observed_body,
            noise_function=lambda time: rng.normal(0.0, 50.0, 3),
        )
    )
    observations = observation_data.simulate_observation_dataset(
        settings_list, estimator.observation_simulators, bodies
    )
    return estimator, parameter_set, observations, model_settings


def test_orbit_estimation_with_variance_components():
    """Range and position sigmas, started 10x off, are recovered by the estimator method."""
    from tudatpy.estimation.observations import ObservationSelectionCondition as Condition

    estimator, parameter_set, observations, model_settings = _orbit_estimation_setup()

    groups = [
        Condition.observable_type(model_settings.one_way_range_type),
        Condition.observable_type(model_settings.position_observable_type),
    ]
    true_sigmas = [1.0, 50.0]
    initial_sigmas = [10.0, 5.0]
    for group, sigma in zip(groups, initial_sigmas):
        observations.set_constant_single_observation_scalar_weight(group, 1.0 / sigma**2)

    true_parameters = parameter_set.parameter_vector
    perturbed = true_parameters.copy()
    perturbed[:3] += 100.0
    perturbed[3:] += 0.1
    parameter_set.parameter_vector = perturbed

    estimation_input = estimation_analysis.EstimationInput(
        observations, convergence_checker=estimation_analysis.estimation_convergence_checker(6)
    )
    estimation_input.define_estimation_settings(print_output_to_terminal=False)
    settings = estimation_analysis.variance_component_settings(
        groups, tolerance=0.01, maximum_iterations=15
    )

    output = estimator.perform_estimation_with_variance_components(estimation_input, settings)

    assert output.is_converged
    assert len(output.iteration_results) == len(output.total_variance_factor_history)
    assert len(output.least_squares_best_iterations) == len(output.iteration_results)
    estimated_sigmas = np.array(initial_sigmas) * np.sqrt(output.total_variance_factors)
    assert estimated_sigmas[0] == pytest.approx(true_sigmas[0], rel=0.08)
    assert estimated_sigmas[1] == pytest.approx(true_sigmas[1], rel=0.3)

    # The final weights are left in the observation collection.
    range_weights = observations.get_weight_diagonal(groups[0])
    assert range_weights[0] == pytest.approx(1.0 / estimated_sigmas[0] ** 2, rel=1.0e-12)

    # The standalone function gives the same factors on the final output (all close to 1).
    estimate = estimation_analysis.compute_variance_components(
        output.estimation_output, observations, groups
    )
    assert np.allclose(estimate.variance_factors, output.iteration_results[-1].variance_factors)
    assert np.all(np.abs(estimate.variance_factors - 1.0) < 0.01)
    assert estimate.group_results[0].number_of_observations == 2 * len(
        np.arange(60.0, 4.0 * 3600.0 - 60.0, 30.0)
    )
