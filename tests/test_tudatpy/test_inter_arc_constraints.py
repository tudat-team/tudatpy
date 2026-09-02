import numpy as np
import pytest

from tudatpy.estimation import estimation_analysis


def test_inter_arc_constraint_factories_accept_expected_inputs():
    """Verify all four Python factories accept multi-body scalar, vector, dense, and broadcast weights."""

    # Different epoch counts per body exercise body-keyed lookup and one-matrix-to-many-epochs broadcasting.
    bodies = ["Earth", "Mars"]
    epochs = {"Earth": [100.0, 200.0], "Mars": [150.0]}

    # Full-state construction mixes component-wise Earth weights with isotropic Mars weights.
    full = estimation_analysis.full_state_continuity(
        bodies,
        epochs,
        position_weights={"Earth": [1.0, 2.0, 3.0], "Mars": 4.0},
        velocity_weights={"Earth": [0.1, 0.2, 0.3], "Mars": 0.4},
        constraint_scaling_factor=2.0,
    )
    # The two masked factories independently cover vector and scalar conversion for their active state blocks.
    position = estimation_analysis.position_only_continuity(
        bodies,
        epochs,
        position_weights={"Earth": [1.0, 2.0, 3.0], "Mars": [4.0, 5.0, 6.0]},
    )
    velocity = estimation_analysis.velocity_only_continuity(
        bodies,
        epochs,
        velocity_weights={"Earth": [0.1, 0.2, 0.3], "Mars": 0.4},
    )

    # A dense rank-one matrix proves the general factory does not require diagonal or full-rank input; Mars' one
    # identity matrix simultaneously exercises broadcasting to its connection list.
    dense_rank_deficient = np.ones((6, 1)) @ np.ones((1, 6))
    general = estimation_analysis.general_continuity(
        bodies,
        epochs,
        {"Earth": [dense_rank_deficient], "Mars": np.eye(6)},
        constraint_scaling_factor=3.0,
    )

    # Every factory must return the common settings type expected by covariance and estimation inputs.
    assert type(full).__name__ == "InterArcStateContinuityConstraintSettings"
    assert type(position).__name__ == "InterArcStateContinuityConstraintSettings"
    assert type(velocity).__name__ == "InterArcStateContinuityConstraintSettings"
    assert type(general).__name__ == "InterArcStateContinuityConstraintSettings"
    # Read-only properties must preserve user ordering and values across the C++/Python boundary.
    assert general.bodies == bodies
    assert general.connection_epochs == epochs
    assert general.constraint_scaling_factor == 3.0


@pytest.mark.parametrize("scalar", [np.float32(1.0), np.float64(1.0), np.int32(1), np.int64(1)])
def test_inter_arc_constraint_factories_accept_numpy_weight_scalars(scalar):
    """Verify NumPy floating-point and integer scalar classes follow the isotropic-weight conversion path."""

    # Parameterization covers common NumPy scalar dtypes that are not exact Python float instances.
    settings = estimation_analysis.full_state_continuity(
        ["Vehicle"],
        {"Vehicle": [100.0]},
        position_weights={"Vehicle": scalar},
        velocity_weights={"Vehicle": scalar},
    )

    # Successful construction confirms the binding converted the scalar before C++ settings validation.
    assert type(settings).__name__ == "InterArcStateContinuityConstraintSettings"


def test_inter_arc_constraint_factories_accept_numpy_weight_arrays_as_sequences():
    """Verify a one-dimensional NumPy array converts to the three Cartesian component weights."""

    settings = estimation_analysis.position_only_continuity(
        ["Vehicle"],
        {"Vehicle": [100.0]},
        position_weights={"Vehicle": np.array([1.0, 2.0, 3.0])},
    )

    # The common settings result proves the array reached the component-wise overload successfully.
    assert type(settings).__name__ == "InterArcStateContinuityConstraintSettings"


def test_inter_arc_constraint_factories_validate_inputs():
    """Verify invalid scaling, map coverage, weight shape, definiteness, and broadcast length are rejected."""

    # A zero penalty denominator has no mathematical interpretation.
    with pytest.raises(RuntimeError):
        estimation_analysis.position_only_continuity(
            ["Vehicle"],
            {"Vehicle": [100.0]},
            position_weights={"Vehicle": 1.0},
            constraint_scaling_factor=0.0,
        )

    # Cartesian component weights must contain exactly three entries.
    with pytest.raises(RuntimeError):
        estimation_analysis.full_state_continuity(
            ["Vehicle"],
            {"Vehicle": [100.0]},
            position_weights={"Vehicle": [1.0, 2.0]},
            velocity_weights={"Vehicle": 1.0},
        )

    # Every requested body must have a corresponding weight; map access reports the missing Mars key.
    with pytest.raises(IndexError):
        estimation_analysis.position_only_continuity(
            ["Earth", "Mars"],
            {"Earth": [100.0], "Mars": [200.0]},
            position_weights={"Earth": 1.0},
        )

    # A negative eigenvalue would turn the soft prior into an indefinite objective contribution.
    indefinite = np.diag([-1.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    with pytest.raises(RuntimeError):
        estimation_analysis.general_continuity(
            ["Vehicle"], {"Vehicle": [100.0]}, {"Vehicle": [indefinite]}
        )

    # Three matrices cannot map unambiguously to two connection epochs and are not a broadcast input.
    with pytest.raises(RuntimeError):
        estimation_analysis.general_continuity(
            ["Vehicle"],
            {"Vehicle": [100.0, 200.0]},
            {"Vehicle": [np.eye(6), np.eye(6), np.eye(6)]},
        )

    # The current feature regularizes six-dimensional translational state discrepancies only.
    with pytest.raises(RuntimeError):
        estimation_analysis.general_continuity(
            ["Vehicle"], {"Vehicle": [100.0]}, {"Vehicle": [np.eye(3)]}
        )


def test_inter_arc_constraints_attach_to_input():
    """Verify all settings variants attach to and round-trip through covariance and estimation inputs."""

    # Distinct scaling factors make ordering and object preservation observable after readback.
    constraints = [
        estimation_analysis.full_state_continuity(
            ["Vehicle"],
            {"Vehicle": [100.0]},
            position_weights={"Vehicle": 1.0},
            velocity_weights={"Vehicle": 1.0},
            constraint_scaling_factor=2.0,
        ),
        estimation_analysis.position_only_continuity(
            ["Vehicle"],
            {"Vehicle": [200.0]},
            position_weights={"Vehicle": 1.0},
            constraint_scaling_factor=3.0,
        ),
        estimation_analysis.velocity_only_continuity(
            ["Vehicle"],
            {"Vehicle": [300.0]},
            velocity_weights={"Vehicle": 1.0},
            constraint_scaling_factor=4.0,
        ),
        estimation_analysis.general_continuity(
            ["Vehicle"],
            {"Vehicle": [400.0]},
            {"Vehicle": [np.eye(6)]},
            constraint_scaling_factor=5.0,
        ),
    ]

    # CovarianceAnalysisInput owns the base implementation of the constraint collection.
    covariance_input = estimation_analysis.CovarianceAnalysisInput(None)
    covariance_input.set_inter_arc_continuity_constraints(constraints)
    covariance_readback = covariance_input.inter_arc_continuity_constraints
    # Readback must retain all entries in insertion order and preserve each settings object's scaling.
    assert len(covariance_readback) == 4
    assert [entry.constraint_scaling_factor for entry in covariance_readback] == [
        2.0,
        3.0,
        4.0,
        5.0,
    ]

    # EstimationInput inherits the same API; verify the derived binding exposes identical behavior.
    estimation_input = estimation_analysis.EstimationInput(None)
    estimation_input.set_inter_arc_continuity_constraints(constraints)

    readback = estimation_input.inter_arc_continuity_constraints
    # The inherited collection must not drop or reorder any settings entry.
    assert len(readback) == 4
    assert [entry.constraint_scaling_factor for entry in readback] == [2.0, 3.0, 4.0, 5.0]
