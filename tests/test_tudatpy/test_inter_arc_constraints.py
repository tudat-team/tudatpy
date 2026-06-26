import numpy as np
import pytest

from tudatpy.estimation import estimation_analysis


def test_inter_arc_constraint_factories_accept_expected_inputs():
    bodies = ["Earth", "Mars"]
    epochs = {"Earth": [100.0, 200.0], "Mars": [150.0]}

    full = estimation_analysis.full_state_continuity(
        bodies,
        epochs,
        position_weights={"Earth": [1.0, 2.0, 3.0], "Mars": 4.0},
        velocity_weights={"Earth": [0.1, 0.2, 0.3], "Mars": 0.4},
        constraint_scaling_factor=2.0,
    )
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

    dense_rank_deficient = np.ones((6, 1)) @ np.ones((1, 6))
    general = estimation_analysis.general_continuity(
        bodies,
        epochs,
        {"Earth": [dense_rank_deficient], "Mars": np.eye(6)},
        constraint_scaling_factor=3.0,
    )

    assert type(full).__name__ == "InterArcStateContinuityConstraintSettings"
    assert type(position).__name__ == "InterArcStateContinuityConstraintSettings"
    assert type(velocity).__name__ == "InterArcStateContinuityConstraintSettings"
    assert type(general).__name__ == "InterArcStateContinuityConstraintSettings"
    assert general.bodies == bodies
    assert general.connection_epochs == epochs
    assert general.constraint_scaling_factor == 3.0


def test_inter_arc_constraint_factories_validate_inputs():
    with pytest.raises(RuntimeError):
        estimation_analysis.position_only_continuity(
            ["Vehicle"],
            {"Vehicle": [100.0]},
            position_weights={"Vehicle": 1.0},
            constraint_scaling_factor=0.0,
        )

    with pytest.raises(RuntimeError):
        estimation_analysis.full_state_continuity(
            ["Vehicle"],
            {"Vehicle": [100.0]},
            position_weights={"Vehicle": [1.0, 2.0]},
            velocity_weights={"Vehicle": 1.0},
        )

    with pytest.raises(RuntimeError):
        estimation_analysis.position_only_continuity(
            ["Earth", "Mars"],
            {"Earth": [100.0], "Mars": [200.0]},
            position_weights={"Earth": 1.0},
        )

    with pytest.raises(RuntimeError):
        estimation_analysis.position_only_continuity(
            ["Earth"],
            {"Earth": [100.0]},
            position_weights={"Earth": 1.0, "Mars": 2.0},
        )

    indefinite = np.diag([-1.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    with pytest.raises(RuntimeError):
        estimation_analysis.general_continuity(
            ["Vehicle"], {"Vehicle": [100.0]}, {"Vehicle": [indefinite]}
        )

    with pytest.raises(RuntimeError):
        estimation_analysis.general_continuity(
            ["Vehicle"],
            {"Vehicle": [100.0, 200.0]},
            {"Vehicle": [np.eye(6), np.eye(6), np.eye(6)]},
        )


def test_inter_arc_constraints_attach_to_input():
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

    covariance_input = estimation_analysis.CovarianceAnalysisInput(None)
    covariance_input.set_inter_arc_continuity_constraints(constraints)
    covariance_readback = covariance_input.inter_arc_continuity_constraints
    assert len(covariance_readback) == 4
    assert [entry.constraint_scaling_factor for entry in covariance_readback] == [
        2.0,
        3.0,
        4.0,
        5.0,
    ]

    estimation_input = estimation_analysis.EstimationInput(None)
    estimation_input.set_inter_arc_continuity_constraints(constraints)

    readback = estimation_input.inter_arc_continuity_constraints
    assert len(readback) == 4
    assert [entry.constraint_scaling_factor for entry in readback] == [2.0, 3.0, 4.0, 5.0]
