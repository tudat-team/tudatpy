import numpy as np
import pytest

from tudatpy.estimation import estimation_analysis


def test_inter_arc_constraint_factories_accept_expected_inputs():
    epochs = [100.0, 200.0]

    full = estimation_analysis.full_state_continuity(
        "Vehicle",
        epochs,
        position_weight=[1.0, 2.0, 3.0],
        velocity_weight=[0.1, 0.2, 0.3],
        constraint_scaling_factor=2.0,
    )
    position = estimation_analysis.position_only_continuity(
        "Vehicle", epochs, position_weight=[1.0, 2.0, 3.0]
    )
    velocity = estimation_analysis.velocity_only_continuity(
        "Vehicle", epochs, velocity_weight=[0.1, 0.2, 0.3]
    )

    dense_rank_deficient = np.ones((6, 1)) @ np.ones((1, 6))
    general = estimation_analysis.general_continuity(
        "Vehicle", epochs, [dense_rank_deficient], constraint_scaling_factor=[3.0, 4.0]
    )

    assert type(full).__name__ == "InterArcStateContinuityConstraintSettings"
    assert type(position).__name__ == "InterArcStateContinuityConstraintSettings"
    assert type(velocity).__name__ == "InterArcStateContinuityConstraintSettings"
    assert type(general).__name__ == "InterArcStateContinuityConstraintSettings"
    assert general.body == "Vehicle"
    assert general.connection_epochs == epochs
    assert general.constraint_scaling_factors == [3.0, 4.0]


def test_inter_arc_constraint_factories_validate_inputs():
    with pytest.raises(RuntimeError):
        estimation_analysis.position_only_continuity(
            "Vehicle", [100.0], constraint_scaling_factor=0.0
        )

    with pytest.raises(RuntimeError):
        estimation_analysis.full_state_continuity("Vehicle", [100.0], position_weight=[1.0, 2.0])

    indefinite = np.diag([-1.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    with pytest.raises(RuntimeError):
        estimation_analysis.general_continuity("Vehicle", [100.0], [indefinite])

    with pytest.raises(RuntimeError):
        estimation_analysis.general_continuity(
            "Vehicle", [100.0, 200.0], [np.eye(6)], constraint_scaling_factor=[1.0, 2.0, 3.0]
        )


def test_inter_arc_constraints_attach_to_input():
    constraints = [
        estimation_analysis.full_state_continuity(
            "Vehicle", [100.0], constraint_scaling_factor=2.0
        ),
        estimation_analysis.position_only_continuity(
            "Vehicle", [200.0], constraint_scaling_factor=[3.0]
        ),
        estimation_analysis.velocity_only_continuity(
            "Vehicle", [300.0], constraint_scaling_factor=4.0
        ),
        estimation_analysis.general_continuity(
            "Vehicle", [400.0], [np.eye(6)], constraint_scaling_factor=[5.0]
        ),
    ]

    covariance_input = estimation_analysis.CovarianceAnalysisInput(None)
    covariance_input.set_inter_arc_continuity_constraints(constraints)
    covariance_readback = covariance_input.inter_arc_continuity_constraints
    assert len(covariance_readback) == 4
    assert [entry.mu_values for entry in covariance_readback] == [[2.0], [3.0], [4.0], [5.0]]

    estimation_input = estimation_analysis.EstimationInput(None)
    estimation_input.set_inter_arc_continuity_constraints(constraints)

    readback = estimation_input.inter_arc_continuity_constraints
    assert len(readback) == 4
    assert [entry.constraint_scaling_factors for entry in readback] == [[2.0], [3.0], [4.0], [5.0]]
