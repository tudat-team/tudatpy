"""Short native multi-arc example; run with a freshly built tudatpy on PYTHONPATH."""

import numpy as np

from tudatpy.dynamics import simulator
from tudatpy.data_input.environment_data import spice


def main():
    spice.load_standard_kernels()
    model = simulator.MarsGravityArcModel(degree=20)
    radius = model.reference_radius + 300_000.0
    speed = np.sqrt(model.gravitational_parameter / radius)
    inclination = np.deg2rad(93.0)
    states = np.array(
        [[radius, 0, 0, 0, speed * np.cos(inclination), speed * np.sin(inclination)]] * 2
    )
    states[1, 1] = 100.0
    simulation = simulator.IndependentArcVariationalSimulator(
        model=model,
        arc_start_times=[778_118_400.0, 778_120_200.0],
        initial_states=states,
        arc_duration=600.0,
        step_size=30.0,
    )
    simulation.set_independent_arc_execution(number_of_workers=2)
    simulation.integrate_full_equations()
    results = simulation.variational_propagation_results
    print("Native independent execution:", simulation.last_integration_used_independent_arcs)
    for index, arc in enumerate(simulation.propagation_results.single_arc_results):
        print(f"Arc {index}: {len(arc.state_history)} saved states")
    epoch = 778_118_700.0
    matrix = simulation.state_transition_interface.full_state_transition_sensitivity_at_epoch(epoch)
    print("Global state-transition/sensitivity matrix:", matrix.shape)
    print(
        "Shared Mars-mu sensitivity:",
        results.single_arc_results[0].sensitivity_matrix_history[epoch].ravel(),
    )


if __name__ == "__main__":
    main()
