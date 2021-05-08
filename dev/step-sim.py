class DynamicsSimulator:
    pass


# # Create simulation object and propagate dynamics.
dynamics_simulator = propagation_setup.SingleArcDynamicsSimulator(
    bodies,
    integrator_settings,
    propagator_settings
)


def guidance_algorithm(r1, r2, v1, v2):
    # do some fancy magic here
    pass

terminated = False

while not terminated:
    # standard step for retrieving the three pillars of tudat simulations
    environment, propagation, estimation = dynamics_simulator.step()

    # check if state is terminal according to propagator
    terminated = propagation.propagator.state

    # retrieve algorithm dependent variables from environment
    r1, v1 = environment.bodies.asteroid.state_vectors
    r2, v2 = environment.bodies.spacecraft.state_vectors

    # retrieve current guidance according to environment state
    direction, magnitude = guidance_algorithm(r1, r2, v1, v2)

    # set newly updated thrust settings
    dynamics_simulator.propagation.bodies.spacecraft.set_thrust(
        direction=direction,
        magnitude=magnitude)
