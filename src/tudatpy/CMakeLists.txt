# Update version
configure_file(
    "${TUDATPY_SOURCE_DIR}/_version.py.in"
    "${TUDATPY_SOURCE_DIR}/src/tudatpy/__init__.py"
    @ONLY
)

# Include libraries and add_extension macro
include(TudatpyLibraries)
include(AddExtension)

# Add extensions
##################################

# Interface - Spice
add_extension(tudatpy.interface.spice)

# IO
add_extension(tudatpy.data.io)

# Astro
add_extension(tudatpy.astro.element_conversion)
add_extension(tudatpy.astro.frame_conversion)
add_extension(tudatpy.astro.fundamentals)
add_extension(tudatpy.astro.gravitation)
add_extension(tudatpy.astro.polyhedron_utilities)
add_extension(tudatpy.astro.time_conversion)
add_extension(tudatpy.astro.two_body_dynamics)

# Constants
add_extension(tudatpy.constants)

# Math
add_extension(tudatpy.math.geometry)
add_extension(tudatpy.math.interpolators)
add_extension(tudatpy.math.numerical_integrators)
add_extension(tudatpy.math.root_finders)

# Numerical simulation
add_extension(tudatpy.numerical_simulation)
add_extension(tudatpy.numerical_simulation.propagation)
add_extension(tudatpy.numerical_simulation.environment)
add_extension(tudatpy.numerical_simulation.estimation)

# Numerical simulation - Environment setup
add_extension(tudatpy.numerical_simulation.environment_setup)
add_extension(tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients)
add_extension(tudatpy.numerical_simulation.environment_setup.atmosphere)
add_extension(tudatpy.numerical_simulation.environment_setup.ephemeris)
add_extension(tudatpy.numerical_simulation.environment_setup.gravity_field)
add_extension(tudatpy.numerical_simulation.environment_setup.gravity_field_variation)
add_extension(tudatpy.numerical_simulation.environment_setup.ground_station)
add_extension(tudatpy.numerical_simulation.environment_setup.radiation_pressure)
add_extension(tudatpy.numerical_simulation.environment_setup.rigid_body)
add_extension(tudatpy.numerical_simulation.environment_setup.rotation_model)
add_extension(tudatpy.numerical_simulation.environment_setup.shape)
add_extension(tudatpy.numerical_simulation.environment_setup.shape_deformation)
add_extension(tudatpy.numerical_simulation.environment_setup.vehicle_systems)

# Numerical simulation - Propagation setup
add_extension(tudatpy.numerical_simulation.propagation_setup)
add_extension(tudatpy.numerical_simulation.propagation_setup.acceleration)
add_extension(tudatpy.numerical_simulation.propagation_setup.dependent_variable)
add_extension(tudatpy.numerical_simulation.propagation_setup.integrator)
add_extension(tudatpy.numerical_simulation.propagation_setup.mass_rate)
add_extension(tudatpy.numerical_simulation.propagation_setup.propagator)
add_extension(tudatpy.numerical_simulation.propagation_setup.torque)
add_extension(tudatpy.numerical_simulation.propagation_setup.thrust)

# Numerical simulation - Estimation setup
add_extension(tudatpy.numerical_simulation.estimation_setup)
add_extension(tudatpy.numerical_simulation.estimation_setup.observation)
add_extension(tudatpy.numerical_simulation.estimation_setup.parameter)

# Trajectory design
add_extension(tudatpy.trajectory_design.shape_based_thrust)
add_extension(tudatpy.trajectory_design.transfer_trajectory)

# Installation
message(STATUS "Installing tudatpy")
install(
    DIRECTORY "${TUDATPY_SOURCE_DIR}/src/tudatpy"
    DESTINATION "${Python_SITELIB}"
)
