# Update version
configure_file(
    "${TUDATPY_SOURCE_DIR}/../../__init__.py.in"
    "${TUDATPY_SOURCE_DIR}/__init__.py"
    @ONLY
)

# include(YacmaExtension)
include(PythonExtension)

# Setup the installation path.
set(TUDATPY_INSTALL_PATH "${YACMA_PYTHON_MODULES_INSTALL_PATH}/tudatpy")
# Add extensions
##################################

# add_extension(numerical_simulation)
# message(FATAL_ERROR "FINISHED VOLUNTARILY")

# Interface - Spice
add_extension(interface.spice)

# Data
add_extension(data)

# Astro
add_extension(astro.element_conversion)
add_extension(astro.ephemerides)
add_extension(astro.frame_conversion)
add_extension(astro.fundamentals)
add_extension(astro.gravitation)
add_extension(astro.polyhedron_utilities)
add_extension(astro.time_conversion)
add_extension(astro.two_body_dynamics)

# Constants
add_extension(constants)

# Math
add_extension(math.geometry)
add_extension(math.interpolators)
add_extension(math.numerical_integrators)
add_extension(math.root_finders)
add_extension(math.statistics)

# Numerical simulation
add_extension(numerical_simulation)
add_extension(numerical_simulation.propagation)
add_extension(numerical_simulation.environment)
add_extension(numerical_simulation.estimation)

# Numerical simulation - Environment setup
add_extension(numerical_simulation.environment_setup)
add_extension(numerical_simulation.environment_setup.aerodynamic_coefficients)
add_extension(numerical_simulation.environment_setup.atmosphere)
add_extension(numerical_simulation.environment_setup.ephemeris)
add_extension(numerical_simulation.environment_setup.gravity_field)
add_extension(numerical_simulation.environment_setup.gravity_field_variation)
add_extension(numerical_simulation.environment_setup.ground_station)
add_extension(numerical_simulation.environment_setup.radiation_pressure)
add_extension(numerical_simulation.environment_setup.rigid_body)
add_extension(numerical_simulation.environment_setup.rotation_model)
add_extension(numerical_simulation.environment_setup.shape)
add_extension(numerical_simulation.environment_setup.shape_deformation)
add_extension(numerical_simulation.environment_setup.vehicle_systems)

# Numerical simulation - Propagation setup
add_extension(numerical_simulation.propagation_setup)
add_extension(numerical_simulation.propagation_setup.acceleration)
add_extension(numerical_simulation.propagation_setup.dependent_variable)
add_extension(numerical_simulation.propagation_setup.integrator)
add_extension(numerical_simulation.propagation_setup.mass_rate)
add_extension(numerical_simulation.propagation_setup.propagator)
add_extension(numerical_simulation.propagation_setup.thrust)
add_extension(numerical_simulation.propagation_setup.torque)

# Numerical simulation - Estimation setup
add_extension(numerical_simulation.estimation_setup)
add_extension(numerical_simulation.estimation_setup.observation)
add_extension(numerical_simulation.estimation_setup.parameter)

# Trajectory design
# add_extension(trajectory_design.shape_based_thrust)
# add_extension(trajectory_design.transfer_trajectory)
add_extension(trajectory_design 1)

# Install
######################################################
install(
    DIRECTORY ${TUDATPY_SOURCE_DIR}
    DESTINATION ${YACMA_PYTHON_MODULES_INSTALL_PATH}
)
