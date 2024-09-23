# Update version
configure_file(
    "${TUDATPY_SOURCE_DIR}/../../_version.py.in"
    "${TUDATPY_SOURCE_DIR}/__init__.py"
    @ONLY
)

include(YacmaExtension)

# Setup the installation path.
set(TUDATPY_INSTALL_PATH "${YACMA_PYTHON_MODULES_INSTALL_PATH}/tudatpy")
# Add extensions
##################################

# Interface - Spice
yacma_extension(tudatpy.interface.spice)

# Data
yacma_extension(tudatpy.data)

# Astro
yacma_extension(tudatpy.astro.element_conversion)
yacma_extension(tudatpy.astro.frame_conversion)
yacma_extension(tudatpy.astro.fundamentals)
yacma_extension(tudatpy.astro.gravitation)
yacma_extension(tudatpy.astro.polyhedron_utilities)
yacma_extension(tudatpy.astro.time_conversion)
yacma_extension(tudatpy.astro.two_body_dynamics)

# Constants
yacma_extension(tudatpy.constants)

# Math
yacma_extension(tudatpy.math.geometry)
yacma_extension(tudatpy.math.interpolators)
yacma_extension(tudatpy.math.numerical_integrators)
yacma_extension(tudatpy.math.root_finders)

# Numerical simulation
yacma_extension(tudatpy.numerical_simulation)
yacma_extension(tudatpy.numerical_simulation.propagation)
yacma_extension(tudatpy.numerical_simulation.environment)
yacma_extension(tudatpy.numerical_simulation.estimation)

# Numerical simulation - Environment setup
yacma_extension(tudatpy.numerical_simulation.environment_setup)
yacma_extension(tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients)
yacma_extension(tudatpy.numerical_simulation.environment_setup.atmosphere)
yacma_extension(tudatpy.numerical_simulation.environment_setup.ephemeris)
yacma_extension(tudatpy.numerical_simulation.environment_setup.gravity_field)
yacma_extension(tudatpy.numerical_simulation.environment_setup.gravity_field_variation)
yacma_extension(tudatpy.numerical_simulation.environment_setup.ground_station)
yacma_extension(tudatpy.numerical_simulation.environment_setup.radiation_pressure)
yacma_extension(tudatpy.numerical_simulation.environment_setup.rigid_body)
yacma_extension(tudatpy.numerical_simulation.environment_setup.rotation_model)
yacma_extension(tudatpy.numerical_simulation.environment_setup.shape)
yacma_extension(tudatpy.numerical_simulation.environment_setup.shape_deformation)
yacma_extension(tudatpy.numerical_simulation.environment_setup.vehicle_systems)

# Numerical simulation - Propagation setup
yacma_extension(tudatpy.numerical_simulation.propagation_setup)
yacma_extension(tudatpy.numerical_simulation.propagation_setup.acceleration)
yacma_extension(tudatpy.numerical_simulation.propagation_setup.dependent_variable)
yacma_extension(tudatpy.numerical_simulation.propagation_setup.integrator)
yacma_extension(tudatpy.numerical_simulation.propagation_setup.mass_rate)
yacma_extension(tudatpy.numerical_simulation.propagation_setup.propagator)
yacma_extension(tudatpy.numerical_simulation.propagation_setup.torque)
yacma_extension(tudatpy.numerical_simulation.propagation_setup.thrust)

# Numerical simulation - Estimation setup
yacma_extension(tudatpy.numerical_simulation.estimation_setup)
yacma_extension(tudatpy.numerical_simulation.estimation_setup.observation)
yacma_extension(tudatpy.numerical_simulation.estimation_setup.parameter)

# Trajectory design
yacma_extension(tudatpy.trajectory_design.shape_based_thrust)
yacma_extension(tudatpy.trajectory_design.transfer_trajectory)

# Install
######################################################
install(
    DIRECTORY ${TUDATPY_SOURCE_DIR}
    DESTINATION ${YACMA_PYTHON_MODULES_INSTALL_PATH}
)
