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
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${TUDATPY_SOURCE_DIR})

# Get sources for extension
set(sources ${TUDATPY_SOURCE_DIR}/kernel.cpp)
update_sources(math sources)
update_sources(astro sources)
update_sources(trajectory_design sources)
update_sources(constants sources)
update_sources(interface sources)
update_sources(data sources)
update_sources(numerical_simulation sources)

YACMA_PYTHON_MODULE(kernel ${sources})

# YACMA_PYTHON_MODULE(
#     kernel

#     ${TUDATPY_SOURCE_DIR}/kernel.cpp

#     # Astro
#     ${TUDATPY_SOURCE_DIR}/astro/astro.cpp
#     ${TUDATPY_SOURCE_DIR}/astro/gravitation/gravitation.cpp
#     ${TUDATPY_SOURCE_DIR}/astro/element_conversion/element_conversion.cpp
#     ${TUDATPY_SOURCE_DIR}/astro/frame_conversion/frame_conversion.cpp
#     ${TUDATPY_SOURCE_DIR}/astro/time_conversion/time_conversion.cpp
#     ${TUDATPY_SOURCE_DIR}/astro/two_body_dynamics/two_body_dynamics.cpp
#     ${TUDATPY_SOURCE_DIR}/astro/fundamentals/fundamentals.cpp
#     ${TUDATPY_SOURCE_DIR}/astro/polyhedron_utilities/polyhedron_utilities.cpp

#     # Trajectory design
#     ${TUDATPY_SOURCE_DIR}/trajectory_design/trajectory_design.cpp
#     ${TUDATPY_SOURCE_DIR}/trajectory_design/transfer_trajectory/transfer_trajectory.cpp
#     ${TUDATPY_SOURCE_DIR}/trajectory_design/shape_based_thrust/shape_based_thrust.cpp

#     # Constants
#     ${TUDATPY_SOURCE_DIR}/constants/constants.cpp

#     # Interface
#     ${TUDATPY_SOURCE_DIR}/interface/interface.cpp
#     ${TUDATPY_SOURCE_DIR}/interface/spice/spice.cpp

#     # Data
#     ${TUDATPY_SOURCE_DIR}/data/data.cpp

#     # Math
#     ${TUDATPY_SOURCE_DIR}/math/math.cpp
#     ${TUDATPY_SOURCE_DIR}/math/interpolators/interpolators.cpp
#     ${TUDATPY_SOURCE_DIR}/math/numerical_integrators/numerical_integrators.cpp
#     ${TUDATPY_SOURCE_DIR}/math/root_finders/root_finders.cpp
#     ${TUDATPY_SOURCE_DIR}/math/geometry/geometry.cpp
#     ${TUDATPY_SOURCE_DIR}/math/statistics/statistics.cpp

#     # Numerical simulation
#     ${TUDATPY_SOURCE_DIR}/numerical_simulation/propagation/propagation.cpp
#     ${TUDATPY_SOURCE_DIR}/numerical_simulation/environment/environment.cpp
#     ${TUDATPY_SOURCE_DIR}/numerical_simulation/estimation/estimation.cpp
#     ${TUDATPY_SOURCE_DIR}/numerical_simulation/estimation/estimation_filter_parser.cpp
#     ${TUDATPY_SOURCE_DIR}/numerical_simulation/estimation/estimation_observation_collection.cpp
#     ${TUDATPY_SOURCE_DIR}/numerical_simulation/estimation/estimation_propagated_covariance.cpp
#     ${TUDATPY_SOURCE_DIR}/numerical_simulation/estimation/estimation_single_observation_set.cpp

#     ${TUDATPY_SOURCE_DIR}/numerical_simulation/environment_setup/environment_setup.cpp
#     ${TUDATPY_SOURCE_DIR}/numerical_simulation/environment_setup/aerodynamic_coefficients/aerodynamic_coefficients.cpp
#     ${TUDATPY_SOURCE_DIR}/numerical_simulation/environment_setup/atmosphere/atmosphere.cpp
#     ${TUDATPY_SOURCE_DIR}/numerical_simulation/environment_setup/ephemeris/ephemeris.cpp
#     ${TUDATPY_SOURCE_DIR}/numerical_simulation/environment_setup/gravity_field/gravity_field.cpp
#     ${TUDATPY_SOURCE_DIR}/numerical_simulation/environment_setup/gravity_field_variation/gravity_field_variation.cpp
#     ${TUDATPY_SOURCE_DIR}/numerical_simulation/environment_setup/radiation_pressure/radiation_pressure.cpp
#     ${TUDATPY_SOURCE_DIR}/numerical_simulation/environment_setup/shape_deformation/shape_deformation.cpp
#     ${TUDATPY_SOURCE_DIR}/numerical_simulation/environment_setup/rotation_model/rotation_model.cpp
#     ${TUDATPY_SOURCE_DIR}/numerical_simulation/environment_setup/shape/shape.cpp
#     ${TUDATPY_SOURCE_DIR}/numerical_simulation/environment_setup/ground_station/ground_station.cpp
#     ${TUDATPY_SOURCE_DIR}/numerical_simulation/environment_setup/rigid_body/rigid_body.cpp
#     ${TUDATPY_SOURCE_DIR}/numerical_simulation/environment_setup/vehicle_systems/vehicle_systems.cpp

#     ${TUDATPY_SOURCE_DIR}/numerical_simulation/propagation_setup/propagation_setup.cpp
#     ${TUDATPY_SOURCE_DIR}/numerical_simulation/propagation_setup/acceleration/acceleration.cpp
#     ${TUDATPY_SOURCE_DIR}/numerical_simulation/propagation_setup/dependent_variable/dependent_variable.cpp
#     ${TUDATPY_SOURCE_DIR}/numerical_simulation/propagation_setup/integrator/integrator.cpp
#     ${TUDATPY_SOURCE_DIR}/numerical_simulation/propagation_setup/mass_rate/mass_rate.cpp
#     ${TUDATPY_SOURCE_DIR}/numerical_simulation/propagation_setup/propagator/propagator.cpp
#     ${TUDATPY_SOURCE_DIR}/numerical_simulation/propagation_setup/torque/torque.cpp
#     ${TUDATPY_SOURCE_DIR}/numerical_simulation/propagation_setup/thrust/thrust.cpp

#     ${TUDATPY_SOURCE_DIR}/numerical_simulation/estimation_setup/estimation_setup.cpp
#     ${TUDATPY_SOURCE_DIR}/numerical_simulation/estimation_setup/observation/observation.cpp
#     ${TUDATPY_SOURCE_DIR}/numerical_simulation/estimation_setup/parameter/parameter.cpp

#     ${TUDATPY_SOURCE_DIR}/numerical_simulation/numerical_simulation.cpp
#     ${TUDATPY_SOURCE_DIR}/numerical_simulation/numerical_simulation_simulator.cpp
#     ${TUDATPY_SOURCE_DIR}/numerical_simulation/numerical_simulation_variational.cpp
#     ${TUDATPY_SOURCE_DIR}/numerical_simulation/numerical_simulation_estimator.cpp

# )


# Link against libraries
target_link_libraries(kernel PRIVATE
    ${Boost_LIBRARIES}
    ${Boost_SYSTEM_LIBRARY}
    ${Tudat_PROPAGATION_LIBRARIES}
    ${Tudat_ESTIMATION_LIBRARIES}
)

# Include headers
target_include_directories(kernel PUBLIC
    $<BUILD_INTERFACE:${Boost_INCLUDE_DIRS}>
    $<BUILD_INTERFACE:${Tudat_INCLUDE_DIRS}>
    $<BUILD_INTERFACE:${EIGEN3_INCLUDE_DIRS}>
    $<BUILD_INTERFACE:${CSpice_INCLUDE_DIRS}>
    $<BUILD_INTERFACE:${Sofa_INCLUDE_DIRS}>
    $<INSTALL_INTERFACE:include>)
    target_include_directories(kernel SYSTEM PRIVATE "${pybind11_INCLUDE_DIR}")
    target_include_directories(kernel SYSTEM PRIVATE "${EIGEN3_INCLUDE_DIRS}")
    target_include_directories(kernel SYSTEM PRIVATE "${CSpice_INCLUDE_DIRS}")
    target_include_directories(kernel SYSTEM PRIVATE "${Sofa_INCLUDE_DIRS}")
    target_include_directories(kernel SYSTEM PRIVATE "${Tudat_INCLUDE_DIRS}")
    target_compile_definitions(kernel PRIVATE "${pybind11_DEFINITIONS}")
    set_target_properties(kernel PROPERTIES CXX_VISIBILITY_PRESET hidden)
    set_target_properties(kernel PROPERTIES VISIBILITY_INLINES_HIDDEN TRUE)

unset(CMAKE_LIBRARY_OUTPUT_DIRECTORY)

install(
    TARGETS kernel
    RUNTIME DESTINATION ${TUDATPY_INSTALL_PATH}
    LIBRARY DESTINATION ${TUDATPY_INSTALL_PATH}
)
# # Data
# add_extension(data)

# # Astro
# add_extension(astro 1)
# # add_extension(astro.element_conversion)
# # add_extension(astro.ephemerides)
# # add_extension(astro.frame_conversion)
# # add_extension(astro.fundamentals)
# # add_extension(astro.gravitation)
# # add_extension(astro.polyhedron_utilities)
# # add_extension(astro.time_conversion)
# # add_extension(astro.two_body_dynamics)

# # Constants
# add_extension(constants)

# # Math
# add_extension(math 1)
# # add_extension(math.geometry)
# # add_extension(math.interpolators)
# # add_extension(math.numerical_integrators)
# # add_extension(math.root_finders)
# # add_extension(math.statistics)

# # Numerical simulation
# add_extension(numerical_simulation 1)
# # add_extension(numerical_simulation.propagation)
# # add_extension(numerical_simulation.environment)
# # add_extension(numerical_simulation.estimation)

# # # Numerical simulation - Environment setup
# # add_extension(numerical_simulation.environment_setup)
# # add_extension(numerical_simulation.environment_setup.aerodynamic_coefficients)
# # add_extension(numerical_simulation.environment_setup.atmosphere)
# # add_extension(numerical_simulation.environment_setup.ephemeris)
# # add_extension(numerical_simulation.environment_setup.gravity_field)
# # add_extension(numerical_simulation.environment_setup.gravity_field_variation)
# # add_extension(numerical_simulation.environment_setup.ground_station)
# # add_extension(numerical_simulation.environment_setup.radiation_pressure)
# # add_extension(numerical_simulation.environment_setup.rigid_body)
# # add_extension(numerical_simulation.environment_setup.rotation_model)
# # add_extension(numerical_simulation.environment_setup.shape)
# # add_extension(numerical_simulation.environment_setup.shape_deformation)
# # add_extension(numerical_simulation.environment_setup.vehicle_systems)

# # # Numerical simulation - Propagation setup
# # add_extension(numerical_simulation.propagation_setup)
# # add_extension(numerical_simulation.propagation_setup.acceleration)
# # add_extension(numerical_simulation.propagation_setup.dependent_variable)
# # add_extension(numerical_simulation.propagation_setup.integrator)
# # add_extension(numerical_simulation.propagation_setup.mass_rate)
# # add_extension(numerical_simulation.propagation_setup.propagator)
# # add_extension(numerical_simulation.propagation_setup.thrust)
# # add_extension(numerical_simulation.propagation_setup.torque)

# # # Numerical simulation - Estimation setup
# # add_extension(numerical_simulation.estimation_setup)
# # add_extension(numerical_simulation.estimation_setup.observation)
# # add_extension(numerical_simulation.estimation_setup.parameter)

# # Trajectory design
# # add_extension(trajectory_design.shape_based_thrust)
# # add_extension(trajectory_design.transfer_trajectory)
# add_extension(trajectory_design 1)

# Install
######################################################
install(
    DIRECTORY ${TUDATPY_SOURCE_DIR}
    DESTINATION ${YACMA_PYTHON_MODULES_INSTALL_PATH}
)
