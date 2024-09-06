# Update version
configure_file(
    "${TUDATPY_SOURCE_DIR}/_version.py.in"
    "${TUDATPY_SOURCE_DIR}/src/tudatpy/__init__.py"
    @ONLY
)

# Include libraries and add_extension macro
include(TudatpyLibraries)
include(AddExtension)

# # The list of pygmo's Python files.
# set(TUDATPY_PYTHON_FILES
#         __init__.py
# #        elements.py
#         prototype.py
# #        io.py
#         )

# # Copy the python files in the current binary dir,
# # so that we can import pygmo from the build dir.
# # NOTE: importing from the build dir will work
# # only on single-configuration generators.
# foreach (TUDATPY_PYTHON_FILE ${TUDATPY_PYTHON_FILES})
#     configure_file("${CMAKE_CURRENT_SOURCE_DIR}/${TUDATPY_PYTHON_FILE}"
#             "${CMAKE_CURRENT_BINARY_DIR}/${TUDATPY_PYTHON_FILE}" COPYONLY)
# endforeach ()

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

# # Core module.
# pybind11_add_module(kernel
#         # kernel
#         kernel/kernel.cpp

#         # kernel/utils
#         kernel/expose_utils.cpp
#         kernel/expose_utils/expose_data.cpp


#         # kernel/astro
#         kernel/expose_astro.cpp
#         kernel/expose_astro/expose_gravitation.cpp
#         kernel/expose_astro/expose_element_conversion.cpp
#         kernel/expose_astro/expose_frame_conversion.cpp
#         kernel/expose_astro/expose_time_conversion.cpp
#         kernel/expose_astro/expose_two_body_dynamics.cpp
#         kernel/expose_astro/expose_fundamentals.cpp
#         kernel/expose_astro/expose_polyhedron_utilities.cpp

#         # kernel/trajectory_design
#         kernel/expose_trajectory_design.cpp
#         kernel/expose_trajectory_design/expose_transfer_trajectory.cpp
#         kernel/expose_trajectory_design/expose_shape_based_thrust.cpp

#         # kernel/constants
#         kernel/expose_constants.cpp

#         # kernel/interface
#         kernel/expose_interface.cpp
#         kernel/expose_interface/expose_spice.cpp

#         #        kernel/expose_interface/expose_json_interface.cpp
#         #        kernel/expose_interface/expose_sofa_interface.cpp

#         # kernel/io
#         kernel/expose_io.cpp

#         # kernel/math
#         kernel/expose_math.cpp
#         kernel/expose_math/expose_interpolators.cpp
#         kernel/expose_math/expose_numerical_integrators.cpp
#         kernel/expose_math/expose_root_finders.cpp
#         kernel/expose_math/expose_geometry.cpp

#         kernel/expose_numerical_simulation.cpp

#         kernel/expose_numerical_simulation/expose_propagation.cpp
#         kernel/expose_numerical_simulation/expose_environment.cpp
#         kernel/expose_numerical_simulation/expose_estimation.cpp

#         kernel/expose_numerical_simulation/expose_environment_setup.cpp
#         kernel/expose_numerical_simulation/expose_environment_setup/expose_aerodynamic_coefficient_setup.cpp
#         kernel/expose_numerical_simulation/expose_environment_setup/expose_atmosphere_setup.cpp
#         kernel/expose_numerical_simulation/expose_environment_setup/expose_ephemeris_setup.cpp
#         kernel/expose_numerical_simulation/expose_environment_setup/expose_gravity_field_setup.cpp
#         kernel/expose_numerical_simulation/expose_environment_setup/expose_gravity_field_variation_setup.cpp
#         kernel/expose_numerical_simulation/expose_environment_setup/expose_radiation_pressure_setup.cpp
#         kernel/expose_numerical_simulation/expose_environment_setup/expose_shape_deformation_setup.cpp
#         kernel/expose_numerical_simulation/expose_environment_setup/expose_rotation_model_setup.cpp
#         kernel/expose_numerical_simulation/expose_environment_setup/expose_shape_setup.cpp
#         kernel/expose_numerical_simulation/expose_environment_setup/expose_ground_station_setup.cpp
#         kernel/expose_numerical_simulation/expose_environment_setup/expose_rigid_body_setup.cpp
#         kernel/expose_numerical_simulation/expose_environment_setup/expose_vehicle_systems_setup.cpp


#         kernel/expose_numerical_simulation/expose_propagation_setup.cpp
#         kernel/expose_numerical_simulation/expose_propagation_setup/expose_acceleration_setup.cpp
#         kernel/expose_numerical_simulation/expose_propagation_setup/expose_dependent_variable_setup.cpp
#         kernel/expose_numerical_simulation/expose_propagation_setup/expose_integrator_setup.cpp
#         kernel/expose_numerical_simulation/expose_propagation_setup/expose_mass_rate_setup.cpp
#         kernel/expose_numerical_simulation/expose_propagation_setup/expose_propagator_setup.cpp
#         kernel/expose_numerical_simulation/expose_propagation_setup/expose_torque_setup.cpp
#         kernel/expose_numerical_simulation/expose_propagation_setup/expose_thrust_setup.cpp

#         kernel/expose_numerical_simulation/expose_estimation_setup.cpp
#         kernel/expose_numerical_simulation/expose_estimation_setup/expose_observation_setup.cpp
#         kernel/expose_numerical_simulation/expose_estimation_setup/expose_estimated_parameter_setup.cpp

# #        kernel/expose_unit_tests.cpp

# #        kernel/expose_example.cpp
#         )

# target_link_libraries(kernel PRIVATE
#     ${TUDATPY_INTERNAL_LIBRARIES}
#     ${TUDATPY_EXTERNAL_LIBRARIES}
# )

# target_link_libraries(kernel PRIVATE
#         ${Boost_LIBRARIES}
#         ${Boost_SYSTEM_LIBRARY}
#         ${Tudat_PROPAGATION_LIBRARIES}
#         ${Tudat_ESTIMATION_LIBRARIES}
#         )

# target_include_directories(kernel PUBLIC
#         $<BUILD_INTERFACE:${Boost_INCLUDE_DIRS}>
#         $<BUILD_INTERFACE:${Tudat_INCLUDE_DIRS}>
#         $<BUILD_INTERFACE:${EIGEN3_INCLUDE_DIRS}>
#         $<BUILD_INTERFACE:${CSpice_INCLUDE_DIRS}>
#         $<BUILD_INTERFACE:${Sofa_INCLUDE_DIRS}>
#         $<INSTALL_INTERFACE:include>)

# target_include_directories(kernel SYSTEM PRIVATE "${pybind11_INCLUDE_DIR}")
# target_include_directories(kernel SYSTEM PRIVATE "${EIGEN3_INCLUDE_DIRS}")
# target_include_directories(kernel SYSTEM PRIVATE "${CSpice_INCLUDE_DIRS}")
# target_include_directories(kernel SYSTEM PRIVATE "${Sofa_INCLUDE_DIRS}")
# target_include_directories(kernel SYSTEM PRIVATE "${Tudat_INCLUDE_DIRS}")
# target_compile_definitions(kernel PRIVATE "${pybind11_DEFINITIONS}")
# set_target_properties(kernel PROPERTIES CXX_VISIBILITY_PRESET hidden)
# set_target_properties(kernel PROPERTIES VISIBILITY_INLINES_HIDDEN TRUE)

# Installation
set(TUDATPY_INSTALL_PATH "${CONDA_PREFIX}")
install(
    DIRECTORY "${TUDATPY_SOURCE_DIR}/src/tudatpy"
    DESTINATION "${TUDATPY_INSTALL_PATH}/tudatpy"
)

# # Install the kernel module.
# install(TARGETS kernel
#         RUNTIME DESTINATION ${TUDATPY_INSTALL_PATH}
#         LIBRARY DESTINATION ${TUDATPY_INSTALL_PATH}
#         )

# # Install Python files.
# install(FILES ${TUDATPY_PYTHON_FILES} "${CMAKE_CURRENT_BINARY_DIR}/_version.py"
#         DESTINATION ${TUDATPY_INSTALL_PATH})

# # Install hybrid C++/Python modules
# macro(INSTALL_HYBRID_CPP_PYTHON_MODULE module_dir_name)
#     install(DIRECTORY ${module_dir_name}/ DESTINATION "${TUDATPY_INSTALL_PATH}/${module_dir_name}")
#     file(COPY ${module_dir_name}/ DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/${module_dir_name}/)
# endmacro()

# # Hybrid modules
# INSTALL_HYBRID_CPP_PYTHON_MODULE(astro)
# INSTALL_HYBRID_CPP_PYTHON_MODULE(trajectory_design)
# INSTALL_HYBRID_CPP_PYTHON_MODULE(constants)
# INSTALL_HYBRID_CPP_PYTHON_MODULE(interface)
# INSTALL_HYBRID_CPP_PYTHON_MODULE(math)
# INSTALL_HYBRID_CPP_PYTHON_MODULE(numerical_simulation)
# # Pure Python modules
# INSTALL_HYBRID_CPP_PYTHON_MODULE(apps)
# INSTALL_HYBRID_CPP_PYTHON_MODULE(bodies)
# INSTALL_HYBRID_CPP_PYTHON_MODULE(plotting)
# INSTALL_HYBRID_CPP_PYTHON_MODULE(data)
# INSTALL_HYBRID_CPP_PYTHON_MODULE(util)
# INSTALL_HYBRID_CPP_PYTHON_MODULE(io)

#install(FILES test.py "${CMAKE_CURRENT_BINARY_DIR}/_version.py")
file(COPY ../../tests DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../../)

# collect all example-X.py files
file(GLOB EXAMPLE_SCRIPTS "../../examples/example*.py")

# copy them into the root build tree
file(COPY ${EXAMPLE_SCRIPTS} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../../)
