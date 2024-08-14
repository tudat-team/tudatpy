# Include list of libraries and add_extension macro
include(TudatpyLibraries)
include(AddExtension)

# # Generate _version.py
# configure_file(
#     "${CMAKE_CURRENT_SOURCE_DIR}/../../_version.py.in"
#     "${CMAKE_CURRENT_SOURCE_DIR}/__init__.py"
#     @ONLY
# )

# Add extensions
##################################

# Interface - Spice
# add_extension(tudatpy.interface.spice)

# # IO
# add_extension(tudatpy.io)

# # Astro
# add_extension(tudatpy.astro.element_conversion)
# add_extension(tudatpy.astro.frame_conversion)
# add_extension(tudatpy.astro.fundamentals)
# add_extension(tudatpy.astro.gravitation)
# add_extension(tudatpy.astro.polyhedron_utilities)
# add_extension(tudatpy.astro.time_conversion)
# add_extension(tudatpy.astro.two_body_dynamics)

# # Constants
# add_extension(tudatpy.constants)

# # Math
# add_extension(tudatpy.math.geometry)
# add_extension(tudatpy.math.interpolators)
# add_extension(tudatpy.math.numerical_integrators)
# add_extension(tudatpy.math.root_finders)

# # Numerical simulation
# add_extension(tudatpy.numerical_simulation)
# add_extension(tudatpy.numerical_simulation.propagation)
# add_extension(tudatpy.numerical_simulation.environment)
# add_extension(tudatpy.numerical_simulation.estimation)

# # Numerical simulation - Environment setup
# add_extension(tudatpy.numerical_simulation.environment_setup)
# add_extension(tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients)
# add_extension(tudatpy.numerical_simulation.environment_setup.atmosphere)
# add_extension(tudatpy.numerical_simulation.environment_setup.ephemeris)
# add_extension(tudatpy.numerical_simulation.environment_setup.gravity_field)
# add_extension(tudatpy.numerical_simulation.environment_setup.gravity_field_variation)
# add_extension(tudatpy.numerical_simulation.environment_setup.ground_station)
# add_extension(tudatpy.numerical_simulation.environment_setup.radiation_pressure)
# add_extension(tudatpy.numerical_simulation.environment_setup.rigid_body)
# add_extension(tudatpy.numerical_simulation.environment_setup.rotation_model)
# add_extension(tudatpy.numerical_simulation.environment_setup.shape)
# add_extension(tudatpy.numerical_simulation.environment_setup.shape_deformation)
# add_extension(tudatpy.numerical_simulation.environment_setup.vehicle_systems)

# # Numerical simulation - Propagation setup
# add_extension(tudatpy.numerical_simulation.propagation_setup)
# add_extension(tudatpy.numerical_simulation.propagation_setup.acceleration)
# add_extension(tudatpy.numerical_simulation.propagation_setup.dependent_variable)
# add_extension(tudatpy.numerical_simulation.propagation_setup.integrator)
# add_extension(tudatpy.numerical_simulation.propagation_setup.mass_rate)
# add_extension(tudatpy.numerical_simulation.propagation_setup.propagator)
# add_extension(tudatpy.numerical_simulation.propagation_setup.torque)
# add_extension(tudatpy.numerical_simulation.propagation_setup.thrust)

# # Numerical simulation - Estimation setup
# add_extension(tudatpy.numerical_simulation.estimation_setup)
# add_extension(tudatpy.numerical_simulation.estimation_setup.observation)
# add_extension(tudatpy.numerical_simulation.estimation_setup.parameter)

# # Trajectory design
# add_extension(tudatpy.trajectory_design.shape_based_thrust)
# add_extension(tudatpy.trajectory_design.transfer_trajectory)

# Generate stubs
##################################

if (EXISTS ${Python_SITEARCH}/${PROJECT_NAME})
    message(STATUS "Tudatpy is already installed")
else()
    message(STATUS "Installing tudatpy")
    file(
        CREATE_LINK
        ${TUDATPY_SOURCE_DIR}/src/${PROJECT_NAME}
        ${Python_SITEARCH}/${PROJECT_NAME}
        SYMBOLIC
    )
endif()

# Manual generation of stubs for __init__.py files
file(GLOB level_zero "${TUDATPY_SOURCE_DIR}/src/tudatpy/__init__.py")
file(GLOB level_one "${TUDATPY_SOURCE_DIR}/src/tudatpy/*/__init__.py")
file(GLOB level_two "${TUDATPY_SOURCE_DIR}/src/tudatpy/*/*/__init__.py")
file(GLOB level_three "${TUDATPY_SOURCE_DIR}/src/tudatpy/*/*/*/__init__.py")
file(GLOB level_four "${TUDATPY_SOURCE_DIR}/src/tudatpy/*/*/*/*/__init__.py")
list(
    APPEND init_files
    ${level_zero} ${level_one} ${level_two} ${level_three} ${level_four}
)
foreach(init_file ${init_files})
    path(RELATIVE_PATH relative_path ${TUDATPY_SOURCE_DIR}/src/tudatpy ${init_file})
    set(stub_path "${TUDATPY_SOURCE_DIR}/src/tudatpy-stubs/${relative_path}i")
    message(STATUS "Generating stubs for ${init_file}")
endforeach()
# if (REMOVE_AFTER_BUILD)
#     message(STATUS "Removing tudatpy symbolic link after build")
#     add_custom_command(
#         TARGET expose_spice POST_BUILD
#         COMMAND ${CMAKE_COMMAND} -E remove ${Python_SITEARCH}/${PROJECT_NAME}
#         COMMENT "Removing symbolic link..."
#     )
# endif()



# # Move stubs into tudatpy
# # - TODO: replace by install directory just like tudatpy
# # - TODO: turn into macro that can becalled to include stubs from particular kernel submodules, just like tudatpy python modules are selected for install
# add_custom_command(TARGET kernel POST_BUILD
#     COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_BINARY_DIR}/kernel ${CMAKE_CURRENT_BINARY_DIR} > /dev/null 2>&1
#     COMMENT "Moving stubs to tudatpy..."
# )

# # Remove tudatpy/__init__.pyi
# add_custom_command(TARGET kernel POST_BUILD
#     COMMAND ${CMAKE_COMMAND} -E remove "${CMAKE_CURRENT_BINARY_DIR}/__init__.pyi" > /dev/null 2>&1
#     COMMENT "Removing __init__.pyi file..."
# )

# # Remove directory (tudatpy/kernel) where stubgen generated the stubs
# add_custom_command(TARGET kernel POST_BUILD
#     COMMAND ${CMAKE_COMMAND} -E remove_directory "${CMAKE_CURRENT_BINARY_DIR}/kernel" > /dev/null 2>&1
#     COMMENT "Removing stub directory..."
# )
#################################################

# # collect all example-X.py files
# file(GLOB EXAMPLE_SCRIPTS "../examples/example*.py")

# # copy them into the root build tree
# file(COPY ${EXAMPLE_SCRIPTS} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../)
