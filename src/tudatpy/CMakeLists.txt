# Configure the version file.
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/_version.py.in" "${CMAKE_CURRENT_BINARY_DIR}/_version.py" @ONLY)

# Configure the files needed to make the python wheels (for PyPi packages)
#if(MINGW OR ${CMAKE_SYSTEM_NAME} MATCHES "Linux")
#    message(STATUS "Creating the files for the generation of a binary wheel.")
#    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/../tools/wheel_setup.py" "${CMAKE_CURRENT_BINARY_DIR}/../wheel/setup.py" @ONLY)
#    if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
#        # NOTE: this is necessary on linux but harmful on mingw.
#        configure_file("${CMAKE_CURRENT_SOURCE_DIR}/../tools/wheel_setup.cfg" "${CMAKE_CURRENT_BINARY_DIR}/../wheel/setup.cfg" @ONLY)
#    endif()
#    if(MINGW)
#        configure_file("${CMAKE_CURRENT_SOURCE_DIR}/../tools/mingw_wheel_libs_python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}.txt" "${CMAKE_CURRENT_BINARY_DIR}/wheel/mingw_wheel_libs_python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}.txt" @ONLY)
#    endif()
#endif()

# Offer the user the choice of overriding the installation directories
#set(INSTALL_LIB_DIR ${CMAKE_INSTALL_PREFIX}/lib CACHE PATH "Installation directory for libraries")
#set(INSTALL_BIN_DIR ${CMAKE_INSTALL_PREFIX}/bin CACHE PATH "Installation directory for executables")
#set(INSTALL_INCLUDE_DIR ${CMAKE_INSTALL_PREFIX}/include CACHE PATH "Installation directory for header files")
#set(INSTALL_DATA_DIR ${CMAKE_INSTALL_PREFIX}/data/ CACHE PATH "Installation directory for data files")
#set(INSTALL_MAN_DIR ${CMAKE_INSTALL_PREFIX} CACHE PATH "Installation directory for man documentation")
#set(INSTALL_CMAKE_DIR ${INSTALL_LIB_DIR}/cmake/cspice CACHE PATH "Installation directory for cmake config files")


# The list of pygmo's Python files.
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

macro (add_module_custom name location)

    # Update installation prefix
    set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/tudatpy/${location})

    # Add module
    pybind11_add_module(${name} MODULE ${ARGN})

    # Add headers
    target_include_directories(${name} PUBLIC
        $<BUILD_INTERFACE:${Boost_INCLUDE_DIRS}>
        $<BUILD_INTERFACE:${Tudat_INCLUDE_DIRS}>
        $<BUILD_INTERFACE:${EIGEN3_INCLUDE_DIRS}>
        $<BUILD_INTERFACE:${CSpice_INCLUDE_DIRS}>
        $<BUILD_INTERFACE:${Sofa_INCLUDE_DIRS}>
        $<INSTALL_INTERFACE:include>
    )
    target_include_directories(${name} SYSTEM PRIVATE
        "${pybind11_INCLUDE_DIR}"
        "${EIGEN3_INCLUDE_DIRS}"
        "${CSpice_INCLUDE_DIRS}"
        "${Sofa_INCLUDE_DIRS}"
        "${Tudat_INCLUDE_DIRS}"
    )

    # Link libraries
    target_link_libraries(${name} PRIVATE
        ${Boost_LIBRARIES}
        ${Boost_SYSTEM_LIBRARY}
        ${Tudat_PROPAGATION_LIBRARIES}
        ${Tudat_ESTIMATION_LIBRARIES}
    )

    # Additional properties
    target_compile_definitions(${name} PRIVATE "${pybind11_DEFINITIONS}")
    set_target_properties(${name} PROPERTIES CXX_VISIBILITY_PRESET hidden)
    set_target_properties(${name} PROPERTIES VISIBILITY_INLINES_HIDDEN TRUE)

    # Restore installation prefix
    unset (CMAKE_LIBRARY_OUTPUT_DIRECTORY)

endmacro()

# # Utils
# add_module_custom(
#     expose_utils
#     src/tudatpy/util
#     util/expose_utils.cpp
# )

# IO
add_module_custom(
    expose_io
    src/tudatpy/io
    io/expose_io.cpp
)

# Interface - Spice
add_module_custom(
    expose_spice
    src/tudatpy/interface/spice
    interface/spice/expose_spice.cpp
)

# Astro - Element conversion
add_module_custom(
    expose_element_conversion
    src/tudatpy/astro/element_conversion
    astro/element_conversion/expose_element_conversion.cpp
)

# Astro - Element conversion
add_module_custom(
    expose_frame_conversion
    src/tudatpy/astro/frame_conversion
    astro/frame_conversion/expose_frame_conversion.cpp
)

# Astro - Fundamentals
add_module_custom(
    expose_fundamentals
    src/tudatpy/astro/fundamentals
    astro/fundamentals/expose_fundamentals.cpp
)

# Astro - Gravitation
add_module_custom(
    expose_gravitation
    src/tudatpy/astro/gravitation
    astro/gravitation/expose_gravitation.cpp
)

# Astro - Polyhedron utilities
add_module_custom(
    expose_polyhedron_utilities
    src/tudatpy/astro/polyhedron_utilities
    astro/polyhedron_utilities/expose_polyhedron_utilities.cpp
)

# Astro - Time conversion
add_module_custom(
    expose_time_conversion
    src/tudatpy/astro/time_conversion
    astro/time_conversion/expose_time_conversion.cpp
)

# Astro - Two body dynamics
add_module_custom(
    expose_two_body_dynamics
    src/tudatpy/astro/two_body_dynamics
    "math/root_finders/expose_root_finders.cpp"
    astro/two_body_dynamics/expose_two_body_dynamics.cpp
)

# Trajectory design - Shape based thrust
add_module_custom(
    expose_shape_based_thrust
    src/tudatpy/trajectory_design/shape_based_thrust
    trajectory_design/shape_based_thrust/expose_shape_based_thrust.cpp
)

# Trajectory design - Transfer trajectory
add_module_custom(
    expose_transfer_trajectory
    src/tudatpy/trajectory_design/transfer_trajectory
    trajectory_design/transfer_trajectory/expose_transfer_trajectory.cpp
)

# Constants
add_module_custom(
    expose_constants
    src/tudatpy/constants
    constants/expose_constants.cpp
)

# Math - Geometry
add_module_custom(
    "expose_geometry"
    "src/tudatpy/math/geometry"
    "math/geometry/expose_geometry.cpp"
)

# Math - Interpolators
add_module_custom(
    "expose_interpolators"
    "src/tudatpy/math/interpolators"
    "math/interpolators/expose_interpolators.cpp"
)

# Math - Numerical integrators
add_module_custom(
    "expose_numerical_integrators"
    "src/tudatpy/math/numerical_integrators"
    "math/numerical_integrators/expose_numerical_integrators.cpp"
)

# Math - Root finders
add_module_custom(
    "expose_root_finders"
    "src/tudatpy/math/root_finders"
    "math/root_finders/expose_root_finders.cpp"
)

# Numerical simulation
add_module_custom(
    "expose_numerical_simulation"
    "src/tudatpy/numerical_simulation"
    "numerical_simulation/expose_numerical_simulation.cpp"
)

# Numerical simulation - Propagation
add_module_custom(
    "expose_propagation"
    "src/tudatpy/numerical_simulation/propagation"
    "numerical_simulation/propagation/expose_propagation.cpp"
)

# Numerical simulation - Environment
add_module_custom(
    "expose_environment"
    "src/tudatpy/numerical_simulation/environment"
    "numerical_simulation/environment/expose_environment.cpp"
)

# Numerical simulation - Estimation
add_module_custom(
    "expose_estimation"
    "src/tudatpy/numerical_simulation/estimation"
    "numerical_simulation/estimation/expose_estimation.cpp"
)

# Numerical simulation - Environment setup
add_module_custom(
    "expose_environment_setup"
    "src/tudatpy/numerical_simulation/environment_setup"
    "numerical_simulation/environment_setup/expose_environment_setup.cpp"
)

# Numerical simulation - Environment setup - Aerodynamic coefficients
add_module_custom(
    "expose_aerodynamic_coefficients"
    "src/tudatpy/numerical_simulation/environment_setup/aerodynamic_coefficients"
    "numerical_simulation/environment_setup/aerodynamic_coefficients/expose_aerodynamic_coefficients.cpp"
)

# Numerical simulation - Environment setup - Atmosphere
add_module_custom(
    "expose_atmosphere"
    "src/tudatpy/numerical_simulation/environment_setup/atmosphere"
    "numerical_simulation/environment_setup/atmosphere/expose_atmosphere.cpp"
)

# Numerical simulation - Environment setup - Ephemeris
add_module_custom(
    "expose_ephemeris"
    "src/tudatpy/numerical_simulation/environment_setup/ephemeris"
    "numerical_simulation/environment_setup/ephemeris/expose_ephemeris.cpp"
)

# Numerical simulation - Environment setup - Gravity field
add_module_custom(
    "expose_gravity_field"
    "src/tudatpy/numerical_simulation/environment_setup/gravity_field"
    "numerical_simulation/environment_setup/gravity_field/expose_gravity_field.cpp"
)

# Numerical simulation - Environment setup - Gravity field variation
add_module_custom(
    "expose_gravity_field_variation"
    "src/tudatpy/numerical_simulation/environment_setup/gravity_field_variation"
    "numerical_simulation/environment_setup/gravity_field_variation/expose_gravity_field_variation.cpp"
)

# Numerical simulation - Environment setup - Ground station
add_module_custom(
    "expose_ground_station"
    "src/tudatpy/numerical_simulation/environment_setup/ground_station"
    "numerical_simulation/environment_setup/ground_station/expose_ground_station.cpp"
)

# Numerical simulation - Environment setup - Radiation pressure
add_module_custom(
    "expose_radiation_pressure"
    "src/tudatpy/numerical_simulation/environment_setup/radiation_pressure"
    "numerical_simulation/environment_setup/radiation_pressure/expose_radiation_pressure.cpp"
)

# Numerical simulation - Environment setup - Rigid body
add_module_custom(
    "expose_rigid_body"
    "src/tudatpy/numerical_simulation/environment_setup/rigid_body"
    "numerical_simulation/environment_setup/rigid_body/expose_rigid_body.cpp"
)

# Numerical simulation - Environment setup - Rotation model
add_module_custom(
    "expose_rotation_model"
    "src/tudatpy/numerical_simulation/environment_setup/rotation_model"
    "numerical_simulation/environment_setup/rotation_model/expose_rotation_model.cpp"
)

# Numerical simulation - Environment setup - Shape
add_module_custom(
    "expose_shape"
    "src/tudatpy/numerical_simulation/environment_setup/shape"
    "numerical_simulation/environment_setup/shape/expose_shape.cpp"
)

# Numerical simulation - Environment setup - Shape deformation
add_module_custom(
    "expose_shape_deformation"
    "src/tudatpy/numerical_simulation/environment_setup/shape_deformation"
    "numerical_simulation/environment_setup/shape_deformation/expose_shape_deformation.cpp"
)

# Numerical simulation - Environment setup - Vehicle systems
add_module_custom(
    "expose_vehicle_systems"
    "src/tudatpy/numerical_simulation/environment_setup/vehicle_systems"
    "numerical_simulation/environment_setup/vehicle_systems/expose_vehicle_systems.cpp"
)

# Numerical simulation - Propagation setup
add_module_custom(
    "expose_propagation_setup"
    "src/tudatpy/numerical_simulation/propagation_setup"
    "numerical_simulation/propagation_setup/expose_propagation_setup.cpp"
)

# Numerical simulation - Propagation setup - Acceleration
add_module_custom(
    "expose_acceleration"
    "src/tudatpy/numerical_simulation/propagation_setup/acceleration"
    "numerical_simulation/propagation_setup/acceleration/expose_acceleration.cpp"
)

# Numerical simulation - Propagation setup - Dependent variable
add_module_custom(
    "expose_dependent_variable"
    "src/tudatpy/numerical_simulation/propagation_setup/dependent_variable"
    "numerical_simulation/propagation_setup/dependent_variable/expose_dependent_variable.cpp"
)

# Numerical simulation - Propagation setup - Integrator
add_module_custom(
    "expose_integrator"
    "src/tudatpy/numerical_simulation/propagation_setup/integrator"
    "numerical_simulation/propagation_setup/integrator/expose_integrator.cpp"
)

# Numerical simulation - Propagation setup - Mass rate
add_module_custom(
    "expose_mass_rate"
    "src/tudatpy/numerical_simulation/propagation_setup/mass_rate"
    "numerical_simulation/propagation_setup/mass_rate/expose_mass_rate.cpp"
)

# Numerical simulation - Propagation setup - Propagator
add_module_custom(
    "expose_propagator"
    "src/tudatpy/numerical_simulation/propagation_setup/propagator"
    "numerical_simulation/propagation_setup/propagator/expose_propagator.cpp"
)

# Numerical simulation - Propagation setup - Torque
add_module_custom(
    "expose_torque"
    "src/tudatpy/numerical_simulation/propagation_setup/torque"
    "numerical_simulation/propagation_setup/torque/expose_torque.cpp"
)

# Numerical simulation - Propagation setup - Thrust
add_module_custom(
    "expose_thrust"
    "src/tudatpy/numerical_simulation/propagation_setup/thrust"
    "numerical_simulation/propagation_setup/thrust/expose_thrust.cpp"
)

# Numerical simulation - Estimation setup
add_module_custom(
    "expose_estimation_setup"
    "src/tudatpy/numerical_simulation/estimation_setup"
    "numerical_simulation/estimation_setup/expose_estimation_setup.cpp"
)

# Numerical simulation - Estimation setup - Parameter
add_module_custom(
    "expose_parameter"
    "src/tudatpy/numerical_simulation/estimation_setup/parameter"
    "numerical_simulation/estimation_setup/parameter/expose_parameter.cpp"
)

# Numerical simulation - Estimation setup - Observation
add_module_custom(
    "expose_observation"
    "src/tudatpy/numerical_simulation/estimation_setup/observation"
    "numerical_simulation/estimation_setup/observation/expose_observation.cpp"
)

# # Kernel
# add_module_custom(kernel kernel
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

# # target_link_libraries(kernel PRIVATE
# #         ${Boost_LIBRARIES}
# #         ${Boost_SYSTEM_LIBRARY}
# #         ${Tudat_PROPAGATION_LIBRARIES}
# #         ${Tudat_ESTIMATION_LIBRARIES}
# #         )

# # target_include_directories(kernel PUBLIC
# #         $<BUILD_INTERFACE:${Boost_INCLUDE_DIRS}>
# #         $<BUILD_INTERFACE:${Tudat_INCLUDE_DIRS}>
# #         $<BUILD_INTERFACE:${EIGEN3_INCLUDE_DIRS}>
# #         $<BUILD_INTERFACE:${CSpice_INCLUDE_DIRS}>
# #         $<BUILD_INTERFACE:${Sofa_INCLUDE_DIRS}>
# #         $<INSTALL_INTERFACE:include>)

# # target_include_directories(kernel SYSTEM PRIVATE "${pybind11_INCLUDE_DIR}")
# # target_include_directories(kernel SYSTEM PRIVATE "${EIGEN3_INCLUDE_DIRS}")
# # target_include_directories(kernel SYSTEM PRIVATE "${CSpice_INCLUDE_DIRS}")
# # target_include_directories(kernel SYSTEM PRIVATE "${Sofa_INCLUDE_DIRS}")
# # target_include_directories(kernel SYSTEM PRIVATE "${Tudat_INCLUDE_DIRS}")
# # target_compile_definitions(kernel PRIVATE "${pybind11_DEFINITIONS}")
# # set_target_properties(kernel PROPERTIES CXX_VISIBILITY_PRESET hidden)
# # set_target_properties(kernel PROPERTIES VISIBILITY_INLINES_HIDDEN TRUE)

# # # Setup the installation path.
# # set(TUDATPY_INSTALL_PATH "${YACMA_PYTHON_MODULES_INSTALL_PATH}/tudatpy")

# # # Install the kernel module.
# # install(TARGETS kernel
# #         RUNTIME DESTINATION ${TUDATPY_INSTALL_PATH}
# #         LIBRARY DESTINATION ${TUDATPY_INSTALL_PATH}
# #         )

# # # Install Python files.
# # install(FILES ${TUDATPY_PYTHON_FILES} "${CMAKE_CURRENT_BINARY_DIR}/_version.py"
# #         DESTINATION ${TUDATPY_INSTALL_PATH})

# # # Install hybrid C++/Python modules
# # macro(INSTALL_HYBRID_CPP_PYTHON_MODULE module_dir_name)
# #     install(DIRECTORY ${module_dir_name}/ DESTINATION "${TUDATPY_INSTALL_PATH}/${module_dir_name}")
# #     file(COPY ${module_dir_name}/ DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/${module_dir_name}/)
# # endmacro()

# # # Hybrid modules
# # INSTALL_HYBRID_CPP_PYTHON_MODULE(astro)
# # INSTALL_HYBRID_CPP_PYTHON_MODULE(trajectory_design)
# # INSTALL_HYBRID_CPP_PYTHON_MODULE(constants)
# # INSTALL_HYBRID_CPP_PYTHON_MODULE(interface)
# # INSTALL_HYBRID_CPP_PYTHON_MODULE(math)
# # INSTALL_HYBRID_CPP_PYTHON_MODULE(numerical_simulation)
# # # Pure Python modules
# # INSTALL_HYBRID_CPP_PYTHON_MODULE(apps)
# # INSTALL_HYBRID_CPP_PYTHON_MODULE(bodies)
# # INSTALL_HYBRID_CPP_PYTHON_MODULE(plotting)
# # INSTALL_HYBRID_CPP_PYTHON_MODULE(data)
# # INSTALL_HYBRID_CPP_PYTHON_MODULE(util)
# # INSTALL_HYBRID_CPP_PYTHON_MODULE(io)

# # #install(FILES test.py "${CMAKE_CURRENT_BINARY_DIR}/_version.py")
# # file(COPY ../tests DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../)

# # # collect all example-X.py files
# # file(GLOB EXAMPLE_SCRIPTS "../examples/example*.py")

# # # copy them into the root build tree
# # file(COPY ${EXAMPLE_SCRIPTS} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../)
