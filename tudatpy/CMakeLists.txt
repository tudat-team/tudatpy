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
set(TUDATPY_PYTHON_FILES
        __init__.py
#        elements.py
        prototype.py
#        io.py
        )

# Copy the python files in the current binary dir,
# so that we can import pygmo from the build dir.
# NOTE: importing from the build dir will work
# only on single-configuration generators.
foreach (TUDATPY_PYTHON_FILE ${TUDATPY_PYTHON_FILES})
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/${TUDATPY_PYTHON_FILE}"
            "${CMAKE_CURRENT_BINARY_DIR}/${TUDATPY_PYTHON_FILE}" COPYONLY)
endforeach ()


# Core module.
YACMA_PYTHON_MODULE(kernel
        # kernel
        kernel/kernel.cpp

        # kernel/utils
        kernel/expose_utils.cpp
        kernel/expose_utils/expose_data.cpp


        # kernel/astro
        kernel/expose_astro.cpp
        kernel/expose_astro/expose_gravitation.cpp
        kernel/expose_astro/expose_element_conversion.cpp
        kernel/expose_astro/expose_frame_conversion.cpp
        kernel/expose_astro/expose_time_conversion.cpp
        kernel/expose_astro/expose_two_body_dynamics.cpp
        kernel/expose_astro/expose_fundamentals.cpp
        kernel/expose_astro/expose_polyhedron_utilities.cpp

        # kernel/trajectory_design
        kernel/expose_trajectory_design.cpp
        kernel/expose_trajectory_design/expose_transfer_trajectory.cpp
        kernel/expose_trajectory_design/expose_shape_based_thrust.cpp

        # kernel/constants
        kernel/expose_constants.cpp

        # kernel/interface
        kernel/expose_interface.cpp
        kernel/expose_interface/expose_spice.cpp

        #        kernel/expose_interface/expose_json_interface.cpp
        #        kernel/expose_interface/expose_sofa_interface.cpp

        # kernel/io
        kernel/expose_io.cpp

        # kernel/math
        kernel/expose_math.cpp
        kernel/expose_math/expose_interpolators.cpp
        kernel/expose_math/expose_numerical_integrators.cpp
        kernel/expose_math/expose_root_finders.cpp
        kernel/expose_math/expose_geometry.cpp

        kernel/expose_numerical_simulation.cpp

        kernel/expose_numerical_simulation/expose_propagation.cpp
        kernel/expose_numerical_simulation/expose_environment.cpp
        kernel/expose_numerical_simulation/expose_estimation.cpp

        kernel/expose_numerical_simulation/expose_environment_setup.cpp
        kernel/expose_numerical_simulation/expose_environment_setup/expose_aerodynamic_coefficient_setup.cpp
        kernel/expose_numerical_simulation/expose_environment_setup/expose_atmosphere_setup.cpp
        kernel/expose_numerical_simulation/expose_environment_setup/expose_ephemeris_setup.cpp
        kernel/expose_numerical_simulation/expose_environment_setup/expose_gravity_field_setup.cpp
        kernel/expose_numerical_simulation/expose_environment_setup/expose_gravity_field_variation_setup.cpp
        kernel/expose_numerical_simulation/expose_environment_setup/expose_radiation_pressure_setup.cpp
        kernel/expose_numerical_simulation/expose_environment_setup/expose_shape_deformation_setup.cpp
        kernel/expose_numerical_simulation/expose_environment_setup/expose_rotation_model_setup.cpp
        kernel/expose_numerical_simulation/expose_environment_setup/expose_shape_setup.cpp
        kernel/expose_numerical_simulation/expose_environment_setup/expose_ground_station_setup.cpp
        kernel/expose_numerical_simulation/expose_environment_setup/expose_rigid_body_setup.cpp
        kernel/expose_numerical_simulation/expose_environment_setup/expose_vehicle_systems_setup.cpp


        kernel/expose_numerical_simulation/expose_propagation_setup.cpp
        kernel/expose_numerical_simulation/expose_propagation_setup/expose_acceleration_setup.cpp
        kernel/expose_numerical_simulation/expose_propagation_setup/expose_dependent_variable_setup.cpp
        kernel/expose_numerical_simulation/expose_propagation_setup/expose_integrator_setup.cpp
        kernel/expose_numerical_simulation/expose_propagation_setup/expose_mass_rate_setup.cpp
        kernel/expose_numerical_simulation/expose_propagation_setup/expose_propagator_setup.cpp
        kernel/expose_numerical_simulation/expose_propagation_setup/expose_torque_setup.cpp
        kernel/expose_numerical_simulation/expose_propagation_setup/expose_thrust_setup.cpp

        kernel/expose_numerical_simulation/expose_estimation_setup.cpp
        kernel/expose_numerical_simulation/expose_estimation_setup/expose_observation_setup.cpp
        kernel/expose_numerical_simulation/expose_estimation_setup/expose_estimated_parameter_setup.cpp

#        kernel/expose_unit_tests.cpp

#        kernel/expose_example.cpp
        )

target_link_libraries(kernel PRIVATE
        ${Boost_LIBRARIES}
        ${Boost_SYSTEM_LIBRARY}
        ${Tudat_PROPAGATION_LIBRARIES}
        ${Tudat_ESTIMATION_LIBRARIES}
        )

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

# Setup the installation path.
set(TUDATPY_INSTALL_PATH "${YACMA_PYTHON_MODULES_INSTALL_PATH}/tudatpy")

# Install the kernel module.
install(TARGETS kernel
        RUNTIME DESTINATION ${TUDATPY_INSTALL_PATH}
        LIBRARY DESTINATION ${TUDATPY_INSTALL_PATH}
        )

# Install Python files.
install(FILES ${TUDATPY_PYTHON_FILES} "${CMAKE_CURRENT_BINARY_DIR}/_version.py"
        DESTINATION ${TUDATPY_INSTALL_PATH})

# Install hybrid C++/Python modules
macro(INSTALL_HYBRID_CPP_PYTHON_MODULE module_dir_name)
    install(DIRECTORY ${module_dir_name}/ DESTINATION "${TUDATPY_INSTALL_PATH}/${module_dir_name}")
    file(COPY ${module_dir_name}/ DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/${module_dir_name}/)
endmacro()

# Hybrid modules
INSTALL_HYBRID_CPP_PYTHON_MODULE(astro)
INSTALL_HYBRID_CPP_PYTHON_MODULE(trajectory_design)
INSTALL_HYBRID_CPP_PYTHON_MODULE(constants)
INSTALL_HYBRID_CPP_PYTHON_MODULE(interface)
INSTALL_HYBRID_CPP_PYTHON_MODULE(math)
INSTALL_HYBRID_CPP_PYTHON_MODULE(numerical_simulation)
# Pure Python modules
INSTALL_HYBRID_CPP_PYTHON_MODULE(apps)
INSTALL_HYBRID_CPP_PYTHON_MODULE(bodies)
INSTALL_HYBRID_CPP_PYTHON_MODULE(plotting)
INSTALL_HYBRID_CPP_PYTHON_MODULE(data)
INSTALL_HYBRID_CPP_PYTHON_MODULE(util)
INSTALL_HYBRID_CPP_PYTHON_MODULE(io)

#install(FILES test.py "${CMAKE_CURRENT_BINARY_DIR}/_version.py")
file(COPY ../../tests DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../../)

# collect all example-X.py files
file(GLOB EXAMPLE_SCRIPTS "../../examples/example*.py")

# copy them into the root build tree
file(COPY ${EXAMPLE_SCRIPTS} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../../)
