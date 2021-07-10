## Configure the version file.
#configure_file("${CMAKE_CURRENT_SOURCE_DIR}/_version.py.in" "${CMAKE_CURRENT_BINARY_DIR}/_version.py" @ONLY)

if (MINGW OR ${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    message(STATUS "Creating the files for the installation from build/tudatpy.")
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/../tools/wheel_setup.py" "${CMAKE_CURRENT_BINARY_DIR}/../setup.py" @ONLY)
    if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
        # NOTE: this is necessary on linux but harmful on mingw.
        #        configure_file("${CMAKE_CURRENT_SOURCE_DIR}/../tools/wheel_setup.cfg" "${CMAKE_CURRENT_BINARY_DIR}/../setup.cfg" @ONLY)
    endif ()
    if (MINGW)
        # configure_file("${CMAKE_CURRENT_SOURCE_DIR}/../tools/mingw_wheel_libs_python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}.txt" "${CMAKE_CURRENT_BINARY_DIR}/wheel/mingw_wheel_libs_python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}.txt" @ONLY)
    endif ()
endif ()


# The list of pygmo's Python files.
set(TUDATPY_PYTHON_FILES
        __init__.py
        #_layer_propagation_setup.py
        elements.py
        prototype.py
        io.py
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

        # kernel/astro
        kernel/expose_astro.cpp
        kernel/expose_astro/expose_two_body_dynamics.cpp
        kernel/expose_astro/expose_fundamentals.cpp
        kernel/expose_astro/expose_gravitation.cpp
        kernel/expose_astro/expose_conversion.cpp
        kernel/expose_astro/expose_ephemerides.cpp
        kernel/expose_astro/expose_aerodynamics.cpp
        kernel/expose_astro/expose_trajectory_design.cpp
        kernel/expose_astro/expose_propagators.cpp
        kernel/expose_astro/expose_reference_frames.cpp
        kernel/expose_astro/expose_shape.cpp
        kernel/expose_astro/expose_observations.cpp
        kernel/expose_astro/expose_ground_stations.cpp


        # kernel/constants
        kernel/expose_constants.cpp

        # kernel/interface
        kernel/expose_interface.cpp
        kernel/expose_interface/expose_spice_interface.cpp
        #        kernel/expose_interface/expose_json_interface.cpp
        #        kernel/expose_interface/expose_sofa_interface.cpp

        # kernel/io
        kernel/expose_io.cpp

        # kernel/math
        kernel/expose_math.cpp
        kernel/expose_math/expose_basic_functions.cpp
        kernel/expose_math/expose_interpolators.cpp
        kernel/expose_math/expose_numerical_integrators.cpp
        kernel/expose_math/expose_root_finders.cpp
        kernel/expose_math/expose_geometry.cpp

        # kernel/simulation
        kernel/expose_simulation.cpp
        kernel/expose_simulation/expose_astro_setup.cpp
        kernel/expose_simulation/expose_environment_setup.cpp
        kernel/expose_simulation/expose_propagation_setup.cpp
        kernel/expose_simulation/expose_estimation_setup.cpp
        kernel/expose_simulation/expose_shape_based_thrust.cpp

        kernel/expose_unit_tests.cpp
#        kernel/docstrings.cpp

        kernel/expose_example.cpp
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

# Install the Python files.
install(FILES ${TUDATPY_PYTHON_FILES} "${CMAKE_CURRENT_BINARY_DIR}/_version.py"
        DESTINATION ${TUDATPY_INSTALL_PATH})

macro(INSTALL_PURE_PYTHON_MODULE_DIR module_dir_name)
    install(DIRECTORY ${module_dir_name}/ DESTINATION "${TUDATPY_INSTALL_PATH}/${module_dir_name}")
    file(COPY ${module_dir_name}/ DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/${module_dir_name}/)
endmacro()

INSTALL_PURE_PYTHON_MODULE_DIR(apps)
INSTALL_PURE_PYTHON_MODULE_DIR(bodies)
INSTALL_PURE_PYTHON_MODULE_DIR(plotting)
INSTALL_PURE_PYTHON_MODULE_DIR(util)

#install(FILES test.py "${CMAKE_CURRENT_BINARY_DIR}/_version.py")
file(COPY ../tests DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../)
file(COPY test.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../)
file(GLOB EXAMPLE_SCRIPTS "../examples/example*.py")
file(COPY ${EXAMPLE_SCRIPTS} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/../)
