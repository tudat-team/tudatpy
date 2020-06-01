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

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")

if (MINGW OR ${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    message(STATUS "Creating the files for the installation from build/tudatpy.")
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/../tools/wheel_setup.py" "${CMAKE_CURRENT_BINARY_DIR}/../setup.py" @ONLY)
    if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
        # NOTE: this is necessary on linux but harmful on mingw.
        #        configure_file("${CMAKE_CURRENT_SOURCE_DIR}/../tools/wheel_setup.cfg" "${CMAKE_CURRENT_BINARY_DIR}/../setup.cfg" @ONLY)
    endif ()
    if (MINGW)
        configure_file("${CMAKE_CURRENT_SOURCE_DIR}/../tools/mingw_wheel_libs_python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}.txt" "${CMAKE_CURRENT_BINARY_DIR}/wheel/mingw_wheel_libs_python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}.txt" @ONLY)
    endif ()
endif ()


# The list of pygmo's Python files.
set(TUDATPY_PYTHON_FILES
        __init__.py
        _layer_simulation_setup.py
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
YACMA_PYTHON_MODULE(core
        kernel/core.cpp
        kernel/expose_constants.cpp
        kernel/expose_root_finders.cpp
        kernel/expose_interpolators.cpp
        kernel/expose_spice_interface.cpp
        kernel/expose_ephemerides.cpp
        kernel/expose_reference_frames.cpp
        kernel/expose_aerodynamics.cpp
        kernel/expose_basic_astrodynamics.cpp
        kernel/expose_gravitation.cpp
        kernel/expose_numerical_integrators.cpp
        kernel/expose_orbital_element_conversions.cpp
        kernel/expose_propagators.cpp
        kernel/expose_simulation_setup.cpp
        kernel/expose_unit_tests.cpp
        #        expose_mathematics.cpp
        #        trampoline_classes.cpp
        kernel/docstrings.cpp
        #        expose_ephemeris.cpp
        )

target_link_libraries(core PRIVATE
        ${Boost_LIBRARIES}
        ${Boost_SYSTEM_LIBRARY}
        ${PYTHON_LIBRARIES}
        ${Tudat_PROPAGATION_LIBRARIES})

target_include_directories(core PUBLIC
        $<BUILD_INTERFACE:${pybind11_INCLUDE_DIR}>
        $<INSTALL_INTERFACE:include>)

target_include_directories(core PUBLIC
        $<BUILD_INTERFACE:${Tudat_INCLUDE_DIRS}>
        $<INSTALL_INTERFACE:include>)

target_include_directories(core PUBLIC
        $<BUILD_INTERFACE:${EIGEN3_INCLUDE_DIRS}>
        $<INSTALL_INTERFACE:include>)

target_include_directories(core PUBLIC
        $<BUILD_INTERFACE:${CSpice_INCLUDE_DIRS}>
        $<INSTALL_INTERFACE:include>)

target_include_directories(core PUBLIC
        $<BUILD_INTERFACE:${Sofa_INCLUDE_DIRS}>
        $<INSTALL_INTERFACE:include>)

target_include_directories(core SYSTEM PRIVATE "${pybind11_INCLUDE_DIR}")
target_compile_definitions(core PRIVATE "${pybind11_DEFINITIONS}")
target_compile_definitions(core PRIVATE "${Tudat_DEFINITIONS}")
set_target_properties(core PROPERTIES CXX_VISIBILITY_PRESET hidden)
set_target_properties(core PROPERTIES VISIBILITY_INLINES_HIDDEN TRUE)

# Setup the installation path.
set(TUDATPY_INSTALL_PATH "${YACMA_PYTHON_MODULES_INSTALL_PATH}/tudatpy")

# Install the core module.
install(TARGETS core
        RUNTIME DESTINATION ${TUDATPY_INSTALL_PATH}
        LIBRARY DESTINATION ${TUDATPY_INSTALL_PATH}
        )

# Install the Python files.
install(FILES ${TUDATPY_PYTHON_FILES} "${CMAKE_CURRENT_BINARY_DIR}/_version.py"
        DESTINATION ${TUDATPY_INSTALL_PATH})


#install(FILES test.py "${CMAKE_CURRENT_BINARY_DIR}/_version.py")
file(COPY test.py DESTINATION ${CMAKE_BINARY_DIR})