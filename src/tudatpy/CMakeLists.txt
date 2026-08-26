# Update version
# This should be replaced with setuptools_scm based versioning
configure_file(
    "${PROJECT_SOURCE_DIR}/configuration/tudatpy/__init__.py.in"
    "${PROJECT_SOURCE_DIR}/src/tudatpy/__init__.py"
    @ONLY
)

# include(YacmaExtension)
include(PythonExtension)

copy_python_in_build()

# Setup the installation path.
set(TUDATPY_INSTALL_PATH "${YACMA_PYTHON_MODULES_INSTALL_PATH}/tudatpy")
# set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${TUDATPY_SOURCE_DIR})

# Get sources for extension
set(sources
        ${TUDATPY_SOURCE_DIR}/kernel.cpp
        ${TUDATPY_SOURCE_DIR}/add_math_to_kernel.cpp
        ${TUDATPY_SOURCE_DIR}/add_astro_to_kernel.cpp
        ${TUDATPY_SOURCE_DIR}/add_trajectory_design_to_kernel.cpp
        ${TUDATPY_SOURCE_DIR}/add_constants_to_kernel.cpp
        ${TUDATPY_SOURCE_DIR}/add_interface_to_kernel.cpp
        ${TUDATPY_SOURCE_DIR}/add_data_to_kernel.cpp
        ${TUDATPY_SOURCE_DIR}/add_dynamics_to_kernel.cpp
        ${TUDATPY_SOURCE_DIR}/add_estimation_to_kernel.cpp
        ${TUDATPY_SOURCE_DIR}/add_exceptions_to_kernel.cpp
)
update_sources(math sources)
update_sources(astro sources)
update_sources(trajectory_design sources)
update_sources(constants sources)
update_sources(interface sources)
update_sources(data sources)
update_sources(numerical_simulation sources)
update_sources(dynamics sources)
update_sources(estimation sources)
update_sources(exceptions sources)

YACMA_PYTHON_MODULE(kernel ${sources})

if(TUDAT_BUILD_WITH_MCD_INTERFACE)
    target_link_libraries(kernel PRIVATE
        mcd_fortran
        tudat_mcd_interface
        ${NetCDF_Fortran_LIBRARIES}
        ${NetCDF_C_LIBRARIES}
    )
    target_include_directories(kernel SYSTEM PRIVATE 
        "${NetCDF_INCLUDE_DIRS}"
        "${CMAKE_SOURCE_DIR}/third_parties/mcd")
    message(STATUS "Tudatpy: MCD support linked to kernel")
endif()

# Link against libraries
target_link_libraries(kernel PRIVATE
    ${Boost_LIBRARIES}
    ${Boost_SYSTEM_LIBRARY}
    ${Tudat_PROPAGATION_LIBRARIES}
    ${Tudat_ESTIMATION_LIBRARIES}
)

# ---- Precompiled header for pybind11 bindings ----
# These headers are included by virtually every binding .cpp file.
# PCH eliminates redundant parsing across all 86+ translation units.
if (TUDAT_BUILD_WITH_PCH)
    target_precompile_headers(kernel PRIVATE
        <pybind11/pybind11.h>
        <pybind11/stl.h>
        <pybind11/eigen.h>
        <pybind11/functional.h>
        <pybind11/numpy.h>
        <pybind11/chrono.h>
        <Eigen/Core>
        <Eigen/Geometry>
        <memory>
        <vector>
        <string>
        <map>
        <functional>
    )
endif ()

# Disable unity build for binding target (each file exposes unique symbols)
set_target_properties(kernel PROPERTIES UNITY_BUILD OFF)

# Include headers
target_include_directories(kernel PUBLIC
    $<BUILD_INTERFACE:${Tudat_INCLUDE_DIRS}>
    $<INSTALL_INTERFACE:include>)

if (TUDATPY_BUILD_WITH_DETAILED_PYBIND11_ERRORS)
    target_compile_definitions(kernel PRIVATE TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS=1)
else ()
    target_compile_definitions(kernel PRIVATE TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS=0)
endif ()

    target_include_directories(kernel SYSTEM PRIVATE "${pybind11_INCLUDE_DIR}")
    target_include_directories(kernel SYSTEM PRIVATE "${Boost_INCLUDE_DIRS}")
    target_include_directories(kernel SYSTEM PRIVATE "${EIGEN3_INCLUDE_DIRS}")
    target_include_directories(kernel SYSTEM PRIVATE "${CSpice_INCLUDE_DIRS}")
    target_include_directories(kernel SYSTEM PRIVATE "${Sofa_INCLUDE_DIRS}")
    target_include_directories(kernel SYSTEM PRIVATE "${Tudat_INCLUDE_DIRS}")
    target_compile_definitions(kernel PRIVATE "${pybind11_DEFINITIONS}")
    if (APPLE AND TUDAT_BUILD_STATIC_LIBRARY AND TUDAT_BUILD_WITH_SERIALIZATION)
        set_target_properties(kernel PROPERTIES CXX_VISIBILITY_PRESET default)
        set_target_properties(kernel PROPERTIES VISIBILITY_INLINES_HIDDEN FALSE)
    else ()
        set_target_properties(kernel PROPERTIES CXX_VISIBILITY_PRESET hidden)
        set_target_properties(kernel PROPERTIES VISIBILITY_INLINES_HIDDEN TRUE)
    endif ()

# Install
#################################################################
# This is a preliminary solution!
install(
    TARGETS kernel
    RUNTIME DESTINATION ${TUDATPY_INSTALL_PATH}
    LIBRARY DESTINATION ${TUDATPY_INSTALL_PATH}
)
install(
    DIRECTORY ${TUDATPY_SOURCE_DIR}
    DESTINATION ${YACMA_PYTHON_MODULES_INSTALL_PATH}
)
