# Update version and build capabilities
# This should be replaced with setuptools_scm based versioning
if(TUDATPY_BUILD_WITH_QUAD_PRECISION_EXPOSURE)
    set(TUDATPY_QUAD_PRECISION_AVAILABLE True)
else()
    set(TUDATPY_QUAD_PRECISION_AVAILABLE False)
endif()
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
        ${TUDATPY_SOURCE_DIR}/pythonStateScalar.cpp
        ${TUDATPY_SOURCE_DIR}/expose_state_scalar_bindings.cpp
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

set(TUDATPY_STATE_SCALAR_BINDING_SOURCES
        dynamics/environment/expose_environment_state_scalar.cpp
        dynamics/environment_setup/ephemeris/expose_ephemeris_state_scalar.cpp
        dynamics/environment_setup/expose_environment_setup_state_scalar.cpp
        dynamics/parameters/expose_parameters_state_scalar.cpp
        dynamics/parameters_setup/expose_parameters_setup_state_scalar.cpp
        dynamics/propagation/expose_propagation_results_state_scalar.cpp
        dynamics/propagation/expose_propagation_state_utilities_state_scalar.cpp
        dynamics/propagation_setup/propagator/expose_propagator_state_scalar.cpp
        dynamics/simulator/expose_simulator_dynamics_state_scalar.cpp
        dynamics/simulator/expose_simulator_variational_state_scalar.cpp
        estimation/estimation_analysis/expose_estimation_analysis_state_scalar.cpp
        estimation/estimation_analysis/expose_estimation_analysis_ephemeris_fit_state_scalar.cpp
        estimation/estimation_analysis/expose_estimation_analysis_estimator_state_scalar.cpp
        estimation/observable_models/observables_simulation/expose_observables_simulation_state_scalar.cpp
        estimation/observations/expose_observations_state_scalar.cpp
        estimation/observations_setup/observations_simulation_settings/expose_observations_simulation_settings_factory_state_scalar.cpp
        estimation/observations_setup/observations_wrapper/expose_observations_wrapper_io_state_scalar.cpp
        estimation/observations_setup/observations_wrapper/expose_observations_wrapper_simulation_state_scalar.cpp
        math/interpolators/expose_interpolators_state_scalar.cpp
)

list(REMOVE_ITEM sources ${TUDATPY_STATE_SCALAR_BINDING_SOURCES})

YACMA_PYTHON_MODULE(kernel ${sources})

target_compile_definitions(kernel PRIVATE
    TUDATPY_BUILD_WITH_QUAD_PRECISION_EXPOSURE=$<BOOL:${TUDATPY_BUILD_WITH_QUAD_PRECISION_EXPOSURE}>
)

function(add_tudatpy_state_scalar_bindings target_name scalar_type binding_namespace supports_serialization)
    add_library(${target_name} OBJECT
        ${TUDATPY_STATE_SCALAR_BINDING_SOURCES}
    )
    set_target_properties(${target_name} PROPERTIES
        POSITION_INDEPENDENT_CODE ON
        UNITY_BUILD OFF
    )
    target_compile_definitions(${target_name} PRIVATE
        TUDATPY_BUILD_WITH_QUAD_PRECISION_EXPOSURE=$<BOOL:${TUDATPY_BUILD_WITH_QUAD_PRECISION_EXPOSURE}>
        TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS=$<BOOL:${TUDATPY_BUILD_WITH_DETAILED_PYBIND11_ERRORS}>
        TUDATPY_STATE_SCALAR_BINDING_NAMESPACE=${binding_namespace}
        TUDATPY_STATE_SCALAR_SUPPORTS_SERIALIZATION=${supports_serialization}
        STATE_SCALAR_TYPE=${scalar_type}
        "${pybind11_DEFINITIONS}"
    )
    target_include_directories(${target_name} PRIVATE
        ${Tudat_INCLUDE_DIRS}
        ${Boost_INCLUDE_DIRS}
        ${EIGEN3_INCLUDE_DIRS}
        ${CSpice_INCLUDE_DIRS}
        ${Sofa_INCLUDE_DIRS}
    )
    target_include_directories(${target_name} SYSTEM PRIVATE
        "${pybind11_INCLUDE_DIR}"
        "${YACMA_PYTHON_INCLUDE_DIR}"
        "${Python_INCLUDE_DIRS}"
        "${Boost_INCLUDE_DIRS}"
        "${EIGEN3_INCLUDE_DIRS}"
        "${CSpice_INCLUDE_DIRS}"
        "${Sofa_INCLUDE_DIRS}"
        "${Tudat_INCLUDE_DIRS}"
    )
    target_link_libraries(${target_name} PRIVATE
        ${Boost_LIBRARIES}
        ${Boost_SYSTEM_LIBRARY}
        ${Tudat_PROPAGATION_LIBRARIES}
        ${Tudat_ESTIMATION_LIBRARIES}
    )
endfunction()

add_tudatpy_state_scalar_bindings(
    kernel_double_state_scalar_bindings
    double
    double_precision
    1
)
target_sources(kernel PRIVATE
    $<TARGET_OBJECTS:kernel_double_state_scalar_bindings>
)

if(TUDATPY_BUILD_WITH_QUAD_PRECISION_EXPOSURE)
    add_tudatpy_state_scalar_bindings(
        kernel_quad_state_scalar_bindings
        tudat::HighPrecisionStateScalar
        quad_precision
        0
    )
    target_sources(kernel PRIVATE
        $<TARGET_OBJECTS:kernel_quad_state_scalar_bindings>
    )
endif()

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
    foreach(state_scalar_target
            kernel_double_state_scalar_bindings
            kernel_quad_state_scalar_bindings)
        if(TARGET ${state_scalar_target})
            target_include_directories(${state_scalar_target} SYSTEM PRIVATE
            "${NetCDF_INCLUDE_DIRS}"
            "${CMAKE_SOURCE_DIR}/third_parties/mcd"
            )
        endif()
    endforeach()
endif()

# Link against libraries
target_link_libraries(kernel PRIVATE
    Boost::boost
    ${Tudat_PROPAGATION_LIBRARIES}
    ${Tudat_ESTIMATION_LIBRARIES}
)

# ---- Precompiled header for pybind11 bindings ----
# These headers are included by virtually every binding .cpp file.
# PCH eliminates redundant parsing across all 86+ translation units.
if (TUDAT_BUILD_WITH_PCH)
    # clang-cl can emit invalid exception copy metadata for types loaded from a
    # PCH, causing pybind11 exception translation to corrupt the Windows heap.
    # See https://github.com/llvm/llvm-project/issues/118276.
    if (WIN32 AND CMAKE_CXX_COMPILER_ID STREQUAL "Clang"
            AND CMAKE_CXX_SIMULATE_ID STREQUAL "MSVC")
        message(STATUS "Disabling TudatPy binding PCH for clang-cl on Windows")
    else ()
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
        foreach(state_scalar_target
                kernel_double_state_scalar_bindings
                kernel_quad_state_scalar_bindings)
            if(TARGET ${state_scalar_target})
                target_precompile_headers(${state_scalar_target} PRIVATE
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
            endif()
        endforeach()
    endif ()
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
    foreach(tudatpy_binding_target
            kernel
            kernel_double_state_scalar_bindings
            kernel_quad_state_scalar_bindings)
        if(TARGET ${tudatpy_binding_target})
            if (APPLE AND TUDAT_BUILD_STATIC_LIBRARY AND TUDAT_BUILD_WITH_SERIALIZATION)
                set_target_properties(${tudatpy_binding_target} PROPERTIES
                    CXX_VISIBILITY_PRESET default
                    VISIBILITY_INLINES_HIDDEN FALSE)
            else ()
                set_target_properties(${tudatpy_binding_target} PROPERTIES
                    CXX_VISIBILITY_PRESET hidden
                    VISIBILITY_INLINES_HIDDEN TRUE)
            endif ()
        endif()
    endforeach()

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
