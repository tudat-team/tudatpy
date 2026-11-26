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
        ${TUDATPY_SOURCE_DIR}/add_data_input_to_kernel.cpp
        ${TUDATPY_SOURCE_DIR}/add_util_to_kernel.cpp
        ${TUDATPY_SOURCE_DIR}/add_dynamics_to_kernel.cpp
        ${TUDATPY_SOURCE_DIR}/add_estimation_to_kernel.cpp
        ${TUDATPY_SOURCE_DIR}/add_exceptions_to_kernel.cpp
)
update_sources(math sources)
update_sources(astro sources)
update_sources(trajectory_design sources)
update_sources(constants sources)
update_sources(interface sources)
update_sources(data_input sources)
update_sources(util sources)
update_sources(numerical_simulation sources)
update_sources(dynamics sources)
update_sources(estimation sources)
update_sources(exceptions sources)

# Explicit instantiations of the heavy <STATE_SCALAR_TYPE, TIME_TYPE> propagation/estimation
# templates, each compiled once; the other kernel sources declare them extern via the per-area headers
# include/tudatpy/kernelExternTemplates*.h. Each group below is compiled as one unity source; the groups
# are chosen so that their members need the same Tudat headers (orbit determination manager and
# estimated parameter creation on their own).
set(kernel_instantiation_groups simulators propagation_functions variational observations estimation parameters orbit_determination)
set(kernel_instantiation_simulators
        ${TUDATPY_SOURCE_DIR}/kernelInstantiationDynamicsSimulator.cpp
        ${TUDATPY_SOURCE_DIR}/kernelInstantiationSingleArcDynamicsSimulator.cpp
        ${TUDATPY_SOURCE_DIR}/kernelInstantiationMultiArcDynamicsSimulator.cpp
        ${TUDATPY_SOURCE_DIR}/kernelInstantiationHybridArcDynamicsSimulator.cpp)
set(kernel_instantiation_propagation_functions
        ${TUDATPY_SOURCE_DIR}/kernelInstantiationPropagationFunctions.cpp
        ${TUDATPY_SOURCE_DIR}/kernelInstantiationStateDerivativeModels.cpp)
set(kernel_instantiation_variational
        ${TUDATPY_SOURCE_DIR}/kernelInstantiationVariationalEquationsSolver.cpp
        ${TUDATPY_SOURCE_DIR}/kernelInstantiationSingleArcVariationalEquationsSolver.cpp
        ${TUDATPY_SOURCE_DIR}/kernelInstantiationMultiArcVariationalEquationsSolver.cpp
        ${TUDATPY_SOURCE_DIR}/kernelInstantiationHybridArcVariationalEquationsSolver.cpp
        ${TUDATPY_SOURCE_DIR}/kernelInstantiationStateDerivativePartials.cpp)
set(kernel_instantiation_observations
        ${TUDATPY_SOURCE_DIR}/kernelInstantiationObservationManager1.cpp
        ${TUDATPY_SOURCE_DIR}/kernelInstantiationObservationManager2.cpp
        ${TUDATPY_SOURCE_DIR}/kernelInstantiationObservationManager3.cpp
        ${TUDATPY_SOURCE_DIR}/kernelInstantiationObservationSimulators.cpp)
set(kernel_instantiation_estimation
        ${TUDATPY_SOURCE_DIR}/kernelInstantiationEstimationFunctions.cpp)
set(kernel_instantiation_parameters
        ${TUDATPY_SOURCE_DIR}/kernelInstantiationParametersToEstimate.cpp)
set(kernel_instantiation_orbit_determination
        ${TUDATPY_SOURCE_DIR}/kernelInstantiationOrbitDeterminationManager.cpp)
foreach (group ${kernel_instantiation_groups})
    list(APPEND sources ${kernel_instantiation_${group}})
endforeach ()

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
        message(STATUS "Disabling TudatPy kernel PCH for clang-cl on Windows")
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
    endif ()
endif ()

# These binding files only use the extern <STATE_SCALAR_TYPE, TIME_TYPE> instantiations; compiling
# each list as one unity source parses its propagation/estimation headers once instead of once per file.
# The orbit determination manager bindings form their own group, since their headers are the heaviest.
# Other binding sources are not merged.
set(kernel_shared_propagation_sources
    ${TUDATPY_SOURCE_DIR}/dynamics/parameters_setup/expose_parameters_setup.cpp
    ${TUDATPY_SOURCE_DIR}/dynamics/propagation/expose_propagation_state_utilities.cpp
    ${TUDATPY_SOURCE_DIR}/dynamics/simulator/expose_simulator_dynamics.cpp
    ${TUDATPY_SOURCE_DIR}/dynamics/simulator/expose_simulator_variational.cpp
    ${TUDATPY_SOURCE_DIR}/estimation/observations/expose_observations_simulation.cpp
    ${TUDATPY_SOURCE_DIR}/estimation/observations_setup/observations_simulation_settings/expose_observations_simulation_settings_factory.cpp
)
set(kernel_shared_estimation_sources
    ${TUDATPY_SOURCE_DIR}/estimation/estimation_analysis/expose_estimation_analysis.cpp
    ${TUDATPY_SOURCE_DIR}/estimation/estimation_analysis/expose_estimation_analysis_estimator.cpp
    ${TUDATPY_SOURCE_DIR}/estimation/estimation_analysis/expose_estimation_analysis_ephemeris_fit.cpp
)
if (TUDAT_BUILD_WITH_UNITY_GROUPS)
    set_target_properties(kernel PROPERTIES UNITY_BUILD ON UNITY_BUILD_MODE GROUP)
    set_source_files_properties(${kernel_shared_propagation_sources}
        PROPERTIES UNITY_GROUP shared_propagation)
    set_source_files_properties(${kernel_shared_estimation_sources}
        PROPERTIES UNITY_GROUP shared_estimation)
    foreach (group ${kernel_instantiation_groups})
        set_source_files_properties(${kernel_instantiation_${group}}
            PROPERTIES UNITY_GROUP instantiation_${group})
    endforeach ()
endif ()

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

# On macOS, debug info stays in the object files, bundle it into a .dSYM so profilers can symbolize after the objects are cleaned.
if (APPLE AND TUDAT_BUILD_WITH_PROFILING_SYMBOLS)
    find_program(DSYMUTIL_EXECUTABLE dsymutil)
    if (DSYMUTIL_EXECUTABLE)
        add_custom_command(TARGET kernel POST_BUILD
            COMMAND ${DSYMUTIL_EXECUTABLE} $<TARGET_FILE:kernel>
            COMMENT "Generating dSYM for kernel"
            VERBATIM
        )
    else ()
        message(WARNING "dsymutil not found: kernel.dSYM will not be generated.")
    endif ()
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
