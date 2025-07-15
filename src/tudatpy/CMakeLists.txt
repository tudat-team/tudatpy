# Update version
configure_file(
    "${TUDATPY_SOURCE_DIR}/../../__init__.py.in"
    "${TUDATPY_SOURCE_DIR}/__init__.py"
    @ONLY
)

# include(YacmaExtension)
include(PythonExtension)

copy_python_in_build()

# Setup the installation path.
set(TUDATPY_INSTALL_PATH "${YACMA_PYTHON_MODULES_INSTALL_PATH}/tudatpy")
# set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${TUDATPY_SOURCE_DIR})

# Get sources for extension
set(sources ${TUDATPY_SOURCE_DIR}/kernel.cpp)
update_sources(math sources)
update_sources(astro sources)
update_sources(trajectory_design sources)
update_sources(constants sources)
update_sources(interface sources)
update_sources(data sources)
update_sources(numerical_simulation sources)
update_sources(dynamics sources)
update_sources(estimation_refactoring sources)
update_sources(exceptions sources)

YACMA_PYTHON_MODULE(kernel ${sources})

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
