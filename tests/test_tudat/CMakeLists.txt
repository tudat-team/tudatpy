include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
# Only the tests see the extern template declarations, so that switching the tests on or off does not recompile the
# libraries and the kernel; the libraries compile these instantiations in their *Instantiation.cpp sources
if (TUDAT_BUILD_TESTS_WITH_EXPLICIT_INSTANTIATIONS)
    add_compile_definitions(TUDAT_BUILD_TESTS_WITH_EXPLICIT_INSTANTIATIONS=1)
endif ()
add_subdirectory(src/astro)
add_subdirectory(src/basics)
add_subdirectory(src/interface)
add_subdirectory(src/math)
add_subdirectory(src/simulation)
add_subdirectory(src/io)
add_subdirectory(src/support)

# Copy test data to build tree (tests may also write output to this directory).
add_custom_target( copy_test_data ALL)
add_custom_command(
    TARGET copy_test_data PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${PROJECT_BINARY_DIR}/tests/data/
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/data ${PROJECT_BINARY_DIR}/tests/data
    COMMENT "Copying test data"
    VERBATIM
)
