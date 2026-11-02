include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
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
