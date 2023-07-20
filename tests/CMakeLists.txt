add_subdirectory(src/astro)
#add_subdirectory(src/basics)
#add_subdirectory(src/interface)
#add_subdirectory(src/math)
#add_subdirectory(src/simulation)
#add_subdirectory(src/io)
#add_subdirectory(src/support)

# Unzip compressed data and place in build tree.
add_custom_target( unpack_test_data ALL)
add_custom_command(TARGET unpack_test_data PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E remove_directory ${PROJECT_BINARY_DIR}/tests/data/
        COMMAND ${CMAKE_COMMAND} -E tar xzf ${CMAKE_CURRENT_SOURCE_DIR}/data.zip
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/data.zip
        COMMENT "Unpacking data.zip"
        VERBATIM)
