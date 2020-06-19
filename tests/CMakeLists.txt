
include_directories(include)
add_subdirectory(src/astro)
add_subdirectory(src/basics)
add_subdirectory(src/interface)
add_subdirectory(src/math)
add_subdirectory(src/simulation)
add_subdirectory(src/io)

#if(TUDAT_INSTALL_TESTS)
install(DIRECTORY data/ DESTINATION "${INSTALL_BIN_DIR}/tudat/tests/data")
file(COPY data/ DESTINATION ${PROJECT_BINARY_DIR}/tests/data/)

#endif()