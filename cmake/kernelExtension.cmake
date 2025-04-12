macro (add_module extension import_path)

    # Get module name and path from import path
    string(REPLACE "." ";" extension_path_list ${import_path})
    list(GET extension_path_list -1 extension_name)
    set(extension_name expose_${extension_name})
    list(JOIN extension_path_list "/" extension_path)
    set(source_dir ${TUDATPY_SOURCE_DIR}/${extension_path})

    # Find all source files for extension [All .cpp files in source_dir]
    file(
        GLOB_RECURSE sources
        RELATIVE "${TUDATPY_SOURCE_DIR}"
        "${source_dir}/*.cpp"
    )
    # if (${ARGN})
    #     file(GLOB_RECURSE sources RELATIVE "${TUDATPY_SOURCE_DIR}" "${source_dir}/*.cpp" ${ARGN})
    #     message(STATUS "SOURCES FOR ${extension_name}: ${sources}")
    # else()
    #     file(GLOB sources RELATIVE "${TUDATPY_SOURCE_DIR}" "${source_dir}/*.cpp")
    # endif()
    # # string(REPLACE ";" " " sources "${sources}")

    # Update output directory for extension
    # set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${source_dir})

    # Create YACMA extension
    YACMA_PYTHON_MODULE(${extension} ${sources})

    # message(STATUS "sources: ${sources}")

    # message(STATUS "TUDATPY_SOURCE_DIR: ${TUDATPY_SOURCE_DIR}")
    # message(STATUS "extension_path: ${extension_path}")
    # message(STATUS "")
    # message(FATAL_ERROR "Import path: ${import_path}")


    # set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${TUDATPY_SOURCE_DIR}/${extension_path})

    # # Add extension
    # YACMA_PYTHON_MODULE(
    #     ${extension_name}
    #     ${extension_path}/${extension_name}.cpp
    # )

    list(APPEND tudat_libraries ${Tudat_PROPAGATION_LIBRARIES})
    list(APPEND tudat_libraries ${Tudat_ESTIMATION_LIBRARIES})
    list(REMOVE_DUPLICATES tudat_libraries)

    target_link_libraries(${extension} PRIVATE
        ${Boost_LIBRARIES}
        ${Boost_SYSTEM_LIBRARY}
        ${tudat_libraries}
    )

    target_include_directories(${extension} PUBLIC
        $<BUILD_INTERFACE:${Boost_INCLUDE_DIRS}>
        $<BUILD_INTERFACE:${Tudat_INCLUDE_DIRS}>
        $<BUILD_INTERFACE:${EIGEN3_INCLUDE_DIRS}>
        $<BUILD_INTERFACE:${CSpice_INCLUDE_DIRS}>
        $<BUILD_INTERFACE:${Sofa_INCLUDE_DIRS}>
        $<INSTALL_INTERFACE:include>)

    target_include_directories(${extension} SYSTEM PRIVATE "${pybind11_INCLUDE_DIR}")
    target_include_directories(${extension} SYSTEM PRIVATE "${EIGEN3_INCLUDE_DIRS}")
    target_include_directories(${extension} SYSTEM PRIVATE "${CSpice_INCLUDE_DIRS}")
    target_include_directories(${extension} SYSTEM PRIVATE "${Sofa_INCLUDE_DIRS}")
    target_include_directories(${extension} SYSTEM PRIVATE "${Tudat_INCLUDE_DIRS}")
    target_compile_definitions(${extension} PRIVATE "${pybind11_DEFINITIONS}")
    set_target_properties(${extension} PROPERTIES CXX_VISIBILITY_PRESET hidden)
    set_target_properties(${extension} PROPERTIES VISIBILITY_INLINES_HIDDEN TRUE)

    # unset(CMAKE_LIBRARY_OUTPUT_DIRECTORY)

    # Update kernel target with extension
    # add_dependencies(kernel ${extension})

    # Install
    install(
        TARGETS ${extension}
        RUNTIME DESTINATION ${TUDATPY_INSTALL_PATH}/${extension_path}
        LIBRARY DESTINATION ${TUDATPY_INSTALL_PATH}/${extension_path}
    )

endmacro()
