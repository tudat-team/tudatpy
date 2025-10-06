macro (add_extension import_path)

    # Get module name and path from import path
    string(REPLACE "." ";" extension_path_list ${import_path})
    list(GET extension_path_list -1 extension_name)
    set(extension_name expose_${extension_name})
    list(JOIN extension_path_list "/" extension_path)
    set(source_dir ${TUDATPY_SOURCE_DIR}/${extension_path})

    # Find all source files for extension [All .cpp files in source_dir]
    if (${ARGN})
        file(GLOB_RECURSE sources RELATIVE "${TUDATPY_SOURCE_DIR}" "${source_dir}/*.cpp" ${ARGN})
        message(STATUS "SOURCES FOR ${extension_name}: ${sources}")
    else()
        file(GLOB sources RELATIVE "${TUDATPY_SOURCE_DIR}" "${source_dir}/*.cpp")
    endif()
    # string(REPLACE ";" " " sources "${sources}")

    # Update output directory for extension
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${source_dir})

    # Create YACMA extension
    YACMA_PYTHON_MODULE(${extension_name} ${sources})

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

    target_link_libraries(${extension_name} PRIVATE
        ${Boost_LIBRARIES}
        ${Boost_SYSTEM_LIBRARY}
        ${Tudat_PROPAGATION_LIBRARIES}
        ${Tudat_ESTIMATION_LIBRARIES}
    )

    target_include_directories(${extension_name} PUBLIC
        $<BUILD_INTERFACE:${Boost_INCLUDE_DIRS}>
        $<BUILD_INTERFACE:${Tudat_INCLUDE_DIRS}>
        $<BUILD_INTERFACE:${EIGEN3_INCLUDE_DIRS}>
        $<BUILD_INTERFACE:${CSpice_INCLUDE_DIRS}>
        $<BUILD_INTERFACE:${Sofa_INCLUDE_DIRS}>
        $<INSTALL_INTERFACE:include>)

    target_include_directories(${extension_name} SYSTEM PRIVATE "${pybind11_INCLUDE_DIR}")
    target_include_directories(${extension_name} SYSTEM PRIVATE "${EIGEN3_INCLUDE_DIRS}")
    target_include_directories(${extension_name} SYSTEM PRIVATE "${CSpice_INCLUDE_DIRS}")
    target_include_directories(${extension_name} SYSTEM PRIVATE "${Sofa_INCLUDE_DIRS}")
    target_include_directories(${extension_name} SYSTEM PRIVATE "${Tudat_INCLUDE_DIRS}")
    target_compile_definitions(${extension_name} PRIVATE "${pybind11_DEFINITIONS}")
    set_target_properties(${extension_name} PROPERTIES CXX_VISIBILITY_PRESET hidden)
    set_target_properties(${extension_name} PROPERTIES VISIBILITY_INLINES_HIDDEN TRUE)

    unset(CMAKE_LIBRARY_OUTPUT_DIRECTORY)

    # Update kernel target with extension
    add_dependencies(kernel ${extension_name})

    # Install
    install(
        TARGETS ${extension_name}
        RUNTIME DESTINATION ${TUDATPY_INSTALL_PATH}/${extension_path}
        LIBRARY DESTINATION ${TUDATPY_INSTALL_PATH}/${extension_path}
    )

endmacro()

macro (update_sources import_path sources)

    # Get module name and path from import path
    string(REPLACE "." ";" extension_path_list ${import_path})
    list(GET extension_path_list -1 extension_name)
    set(extension_name expose_${extension_name})
    list(JOIN extension_path_list "/" extension_path)
    set(source_dir ${TUDATPY_SOURCE_DIR}/${extension_path})

    # Find all source files for extension [All .cpp files in source_dir]
    file(
        GLOB_RECURSE new_sources
        RELATIVE "${TUDATPY_SOURCE_DIR}"
        "${source_dir}/*.cpp"
    )

    list(APPEND sources ${new_sources})

endmacro()

macro (copy_python_in_build)

    # Copy python files to build directory
    file(GLOB_RECURSE py_files "${TUDATPY_SOURCE_DIR}/**/*.py")
    foreach(py_file ${py_files})
        file(RELATIVE_PATH py_file_name ${TUDATPY_SOURCE_DIR} ${py_file})
        get_filename_component(parents ${py_file_name} DIRECTORY)
        file(COPY ${py_file_name} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/${parents})
    endforeach()

    # Copy __init__.py in base directory
    file(
        COPY ${TUDATPY_SOURCE_DIR}/__init__.py
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
    )

endmacro()
