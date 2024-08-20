macro (add_extension import_path)

    # Get source from import path
    string(REPLACE "." ";" import_path_list ${import_path})
    list(GET import_path_list -1 extension_name)
    string(REPLACE "." "/" extension_path ${import_path})
    set(
        src_path
        "${TUDATPY_SOURCE_DIR}/src/${extension_path}/expose_${extension_name}.cpp"
    )

    # Update installation prefix
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${TUDATPY_SOURCE_DIR}/src/${extension_path}")

    # Add extension
    pybind11_add_module("expose_${extension_name}" MODULE ${src_path})

    # Headers
    target_include_directories("expose_${extension_name}" PUBLIC
        $<BUILD_INTERFACE:${Boost_INCLUDE_DIRS}>
        $<BUILD_INTERFACE:${Tudat_INCLUDE_DIRS}>
        $<BUILD_INTERFACE:${EIGEN3_INCLUDE_DIRS}>
        $<BUILD_INTERFACE:${CSpice_INCLUDE_DIRS}>
        $<BUILD_INTERFACE:${Sofa_INCLUDE_DIRS}>
        $<INSTALL_INTERFACE:include>)

    target_include_directories(
        "expose_${extension_name}" SYSTEM PRIVATE "${pybind11_INCLUDE_DIR}"
    )
    target_include_directories(
        "expose_${extension_name}" SYSTEM PRIVATE "${EIGEN3_INCLUDE_DIRS}")
    target_include_directories(
        "expose_${extension_name}" SYSTEM PRIVATE "${CSpice_INCLUDE_DIRS}"
    )
    target_include_directories(
        "expose_${extension_name}" SYSTEM PRIVATE "${Sofa_INCLUDE_DIRS}"
    )
    target_include_directories(
        "expose_${extension_name}" SYSTEM PRIVATE "${Tudat_INCLUDE_DIRS}"
    )
    target_compile_definitions(
        "expose_${extension_name}" PRIVATE "${pybind11_DEFINITIONS}"
    )
    set_target_properties(
        "expose_${extension_name}" PROPERTIES CXX_VISIBILITY_PRESET hidden
    )
    set_target_properties(
        "expose_${extension_name}" PROPERTIES VISIBILITY_INLINES_HIDDEN TRUE
    )

    # Link
    target_link_libraries("expose_${extension_name}" PRIVATE
        ${TUDATPY_INTERNAL_LIBRARIES}
        ${TUDATPY_EXTERNAL_LIBRARIES}
    )

    # Restore installation prefix
    unset(CMAKE_LIBRARY_OUTPUT_DIRECTORY)

endmacro()
