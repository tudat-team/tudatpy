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

    # Link
    target_link_libraries("expose_${extension_name}" PRIVATE
        ${TUDATPY_INTERNAL_LIBRARIES}
        ${TUDATPY_EXTERNAL_LIBRARIES}
    )

    # Restore installation prefix
    unset(CMAKE_LIBRARY_OUTPUT_DIRECTORY)

    # Generate stubs
    add_custom_command(
        TARGET "expose_${extension_name}"
        POST_BUILD
        COMMAND stubgen -p ${import_path} -o . > /dev/null 2>&1
        WORKING_DIRECTORY "${TUDATPY_SOURCE_DIR}/src"
        COMMENT "Generating stubs for ${extension_name}..."
    )

endmacro()
