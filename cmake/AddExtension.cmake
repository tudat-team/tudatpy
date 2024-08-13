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

endmacro()


macro (generate_stubs import_path)

    string(REPLACE "." ";" import_path_list ${import_path})
    list(GET import_path_list -1 extension_name)

    # Generate stubs
    add_custom_target(${extension_name})
    add_custom_command(
        TARGET ${extension_name} POST_BUILD
        COMMAND pybind11-stubgen ${import_path} -o . --root-suffix=-stubs
        WORKING_DIRECTORY "${TUDATPY_SOURCE_DIR}/src"
        COMMENT "Generating stubs for ${import_path}..."
    )
endmacro()
