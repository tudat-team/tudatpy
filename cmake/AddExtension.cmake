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

    # Generate stubs for extension
    add_custom_command(
        TARGET "expose_${extension_name}" POST_BUILD
        COMMAND pybind11-stubgen "${import_path}.expose_${extension_name}" -o . --root-suffix=-stubs --numpy-array-remove-parameters --print-invalid-expressions-as-is
        WORKING_DIRECTORY "${TUDATPY_SOURCE_DIR}/src"
        COMMENT "Generating stubs for ${import_path}..."
    )

    # # Generate __init__.pyi for extension
    # if (EXISTS "${TUDATPY_SOURCE_DIR}/src/tudatpy-stubs")
    # else()
    #     file(MAKE_DIRECTORY "${TUDATPY_SOURCE_DIR}/src/tudatpy-stubs")
    # endif()
    # list(REMOVE_AT import_path_list 0)
    # list(JOIN import_path_list "/" stub_path)
    # set(stub_path "${TUDATPY_SOURCE_DIR}/src/tudatpy-stubs/${stub_path}")
    # configure_file(
    #     "${TUDATPY_SOURCE_DIR}/src/${extension_path}/__init__.py"
    #     "${stub_path}/__init__.pyi"
    #     COPYONLY
    # )

endmacro()
