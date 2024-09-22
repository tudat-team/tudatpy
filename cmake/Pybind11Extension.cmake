macro (add_extension import_path)

    # Get module name and path from import path
    string(REPLACE "." ";" extension_path_list ${import_path})
    list(GET extension_path_list -1 extension_name)
    set(extension_name expose_${extension_name})
    list(SUBLIST extension_path_list 1 -1 extension_path_list)
    list(JOIN extension_path_list "/" extension_path)

    # Add extension
    pybind11_add_module(
        ${extension_name} MODULE
        ${TUDATPY_SOURCE_DIR}/${extension_path}/${extension_name}.cpp
    )

    # Link against libraries
    target_link_directories(${extension_name} PRIVATE "${Tudat_INCLUDE_DIRS}/../lib")
    target_link_libraries(
        ${extension_name} PRIVATE
        ${TUDATPY_INTERNAL_LIBRARIES}
        ${TUDATPY_EXTERNAL_LIBRARIES}
    )

    # Headers
    target_include_directories(${extension_name} PUBLIC
        $<BUILD_INTERFACE:${Boost_INCLUDE_DIRS}>
        $<BUILD_INTERFACE:${Tudat_INCLUDE_DIRS}>
        $<BUILD_INTERFACE:${EIGEN3_INCLUDE_DIRS}>
        $<BUILD_INTERFACE:${CSpice_INCLUDE_DIRS}>
        $<BUILD_INTERFACE:${Sofa_INCLUDE_DIRS}>
        $<INSTALL_INTERFACE:include>)

    target_include_directories(
        ${extension_name} SYSTEM PRIVATE ${pybind11_INCLUDE_DIR}
    )
    target_include_directories(
        ${extension_name} SYSTEM PRIVATE ${EIGEN3_INCLUDE_DIRS})
    target_include_directories(
        ${extension_name} SYSTEM PRIVATE ${CSpice_INCLUDE_DIRS}
    )
    target_include_directories(
        ${extension_name} SYSTEM PRIVATE ${Sofa_INCLUDE_DIRS}
    )
    target_include_directories(
        ${extension_name} SYSTEM PRIVATE ${Tudat_INCLUDE_DIRS}
    )
    target_compile_definitions(
        ${extension_name} PRIVATE ${pybind11_DEFINITIONS}
    )
    set_target_properties(
        ${extension_name} PROPERTIES CXX_VISIBILITY_PRESET hidden
    )
    set_target_properties(
        ${extension_name} PROPERTIES VISIBILITY_INLINES_HIDDEN TRUE
    )

    # Add as dependency for kernel target
    add_dependencies(kernel ${extension_name})

    # Install extension
    install(
        TARGETS ${extension_name}
        DESTINATION tudatpy/${extension_path}
    )

endmacro()
