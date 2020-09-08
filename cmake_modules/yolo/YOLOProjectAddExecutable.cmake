function(TUDAT_ADD_EXECUTABLE arg1 arg2)
    # arg1: executable name
    # arg2: sources

    # Create target name.
    set(target_name "${arg1}")

    # Add executable.
    add_executable(${target_name} ${arg2})

    #==========================================================================
    # TARGET-CONFIGURATION.
    #==========================================================================
    # NOTE: make sure the include directories from the current build
    # are included first, so that if there is already a pagmo installation
    # in the prefix path we don't risk including the headers from that
    # one instead.
    target_include_directories("${target_name}" PUBLIC
            $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
            $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
            $<INSTALL_INTERFACE:include>)

    target_include_directories("${target_name}"
            SYSTEM PRIVATE "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}" "${CSpice_INCLUDE_DIRS}" "${Sofa_INCLUDE_DIRS}" "${TudatResources_INCLUDE_DIRS}"
            )

    message(STATUS LIB:BOOST:${Boost_LIBRARIES})
    target_link_libraries("${target_name}"
            PUBLIC    ${ARGN}
            PRIVATE   "${Boost_LIBRARIES}"
            )

    #==========================================================================
    # BUILD-TREE.
    #==========================================================================
    set_target_properties(${target_name}
            PROPERTIES
            LINKER_LANGUAGE CXX
            ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib"
            LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib"
            RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin"
            )

    # Let's setup the target C++ standard, but only if the user did not provide it manually.
    if (NOT CMAKE_CXX_STANDARD)
        set_property(TARGET ${target_name} PROPERTY CXX_STANDARD 17)
    endif ()
    set_property(TARGET ${target_name} PROPERTY CXX_STANDARD_REQUIRED YES)
    set_property(TARGET ${target_name} PROPERTY CXX_EXTENSIONS NO)
#    add_test(${target_name} "${PROJECT_BINARY_DIR}/tests/${target_name}")

    #==========================================================================
    # INSTALL-TREE.
    #==========================================================================
    install(TARGETS "${target_name}"
            EXPORT tudat_export
            LIBRARY DESTINATION "${INSTALL_LIB_DIR}"
            ARCHIVE DESTINATION "${INSTALL_LIB_DIR}"
            RUNTIME DESTINATION "${INSTALL_BIN_DIR}/tudat"
            )

    # Clean up set variables.
    unset(target_name)

endfunction()