include(CMakeParseArguments)

function(TUDAT_ADD_LIBRARY arg1 arg2 arg3)
    # arg1: library name
    # arg3: sources
    # arg4: headers
    # Design function parser.
    cmake_parse_arguments(
            PARSED_ARGS
            ""
            ""
            "PUBLIC_LINKS;PRIVATE_LINKS;INTERFACE_LINKS;PRIVATE_INCLUDES"
            ${ARGN})

    set(target_name "tudat_${arg1}")
    # Setup of library.
    if (TUDAT_BUILD_STATIC_LIBRARY)
        add_library(${target_name} STATIC "${arg2}")
    else ()
        add_library(${target_name} SHARED "${arg2}")
    endif ()
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
            SYSTEM PRIVATE ${PARSED_ARGS_PRIVATE_INCLUDES} ${Boost_INCLUDE_DIRS} ${CSpice_INCLUDE_DIRS} ${Sofa_INCLUDE_DIRS} ${TudatResources_INCLUDE_DIRS}
            )

    target_link_libraries("${target_name}"
            PUBLIC    ${PARSED_ARGS_PUBLIC_LINKS}
            PRIVATE   ${PARSED_ARGS_PRIVATE_LINKS}
            INTERFACE ${PARSED_ARGS_INTERFACE_LINKS}
            )
    #==========================================================================
    # BUILD-TREE.
    #==========================================================================
    set_target_properties("${target_name}"
            PROPERTIES
            LINKER_LANGUAGE CXX
            ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
            LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
            )

    add_library(${PROJECT_NAME}::${target_name} ALIAS "${target_name}")
    #==========================================================================
    # INSTALL-TREE.
    #==========================================================================
    install(TARGETS              "${target_name}"
            EXPORT               tudat_export
            LIBRARY DESTINATION  "${INSTALL_LIB_DIR}"
            ARCHIVE DESTINATION  "${INSTALL_LIB_DIR}"
            RUNTIME DESTINATION  "${INSTALL_BIN_DIR}"
            )
    unset(target_name)
endfunction()

function(TUDAT_ADD_TEST_LIBRARY arg1 arg2 arg3)
    # arg1: library name
    # arg3: sources
    # arg4: headers
    # Design function parser.
    cmake_parse_arguments(
            PARSED_ARGS
            ""
            ""
            "PUBLIC_LINKS;PRIVATE_LINKS;INTERFACE_LINKS;PRIVATE_INCLUDES"
            ${ARGN})

    set(target_name "tudat_test_${arg1}")
    # Setup of library.
    if (TUDAT_BUILD_STATIC_LIBRARY)
        add_library(${target_name} STATIC "${arg2}")
    else ()
        add_library(${target_name} SHARED "${arg2}")
    endif ()
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
            $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/tests/include>
            $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/tests/include>
            $<INSTALL_INTERFACE:include>)

    target_include_directories("${target_name}"
            SYSTEM PRIVATE ${PARSED_ARGS_PRIVATE_INCLUDES} ${Boost_INCLUDE_DIRS} ${CSpice_INCLUDE_DIRS} ${Sofa_INCLUDE_DIRS} ${TudatResources_INCLUDE_DIRS}
            )

    target_link_libraries("${target_name}"
            PUBLIC    ${PARSED_ARGS_PUBLIC_LINKS}
            PRIVATE   ${PARSED_ARGS_PRIVATE_LINKS}
            INTERFACE ${PARSED_ARGS_INTERFACE_LINKS}
            )
    #==========================================================================
    # BUILD-TREE.
    #==========================================================================
    set_target_properties("${target_name}"
            PROPERTIES
            LINKER_LANGUAGE CXX
            ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
            LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
            )

    add_library(${PROJECT_NAME}::${target_name} ALIAS "${target_name}")
    #==========================================================================
    # INSTALL-TREE.
    #==========================================================================
    install(TARGETS              "${target_name}"
            EXPORT               tudat_export
            LIBRARY DESTINATION  "${INSTALL_LIB_DIR}"
            ARCHIVE DESTINATION  "${INSTALL_LIB_DIR}"
            RUNTIME DESTINATION  "${INSTALL_BIN_DIR}"
            )
    unset(target_name)
endfunction()
