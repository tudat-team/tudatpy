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
    target_include_directories("${target_name}"
            PUBLIC
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
            $<INSTALL_INTERFACE:include>
            )

    target_include_directories("${target_name}"
            SYSTEM PRIVATE ${PARSED_ARGS_PRIVATE_INCLUDES}
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


#include(CMakeParseArguments)
#
## Get upper project name.
#string(TOUPPER ${PROJECT_NAME} _YOLO_UPPER_PROJECT_NAME)
#string(TOLOWER ${PROJECT_NAME} _YOLO_LOWER_PROJECT_NAME)
#
#function("${_YOLO_UPPER_PROJECT_NAME}_ADD_LIBRARY" arg1 arg2 arg3)
#    # arg1: library name
#    # arg3: sources
#    # arg4: headers
#
#    # Design function parser.
#    # https://cmake.org/cmake/help/v3.0/module/CMakeParseArguments.html
#    set(options BUILD_STATIC)
#    set(oneValueArgs PREFIX EXPORT)
#    set(multiValueArgs LINK_ITEMS)
#    cmake_parse_arguments("_LIBRARY" "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
#
#    # TODO: Detect mode. Parser or Standard.
#
#    # Create target name.
#    set(target_name "${arg1}")
#    set(target_name_no_prefix "${arg1}")
#
#    # Prepend custom prefix if desired.
#    if (NOT "${_LIBRARY_PREFIX}" STREQUAL "")
#        string(PREPEND target_name "${_LIBRARY_PREFIX}_")
#    elseif (NOT "${_${_YOLO_UPPER_PROJECT_NAME}_LIBRARY_PREFIX}" STREQUAL "")
#        string(PREPEND target_name "${_${_YOLO_UPPER_PROJECT_NAME}_LIBRARY_PREFIX}_")
#    endif ()
#
#    # Check if parsed argument exists (takes precedence).
#    if (DEFINED ${_LIBRARY_BUILD_STATIC})
#        if (${_LIBRARY_BUILD_STATIC})
#            add_library("${target_name}" STATIC ${arg2})
#        else ()
#            add_library("${target_name}" SHARED ${arg2})
#        endif ()
#        # Check if global argument exists.
#    elseif (DEFINED ${_${_YOLO_UPPER_PROJECT_NAME}_LIBRARY_BUILD_STATIC})
#        if (${_${_YOLO_UPPER_PROJECT_NAME}_LIBRARY_BUILD_STATIC})
#            add_library("${target_name}" STATIC ${arg2})
#        else ()
#            add_library("${target_name}" SHARED ${arg2})
#        endif ()
#    else ()
#        # Shared by default, if nothing is defined.
#        add_library("${target_name}" SHARED ${arg2})
#    endif ()
#
#    macro(PARSE_LIBRARY_ARGS arg default)
#        if (DEFINED ${_LIBRARY_${arg}})
#            set("_${arg}" ${_LIBRARY_${arg}})
#        elseif (DEFINED "${_${_YOLO_UPPER_PROJECT_NAME}_LIBRARY_${arg}}")
#            set("_${arg}" ${_${_YOLO_UPPER_PROJECT_NAME}_LIBRARY_${arg}})
#        else ()
#            set("_${arg}" "${default}")
#        endif ()
#    endmacro()
#
#    # BUILD-TREE : ARCHIVE_OUTPUT_DIRECTORY
#
##    if (DEFINED ${_LIBRARY_ARCHIVE_OUTPUT_DIRECTORY})
##        set(_ARCHIVE_OUTPUT_DIRECTORY ${_LIBRARY_ARCHIVE_OUTPUT_DIRECTORY})
##    elseif (DEFINED ${_${_YOLO_UPPER_PROJECT_NAME}_LIBRARY_ARCHIVE_OUTPUT_DIRECTORY})
##        set(_ARCHIVE_OUTPUT_DIRECTORY ${_${_YOLO_UPPER_PROJECT_NAME}_LIBRARY_ARCHIVE_OUTPUT_DIRECTORY})
##    else ()
##        set(_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
##    endif ()
#
#    # BUILD-TREE : LIBRARY_OUTPUT_DIRECTORY
##    if (DEFINED ${_LIBRARY_LIBRARY_OUTPUT_DIRECTORY})
##        set(_LIBRARY_OUTPUT_DIRECTORY ${_LIBRARY_LIBRARY_OUTPUT_DIRECTORY})
##    elseif (DEFINED ${_${_YOLO_UPPER_PROJECT_NAME}_LIBRARY_LIBRARY_OUTPUT_DIRECTORY})
##        set(_LIBRARY_OUTPUT_DIRECTORY ${_${_YOLO_UPPER_PROJECT_NAME}_LIBRARY_LIBRARY_OUTPUT_DIRECTORY})
##    else ()
##        set(_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
##    endif ()
#
#    # BUILD-TREE : RUNTIME_OUTPUT_DIRECTORY
##    if (DEFINED ${_LIBRARY_RUNTIME_OUTPUT_DIRECTORY})
##        set(_RUNTIME_OUTPUT_DIRECTORY ${_LIBRARY_RUNTIME_OUTPUT_DIRECTORY})
##    elseif (DEFINED ${_${_YOLO_UPPER_PROJECT_NAME}_LIBRARY_RUNTIME_OUTPUT_DIRECTORY})
##        set(_RUNTIME_OUTPUT_DIRECTORY ${_${_YOLO_UPPER_PROJECT_NAME}_LIBRARY_RUNTIME_OUTPUT_DIRECTORY})
##    else ()
##        set(_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
##    endif ()
#
#    # BUILD-TREE : COMPILE_FLAGS
#
##    if (DEFINED ${_LIBRARY_COMPILE_FLAGS})
##        set(_COMPILE_FLAGS ${_LIBRARY_COMPILE_FLAGS})
##    elseif (DEFINED ${_${_YOLO_UPPER_PROJECT_NAME}_LIBRARY_COMPILE_FLAGS})
##        set(_COMPILE_FLAGS ${_${_YOLO_UPPER_PROJECT_NAME}_LIBRARY_COMPILE_FLAGS})
##    else ()
##        set(_COMPILE_FLAGS "${CMAKE_CXX_FLAGS}")
##    endif ()
#
#    # BUILD-TREE.
#    PARSE_LIBRARY_ARGS("ARCHIVE_OUTPUT_DIRECTORY"  "${CMAKE_BINARY_DIR}/lib")
#    PARSE_LIBRARY_ARGS("LIBRARY_OUTPUT_DIRECTORY"  "${CMAKE_BINARY_DIR}/lib")
#    PARSE_LIBRARY_ARGS("RUNTIME_OUTPUT_DIRECTORY"  "${CMAKE_BINARY_DIR}/bin")
#    PARSE_LIBRARY_ARGS("COMPILE_FLAGS" "")
#    set_target_properties("${target_name}"
#            PROPERTIES
#            LINKER_LANGUAGE CXX
#            ARCHIVE_OUTPUT_DIRECTORY "${_ARCHIVE_OUTPUT_DIRECTORY}"
#            LIBRARY_OUTPUT_DIRECTORY "${_LIBRARY_OUTPUT_DIRECTORY}"
#            RUNTIME_OUTPUT_DIRECTORY "${_RUNTIME_OUTPUT_DIRECTORY}"
##            COMPILE_FLAGS            "${CMAKE_CXX_FLAGS} ${_COMPILE_FLAGS}"
#            )
#
#    add_library(${PROJECT_NAME}::${target_name} ALIAS "${target_name}")
#
#    # Link libraries if extra arguments provided.
#    target_link_libraries("${target_name}" PUBLIC
#            "${ARGN}"
#            "${_${_YOLO_UPPER_PROJECT_NAME}_TEST_CASE_LINK_ITEMS}"
#            )
#
#            "${_${_YOLO_UPPER_PROJECT_NAME}_TEST_CASE_LINK_ITEMS}"
#    # INSTALL TREE.
#    PARSE_LIBRARY_ARGS("LIBRARY_DESTINATION"  "${CMAKE_INSTALL_PREFIX}/lib")
#    PARSE_LIBRARY_ARGS("ARCHIVE_DESTINATION"  "${CMAKE_INSTALL_PREFIX}/lib")
#    PARSE_LIBRARY_ARGS("RUNTIME_DESTINATION"  "${CMAKE_INSTALL_PREFIX}/bin")
##    PARSE_LIBRARY_ARGS("INCLUDES_DESTINATION" "${CMAKE_INSTALL_PREFIX}/include")
#    PARSE_LIBRARY_ARGS("EXPORT"               "${_YOLO_LOWER_PROJECT_NAME}_export")
#    install(TARGETS              "${target_name}"
#            EXPORT               "${_EXPORT}"
#            LIBRARY DESTINATION  "${_LIBRARY_DESTINATION}"
#            ARCHIVE DESTINATION  "${_ARCHIVE_DESTINATION}"
#            RUNTIME DESTINATION  "${_RUNTIME_DESTINATION}"
##            INCLUDES DESTINATION "${_INCLUDES_DESTINATION}"
##            COMPONENT Development
#            )
#
#
#
#    # Cleanup.
#    unset(_ARCHIVE_OUTPUT_DIRECTORY)
#    unset(_LIBRARY_OUTPUT_DIRECTORY)
#    unset(_RUNTIME_OUTPUT_DIRECTORY)
#    unset(_EXPORT)
#    unset(_LIBRARY_DESTINATION)
#    unset(_ARCHIVE_DESTINATION)
#    unset(_RUNTIME_DESTINATION)
#    unset(_INCLUDES_DESTINATION)
#    unset(_COMPILE_FLAGS)
#endfunction()
