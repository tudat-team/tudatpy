include(CMakeParseArguments)

function("TUDAT_ADD_TEST_CASE" arg1)
    if (TUDAT_BUILD_TESTS)

        # arg1 : Test name. Will add source file ${CMAKE_CURRENT_SOURCE_DIR}/tests/unitTest${arg1}.cpp
        # _${PROJECT_NAME}_TEST_CASE_ITEMS : Global dependencies to link to all.
        # ADD_DIRNAME : (bool) Adds the current dirname as prefix to test.
        # INSTALL : (bool) Install the test in the install-tree.
        cmake_parse_arguments(
                PARSED_ARGS
                ""
                ""
                "SOURCES;PRIVATE_LINKS"
                ${ARGN})

        # Create target name.
        get_filename_component(dirname ${CMAKE_CURRENT_SOURCE_DIR} NAME)
        set(target_name "test_${dirname}_${arg1}")

        # Add executable.
        add_executable(${target_name} ${CMAKE_CURRENT_SOURCE_DIR}/tests/unitTest${arg1}.cpp ${PARSED_ARGS_SOURCES})

        #==========================================================================
        # TARGET-CONFIGURATION.
        #==========================================================================
        target_include_directories("${target_name}"
                PUBLIC
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                $<INSTALL_INTERFACE:include>
                )

        target_include_directories("${target_name}"
                SYSTEM PRIVATE "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}" "${CSpice_INCLUDE_DIRS}" "${Sofa_INCLUDE_DIRS}"
                )

        target_link_libraries("${target_name}"
                PUBLIC ${PARSED_ARGS_PRIVATE_LINKS}
                PRIVATE "${Boost_LIBRARIES}"
                )

        #==========================================================================
        # BUILD-TREE.
        #==========================================================================
        set_target_properties(${target_name}
                PROPERTIES
                LINKER_LANGUAGE CXX
                ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
                LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
                RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests"
                )

        # Let's setup the target C++ standard, but only if the user did not provide it manually.
        if (NOT CMAKE_CXX_STANDARD)
            set_property(TARGET ${target_name} PROPERTY CXX_STANDARD 17)
        endif ()
        set_property(TARGET ${target_name} PROPERTY CXX_STANDARD_REQUIRED YES)
        set_property(TARGET ${target_name} PROPERTY CXX_EXTENSIONS NO)
        add_test(${target_name} "${CMAKE_BINARY_DIR}/tests/${target_name}")

        # Clean up set variables.
        unset(target_name)
        unset(dirname)
    endif ()
endfunction()

# Clean up.
#unset(_YOLO_UPPER_PROJECT_NAME)

#function(TUDAT_ADD_TESTCASE arg1)
#    if (TUDAT_BUILD_TESTS)
#        # arg1: test name
#        # arg2: test path
#        #        get_filename_component(name_top "${CMAKE_CURRENT_SOURCE_DIR}/.." NAME)
#        get_filename_component(name_bot ${CMAKE_CURRENT_SOURCE_DIR} NAME)
#        set(executable_name "test_${name_bot}_${arg1}")
#        add_executable(${executable_name} ${CMAKE_CURRENT_SOURCE_DIR}/tests/unitTest${arg1}.cpp)
#        if (CMAKE_BUILD_TYPE MATCHES DEBUG)
#            if (${arg1} MATCHES "tudat_sofa_interface")
#                message("Tudat TESTCASE: ${arg1}")
#                message("ARGN: ${ARGN}")
#            endif ()
#        endif ()
#        set_target_properties(${executable_name}
#                PROPERTIES
#                LINKER_LANGUAGE CXX
#                ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
#                LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
#                RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests"
#                COMPILE_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_USE_CXX11_ABI=0"
#                )
#        # All tests require Boost_LIBRARIES, so they are kept constant here.
#        target_link_libraries(${executable_name} PUBLIC
#                "${ARGN}"
#                "${Boost_LIBRARIES}"
#                )
#        # Has to do with a fundamental change in the way that vectors, maps
#        # etc work compared to the unit-test-framework. Not the best link for
#        # resource but:
#        # https://stackoverflow.com/questions/33644088/linker-error-while-building-unit-tests-with-boost
#
#        #        target_compile_options(${executable_name} PRIVATE
#        #                "$<$<CONFIG:Debug>:${TUDAT_CXX_FLAGS_DEBUG}>"
#        #                "$<$<CONFIG:Release>:${TUDAT_CXX_FLAGS_RELEASE}>"
#        #                "$<$<CONFIG:RelWithDebInfo>:${TUDAT_CXX_FLAGS_RELEASE}>"
#        #                "$<$<CONFIG:MinSizeRel>:${TUDAT_CXX_FLAGS_RELEASE}>"
#        #                )
#        # Let's setup the target C++ standard, but only if the user did not provide it manually.
#        if (NOT CMAKE_CXX_STANDARD)
#            set_property(TARGET ${executable_name} PROPERTY CXX_STANDARD 17)
#        endif ()
#        set_property(TARGET ${executable_name} PROPERTY CXX_STANDARD_REQUIRED YES)
#        set_property(TARGET ${executable_name} PROPERTY CXX_EXTENSIONS NO)
#        add_test(${executable_name} "${CMAKE_BINARY_DIR}/tests/${executable_name}")
#    endif ()
#endfunction()