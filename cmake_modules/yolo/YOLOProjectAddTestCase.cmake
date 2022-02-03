include(CMakeParseArguments)

# TODO: remove this: all tests shoulld be run, if some are omitted, a warning/error should at least be printed
if (TUDAT_SKIP_JSON_TESTS)
    # https://github.com/tudat-team/tudat/issues/8
    set(TEST_TO_BE_SKIPPED
            test_json_Acceleration
            test_json_Aerodynamics
            test_json_Body
            test_json_Ephemeris
            test_json_GroundStation
            test_json_Interpolation
            test_json_Propagator
            test_json_SimulationSingleSatellite
            test_json_SimulationSinglePerturbedSatellite
            test_json_SimulationInnerSolarSystem
            test_json_SimulationGalileoConstellation
            test_json_SimulationThrustAlongVelocityVector
            test_json_SimulationThrustAccelerationFromFile
            test_json_State
            test_json_Thrust
            )
else ()
    set(TEST_TO_BE_SKIPPED
            )
endif ()

if (TUDAT_SKIP_BROKEN_MSVC_CLANG_PRECISION_TESTS)
    # https://github.com/tudat-team/tudat/issues/7
    set(TEST_TO_BE_SKIPPED ${TEST_TO_BE_SKIPPED})
endif ()


function("TUDAT_ADD_TEST_CASE" arg1)
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

    if (${target_name} IN_LIST TEST_TO_BE_SKIPPED)

    else ()
        # Add executable.
        add_executable(${target_name} ${CMAKE_CURRENT_SOURCE_DIR}/unitTest${arg1}.cpp ${PARSED_ARGS_SOURCES})

        #==========================================================================
        # TARGET-CONFIGURATION.
        #==========================================================================
        # NOTE: make sure the include directories from the current build
        # are included first, so that if there is already a pagmo installation
        # in the prefix path we don't risk including the headers from that
        # one instead.
        target_include_directories("${target_name}" PUBLIC
                $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>  # Configured test headers
                $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/tests/include>  # Test specific headers
                $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>        # Project headers
                $<INSTALL_INTERFACE:include>
                )                           # Installed headers

        target_include_directories("${target_name}"
                SYSTEM PRIVATE
                "${EIGEN3_INCLUDE_DIRS}"
                "${Boost_INCLUDE_DIRS}"
                "${CSpice_INCLUDE_DIRS}"
                "${Sofa_INCLUDE_DIRS}"
                "${TudatResources_INCLUDE_DIRS}"
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
                ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib"
                LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib"
                RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/tests"
                )

        # Let's setup the target C++ standard, but only if the user did not provide it manually.
        if (NOT CMAKE_CXX_STANDARD)
            set_property(TARGET ${target_name} PROPERTY CXX_STANDARD 17)
        endif ()
        set_property(TARGET ${target_name} PROPERTY CXX_STANDARD_REQUIRED YES)
        set_property(TARGET ${target_name} PROPERTY CXX_EXTENSIONS NO)
        add_test(${target_name} "${PROJECT_BINARY_DIR}/tests/${target_name}")

        if (TUDAT_INSTALL_TESTS)
            #==========================================================================
            # INSTALL-TREE.
            #==========================================================================
            install(TARGETS "${target_name}"
                    EXPORT tudat_export
                    LIBRARY DESTINATION "${INSTALL_LIB_DIR}"
                    ARCHIVE DESTINATION "${INSTALL_LIB_DIR}"
                    RUNTIME DESTINATION "${INSTALL_BIN_DIR}/tudat/tests"
                    )
        endif (TUDAT_INSTALL_TESTS)

        # Clean up set variables.
        unset(target_name)
        unset(dirname)
    endif ()
endfunction()
