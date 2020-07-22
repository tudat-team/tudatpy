include(CMakeParseArguments)

# TODO: remove this: all tests shoulld be run, if some are omitted, a warning/error should at least be printed
set(TESTS_REQUIRING_TEST_DATA
        # The mysterious failures
#        test_basic_astro_EmpiricalAcceleration
#        test_earth_orientation_EarthOrientationCalculator
#        test_earth_orientation_EopReader
#        test_earth_orientation_PolarMotionCalculator
#        test_earth_orientation_TimeScaleConverter
#        test_earth_orientation_ShortPeriodEopCorrections
#        test_electromagnetism_PanelledRadiationPressure

        # SPICE(FILENAMETOOLONG) --
        #
        #  Input file name </home/ggarrett/miniconda3/conda-bld/tudat_1592845536964/_h_env
        #  _placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold
        #  _placehold_placehold_placehold_placehold_placehold_placehold_placehold_placehold
        #  _placehold_placehold_placehold_pl/resource/spice_kernels/pck00010.tpc> has
        #  length 291 characters. The limit on the length of file names stored by FURNSH
        #  is 255 characters.
        #
        #  A traceback follows.  The name of the highest level module is first.
        #  furnsh_c --> FURNSH
#        test_aerodynamics_AerodynamicMomentAndAerodynamicForce
#        test_aerodynamics_ControlSurfaceIncrements
#        test_aerodynamics_AerodynamicCoefficientsFromFile
#        test_aerodynamics_WindModel
#        test_spice_SpiceInterface
#        test_simulation_EnvironmentModelSetup
#        test_simulation_AccelerationModelSetup

        # Fails remotely
        #        test_aerodynamics_AerodynamicMomentAndAerodynamicForce
        #        test_aerodynamics_TabulatedAtmosphere
        #        test_aerodynamics_ControlSurfaceIncrements
        #        test_aerodynamics_WindModel
        #        test_ephemerides_ApproximatePlanetPositions
        #        test_ephemerides_TabulatedEphemeris
        #        test_spice_SpiceInterface
        #        test_simulation_EnvironmentModelSetup
        #        test_simulation_AccelerationModelSetup
        #        test_io_BasicInputOutput
        #        # Fails locally
        #        test_aerodynamics_AerodynamicCoefficientsFromFile
        #        test_basic_astro_EmpiricalAcceleration
        #        test_earth_orientation_EarthOrientationCalculator
        #        test_earth_orientation_EopReader
        #        test_earth_orientation_PolarMotionCalculator
        #        test_earth_orientation_TimeScaleConverter
        #        test_earth_orientation_ShortPeriodEopCorrections
        #        test_electromagnetism_PanelledRadiationPressure
        #        test_interpolators_CubicSplineInterpolator
        #        test_interpolators_LinearInterpolator
        #        test_interpolators_MultiLinearInterpolator
        #        test_integrators_EulerIntegrator
        #        test_integrators_RungeKutta4Integrator
        #        test_integrators_RungeKuttaFehlberg45Integrator
        #        test_integrators_RungeKuttaFehlberg78Integrator
        #        test_integrators_RungeKutta87DormandPrinceIntegrator
        #        test_quadrature_GaussianQuadrature
        #        test_io_MapTextFileReader
        #        test_io_MatrixTextFileReader

        # HERE
        #        test_io_TwoLineElementsTextFileReader
        #        test_io_MissileDatcomReader
        #        test_io_MissileDatcomData
        #        test_io_DictionaryInputSystem
        #        test_io_BasicInputOutput

        #        test_io_MultiArrayReader
        #        test_io_MultiArrayWriter
        #        test_io_AerodynamicCoefficientReader
        )

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

    if (${target_name} IN_LIST TESTS_REQUIRING_TEST_DATA)

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
