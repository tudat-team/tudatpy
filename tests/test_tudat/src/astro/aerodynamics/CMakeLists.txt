#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

if (TUDAT_BUILD_WITH_PCH AND TUDAT_BUILD_WITH_TEST_PCH)
    add_library(tudat_aerodynamics_test_case_support
        aerodynamicsTestCaseSupport.cpp
    )

    target_include_directories(tudat_aerodynamics_test_case_support PUBLIC
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/tests/test_tudat/include>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    )

    target_include_directories(tudat_aerodynamics_test_case_support SYSTEM PRIVATE
        "${EIGEN3_INCLUDE_DIRS}"
        "${Boost_INCLUDE_DIRS}"
        "${CSpice_INCLUDE_DIRS}"
        "${Sofa_INCLUDE_DIRS}"
        "${TudatResources_INCLUDE_DIRS}"
    )

    target_link_libraries(tudat_aerodynamics_test_case_support PUBLIC
        ${Tudat_PROPAGATION_LIBRARIES}
    )
    target_link_libraries(tudat_aerodynamics_test_case_support PRIVATE
        "${Boost_LIBRARIES}"
    )

    set_target_properties(tudat_aerodynamics_test_case_support PROPERTIES POSITION_INDEPENDENT_CODE OFF)
    target_compile_options(tudat_aerodynamics_test_case_support PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-fPIE>
    )

    target_precompile_headers(tudat_aerodynamics_test_case_support PRIVATE
        <tudat/simulation/environment_setup/body.h>
        <tudat/simulation/environment_setup/createBodiesFactory.h>
        <tudat/simulation/environment_setup/defaultBodies.h>
        <tudat/simulation/propagation_setup/createAccelerationModels.h>
        <tudat/simulation/propagation_setup/singleArcDynamicsSimulator.h>
        <tudat/simulation/estimation_setup/createNumericalSimulator.h>
        <tudat/interface/spice/spiceInterface.h>
        <tudat/io/basicInputOutput.h>
        <Eigen/Core>
        <Eigen/Geometry>
        <limits>
        <map>
        <memory>
        <string>
        <vector>
    )

    set(TUDAT_TEST_PCH_REUSE_TARGET tudat_aerodynamics_test_case_support)
endif ()

TUDAT_ADD_TEST_CASE(AerodynamicMomentAndAerodynamicForce
        PRIVATE_LINKS
        ${Tudat_PROPAGATION_LIBRARIES}
        Tudat::tudat_aerodynamics
        )

TUDAT_ADD_TEST_CASE(AerodynamicsNamespace
        PRIVATE_LINKS
        Tudat::tudat_aerodynamics
        Tudat::tudat_basic_mathematics
        Tudat::tudat_statistics
        )

TUDAT_ADD_TEST_CASE(CoefficientGenerator
        PRIVATE_LINKS
        Tudat::tudat_aerodynamics
        Tudat::tudat_geometric_shapes
        Tudat::tudat_interpolators
        Tudat::tudat_basic_mathematics
        Tudat::tudat_statistics
        )

TUDAT_ADD_TEST_CASE(ExponentialAtmosphere
        PRIVATE_LINKS
        Tudat::tudat_aerodynamics
        Tudat::tudat_basic_mathematics
        Tudat::tudat_statistics
        )

TUDAT_ADD_TEST_CASE(CustomConstantTemperatureAtmosphere
        PRIVATE_LINKS
        Tudat::tudat_aerodynamics
        Tudat::tudat_basic_mathematics
        Tudat::tudat_statistics
        )

 TUDAT_ADD_TEST_CASE(TabulatedAtmosphere
         PRIVATE_LINKS
         Tudat::tudat_aerodynamics
         Tudat::tudat_input_output
         Tudat::tudat_interpolators
         Tudat::tudat_basic_mathematics
         Tudat::tudat_statistics
         Tudat::tudat_basic_astrodynamics
         Tudat::tudat_statistics
         )

TUDAT_ADD_TEST_CASE(HeatTransfer
        PRIVATE_LINKS
        Tudat::tudat_aerodynamics
        Tudat::tudat_root_finders
        Tudat::tudat_basic_mathematics
        Tudat::tudat_statistics
        Tudat::tudat_basic_astrodynamics
        Tudat::tudat_statistics
        )

TUDAT_ADD_TEST_CASE(ControlSurfaceIncrements
        PRIVATE_LINKS
        ${Tudat_PROPAGATION_LIBRARIES}
        )

TUDAT_ADD_TEST_CASE(AerodynamicCoefficientsFromFile
        PRIVATE_LINKS
        ${Tudat_PROPAGATION_LIBRARIES}
        )

#TUDAT_ADD_TEST_CASE(RarefieldFlowGenerator
#        PRIVATE_LINKS
#        ${Tudat_PROPAGATION_LIBRARIES}
#)

TUDAT_ADD_TEST_CASE(WindModel
        PRIVATE_LINKS
        ${Tudat_PROPAGATION_LIBRARIES}
        )

TUDAT_ADD_TEST_CASE(MarsDtmAtmosphere
        PRIVATE_LINKS
        ${Tudat_PROPAGATION_LIBRARIES}
)

TUDAT_ADD_TEST_CASE(NRLMSISE00Atmosphere
        PRIVATE_LINKS
        ${Tudat_PROPAGATION_LIBRARIES}
        )

if(TUDAT_BUILD_WITH_MCD_INTERFACE)
    TUDAT_ADD_TEST_CASE(McdAtmosphere
                PRIVATE_LINKS
                ${Tudat_PROPAGATION_LIBRARIES}
                Tudat::tudat_estimatable_parameters
                Tudat::tudat_aerodynamics
        )

    if(WIN32 AND TARGET test_aerodynamics_McdAtmosphere)
        # The MCD Fortran routines use substantial local storage while reading
        # climatology data. The default Windows executable stack can overflow
        # before the first MCD query returns.
        target_link_options(test_aerodynamics_McdAtmosphere PRIVATE
                $<$<CXX_COMPILER_ID:MSVC>:/STACK:67108864>
                )
    endif()
endif()

TUDAT_ADD_TEST_CASE(ComaModel
        PRIVATE_LINKS
        ${Tudat_PROPAGATION_LIBRARIES}
)

TUDAT_ADD_TEST_CASE(ComaWindModel
        PRIVATE_LINKS
        ${Tudat_PROPAGATION_LIBRARIES}
)
