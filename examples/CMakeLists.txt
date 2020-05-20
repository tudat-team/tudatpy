#
#add_subdirectory( "${CMAKE_CURRENT_SOURCE_DIR}/satellitePropagatorExamples/SatellitePropagatorExamples/")
#
message("HI(${TUDAT_BUILD_TUDAT_TUTORIALS}")
message("HI(${TUDAT_BUILD_PAGMO_TUTORIALS}")
if (TUDAT_BUILD_TUDAT_TUTORIALS)
    message("Building Tudat tutorials.")
    add_subdirectory(tudat)
endif ()

if (TUDAT_BUILD_PAGMO_TUTORIALS)
    message("Building pagmo tutorials.")
    add_subdirectory(pagmo)
endif ()
