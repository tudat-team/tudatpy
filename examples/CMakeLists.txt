


#if (TUDAT_BUILD_TUDAT_TUTORIALS)
#    message("Building Tudat tutorials.")
#    add_subdirectory(tudat)
#endif ()

#if (TUDAT_BUILD_WITH_PAGMO)
#    message("Building pagmo tutorials.")
#    add_subdirectory(pagmo)
#endif ()
