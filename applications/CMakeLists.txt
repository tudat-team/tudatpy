cmake_minimum_required(VERSION 3.15.0)

TUDAT_ADD_LIBRARY(liblro
        "src/liblro.cpp"
        "src/liblro.h"
        PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_EXECUTABLE(application_lro_baseline
        "src/lro_baseline.cpp"
        ${Tudat_PROPAGATION_LIBRARIES} tudat_liblro
        )

TUDAT_ADD_EXECUTABLE(application_lro_json
        "src/lro_json.cpp"
        ${Tudat_PROPAGATION_LIBRARIES} tudat_liblro
        )

TUDAT_ADD_EXECUTABLE(application_lro_analysis
        "src/lro_analysis.cpp"
        ${Tudat_PROPAGATION_LIBRARIES}
        )

TUDAT_ADD_EXECUTABLE(application_spice
        "src/spice_exe.cpp"
        ${Tudat_PROPAGATION_LIBRARIES}
        )

TUDAT_ADD_EXECUTABLE(application_shadow_function_estimation
        "src/shadow_function_estimation.cpp"
        ${Tudat_PROPAGATION_LIBRARIES}
        )
