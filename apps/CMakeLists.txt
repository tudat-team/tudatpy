cmake_minimum_required(VERSION 3.15.0)

TUDAT_ADD_EXECUTABLE(application_lro_baseline
        "src/lro_baseline.cpp"
        ${Tudat_PROPAGATION_LIBRARIES}
        )

TUDAT_ADD_EXECUTABLE(application_paneled_srp
        "src/paneled_srp.cpp"
        ${Tudat_PROPAGATION_LIBRARIES}
        )

TUDAT_ADD_EXECUTABLE(application_cannonball_srp
        "src/cannonball_srp.cpp"
        ${Tudat_PROPAGATION_LIBRARIES}
        )

TUDAT_ADD_EXECUTABLE(application_spice
        "src/spice_exe.cpp"
        ${Tudat_PROPAGATION_LIBRARIES}
        )
