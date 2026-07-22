# Polymorphic registrations are instantiated once per domain here instead of in every consumer
# translation unit. Non-template archive/file-I/O definitions remain in their owning model
# libraries so direct per-library consumers retain the existing link interface.

set(serialization_SOURCES
        registrationsAcceleration.cpp
        registrationsEnvironment.cpp
        registrationsMath.cpp
        registrationsPropagation.cpp
        )

set(serialization_PRIVATE_LINKS
        tudat_interpolators
        tudat_environment_setup
        tudat_propagation_setup
        )

if(TUDAT_BUILD_WITH_ESTIMATION_TOOLS)
    list(APPEND serialization_SOURCES
            registrationsEstimation.cpp
            )
    list(APPEND serialization_PRIVATE_LINKS
            tudat_observation_models
            tudat_estimation_setup
            tudat_orbit_determination
            )
endif()

TUDAT_ADD_LIBRARY("serialization"
        "${serialization_SOURCES}"
        ""
        PUBLIC_SERIALIZATION
        PRIVATE_LINKS ${serialization_PRIVATE_LINKS}
        )
