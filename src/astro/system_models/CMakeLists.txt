#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

# Set the header files.
set(system_models_HEADERS
  "engineModel.h"
  "vehicleSystems.h"
)

set(system_models_SOURCES
  "dummySystemModel.cpp"
)

TUDAT_ADD_LIBRARY("system_models"
     "${system_models_SOURCES}"
     "${system_models_HEADERS}"
     )

# Add static libraries.
#add_library(tudat_system_models STATIC ${SYSTEMMODELS_SOURCES} ${SYSTEMMODELS_HEADERS})
#setup_tudat_library_target(tudat_system_models "${SRCROOT}${SYSTEMMODELSDIR}")
#set_target_properties(tudat_system_models PROPERTIES LINKER_LANGUAGE CXX)

