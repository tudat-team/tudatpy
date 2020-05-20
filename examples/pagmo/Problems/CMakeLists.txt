 #    Copyright (c) 2010-2019, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.

# Set the header files.
set(MY_PAGMO_PROBLEMS_HEADERS
  "${CMAKE_CURRENT_SOURCE_DIR}/earthMarsTransfer.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/propagationTargeting.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/multipleGravityAssist.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/himmelblau.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/getAlgorithm.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/saveOptimizationResults.h"
)

set(MY_PAGMO_PROBLEMS_SOURCES
  "${CMAKE_CURRENT_SOURCE_DIR}/earthMarsTransfer.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/propagationTargeting.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/multipleGravityAssist.cpp"
)

add_library(pagmo2_library_example_problems STATIC ${MY_PAGMO_PROBLEMS_SOURCES} ${MY_PAGMO_PROBLEMS_HEADERS})
setup_tudat_library_target(pagmo2_library_example_problems "${CMAKE_CURRENT_SOURCE_DIR}")


#target_compile_options(pagmo2_library_example_tudat_problems PRIVATE "$<$<CONFIG:DEBUG>:${PAGMO_CXX_FLAGS_DEBUG}>" "$<$<CONFIG:RELEASE>:${PAGMO_CXX_FLAGS_RELEASE}>")
#set_property(TARGET pagmo2_library_example_tudat_problems PROPERTY CXX_STANDARD 11)
#set_property(TARGET pagmo2_library_example_tudat_problems PROPERTY CXX_STANDARD_REQUIRED YES)
#set_property(TARGET pagmo2_library_example_tudat_problems PROPERTY CXX_EXTENSIONS NO)
#target_link_libraries( pagmo2_library_example_tudat_problems ${TUDAT_APPLICATION_PROPAGATION_LIBRARIES}
#    ${Boost_LIBRARIES} pthread )
