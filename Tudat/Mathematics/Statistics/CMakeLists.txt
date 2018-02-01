 #    Copyright (c) 2010-2018, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.
 #

# Add source files.
set(STATISTICS_SOURCES
  "${SRCROOT}${MATHEMATICSDIR}/Statistics/basicStatistics.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/Statistics/simpleLinearRegression.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/Statistics/multiVariateGaussianProbabilityDistributions.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/Statistics/continuousProbabilityDistributions.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/Statistics/boostProbabilityDistributions.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/Statistics/kernelDensityDistribution.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/Statistics/randomSampling.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/Statistics/randomVariableGenerator.cpp"
)

# Add header files.
set(STATISTICS_HEADERS 
  "${SRCROOT}${MATHEMATICSDIR}/Statistics/basicStatistics.h"
  "${SRCROOT}${MATHEMATICSDIR}/Statistics/simpleLinearRegression.h"
  "${SRCROOT}${MATHEMATICSDIR}/Statistics/multiVariateGaussianProbabilityDistributions.h"
  "${SRCROOT}${MATHEMATICSDIR}/Statistics/continuousProbabilityDistributions.h"
  "${SRCROOT}${MATHEMATICSDIR}/Statistics/boostProbabilityDistributions.h"
  "${SRCROOT}${MATHEMATICSDIR}/Statistics/kernelDensityDistribution.h"
  "${SRCROOT}${MATHEMATICSDIR}/Statistics/randomSampling.h"
  "${SRCROOT}${MATHEMATICSDIR}/Statistics/randomVariableGenerator.h"
)

# Add static libraries.
add_library(tudat_statistics STATIC ${STATISTICS_SOURCES} ${STATISTICS_HEADERS})
setup_tudat_library_target(tudat_statistics "${SRCROOT}${MATHEMATICSDIR}/Statistics")

# Add unit tests.
add_executable(test_SimpleLinearRegression "${SRCROOT}${MATHEMATICSDIR}/Statistics/UnitTests/unitTestSimpleLinearRegression.cpp")
setup_custom_test_program(test_SimpleLinearRegression "${SRCROOT}${MATHEMATICSDIR}/Statistics")
target_link_libraries(test_SimpleLinearRegression tudat_statistics ${Boost_LIBRARIES})

add_executable(test_BasicStatistics "${SRCROOT}${MATHEMATICSDIR}/Statistics/UnitTests/unitTestBasicStatistics.cpp")
setup_custom_test_program(test_BasicStatistics "${SRCROOT}${MATHEMATICSDIR}/Statistics")
target_link_libraries(test_BasicStatistics tudat_statistics ${Boost_LIBRARIES})

add_executable(test_MultiVariateProbabilityDistributions "${SRCROOT}${MATHEMATICSDIR}/Statistics/UnitTests/unitTestMultiVariateGaussianProbabilityDistributions.cpp")
setup_custom_test_program(test_MultiVariateProbabilityDistributions "${SRCROOT}${MATHEMATICSDIR}/Statistics")
target_link_libraries(test_MultiVariateProbabilityDistributions tudat_statistics ${Boost_LIBRARIES})

add_executable(test_BoostDistributions "${SRCROOT}${MATHEMATICSDIR}/Statistics/UnitTests/unitTestBoostRandomVariables.cpp")
setup_custom_test_program(test_BoostDistributions "${SRCROOT}${MATHEMATICSDIR}/Statistics")
target_link_libraries(test_BoostDistributions tudat_statistics ${Boost_LIBRARIES})

add_executable(test_KernelDensityDistribution "${SRCROOT}${MATHEMATICSDIR}/Statistics/UnitTests/unitTestKernelDensityDistribution.cpp")
setup_custom_test_program(test_KernelDensityDistribution "${SRCROOT}${MATHEMATICSDIR}/Statistics")
target_link_libraries(test_KernelDensityDistribution tudat_statistics ${Boost_LIBRARIES})

add_executable(test_RandomSampling "${SRCROOT}${MATHEMATICSDIR}/Statistics/UnitTests/unitTestRandomSampling.cpp")
setup_custom_test_program(test_RandomSampling "${SRCROOT}${MATHEMATICSDIR}/Statistics")
target_link_libraries(test_RandomSampling tudat_statistics ${TUDAT_EXTERNAL_LIBRARIES} ${Boost_LIBRARIES})
