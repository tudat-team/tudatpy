 #    Copyright (c) 2010-2015, Delft University of Technology
 #    All rights reserved.
 #
 #    Redistribution and use in source and binary forms, with or without modification, are
 #    permitted provided that the following conditions are met:
 #      - Redistributions of source code must retain the above copyright notice, this list of
 #        conditions and the following disclaimer.
 #      - Redistributions in binary form must reproduce the above copyright notice, this list of
 #        conditions and the following disclaimer in the documentation and/or other materials
 #        provided with the distribution.
 #      - Neither the name of the Delft University of Technology nor the names of its contributors
 #        may be used to endorse or promote products derived from this software without specific
 #        prior written permission.
 #
 #    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 #    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 #    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 #    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 #    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 #    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 #    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 #    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 #    OF THE POSSIBILITY OF SUCH DAMAGE.
 #
 #    Changelog
 #      YYMMDD    Author            Comment
 #      110820    S.M. Persson      File created.
 #      111025    K. Kumar          Adapted file to work with Revision 194.
 #      111026    K. Kumar          Adapted file so all headers show in project tree in Qt Creator.
 #
 #    References
 #
 #    Notes
 #

# Add source files.
set(PROPAGATORS_SOURCES
  "${SRCROOT}${PROPAGATORSDIR}/nBodyStateDerivative.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/nBodyEnckeStateDerivative.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/variationalEquations.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/stateTransitionMatrixInterface.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/environmentUpdateTypes.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/singleStateTypeDerivative.cpp"
)

# Add header files.
set(PROPAGATORS_HEADERS
  "${SRCROOT}${PROPAGATORSDIR}/centralBodyData.h"
  "${SRCROOT}${PROPAGATORSDIR}/nBodyStateDerivative.h"
  "${SRCROOT}${PROPAGATORSDIR}/nBodyCowellStateDerivative.h"
  "${SRCROOT}${PROPAGATORSDIR}/nBodyEnckeStateDerivative.h"
  "${SRCROOT}${PROPAGATORSDIR}/dynamicsStateDerivativeModel.h"
  "${SRCROOT}${PROPAGATORSDIR}/singleStateTypeDerivative.h"
  "${SRCROOT}${PROPAGATORSDIR}/integrateEquations.h"
  "${SRCROOT}${PROPAGATORSDIR}/bodyMassStateDerivative.h"
  "${SRCROOT}${PROPAGATORSDIR}/variationalEquations.h"
  "${SRCROOT}${PROPAGATORSDIR}/stateTransitionMatrixInterface.h"
  "${SRCROOT}${PROPAGATORSDIR}/environmentUpdateTypes.h"
)

# Add static libraries.
add_library(tudat_propagators STATIC ${PROPAGATORS_SOURCES} ${PROPAGATORS_HEADERS})
setup_tudat_library_target(tudat_propagators "${SRCROOT}${PROPAGATORSDIR}")

# Add unit tests.
add_executable(test_CentralBodyData "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestCentralBodyData.cpp")
setup_custom_test_program(test_CentralBodyData "${SRCROOT}${PROPAGATORSDIR}")
target_link_libraries(test_CentralBodyData tudat_propagators tudat_basic_mathematics ${Boost_LIBRARIES})

if(USE_CSPICE)

add_executable(test_CowellStateDerivative "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestCowellStateDerivative.cpp")
setup_custom_test_program(test_CowellStateDerivative "${SRCROOT}${PROPAGATORSDIR}")
target_link_libraries(test_CowellStateDerivative ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_EnckeStateDerivative "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestEnckeStateDerivative.cpp")
setup_custom_test_program(test_EnckeStateDerivative "${SRCROOT}${PROPAGATORSDIR}")
target_link_libraries(test_EnckeStateDerivative ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})


add_executable(test_SequentialVariationEquationIntegration "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestSequentialVariationalEquationIntegration.cpp")
setup_custom_test_program(test_SequentialVariationEquationIntegration "${SRCROOT}${PROPAGATORSDIR}")
target_link_libraries(test_SequentialVariationEquationIntegration ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_VariationalEquations "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestVariationalEquationPropagation.cpp")
setup_custom_test_program(test_VariationalEquations "${SRCROOT}${PROPAGATORSDIR}")
target_link_libraries(test_VariationalEquations ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_EnvironmentModelUpdater "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestEnvironmentUpdater.cpp")
setup_custom_test_program(test_EnvironmentModelUpdater "${SRCROOT}${PROPAGATORSDIR}")
target_link_libraries(test_EnvironmentModelUpdater ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})


add_executable(test_BodyMassPropagation "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestBodyMassPropagation.cpp")
setup_custom_test_program(test_BodyMassPropagation "${SRCROOT}${PROPAGATORSDIR}")
target_link_libraries(test_BodyMassPropagation ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})


add_executable(test_MultiTypeStatePropagation "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestMultiTypeStatePropagation.cpp")
setup_custom_test_program(test_MultiTypeStatePropagation "${SRCROOT}${PROPAGATORSDIR}")
target_link_libraries(test_MultiTypeStatePropagation ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})


add_executable(test_DependentVariableOutput "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestDependentVariableOutput.cpp")
setup_custom_test_program(test_DependentVariableOutput "${SRCROOT}${PROPAGATORSDIR}")
target_link_libraries(test_DependentVariableOutput ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_StoppingConditions "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestStoppingConditions.cpp")
setup_custom_test_program(test_StoppingConditions "${SRCROOT}${PROPAGATORSDIR}")
target_link_libraries(test_StoppingConditions ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

endif( )
