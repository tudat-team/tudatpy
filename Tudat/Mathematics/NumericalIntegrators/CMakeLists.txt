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
 #      120127    B. Tong Minh      File created.
 #      121219    K. Kumar          Removed Euler integrator files (migrated to Tudat Core).
 #      130916    K. Kumar          Reformatted unit test entries.
 #
 #    References
 #
 #    Notes
 #

# Add source files.
set(NUMERICALINTEGRATORS_SOURCES
  "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/rungeKuttaCoefficients.cpp"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/numericalIntegrator.cpp"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/UnitTests/numericalIntegratorTests.cpp"
)

# Add header files.
set(NUMERICALINTEGRATORS_HEADERS 
  "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/rungeKuttaVariableStepSizeIntegrator.h"
  "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/rungeKuttaCoefficients.h"
  "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/UnitTests/burdenAndFairesNumericalIntegratorTest.h"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/euler.h"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/numericalIntegrator.h"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/reinitializableNumericalIntegrator.h"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/rungeKutta4Integrator.h"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/createNumericalIntegrator.h"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/UnitTests/numericalIntegratorTests.h"
  "${SRCROOT}${NUMERICALINTEGRATORSDIR}/UnitTests/numericalIntegratorTestFunctions.h"
)

# Add static libraries.
add_library(tudat_numerical_integrators STATIC ${NUMERICALINTEGRATORS_SOURCES} ${NUMERICALINTEGRATORS_HEADERS})
setup_tudat_library_target(tudat_numerical_integrators "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators")

# Add unit tests.
add_executable(test_EulerIntegrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/UnitTests/unitTestEulerIntegrator.cpp")
setup_custom_test_program(test_EulerIntegrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators")
target_link_libraries(test_EulerIntegrator tudat_numerical_integrators tudat_input_output ${Boost_LIBRARIES})

add_executable(test_NumericalIntegrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/UnitTests/unitTestNumericalIntegrator.cpp")
setup_custom_test_program(test_NumericalIntegrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators")
target_link_libraries(test_NumericalIntegrator tudat_numerical_integrators ${Boost_LIBRARIES})

add_executable(test_RungeKutta4Integrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/UnitTests/unitTestRungeKutta4Integrator.cpp")
setup_custom_test_program(test_RungeKutta4Integrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators")
target_link_libraries(test_RungeKutta4Integrator tudat_numerical_integrators tudat_input_output ${Boost_LIBRARIES})

add_executable(test_RungeKuttaVariableStepSizeIntegrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/UnitTests/unitTestRungeKuttaVariableStepSizeIntegrator.cpp")
setup_custom_test_program(test_RungeKuttaVariableStepSizeIntegrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators")
target_link_libraries(test_RungeKuttaVariableStepSizeIntegrator tudat_numerical_integrators ${Boost_LIBRARIES})

add_executable(test_RungeKuttaCoefficients "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/UnitTests/unitTestRungeKuttaCoefficients.cpp")
setup_custom_test_program(test_RungeKuttaCoefficients "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators")
target_link_libraries(test_RungeKuttaCoefficients tudat_numerical_integrators ${Boost_LIBRARIES})

add_executable(test_RungeKuttaFehlberg45Integrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/UnitTests/unitTestRungeKuttaFehlberg45Integrator.cpp")
setup_custom_test_program(test_RungeKuttaFehlberg45Integrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators")
target_link_libraries(test_RungeKuttaFehlberg45Integrator tudat_numerical_integrators tudat_input_output ${Boost_LIBRARIES})

add_executable(test_RungeKuttaFehlberg56Integrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/UnitTests/unitTestRungeKuttaFehlberg56Integrator.cpp")
setup_custom_test_program(test_RungeKuttaFehlberg56Integrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators")
target_link_libraries(test_RungeKuttaFehlberg56Integrator tudat_numerical_integrators tudat_input_output ${Boost_LIBRARIES})

add_executable(test_RungeKuttaFehlberg78Integrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/UnitTests/unitTestRungeKuttaFehlberg78Integrator.cpp")
setup_custom_test_program(test_RungeKuttaFehlberg78Integrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators")
target_link_libraries(test_RungeKuttaFehlberg78Integrator tudat_numerical_integrators tudat_input_output ${Boost_LIBRARIES})

add_executable(test_RungeKutta87DormandPrinceIntegrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/UnitTests/unitTestRungeKutta87DormandPrinceIntegrator.cpp")
setup_custom_test_program(test_RungeKutta87DormandPrinceIntegrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators")
target_link_libraries(test_RungeKutta87DormandPrinceIntegrator tudat_numerical_integrators tudat_input_output ${Boost_LIBRARIES})
