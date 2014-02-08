 #    Copyright (c) 2010-2013, Delft University of Technology
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
 #      121123    D. Dirkx          Adapted for ElectroMagnetism sub-directory.
 #      130124    K. Kumar          Updated file header.
 #
 #    References
 #
 #    Notes
 #

# Set the source files.
set(ELECTROMAGNETISM_SOURCES
  "${SRCROOT}${ELECTROMAGNETISMDIR}/cannonBallRadiationPressureAcceleration.h"
  "${SRCROOT}${ELECTROMAGNETISMDIR}/cannonBallRadiationPressureForce.h"
  "${SRCROOT}${ELECTROMAGNETISMDIR}/lorentzStaticMagneticForce.h"
  "${SRCROOT}${ELECTROMAGNETISMDIR}/lorentzStaticMagneticAcceleration.h"
)

# Set the header files.
set(ELECTROMAGNETISM_HEADERS 
  "${SRCROOT}${ELECTROMAGNETISMDIR}/cannonBallRadiationPressureAcceleration.cpp"
  "${SRCROOT}${ELECTROMAGNETISMDIR}/cannonBallRadiationPressureForce.cpp"
  "${SRCROOT}${ELECTROMAGNETISMDIR}/lorentzStaticMagneticForce.cpp"
  "${SRCROOT}${ELECTROMAGNETISMDIR}/lorentzStaticMagneticAcceleration.cpp"
)

# Add static libraries.
add_library(tudat_electro_magnetism STATIC ${ELECTROMAGNETISM_SOURCES} ${ELECTROMAGNETISM_HEADERS})
setup_tudat_library_target(tudat_electro_magnetism "${SRCROOT}${ELECTROMAGNETISMDIR}")

# Add unit tests.
add_executable(test_CannonBallRadiationPressureAccelerationAndForce "${SRCROOT}${ELECTROMAGNETISMDIR}/UnitTests/unitTestCannonBallRadiationPressureAccelerationAndForce.cpp")
setup_custom_test_program(test_CannonBallRadiationPressureAccelerationAndForce "${SRCROOT}${ELECTROMAGNETISMDIR}")
target_link_libraries(test_CannonBallRadiationPressureAccelerationAndForce tudat_electro_magnetism ${Boost_LIBRARIES})

add_executable(test_LorentzStaticMagneticAccelerationAndForce "${SRCROOT}${ELECTROMAGNETISMDIR}/UnitTests/unitTestLorentzStaticMagneticAccelerationAndForce.cpp")
setup_custom_test_program(test_LorentzStaticMagneticAccelerationAndForce "${SRCROOT}${ELECTROMAGNETISMDIR}")
target_link_libraries(test_LorentzStaticMagneticAccelerationAndForce tudat_electro_magnetism ${Boost_LIBRARIES})
