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
 #      120717    D. Dirkx          Modified for SpiceInterface.
 #
 #    References
 #
 #    Notes
 #

#Set the source files.
set(SPICEINTERFACE_SOURCES
  "${SRCROOT}${EXTERNALDIR}/SpiceInterface/spiceEphemeris.cpp"
  "${SRCROOT}${EXTERNALDIR}/SpiceInterface/spiceRotationalEphemeris.cpp"
  "${SRCROOT}${EXTERNALDIR}/SpiceInterface/spiceInterface.cpp"
)

# Set the header files.
set(SPICEINTERFACE_HEADERS
  "${SRCROOT}${EXTERNALDIR}/SpiceInterface/spiceEphemeris.h"
  "${SRCROOT}${EXTERNALDIR}/SpiceInterface/spiceRotationalEphemeris.h"
  "${SRCROOT}${EXTERNALDIR}/SpiceInterface/spiceInterface.h"
)

# Add static libraries.
add_library(tudat_spice_interface STATIC ${SPICEINTERFACE_SOURCES} ${SPICEINTERFACE_HEADERS})
setup_tudat_library_target(tudat_spice_interface "${SRCROOT}${SPICEINTERFACEDIR}")

# Add unit tests.
add_executable(test_SpiceInterface "${SRCROOT}${EXTERNALDIR}/SpiceInterface/UnitTests/unitTestSpiceInterface.cpp")
setup_custom_test_program(test_SpiceInterface "${SRCROOT}${EXTERNALDIR}/SpiceInterface")
target_link_libraries(test_SpiceInterface tudat_ephemerides tudat_basic_mathematics tudat_spice_interface tudat_basic_astrodynamics ${SPICE_LIBRARIES} ${Boost_LIBRARIES})
