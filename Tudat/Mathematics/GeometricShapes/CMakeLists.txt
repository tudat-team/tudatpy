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
 #      111025    K. Kumar          Adapted file to work with Revision 194.
 #      111026    K. Kumar          Adapted file so all headers show in project tree in Qt Creator.
 #
 #    References
 #
 #    Notes
 #

# Add source files.
set(GEOMETRICSHAPES_SOURCES
  "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes/capsule.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes/compositeSurfaceGeometry.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes/conicalFrustum.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes/lawgsPartGeometry.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes/quadrilateralMeshedSurfaceGeometry.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes/singleSurfaceGeometry.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes/sphereSegment.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes/torus.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes/geometricShapesToFile.cpp"
)

# Add header files.
set(GEOMETRICSHAPES_HEADERS 
  "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes/capsule.h"
  "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes/compositeSurfaceGeometry.h"
  "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes/conicalFrustum.h"
  "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes/lawgsPartGeometry.h"
  "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes/quadrilateralMeshedSurfaceGeometry.h"
  "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes/singleSurfaceGeometry.h"
  "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes/sphereSegment.h"
  "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes/surfaceGeometry.h"
  "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes/torus.h"
  "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes/geometricShapesToFile.h"
)

# Add static libraries.
add_library(tudat_geometric_shapes STATIC ${GEOMETRICSHAPES_SOURCES})
setup_tudat_library_target(tudat_geometric_shapes "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes")

# Add unit tests.
add_executable(test_LawgsSurfaceGeometry "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes/UnitTests/unitTestLawgsSurfaceGeometry.cpp")
setup_custom_test_program(test_LawgsSurfaceGeometry "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes")
target_link_libraries(test_LawgsSurfaceGeometry tudat_geometric_shapes tudat_basic_mathematics ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})
