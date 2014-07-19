 #    Copyright (c) 2010-2014, Delft University of Technology
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
 #      111027    K. Kumar          Removed RandomNumberGenerators as sub-directory.
 #
 #    References
 #
 #    Notes
 #

# Add sub-directories.
add_subdirectory("${SRCROOT}${MATHEMATICSDIR}/BasicMathematics")
add_subdirectory("${SRCROOT}${MATHEMATICSDIR}/GeometricShapes")
add_subdirectory("${SRCROOT}${MATHEMATICSDIR}/Interpolators")
add_subdirectory("${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators")
add_subdirectory("${SRCROOT}${MATHEMATICSDIR}/RootFinders")
add_subdirectory("${SRCROOT}${MATHEMATICSDIR}/Statistics")

# Get target properties for static libraries.
get_target_property(BASICMATHEMATICS tudat_basic_mathematics SOURCES)
get_target_property(GEOMETRICSHAPESSOURCES tudat_geometric_shapes SOURCES)
get_target_property(INTERPOLATORSSOURCES tudat_interpolators SOURCES)
get_target_property(NUMERICALINTEGRATORSSOURCES tudat_numerical_integrators SOURCES)
get_target_property(ROOTFINDERSSOURCES tudat_root_finders SOURCES)
get_target_property(STATISTICSSOURCES tudat_statistics SOURCES)
