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
 #      12xxxx    B. Tong Minh      File created based on FindEigen3.cmake.
 #      12xxxx    P. Musegaas
 #      12xxxx    D. Dirkx          Adapted to detect the SPICE library.
 #      140127    D. Dirkx          Adapted for custom Spice kernel folder.
 #
 #
 #    References
 #      FindEigen3.cmake.
 #
 #    Notes
 #      This script tries to find SPICE library.
 #
 #      Original copyright statements (from FindEigen3.cmake:
 #          Copyright (c) 2006, 2007 Montel Laurent, <montel@kde.org>
 #          Copyright (c) 2008, 2009 Gael Guennebaud, <g.gael@free.fr>
 #          Copyright (c) 2009 Benoit Jacob <jacob.benoit.1@gmail.com>
 #
 #      FindEigen3.cmake states that redistribution and use is allowed according to the terms of
 #      the 2-clause BSD license.

if(NOT SPICE_BASE_PATH)
    find_path(SPICE_BASE_PATH NAMES SpiceUsr.h
        PATHS
            ${PROJECT_SOURCE_DIR}/External
        PATH_SUFFIXES cspice/include
)
endif(NOT SPICE_BASE_PATH)

if(NOT SPICE_BASE_PATH)
  message(STATUS "WARNING: SPICE not found! Make sure SPICE is installed or set USE_CSPICE to 'false'.")
endif(NOT SPICE_BASE_PATH)

set(SPICE_BASE_PATH ${SPICE_BASE_PATH}/..)
set(SPICE_INCLUDE_DIR ${SPICE_BASE_PATH}/..)
set(SPICE_LIBRARIES_DIR ${SPICE_BASE_PATH}/lib)
set(SPICE_LIBRARIES "cspice")
link_directories(${SPICE_LIBRARIES_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SPICE DEFAULT_MSG SPICE_INCLUDE_DIR)

mark_as_advanced(SPICE_INCLUDE_DIR)
