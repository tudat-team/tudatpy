#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#
#    Note that this license only applies to this CMakeLists.txt file, not the
#    sofa library itself
#

# Specify minimum CMake version required to compile Sofa.
cmake_minimum_required(VERSION 3.3)

# Trick for recompiling changes for BUILD-TREE.
if (EXISTS "${CMAKE_BINARY_DIR}/CMakeCache.txt")
    file(REMOVE "${CMAKE_BINARY_DIR}/CMakeCache.txt")
endif ()

# Specific project name.
project(Sofa VERSION 1.0.0 LANGUAGES C)

# Get lower string name of PROJECT_NAME for directory configuration.
string(TOLOWER ${PROJECT_NAME} PROJECT_NAME_LOWER)
string(TOUPPER ${PROJECT_NAME} PROJECT_NAME_UPPER)

# +============================================================================
# BUILD OPTIONS
#  Offer the user the choice of defining the build variation.
# +============================================================================
# Build static library instead of dynamic.
option(SOFA_BUILD_STATIC_LIBRARY "Build sofa as a static library, instead of dynamic." OFF)

# +============================================================================
# INSTALL TREE CONFIGURATION
#  Offer the user the choice of overriding the installation directories.
# +============================================================================
set(${PROJECT_NAME_UPPER}_INSTALL_LIB_DIR     "${CMAKE_INSTALL_PREFIX}/lib/"                       CACHE PATH "Installation directory for libraries")
set(${PROJECT_NAME_UPPER}_INSTALL_BIN_DIR     "${CMAKE_INSTALL_PREFIX}/bin/"                       CACHE PATH "Installation directory for executables")
set(${PROJECT_NAME_UPPER}_INSTALL_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}/include/"                   CACHE PATH "Installation directory for header files")
set(${PROJECT_NAME_UPPER}_INSTALL_DATA_DIR    "${CMAKE_INSTALL_PREFIX}/data/${PROJECT_NAME_LOWER}" CACHE PATH "Installation directory for data files")
set(${PROJECT_NAME_UPPER}_INSTALL_MAN_DIR     "${CMAKE_INSTALL_PREFIX}/"                           CACHE PATH "Installation directory for man documentation")
set(${PROJECT_NAME_UPPER}_INSTALL_TESTS_DIR   "${CMAKE_BINARY_DIR}/tests/"                         CACHE PATH "Installation directory for tests (default=OFF)")
set(${PROJECT_NAME_UPPER}_INSTALL_CMAKE_DIR   "${${PROJECT_NAME_UPPER}_INSTALL_LIB_DIR}/cmake/${PROJECT_NAME_LOWER}"     CACHE PATH "Installation directory for cmake config files")

# Make relative paths absolute (needed later on)
foreach (p LIB BIN INCLUDE DATA CMAKE)
    set(var ${PROJECT_NAME_UPPER}_INSTALL_${p}_DIR)
    set(${PROJECT_NAME_UPPER}_RELATIVE_INSTALL_${p}_DIR ${${PROJECT_NAME_UPPER}_INSTALL_${p}_DIR})
    if (NOT IS_ABSOLUTE "${${var}}")
        set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
    endif ()
endforeach ()

# Find relative paths for cmake-config file. (for find_package)
set(${PROJECT_NAME}_CMAKE_DIR       "${${PROJECT_NAME_UPPER}_INSTALL_CMAKE_DIR}")
file(RELATIVE_PATH rel_include_dir  "${${PROJECT_NAME}_CMAKE_DIR}" "${${PROJECT_NAME_UPPER}_INSTALL_INCLUDE_DIR}")
file(RELATIVE_PATH rel_lib_dir      "${${PROJECT_NAME}_CMAKE_DIR}" "${${PROJECT_NAME_UPPER}_INSTALL_LIB_DIR}")
file(RELATIVE_PATH rel_data_dir     "${${PROJECT_NAME}_CMAKE_DIR}" "${${PROJECT_NAME_UPPER}_INSTALL_DATA_DIR}")
file(RELATIVE_PATH rel_bin_dir      "${${PROJECT_NAME}_CMAKE_DIR}" "${${PROJECT_NAME_UPPER}_INSTALL_BIN_DIR}")
file(RELATIVE_PATH rel_man_dir      "${${PROJECT_NAME}_CMAKE_DIR}" "${${PROJECT_NAME_UPPER}_INSTALL_MAN_DIR}")

# Set relative paths for config.cmake.
foreach (p include lib data bin man)
    string(TOUPPER ${p} P)
    set(${PROJECT_NAME}_RELATIVE_INSTALL_${P}_DIR ${rel_${p}_dir})
endforeach ()

# +============================================================================
# CMAKE BUILD TREE CONFIGURATION
# +============================================================================
set(${PROJECT_NAME_UPPER}_BINARY_LIB_DIR      "${CMAKE_BINARY_DIR}/lib/"                            CACHE PATH "Build directory for libraries")
set(${PROJECT_NAME_UPPER}_BINARY_BIN_DIR      "${CMAKE_BINARY_DIR}/bin/"                            CACHE PATH "Build directory for executables")
set(${PROJECT_NAME_UPPER}_BINARY_INCLUDE_DIR  "${CMAKE_BINARY_DIR}/include/"                        CACHE PATH "Build directory for header files")
set(${PROJECT_NAME_UPPER}_BINARY_DATA_DIR     "${CMAKE_INSTALL_PREFIX}/data/${PROJECT_NAME_LOWER}/" CACHE PATH "Build directory for data files (default=OFF)")
set(${PROJECT_NAME_UPPER}_BINARY_MAN_DIR      "${CMAKE_BINARY_DIR}/"                                CACHE PATH "Build directory for man documentation")
set(${PROJECT_NAME_UPPER}_BINARY_TESTS_DIR    "${CMAKE_BINARY_DIR}/tests/"                          CACHE PATH "Build directory for tests")
set(${PROJECT_NAME_UPPER}_BINARY_CMAKE_DIR    "${BINARY_LIB_DIR}/cmake/${PROJECT_NAME_LOWER}"       CACHE PATH "Build directory for cmake config files")

# Set the source files.
set(SOFA_SOURCE_FILES
        src/ut1tt.c
        src/epb.c
        src/tttai.c
        src/atio13.c
        src/ir.c
        src/fk45z.c
        src/epb2jd.c
        src/epj.c
        src/c2ixys.c
        src/ldn.c
        src/seps.c
        src/tpors.c
        src/xys00b.c
        src/fw2xy.c
        src/pmsafe.c
        src/falp03.c
        src/apco.c
        src/sp00.c
        src/s2p.c
        src/pvppv.c
        src/pn.c
        src/fane03.c
        src/gst06.c
        src/prec76.c
        src/s06.c
        src/pdp.c
        src/pvtob.c
        src/apci.c
        src/bpn2xy.c
        src/obl80.c
        src/pr00.c
        src/gd2gc.c
        src/rxr.c
        src/num06a.c
        src/pn00a.c
        src/pn00b.c
        src/pap.c
        src/hfk5z.c
        src/pnm06a.c
        src/pb06.c
        src/num00a.c
        src/tpstv.c
        src/p2pv.c
        src/ee00.c
        src/anp.c
        src/pvmpv.c
        src/ld.c
        src/aticqn.c
        src/gc2gde.c
        src/pvstar.c
        src/pvm.c
        src/a2af.c
        src/cr.c
        src/faur03.c
        src/atoc13.c
        src/gmst82.c
        src/fame03.c
        src/gd2gce.c
        src/aticq.c
        src/ppsp.c
        src/c2t00b.c
        src/faf03.c
        src/atoiq.c
        src/xy06.c
        src/d2tf.c
        src/ltpb.c
        src/bp00.c
        src/tttcg.c
        src/obl06.c
        src/pn06.c
        src/s00a.c
        src/a2tf.c
        src/cpv.c
        src/gmst06.c
        src/faom03.c
        src/numat.c
        src/pom00.c
        src/ltpecl.c
        src/pnm00b.c
        src/c2ibpn.c
        src/eqec06.c
        src/nutm80.c
        src/s00.c
        src/pvdpv.c
        src/num00b.c
        src/atic13.c
        src/tcbtdb.c
        src/atciqz.c
        src/c2t00a.c
        src/sxpv.c
        src/zp.c
        src/rv2m.c
        src/pmat06.c
        src/fk54z.c
        src/pvxpv.c
        src/s2xpv.c
        src/fk524.c
        src/tpxes.c
        src/c2i00a.c
        src/pn00.c
        src/tporv.c
        src/c2ixy.c
        src/pvu.c
        src/fk5hip.c
        src/pv2s.c
        src/refco.c
        src/hd2pa.c
        src/pmpx.c
        src/fad03.c
        src/ltecm.c
        src/c2i00b.c
        src/rxp.c
        src/tdbtt.c
        src/fk5hz.c
        src/xys06a.c
        src/pnm80.c
        src/ee00b.c
        src/apci13.c
        src/c2t06a.c
        src/dat.c
        src/tcgtt.c
        src/ae2hd.c
        src/era00.c
        src/tf2a.c
        src/p06e.c
        src/apcg13.c
        src/rz.c
        src/apio.c
        src/lteceq.c
        src/pmat00.c
        src/utcut1.c
        src/tpxev.c
        src/zpv.c
        src/ltp.c
        src/atco13.c
        src/s00b.c
        src/rm2v.c
        src/jd2cal.c
        src/eors.c
        src/ut1utc.c
        src/tdbtcb.c
        src/hd2ae.c
        src/fw2m.c
        src/trxp.c
        src/fal03.c
        src/lteqec.c
        src/taiutc.c
        src/pnm00a.c
        src/pv2p.c
        src/gst94.c
        src/af2a.c
        src/bp06.c
        src/ee00a.c
        src/gst06a.c
        src/atciq.c
        src/p2s.c
        src/c2s.c
        src/taiut1.c
        src/g2icrs.c
        src/gmst00.c
        src/rxpv.c
        src/starpm.c
        src/faju03.c
        src/c2i06a.c
        src/atoi13.c
        src/dtf2d.c
        src/tr.c
        src/anpm.c
        src/nut06a.c
        src/apcs13.c
        src/ry.c
        src/sepp.c
        src/ppp.c
        src/starpv.c
        src/cp.c
        src/c2tpe.c
        src/gst00b.c
        src/pm.c
        src/atci13.c
        src/c2txy.c
        src/dtdb.c
        src/eect00.c
        src/h2fk5.c
        src/fapa03.c
        src/sxp.c
        src/pvup.c
        src/ecm06.c
        src/tttdb.c
        src/tpsts.c
        src/ldsun.c
        src/aper.c
        src/icrs2g.c
        src/ttut1.c
        src/apio13.c
        src/zr.c
        src/ee06a.c
        src/utctai.c
        src/atciqn.c
        src/apcg.c
        src/taitt.c
        src/bi00.c
        src/ab.c
        src/eqeq94.c
        src/apcs.c
        src/pn06a.c
        src/pfw06.c
        src/cal2jd.c
        src/gst00a.c
        src/c2tcio.c
        src/fasa03.c
        src/nut80.c
        src/apco13.c
        src/c2teqx.c
        src/epv00.c
        src/pxp.c
        src/epj2jd.c
        src/fk52h.c
        src/aper13.c
        src/fae03.c
        src/ltpequ.c
        src/eo06a.c
        src/pmat76.c
        src/jdcalf.c
        src/s2c.c
        src/nut00a.c
        src/gc2gd.c
        src/tf2d.c
        src/nut00b.c
        src/eceq06.c
        src/pas.c
        src/trxpv.c
        src/fama03.c
        src/fk425.c
        src/s06a.c
        src/pmp.c
        src/xys00a.c
        src/atioq.c
        src/rx.c
        src/d2dtf.c
        src/eform.c
        src/fave03.c
        src/s2pv.c
        src/ut1tai.c
        src/plan94.c
        )

# Setup of library.
if (SOFA_BUILD_STATIC_LIBRARY)
    add_library(sofa STATIC "${SOFA_SOURCE_FILES}")
else ()
    add_library(sofa SHARED "${SOFA_SOURCE_FILES}")
    # Link to .m lib
    target_link_libraries(sofa PRIVATE m)
endif ()

# Add headers to library target.
target_include_directories(sofa PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
        )

# Setup of namespaced library.
add_library(Sofa::sofa ALIAS sofa)

# Testing.
enable_testing()
add_executable(sofa-test tests/t_sofa_c.c)
target_link_libraries(sofa-test Sofa::sofa)
add_test(sofa-test sofa-test)

# Installation of the library.
install(TARGETS sofa
        EXPORT sofa_export
        LIBRARY DESTINATION ${${PROJECT_NAME_UPPER}_RELATIVE_INSTALL_LIB_DIR}
        ARCHIVE DESTINATION ${${PROJECT_NAME_UPPER}_RELATIVE_INSTALL_LIB_DIR}
        RUNTIME DESTINATION ${${PROJECT_NAME_UPPER}_RELATIVE_INSTALL_BIN_DIR}
        )

# Set <project>_LIBRARIES for <project>-config.cmake.
set(Sofa_LIBRARIES "Sofa::sofa")

# Installation of the header files.
install(DIRECTORY include/ DESTINATION "${${PROJECT_NAME_UPPER}_RELATIVE_INSTALL_INCLUDE_DIR}")

# Installation of the data files.
#install(DIRECTORY data/ DESTINATION "${${PROJECT_NAME_UPPER}_RELATIVE_INSTALL_DATA_DIR}")

# Installation of the export modules.
install(EXPORT sofa_export NAMESPACE Sofa:: DESTINATION "${${PROJECT_NAME_UPPER}_RELATIVE_INSTALL_CMAKE_DIR}")

# Add all targets to the build-tree export set
#export(TARGETS ${Sofa_LIBRARIES} NAMESPACE Sofa:: FILE "${${PROJECT_NAME_UPPER}_BINARY_CMAKE_DIR}/sofa_export.cmake")

# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
#set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)
#export(PACKAGE Sofa)

# Configure cmake files.
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/sofa-config.cmake.in" "${CMAKE_CURRENT_BINARY_DIR}/sofa-config.cmake" @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/sofa-config.cmake" DESTINATION "${${PROJECT_NAME_UPPER}_INSTALL_CMAKE_DIR}")

# Take care of versioning.
include(CMakePackageConfigHelpers)
# NOTE: SameMinorVersion available only
# since CMake 3.11.
if (${CMAKE_VERSION} VERSION_LESS "3.11.0")
    write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/sofa-config-version.cmake" VERSION ${Sofa_VERSION}
            COMPATIBILITY SameMajorVersion)
else ()
    write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/sofa-config-version.cmake" VERSION ${Sofa_VERSION}
            COMPATIBILITY SameMinorVersion)
endif ()
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/sofa-config-version.cmake" DESTINATION "${${PROJECT_NAME_UPPER}_RELATIVE_INSTALL_CMAKE_DIR}")


# End of file.
