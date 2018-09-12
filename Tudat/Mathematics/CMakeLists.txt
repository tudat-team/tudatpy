 #    Copyright (c) 2010-2018, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.
 #
 #    Notes
 #

# Define the sub-directories.
set(BASICMATHEMATICSDIR "${MATHEMATICSDIR}/BasicMathematics")
set(GEOMETRICSHAPESDIR "${MATHEMATICSDIR}/GeometricShapes")
set(INTERPOLATORSDIR "${MATHEMATICSDIR}/Interpolators")
set(FILTERSDIR "${MATHEMATICSDIR}/Filters")
set(NUMERICALINTEGRATORSDIR "${MATHEMATICSDIR}/NumericalIntegrators")
set(ROOTFINDERSDIR "${MATHEMATICSDIR}/RootFinders")
set(STATISTICSDIR "${MATHEMATICSDIR}/Statistics")
set(NUMERICALQUADRATUREDIR "${MATHEMATICSDIR}/NumericalQuadrature")

# Add sub-directories.
add_subdirectory("${SRCROOT}${MATHEMATICSDIR}/BasicMathematics")
add_subdirectory("${SRCROOT}${MATHEMATICSDIR}/GeometricShapes")
add_subdirectory("${SRCROOT}${MATHEMATICSDIR}/Interpolators")
add_subdirectory("${SRCROOT}${MATHEMATICSDIR}/Filters")
add_subdirectory("${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators")
add_subdirectory("${SRCROOT}${MATHEMATICSDIR}/RootFinders")
add_subdirectory("${SRCROOT}${MATHEMATICSDIR}/Statistics")
add_subdirectory("${SRCROOT}${MATHEMATICSDIR}/NumericalQuadrature")

# Get target properties for static libraries.
get_target_property(BASICMATHEMATICS tudat_basic_mathematics SOURCES)
get_target_property(GEOMETRICSHAPESSOURCES tudat_geometric_shapes SOURCES)
get_target_property(INTERPOLATORSSOURCES tudat_interpolators SOURCES)
get_target_property(FILTERSSOURCES tudat_filters SOURCES)
get_target_property(NUMERICALINTEGRATORSSOURCES tudat_numerical_integrators SOURCES)
get_target_property(ROOTFINDERSSOURCES tudat_root_finders SOURCES)
get_target_property(STATISTICSSOURCES tudat_statistics SOURCES)
get_target_property(NUMERICALQUADRATURESOURCES tudat_numerical_quadrature SOURCES)
