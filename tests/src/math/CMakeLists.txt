#    Copyright (c) 2010-2019, Delft University of Technology
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

add_subdirectory("basic")
add_subdirectory("geometric")
add_subdirectory("interpolators")
add_subdirectory("integrators")
add_subdirectory("root_finders")
add_subdirectory("statistics")
add_subdirectory("quadrature")
if(TUDAT_BUILD_WITH_FILTERS)
  add_subdirectory("filters")
endif()
