#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

# Add source files.
set(data_HEADERS downloadFile.h)

# Add header files.
set(data_SOURCES downloadFile.cpp)

# Add library.
TUDAT_ADD_LIBRARY("data"
        "${data_SOURCES}"
        "${data_HEADERS}"
        PUBLIC_LINKS "${CURL_LIBRARIES}"
        )


