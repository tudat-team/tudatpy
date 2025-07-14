#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

## Add tests
TUDAT_ADD_TEST_CASE(Bisection
    PRIVATE_LINKS
    tudat_root_finders
    )

TUDAT_ADD_TEST_CASE( HalleyRootFinder
    PRIVATE_LINKS
    tudat_root_finders
    )

TUDAT_ADD_TEST_CASE(NewtonRaphson
    PRIVATE_LINKS
    tudat_root_finders
    )

TUDAT_ADD_TEST_CASE(SecantRootFinder
    PRIVATE_LINKS
    tudat_root_finders
    )
