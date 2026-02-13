#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.

TUDAT_ADD_TEST_CASE(MetricPartials
        PRIVATE_LINKS
        ${Tudat_ESTIMATION_LIBRARIES}
        Tudat::tudat_electromagnetism
        )

TUDAT_ADD_TEST_CASE(SolarSystemMetricPartials
        PRIVATE_LINKS
        ${Tudat_ESTIMATION_LIBRARIES}
        Tudat::tudat_electromagnetism
        )
