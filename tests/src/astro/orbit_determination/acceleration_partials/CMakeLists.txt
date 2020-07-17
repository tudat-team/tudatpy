#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.

if(TUDAT_BUILD_WITH_SPICE_INTERFACE)

    TUDAT_ADD_TEST_CASE(AccelerationPartials
            PRIVATE_LINKS
            ${Tudat_PROPAGATION_LIBRARIES}
            ${CSpice_LIBRARIES}
            ${Sofa_LIBRARIES}
            )

    TUDAT_ADD_TEST_CASE(SphericalHarmonicPartials
            PRIVATE_LINKS
#            Tudat_ESTIMATION_LIBRARIES
            )

    TUDAT_ADD_TEST_CASE(MutualSphericalHarmonicPartials
            PRIVATE_LINKS
#            Tudat_ESTIMATION_LIBRARIES
            )

endif( )


