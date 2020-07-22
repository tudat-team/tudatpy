#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#


TUDAT_ADD_TEST_CASE(ReferenceFrameTransformations PRIVATE_LINKS tudat_reference_frames)

TUDAT_ADD_TEST_CASE(ApparentAccelerationModel PRIVATE_LINKS tudat_reference_frames tudat_basic_astrodynamics)

TUDAT_ADD_TEST_CASE(AerodynamicAngleCalculator PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})
