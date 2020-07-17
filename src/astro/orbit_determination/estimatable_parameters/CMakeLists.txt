 #    Copyright (c) 2010-2019, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.

# Set the source files.
set(estimatable_parameters_SOURCES
  "estimatableParameter.cpp"
  "sphericalHarmonicSineCoefficients.cpp"
  "sphericalHarmonicCosineCoefficients.cpp"
  "tidalLoveNumber.cpp"
)

# Set the header files.
set(estimatable_parameters_HEADERS
  "estimatableParameter.h"
  "gravitationalParameter.h"
  "initialTranslationalState.h"
  "initialRotationalState.h"
  "radiationPressureCoefficient.h"
  "sphericalHarmonicSineCoefficients.h"
  "sphericalHarmonicCosineCoefficients.h"
  "constantRotationRate.h"
  "constantDragCoefficient.h"
  "constantRotationalOrientation.h"
  "observationBiasParameter.h"
  "ppnParameters.h"
  "groundStationPosition.h"
  "equivalencePrincipleViolationParameter.h"
  "empiricalAccelerationCoefficients.h"
  "tidalLoveNumber.h"
  "directTidalTimeLag.h"
  "meanMomentOfInertiaParameter.h"
  "desaturationDeltaV.h"
  "periodicSpinVariation.h"
  "polarMotionAmplitude.h"
  "coreFactor.h"
  "freeCoreNutationRate.h"
)


TUDAT_ADD_LIBRARY("estimatable_parameters"
        "${estimatable_parameters_SOURCES}"
        "${estimatable_parameters_HEADERS}"
#        PUBLIC_LINKS
#        tudat_input_output
#        tudat_geometric_shapes
#        tudat_reference_frames
#        tudat_interpolators
        PRIVATE_LINKS "${Boost_LIBRARIES}"
        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}"
        )

