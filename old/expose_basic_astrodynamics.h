//
// Created by ggarrett on 29-04-20.
//

#ifndef TUDATBUNDLE_EXPOSE_BASIC_ASTRODYNAMICS_H
#define TUDATBUNDLE_EXPOSE_BASIC_ASTRODYNAMICS_H

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include "tudat/astro/basic_astro/accelerationModelTypes.h"
#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/astro/basic_astro/attitudeElementConversions.h"
#include "tudat/astro/basic_astro/celestialBodyConstants.h"
#include "tudat/astro/basic_astro/convertMeanToEccentricAnomalies.h"
#include "tudat/astro/basic_astro/clohessyWiltshirePropagator.h"
#include "tudat/astro/basic_astro/customTorque.h"
#include "tudat/astro/basic_astro/geodeticCoordinateConversions.h"
#include "tudat/astro/basic_astro/keplerPropagator.h"
#include "tudat/astro/basic_astro/missionGeometry.h"
#include "tudat/astro/basic_astro/modifiedEquinoctialElementConversions.h"
#include "tudat/astro/basic_astro/stateVectorIndices.h"
#include "tudat/astro/basic_astro/timeConversions.h"
//#include "Tudat/Astrodynamics/BasicAstrodynamics/UnitTests/testAccelerationModels.h"
//#include "Tudat/Astrodynamics/BasicAstrodynamics/UnitTests/testBody.h"
//#include "Tudat/Astrodynamics/BasicAstrodynamics/UnitTests/keplerPropagatorTestData.h"
#include "tudat/astro/basic_astro/astrodynamicsFunctions.h"
#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/basic_astro/bodyShapeModel.h"
#include "tudat/astro/basic_astro/oblateSpheroidBodyShapeModel.h"
#include "tudat/astro/basic_astro/sphericalBodyShapeModel.h"
#include "tudat/astro/basic_astro/massRateModel.h"
#include "tudat/astro/basic_astro/sphericalStateConversions.h"
#include "tudat/astro/basic_astro/unifiedStateModelQuaternionElementConversions.h"
#include "tudat/astro/basic_astro/unifiedStateModelModifiedRodriguesParameterElementConversions.h"
#include "tudat/astro/basic_astro/unifiedStateModelExponentialMapElementConversions.h"
#include "tudat/astro/basic_astro/stateRepresentationConversions.h"
#include "tudat/astro/basic_astro/empiricalAcceleration.h"
#include "tudat/astro/basic_astro/torqueModel.h"
#include "tudat/astro/basic_astro/torqueModelTypes.h"
#include "tudat/astro/basic_astro/dissipativeTorqueModel.h"

namespace py = pybind11;

namespace tudatpy {
    void expose_basic_astrodynamics(py::module &m);
}

#endif //TUDATBUNDLE_EXPOSE_BASIC_ASTRODYNAMICS_H
