//
// Created by ggarrett on 29-04-20.
//

#ifndef TUDATBUNDLE_EXPOSE_BASIC_ASTRODYNAMICS_H
#define TUDATBUNDLE_EXPOSE_BASIC_ASTRODYNAMICS_H

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModelTypes.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/attitudeElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/celestialBodyConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/convertMeanToEccentricAnomalies.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/clohessyWiltshirePropagator.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/customTorque.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/geodeticCoordinateConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/keplerPropagator.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/missionGeometry.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/modifiedEquinoctialElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
//#include "Tudat/Astrodynamics/BasicAstrodynamics/UnitTests/testAccelerationModels.h"
//#include "Tudat/Astrodynamics/BasicAstrodynamics/UnitTests/testBody.h"
//#include "Tudat/Astrodynamics/BasicAstrodynamics/UnitTests/keplerPropagatorTestData.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/astrodynamicsFunctions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/bodyShapeModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/oblateSpheroidBodyShapeModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/sphericalBodyShapeModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/massRateModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/sphericalStateConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unifiedStateModelQuaternionElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unifiedStateModelModifiedRodriguesParameterElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unifiedStateModelExponentialMapElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/stateRepresentationConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/empiricalAcceleration.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/torqueModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/torqueModelTypes.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/dissipativeTorqueModel.h"

namespace py = pybind11;

namespace tudatpy {
    void expose_basic_astrodynamics(py::module &m);
}

#endif //TUDATBUNDLE_EXPOSE_BASIC_ASTRODYNAMICS_H
