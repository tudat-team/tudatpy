//
// Created by ggarrett on 29-04-20.
//

#ifndef TUDATBUNDLE_EXPOSE_GRAVITATION_H
#define TUDATBUNDLE_EXPOSE_GRAVITATION_H

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>

#include "Tudat/Astrodynamics/Gravitation/basicSolidBodyTideGravityFieldVariations.h"
#include "Tudat/Astrodynamics/Gravitation/gravityFieldVariations.h"
#include "Tudat/Astrodynamics/Gravitation/centralGravityModel.h"
#include "Tudat/Astrodynamics/Gravitation/centralJ2GravityModel.h"
#include "Tudat/Astrodynamics/Gravitation/centralJ2J3GravityModel.h"
#include "Tudat/Astrodynamics/Gravitation/centralJ2J3J4GravityModel.h"
#include "Tudat/Astrodynamics/Gravitation/gravityFieldModel.h"
#include "Tudat/Astrodynamics/Gravitation/jacobiEnergy.h"
#include "Tudat/Astrodynamics/Gravitation/librationPoint.h"
#include "Tudat/Astrodynamics/Gravitation/sphericalHarmonicsGravityModel.h"
#include "Tudat/Astrodynamics/Gravitation/sphericalHarmonicsGravityModelBase.h"
#include "Tudat/Astrodynamics/Gravitation/sphericalHarmonicsGravityField.h"
#include "Tudat/Astrodynamics/Gravitation/thirdBodyPerturbation.h"
#include "Tudat/Astrodynamics/Gravitation/timeDependentSphericalHarmonicsGravityField.h"
#include "Tudat/Astrodynamics/Gravitation/unitConversionsCircularRestrictedThreeBodyProblem.h"
#include "Tudat/Astrodynamics/Gravitation/UnitTests/planetTestData.h"
#include "Tudat/Astrodynamics/Gravitation/triAxialEllipsoidGravity.h"
#include "Tudat/Astrodynamics/Gravitation/tabulatedGravityFieldVariations.h"
#include "Tudat/Astrodynamics/Gravitation/mutualSphericalHarmonicGravityModel.h"
#include "Tudat/Astrodynamics/Gravitation/secondDegreeGravitationalTorque.h"
#include "Tudat/Astrodynamics/Gravitation/directTidalDissipationAcceleration.h"
#include "Tudat/Astrodynamics/Gravitation/sphericalHarmonicGravitationalTorque.h"

namespace py = pybind11;

namespace tudatpy{

    void expose_gravitation(py::module &m);

}



#endif //TUDATBUNDLE_EXPOSE_GRAVITATION_H
