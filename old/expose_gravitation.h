//
// Created by ggarrett on 29-04-20.
//

#ifndef TUDATBUNDLE_EXPOSE_GRAVITATION_H
#define TUDATBUNDLE_EXPOSE_GRAVITATION_H

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>

#include "tudat/astro/gravitation/basicSolidBodyTideGravityFieldVariations.h"
#include "tudat/astro/gravitation/gravityFieldVariations.h"
#include "tudat/astro/gravitation/centralGravityModel.h"
#include "tudat/astro/gravitation/centralJ2GravityModel.h"
#include "tudat/astro/gravitation/centralJ2J3GravityModel.h"
#include "tudat/astro/gravitation/centralJ2J3J4GravityModel.h"
#include "tudat/astro/gravitation/gravityFieldModel.h"
#include "tudat/astro/gravitation/jacobiEnergy.h"
#include "tudat/astro/gravitation/librationPoint.h"
#include "tudat/astro/gravitation/sphericalHarmonicsGravityModel.h"
#include "tudat/astro/gravitation/sphericalHarmonicsGravityModelBase.h"
#include "tudat/astro/gravitation/sphericalHarmonicsGravityField.h"
#include "tudat/astro/gravitation/thirdBodyPerturbation.h"
#include "tudat/astro/gravitation/timeDependentSphericalHarmonicsGravityField.h"
#include "tudat/astro/gravitation/unitConversionsCircularRestrictedThreeBodyProblem.h"
#include "tudat/astro/gravitation/tests/planetTestData.h"   // <---- ?
#include "tudat/astro/gravitation/triAxialEllipsoidGravity.h"
#include "tudat/astro/gravitation/tabulatedGravityFieldVariations.h"
#include "tudat/astro/gravitation/mutualSphericalHarmonicGravityModel.h"
#include "tudat/astro/gravitation/secondDegreeGravitationalTorque.h"
#include "tudat/astro/gravitation/directTidalDissipationAcceleration.h"
#include "tudat/astro/gravitation/sphericalHarmonicGravitationalTorque.h"

namespace py = pybind11;

namespace tudatpy{

    void expose_gravitation(py::module &m);

}



#endif //TUDATBUNDLE_EXPOSE_GRAVITATION_H
