//
// Created by ggarrett on 26-04-20.
//

#ifndef TUDATBUNDLE_EXPOSE_PROPAGATORS_H
#define TUDATBUNDLE_EXPOSE_PROPAGATORS_H

#include <pybind11/pybind11.h>

//#include "Tudat/Astrodynamics/Propagators/centralBodyData.h"
#include "tudat/astro/propagators/nBodyStateDerivative.h"
#include "tudat/astro/propagators/nBodyCowellStateDerivative.h"
//#include "Tudat/Astrodynamics/Propagators/nBodyEnckeStateDerivative.h"
//#include "Tudat/Astrodynamics/Propagators/nBodyGaussKeplerStateDerivative.h"
//#include "Tudat/Astrodynamics/Propagators/nBodyGaussModifiedEquinoctialStateDerivative.h"
//#include "Tudat/Astrodynamics/Propagators/nBodyUnifiedStateModelQuaternionsStateDerivative.h"
//#include "Tudat/Astrodynamics/Propagators/nBodyUnifiedStateModelModifiedRodriguesParametersStateDerivative.h"
//#include "Tudat/Astrodynamics/Propagators/nBodyUnifiedStateModelExponentialMapStateDerivative.h"
//#include "Tudat/Astrodynamics/Propagators/dynamicsStateDerivativeModel.h"
#include "tudat/astro/propagators/singleStateTypeDerivative.h"
//#include "Tudat/Astrodynamics/Propagators/integrateEquations.h"
//#include "Tudat/Astrodynamics/Propagators/bodyMassStateDerivative.h"
//#include "Tudat/Astrodynamics/Propagators/variationalEquations.h"
//#include "Tudat/Astrodynamics/Propagators/stateTransitionMatrixInterface.h"
//#include "Tudat/Astrodynamics/Propagators/environmentUpdateTypes.h"
//#include "Tudat/Astrodynamics/Propagators/customStateDerivative.h"
//#include "Tudat/Astrodynamics/Propagators/rotationalMotionStateDerivative.h"
//#include "Tudat/Astrodynamics/Propagators/rotationalMotionQuaternionsStateDerivative.h"
//#include "Tudat/Astrodynamics/Propagators/rotationalMotionModifiedRodriguesParametersStateDerivative.h"
//#include "Tudat/Astrodynamics/Propagators/rotationalMotionExponentialMapStateDerivative.h"
//#include "Tudat/Astrodynamics/Propagators/stateDerivativeCircularRestrictedThreeBodyProblem.h"
//#include "Tudat/Astrodynamics/Propagators/getZeroProperModeRotationalInitialState.h"
//#include "Tudat/Astrodynamics/Propagators/propagateCovariance.h"

#include "tudat/simulation/propagation/dynamicsSimulator.h"


namespace py = pybind11;
namespace tudatpy {
    void expose_propagators(py::module &m);

}

#endif //TUDATBUNDLE_EXPOSE_PROPAGATORS_H
