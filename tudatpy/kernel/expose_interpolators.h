//
// Created by ggarrett on 29-04-20.
//

#ifndef TUDAT_EXPOSE_INTERPOLATORS_H
#define TUDAT_EXPOSE_INTERPOLATORS_H

#include <pybind11/pybind11.h>

#include "tudat/math/interpolators/cubicSplineInterpolator.h"
#include "tudat/math/interpolators/hermiteCubicSplineInterpolator.h"
#include "tudat/math/interpolators/linearInterpolator.h"
#include "tudat/math/interpolators/lagrangeInterpolator.h"
#include "tudat/math/interpolators/interpolator.h"
#include "tudat/math/interpolators/lookupScheme.h"
#include "tudat/math/interpolators/multiDimensionalInterpolator.h"
#include "tudat/math/interpolators/oneDimensionalInterpolator.h"
#include "tudat/math/interpolators/multiLinearInterpolator.h"
#include "tudat/math/interpolators/piecewiseConstantInterpolator.h"
#include "tudat/math/interpolators/jumpDataLinearInterpolator.h"
#include "tudat/math/interpolators/createInterpolator.h"

namespace py = pybind11;

namespace tudatpy {

    void expose_interpolators(py::module &m);

}


#endif //TUDATBUNDLE_EXPOSE_INTERPOLATORS_H
