//
// Created by ggarrett on 29-04-20.
//

#ifndef TUDATBUNDLE_EXPOSE_INTERPOLATORS_H
#define TUDATBUNDLE_EXPOSE_INTERPOLATORS_H

#include <pybind11/pybind11.h>

#include "Tudat/Mathematics/Interpolators/cubicSplineInterpolator.h"
#include "Tudat/Mathematics/Interpolators/hermiteCubicSplineInterpolator.h"
#include "Tudat/Mathematics/Interpolators/linearInterpolator.h"
#include "Tudat/Mathematics/Interpolators/lagrangeInterpolator.h"
#include "Tudat/Mathematics/Interpolators/interpolator.h"
#include "Tudat/Mathematics/Interpolators/lookupScheme.h"
#include "Tudat/Mathematics/Interpolators/multiDimensionalInterpolator.h"
#include "Tudat/Mathematics/Interpolators/oneDimensionalInterpolator.h"
#include "Tudat/Mathematics/Interpolators/multiLinearInterpolator.h"
#include "Tudat/Mathematics/Interpolators/piecewiseConstantInterpolator.h"
#include "Tudat/Mathematics/Interpolators/jumpDataLinearInterpolator.h"
#include "Tudat/Mathematics/Interpolators/createInterpolator.h"

namespace py = pybind11;

namespace tudatpy {

    void expose_interpolators(py::module &m);

}


#endif //TUDATBUNDLE_EXPOSE_INTERPOLATORS_H
