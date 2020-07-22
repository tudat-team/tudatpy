//
// Created by elmar on 13-07-20.
//

#ifndef TUDATPY_EXPOSE_LOW_THRUST_H
#define TUDATPY_EXPOSE_LOW_THRUST_H

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>

#include <tudat/astro/low_thrust/lowThrustLeg.h>
#include <tudat/astro/low_thrust/lowThrustLegSettings.h>
#include <tudat/astro/low_thrust/shape_based/sphericalShaping.h>
#include <tudat/math/root_finders/createRootFinder.h>

namespace py = pybind11;

namespace tudatpy {

	void expose_low_thrust(py::module &m);

};

#endif //TUDATPY_EXPOSE_LOW_THRUST_H
