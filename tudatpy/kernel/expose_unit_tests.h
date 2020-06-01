//
// Created by ggarrett on 02-05-20.
//

#ifndef TUDATBUNDLE_EXPOSE_UNIT_TESTS_H
#define TUDATBUNDLE_EXPOSE_UNIT_TESTS_H

#include <pybind11/pybind11.h>
//#include "Tudat/Astrodynamics/Aerodynamics/UnitTests/testApolloCapsuleCoefficients.h"
#include "tudat/astro/aerodynamics/hypersonicLocalInclinationAnalysis.h"

namespace py = pybind11;

namespace tudatpy {
    void expose_unit_tests(py::module &m);

}

#endif //TUDATBUNDLE_EXPOSE_UNIT_TESTS_H
