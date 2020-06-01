//
// Created by ggarrett on 02-05-20.
//

#ifndef TUDATBUNDLE_EXPOSE_AERODYNAMICS_H
#define TUDATBUNDLE_EXPOSE_AERODYNAMICS_H

#include <pybind11/pybind11.h>
#include "tudat/astro/aerodynamics/hypersonicLocalInclinationAnalysis.h"
#include "tudat/astro/aerodynamics/flightConditions.h"
#include "tudat/astro/aerodynamics/aerodynamicCoefficientInterface.h"

namespace py = pybind11;

namespace tudatpy {

    void expose_aerodynamics(py::module &m);

};


#endif //TUDATBUNDLE_EXPOSE_AERODYNAMICS_H
