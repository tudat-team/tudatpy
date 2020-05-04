//
// Created by ggarrett on 02-05-20.
//

#ifndef TUDATBUNDLE_EXPOSE_AERODYNAMICS_H
#define TUDATBUNDLE_EXPOSE_AERODYNAMICS_H

#include <pybind11/pybind11.h>
#include "Tudat/Astrodynamics/Aerodynamics/hypersonicLocalInclinationAnalysis.h"
#include "Tudat/Astrodynamics/Aerodynamics/flightConditions.h"
#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicCoefficientInterface.h"

namespace py = pybind11;

namespace tudatpy {

    void expose_aerodynamics(py::module &m);

};


#endif //TUDATBUNDLE_EXPOSE_AERODYNAMICS_H
