//
// Created by ggarrett on 29-04-20.
//

#ifndef TUDATBUNDLE_EXPOSE_SPICE_INTERFACE_H
#define TUDATBUNDLE_EXPOSE_SPICE_INTERFACE_H

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "Tudat/External/SpiceInterface/spiceInterface.h"

namespace py = pybind11;
namespace tsi = tudat::spice_interface;

namespace tudatpy {

    void expose_spice_interface(py::module &m);
}


#endif //TUDATBUNDLE_EXPOSE_SPICE_INTERFACE_H
