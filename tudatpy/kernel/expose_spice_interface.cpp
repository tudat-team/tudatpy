//
// Created by ggarrett on 29-04-20.
//

#include "expose_spice_interface.h"
#include "docstrings.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace tsi = tudat::spice_interface;

namespace tudatpy {

    void expose_spice_interface(py::module &m) {

        m.def("load_standard_spice_kernels",
              py::overload_cast<
                      const std::vector<std::string>,
                      const std::string
                      >(&tudat::spice_interface::loadStandardSpiceKernels),
              py::arg("alternative_kernels") = std::vector<std::string>(),
#ifdef TUDAT_INSTALL
              py::arg("state") = "install",
#else
              py::arg("state") = "build",
#endif
              tudatpy::load_standard_spice_kernels_docstring().c_str());

        m.def("load_standard_spice_kernels",
              py::overload_cast<
                      const std::vector<std::string>
                      >(&tudat::spice_interface::loadStandardSpiceKernels),
              py::arg("alternative_kernels") = std::vector<std::string>(),
              tudatpy::load_standard_spice_kernels_docstring().c_str());

        m.def("clear_spice_kernels",
              &tudat::spice_interface::clearSpiceKernels,
              tudatpy::clear_spice_kernels_docstring().c_str());

        m.def("get_average_radius",
              &tudat::spice_interface::getAverageRadius,
              "<no_doc>");


    };
}
