//
// Created by ggarrett on 29-04-20.
//

#include "expose_spice_interface.h"
#include "docstrings.h"

#include <pybind11/pybind11.h>

namespace py = pybind11;
namespace tsi = tudat::spice_interface;

namespace tudatpy {

    void expose_spice_interface(py::module &m) {

        m.def("load_standard_kernels",
              &tudat::spice_interface::loadStandardSpiceKernels,
              py::arg("alternative_kernels") = std::vector<std::string>(),
              tudatpy::load_standard_kernels_docstring().c_str());

        m.def("load_kernel",
              &tudat::spice_interface::loadSpiceKernelInTudat,
              "<no_doc>");

        m.def("clear_kernels",
              &tudat::spice_interface::clearSpiceKernels,
              tudatpy::clear_spice_kernels_docstring().c_str());

        m.def("get_average_radius",
              &tudat::spice_interface::getAverageRadius,
              "<no_doc>");

        m.def("get_body_gravitational_parameter",
			  &tudat::spice_interface::getBodyGravitationalParameter,
			  "<no_doc>");

        m.def("get_body_cartesian_state_at_epoch",
			  &tudat::spice_interface::getBodyCartesianStateAtEpoch,
			  py::arg("target_body_name"),
			  py::arg("observer_body_name"),
			  py::arg("reference_frame_name"),
			  py::arg("aberration_corrections"),
			  py::arg("ephemeris_time"));

    };
}
