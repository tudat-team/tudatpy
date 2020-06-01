//
// Created by ggarrett on 02-05-20.
//

#include "expose_reference_frames.h"


namespace trf = tudat::reference_frames;

namespace py = pybind11;
namespace tudatpy {

    void expose_reference_frames(py::module &m) {

        py::class_<trf::AerodynamicAngleCalculator,
                std::shared_ptr<trf::AerodynamicAngleCalculator>>(m, "AerodynamicAngleCalculator")
                .def("set_orientation_angle_functions",
                     py::overload_cast<
                             const std::function<double()>,
                             const std::function<double()>,
                             const std::function<double()>,
                             const std::function<void(const double)>
                     >(&trf::AerodynamicAngleCalculator::setOrientationAngleFunctions),
                     py::arg("angle_of_attack_function") = std::function<double()>(),
                     py::arg("angle_of_sideslip_function") = std::function<double()>(),
                     py::arg("bank_angle_function") = std::function<double()>(),
                     py::arg("angle_update_function") = std::function<void(const double)>(),
                     "<no_doc>"
                )
                .def("set_orientation_angle_functions",
                     py::overload_cast<
                             const double,
                             const double,
                             const double
                     >(&trf::AerodynamicAngleCalculator::setOrientationAngleFunctions),
                     py::arg("angle_of_attack") = TUDAT_NAN,
                     py::arg("angle_of_sideslip") = TUDAT_NAN,
                     py::arg("bank_angle") = TUDAT_NAN,
                     "<no_doc>"
                );


    }
}