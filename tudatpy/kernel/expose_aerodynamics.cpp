//
// Created by ggarrett on 02-05-20.
//

#include "expose_aerodynamics.h"


namespace py = pybind11;
namespace ta = tudat::aerodynamics;
namespace tudatpy {

    void expose_aerodynamics(py::module &m) {

        py::class_<ta::AerodynamicCoefficientInterface,
                std::shared_ptr<ta::AerodynamicCoefficientInterface>
        > AerodynamicCoefficientInterface_(m, "AerodynamicCoefficientInterface", "<no_doc, only_dec>");

        py::class_<ta::AerodynamicCoefficientGenerator<3, 6>,
                std::shared_ptr<ta::AerodynamicCoefficientGenerator<3, 6>>,
                ta::AerodynamicCoefficientInterface
        > AerodynamicCoefficientGenerator36_(m, "AerodynamicCoefficientGenerator36", "<no_doc, only_dec>");

        py::class_<ta::HypersonicLocalInclinationAnalysis,
                std::shared_ptr<ta::HypersonicLocalInclinationAnalysis>,
                ta::AerodynamicCoefficientGenerator<3, 6>
        > HypersonicLocalInclinationAnalysis_(m, "HypersonicLocalInclinationAnalysis");

        py::class_<ta::FlightConditions,
                std::shared_ptr<ta::FlightConditions>>(m, "FlightConditions")
                .def(py::init<
                             const std::shared_ptr<tudat::basic_astrodynamics::BodyShapeModel>,
                             const std::shared_ptr<tudat::reference_frames::AerodynamicAngleCalculator>
                     >(),
                     py::arg("shape_model"),
                     py::arg("aerodynamic_angle_calculator") = std::shared_ptr<tudat::reference_frames::AerodynamicAngleCalculator>()
                )
                .def("get_aerodynamic_angle_calculator", &ta::FlightConditions::getAerodynamicAngleCalculator);

    };

};
