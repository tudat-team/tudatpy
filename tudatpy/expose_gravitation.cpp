//
// Created by ggarrett on 29-04-20.
//

#include "expose_gravitation.h"

namespace py = pybind11;
namespace tg = tudat::gravitation;

namespace tudatpy {

    void expose_gravitation(py::module &m){

        py::class_ < tg::GravityFieldModel,
                std::shared_ptr < tg::GravityFieldModel >> (m, "GravityFieldModel")
                .def(py::init<
                             const double,
                             const std::function<void()>
                     >(),
                     py::arg("gravitational_parameter"),
                     py::arg("update_inertia_tensor") = std::function<void()>())
                .def("get_gravitational_parameter", &tg::GravityFieldModel::getGravitationalParameter)
                .def_property("gravitational_parameter", &tg::GravityFieldModel::getGravitationalParameter,
                              &tg::GravityFieldModel::resetGravitationalParameter);


    };

}
