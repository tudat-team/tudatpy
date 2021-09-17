/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#include <tudat/astro/aerodynamics/aerodynamicGuidance.h>
#include <tudat/astro/basic_astro.h>
#include <tudat/astro/propagators.h>


#include "expose_propagation.h"

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>


namespace py = pybind11;

namespace ta = tudat::aerodynamics;
namespace tp = tudat::propagators;
namespace tba = tudat::basic_astrodynamics;

namespace tudat
{

namespace aerodynamics
{

class PyAerodynamicGuidance : public ta::AerodynamicGuidance {
public:
    /* Inherit the constructors */
    using AerodynamicGuidance::AerodynamicGuidance;
    using AerodynamicGuidance::currentAngleOfAttack_;
    using AerodynamicGuidance::currentAngleOfSideslip_;
    using AerodynamicGuidance::currentBankAngle_;

    void updateGuidance( const double currentTime ) override {
        PYBIND11_OVERLOAD_PURE(void, AerodynamicGuidance, updateGuidance, currentTime ); }
};

}

}


namespace tudatpy {
namespace simulation {
namespace propagation {


void expose_propagation(py::module &m) {



    py::class_<ta::AerodynamicGuidance, ta::PyAerodynamicGuidance,
            std::shared_ptr< ta::AerodynamicGuidance > >(m, "AerodynamicGuidance")
            .def(py::init<>())
            .def("updateGuidance", &ta::AerodynamicGuidance::updateGuidance, py::arg("current_time") )
            .def_readwrite("angle_of_attack", &ta::PyAerodynamicGuidance::currentAngleOfAttack_)
            .def_readwrite("bank_angle", &ta::PyAerodynamicGuidance::currentBankAngle_)
            .def_readwrite("sideslip_angle", &ta::PyAerodynamicGuidance::currentAngleOfSideslip_);

    py::class_<
        tba::AccelerationModel<Eigen::Vector3d>,
        std::shared_ptr<tba::AccelerationModel<Eigen::Vector3d>>>(m, "AccelerationModel");

    py::class_<
        tba::TorqueModel,
        std::shared_ptr<tba::TorqueModel> >
            (m, "TorqueModel");


    m.def("get_single_integration_size",
          &tp::getSingleIntegrationSize,
          py::arg("state_type"));

    m.def("get_single_integration_differential_equation_order",
          &tp::getSingleIntegrationDifferentialEquationOrder,
          py::arg("state_type"));

    m.def("get_generalized_acceleration_size",
          &tp::getGeneralizedAccelerationSize,
          py::arg("state_type"));


    py::class_<
        tp::NBodyStateDerivative<double, double>,
        std::shared_ptr<tp::NBodyStateDerivative<double, double>>>
        NBodyStateDerivative_(m, "NBodyStateDerivative");

    py::class_<
        tp::NBodyCowellStateDerivative<double, double>,
        std::shared_ptr<tp::NBodyCowellStateDerivative<double, double>>>(m, "NBodyCowellStateDerivative")
        .def(py::init<
                 const tudat::basic_astrodynamics::AccelerationMap &,
                 const std::shared_ptr<tp::CentralBodyData<double, double>>,
                 const std::vector<std::string> &>(),
             py::arg("acceleration_models_per_body"),
             py::arg("central_body_data"),
             py::arg("bodies_to_integrate"));


}
}// namespace propagation
}// namespace simulation
}// namespace tudatpy
