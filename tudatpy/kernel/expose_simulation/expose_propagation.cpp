/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudatpy/docstrings.h"

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
namespace tss = tudat::simulation_setup;
namespace tni = tudat::numerical_integrators;

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
    // FREE FUNCTIONS

    // NOTE: the following 5 functions do not strictly belong to the actual "propagator settings".

    m.def("get_initial_state_of_bodies",
          py::overload_cast<const std::vector<std::string> &,
                  const std::vector<std::string> &,
                  const tss::SystemOfBodies &,
                  const double>(
                  &tp::getInitialStatesOfBodies<>),
          py::arg("bodies_to_propagate"),
          py::arg("central_bodies"),
          py::arg("body_system"),
          py::arg("initial_time"),
          get_docstring("get_initial_state_of_bodies").c_str());

    m.def("get_initial_state_of_body",// overload [2/2]
          py::overload_cast<const std::string&,
          const std::string&,
          const tss::SystemOfBodies&,
          const double>(
              &tp::getInitialStateOfBody<>),
          py::arg("body_to_propagate"),
          py::arg("central_body"),
          py::arg("bodies"),
          py::arg("initial_time"));

    m.def("get_initial_rotational_state_of_body",
          py::overload_cast<const std::string&,
          const std::string&,
          const tss::SystemOfBodies&,
          const double>(
              &tp::getInitialRotationalStateOfBody<>),
          py::arg("body_to_propagate"),
          py::arg("base_orientation"),
          py::arg("bodies"),
          py::arg("initial_time"));

    m.def("get_zero_proper_mode_rotational_state",
          py::overload_cast<
          const tss::SystemOfBodies&,
          const std::shared_ptr< tni::IntegratorSettings< double > >,
          const std::shared_ptr< tp::SingleArcPropagatorSettings< double > >,
          const double,
          const std::vector< double >,
          const bool >( &tp::getZeroProperModeRotationalState< > ),
          py::arg("bodies"),
          py::arg("integrator_settings"),
          py::arg("propagator_settings"),
          py::arg("body_mean_rotational_rate"),
          py::arg("dissipation_times"),
          py::arg("propagate_undamped") = true );

    m.def("combine_initial_states",
          &tp::createCombinedInitialState<double>,
          py::arg("propagator_settings_per_type"),
          get_docstring("combine_initial_states").c_str());

    py::class_<
        tba::AccelerationModel<Eigen::Vector3d>,
        std::shared_ptr<tba::AccelerationModel<Eigen::Vector3d>>>(m, "AccelerationModel");

    // First overload
    m.def("create_acceleration_models",
          py::overload_cast<const tss::SystemOfBodies &,
                  const tss::SelectedAccelerationMap &,
                  const std::map<std::string, std::string> &>(
                  &tss::createAccelerationModelsMap),
          py::arg("body_system"),
          py::arg("selected_acceleration_per_body"),
          py::arg("central_bodies"),
          get_docstring("create_acceleration_models", 0).c_str());

    // Second overload
    m.def("create_acceleration_models",
          py::overload_cast<const tss::SystemOfBodies &,
                  const tss::SelectedAccelerationMap &,
                  const std::vector<std::string> &,
                  const std::vector<std::string> &>(
                  &tss::createAccelerationModelsMap),
          py::arg("body_system"),
          py::arg("selected_acceleration_per_body"),
          py::arg("bodies_to_propagate"),
          py::arg("central_bodies"),
          get_docstring("create_acceleration_models", 1).c_str());

    m.def("create_torque_models",
          &tss::createTorqueModelsMap,
          py::arg("body_system"),
          py::arg("selected_torque_per_body"),
          py::arg("bodies_to_propagate"),
          get_docstring("create_torque_models").c_str());


}
}// namespace propagation
}// namespace simulation
}// namespace tudatpy
