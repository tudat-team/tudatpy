/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_mass_rate_setup.h"

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;
namespace tba = tudat::basic_astrodynamics;
namespace tss = tudat::simulation_setup;
namespace tp = tudat::propagators;
namespace tinterp = tudat::interpolators;
namespace te = tudat::ephemerides;
namespace tni = tudat::numerical_integrators;
namespace trf = tudat::reference_frames;
namespace tmrf = tudat::root_finders;

namespace tudatpy {
namespace simulation {
namespace propagation_setup {

    void expose_mass_rate_setup(py::module &m) {

        // Enums

        py::enum_<tba::AvailableMassRateModels>(m, "AvailableMassRateModels")
                .value("undefined_mass_rate_type", tba::AvailableMassRateModels::undefined_mass_rate_model)
                .value("custom_mass_rate_type", tba::AvailableMassRateModels::custom_mass_rate_model)
                .value("from_thrust_mass_rate_type", tba::AvailableMassRateModels::from_thrust_mass_rate_model)
                .export_values();

        // Classes

        py::class_<tss::MassRateModelSettings,
                std::shared_ptr<tss::MassRateModelSettings>>(m, "MassRateModelSettings")
//                .def(py::init<const tudat::basic_astrodynamics::AvailableMassRateModels>(),
//                     py::arg("mass_rate_type"));

        py::class_<tss::FromThrustMassModelSettings,
                std::shared_ptr<tss::FromThrustMassModelSettings>,
                tss::MassRateModelSettings>(m, "FromThrustMassModelSettings")
//                .def(py::init<const bool, const std::string &>(),
//                     py::arg("use_all_thrust_models") = 1,
//                     py::arg("associated_thrust_source") = "");

        py::class_<tss::CustomMassRateModelSettings,
                std::shared_ptr<tss::CustomMassRateModelSettings>,
                tss::CustomMassRateModelSettings>(m, "CustomMassRateModelSettings")

        // Factory functions

        m.def("from_thrust", &tss::fromThrustMassRate,
              py::arg("use_all_thrust_models") = 1,
              py::arg("associated_thrust_source") = "");

        m.def("custom", &tss::customMassRate,
              py::arg("mass_rate_function"));

    }

}// namespace propagation_setup
}// namespace simulation
}// namespace tudatpy