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
#include <tudat/basics/deprecationWarnings.h>

#include "docstrings.h"
#include <tudat/simulation/propagation_setup.h>

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


namespace tudat
{
namespace simulation_setup
{

inline std::shared_ptr< MassRateModelSettings > customMassRateDeprecated(
        const std::function< double( const double ) > massRateFunction )

{
    static bool isWarningPrinted = false;
    if( isWarningPrinted == false )
    {
        tudat::utilities::printDeprecationWarning( "tudatpy.numerical_simulation.propagation_setup.mass_rate.custom",
                             "tudatpy.numerical_simulation.propagation_setup.mass_rate.custom_mass_rate");
        isWarningPrinted = true;
    }

    return customMassRate( massRateFunction );

}
}
}

namespace tudatpy {
namespace numerical_simulation {
namespace propagation_setup {
namespace mass_rate {

    void expose_mass_rate_setup(py::module &m) {

        // Enums

        py::enum_<tba::AvailableMassRateModels>(m, "AvailableMassRateModels",
                                                get_docstring("AvailableMassRateModels").c_str())
                .value("undefined_mass_rate_type", tba::AvailableMassRateModels::undefined_mass_rate_model,
                       get_docstring("AvailableMassRateModels.undefined_mass_rate_type").c_str())
                .value("custom_mass_rate_type", tba::AvailableMassRateModels::custom_mass_rate_model,
                       get_docstring("AvailableMassRateModels.custom_mass_rate_type").c_str())
                .value("from_thrust_mass_rate_type", tba::AvailableMassRateModels::from_thrust_mass_rate_model,
                       get_docstring("AvailableMassRateModels.from_thrust_mass_rate_type").c_str())
                .export_values();

        // Classes

        py::class_<tss::MassRateModelSettings,
                std::shared_ptr<tss::MassRateModelSettings>>(m, "MassRateModelSettings",
                        get_docstring("MassRateModelSettings").c_str());
//                .def(py::init<const tudat::basic_astrodynamics::AvailableMassRateModels>(),
//                     py::arg("mass_rate_type"));

        py::class_<tss::FromThrustMassRateSettings,
                std::shared_ptr<tss::FromThrustMassRateSettings>,
                tss::MassRateModelSettings>(m, "FromThrustMassRateSettings",
                                            get_docstring("FromThrustMassRateSettings").c_str());
//                .def(py::init<const bool, const std::string &>(),
//                     py::arg("use_all_thrust_models") = 1,
//                     py::arg("associated_thrust_source") = "");

        py::class_<tss::CustomMassRateSettings,
                std::shared_ptr<tss::CustomMassRateSettings>,
                tss::MassRateModelSettings>(m, "CustomMassRateSettings",
                                            get_docstring("CustomMassRateSettings").c_str());


        // Factory functions

        m.def("from_thrust", &tss::fromThrustMassRate,
              py::arg("use_all_thrust_models") = 1,
              py::arg("associated_thrust_source") = "",
              get_docstring("from_thrust").c_str());

        m.def("custom", &tss::customMassRateDeprecated,
              py::arg("mass_rate_function") );

        m.def("custom_mass_rate", &tss::customMassRate,
              py::arg("mass_rate_function"),
              get_docstring("custom_mass_rate").c_str());

    }

}// namespace mass
}// namespace propagation_setup
}// namespace numerical_simulation
}// namespace tudatpy
