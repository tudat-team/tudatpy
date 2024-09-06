/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <tudat/basics/deprecationWarnings.h>
#include <tudat/simulation/propagation_setup.h>


namespace py = pybind11;
namespace tba = tudat::basic_astrodynamics;
namespace tss = tudat::simulation_setup;
namespace tp = tudat::propagators;
namespace tinterp = tudat::interpolators;
namespace te = tudat::ephemerides;
namespace tni = tudat::numerical_integrators;
namespace trf = tudat::reference_frames;
namespace tmrf = tudat::root_finders;


namespace tudat {
    namespace simulation_setup {

        inline std::shared_ptr<MassRateModelSettings> customMassRateDeprecated(
            const std::function<double(const double)> massRateFunction)

        {
            static bool isWarningPrinted = false;
            if(isWarningPrinted == false) {
                tudat::utilities::printDeprecationWarning(
                    "tudatpy.numerical_simulation.propagation_setup.mass_rate."
                    "custom",
                    "tudatpy.numerical_simulation.propagation_setup.mass_rate."
                    "custom_mass_rate");
                isWarningPrinted = true;
            }

            return customMassRate(massRateFunction);
        }
    }  // namespace simulation_setup
}  // namespace tudat

namespace tudatpy {
    namespace numerical_simulation {
        namespace propagation_setup {
            namespace mass_rate {

                PYBIND11_MODULE(expose_mass_rate, m) {
                    // Enums

                    py::enum_<tba::AvailableMassRateModels>(
                        m, "AvailableMassRateModels",
R"doc(Enumeration of available mass rate models.

	Enumeration of mass rate models supported by tudat.


	:member undefined_mass_rate_type:
	:member custom_mass_rate_type:
	:member from_thrust_mass_rate_type:
)doc")
                        .value("undefined_mass_rate_type",
                               tba::AvailableMassRateModels::
                                   undefined_mass_rate_model,
"")
                        .value(
                            "custom_mass_rate_type",
                            tba::AvailableMassRateModels::
                                custom_mass_rate_model,
"")
                        .value("from_thrust_mass_rate_type",
                               tba::AvailableMassRateModels::
                                   from_thrust_mass_rate_model,
"")
                        .export_values();

                    // Classes

                    py::class_<tss::MassRateModelSettings,
                               std::shared_ptr<tss::MassRateModelSettings>>(
                        m, "MassRateModelSettings",
R"doc(Functional base class to define settings for mass rates.

	Base class for providing settings for a mass rate model, that defines the models to be used to numerically propagate the
	mass of a body during a simulation. If any additional information (besides the type of the mass rate model) is required,
	these must be implemented in a derived class (dedicated for each mass rate model type).

)doc");
                    //                .def(py::init<const
                    //                tudat::basic_astrodynamics::AvailableMassRateModels>(),
                    //                     py::arg("mass_rate_type"));

                    py::class_<tss::FromThrustMassRateSettings,
                               std::shared_ptr<tss::FromThrustMassRateSettings>,
                               tss::MassRateModelSettings>(
                        m, "FromThrustMassRateSettings",
R"doc(`MassRateModelSettings`-derived class to define settings for a mass rate model derived from a thrust model.

	`MassRateModelSettings`-derived class to define settings for a mass rate model derived from a thrust model.

)doc");
                    //                .def(py::init<const bool, const
                    //                std::string &>(),
                    //                     py::arg("use_all_thrust_models") = 1,
                    //                     py::arg("associated_thrust_source") =
                    //                     "");

                    py::class_<tss::CustomMassRateSettings,
                               std::shared_ptr<tss::CustomMassRateSettings>,
                               tss::MassRateModelSettings>(
                        m, "CustomMassRateSettings",
R"doc(`MassRateModelSettings`-derived class to define settings for a custom mass rate model.

	`MassRateModelSettings`-derived class to define settings for a custom mass rate model.

)doc");


                    // Factory functions

                    m.def("from_thrust", &tss::fromThrustMassRate,
                          py::arg("use_all_thrust_models") = 1,
                          py::arg("associated_thrust_source") = "",
R"doc(Creates the settings for a mass rate model defined from a thrust model.

	Creates the settings for a mass rate model defined from a thrust model. The mass rate model is derived from
	the associated body's engine model. It is possible to consider only a specific engine or all engines.


	:param use_all_thrust_models:
		Denotes whether all engines of the associated body are to be combined into a single thrust model.
	:param associated_thrust_source:
		Name of engine model from which thrust is to be derived (must be empty if the first argument is set to true).
	:return:
		From thrust mass rate settings object.
)doc");

                    m.def("custom", &tss::customMassRateDeprecated,
                          py::arg("mass_rate_function"));

                    m.def("custom_mass_rate", &tss::customMassRate,
                          py::arg("mass_rate_function"),
R"doc(Creates the settings for a mass rate model defined from a thrust model.

	Creates the settings for a custom mass rate model defined through a mass rate function. The function must have
	time as an independent variable.


	:param mass_rate_function:
		Function of time defining the custom mass rate.
	:return:
		Custom mass rate settings object.
)doc");
                }

            }  // namespace mass_rate
        }  // namespace propagation_setup
    }  // namespace numerical_simulation
}  // namespace tudatpy
