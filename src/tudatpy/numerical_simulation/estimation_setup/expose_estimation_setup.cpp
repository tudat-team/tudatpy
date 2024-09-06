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
#include <tudat/simulation/estimation_setup.h>

#include "tudatpy/scalarTypes.h"

namespace py = pybind11;
namespace tss = tudat::simulation_setup;
namespace tp = tudat::propagators;
namespace tom = tudat::observation_models;
namespace tep = tudat::estimatable_parameters;

namespace tudatpy {
    namespace numerical_simulation {
        namespace estimation_setup {


            PYBIND11_MODULE(expose_estimation_setup, m) {
                // *************** PARAMETER ***************

                m.def("print_parameter_names",
                      &tep::printEstimatableParameterEntries<double>,
                      py::arg("parameter_set"));


                // # EstimatableParameterSettings --> EstimatableParameterSet #
                m.def(
                    "create_parameter_set",
                    &tss::createParametersToEstimate<double>,
                    py::arg("parameter_settings"), py::arg("bodies"),
                    py::arg("propagator_settings") = nullptr,
                    py::arg("consider_parameters_names") = std::vector<
                        std::shared_ptr<tep::EstimatableParameterSettings> >(),
R"doc(Function for creating a consolidated set of estimatable parameters.

	Function for creating a consolidated parameter from the given estimatable parameter settings.
	The function checks for consistency between the parameter settings and the models contained in the simulation setup (given by the `bodies` & and `propagator_settings` parameters).


	:param parameter_settings:
		List of objects that define the settings for the parameters that are to be created. Each entry in this list is typically created by a call to a factory function in the :ref:`\`\`parameter\`\`` module

	:param bodies:
		Object consolidating all bodies and environment models, including ground station models, that constitute the physical environment.

	:param propagator_settings:
		Object containing the consolidated propagation settings of the simulation.

	:return:
		Instance of :class:`~tudatpy.numerical_simulation.estimation.EstimatableParameterSet` class, consolidating all estimatable parameters and simulation models.

)doc");


                // ************** OBSERVATION ***************


                // #   Observation Model Settings --> Observation Simulator #
                m.def("create_observation_simulators",
                      py::overload_cast<const std::vector<std::shared_ptr<
                                            tom::ObservationModelSettings> >&,
                                        const tss::SystemOfBodies&>(
                          &tom::createObservationSimulators<double, TIME_TYPE>),
                      py::arg("observation_settings"), py::arg("bodies"),
R"doc(Function for creating observation simulator objects.

	Factory function for creating observation simulator objects from observation settings.
	Note that each observation (i.e. combination of observable and link geometry) requires its own observation simulator object.


	:param observation_settings:
		List of settings objects, each object defining the observation model settings for one combination of observable and link geometry that is to be simulated.

	:param bodies:
		Object consolidating all bodies and environment models, including ground station models, that constitute the physical environment.

	:return:
		List of :class:`~tudatpy.numerical_simulation.estimation.ObservationSimulator` objects, each object hosting the functionality for simulating one combination of observable type and link geometry.
)doc");


                m.def(
                    "single_type_observation_collection",
                    py::overload_cast<
                        const tom::ObservableType, const tom::LinkDefinition&,
                        const std::vector<
                            Eigen::Matrix<double, Eigen::Dynamic, 1> >&,
                        const std::vector<TIME_TYPE>, const tom::LinkEndType,
                        const std::shared_ptr<
                            tom::ObservationAncilliarySimulationSettings> >(
                        &tom::createManualObservationCollection<double,
                                                                TIME_TYPE>),
                    py::arg("observable_type"), py::arg("link_ends"),
                    py::arg("observations_list"), py::arg("times_list"),
                    py::arg("reference_link_end"),
                    py::arg("ancilliary_settings") = nullptr,
"");
            }

        }  // namespace estimation_setup
    }  // namespace numerical_simulation
}  // namespace tudatpy
