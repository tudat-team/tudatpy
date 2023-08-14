/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_estimation_setup.h"
#include "expose_estimation_setup/expose_estimated_parameter_setup.h"
#include "expose_estimation_setup/expose_observation_setup.h"

#include "tudatpy/docstrings.h"
#include "tudatpy/scalarTypes.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>

namespace py = pybind11;
namespace tss = tudat::simulation_setup;
namespace tp = tudat::propagators;
namespace tom = tudat::observation_models;
namespace tep = tudat::estimatable_parameters;

namespace tudatpy {
namespace numerical_simulation {
namespace estimation_setup {


void expose_estimation_setup(py::module &m) {

    // *************** PARAMETER ***************

    m.def("print_parameter_names",
          &tep::printEstimatableParameterEntries< double >,
          py::arg("parameter_set") );


    // # EstimatableParameterSettings --> EstimatableParameterSet #
    m.def("create_parameter_set",
          &tss::createParametersToEstimate< double >,
          py::arg("parameter_settings"),
          py::arg("bodies"),
          py::arg("propagator_settings") = nullptr,
          py::arg("consider_parameters_names") = std::vector< std::shared_ptr< tep::EstimatableParameterSettings > >( ),
          get_docstring("create_parameter_set").c_str() );


    // ************** OBSERVATION ***************



    // #   Observation Model Settings --> Observation Simulator #
    m.def("create_observation_simulators",
          py::overload_cast< const std::vector< std::shared_ptr< tom::ObservationModelSettings > >&, const tss::SystemOfBodies& >(
                  &tom::createObservationSimulators< double, TIME_TYPE > ),
          py::arg( "observation_settings" ),
          py::arg( "bodies" ),
          get_docstring("create_observation_simulators").c_str() );


    m.def("single_type_observation_collection",
          py::overload_cast<
            const tom::ObservableType,
            const tom::LinkDefinition&,
            const std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >&,
            const std::vector< TIME_TYPE >,
            const tom::LinkEndType,
            const std::shared_ptr< tom::ObservationAncilliarySimulationSettings < TIME_TYPE > > >(
        &tom::createManualObservationCollection< double, TIME_TYPE > ),
          py::arg("observable_type"),
          py::arg("link_ends"),
          py::arg("observations_list"),
          py::arg("times_list"),
          py::arg("reference_link_end" ),
          py::arg("ancilliary_settings" ) = nullptr,
          get_docstring("single_type_observaion_collection").c_str() );



    // ************** Modules ***************

    auto parameter_setup = m.def_submodule("parameter");
    parameter::expose_estimated_parameter_setup(parameter_setup);

    auto observation_setup = m.def_submodule("observation");
    observation::expose_observation_setup(observation_setup);
}

}
}
}// namespace tudatpy
