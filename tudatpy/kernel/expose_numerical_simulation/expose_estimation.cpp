/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_estimation.h"

#include "tudat/simulation/estimation_setup/fitOrbitToEphemeris.h"
#include "tudat/astro/propagators/propagateCovariance.h"
#include "tudat/basics/utilities.h"

#include "docstrings.h"
#include "scalarTypes.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>

namespace py = pybind11;
namespace tep = tudat::estimatable_parameters;
namespace tp = tudat::propagators;
namespace ts = tudat::statistics;
namespace tss = tudat::simulation_setup;
namespace tni = tudat::numerical_integrators;
namespace tom = tudat::observation_models;
namespace trf = tudat::reference_frames;
namespace te = tudat::ephemerides;



namespace tudatpy {
namespace numerical_simulation {
namespace estimation {



void expose_estimation(py::module &m) {



    py::class_<tom::ObservationViabilityCalculator,
            std::shared_ptr<tom::ObservationViabilityCalculator>>(m, "ObservationViabilityCalculator",
                                                                  get_docstring("ObservationViabilityCalculator").c_str() )
            .def("is_observation_viable", &tom::ObservationViabilityCalculator::isObservationViable,
                 py::arg( "link_end_states" ),
                 py::arg( "link_end_times" ),
                 get_docstring("ObservationViabilityCalculator.is_observation_viable").c_str() );

    py::class_<tom::ObservationSimulatorBase<STATE_SCALAR_TYPE, TIME_TYPE>,
            std::shared_ptr<tom::ObservationSimulatorBase<STATE_SCALAR_TYPE, TIME_TYPE>>>(m, "ObservationSimulator",
                                                                           get_docstring("ObservationSimulator").c_str() );

    py::class_<tom::ObservationSimulator<1,STATE_SCALAR_TYPE, TIME_TYPE>,
            std::shared_ptr<tom::ObservationSimulator<1,STATE_SCALAR_TYPE, TIME_TYPE>>,
            tom::ObservationSimulatorBase<STATE_SCALAR_TYPE, TIME_TYPE>>(m, "ObservationSimulator_1",
                                                          get_docstring("ObservationSimulator_1").c_str() );

    py::class_<tom::ObservationSimulator<2,STATE_SCALAR_TYPE, TIME_TYPE>,
            std::shared_ptr<tom::ObservationSimulator<2,STATE_SCALAR_TYPE, TIME_TYPE>>,
            tom::ObservationSimulatorBase<STATE_SCALAR_TYPE, TIME_TYPE>>(m, "ObservationSimulator_2",
                                                          get_docstring("ObservationSimulator_2").c_str() );

    py::class_<tom::ObservationSimulator<3,STATE_SCALAR_TYPE, TIME_TYPE>,
            std::shared_ptr<tom::ObservationSimulator<3,STATE_SCALAR_TYPE, TIME_TYPE>>,
            tom::ObservationSimulatorBase<STATE_SCALAR_TYPE, TIME_TYPE>>(m, "ObservationSimulator_3",
                                                          get_docstring("ObservationSimulator_3").c_str() );

    py::class_<tom::ObservationSimulator<6,STATE_SCALAR_TYPE, TIME_TYPE>,
            std::shared_ptr<tom::ObservationSimulator<6,STATE_SCALAR_TYPE, TIME_TYPE>>,
            tom::ObservationSimulatorBase<STATE_SCALAR_TYPE, TIME_TYPE>>(m, "ObservationSimulator_6",
                                                          get_docstring("ObservationSimulator_6").c_str() );

    m.def("simulate_observations",
          &tss::simulateObservations<STATE_SCALAR_TYPE, TIME_TYPE>,
          py::arg("simulation_settings"),
          py::arg("observation_simulators" ),
          py::arg("bodies"),
          get_docstring("simulate_observations").c_str() );

    m.def("compute_residuals_and_dependent_variables",
          &tss::computeResidualsAndDependentVariables<STATE_SCALAR_TYPE, TIME_TYPE>,
          py::arg("observation_collection"),
          py::arg("observation_simulators" ),
          py::arg("bodies"),
          get_docstring("compute_and_set_residuals").c_str() );


    m.def("create_pseudo_observations_and_models",
          &tss::simulatePseudoObservations<TIME_TYPE, STATE_SCALAR_TYPE>,
          py::arg("bodies"),
          py::arg("observed_bodies" ),
          py::arg("central_bodies" ),
          py::arg("initial_time"),
          py::arg("final_time"),
          py::arg("time_step"),
          get_docstring("create_pseudo_observations_and_models").c_str() );

    m.def("create_best_fit_to_ephemeris",
          &tss::createBestFitToCurrentEphemeris<TIME_TYPE, STATE_SCALAR_TYPE>,
          py::arg("bodies"),
          py::arg("acceleration_models"),
          py::arg("observed_bodies" ),
          py::arg("central_bodies" ),
          py::arg("integrator_settings" ),
          py::arg("initial_time"),
          py::arg("final_time"),
          py::arg("data_point_interval"),
          py::arg("additional_parameter_names") = std::vector< std::shared_ptr< tep::EstimatableParameterSettings > >( ),
          py::arg("number_of_iterations") = 3,
          py::arg("reintegrate_variational_equations") = true,
          py::arg("results_print_frequency") = 0.0,
          get_docstring("create_best_fit_to_ephemeris").c_str() );

    m.def("set_existing_observations",
          &tss::setExistingObservations<STATE_SCALAR_TYPE, TIME_TYPE>,
          py::arg("observations"),
          py::arg("reference_link_end" ),
          py::arg("ancilliary_settings_per_observatble" ) = std::map< tom::ObservableType,
          std::shared_ptr< tom::ObservationAncilliarySimulationSettings > >( ) );

    m.def("compute_target_angles_and_range",
          &tss::getTargetAnglesAndRange,
          py::arg("bodies"),
          py::arg("station_id" ),
          py::arg("target_body" ),
          py::arg("observation_times"),
          py::arg("is_station_transmitting"),
          get_docstring("compute_target_angles_and_range").c_str() );

    
    m.def("create_filtered_observation_collection",
          py::overload_cast<
                  const std::shared_ptr< tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > >,
                  const std::map< std::shared_ptr< tom::ObservationCollectionParser >, std::shared_ptr< tom::ObservationFilterBase > >& >(
          &tom::filterObservations< STATE_SCALAR_TYPE, TIME_TYPE > ),
          py::arg("original_observation_collection"),
          py::arg("observation_filters_map"),
          get_docstring("create_filtered_observation_collection").c_str() );

    m.def("create_filtered_observation_collection",
          py::overload_cast<
                  const std::shared_ptr< tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > >,
                  const std::shared_ptr< tom::ObservationFilterBase >,
                  const std::shared_ptr< tom::ObservationCollectionParser > >(
          &tom::filterObservations< STATE_SCALAR_TYPE, TIME_TYPE > ),
          py::arg("original_observation_collection"),
          py::arg("observation_filter"),
          py::arg("observation_parser") = std::make_shared< tom::ObservationCollectionParser >( ),
          get_docstring("create_filtered_observation_collection").c_str() );

    m.def("split_observation_collection",
          &tom::splitObservationSets< STATE_SCALAR_TYPE, TIME_TYPE >,
          py::arg("original_observation_collection"),
          py::arg("observation_set_splitter"),
          py::arg("observation_parser") = std::make_shared< tom::ObservationCollectionParser >( ),
          get_docstring("split_observation_collection").c_str() );

    m.def("create_new_observation_collection",
         &tom::createNewObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >,
         py::arg("original_observation_collection"),
         py::arg("observation_parser") = std::make_shared< tom::ObservationCollectionParser >( ),
         get_docstring("create_new_observation_collection ").c_str() );





}

}
}
}// namespace tudatpy
