/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_observations.h"

#include <tudat/astro/observation_models.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace tom = tudat::observation_models;

namespace tudatpy {

void expose_observations(py::module &m) {


    py::class_<tom::ObservationViabilityCalculator,
               std::shared_ptr<tom::ObservationViabilityCalculator>>(m, "ObservationViabilityCalculator")
            .def("is_observation_viable", &tom::ObservationViabilityCalculator::isObservationViable,
                 py::arg( "link_end_states" ),
                 py::arg( "link_end_times" ) );

    py::class_<tom::ObservationSimulatorBase<double,double>,
               std::shared_ptr<tom::ObservationSimulatorBase<double,double>>>(m, "ObservationSimulator");

    py::class_<tom::ObservationSimulator<1,double,double>,
               std::shared_ptr<tom::ObservationSimulator<1,double,double>>,
               tom::ObservationSimulatorBase<double,double>>(m, "ObservationSimulator_1");

    py::class_<tom::ObservationSimulator<2,double,double>,
               std::shared_ptr<tom::ObservationSimulator<2,double,double>>,
               tom::ObservationSimulatorBase<double,double>>(m, "ObservationSimulator_2");

    py::class_<tom::ObservationSimulator<3,double,double>,
               std::shared_ptr<tom::ObservationSimulator<3,double,double>>,
               tom::ObservationSimulatorBase<double,double>>(m, "ObservationSimulator_3");

    py::class_<tom::ObservationSimulator<6,double,double>,
               std::shared_ptr<tom::ObservationSimulator<6,double,double>>,
               tom::ObservationSimulatorBase<double,double>>(m, "ObservationSimulator_6");

    py::class_<tom::ObservationSimulationTimeSettings<double>,
               std::shared_ptr<tom::ObservationSimulationTimeSettings<double>>>(m, "ObservationSimulationTimeSettings");

    py::class_<tom::TabulatedObservationSimulationTimeSettings<double>,
               std::shared_ptr<tom::TabulatedObservationSimulationTimeSettings<double>>,
               tom::ObservationSimulationTimeSettings<double> >(m, "TabulatedObservationSimulationTimeSettings")
            .def(py::init<
                 const tom::LinkEndType, const std::vector< double > >(),
                 py::arg("reference_link_end"),
                 py::arg("observation_times") );

    py::class_<tom::ArcLimitedObservationSimulationTimeSettings<double>,
               std::shared_ptr<tom::ArcLimitedObservationSimulationTimeSettings<double>>,
               tom::ObservationSimulationTimeSettings<double> >(m, "ArcLimitedObservationSimulationTimeSettings")
            .def(py::init<
                 const tom::LinkEndType,
                 const double,
                 const double,
                 const double,
                 const double,
                 const int >(),
                 py::arg("reference_link_end"),
                 py::arg("start_time"),
                 py::arg("end_time"),
                 py::arg("observation_interval"),
                 py::arg("arc_duration"),
                 py::arg("observations_per_arc") );

    m.def("simulate_observations",
          py::overload_cast<
          const std::map< tom::ObservableType, std::map< tom::LinkEnds, std::pair< std::vector< double >, tom::LinkEndType > > >&,
          const std::map< tom::ObservableType, std::shared_ptr< tom::ObservationSimulatorBase< double, double > > >&,
          const tom::PerObservableObservationViabilityCalculatorList >(
              &tom::simulateObservations< double, double > ),
          py::arg("observation_to_simulate"),
          py::arg("observation_simulator"),
          py::arg("observation_viability_calculators") );

    m.def("simulate_observations",
          py::overload_cast<
          const std::map< tom::ObservableType, std::map< tom::LinkEnds, std::shared_ptr< tom::ObservationSimulationTimeSettings< double > > > >&,
          const std::map< tom::ObservableType, std::shared_ptr< tom::ObservationSimulatorBase< double, double > > >&,
          const tom::PerObservableObservationViabilityCalculatorList >(
              &tom::simulateObservations< double, double > ),
          py::arg("observation_to_simulate"),
          py::arg("observation_simulator"),
          py::arg("observation_viability_calculators") );
}

}// namespace tudatpy
