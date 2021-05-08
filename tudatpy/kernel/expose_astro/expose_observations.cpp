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

#include <stdio.h>
#include <time.h>

#include <tudat/astro/observation_models.h>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>

namespace py = pybind11;
namespace tom = tudat::observation_models;

namespace tudat
{

namespace observation_models
{

LinkEndType getDefaultReferenceLinkEndType(
        const ObservableType observableType )
{
    LinkEndType referenceLinkEndType;
    switch( observableType )
    {
    case one_way_range:
        referenceLinkEndType = receiver;
        break;
    case angular_position:
        referenceLinkEndType = receiver;
        break;
    case one_way_doppler:
        referenceLinkEndType = receiver;
        break;
    case one_way_differenced_range:
        referenceLinkEndType = receiver;
        break;
    case n_way_range:
        referenceLinkEndType = receiver;
        break;
    case two_way_doppler:
        referenceLinkEndType = receiver;
        break;
    case position_observable:
        referenceLinkEndType = observed_body;
        break;
    case velocity_observable:
        referenceLinkEndType = observed_body;
        break;
    case euler_angle_313_observable:
        referenceLinkEndType = observed_body;
        break;
    default:
        throw std::runtime_error( "Error, default reference link end not defined for observable " +
                                  std::to_string( observableType ) );
    }
    return referenceLinkEndType;
}

template< typename ObservationScalarType = double, typename TimeType = double >
std::map< ObservableType, std::map< LinkEnds, std::pair< std::vector< TimeType >, LinkEndType > > > getObservationSettingsMapFormat(
        const std::vector< std::tuple< tom::ObservableType, LinkEnds, std::vector< TimeType > > >& observationsToSimulate )
{
    std::map< ObservableType, std::map< LinkEnds, std::pair< std::vector< TimeType >, LinkEndType > > > sortedObservationsToSimulate;
    for( int i = 0; i < observationsToSimulate.size( ); i++ )
    {
        ObservableType currentObservableType = std::get< 0 >( observationsToSimulate.at( i ) );
        LinkEnds currentLinkEnds = std::get< 1 >( observationsToSimulate.at( i ) );
        std::vector< TimeType >  currentLinkEndTimes = std::get< 2 >( observationsToSimulate.at( i ) );
        LinkEndType currentReferenceLinkEndType = getDefaultReferenceLinkEndType( currentObservableType );

        sortedObservationsToSimulate[ currentObservableType ][ currentLinkEnds ] = std::make_pair(
                    currentLinkEndTimes, currentReferenceLinkEndType );
    }
    return sortedObservationsToSimulate;
}

template< typename ObservationScalarType = double, typename TimeType = double >
std::vector< std::tuple< ObservableType, LinkEnds, Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >,
std::vector< TimeType >, LinkEndType > > getObservationsVectorFormat(
        const std::map< ObservableType, std::map< LinkEnds, std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >,
                std::pair< std::vector< TimeType >, LinkEndType > > > >& sortedObservations )
{
    std::vector< std::tuple< ObservableType, LinkEnds, Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >,
    std::vector< TimeType >, LinkEndType > > observationsOutput;

    for( auto observableIt : sortedObservations )
    {
        for( auto linkEndsIt : observableIt.second )
        {
            ObservableType currentObservableType = observableIt.first;
            LinkEnds currentLinkEnds = linkEndsIt.first;
            Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > currentObservable = linkEndsIt.second.first;
            std::vector< TimeType > currentTimes = linkEndsIt.second.second.first;
            LinkEndType currentReferenceLinkEndsType = linkEndsIt.second.second.second;

            observationsOutput.push_back(
                        std::make_tuple(
                            currentObservableType, currentLinkEnds, currentObservable,
                            currentTimes, currentReferenceLinkEndsType ) );
        }
    }
    return observationsOutput;
}


template< typename ObservationScalarType = double, typename TimeType = double >
std::vector< std::tuple< ObservableType, LinkEnds, Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >,
std::vector< TimeType >, LinkEndType > >
simulateObservations(
        const std::vector< std::tuple< tom::ObservableType, LinkEnds, std::vector< TimeType > > >& observationsToSimulate,
        const std::map< ObservableType, std::shared_ptr< ObservationSimulatorBase< ObservationScalarType, TimeType > > >&
        observationSimulators,
        const PerObservableObservationViabilityCalculatorList viabilityCalculatorList =
        PerObservableObservationViabilityCalculatorList( ) )
{
    std::map< ObservableType, std::map< LinkEnds, std::pair< std::vector< TimeType >, LinkEndType > > > sortedObservationsToSimulate =
        getObservationSettingsMapFormat( observationsToSimulate );

    std::map< ObservableType, std::map< LinkEnds, std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >,
            std::pair< std::vector< TimeType >, LinkEndType > > > > sortedObservations = simulateObservations(
                sortedObservationsToSimulate, observationSimulators, viabilityCalculatorList );

    return getObservationsVectorFormat( sortedObservations );
}

template< typename ObservationScalarType = double, typename TimeType = double >
std::vector< std::tuple< ObservableType, LinkEnds, Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >,
std::vector< TimeType >, LinkEndType > >
simulateObservationsWithNoise(
        const std::vector< std::tuple< tom::ObservableType, LinkEnds, std::vector< TimeType > > >& observationsToSimulate,
        const std::map< ObservableType, std::shared_ptr< ObservationSimulatorBase< ObservationScalarType, TimeType > > >&
        observationSimulators,
        const std::map< ObservableType, std::function< double( const double ) > >& noiseFunctions,
        const PerObservableObservationViabilityCalculatorList viabilityCalculatorList =
        PerObservableObservationViabilityCalculatorList( ) )
{
    std::map< ObservableType, std::map< LinkEnds, std::pair< std::vector< TimeType >, LinkEndType > > > sortedObservationsToSimulate =
        getObservationSettingsMapFormat( observationsToSimulate );

    std::map< ObservableType, std::map< LinkEnds, std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >,
            std::pair< std::vector< TimeType >, LinkEndType > > > > sortedObservations = simulateObservationsWithNoise(
                createObservationSimulationTimeSettingsMap(
                        sortedObservationsToSimulate ), observationSimulators, noiseFunctions, viabilityCalculatorList );

    return getObservationsVectorFormat( sortedObservations );
}

}

}

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
          py::arg("observation_viability_calculators") = tom::PerObservableObservationViabilityCalculatorList( ) );

    m.def("simulate_observations",
          py::overload_cast<
          const std::map< tom::ObservableType, std::map< tom::LinkEnds, std::shared_ptr< tom::ObservationSimulationTimeSettings< double > > > >&,
          const std::map< tom::ObservableType, std::shared_ptr< tom::ObservationSimulatorBase< double, double > > >&,
          const tom::PerObservableObservationViabilityCalculatorList >(
              &tom::simulateObservations< double, double > ),
          py::arg("observation_to_simulate"),
          py::arg("observation_simulator"),
          py::arg("observation_viability_calculators") = tom::PerObservableObservationViabilityCalculatorList( ) );

    m.def("simulate_observations",
          py::overload_cast<
          const std::vector< std::tuple< tom::ObservableType, tom::LinkEnds, std::vector< double > > >&,
          const std::map< tom::ObservableType, std::shared_ptr< tom::ObservationSimulatorBase< double, double > > >&,
          const tom::PerObservableObservationViabilityCalculatorList >(
              &tom::simulateObservations< double, double > ),
          py::arg("observation_to_simulate"),
          py::arg("observation_simulator"),
          py::arg("observation_viability_calculators") = tom::PerObservableObservationViabilityCalculatorList( ) );

    m.def("simulate_noisy_observations",
          py::overload_cast<
          const std::vector< std::tuple< tom::ObservableType, tom::LinkEnds, std::vector< double > > >&,
          const std::map< tom::ObservableType, std::shared_ptr< tom::ObservationSimulatorBase< double, double > > >&,
          const std::map< tom::ObservableType, std::function< double( const double ) > >&,
          const tom::PerObservableObservationViabilityCalculatorList >(
              &tom::simulateObservationsWithNoise< double, double > ),
          py::arg("observation_to_simulate"),
          py::arg("observation_simulator"),
          py::arg("noise_functions"),
          py::arg("observation_viability_calculators") = tom::PerObservableObservationViabilityCalculatorList( ) );

    m.def("gaussian_noise_function",
              &tom::getGaussianDistributionNoiseFunction,
          py::arg("standard_deviation"),
          py::arg("mean") = 0.0,
          py::arg("seed") = time(NULL) );

}

}// namespace tudatpy
