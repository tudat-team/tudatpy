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

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "scalarTypes.h"
#include "tudat/astro/propagators/propagateCovariance.h"
#include "tudat/basics/utilities.h"
#include "tudat/simulation/estimation_setup/fitOrbitToEphemeris.h"

namespace py = pybind11;
namespace tep = tudat::estimatable_parameters;
namespace tp = tudat::propagators;
namespace ts = tudat::statistics;
namespace tss = tudat::simulation_setup;
namespace tni = tudat::numerical_integrators;
namespace tom = tudat::observation_models;
namespace trf = tudat::reference_frames;
namespace te = tudat::ephemerides;

namespace tudat
{

namespace simulation_setup
{

std::pair< std::vector< double >, std::vector< Eigen::VectorXd > > getTargetAnglesAndRangeVector(
        const simulation_setup::SystemOfBodies &bodies,
        const std::pair< std::string, std::string > groundStationId,
        const std::string &targetBody,
        const std::vector< double > times,
        const bool transmittingToTarget )
{
    std::map< double, Eigen::VectorXd > targetAnglesAndRange = getTargetAnglesAndRange(
            bodies, groundStationId, targetBody, times, transmittingToTarget );
    return std::make_pair( utilities::createVectorFromMapKeys( targetAnglesAndRange ),
                           utilities::createVectorFromMapValues( targetAnglesAndRange ) );
}

}  // namespace simulation_setup

}  // namespace tudat

namespace tudatpy
{
namespace numerical_simulation
{
namespace estimation
{

PYBIND11_MODULE( expose_estimation, m )
// ( py::module &m )
{
    py::class_< tom::ObservationViabilityCalculator,
                std::shared_ptr< tom::ObservationViabilityCalculator > >(
            m,
            "ObservationViabilityCalculator",
            R"doc(

        Template class for observation viability calculators.

        Template class for classes which conducts viability calculations on simulated observations.
        Instances of the applicable ObservationViabilityCalculators are automatically created from the given :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects during the simulation of observations (:func:`~tudatpy.numerical_simulation.estimation.simulate_observations`).
        The user typically does not interact directly with this class.





     )doc" )
            .def( "is_observation_viable",
                  &tom::ObservationViabilityCalculator::isObservationViable,
                  py::arg( "link_end_states" ),
                  py::arg( "link_end_times" ),
                  R"doc(

        Function to check whether an observation is viable.

        Function to check whether an observation is viable.
        The calculation is performed based on the given times and link end states.
        Note, that this function is called automatically during the simulation of observations.
        Direct calls to this function are generally not required.


        Parameters
        ----------
        link_end_states : List[ numpy.ndarray[numpy.float64[6, 1]] ]
            Vector of states of the link ends involved in the observation.
        link_end_times : List[float]
            Vector of times at the link ends involved in the observation.
        Returns
        -------
        bool
            True if observation is viable, false if not.





    )doc" );

    py::class_< tom::ObservationSimulatorBase< STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tom::ObservationSimulatorBase< STATE_SCALAR_TYPE, TIME_TYPE > > >(
            m,
            "ObservationSimulator",
            R"doc(

        Class hosting the functionality for simulating observations.

        Class hosting the functionality for simulating a given observable over a defined link geometry.
        Instances of this class are automatically created from the given :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSettings` objects upon instantiation of the :class:`~tudatpy.numerical_simulation.Estimator` class.





     )doc" );

    py::class_< tom::ObservationSimulator< 1, STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tom::ObservationSimulator< 1, STATE_SCALAR_TYPE, TIME_TYPE > >,
                tom::ObservationSimulatorBase< STATE_SCALAR_TYPE, TIME_TYPE > >(
            m, "ObservationSimulator_1", R"doc(No documentation found.)doc" );

    py::class_< tom::ObservationSimulator< 2, STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tom::ObservationSimulator< 2, STATE_SCALAR_TYPE, TIME_TYPE > >,
                tom::ObservationSimulatorBase< STATE_SCALAR_TYPE, TIME_TYPE > >(
            m, "ObservationSimulator_2", R"doc(No documentation found.)doc" );

    py::class_< tom::ObservationSimulator< 3, STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tom::ObservationSimulator< 3, STATE_SCALAR_TYPE, TIME_TYPE > >,
                tom::ObservationSimulatorBase< STATE_SCALAR_TYPE, TIME_TYPE > >(
            m, "ObservationSimulator_3", R"doc(No documentation found.)doc" );

    py::class_< tom::ObservationSimulator< 6, STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tom::ObservationSimulator< 6, STATE_SCALAR_TYPE, TIME_TYPE > >,
                tom::ObservationSimulatorBase< STATE_SCALAR_TYPE, TIME_TYPE > >(
            m, "ObservationSimulator_6", R"doc(No documentation found.)doc" );

    m.def( "simulate_observations",
           &tss::simulateObservations< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "simulation_settings" ),
           py::arg( "observation_simulators" ),
           py::arg( "bodies" ),
           R"doc(

Function to simulate observations.

Function to simulate observations from set observation simulators and observation simulator settings.
Automatically iterates over all provided observation simulators, generating the full set of simulated observations.


Parameters
----------
observation_to_simulate : List[ :class:`ObservationSimulationSettings` ]
    List of settings objects, each object providing the observation time settings for simulating one type of observable and link end set.

observation_simulators : List[ :class:`~tudatpy.numerical_simulation.estimation.ObservationSimulator` ]
    List of :class:`~tudatpy.numerical_simulation.estimation.ObservationSimulator` objects, each object hosting the functionality for simulating one type of observable and link end set.

bodies : :class:`~tudatpy.numerical_simulation.environment.SystemOfBodies`
    Object consolidating all bodies and environment models, including ground station models, that constitute the physical environment.

Returns
-------
:class:`~tudatpy.numerical_simulation.estimation.ObservationCollection`
    Object collecting all products of the observation simulation.






    )doc" );

    m.def( "compute_residuals_and_dependent_variables",
           &tss::computeResidualsAndDependentVariables< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "observation_collection" ),
           py::arg( "observation_simulators" ),
           py::arg( "bodies" ),
           R"doc(No documentation found.)doc" );

    m.def( "create_pseudo_observations_and_models",
           &tss::simulatePseudoObservations< TIME_TYPE, STATE_SCALAR_TYPE >,
           py::arg( "bodies" ),
           py::arg( "observed_bodies" ),
           py::arg( "central_bodies" ),
           py::arg( "initial_time" ),
           py::arg( "final_time" ),
           py::arg( "time_step" ),
           R"doc(No documentation found.)doc" );

    m.def( "create_best_fit_to_ephemeris",
           &tss::createBestFitToCurrentEphemeris< TIME_TYPE, STATE_SCALAR_TYPE >,
           py::arg( "bodies" ),
           py::arg( "acceleration_models" ),
           py::arg( "observed_bodies" ),
           py::arg( "central_bodies" ),
           py::arg( "integrator_settings" ),
           py::arg( "initial_time" ),
           py::arg( "final_time" ),
           py::arg( "data_point_interval" ),
           py::arg( "additional_parameter_names" ) =
                   std::vector< std::shared_ptr< tep::EstimatableParameterSettings > >( ),
           py::arg( "number_of_iterations" ) = 3,
           py::arg( "reintegrate_variational_equations" ) = true,
           py::arg( "results_print_frequency" ) = 0.0,
           R"doc(No documentation found.)doc" );

    m.def( "set_existing_observations",
           &tss::setExistingObservations< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "observations" ),
           py::arg( "reference_link_end" ),
           py::arg( "ancilliary_settings_per_observatble" ) =
                   std::map< tom::ObservableType,
                             std::shared_ptr< tom::ObservationAncilliarySimulationSettings > >( ) );

    m.def( "compute_target_angles_and_range",
           &tss::getTargetAnglesAndRange,
           py::arg( "bodies" ),
           py::arg( "station_id" ),
           py::arg( "target_body" ),
           py::arg( "observation_times" ),
           py::arg( "is_station_transmitting" ),
           R"doc(

Function to compute the azimuth angle, elevation angle and range at a ground station.

Function to compute the azimuth angle, elevation angle and range at a ground station. This functions is provided as a function of
convenience, to prevent users having to manually define the relevant settings for this often-needed functionality. This function
takes an observing station and a target body as input, and provides the observed angles and current range (without correction for aberrations, with correction for light time)
as observed at that station


Parameters
----------
bodies : SystemOfBodies
    System of bodies that defines the full physical environment

station_id : tuple[ str, str]
    Identifier for the observing station, as a pair of strings: the body name and the station name.

target_body : str
    Name of body which is observed by ground station

observation_times : list[float]
    List of times at which the ground station observations are to be analyzed

is_station_transmitting : bool
    Boolean defining whether the observation times define times at which the station is transmitting to, or receiving from, the ground station.
    This has an impact on the whether the light-time is computed forward or backward in time from the ground station to the target

Returns
-------
dict[float,numpy.ndarray[numpy.float64[3, 1]]]
    Dictionary with the required output. Key defines the observation time, the value is an array of size three containing entry 0 - elevation angle, entry 1 - azimuth angle, entry 2 - range






    )doc" );

    m.def( "compute_target_angles_and_range_vectors",
           &tss::getTargetAnglesAndRangeVector,
           py::arg( "bodies" ),
           py::arg( "station_id" ),
           py::arg( "target_body" ),
           py::arg( "observation_times" ),
           py::arg( "is_station_transmitting" ),
           R"doc(

Function to compute the azimuth angle, elevation angle and range at a ground station.

Function to compute the azimuth angle, elevation angle and range at a ground station. This functions is provided as a function of
convenience, to prevent users having to manually define the relevant settings for this often-needed functionality. This function
takes an observing station and a target body as input, and provides the observed angles and current range (without correction for aberrations, with correction for light time)
as observed at that station


Parameters
----------
bodies : SystemOfBodies
    System of bodies that defines the full physical environment

station_id : tuple[ str, str]
    Identifier for the observing station, as a pair of strings: the body name and the station name.

target_body : str
    Name of body which is observed by ground station

observation_times : list[float]
    List of times at which the ground station observations are to be analyzed

is_station_transmitting : bool
    Boolean defining whether the observation times define times at which the station is transmitting to, or receiving from, the ground station.
    This has an impact on the whether the light-time is computed forward or backward in time from the ground station to the target

Returns
-------
dict[float,numpy.ndarray[numpy.float64[3, 1]]]
    Dictionary with the required output. Key defines the observation time, the value is an array of size three containing entry 0 - elevation angle, entry 1 - azimuth angle, entry 2 - range






    )doc" );

    m.def( "create_filtered_observation_collection",
           py::overload_cast< const std::shared_ptr<
                                      tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > >,
                              const std::map< std::shared_ptr< tom::ObservationCollectionParser >,
                                              std::shared_ptr< tom::ObservationFilterBase > > & >(
                   &tom::filterObservations< STATE_SCALAR_TYPE, TIME_TYPE > ),
           py::arg( "original_observation_collection" ),
           py::arg( "observation_filters_map" ),
           R"doc(No documentation found.)doc" );

    m.def( "create_filtered_observation_collection",
           py::overload_cast< const std::shared_ptr<
                                      tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > >,
                              const std::shared_ptr< tom::ObservationFilterBase >,
                              const std::shared_ptr< tom::ObservationCollectionParser > >(
                   &tom::filterObservations< STATE_SCALAR_TYPE, TIME_TYPE > ),
           py::arg( "original_observation_collection" ),
           py::arg( "observation_filter" ),
           py::arg( "observation_parser" ) =
                   std::make_shared< tom::ObservationCollectionParser >( ),
           R"doc(No documentation found.)doc" );

    m.def( "split_observation_collection",
           &tom::splitObservationSets< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "original_observation_collection" ),
           py::arg( "observation_set_splitter" ),
           py::arg( "observation_parser" ) =
                   std::make_shared< tom::ObservationCollectionParser >( ),
           R"doc(No documentation found.)doc" );

    m.def( "create_new_observation_collection",
           &tom::createNewObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "original_observation_collection" ),
           py::arg( "observation_parser" ) =
                   std::make_shared< tom::ObservationCollectionParser >( ),
           R"doc(No documentation found.)doc" );

    // Include source code from other files
    expose_estimation_observation_collection( m );
    expose_estimation_propagated_covariance( m );
    expose_estimation_single_observation_set( m );
    expose_estimation_filter_parser( m );
}

}  // namespace estimation
}  // namespace numerical_simulation
}  // namespace tudatpy
