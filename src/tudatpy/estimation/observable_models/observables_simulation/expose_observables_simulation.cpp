/*    Copyright (c) 2010-2021, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif
#include "expose_observables_simulation.h"
#include <pybind11/functional.h>
#include "scalarTypes.h"
#include "tudat/astro/observation_models/observationManager.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"

namespace tom = tudat::observation_models;

namespace tudatpy
{
namespace estimation
{
namespace observable_models
{

namespace observables_simulation
{

void expose_observables_simulation( py::module& m )
{
    py::class_< tom::ObservationManagerBase< STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tom::ObservationManagerBase< STATE_SCALAR_TYPE, TIME_TYPE > > >( m, "ObservationManager" );

    py::class_< tom::ObservationViabilityCalculator, std::shared_ptr< tom::ObservationViabilityCalculator > >(
            m,
            "ObservationViabilityCalculator",
            R"doc(

         Template class for observation viability calculators.

         Template class for classes which conducts viability calculations on simulated observations.
         Instances of the applicable ObservationViabilityCalculators are automatically created from the given :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` objects during the simulation of observations (:func:`~tudatpy.estimation.observations_setup.observations_wrapper.simulate_observations`).
         The user typically does not interact directly with this class.





      )doc" )
            .def( "is_observation_viable",
                  &tom::ObservationViabilityCalculator::isObservationViable,
                  py::arg( "link_end_states" ),
                  py::arg( "link_end_times" ),
                  py::arg( "observation_value" ) = Eigen::VectorXd( ),
                  R"doc(

         Function to check whether an observation is viable.

         Function to check whether an observation is viable.
         The calculation is performed based on the given times and link end states, and the current observation value.
         Note, that this function is called automatically during the simulation of observations.
         Direct calls to this function are generally not required.

         Parameters
         ----------
         link_end_states : List[ numpy.ndarray[numpy.float64[6, 1]] ]
             Vector of states of the link ends involved in the observation.
         link_end_times : list[float]
             Vector of times at the link ends involved in the observation.
         observation_value : numpy.ndarray[numpy.float64[]]
             Current simulated observation value.
         Returns
         -------
         bool
             True if observation is viable, false if not.   

     )doc" );

    py::class_< tom::ObservationSimulatorBase< STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tom::ObservationSimulatorBase< STATE_SCALAR_TYPE, TIME_TYPE > > >( m,
                                                                                                    "ObservationSimulator",
                                                                                                    R"doc(

         Class hosting the functionality for simulating observations.

         Class hosting the functionality for simulating a given observable over a defined link geometry.
         Instances of this class are automatically created from the given :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationModelSettings` objects upon instantiation of the :class:`~tudatpy.estimation.estimation_analysis.Estimator` class.





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
}

}  // namespace observables_simulation
}  // namespace observable_models
}  // namespace estimation
}  // namespace tudatpy
