/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_propagation_setup.h"

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "expose_propagation_setup/expose_acceleration_setup.h"
#include "expose_propagation_setup/expose_dependent_variable_setup.h"
#include "expose_propagation_setup/expose_integrator_setup.h"
#include "expose_propagation_setup/expose_mass_rate_setup.h"
#include "expose_propagation_setup/expose_propagator_setup.h"
#include "expose_propagation_setup/expose_thrust_setup.h"
#include "expose_propagation_setup/expose_torque_setup.h"

namespace py = pybind11;
namespace tba = tudat::basic_astrodynamics;
namespace tss = tudat::simulation_setup;
namespace tp = tudat::propagators;
namespace tinterp = tudat::interpolators;
namespace te = tudat::ephemerides;
namespace tni = tudat::numerical_integrators;
namespace trf = tudat::reference_frames;
namespace tmrf = tudat::root_finders;

namespace tudatpy
{
namespace numerical_simulation
{
namespace propagation_setup
{

void expose_propagation_setup( py::module &m )
{
    auto thrust_setup = m.def_submodule( "thrust" );
    thrust::expose_thrust_setup( thrust_setup );

    auto acceleration_setup = m.def_submodule( "acceleration" );
    acceleration::expose_acceleration_setup( acceleration_setup );

    auto torque_setup = m.def_submodule( "torque" );
    torque::expose_torque_setup( torque_setup );

    auto integrator_setup = m.def_submodule( "integrator" );
    integrator::expose_integrator_setup( integrator_setup );

    auto propagator_setup = m.def_submodule( "propagator" );
    propagator::expose_propagator_setup( propagator_setup );

    auto mass_setup = m.def_submodule( "mass_rate" );
    mass_rate::expose_mass_rate_setup( mass_setup );

    auto dependent_variable_setup = m.def_submodule( "dependent_variable" );
    dependent_variable::expose_dependent_variable_setup( dependent_variable_setup );

    m.def( "create_acceleration_models",
           py::overload_cast< const tss::SystemOfBodies &,
                              const tss::SelectedAccelerationMap &,
                              const std::vector< std::string > &,
                              const std::vector< std::string > & >( &tss::createAccelerationModelsMap ),
           py::arg( "body_system" ),
           py::arg( "selected_acceleration_per_body" ),
           py::arg( "bodies_to_propagate" ),
           py::arg( "central_bodies" ),
           R"doc(

Function to create a set of acceleration models from a dictionary of bodies linked to acceleration model types.

Function to create a set of acceleration models from a map of bodies and acceleration model types. The propagated
bodies and central bodies are provided as two separate lists with the same order.


Parameters
----------
body_system : SystemOfBodies
    System of bodies to be used in the propagation.
selected_acceleration_per_body : Dict[str, Dict[str, List[AccelerationSettings]]]
    Key-value container, with key denoting the body undergoing the acceleration, and the value containing an additional key-value container, with the body exerting acceleration, and list of acceleration settings exerted by this body.
bodies_to_propagate : list
    List of bodies to propagate.
central_bodies : list
    List of central bodies, each referred to each propagated body in the same order.
Returns
-------
AccelerationMap
    Set of accelerations acting on the bodies to propagate, provided as dual key-value container, similar to the acceleration settings input, but now with ``AccelerationModel`` lists as inner value





Examples
--------
In this example, the acceleration model includes a spherical harmonic (degree and order 5) gravitational acceleration
and aerodynamic acceleration exerted by the Earth, as well as point-mass gravity exerted by the Sun and the Moon.
The variable ``accelerations_settings_vehicle`` denotes the list of bodies exerting accelerations and the types of
accelerations, while the variable ``acceleration_settings`` associates this list with the body undergoing the
acceleration (``"Vehicle"``).

.. code-block:: python

   # Define bodies that are propagated.
   bodies_to_propagate = ["Vehicle"]

   # Define central bodies.
   central_bodies = ["Earth"]

   # Define accelerations acting on Vehicle
   accelerations_settings_vehicle = dict(
       Sun=[propagation_setup.acceleration.point_mass_gravity()],
       Moon=[propagation_setup.acceleration.point_mass_gravity()],
       Earth=[propagation_setup.acceleration.spherical_harmonic_gravity(5, 5),
              propagation_setup.acceleration.aerodynamic()]
   )

   # Create global accelerations settings dictionary
   acceleration_settings = {"Vehicle": accelerations_settings_vehicle}

   # Create acceleration models
   acceleration_models = propagation_setup.create_acceleration_models(
       bodies, acceleration_settings,  bodies_to_propagate, central_bodies)


    )doc" );

    m.def( "create_torque_models",
           &tss::createTorqueModelsMap,
           py::arg( "body_system" ),
           py::arg( "selected_torque_per_body" ),
           py::arg( "bodies_to_propagate" ),
           R"doc(

Function to create a set of acceleration models from a dictionary of bodies linked to acceleration model types.

Function to create a set of acceleration models from a map of bodies and acceleration model types. The propagated
bodies is provided as a list.


Parameters
----------
body_system : SystemOfBodies
    System of bodies to be used in the propagation.
selected_torque_per_body : Dict[str, Dict[str, List[TorqueSettings]]]
    Key-value container, with key denoting the body undergoing the torque, and the value containing an additional key-value container, with the body exerting torque, and list of torque settings exerted by this body.
bodies_to_propagate : list
    List of bodies to propagate.
Returns
-------
TorqueModelMap
    Set of torques acting on the bodies to propagate, provided as dual key-value container, similar to the torque settings input, but now with ``TorqueModel`` lists as inner value





Examples
--------

In this example, the following torques are exerted on the vehicle: spherical-harmonic gravitational torque
(up to order 4 and degree 4) and aerodynamic torque exerted by the Earth, second-degree gravitational torque
exerted by the Sun and the Moon.
The variable ``torques_settings_vehicle`` denotes the list of bodies exerting torques and the types of
torques, while the variable ``torque_settings`` associates this list with the body undergoing the
torque.

.. code-block:: python

  # Define bodies that are propagated
  bodies_to_propagate = ["Vehicle"]

  # Define torques per each exerting body
  torque_settings_vehicle = dict(
      Sun=[propagation_setup.torque.second_degree_gravitational()],
      Moon=[propagation_setup.torque.second_degree_gravitational()],
      Earth=[propagation_setup.torque.spherical_harmonic_gravitational(4, 4),
             propagation_setup.torque.aerodynamic()]
  )

  # Create global torque settings dictionary
  torque_settings = {"Vehicle": torque_settings_vehicle}

  # Create torque models
  torque_models = propagation_setup.create_torque_models(
      bodies, torque_settings,  bodies_to_propagate )


    )doc" );

    m.def( "create_mass_rate_models",
           &tss::createMassRateModelsMap,
           py::arg( "body_system" ),
           py::arg( "selected_mass_rates_per_body" ),
           py::arg( "acceleration_models" ) = nullptr,
           R"doc(

Function to create a set of mass-rate models from associated settings.

Function to create a set of mass-rate models from a map of bodies and mass-rate model types.
If the mass-rate depends on any acceleration models (e.g. thrust), the acceleration
models must be provided as an input.


Parameters
----------
body_system : SystemOfBodies
    System of bodies to be used in the propagation.
selected_mass_rates_per_body : Dict[str, List[MassRateModelSettings]]
    Key-value container, with key denoting the body with changing mass, and the value containing a list of mass rate settings (in most cases, this list will have only a single entry)
acceleration_models : AccelerationMap
    Sorted list of acceleration models, as created by :func:`create_acceleration_models`
Returns
-------
MassRateModelMap
    Set of mass-rate models, as key-value container, same as the settings input, with the difference that the rate settings objects have been processed into the associated objects calculating the actual mass-rate changes.






    )doc" );
}
}  // namespace propagation_setup
}  // namespace numerical_simulation
}  // namespace tudatpy
