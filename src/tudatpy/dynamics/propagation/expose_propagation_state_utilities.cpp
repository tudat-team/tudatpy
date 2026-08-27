/*    Copyright (c) 2010-2019, Delft University of Technology
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
#include "expose_propagation_bindings.h"

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <tudat/astro/aerodynamics/aerodynamicGuidance.h>
#include <tudat/astro/basic_astro/accelerationModel.h>
#include <tudat/astro/basic_astro/massRateModel.h>
#include <tudat/astro/basic_astro/torqueModel.h>
#include <tudat/astro/propagators/getZeroProperModeRotationalInitialState.h>
#include <tudat/astro/propagators/singleStateTypeDerivative.h>
#include <tudat/simulation/environment_setup/body.h>
#include <tudat/simulation/propagation_setup/propagationSettings.h>

#include "scalarTypes.h"

namespace py = pybind11;

namespace ta = tudat::aerodynamics;
namespace tp = tudat::propagators;
namespace tba = tudat::basic_astrodynamics;
namespace tss = tudat::simulation_setup;

namespace tudat
{

namespace aerodynamics
{

class PyAerodynamicGuidance : public ta::AerodynamicGuidance
{
public:
    using AerodynamicGuidance::AerodynamicGuidance;

    using AerodynamicGuidance::currentAngleOfAttack_;
    using AerodynamicGuidance::currentAngleOfSideslip_;
    using AerodynamicGuidance::currentBankAngle_;

    void updateGuidance( const double currentTime ) override
    {
        PYBIND11_OVERLOAD_PURE( void, AerodynamicGuidance, updateGuidance, currentTime );
    }
};

}  // namespace aerodynamics

}  // namespace tudat

namespace tudatpy
{
namespace dynamics
{
namespace propagation
{

void expose_propagation_state_utility_types( py::module& m )
{
    py::class_< tba::TorqueModel, std::shared_ptr< tba::TorqueModel > >( m, "TorqueModel" );

    py::class_< tba::AccelerationModel< Eigen::Vector3d >, std::shared_ptr< tba::AccelerationModel< Eigen::Vector3d > > >(
            m, "AccelerationModel" );

    py::class_< tba::MassRateModel, std::shared_ptr< tba::MassRateModel > >( m, "MassRateModel" );
}

void expose_propagation_state_utility_bindings( py::module& m )
{
    py::class_< ta::AerodynamicGuidance, ta::PyAerodynamicGuidance, std::shared_ptr< ta::AerodynamicGuidance > >( m, "AerodynamicGuidance" )
            .def( py::init<>( ) )
            .def( "updateGuidance",
                  &ta::AerodynamicGuidance::updateGuidance,
                  py::arg( "current_time" ) )  // The current_time parameter is now expected to be a Time object
            .def_readwrite( "angle_of_attack", &ta::PyAerodynamicGuidance::currentAngleOfAttack_ )
            .def_readwrite( "bank_angle", &ta::PyAerodynamicGuidance::currentBankAngle_ )
            .def_readwrite( "sideslip_angle", &ta::PyAerodynamicGuidance::currentAngleOfSideslip_ );

    m.def( "get_single_integration_differential_equation_order",
           &tp::getSingleIntegrationDifferentialEquationOrder,
           py::arg( "state_type" ) );

    m.def( "get_generalized_acceleration_size", &tp::getGeneralizedAccelerationSize, py::arg( "state_type" ) );

    m.def( "get_state_of_bodies",
           py::overload_cast< const std::vector< std::string >&,
                              const std::vector< std::string >&,
                              const tss::SystemOfBodies&,
                              const TIME_TYPE >( &tp::getInitialStatesOfBodies< TIME_TYPE, STATE_SCALAR_TYPE > ),
           py::arg( "bodies_to_propagate" ),
           py::arg( "central_bodies" ),
           py::arg( "body_system" ),
           py::arg( "initial_time" ),
           R"doc(

 Function to get the translational states of a set of bodies, with respect to some set of central bodies, at the requested time.

 Function to get the translational states of a set of bodies, with respect to some set of central bodies, at the requested time. This function
 is typically used to extract an initial state for a propagation of a set of bodies, for which the initial state is extracted from the
 existing ephemerides of the bodies.


 Parameters
 ----------
 bodies_to_propagate : list[str]
     List of names of bodies for which the state is to be extracted
 central_bodies : list[str]
     List of central bodies, w.r.t. which the states are to be computed (in the same order as ``bodies_to_propagate``)
 body_system : SystemOfBodies
     System of bodies that define the environment
 initial_time : astro.time_representation.Time
     Time at which the states are to be extracted from the environment (Time object representing seconds since J2000 TDB)
 Returns
 -------
 numpy.ndarray
     Vector of size :math:`6\times N`, with the translational states of each entry of body from
     ``bodies_to_propagate`` w.r.t. the corresponding central body in ``central_bodies``.







     )doc" );

    m.def( "get_initial_state_of_bodies",
           py::overload_cast< const std::vector< std::string >&,
                              const std::vector< std::string >&,
                              const tss::SystemOfBodies&,
                              const TIME_TYPE >( &tp::getInitialStatesOfBodies< TIME_TYPE, STATE_SCALAR_TYPE > ),
           py::arg( "bodies_to_propagate" ),
           py::arg( "central_bodies" ),
           py::arg( "body_system" ),
           py::arg( "initial_time" ) );

    m.def( "get_initial_state_of_body",  // overload [2/2]
           py::overload_cast< const std::string&, const std::string&, const tss::SystemOfBodies&, const TIME_TYPE >(
                   &tp::getInitialStateOfBody< TIME_TYPE, STATE_SCALAR_TYPE > ),
           py::arg( "body_to_propagate" ),
           py::arg( "central_body" ),
           py::arg( "bodies" ),
           py::arg( "initial_time" ) );

    m.def( "get_initial_rotational_state_of_body",
           py::overload_cast< const std::string&, const std::string&, const tss::SystemOfBodies&, const TIME_TYPE >(
                   &tp::getInitialRotationalStateOfBody< TIME_TYPE, STATE_SCALAR_TYPE > ),
           py::arg( "body_to_propagate" ),
           py::arg( "base_orientation" ),
           py::arg( "bodies" ),
           py::arg( "initial_time" ) );

    py::class_< tp::DampedInitialRotationalStateResults< STATE_SCALAR_TYPE >,
                std::shared_ptr< tp::DampedInitialRotationalStateResults< STATE_SCALAR_TYPE > > >( m,
                                                                                                   "RotationalProperModeDampingResults",
                                                                                                   R"doc(

         Object that stores the results of the algorithm to damp the proper mode of rotational dynamics for an initial state,
         as computed by the :func:`~get_damped_proper_mode_initial_rotational_state` function






      )doc" )
            .def_readwrite( "damped_initial_state",
                            &tp::DampedInitialRotationalStateResults< STATE_SCALAR_TYPE >::initialState_,
                            R"doc(

         Initial state produced by the damping algorithm, for which the signature of the proper mode should be
         removed (or at least, substantially reduced). Note that this initial state corresponds to the *full* state vector
         that is provided to the ``get_damped_proper_mode_initial_rotational_state`` function (e.g. is size 7
         for rotational dynamics of a single body, size 13 for coupled orbital-rotational dynamics of a single body, etc.)


         :type: numpy.ndarray
      )doc" )
            .def_readwrite( "forward_backward_states",
                            &tp::DampedInitialRotationalStateResults< STATE_SCALAR_TYPE >::forwardBackwardPropagatedStates_,
                            R"doc(

         Data structure that contains the full state histories used by the damping algorithm. The contents are are as follows:

         * The :math:`i^{th}` entry of the list corresponds to the :math:`i^{th}` iteration of the forward-backward propagation
         * Each tuple in the list contains two dictionaries, the first one corresponding to the forward propagation results, the seconds one to the backward propagation results


         :type: list[tuple[dict[float,numpy.ndarray],dict[float,numpy.ndarray]]]
      )doc" )
            .def_readwrite( "forward_backward_dependent_variables",
                            &tp::DampedInitialRotationalStateResults< STATE_SCALAR_TYPE >::forwardBackwardDependentVariables_,
                            R"doc(

         As ``forward_backward_states``, but for the dependent variables.


         :type: list[tuple[dict[float,numpy.ndarray],dict[float,numpy.ndarray]]]
      )doc" );

    m.def( "get_damped_proper_mode_initial_rotational_state",
           py::overload_cast< const tss::SystemOfBodies&,
                              const std::shared_ptr< tp::SingleArcPropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE > >,
                              const double,
                              const std::vector< double >,
                              const bool >( &tp::getZeroProperModeRotationalStateWithStruct< TIME_TYPE, STATE_SCALAR_TYPE > ),
           py::arg( "bodies" ),
           py::arg( "propagator_settings" ),
           py::arg( "body_mean_rotational_rate" ),
           py::arg( "dissipation_times" ),
           py::arg( "propagate_undamped" ) = true,
           R"doc(

 Function to compute an initial rotational state for which the proper mode of rotation is damped.

 Function to compute an initial rotational state for which the proper mode of rotation is damped, using the algorithm
 used by Rambaux et al. (2010) to compute an initial rotational state for Phobos. This algorithm propagates the
 dynamics of the system a number of times, with the settings specified by the user and a specific modification to
 damp the proper mode. Since a number of propagations are performed by this function, it may take some time to run.
 Specifically, the algorithm works as follows:

 * Introduce a damping torque (see below) to damp the proper mode, with damping time :math:`\tau_{d}`
 * Propagate the dynamics forward in time for a period of :math:`10\tau_{d}`
 * Remove the virtual torque, and propagate the dynamics back to the initial time :math:`t_{0}`
 * Repeat the above for the list of damping times provided by the user

 The state after the final backwards propagation to :math:`t_{0}` is provided as output by this function, to be
 used as damped initial state. The output from this function also provides the user access to the full state history
 and dependent variable history of the forward and backward propagations, to allow a user to track and validate
 the progress of the algorithm.

 The damping torque :math:`\Gamma` is defined as follows:

 .. math::
    \boldsymbol{\Gamma}= -\frac{1}{\tau_{d}}\mathbf{I}\begin{pmatrix}\omega_{x}\\ \omega_{y}\\ \omega_{x}-\omega_{p} \end{pmatrix}

 where :math:\mathbf{I}` is the body's inertia tensor (in its body-fixed frame), :math:`\tau_{d}` the damping time of the
 current propagation, and :math:`\omega_{x}, \omega_{y}, \omega_{z}` the body's current rotation about its
 body-fixed, x-, y- and z-axes, respectively. The damping torque is implemented to damp out all rotations along
 the body-fixed x- and y-axes, and any deviations from constant rotation with frequency :\omega_{p}: about the body-fixed z-axis.

 .. note:: The mean rotation rate of the body :math:`\omega_{p}` is a user-defined input, and must be tuned to the dynamics of the system.


 Parameters
 ----------
 bodies : SystemOfBodies
     Set of body objects that defines the environment
 propagator_settings : SingleArcPropagatorSettings
     Propagator settings for the dynamics of which the initial rotational state is to be damped. These propagator
     settings must be for rotational dynamics only, or for multi-type rotational dynamics that contains rotational
     dynamics for a single body (e.g. translational-rotational dynamics for a single body)

 body_mean_rotational_rate : float
     Mean rotational rate :math:`\omega_{p}` to which the damping algorithm will force the body-fixed rotation about its z-axis.
 dissipation_times : list[ float ]
     List of damping times :math:`\tau_{d}` for which the algorithm is to be run. Note that this list should be organized in ascending order for the algorithm to perform properly
 propagate_undamped : bool, default = True
     Boolean defining whether the first forward/backward propagation performed by the damping algorithm has damping turned off (damping turned off if True, damping turned on if False).
     Propagating without any damping before starting the damping algorithm is useful for verification purposes, but not required for the algorithm itself.

 Returns
 -------
 DampedInitialRotationalStateResults
     Object that contains the results of the damping algorithm (final damped rotational state, and forward/backward propagation results).






     )doc" );

    m.def( "combine_initial_states",
           &tp::createCombinedInitialState< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "propagator_settings_per_type" ),
           R"doc(

 Function to retrieve the initial state for a list of propagator settings.

 Function to retrieve the initial state for a list of propagator settings. This way, the initial state for
 different quantities to be propagated (e.g., translational state, rotational state, mass) are retrieved and
 organized in a single container.


 Parameters
 ----------
 propagator_settings_per_type : dict
     Propagator settings where the type of propagation is reported as key and the respective list of propagator settings as value.
 Returns
 -------
 numpy.ndarray
     Vector of initial states, sorted in order of IntegratedStateType, and then in the order of the vector of SingleArcPropagatorSettings of given type.






     )doc" );
}

}  // namespace propagation
}  // namespace dynamics
}  // namespace tudatpy
