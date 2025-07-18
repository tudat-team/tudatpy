/*    Copyright (c) 2010-2021, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#define PYBIND11_DETAILED_ERROR_MESSAGES
#include "expose_integrator.h"

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <tudat/simulation/propagation_setup.h>

#include "scalarTypes.h"

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
namespace integrator
{

void expose_integrator_setup( py::module &m )
{
    // ENUMS
    py::enum_< tni::MinimumIntegrationTimeStepHandling >(
            m, "MinimumIntegrationTimeStepHandling", R"doc(

Enumeration defining possible behaviours when :math:`\Delta t_{rec}<\Delta t_{\min}`. in step-size control (e.g. recommended time step is smaller than minimum time step)





      )doc" )
            .value( "throw_exception_below_minimum",
                    tni::MinimumIntegrationTimeStepHandling::throw_exception_below_minimum,
                    R"doc(

The propagation is terminated and a :class:`tudatpy.exceptions.MinimumStepSizeViolatedError` is thrown.

)doc" )
            .value( "set_to_minimum_step_silently",
                    tni::MinimumIntegrationTimeStepHandling::set_to_minimum_step_silently,
                    R"doc(

The final time step is set to :math:`\Delta t=\Delta t_{\min}`, violating requirements of step-size control algorithm, without any message to user"

)doc" )
      .value( "set_to_minimum_step_single_warning",
                    tni::MinimumIntegrationTimeStepHandling::set_to_minimum_step_single_warning,
                    R"doc(

The final time step is set to :math:`\Delta t=\Delta t_{\min}`, violating requirements of step-size control algorithm, a warning is printed to the terminal the first time this happens during a propagation"

)doc" )
            .value( "set_to_minimum_step_every_time_warning",
                    tni::MinimumIntegrationTimeStepHandling::set_to_minimum_step_every_time_warning,
                    R"doc(

The final time step is set to :math:`\Delta t=\Delta t_{\min}`, violating requirements of step-size control algorithm, a warning is printed to the terminal every time this happens during a propagation"

)doc" )
            .export_values( );

    py::enum_< tni::AvailableIntegrators >( m, "AvailableIntegrators", R"doc(

         Enumeration of integrators available with tudat.





      )doc" )
            //       .value("euler_type",
            //       tni::AvailableIntegrators::euler)
            //       .value("runge_kutta_4_type",
            //       tni::AvailableIntegrators::rungeKutta4)
            .value( "runge_kutta_fixed_step_size_type",
                    tni::AvailableIntegrators::rungeKuttaFixedStepSize )
            .value( "runge_kutta_variable_step_size_type",
                    tni::AvailableIntegrators::rungeKuttaVariableStepSize )
            .value( "bulirsch_stoer_type", tni::AvailableIntegrators::bulirschStoer )
            .value( "adams_bashforth_moulton_type",
                    tni::AvailableIntegrators::adamsBashforthMoulton )
            .export_values( );

    py::enum_< tni::CoefficientSets >( m,
                                       "CoefficientSets",
                                       R"doc(

Coefficient sets for Runge-Kutta-type integrators.

Coefficient sets for Runge-Kutta-type integrators. The coefficients are defined
in a Butcher Tableau, with an coefficient set yielding an x(y) method yielding an integrator
with global truncation error of :math:`O(\Delta t^{x})`. Some of these coefficients also contain an embedded integrator of :math:`O(\Delta t^{y})`
for step size control.

      )doc" )
            .value( "euler_forward",
                    tni::forwardEuler,
                    R"doc(

Coefficients for the classic forward Euler method

)doc" )
            .value( "rk_4",
                    tni::rungeKutta4Classic,
                    R"doc(

Coefficients for the original Runge-Kutta method of order 4

)doc" )
            .value( "explicit_mid_point",
                    tni::explicitMidPoint,
                    R"doc(

Coefficients for the explicit midpoint method

)doc" )
            .value( "explicit_trapezoid_rule",
                    tni::explicitTrapezoidRule,
                    R"doc(

Coefficients for the explicit trapezoid rule, also called Heun's method or improved Euler's method

)doc" )
            .value( "ralston",
                    tni::ralston,
                    R"doc(

Coefficients for Ralston's method

)doc" )
            .value( "rk_3",
                    tni::rungeKutta3,
                    R"doc(

Coefficients for the Runge-Kutta method of order 3

)doc" )
            .value( "ralston_3",
                    tni::ralston3,
                    R"doc(

Coefficients for Ralston's third-order method

)doc" )
            .value( "SSPRK3",
                    tni::SSPRK3,
                    R"doc(

Coefficients for the Strong Stability Preserving Runge-Kutta third-order method

)doc" )
            .value( "ralston_4",
                    tni::ralston4,
                    R"doc(

Coefficients for Ralston's fourth-order method

)doc" )
            .value( "three_eight_rule_rk_4",
                    tni::threeEighthRuleRK4,
                    R"doc(

Coefficients for the classic Runge Kutta 3/8-rule fourth-order method

)doc" )
            .value( "heun_euler",
                    tni::heunEuler,
                    R"doc(

Coefficients for the Heun's method of order 2 with an embedded Euler method of order 1

)doc" )
            .value( "rkf_12",
                    tni::rungeKuttaFehlberg12,
                    R"doc(

Coefficients for the Runge-Kutta-Fehlberg method of order 2 with an embedded 1st order

)doc" )
            .value( "rkf_45",
                    tni::rungeKuttaFehlberg45,
                    R"doc(

Coefficients for the Runge-Kutta-Fehlberg method of order 5 with an embedded 4th order

)doc" )
            .value( "rkf_56",
                    tni::rungeKuttaFehlberg56,
                    R"doc(

Coefficients for the Runge-Kutta-Fehlberg method of order 6 with an embedded 5th order

)doc" )
            .value( "rkf_78",
                    tni::rungeKuttaFehlberg78,
                    R"doc(

Coefficients for the Runge-Kutta-Fehlberg method of order 8 with an embedded 7th order

)doc" )
            .value( "rkdp_87",
                    tni::rungeKutta87DormandPrince,
                    R"doc(

Coefficients for the Dormand-Prince method of order 7 with an embedded 8th order

)doc" )
            .value( "rkf_89",
                    tni::rungeKuttaFehlberg89,
                    R"doc(

Coefficients for the Runge-Kutta-Fehlberg method of order 9 with an embedded 8th order

)doc" )
            .value( "rkv_89",
                    tni::rungeKuttaVerner89,
                    R"doc(

Coefficients for the Runge-Kutta-Verner method of order 9 with an embedded 8th order

)doc" )
            .value( "rkf_108",
                    tni::rungeKuttaFeagin108,
                    R"doc(

Coefficients for the Runge-Kutta-Feagin method of order 8 with an embedded 10th order

)doc" )
            .value( "rkf_1210",
                    tni::rungeKuttaFeagin1210,
                    R"doc(

Coefficients for the Runge-Kutta-Feagin method of order 10 with an embedded 12ve order

)doc" )
            .value( "rkf_1412",
                    tni::rungeKuttaFeagin1412,
                    R"doc(

Coefficients for the Runge-Kutta-Feagin method of order 12 with an embedded 14th order

)doc" )
            .export_values( );

    py::enum_< tni::RungeKuttaCoefficients::OrderEstimateToIntegrate >( m,
                                                                        "OrderToIntegrate",
                                                                        R"doc(

Enumeration defining Which integrator order needs to be integrated, only used for coefficient sets with an embedded order.





      )doc" )
            .value( "lower",
                    tni::RungeKuttaCoefficients::OrderEstimateToIntegrate::lower,
                    R"doc(

For a method of order :math:`p`, with embedded method of order :math:`q`, the step is taken using the method with order :math:`\min(p,q)`

)doc" )
            .value( "higher",
                    tni::RungeKuttaCoefficients::OrderEstimateToIntegrate::higher,
                    R"doc(

For a method of order :math:`p`, with embedded method of order :math:`q`, the step is taken using the method with order :math:`\max(p,q)`

)doc" )
            .export_values( );

    py::enum_< tni::ExtrapolationMethodStepSequences >( m,
                                                        "ExtrapolationMethodStepSequences",
                                                        R"doc(

Enumeration of available extrapolation method substep sequences, with :math:`n_{j}` defining the number of substeps in iteration :math:`j`.





      )doc" )
            .value( "bulirsch_stoer_sequence",
                    tni::ExtrapolationMethodStepSequences::bulirsch_stoer_sequence,
                    R"doc(

Sequence for which :math:`n_{j}=2n_{j-2}` (2, 4, 6, 8, 12, 16, 24, ....)

)doc" )
            .value( "deufelhard_sequence",
                    tni::ExtrapolationMethodStepSequences::deufelhard_sequence,
                    R"doc(

Sequence for which :math:`n_{j}=2(j+1)` (2, 4, 6, 8, 10, 12, 14, ....)

)doc" )
            .export_values( );

    // CLASSES
    py::class_< tni::IntegratorSettings< TIME_TYPE >,
                std::shared_ptr< tni::IntegratorSettings< TIME_TYPE > > >( m,
                                                                           "IntegratorSettings",
                                                                           R"doc(

         Functional base class to define settings for integrators.

         Class to define settings for numerical integrators, for instance for use in numerical integration of equations of motion/
         variational equations. This class can be used for simple integrators such as fixed step RK and Euler. Integrators that
         require more settings to define have their own derived class.





      )doc" )
            .def_readwrite( "initial_time",
                            &tni::IntegratorSettings< TIME_TYPE >::initialTimeDeprecated_ );

    py::class_< tni::RungeKuttaFixedStepSizeSettings< TIME_TYPE >,
                std::shared_ptr< tni::RungeKuttaFixedStepSizeSettings< TIME_TYPE > >,
                tni::IntegratorSettings< TIME_TYPE > >( m,
                                                        "RungeKuttaFixedStepSizeSettings",
                                                        R"doc(

         `IntegratorSettings`-derived class to define settings for Runge Kutta integrators with a fixed step size





      )doc" );

    py::class_< tni::RungeKuttaVariableStepSizeBaseSettings< TIME_TYPE >,
                std::shared_ptr< tni::RungeKuttaVariableStepSizeBaseSettings< TIME_TYPE > >,
                tni::IntegratorSettings< TIME_TYPE > >(
            m, "RungeKuttaVariableStepSizeBaseSettings", R"doc(No documentation found.)doc" );

    py::class_<
            tni::RungeKuttaVariableStepSizeSettingsVectorTolerances< TIME_TYPE >,
            std::shared_ptr< tni::RungeKuttaVariableStepSizeSettingsVectorTolerances< TIME_TYPE > >,
            tni::RungeKuttaVariableStepSizeBaseSettings< TIME_TYPE > >(
            m,
            "RungeKuttaVariableStepSizeSettingsVectorTolerances",
            R"doc(No documentation found.)doc" );

    py::class_<
            tni::RungeKuttaVariableStepSizeSettingsScalarTolerances< TIME_TYPE >,
            std::shared_ptr< tni::RungeKuttaVariableStepSizeSettingsScalarTolerances< TIME_TYPE > >,
            tni::RungeKuttaVariableStepSizeBaseSettings< TIME_TYPE > >(
            m,
            "RungeKuttaVariableStepSizeSettingsScalarTolerances",
            R"doc(No documentation found.)doc" );

    py::class_< tni::BulirschStoerIntegratorSettings< TIME_TYPE >,
                std::shared_ptr< tni::BulirschStoerIntegratorSettings< TIME_TYPE > >,
                tni::IntegratorSettings< TIME_TYPE > >( m,
                                                        "BulirschStoerIntegratorSettings",
                                                        R"doc(

         `IntegratorSettings`-derived class to define settings for Bulirsch-Stoer integrator settings.





      )doc" );

    py::class_< tni::AdamsBashforthMoultonSettings< TIME_TYPE >,
                std::shared_ptr< tni::AdamsBashforthMoultonSettings< TIME_TYPE > >,
                tni::IntegratorSettings< TIME_TYPE > >( m,
                                                        "AdamsBashforthMoultonSettings",
                                                        R"doc(

         `IntegratorSettings`-derived class to define settings for Adams-Bashforth-Moulton integrator settings.





      )doc" );

    py::class_< tni::IntegratorStepSizeControlSettings,
                std::shared_ptr< tni::IntegratorStepSizeControlSettings > >(
            m,
            "IntegratorStepSizeControlSettings",
            R"doc(

         Base class to define settings for step-size control algorithm.

         Base class to define settings for step-size control algorithm, typically created by one of the functions provided in this module





      )doc" )
            .def_readwrite( "safety_factor",
                            &tni::IntegratorStepSizeControlSettings::safetyFactorForNextStepSize_ )
            .def_readwrite(
                    "minimum_step_decrease",
                    &tni::IntegratorStepSizeControlSettings::minimumFactorDecreaseForNextStepSize_ )
            .def_readwrite( "maximum_step_decrease",
                            &tni::IntegratorStepSizeControlSettings::
                                    maximumFactorDecreaseForNextStepSize_ );

    py::class_< tni::IntegratorStepSizeValidationSettings,
                std::shared_ptr< tni::IntegratorStepSizeValidationSettings > >(
            m,
            "IntegratorStepSizeValidationSettings",
            R"doc(

         Base class to define settings for step-size validation algorithm.

         Base class to define settings for step-size validation algorithm, typically created by one of the functions provided in this module





      )doc" )
            .def_readwrite( "minimum_step",
                            &tni::IntegratorStepSizeValidationSettings::minimumStep_ )
            .def_readwrite( "maximum_step",
                            &tni::IntegratorStepSizeValidationSettings::maximumStep_ )
            .def_readwrite( "minimum_step_handling",
                            &tni::IntegratorStepSizeValidationSettings::
                                    minimumIntegrationTimeStepHandling_ );

    // FACTORY FUNCTIONS
    m.def( "print_butcher_tableau",
           &tni::printButcherTableau,
           py::arg( "coefficient_set" ),
           R"doc(

 Print the Butcher tableau of a given coefficient set.


 Parameters
 ----------
 coefficient_set : CoefficientSets
     Coefficient set of which the Butcher tableau will be printed.





     )doc" );

    m.def( "step_size_validation",
           &tni::stepSizeValidationSettings,
           py::arg( "minimum_step" ),
           py::arg( "maximum_step" ),
           py::arg( "minimum_step_size_handling" ) = tni::throw_exception_below_minimum,
           py::arg( "accept_infinity_step" ) = false,
           py::arg( "accept_nan_step" ) = false,
           R"doc(

 Creates settings step size validation in a variable step-size integrator.

 Function to create settings step size validation in a variable step-size integrator. The validation
 model takes the proposed new step size  :math:`\Delta t_{rec}` as input, and checks if it meets predefined conditions, specifically
 whether the proposed time step falls in a given predefined range :math:`[\Delta t_{\min}, \Delta t_{\max}]`.
 This function also provides the option of handling recommended step sizes below :math:`\Delta t_{\min}` in various ways,
 and control on how to deal with recommend Inf/NaN step sizes.


 Parameters
 ----------
 minimum_step : float
     Value of minimum permitted time step :math:`\Delta t_{\min}`.
 maximum_step : float
     Value of maximum permitted time step :math:`\Delta t_{\max}`.
 minimum_step_size_handling : MinimumIntegrationTimeStepHandling, default = throw_exception_below_minimum
     Entry defining the behaviour when :math:`\Delta t_{rec}<\Delta t_{\min}`.
 accept_infinity_step : bool, default = False
     Entry defining whether to accept a step size of infinity (if False, exception is throw in such cases)
 accept_nan_step : bool, default = False
     Entry defining whether to accept a step size of NaN (if False, exception is throw in such cases)
 Returns
 -------
 IntegratorStepSizeValidationSettings
     Object containing settings for step-size validation.






     )doc" );

    m.def( "step_size_control_elementwise_scalar_tolerance",
           &tni::perElementIntegratorStepSizeControlSettings< double >,
           py::arg( "relative_error_tolerance" ),
           py::arg( "absolute_error_tolerance" ),
           py::arg( "safety_factor" ) = 0.8,
           py::arg( "minimum_factor_increase" ) = 0.1,
           py::arg( "maximum_factor_increase" ) = 4.0,
           R"doc(

 Creates settings for integrator step-size control, using element-wise analysis for the propagated states.

 Function to create settings for integrator step-size control, using element-wise analysis for the propagated states. For a propagated
 state :math:`\mathbf{x}` with entries :math:`x_{i}`, and an estimate :math:`\boldsymbol{\epsilon}` for the current local error, the following
 algorithm is performed per element :math:`i` to calculate the required error :math:`\epsilon_{i,req}` on this element:

 .. math::

    \epsilon_{i,req}=\epsilon_{r}x_{i}+\epsilon_{a}

 A proposed modification to the step size is then computed, using the most constraining of all state elements

 .. math::

    \bar{\Delta t_{rec.}}&=\Delta t\left(\min_{i}\left(\frac{\epsilon_{i,req}}{\epsilon_{i}}\right)\right)^{p}\\
    \Delta t_{rec.}&=K\bar{\Delta t_{rec.}}

 with :math:`p` the order of the local truncation error of the method for which step-size control is being applied,
 :math:`\Delta t_{rec.}` the new, recommended step size, and :math:`\Delta t` the current step size. The factor :math:`K` is a safety factor
 used make the time step slightly smaller than strictly required.

 A minimum and maximum change in time step may be provided by the user, such that if :math:`\Delta t_{rec.}/\Delta t` is too large or too small,
 the proposed increase/decrease to the step size is constrained to this limit value. That is, if :math:`\Delta t_{rec.}/\Delta t` proposed by the algorithm is
 1000, and the ``maximum_factor_increase`` input is equal to 20, the algorithm will use :math:`\Delta t_{rec.}/\Delta t=20` in what follows.

 For cases where :math:`\bar{\Delta t_{rec.}}/\Delta t < 1`, the step is recommended to be recomputed with the new proposed step size (e.g. the current step
 is not accepted, and will be re-attempted with a smaller step size). For cases where :math:`\bar{\Delta t_{rec.}}/\Delta t > 1`, the step is accepted, and
 the next step will be performed with the new, higher, step size.


 Parameters
 ----------
 relative_error_tolerance : float
     Value of relative error tolerance :math:`\epsilon_{r}`.
 absolute_error_tolerance : float
     Value of absolute error tolerance :math:`\epsilon_{a}`.
 safety_factor : float, default = 0.8
     Safety factor :math:`K` for step size control
 minimum_factor_increase : float, default = 0.1
     Minimum permissible value for :math:`\Delta t_{rec.}/\Delta t`
 maximum_factor_increase : float, default = 4.0
     Maximum permissible value for :math:`\Delta t_{rec.}/\Delta t`
 Returns
 -------
 IntegratorStepSizeControlSettings
     Object containing settings for per-element step-size control.






     )doc" );

    m.def( "step_size_control_elementwise_matrix_tolerance",
           &tni::perElementIntegratorStepSizeControlSettings< Eigen::MatrixXd >,
           py::arg( "relative_error_tolerance" ),
           py::arg( "absolute_error_tolerance" ),
           py::arg( "safety_factor" ) = 0.8,
           py::arg( "minimum_factor_increase" ) = 0.1,
           py::arg( "maximum_factor_increase" ) = 4.0,
           R"doc(

 Creates settings for integrator step-size control, using element-wise analysis for the propagated states.

 Function to create settings for integrator step-size control, using element-wise analysis for the propagated states. This function
 is similar to :func:`~tudatpy.numerical_simulation.propagation_setup.integrator.step_size_control_elementwise_scalar_tolerance`,
 with the differences that the tolerances are provided as a vector/matrix (which must be of equal size as the propagates state), such that
 different tolerances can be provided for each state element. The behaviour of the algorithm is then such that
 :math:`\epsilon_{r}\rightarrow\epsilon_{r,i}` and :math:`\epsilon_{a}\rightarrow\epsilon_{a,i}`.

 If the size of the tolerances used as input differ from one another, or differ from the size of the state vector, an exception is thrown


 Parameters
 ----------
 relative_error_tolerance : numpy.ndarray[numpy.float64[m, n]]
     Values of relative error tolerance :math:`\boldsymbol{\epsilon}_{r}`.
 absolute_error_tolerance : numpy.ndarray[numpy.float64[m, n]]
     Values of absolute error tolerance :math:`\boldsymbol{\epsilon}_{a}`.
 safety_factor : float, default = 0.8
     Safety factor :math:`K` for step size control
 minimum_factor_increase : float, default = 0.1
     Minimum permissible value for :math:`\Delta t_{rec.}/\Delta t`
 maximum_factor_increase : float, default = 4.0
     Maximum permissible value for :math:`\Delta t_{rec.}/\Delta t`
 Returns
 -------
 IntegratorStepSizeControlSettings
     Object containing settings for per-element step-size control.






     )doc" );

    m.def( "step_size_control_blockwise_scalar_tolerance",
           &tni::perBlockIntegratorStepSizeControlSettings< double >,
           py::arg( "block_indices" ),
           py::arg( "relative_error_tolerance" ),
           py::arg( "absolute_error_tolerance" ),
           py::arg( "safety_factor" ) = 0.8,
           py::arg( "minimum_factor_increase" ) = 0.1,
           py::arg( "maximum_factor_increase" ) = 4.0,
           R"doc(

 Creates settings for integrator step-size control, using block-wise analysis for the propagated states.

 Function to create settings for integrator step-size control, using block-wise analysis for the propagated states. This function
 is similar to :func:`~tudatpy.numerical_simulation.propagation_setup.integrator.step_size_control_elementwise_scalar_tolerance`,
 with the difference that the error estimation :math:`\boldsymbol{\epsilon}` is not used on an element-by-element basis, but using the norms
 of user defined matrix blocks. This is for instance very useful when propagating Cartesian states, where the tolerances are then
 typically applied twice: once to the norm of the position error, and once to the norm of the velocity error.

 The algorithm is then run, using the modification
 that :math:`\epsilon_{i}\rightarrow||\boldsymbol{\epsilon_{[i,k],[j,l]}}||`. Where the indices on the right-hand side denote start row :math:`i`,
 start column :math:`j`, number of rows :math:`k` and number of columns :math:`l`. over
 which the state error norm is to be taken. For a single Cartesian state vector, the norm is taken on blocks :math:`[0,3],[0,1]` and :math:`[3,3],[0,1]`

 .. note::

     If you would like to create block indices that group the position and velocity elements, take a look at the :func:`~tudatpy.numerical_simulation.propagation_setup.integrator.standard_cartesian_state_element_blocks` function.
     The function will create a list of blocks that can be used as input to the `block_indices` argument of this function.

 Parameters
 ----------
 block_indices : list[tuple[int,int,int,int]]
     List of matrix blocks over which the norms are to be taken (with entries of the tuple denoting :math:`i,j,k,l`, in order.
 relative_error_tolerance : float
     Value of relative error tolerance :math:`\epsilon_{r}`.
 absolute_error_tolerance : float
     Value of absolute error tolerance :math:`\epsilon_{a}`.
 safety_factor : float, default = 0.8
     Safety factor :math:`K` for step size control
 minimum_factor_increase : float, default = 0.1
     Minimum permissible value for :math:`\Delta t_{rec.}/\Delta t`
 maximum_factor_increase : float, default = 4.0
     Maximum permissible value for :math:`\Delta t_{rec.}/\Delta t`
 Returns
 -------
 IntegratorStepSizeControlSettings
     Object containing settings for per-element step-size control.


 Examples
 --------
 In this example, step size control settings are created for a Cartesian state vector, which group the position and velocity elements for the step size validation. Note, these block indices can also be conveniently created using the :func:`~tudatpy.numerical_simulation.propagation_setup.integrator.standard_cartesian_state_element_blocks` function.
 Here we will create them manually for demonstration purposes.
 We would like to create integrator settings with the following settings:
 
 - Relative error tolerance of 1e-12
 - Absolute error tolerance of 1e-12
 - Relative and absolute error tolerances are applied to the position and velocity blocks of a single Cartesian state vector
 
 .. code-block:: python

     from tudatpy.numerical_simulation import propagation_setup
     ...
     # Define integrator step settings
     initial_time_step = 10.0
     minimum_step_size = 1.0e-12
     maximum_step_size = 60.0
     relative_tolerance = 1.0e-12
     absolute_tolerance = 1.0e-12
     """
     Our Cartesian state y is a 2D vector, with dimensions [6, 1].
     # column-index
            0
     y = [[ x ], # row-index 0
          [ y ], # row-index 1
          [ z ], # row-index 2
          [ vx ], # row-index 3
          [ vy ], # row-index 4
          [ vz ]] # row-index 5
     The block indices are denoted as i, j, k, l, where:
     - i: start row index of the block
     - j: start column index of the block
     - k: number of rows in the block
     - l: number of columns in the block
     The corresponding block indices for the position block are therefore
     - i = 0 (start the block at row-index 0)
     - j = 0 (start the block at column-index 0)
     - k = 3 (the block has 3 rows: x, y, z)
     - l = 1 (the block has 1 column)
     which gives us the block indices (i=0, j=0, k=3, l=1).

     For the velocity block, the indices are:
     - i = 3 (start the block at row-index 3)
     - j = 0 (start the block at column-index 0)
     - k = 3 (the block has 3 rows: vx, vy, vz)
     - l = 1 (the block has 1 column)
     which gives us the block indices (i=3, j=0, k=3, l=1).
     """
     # Manually define block indices for position and velocity,
     # which is equivalent to the standard block indices using:
     # block_indices = propagation_setup.integrator.standard_cartesian_state_element_blocks(6, 1)
     block_indices = [(0, 0, 3, 1), (3, 0, 3, 1)]
     step_size_control_settings = (
         propagation_setup.integrator.step_size_control_blockwise_scalar_tolerance(
             block_indices, relative_tolerance, absolute_tolerance
         )
     )
     step_size_validation_settings = propagation_setup.integrator.step_size_validation(
         minimum_step=minimum_step_size, maximum_step=maximum_step_size
     )
     # Retrieve coefficient set
     coefficient_set = propagation_setup.integrator.rkf_78
     variable_step_integrator_settings = (
         propagation_setup.integrator.runge_kutta_variable_step(
             initial_time_step,
             coefficient_set,
             step_size_control_settings,
             step_size_validation_settings,
         )
     )
     ...


     )doc" );

    m.def( "step_size_control_blockwise_matrix_tolerance",
           &tni::perBlockIntegratorStepSizeControlSettings< Eigen::MatrixXd >,
           py::arg( "block_indices" ),
           py::arg( "relative_error_tolerance" ),
           py::arg( "absolute_error_tolerance" ),
           py::arg( "safety_factor" ) = 0.8,
           py::arg( "minimum_factor_increase" ) = 0.1,
           py::arg( "maximum_factor_increase" ) = 4.0,
           R"doc(

 Creates settings for integrator step-size control, using block-wise analysis for the propagated states.

 Function to create settings for integrator step-size control, using block-wise analysis for the propagated states. This function
 is similar to :func:`~tudatpy.numerical_simulation.propagation_setup.integrator.step_size_control_blockwise_scalar_tolerance`,
 with the differences that the tolerances are provided as a list (which must be of equal size as the number of state blocks used), such that
 different tolerances can be provided for each state block.

 If the size of the tolerances used as input differ from one another, or differ from the number of blocks, an exception is thrown


 .. note::

    If you would like to create block indices that group the position and velocity elements, take a look at the :func:`~tudatpy.numerical_simulation.propagation_setup.integrator.standard_cartesian_state_element_blocks` function.
    The function will create a list of blocks that can be used as input to the `block_indices` argument of this function.

 Parameters
 ----------
 block_indices : list[tuple[int,int,int,int]]
     List of matrix blocks over which the norms are to be taken (with entries of the tuple denoting :math:`i,j,k,l`, in order.
 relative_error_tolerance : numpy.ndarray[numpy.float64[m, 1]]
     Values of relative error tolerance :math:`\boldsymbol{\epsilon}_{r}`.
 absolute_error_tolerance : numpy.ndarray[numpy.float64[m, 1]]
     Values of absolute error tolerance :math:`\boldsymbol{\epsilon}_{a}`.
 safety_factor : float, default = 0.8
     Safety factor :math:`K` for step size control
 minimum_factor_increase : float, default = 0.1
     Minimum permissible value for :math:`\Delta t_{rec.}/\Delta t`
 maximum_factor_increase : float, default = 4.0
     Maximum permissible value for :math:`\Delta t_{rec.}/\Delta t`
 Returns
 -------
 IntegratorStepSizeControlSettings
     Object containing settings for per-element step-size control.

 Examples
 --------
 In this example, step size control settings are created for a Cartesian state vector, which group the position and velocity elements for the step size validation. Note, these block indices can also be conveniently created using the :func:`~tudatpy.numerical_simulation.propagation_setup.integrator.standard_cartesian_state_element_blocks` function.
 Here we will create them manually for demonstration purposes.
 We would like to create integrator settings with the following settings:

 - Relative error tolerance of 1e-12 for position and velocity blocks
 - Absolute error tolerance of 1e-9 for position and 1e-12 for velocity blocks
 
 .. code-block:: python

     import numpy as np
     from tudatpy.numerical_simulation import propagation_setup
     ...
     # Define integrator step settings
     initial_time_step = 10.0
     minimum_step_size = 1.0e-12
     maximum_step_size = 60.0
     relative_tolerance_pos = 1.0e-12
     relative_tolerance_vel = 1.0e-12
     absolute_tolerance_pos = 1.0e-9
     absolute_tolerance_vel = 1.0e-12
     """
     Our Cartesian state y is a 2D vector, with dimensions [6, 1].
     # column-index
         0
     y = [[ x ], # row-index 0
         [ y ], # row-index 1
         [ z ], # row-index 2
         [ vx ], # row-index 3
         [ vy ], # row-index 4
         [ vz ]] # row-index 5
     The block indices are denoted as i, j, k, l, where:
     - i: start row index of the block
     - j: start column index of the block
     - k: number of rows in the block
     - l: number of columns in the block
     The corresponding block indices for the position block are therefore
     - i = 0 (start the block at row-index 0)
     - j = 0 (start the block at column-index 0)
     - k = 3 (the block has 3 rows: x, y, z)
     - l = 1 (the block has 1 column)
     which gives us the block indices (i=0, j=0, k=3, l=1).
     For the velocity block, the indices are:
     - i = 3 (start the block at row-index 3)
     - j = 0 (start the block at column-index 0)
     - k = 3 (the block has 3 rows: vx, vy, vz)
     - l = 1 (the block has 1 column)
     which gives us the block indices (i=3, j=0, k=3, l=1).
     """
     # Manually define block indices for position and velocity,
     # which is equivalent to the standard block indices using:
     # block_indices = propagation_setup.integrator.standard_cartesian_state_element_blocks(6, 1)
     block_indices = [(0, 0, 3, 1), (3, 0, 3, 1)]
     # Different from the scalar tolerance, the matrix tolerance is defined as
     # the relative and absolute tolerances for each block.
     relative_tolerances = np.array([[relative_tolerance_pos], [relative_tolerance_vel]])
     absolute_tolerances = np.array([[absolute_tolerance_pos], [absolute_tolerance_vel]])
     step_size_control_settings = (
         propagation_setup.integrator.step_size_control_blockwise_matrix_tolerance(
             block_indices, relative_tolerances, absolute_tolerances
         )
     )
     step_size_validation_settings = propagation_setup.integrator.step_size_validation(
         minimum_step=minimum_step_size, maximum_step=maximum_step_size
     )
     # Retrieve coefficient set
     coefficient_set = propagation_setup.integrator.rkf_78
     variable_step_integrator_settings = (
         propagation_setup.integrator.runge_kutta_variable_step(
             initial_time_step,
             coefficient_set,
             step_size_control_settings,
             step_size_validation_settings,
         )
     )
     ...




     )doc" );

    m.def( "standard_cartesian_state_element_blocks",
           &tni::getStandardCartesianStatesElementsToCheck,
           py::arg( "number_of_rows" ),
           py::arg( "number_of_columns" ),
           R"doc(

 Function to generate step size control blocks on position and velocity elements for numerical integration

 Function to generate step size control blocks on position and velocity elements for numerical integration, typically provided
 to the :func:`~tudatpy.numerical_simulation.propagation_setup.integrator.step_size_control_custom_blockwise_scalar_tolerance` or
 :func:`~tudatpy.numerical_simulation.propagation_setup.integrator.step_size_control_custom_blockwise_matrix_tolerance` function.
 By providing this function to one of these step-size control functions, the final column of the state vector is taken (such that
 it works  both for state-only, and variational equations and state propagation) and combined into :math:`N` blocks of size 3.
 The step-size control is then done on each of these blocks, which will represent the position and velocity blocks.


 Parameters
 ----------
 number_of_rows : int
     Number of rows in state vector
 number_of_columns : int
     Number of columns in state vector
 Returns
 -------
 list[tuple[int,int,int,int]]
     List of matrix blocks over which the step size control is to be done (see ``block_indices_function`` input to :func:`~tudatpy.numerical_simulation.propagation_setup.integrator.step_size_control_custom_blockwise_scalar_tolerance`)






     )doc" );

    m.def( "standard_rotational_state_element_blocks",
           &tni::getStandardRotationalStatesElementsToCheck,
           py::arg( "number_of_rows" ),
           py::arg( "number_of_columns" ) );

    m.def( "step_size_control_custom_blockwise_scalar_tolerance",
           &tni::perBlockFromFunctionIntegratorStepSizeControlSettings< double >,
           py::arg( "block_indices_function" ),
           py::arg( "relative_error_tolerance" ),
           py::arg( "absolute_error_tolerance" ),
           py::arg( "safety_factor" ) = 0.8,
           py::arg( "minimum_factor_increase" ) = 0.1,
           py::arg( "maximum_factor_increase" ) = 4.0,
           R"doc(

 Creates settings for integrator step-size control, using block-wise analysis for the propagated states.

 Function to create settings for integrator step-size control, using block-wise analysis for the propagated states. This function
 is similar to :func:`~tudatpy.numerical_simulation.propagation_setup.integrator.step_size_control_blockwise_scalar_tolerance`,
 but rather than providing the ``block_indices`` directly, a function to determine the block indices, based on the size of the
 propagated state, is provided. For instance, the :func:`~tudatpy.numerical_simulation.propagation_setup.integrator.standard_cartesian_state_element_blocks`
 can be provided to this function (as ``block_indices_function``), which will adapt the block indices depending on the size of the propagated state
 (e.g. regardless of how many bodies are propagated, step size control will always be done on position and velocity element blocks)




 Parameters
 ----------
 block_indices_function : Callable[[int,int],list[tuple[int,int,int,int]]]
     Function returning list of matrix blocks over which the norms are to be taken (with entries of the tuple denoting :math:`i,j,k,l`, in order, with number of rows and columns of propagated state as input.
 relative_error_tolerance : float
     Value of relative error tolerance :math:`\epsilon_{r}`.
 absolute_error_tolerance : float
     Value of absolute error tolerance :math:`\epsilon_{a}`.
 safety_factor : float, default = 0.8
     Safety factor :math:`K` for step size control
 minimum_factor_increase : float, default = 0.1
     Minimum permissible value for :math:`\Delta t_{rec.}/\Delta t`
 maximum_factor_increase : float, default = 4.0
     Maximum permissible value for :math:`\Delta t_{rec.}/\Delta t`
 Returns
 -------
 IntegratorStepSizeControlSettings
     Object containing settings for per-element step-size control.






     )doc" );

    m.def( "step_size_control_custom_blockwise_matrix_tolerance",
           &tni::perBlockFromFunctionIntegratorStepSizeControlSettings< Eigen::MatrixXd >,
           py::arg( "block_indices_function" ),
           py::arg( "relative_error_tolerance" ),
           py::arg( "absolute_error_tolerance" ),
           py::arg( "safety_factor" ) = 0.8,
           py::arg( "minimum_factor_increase" ) = 0.1,
           py::arg( "maximum_factor_increase" ) = 4.0,
           R"doc(

 Creates settings for integrator step-size control, using block-wise analysis for the propagated states.

 Function to create settings for integrator step-size control, using block-wise analysis for the propagated states. This function
 is similar to :func:`~tudatpy.numerical_simulation.propagation_setup.integrator.step_size_control_custom_blockwise_scalar_tolerance`,
 but uses blockwise tolerances (as in :func:`~tudatpy.numerical_simulation.propagation_setup.integrator.step_size_control_blockwise_matrix_tolerance`)



 Parameters
 ----------
 block_indices_function : Callable[[int,int],list[tuple[int,int,int,int]]]
     Function returning list of matrix blocks over which the norms are to be taken (with entries of the tuple denoting :math:`i,j,k,l`, in order, with number of rows and columns of propagated state as input.
 relative_error_tolerance : numpy.ndarray[numpy.float64[m, 1]]
     Values of relative error tolerance :math:`\boldsymbol{\epsilon}_{r}`.
 absolute_error_tolerance : numpy.ndarray[numpy.float64[m, 1]]
     Values of absolute error tolerance :math:`\boldsymbol{\epsilon}_{a}`.
 safety_factor : float, default = 0.8
     Safety factor :math:`K` for step size control
 minimum_factor_increase : float, default = 0.1
     Minimum permissible value for :math:`\Delta t_{rec.}/\Delta t`
 maximum_factor_increase : float, default = 4.0
     Maximum permissible value for :math:`\Delta t_{rec.}/\Delta t`
 Returns
 -------
 IntegratorStepSizeControlSettings
     Object containing settings for per-element step-size control.






     )doc" );

    m.def( "runge_kutta_fixed_step",
           &tni::rungeKuttaFixedStepSettings< TIME_TYPE >,
           py::arg( "time_step" ),
           py::arg( "coefficient_set" ),
           py::arg( "order_to_use" ) = tni::RungeKuttaCoefficients::OrderEstimateToIntegrate::lower,
           py::arg( "assess_termination_on_minor_steps" ) = false,
           R"doc(

 Creates the settings for the Runge-Kutta fixed step size integrator.

 Function to create settings for the Runge-Kutta integrator with a constant step size.
 Different coefficient sets (Butcher's tableau) can be used (see the `CoefficientSets` enum).


 Parameters
 ----------
 time_step : float
     Initial time step to be used.
 coefficient_set : CoefficientSets
     Coefficient set (Butcher's tableau) to be used in the integration.
 order_to_use : OrderToIntegrate, default=OrderToIntegrate.lower
     If the coefficient set is supposed to be for variable step sizes (with an embedded method of a different order),
     this parameter can be used to set the order that will be used.

 assess_termination_on_minor_steps : bool, default=false
     Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the
     integrator (true) or only at the end of each integration step (false).

 Returns
 -------
 IntegratorSettings
     Object containing settings for the integrator.





 Examples
 --------
 In this example, settings for the classical RK4 integrator with 30 second time step are created

 .. code-block:: python

   # Create RK4 settings
   integrator_settings = integrator.runge_kutta_fixed_step(
       time_step = 30.0,
       coefficient_set = integrator.CoefficientSets.rk_4 )

 In this example, settings for fixed-step integration using the higher-order (8th-order) of the two
 embedded propagators of the RKF7(8) method are created, with a time-step of 120 seconds.

 .. code-block:: python

   # Create 8th order RKF settings
   integrator_settings = integrator.runge_kutta_fixed_step(
       time_step = 120.0,
       coefficient_set = integrator.CoefficientSets.rkf_78,
       order_to_use = integrator.OrderToIntegrate.higher )



     )doc" );

    m.def( "runge_kutta_variable_step",
           &tni::multiStageVariableStepSizeSettings< TIME_TYPE >,
           py::arg( "initial_time_step" ),
           py::arg( "coefficient_set" ),
           py::arg( "step_size_control_settings" ),
           py::arg( "step_size_validation_settings" ),
           py::arg( "assess_termination_on_minor_steps" ) = false,
           R"doc(

 Creates the settings for the Runge-Kutta variable step size integrator.

 Function to create settings for the Runge-Kutta variable step size integrator.
 Different coefficient sets (Butcher's tableau) can be used (see the :class:`CoefficientSets` enum).
 The step-size control algorithm is defined by a :class:`~tudatpy.numerical_simulation.propagation_setup.integrator.IntegratorStepSizeControlSettings` and
 :class:`~tudatpy.numerical_simulation.propagation_setup.integrator.IntegratorStepSizeValidationSettings` object, created using one of the functions
 listed above.


 Parameters
 ----------
 initial_time_step : float
     Initial time step to be used.
 coefficient_set : CoefficientSets
     Coefficient set (Butcher's tableau) to be used in the integration.
 step_size_control_settings : IntegratorStepSizeControlSettings
     Object used to control the step size, by computing a new step size :math:`\Delta t_{rec.}`, from the embedded Runge-Kutta integrator pair,
     and recommending whether the steps is to be accepted, or recomputed with a different time step.

 step_size_validation_settings : IntegratorStepSizeValidationSettings
     Object used to validate whether the :math:`\Delta t_{rec.}` provided by model defined by the ``step_size_control_settings`` meets with user-defined
     criteria (minimum, maximum values, etc.)

 assess_termination_on_minor_steps : bool, default=false
     Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the
     integrator (true) or only at the end of each integration step (false).

 Returns
 -------
 IntegratorSettings
     Object containing settings for the integrator.





 Examples
 --------
 In this example, settings for the variable step RK4(5) integrator are created, with the same tolerances (:math:`10^{-10}`)
 applied element-wise on the propagated state. The minimum and maximum time steps are set to 0.001 and 1000 seconds,
 and the initial step is set to 30 seconds. All other inputs are left on their defaults

 .. code-block:: python

   # Create RK4(5) settings
   control_settings = integrator.step_size_control_elementwise_scalar_tolerance( 1.0E-10, 1.0E-10 )
   validation_settings = integrator.step_size_validation( 0.001, 1000.0 )
   integrator_settings = integrator.runge_kutta_variable_step(
       initial_time_step = 30.0,
       coefficient_set = integrator.CoefficientSets.rkf_45,
       step_size_control_settings = control_settings,
       step_size_validation_settings = validation_settings )

 In this example, the above is modified such that step-size control is applied on position and velocity
 element blocks.

 .. code-block:: python

   # Create RK4(5) settings
   control_settings = integrator.step_size_control_custom_blockwise_scalar_tolerance(
       integrator.standard_cartesian_state_element_blocks
       1.0E-10, 1.0E-10 )
   validation_settings = integrator.step_size_validation( 0.001, 1000.0 )
   integrator_settings = integrator.runge_kutta_variable_step(
       initial_time_step = 30.0,
       coefficient_set = integrator.CoefficientSets.rkf_45,
       step_size_control_settings = control_settings,
       step_size_validation_settings = validation_settings )


     )doc" );

    m.def( "bulirsch_stoer_variable_step",
           &tni::bulirschStoerVariableStepIntegratorSettings< TIME_TYPE >,
           py::arg( "initial_time_step" ),
           py::arg( "extrapolation_sequence" ),
           py::arg( "maximum_number_of_steps" ),
           py::arg( "step_size_control_settings" ),
           py::arg( "step_size_validation_settings" ),
           py::arg( "assess_termination_on_minor_steps" ) = false,
           R"doc(

 Creates the settings for the variable time-step Bulirsch-Stoer integrator.

 Function to create settings for the variable time-step Bulirsch-Stoer integrator. This integrator
 works by performing the same (typically very large) step multiple times, using an ever increasing number of substeps.
 Each substep is performed using the modified midpoint method. The successive integrations from :math:`t_{i}` to :math:`t_{i}+\Delta t`
 are (in principle) done using ever increasing accuracy, as the size of the substep decreases. This integrator works
 by extrapolating the behaviour to a substep length of 0 (e.g. an infinite number of substeps), at which the solution should be perfect.
 The number of substeps on the :math:`i^{t}` iteration are done using the number of substeps defined by  entry :math:`i` of the
 ``extrapolation_sequence`` input. The number of iterations for a single step is defined by the ``maximum_number_of_steps`` entry.
 For instance, using the ``bulirsch_stoer_sequence`` sequence, and 5 iterations, the same step is done using 2, 4, 6, 8 and 12 substeps,
 and the results are then extrapolated to an infinite number of steps. Different extrapolation sequences can be used (see the `ExtrapolationMethodStepSequences` enum).

 The step-size control algorithm is defined by a :class:`~tudatpy.numerical_simulation.propagation_setup.integrator.IntegratorStepSizeControlSettings` and
 :class:`~tudatpy.numerical_simulation.propagation_setup.integrator.IntegratorStepSizeValidationSettings` object, created using one of the functions
 listed above. The time step control uses the result from the final, and second to final iteration to generate an error estimate of the current step.


 Parameters
 ----------
 time_step : float
     Initial time step to be used.
 extrapolation_sequence : ExtrapolationMethodStepSequences
     Extrapolation sequence to be used for the integration (defining the number of substeps in iteration :math:`i`).
 maximum_number_of_steps : int
     Number of entries from the sequence to be used (e.g., total number of iterations used for a single extrapolation and time step).
 step_size_control_settings : IntegratorStepSizeControlSettings
     Object used to control the step size, by computing a new step size :math:`\Delta t_{rec.}`, from the embedded Runge-Kutta integrator pair,
     and recommending whether the steps is to be accepted, or recomputed with a different time step.

 step_size_validation_settings : IntegratorStepSizeValidationSettings
     Object used to validate whether the :math:`\Delta t_{rec.}` provided by model defined by the ``step_size_control_settings`` meets with user-defined
     criteria (minimum, maximum values, etc.)

 assess_termination_on_minor_steps : bool, default=false
     Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the integrator (true) or only at the end of each integration step (false).
 Returns
 -------
 IntegratorSettings
     Object containing settings for the integrator.





 Examples
 --------
 In this example, settings for the Bulirsch-Stoer integrator with 600 second initial time step are created, using the typical
 sequence, using 6 iterations of the same step. By using the Bulirsch-Stoer sequence, this means that the same step is done
 using 2, 4, 6, 8, 12 and 16 substeps. The same tolerances (:math:`10^{-10}`)
 applied element-wise on the propagated state. The minimum and maximum time steps are set to 0.1 and 10000 seconds,
 and the initial step is set to 600 seconds. All other inputs are left on their defaults

 .. code-block:: python

   # Create BS settings
   control_settings = integrator.step_size_control_elementwise_scalar_tolerance( 1.0E-10, 1.0E-10 )
   validation_settings = integrator.step_size_validation( 0.1, 10000.0 )
   integrator_settings = integrator.bulirsch_stoer_variable_step(
       initial_time_step = 600.0,
       extrapolation_sequence = integrator.bulirsch_stoer_sequence,
       maximum_number_of_steps = 6,
       step_size_control_settings = control_settings,
       step_size_validation_settings = validation_settings )


     )doc" );

    m.def( "bulirsch_stoer_fixed_step",
           &tni::bulirschStoerFixedStepIntegratorSettings< TIME_TYPE >,
           py::arg( "time_step" ),
           py::arg( "extrapolation_sequence" ),
           py::arg( "maximum_number_of_steps" ),
           py::arg( "assess_termination_on_minor_steps" ) = false,
           R"doc(

 Creates the settings for the fixed time-step Bulirsch-Stoer integrator.

 Function to create settings for the fixed time-step Bulirsch-Stoer integrator. The
 underlying method is the same as :func:`~tudatpy.numerical_simulation.propagation_setup.integrator.bulirsch_stoer_variable_step`,
 but using a fixed, user-defined, time step.


 Parameters
 ----------
 time_step : float
     Time step to be used.
 extrapolation_sequence : ExtrapolationMethodStepSequences
     Extrapolation sequence to be used for the integration (defining the number of substeps in iteration :math:`i`).
 maximum_number_of_steps : int
     Number of entries from the sequence to be used (e.g., total number of iterations used for a single extrapolation and time step).
 assess_termination_on_minor_steps : bool, default=false
     Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the integrator (true) or only at the end of each integration step (false).
 Returns
 -------
 IntegratorSettings
     Object containing settings for the integrator.





 Examples
 --------
 In this example, settings for the Bulirsch-Stoer integrator with 600 second time step are created, using the typical
 sequence, using 6 iterations of the same step. By using the Bulirsch-Stoer sequence, this means that the same step is done
 using 2, 4, 6, 8, 12 and 16 substeps

 .. code-block:: python

   # Create BS settings
   integrator_settings = integrator.bulirsch_stoer_fixed_step(
       time_step = 300.0,
       extrapolation_sequence = integrator.bulirsch_stoer_sequence,
       maximum_number_of_steps = 6 )


     )doc" );

    m.def( "adams_bashforth_moulton",
           &tni::adamsBashforthMoultonSettings< TIME_TYPE >,
           py::arg( "initial_time_step" ),
           py::arg( "minimum_step_size" ),
           py::arg( "maximum_step_size" ),
           py::arg( "relative_error_tolerance" ) = 1.0E-12,
           py::arg( "absolute_error_tolerance" ) = 1.0E-12,
           py::arg( "minimum_order" ) = 6,
           py::arg( "maximum_order" ) = 11,
           py::arg( "assess_termination_on_minor_steps" ) = false,
           py::arg( "bandwidth" ) = 200.0,
           R"doc(

Creates the settings for the Adams-Bashforth-Moulton integrator.

Function to create settings for the Adams-Bashforth-Moulton multistep integrator.
For this integrator, the step size and order are both according to a control algorithm
similar to :func:`~tudatpy.numerical_simulation.propagation_setup.integrator.step_size_control_elementwise_scalar_tolerance`.
The integrator is initialized using an RKF7(8) integrator.

NOTE: this integrator's step-size and order control algorithm work in a method that is overly simplistic,
when increasing/decreasing the order, existing function evaluations are re-used, without any recomputations.
Similarly, when halving or doubling the time-step, the existing interpolating polynomial is evaluated at the relevant points.
This can lead to unwanted behaviour, where the time-step reduces to unrealistically low values. It is strongly
recommended that a reasonable minimum step is provided to this function, to partially mitigate this behaviour.

Parameters
----------
initial_time_step : float
    Initial time step to be used.
minimum_step_size : float
    Minimum time step to be used during the integration.
maximum_step_size : float
    Maximum time step to be used during the integration.
relative_error_tolerance : float, default=1.0E-12
    Relative tolerance to adjust the time step.
absolute_error_tolerance : float, default=1.0E-12
    Relative tolerance to adjust the time step.
minimum_order
    Minimum order of the integrator.
maximum_order
    Maximum order of the integrator.
assess_termination_on_minor_steps : bool, default=false
    Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the integrator (true) or only at the end of each integration step (false).
bandwidth : float, default=200.0
    Maximum error factor for doubling the step size.
Returns
-------
IntegratorSettings
    Object containing settings for the integrator.






     )doc" );

    m.def( "adams_bashforth_moulton_fixed_order",
           &tni::adamsBashforthMoultonSettingsFixedOrder< TIME_TYPE >,
           py::arg( "initial_time_step" ),
           py::arg( "minimum_step_size" ),
           py::arg( "maximum_step_size" ),
           py::arg( "relative_error_tolerance" ) = 1.0E-12,
           py::arg( "absolute_error_tolerance" ) = 1.0E-12,
           py::arg( "order" ) = 6,
           py::arg( "assess_termination_on_minor_steps" ) = false,
           py::arg( "bandwidth" ) = 200.0,
           R"doc(

Creates the settings for the Adams-Bashforth-Moulton integrator of fixed order.

Same as :func:`~tudatpy.numerical_simulation.propagation_setup.integrator.adams_bashforth_moulton`, but
with fixed order and variable step

Parameters
----------
initial_time_step : float
    Initial time step to be used.
minimum_step_size : float
    Minimum time step to be used during the integration.
maximum_step_size : float
    Maximum time step to be used during the integration.
relative_error_tolerance : float, default=1.0E-12
    Relative tolerance to adjust the time step.
absolute_error_tolerance : float, default=1.0E-12
    Relative tolerance to adjust the time step.
order
    Order of the integrator.
assess_termination_on_minor_steps : bool, default=false
    Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the integrator (true) or only at the end of each integration step (false).
bandwidth : float, default=200.0
    Maximum error factor for doubling the step size.
Returns
-------
IntegratorSettings
    Object containing settings for the integrator.






     )doc" );

    m.def( "adams_bashforth_moulton_fixed_step",
           &tni::adamsBashforthMoultonSettingsFixedStep< TIME_TYPE >,
           py::arg( "time_step" ),
           py::arg( "relative_error_tolerance" ) = 1.0E-12,
           py::arg( "absolute_error_tolerance" ) = 1.0E-12,
           py::arg( "minimum_order" ) = 6,
           py::arg( "maximum_order" ) = 11,
           py::arg( "assess_termination_on_minor_steps" ) = false,
           py::arg( "bandwidth" ) = 200.0,
           R"doc(

 Creates the settings for the Adams-Bashforth-Moulton fixed-step integrator.

 Same as :func:`~tudatpy.numerical_simulation.propagation_setup.integrator.adams_bashforth_moulton`, but
 with fixed step and variable order


 Parameters
 ----------
 time_step : float
     Initial time step to be used.
 relative_error_tolerance : float, default=1.0E-12
     Relative tolerance to adjust the time step.
 absolute_error_tolerance : float, default=1.0E-12
     Relative tolerance to adjust the time step.
 minimum_order
     Minimum order of the integrator.
 maximum_order
     Maximum order of the integrator.
 assess_termination_on_minor_steps : bool, default=false
     Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the integrator (true) or only at the end of each integration step (false).
 bandwidth : float, default=200.0
     Maximum error factor for doubling the step size.
 Returns
 -------
 IntegratorSettings
     Object containing settings for the integrator.






     )doc" );

    m.def( "adams_bashforth_moulton_fixed_step_fixed_order",
           &tni::adamsBashforthMoultonSettingsFixedStepFixedOrder< TIME_TYPE >,
           py::arg( "time_step" ),
           py::arg( "order" ) = 6,
           py::arg( "assess_termination_on_minor_steps" ) = false,
           R"doc(

 Creates the settings for the Adams-Bashforth-Moulton fixed-step, fixed-order integrator.

 Same as :func:`~tudatpy.numerical_simulation.propagation_setup.integrator.adams_bashforth_moulton`, but
 with fixed step and fixed order


 Parameters
 ----------
 time_step : float
     Initial time step to be used.
 order
     Order of the integrator.
 assess_termination_on_minor_steps : bool, default=false
     Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the integrator (true) or only at the end of each integration step (false).
 Returns
 -------
 IntegratorSettings
     Object containing settings for the integrator.






     )doc" );

    /*!
     * DEPRECATED
     * -------------------------------------------------------------------------
     *
     */

    m.def( "runge_kutta_variable_step_size_vector_tolerances",
           &tni::rungeKuttaVariableStepSettingsVectorTolerances< TIME_TYPE >,
           py::arg( "initial_time_step" ),
           py::arg( "coefficient_set" ),
           py::arg( "minimum_step_size" ),
           py::arg( "maximum_step_size" ),
           py::arg( "relative_error_tolerance" ),
           py::arg( "absolute_error_tolerance" ),
           py::arg( "assess_termination_on_minor_steps" ) = false,
           py::arg( "safety_factor" ) = 0.8,
           py::arg( "maximum_factor_increase" ) = 4.0,
           py::arg( "minimum_factor_increase" ) = 0.1,
           py::arg( "throw_exception_if_minimum_step_exceeded" ) = true,
           R"doc(

 Creates the settings for the Runge-Kutta variable step size integrator with vector tolerances.

 NOTE: THIS FUNCTION IS DEPRECATED, IT IS RECOMMENDED TO USE THE NEW :func:`~tudatpy.numerical_simulation.propagation_setup.integrator.runge_kutta_variable_step` INTERFACE INSTEAD

 Function to create settings for the Runge-Kutta variable step size integrator with vector tolerances.
 For this integrator, the step size is varied based on the tolerances and safety factor provided.
 The tolerance is composed of an absolute and a relative part.
 Different coefficient sets (Butcher's tableau) can be used (see the `CoefficientSets` enum).


 Parameters
 ----------
 initial_time_step : float
     Initial time step to be used.
 coefficient_set : CoefficientSets
     Coefficient set (Butcher's tableau) to be used in the integration.
 minimum_step_size : float
     Minimum time step to be used during the integration.
 maximum_step_size : float
     Maximum time step to be used during the integration.
 relative_error_tolerance : numpy.ndarray[numpy.float64[m, n]]
     Relative vector tolerance to adjust the time step.
 absolute_error_tolerance : numpy.ndarray[numpy.float64[m, n]]
     Absolute vector tolerance to adjust the time step.
 assess_termination_on_minor_steps : bool, default=false
     Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the
     integrator (true) or only at the end of each integration step (false).

 safety_factor : float, default=0.8
     Safety factor used in the step size control.
 maximum_factor_increase : float, default=4.0
     Maximum increase between consecutive time steps, expressed as the factor between new and old step size.

 minimum_factor_increase : float, default=0.1
     Minimum increase between consecutive time steps, expressed as the factor between new and old step size.

 throw_exception_if_minimum_step_exceeded : bool, default=true
     If set to false, the variable step integrator will use the minimum step size specified when the algorithm
     computes the optimum one to be lower, instead of throwing an exception.

 Returns
 -------
 RungeKuttaVariableStepSizeSettingsVectorTolerances
     RungeKuttaVariableStepSizeSettingsVectorTolerances object.






     )doc" );

    m.def( "runge_kutta_variable_step_size",
           &tni::rungeKuttaVariableStepSettingsScalarTolerances< TIME_TYPE >,
           py::arg( "initial_time_step" ),
           py::arg( "coefficient_set" ),
           py::arg( "minimum_step_size" ),
           py::arg( "maximum_step_size" ),
           py::arg( "relative_error_tolerance" ),
           py::arg( "absolute_error_tolerance" ),
           py::arg( "assess_termination_on_minor_steps" ) = false,
           py::arg( "safety_factor" ) = 0.8,
           py::arg( "maximum_factor_increase" ) = 4.0,
           py::arg( "minimum_factor_increase" ) = 0.1,
           py::arg( "throw_exception_if_minimum_step_exceeded" ) = true,
           R"doc(

 Creates the settings for the Runge-Kutta variable step size integrator with scalar tolerances.

 NOTE: THIS FUNCTION IS DEPRECATED, IT IS RECOMMENDED TO USE THE NEW :func:`~tudatpy.numerical_simulation.propagation_setup.integrator.runge_kutta_variable_step` INTERFACE INSTEAD

 Function to create settings for the Runge-Kutta variable step size integrator with scalar tolerances.
 For this integrator, the step size is varied based on the tolerances and safety factor provided.
 The tolerance is composed of an absolute and a relative part.
 Different coefficient sets (Butcher's tableau) can be used (see the `CoefficientSets` enum).


 Parameters
 ----------
 initial_time_step : float
     Initial time step to be used.
 coefficient_set : CoefficientSets
     Coefficient set (Butcher's tableau) to be used in the integration.
 minimum_step_size : float
     Minimum time step to be used during the integration.
 maximum_step_size : float
     Maximum time step to be used during the integration.
 relative_error_tolerance : numpy.ndarray[numpy.float64[m, n]]
     Relative vector tolerance to adjust the time step.
 absolute_error_tolerance : numpy.ndarray[numpy.float64[m, n]]
     Absolute vector tolerance to adjust the time step.
 assess_termination_on_minor_steps : bool, default=false
     Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the
     integrator (true) or only at the end of each integration step (false).

 safety_factor : float, default=0.8
     Safety factor used in the step size control.
 maximum_factor_increase : float, default=4.0
     Maximum increase between consecutive time steps, expressed as the factor between new and old step size.

 minimum_factor_increase : float, default=0.1
     Minimum increase between consecutive time steps, expressed as the factor between new and old step size.

 throw_exception_if_minimum_step_exceeded : bool, default=true
     If set to false, the variable step integrator will use the minimum step size specified when the algorithm
     computes the optimum one to be lower, instead of throwing an exception.

 Returns
 -------
 RungeKuttaVariableStepSettingsScalarTolerances
     RungeKuttaVariableStepSettingsScalarTolerances object.






     )doc" );

    /*!
     * DEPRECATED UNDOCUMENTED
     * -------------------------------------------------------------------------
     *
     */

    m.def( "euler",
           &tni::eulerSettingsDeprecated< TIME_TYPE >,
           py::arg( "initial_time" ),
           py::arg( "initial_time_step" ),
           py::arg( "assess_termination_on_minor_steps" ) = false );

    m.def( "euler",
           &tni::eulerSettings< TIME_TYPE >,
           py::arg( "initial_time_step" ),
           py::arg( "assess_termination_on_minor_steps" ) = false );

    m.def( "runge_kutta_4",
           &tni::rungeKutta4SettingsDeprecated< TIME_TYPE >,
           py::arg( "initial_time" ),
           py::arg( "initial_time_step" ),
           py::arg( "assess_termination_on_minor_steps" ) = false );

    m.def( "runge_kutta_4",
           &tni::rungeKutta4Settings< TIME_TYPE >,
           py::arg( "initial_time_step" ),
           py::arg( "assess_termination_on_minor_steps" ) = false );

    m.def( "runge_kutta_fixed_step_size",
           &tni::rungeKuttaFixedStepSettingsDeprecated< TIME_TYPE >,
           py::arg( "initial_time" ),
           py::arg( "initial_time_step" ),
           py::arg( "coefficient_set" ),
           py::arg( "order_to_use" ) = tni::RungeKuttaCoefficients::OrderEstimateToIntegrate::lower,
           py::arg( "assess_termination_on_minor_steps" ) = false );

    m.def( "runge_kutta_fixed_step_size",
           &tni::rungeKuttaFixedStepSettings< TIME_TYPE >,
           py::arg( "initial_time_step" ),
           py::arg( "coefficient_set" ),
           py::arg( "order_to_use" ) = tni::RungeKuttaCoefficients::OrderEstimateToIntegrate::lower,
           py::arg( "assess_termination_on_minor_steps" ) = false );

    m.def( "runge_kutta_variable_step_size",
           &tni::rungeKuttaVariableStepSettingsScalarTolerancesDeprecated< TIME_TYPE >,
           py::arg( "initial_time" ),
           py::arg( "initial_time_step" ),
           py::arg( "coefficient_set" ),
           py::arg( "minimum_step_size" ),
           py::arg( "maximum_step_size" ),
           py::arg( "relative_error_tolerance" ),
           py::arg( "absolute_error_tolerance" ),
           py::arg( "assess_termination_on_minor_steps" ) = false,
           py::arg( "safety_factor" ) = 0.8,
           py::arg( "maximum_factor_increase" ) = 4.0,
           py::arg( "minimum_factor_increase" ) = 0.1,
           py::arg( "throw_exception_if_minimum_step_exceeded" ) = true );

    m.def( "runge_kutta_variable_step_size_vector_tolerances",
           &tni::rungeKuttaVariableStepSettingsVectorTolerancesDeprecated< TIME_TYPE >,
           py::arg( "initial_time" ),
           py::arg( "initial_time_step" ),
           py::arg( "coefficient_set" ),
           py::arg( "minimum_step_size" ),
           py::arg( "maximum_step_size" ),
           py::arg( "relative_error_tolerance" ),
           py::arg( "absolute_error_tolerance" ),
           py::arg( "assess_termination_on_minor_steps" ) = false,
           py::arg( "safety_factor" ) = 0.8,
           py::arg( "maximum_factor_increase" ) = 4.0,
           py::arg( "minimum_factor_increase" ) = 0.1,
           py::arg( "throw_exception_if_minimum_step_exceeded" ) = true );

    m.def( "bulirsch_stoer",
           &tni::bulirschStoerIntegratorSettingsDeprecated< TIME_TYPE >,
           py::arg( "initial_time" ),
           py::arg( "initial_time_step" ),
           py::arg( "extrapolation_sequence" ),
           py::arg( "maximum_number_of_steps" ),
           py::arg( "minimum_step_size" ),
           py::arg( "maximum_step_size" ),
           py::arg( "relative_error_tolerance" ) = 1.0E-12,
           py::arg( "absolute_error_tolerance" ) = 1.0E-12,
           py::arg( "assess_termination_on_minor_steps" ) = false,
           py::arg( "safety_factor" ) = 0.7,
           py::arg( "maximum_factor_increase" ) = 10.0,
           py::arg( "minimum_factor_increase" ) = 0.1 );

    m.def( "bulirsch_stoer",
           &tni::bulirschStoerIntegratorSettingsDeprecatedNew< TIME_TYPE >,
           py::arg( "initial_time_step" ),
           py::arg( "extrapolation_sequence" ),
           py::arg( "maximum_number_of_steps" ),
           py::arg( "minimum_step_size" ),
           py::arg( "maximum_step_size" ),
           py::arg( "relative_error_tolerance" ) = 1.0E-12,
           py::arg( "absolute_error_tolerance" ) = 1.0E-12,
           py::arg( "assess_termination_on_minor_steps" ) = false,
           py::arg( "safety_factor" ) = 0.7,
           py::arg( "maximum_factor_increase" ) = 10.0,
           py::arg( "minimum_factor_increase" ) = 0.1,
           R"doc(

 Creates the settings for the Bulirsch-Stoer integrator.


 NOTE: THIS FUNCTION IS DEPRECATED, IT IS RECOMMENDED TO USE THE NEW :func:`~tudatpy.numerical_simulation.propagation_setup.integrator.bulirsch_stoer_variable_step` INTERFACE INSTEAD

 Function to create settings for the Bulirsch-Stoer integrator.
 For this integrator, the step size is varied based on the tolerances and safety factor provided.
 The tolerance is composed of an absolute and a relative part.
 Different extrapolation sequences can be used (see the `ExtrapolationMethodStepSequences` enum).


 Parameters
 ----------
 initial_time_step : float
     Initial time step to be used.
 extrapolation_sequence : ExtrapolationMethodStepSequences
     Extrapolation sequence to be used in the integration.
 maximum_number_of_steps : int
     Number of entries in the sequence (e.g., number of integrations used for a single extrapolation).
 minimum_step_size : float
     Minimum time step to be used during the integration.
 maximum_step_size : float
     Maximum time step to be used during the integration.
 relative_error_tolerance : float, default=1.0E-12
     Relative tolerance to adjust the time step.
 absolute_error_tolerance : float, default=1.0E-12
     Relative tolerance to adjust the time step.
 assess_termination_on_minor_steps : bool, default=false
     Whether the propagation termination conditions should be evaluated during the intermediate sub-steps of the integrator (true) or only at the end of each integration step (false).
 safety_factor : float, default=0.7
     Safety factor used in the step size control.
 maximum_factor_increase : float, default=10.0
     Maximum increase between consecutive time steps, expressed as the factor between new and old step size.
 minimum_factor_increase : float, default=0.1
     Minimum increase between consecutive time steps, expressed as the factor between new and old step size.
 Returns
 -------
 BulirschStoerIntegratorSettings
     BulirschStoerIntegratorSettings object.






     )doc" );

    m.def( "adams_bashforth_moulton",
           &tni::adamsBashforthMoultonSettingsDeprecated< TIME_TYPE >,
           py::arg( "initial_time" ),
           py::arg( "initial_time_step" ),
           py::arg( "minimum_step_size" ),
           py::arg( "maximum_step_size" ),
           py::arg( "relative_error_tolerance" ) = 1.0E-12,
           py::arg( "absolute_error_tolerance" ) = 1.0E-12,
           py::arg( "minimum_order" ) = 6,
           py::arg( "maximum_order" ) = 11,
           py::arg( "assess_termination_on_minor_steps" ) = false,
           py::arg( "bandwidth" ) = 200.0,
           "" );
}

}  // namespace integrator
}  // namespace propagation_setup
}  // namespace numerical_simulation
}  // namespace tudatpy
