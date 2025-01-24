/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_aerodynamic_coefficient_setup.h"

#include <tudat/astro/reference_frames/referenceFrameTransformations.h>
#include <tudat/basics/deprecationWarnings.h>
#include <tudat/simulation/environment_setup.h>

//#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
//#include <pybind11/numpy.h>
#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace tss = tudat::simulation_setup;
namespace ti = tudat::interpolators;
namespace ta = tudat::aerodynamics;

namespace tudat {
    namespace simulation_setup {
        //! @get_docstring(customAerodynamicCoefficientSettings)
        inline std::shared_ptr<AerodynamicCoefficientSettings>
        customAerodynamicCoefficientSettingsDeprecatedPy(
            const std::function<Eigen::Vector3d(const std::vector<double> &)>
                forceCoefficientFunction,
            const double referenceArea,
            const std::vector<
                aerodynamics::AerodynamicCoefficientsIndependentVariables>
                independentVariableNames,
            const bool areCoefficientsInAerodynamicFrame = true,
            const bool areCoefficientsInNegativeAxisDirection = true) {
            static bool isWarningPrinted = false;
            if(isWarningPrinted == false) {
                tudat::utilities::printDeprecationWarning(
                    "tudatpy.numerical_simulation.environment_setup."
                    "aerodynamic_coefficients.custom",
                    "tudatpy.numerical_simulation.environment_setup."
                    "aerodynamic_coefficients.custom_aerodynamic_force_"
                    "coefficients");
                isWarningPrinted = true;
            }

            return customAerodynamicCoefficientSettingsDeprecated(
                forceCoefficientFunction, referenceArea,
                independentVariableNames, areCoefficientsInAerodynamicFrame,
                areCoefficientsInNegativeAxisDirection);
        }


    }  // namespace simulation_setup
}  // namespace tudat

namespace tudatpy {
    namespace numerical_simulation {
        namespace environment_setup {
            namespace aerodynamic_coefficients {

                void expose_aerodynamic_coefficient_setup(py::module &m) {
                    /////////////////////////////////////////////////////////////////////////////
                    // createAerodynamicCoefficientInterface.h
                    /////////////////////////////////////////////////////////////////////////////
                    py::class_<
                        tss::AerodynamicCoefficientSettings,
                        std::shared_ptr<tss::AerodynamicCoefficientSettings>>(
                        m, "AerodynamicCoefficientSettings",
                        R"doc(

        Base class for providing settings for aerodynamic interface model.

        Functional (base) class for settings of aerodynamic interface models that require no
        information in addition to their type.
        Aerodynamic interface model settings requiring additional information must be defined using an object derived from this class.





     )doc")
                        .def_property("add_force_contribution_to_moments",
                                      &tss::AerodynamicCoefficientSettings::
                                          getAddForceContributionToMoments,
                                      &tss::AerodynamicCoefficientSettings::
                                          setAddForceContributionToMoments,
                                      R"doc(

        Variable that toggles whether to add the force contribution to the moment coefficients as:

        .. math::
           \Delta \mathbf{C}_{M} = (\mathbf{r}_{ref}-\mathbf{r}_{com})\times \Delta \mathbf{C}_{F}

        where :math:`(\mathbf{r}_{ref}-\mathbf{r}_{com})` is the vector from the center of mass to the moment reference point, and :math:`\mathbf{C}_{F}` and :math:`\mathbf{C}_{M}` is the vector of forc and moment coefficients. Note that, if the force and moment coefficients afre defined in different frames, the relevant frame conversions are automatically performed.
        By default, his boolean is set to false, implicitly assuming that the moment coefficients are provided w.r.t. the (constant) center of mass.


        :type: bool
     )doc")
                        .def_property("moment_reference_point",
                                      &tss::AerodynamicCoefficientSettings::
                                          getMomentReferencePoint,
                                      &tss::AerodynamicCoefficientSettings::
                                          setMomentReferencePoint,
                                      R"doc(No documentation found.)doc")
                        .def("add_single_control_surface",
                             &tss::AerodynamicCoefficientSettings::
                                 addControlSurfaceSettings,
                             py::arg("control_surface_settings"),
                             py::arg("control_surface_name"),
                             R"doc(

        Function to add settings for a single control surface to the coefficient settings. Note that, in Tudat, the
        control surface aerodynamic database inherits the reference properties (length, area, moment reference point)
        from the ``AerodynamicCoefficientSettings`` to which it is assigned.



        Parameters
        ----------
        control_surface_settings : ControlSurfaceIncrementAerodynamicCoefficientSettings
            Settings for aerodynamic coefficients of a control surface

        control_surface_name : str
            Name by which the control surface will be identified





    )doc");

                    py::class_<tss::ConstantAerodynamicCoefficientSettings,
                               std::shared_ptr<
                                   tss::ConstantAerodynamicCoefficientSettings>,
                               tss::AerodynamicCoefficientSettings>(
                        m, "ConstantAerodynamicCoefficientSettings",
                        R"doc(

        Class for defining model settings from constant aerodynamic coefficients.

        `AerodynamicCoefficientSettings` derived class for aerodynamic interface model settings using only constant aerodynamic coefficients.




     )doc");

                    py::class_<tss::CustomAerodynamicCoefficientSettings,
                               std::shared_ptr<
                                   tss::CustomAerodynamicCoefficientSettings>,
                               tss::AerodynamicCoefficientSettings>(
                        m, "CustomAerodynamicCoefficientSettings",
                        R"doc(No documentation found.)doc");
                    //
                    //        py::class_<tss::TabulatedAerodynamicCoefficientSettings,
                    //                std::shared_ptr<tss::TabulatedAerodynamicCoefficientSettings>,
                    //                tss::AerodynamicCoefficientSettings>(
                    //                        m,
                    //                        "TabulatedAerodynamicCoefficientSettings",
                    //                        get_docstring("TabulatedAerodynamicCoefficientSettings").c_str());

                    py::class_<
                        tss::ScaledAerodynamicCoefficientInterfaceSettings,
                        std::shared_ptr<
                            tss::ScaledAerodynamicCoefficientInterfaceSettings>,
                        tss::AerodynamicCoefficientSettings>(
                        m, "ScaledAerodynamicCoefficientInterfaceSettings",
                        R"doc(No documentation found.)doc");

                    py::class_<
                        tss::
                            ControlSurfaceIncrementAerodynamicCoefficientSettings,
                        std::shared_ptr<
                            tss::
                                ControlSurfaceIncrementAerodynamicCoefficientSettings>>(
                        m,
                        "ControlSurfaceIncrementAerodynamicCoefficientSettings",
                        R"doc(No documentation found.)doc");

                    py::class_<
                        tss::
                            CustomControlSurfaceIncrementAerodynamicCoefficientSettings,
                        std::shared_ptr<
                            tss::
                                CustomControlSurfaceIncrementAerodynamicCoefficientSettings>,
                        tss::
                            ControlSurfaceIncrementAerodynamicCoefficientSettings>(
                        m,
                        "CustomControlSurfaceIncrementAerodynamicCoefficientSet"
                        "tings",
                        R"doc(No documentation found.)doc");

                    m.def("constant",
                          py::overload_cast<
                              const double, const Eigen::Vector3d &,
                              const ta::AerodynamicCoefficientFrames>(
                              &tss::constantAerodynamicCoefficientSettings),
                          py::arg("reference_area"),
                          py::arg("constant_force_coefficient"),
                          py::arg("force_coefficients_frame") =
                              ta::negative_aerodynamic_frame_coefficients,
                          R"doc(

Function for creating aerodynamic interface model settings entirely from constant coefficients.

Function for settings object, defining aerodynamic interface model entirely from constant aerodynamic coefficients,
i.e. coefficients are not a function of any independent variables.


Parameters
----------
reference_area : float
    Reference area with which aerodynamic forces and moments are non-dimensionalized.
constant_force_coefficient : numpy.ndarray
    Constant force coefficients.
force_coefficients_frame : AerodynamicCoefficientFrames, default = negative_aerodynamic_frame_coefficients
    Variable defining the frame in which the force coefficients are defined. By default, this is the negative aerodynamic
    frame, so that the coefficients are for drag, side force and lift

Returns
-------
ConstantAerodynamicCoefficientSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.ConstantAerodynamicCoefficientSettings` class





Examples
--------
In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` for the artificial body "Vehicle", using only constant aerodynamic coefficients:

.. code-block:: python

  # Define the reference area and constant aerodynamic coefficients
  reference_area = 20.0
  drag_coefficient = 1.5
  lift_coefficient = 0.3
  # Create the aerodynamic interface settings
  aero_coefficient_settings = environment_setup.aerodynamic_coefficients.constant(
      reference_area,
      constant_force_coefficient=[drag_coefficient, 0, lift_coefficient],
      force_coefficients_frame=environment.negative_aerodynamic_frame_coefficients,
  )
  # Assign aerodynamic interface to the vehicle
  environment_setup.add_aerodynamic_coefficient_interface(bodies, "Vehicle", aero_coefficient_settings)


    )doc");


                    m.def(
                        "constant_force_and_moment",
                        &tss::
                            constantAerodynamicForceAndMomentCoefficientSettings,
                        py::arg("reference_length"), py::arg("reference_area"),
                        py::arg("moment_reference_point"),
                        py::arg("constant_force_coefficient"),
                        py::arg("constant_moment_coefficient"),
                        py::arg("force_coefficients_frame") =
                            ta::negative_aerodynamic_frame_coefficients,
                        py::arg("moment_coefficients_frame") =
                            ta::body_fixed_frame_coefficients,
                        R"doc(

Function for creating aerodynamic interface model settings entirely from constant coefficients.

Function for settings object, defining aerodynamic interface model entirely from constant aerodynamic coefficients,
i.e. coefficients are not a function of any independent variables.


Parameters
----------
reference_area : float
    Reference area with which aerodynamic forces and moments are non-dimensionalized.
constant_force_coefficient : numpy.ndarray
    Constant force coefficients.
force_coefficients_frame : AerodynamicCoefficientFrames, default = negative_aerodynamic_frame_coefficients
    Variable defining the frame in which the force coefficients are defined. By default, this is the negative aerodynamic
    frame, so that the coefficients are for drag, side force and lift

Returns
-------
ConstantAerodynamicCoefficientSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.ConstantAerodynamicCoefficientSettings` class





Examples
--------
In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` for the artificial body "Vehicle", using only constant aerodynamic coefficients:

.. code-block:: python

  # Define the reference area and constant aerodynamic coefficients
  reference_area = 20.0
  drag_coefficient = 1.5
  lift_coefficient = 0.3
  # Create the aerodynamic interface settings
  aero_coefficient_settings = environment_setup.aerodynamic_coefficients.constant(
      reference_area,
      constant_force_coefficient=[drag_coefficient, 0, lift_coefficient],
      force_coefficients_frame=environment.negative_aerodynamic_frame_coefficients,
  )
  # Assign aerodynamic interface to the vehicle
  environment_setup.add_aerodynamic_coefficient_interface(bodies, "Vehicle", aero_coefficient_settings)


    )doc");

                    m.def("constant",
                          py::overload_cast<
                              const double, const Eigen::Vector3d &,
                              const ta::AerodynamicCoefficientFrames>(
                              &tss::constantAerodynamicCoefficientSettings),
                          py::arg("reference_area"),
                          py::arg("constant_force_coefficient"),
                          py::arg("force_coefficients_frame") =
                              ta::negative_aerodynamic_frame_coefficients,
                          R"doc(

Function for creating aerodynamic interface model settings entirely from constant coefficients.

Function for settings object, defining aerodynamic interface model entirely from constant aerodynamic coefficients,
i.e. coefficients are not a function of any independent variables.


Parameters
----------
reference_area : float
    Reference area with which aerodynamic forces and moments are non-dimensionalized.
constant_force_coefficient : numpy.ndarray
    Constant force coefficients.
force_coefficients_frame : AerodynamicCoefficientFrames, default = negative_aerodynamic_frame_coefficients
    Variable defining the frame in which the force coefficients are defined. By default, this is the negative aerodynamic
    frame, so that the coefficients are for drag, side force and lift

Returns
-------
ConstantAerodynamicCoefficientSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.ConstantAerodynamicCoefficientSettings` class





Examples
--------
In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` for the artificial body "Vehicle", using only constant aerodynamic coefficients:

.. code-block:: python

  # Define the reference area and constant aerodynamic coefficients
  reference_area = 20.0
  drag_coefficient = 1.5
  lift_coefficient = 0.3
  # Create the aerodynamic interface settings
  aero_coefficient_settings = environment_setup.aerodynamic_coefficients.constant(
      reference_area,
      constant_force_coefficient=[drag_coefficient, 0, lift_coefficient],
      force_coefficients_frame=environment.negative_aerodynamic_frame_coefficients,
  )
  # Assign aerodynamic interface to the vehicle
  environment_setup.add_aerodynamic_coefficient_interface(bodies, "Vehicle", aero_coefficient_settings)


    )doc");


                    m.def(
                        "custom_aerodynamic_force_coefficients",
                        py::overload_cast<
                            const std::function<Eigen::Vector3d(
                                const std::vector<double> &)>,
                            const double,
                            const std::vector<
                                ta::AerodynamicCoefficientsIndependentVariables>,
                            const ta::AerodynamicCoefficientFrames>(
                            &tss::customAerodynamicCoefficientSettings),
                        py::arg("force_coefficient_function"),
                        py::arg("reference_area"),
                        py::arg("independent_variable_names"),
                        py::arg("force_coefficients_frame") =
                            ta::negative_aerodynamic_frame_coefficients,
                        R"doc(

Function for creating aerodynamic interface model settings from custom coefficients.

Function for settings object, defining aerodynamic interface model via a custom force coefficient function
(function of independent variable).


Parameters
----------
force_coefficient_function : callable[[list[float]], numpy.ndarray[numpy.float64[3, 1]]]
    Function that is defining the aerodynamic coefficients as function of an independent variable (see arg independent_variable_names).
reference_area : float
    Reference area with which aerodynamic forces and moments are non-dimensionalized.
independent_variable_name : list[environment.AerodynamicCoefficientsIndependentVariables]
    Vector with identifiers for the independent variable w.r.t. which the aerodynamic coefficients are defined.
force_coefficients_frame : AerodynamicCoefficientFrames, default = negative_aerodynamic_frame_coefficients
    Variable defining the frame in which the force coefficients are defined. By default, this is the negative aerodynamic
    frame, so that the coefficients are for drag, side force and lift (:math:`C_{D}, C_{S}, C_{L}`)

Returns
-------
CustomAerodynamicCoefficientSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.CustomAerodynamicCoefficientSettings` class





Examples
--------
In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` for the artificial body "Vehicle", using a function based on the mach number:

.. code-block:: python

  def force_coefficients(variables_list):
    # Extract the mach number
    mach_number = variables_list[0]
    # If the mach number is below 3, use fixed coefficients
    if mach_number <= 3:
        return [0.99, 0, 1.08]
    # Same if the mach number is above 10
    elif mach_number >= 10:
        return [0.82, 0, 0.88]
    # Otherwise, vary linearly between the ones at M=3 and M=10
    CD = 1.0667-0.02457*mach_number
    CL = 1.1636-0.02786*mach_number
    return [CD, 0, CL]
  # Create the aerodynamic interface settings
  aero_coefficient_settings = environment_setup.aerodynamic_coefficients.custom(
      force_coefficients,
      reference_area=1.50,
      independent_variable_names=[environment.AerodynamicCoefficientsIndependentVariables.mach_number_dependent]
  )
  # Assign the aerodynamic coefficient interface to the vehicle
  environment_setup.add_aerodynamic_coefficient_interface(bodies, "Vehicle", aero_coefficient_settings)


    )doc");

                    m.def(
                        "custom_aerodynamic_force_and_moment_coefficients",
                        py::overload_cast<
                            const std::function<Eigen::Vector3d(
                                const std::vector<double> &)>,
                            const std::function<Eigen::Vector3d(
                                const std::vector<double> &)>,
                            const double, const double,
                            const std::vector<
                                ta::AerodynamicCoefficientsIndependentVariables>,
                            const ta::AerodynamicCoefficientFrames,
                            const ta::AerodynamicCoefficientFrames,
                            const Eigen::Vector3d &>(
                            &tss::customAerodynamicCoefficientSettings),
                        py::arg("force_coefficient_function"),
                        py::arg("moment_coefficient_function"),
                        py::arg("reference_length"), py::arg("reference_area"),
                        py::arg("independent_variable_names"),
                        py::arg("force_coefficients_frame") =
                            ta::negative_aerodynamic_frame_coefficients,
                        py::arg("moment_coefficients_frame") =
                            ta::body_fixed_frame_coefficients,
                        py::arg("moment_reference_point") =
                            Eigen::Vector3d::Constant(TUDAT_NAN),
                        R"doc(

Function for creating aerodynamic interface model settings from custom coefficients.

Function for settings object, defining aerodynamic interface model via a custom force and moment coefficient function
(function of independent variable).


Parameters
----------
force_coefficient_function : callable[[list[float]], numpy.ndarray[numpy.float64[3, 1]]]
    Function that is defining the aerodynamic force coefficients as function of an independent variable (see arg independent_variable_names).
moment_coefficient_function : callable[[list[float]], numpy.ndarray[numpy.float64[3, 1]]]
    Function that is defining the aerodynamic moment coefficients as function of an independent variable (see arg independent_variable_names).
reference_area : float
    Reference area with which aerodynamic forces and moments are non-dimensionalized.
reference_length : float
    Reference length with which aerodynamic moments are non-dimensionalized.
independent_variable_name : list[environment.AerodynamicCoefficientsIndependentVariables]
    Vector with identifiers for the independent variable w.r.t. which the aerodynamic coefficients are defined.
force_coefficients_frame : AerodynamicCoefficientFrames, default = negative_aerodynamic_frame_coefficients
    Variable defining the frame in which the force coefficients are defined. By default, this is the negative aerodynamic
    frame, so that the coefficients are for drag, side force and lift (:math:`C_{D}, C_{S}, C_{L}`)

moment_coefficients_frame : AerodynamicCoefficientFrames, default = positive_body_frame_coefficients
    Variable defining the frame in which the moment coefficients are defined. By default, this is the positive body
    frame, so that the coefficients are roll, pitch and yaw (:math:`C_{l}, C_{m}, C_{n}`)

moment_reference_point : numpy.ndarray[numpy.float64[3, 1]] = np.full([3, 1], np.nan)
    Point w.r.t. aerodynamic moment coefficients are defined. This variable is used to calculate the contribution of the aerodynamic
    force coefficients to the effective moment coefficients. See the ``add_force_contribution_to_moments`` attribute of the
    :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` for more details.
    If the present input is set to NaN (as is the default), the reference point is left undefined, and the aerodynamic moments are computed
    without computing any force coefficient contribution to the moment coefficients.

Returns
-------
CustomAerodynamicCoefficientSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.CustomAerodynamicCoefficientSettings` class






    )doc");


                    m.def(
                        "tabulated",
                        py::overload_cast<
                            const std::vector<double>,
                            const std::vector<Eigen::Vector3d>,
                            const std::vector<Eigen::Vector3d>, const double,
                            const double,
                            const ta::
                                AerodynamicCoefficientsIndependentVariables,
                            const ta::AerodynamicCoefficientFrames,
                            const ta::AerodynamicCoefficientFrames,
                            const Eigen::Vector3d &,
                            const std::shared_ptr<ti::InterpolatorSettings>>(
                            &tss::
                                oneDimensionalTabulatedAerodynamicCoefficientSettings),
                        py::arg("independent_variables"),
                        py::arg("force_coefficients"),
                        py::arg("moment_coefficients"),
                        py::arg("reference_length"), py::arg("reference_area"),
                        py::arg("independent_variable_name"),
                        py::arg("force_coefficients_frame") =
                            ta::negative_aerodynamic_frame_coefficients,
                        py::arg("moment_coefficients_frame") =
                            ta::body_fixed_frame_coefficients,
                        py::arg("moment_reference_point") =
                            Eigen::Vector3d::Constant(TUDAT_NAN),
                        py::arg("interpolator_settings") = nullptr,
                        R"doc(

Function for creating aerodynamic interface model settings from user-defined, 1-d tabulated coefficients.

Function for settings object, defining aerodynamic interface model via user-defined, 1-dimensional, tabulated aerodynamic force and moment coefficients
(tabulated w.r.t. independent variable).


Parameters
----------
independent_variables : list[float]
    Values of independent variables at which the coefficients in the input multi vector are defined (size 1).
force_coefficients : list[numpy.ndarray[numpy.float64[3, 1]]]
    Values of force coefficients at independent variables defined by independent_variables.
moment_coefficients : list[numpy.ndarray[numpy.float64[3, 1]]]
    Values of moment coefficients at independent variables defined by independent_variables.
reference_length : float
    Reference length with which aerodynamic moments about x- and z- axes are non-dimensionalized.
reference_area : float
    Reference area with which aerodynamic forces and moments are non-dimensionalized.
independent_variable_name : list[environment.AerodynamicCoefficientsIndependentVariables]
    Vector with identifiers for the independent variable w.r.t. which the aerodynamic coefficients are defined.
force_coefficients_frame : AerodynamicCoefficientFrames, default = negative_aerodynamic_frame_coefficients
    Variable defining the frame in which the force coefficients are defined. By default, this is the negative aerodynamic
    frame, so that the coefficients are for drag, side force and lift (:math:`C_{D}, C_{S}, C_{L}`)

moment_coefficients_frame : AerodynamicCoefficientFrames, default = positive_body_frame_coefficients
    Variable defining the frame in which the moment coefficients are defined. By default, this is the positive body
    frame, so that the coefficients are roll, pitch yaw (:math:`C_{l}, C_{m}, C_{n}`)

moment_reference_point : numpy.ndarray[numpy.float64[3, 1]] = np.full([3, 1], np.nan)
    Point w.r.t. aerodynamic moment coefficients are defined. This variable is used to calculate the contribution of the aerodynamic
    force coefficients to the effective moment coefficients. See the ``add_force_contribution_to_moments`` attribute of the
    :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` for more details.
    If the present input is set to NaN (as is the default), the reference point is left undefined, and the aerodynamic moments are computed
    without computing any force coefficient contribution to the moment coefficients.

interpolator_settings : math.interpolators.InterpolatorSettings, default = None
    Interpolator settings object, where the conditions for interpolation of tabulated inputs are saved.
Returns
-------
TabulatedAerodynamicCoefficientSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.TabulatedAerodynamicCoefficientSettings` class (via :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.TabulatedAerodynamicCoefficientSettingsBase` class)





Examples
--------
In this example, aerodynamic force and moment coefficients are defined as multi-dimensional arrays.
The values for the aerodynamic coefficients vary with Mach number, and are defined for Mach numbers of 3, 5, 10, and 15.
This example also shows how to set the required reference point, lengths, and area.

.. code-block:: python

  # Define the aerodynamic force coefficients [CD, CS, CL] for different mach numbers
  aero_coefficients_array_force = [
      [0.7647, 0, 0.9722],
      [0.6729, 0, 0.8461],
      [0.6240, 0, 0.7838],
      [0.6246, 0, 0.7841]
  ]
  # Define the aerodynamic moment coefficients for different mach numbers
  aero_coefficients_array_moment = [
      [0.45, 0, 0],
      [0.50, 0, 0],
      [0.53, 0, 0],
      [0.55, 0, 0]
  ]
  # Create the aerodynamic interface settings
  aero_coefficient_settings = environment_setup.aerodynamic_coefficients.tabulated(
      independent_variables=[3, 5, 10, 15],       # Mach number at which the coefficients are defined
      force_coefficients=aero_coefficients_array_force,
      moment_coefficients=aero_coefficients_array_moment,
      reference_length=0.25,
      reference_area=1.50,
      independent_variable_name=environment.AerodynamicCoefficientsIndependentVariables.mach_number_dependent
  )
  # Assign the aerodynamic coefficient interface to the vehicle
  environment_setup.add_aerodynamic_coefficient_interface(bodies, "Vehicle", aero_coefficient_settings)


    )doc");

                    m.def(
                        "tabulated_force_only",
                        py::overload_cast<
                            const std::vector<double>,
                            const std::vector<Eigen::Vector3d>, const double,
                            const ta::
                                AerodynamicCoefficientsIndependentVariables,
                            const ta::AerodynamicCoefficientFrames,
                            const std::shared_ptr<ti::InterpolatorSettings>>(
                            &tss::
                                oneDimensionalTabulatedAerodynamicCoefficientSettings),
                        py::arg("independent_variables"),
                        py::arg("force_coefficients"),
                        py::arg("reference_area"),
                        py::arg("independent_variable_name"),
                        py::arg("force_coefficients_frame") =
                            ta::negative_aerodynamic_frame_coefficients,
                        py::arg("interpolator_settings") = nullptr,
                        R"doc(

Function for creating aerodynamic interface model settings from user-defined, 1-d tabulated force coefficients.

Function for settings object, defining aerodynamic interface model via user-defined, 1-dimensional, tabulated aerodynamic force coefficients
(tabulated w.r.t. independent variable).


Parameters
----------
independent_variables : list[float]
    Values of independent variables at which the coefficients in the input multi vector are defined (size 1)
force_coefficients : list[numpy.ndarray[numpy.float64[3, 1]]]
    Values of force coefficients at independent variables defined by independent_variables.
reference_area : float
    Reference area with which aerodynamic forces and moments are non-dimensionalized.
independent_variable_name : environment.AerodynamicCoefficientsIndependentVariables
    Identifier of the independent variable w.r.t. which the aerodynamic coefficients are defined.
force_coefficients_frame : AerodynamicCoefficientFrames, default = negative_aerodynamic_frame_coefficients
    Variable defining the frame in which the force coefficients are defined. By default, this is the negative aerodynamic
    frame, so that the coefficients are for drag, side force and lift (:math:`C_{D}, C_{S}, C_{L}`)

interpolator_settings : math.interpolators.InterpolatorSettings, default = None
    Interpolator settings object, where the conditions for interpolation of tabulated inputs are saved.
    Pointer to an interpolator settings object where the conditions for interpolation of tabulated inputs are saved.

Returns
-------
TabulatedAerodynamicCoefficientSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.TabulatedAerodynamicCoefficientSettings` class





Examples
--------
In this example, aerodynamic force coefficients are defined as a multi-dimensional array.
The values for the force coefficients vary with Mach number, and are defined for Mach numbers of 3, 5, 10, and 15.

.. code-block:: python

  # Define the aerodynamic coefficients [CD, CS, CL] for different mach numbers
  aero_coefficients_array = [
      [0.7647, 0, 0.9722],
      [0.6729, 0, 0.8461],
      [0.6240, 0, 0.7838],
      [0.6246, 0, 0.7841]
  ]
  # Create the aerodynamic interface settings
  aero_coefficient_settings = environment_setup.aerodynamic_coefficients.tabulated_force_only(
      independent_variables=[3.0, 5.0, 10.0, 15.0],       # Mach number at which the coefficients are defined
      force_coefficients=aero_coefficients_array,
      reference_area=1.50,
      independent_variable_name=environment.AerodynamicCoefficientsIndependentVariables.mach_number_dependent
  )
  # Assign the aerodynamic coefficient interface to the vehicle
  environment_setup.add_aerodynamic_coefficient_interface(bodies, "Vehicle", aero_coefficient_settings)


    )doc");

                    m.def(
                        "tabulated_force_only_from_files",
                        py::overload_cast<
                            const std::map<int, std::string>, const double,
                            const std::vector<
                                ta::AerodynamicCoefficientsIndependentVariables>,
                            const ta::AerodynamicCoefficientFrames,
                            const std::shared_ptr<ti::InterpolatorSettings>>(
                            &tss::
                                readTabulatedAerodynamicCoefficientsFromFiles),
                        py::arg("force_coefficient_files"),
                        py::arg("reference_area"),
                        py::arg("independent_variable_names"),
                        py::arg("force_coefficients_frame") =
                            ta::negative_aerodynamic_frame_coefficients,
                        py::arg("interpolator_settings") = nullptr,
                        R"doc(

Function for creating aerodynamic interface model settings from tabulated force coefficients from files.

Function for settings object, defining aerodynamic interface model via user-defined, tabulated aerodynamic force coefficients
(tabulated w.r.t. independent variable), obtained from data files.


Parameters
----------
force_coefficient_files : Dict[int, str]
    Path of the aerodynamic coefficient files corresponding to the force coefficient of the given dict key.
reference_area : float
    Reference area with which aerodynamic forces and moments are non-dimensionalized.
independent_variable_names : list[environment.AerodynamicCoefficientsIndependentVariables]
    Vector with identifiers for the independent variable w.r.t. which the aerodynamic coefficients are defined.
force_coefficients_frame : AerodynamicCoefficientFrames, default = negative_aerodynamic_frame_coefficients
    Variable defining the frame in which the force coefficients are defined. By default, this is the negative aerodynamic
    frame, so that the coefficients are for drag, side force and lift (:math:`C_{D}, C_{S}, C_{L}`)

interpolator_settings : math.interpolators.InterpolatorSettings, default = None
    Interpolator settings object, where the conditions for interpolation of tabulated inputs are saved.

Returns
-------
TabulatedAerodynamicCoefficientSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.TabulatedAerodynamicCoefficientSettings` class





Examples
--------
In this example, the drag and lift coefficients of the Space Transport System are defined from two data files.
Both of these data files contain coefficient values dependent on both the angle of attack and the mach number,
as shown in the example in the `independent_variable_names` input.
This example is taken from the `reentry trajectory example <https://github.com/tudat-team/tudatpy-examples/blob/1f8180b0064226175bbe66e3eaf044f229a897f6/propagation/reentry_trajectory.py>`_.

.. code-block:: python

  # Define the aerodynamic coefficient files (leave C_S empty)
  aero_coefficients_files = {0: "input/STS_CD.dat", 2:"input/STS_CL.dat"}
  # Setup the aerodynamic coefficients settings tabulated from the files
  coefficient_settings = environment_setup.aerodynamic_coefficients.tabulated_force_only_from_files(
      force_coefficient_files=aero_coefficients_files,
      reference_area=2690.0*0.3048*0.3048,
      independent_variable_names=[environment.angle_of_attack_dependent, environment.mach_number_dependent]
  )
  # Add predefined aerodynamic coefficients database to the body
  environment_setup.add_aerodynamic_coefficient_interface(bodies, "STS", coefficient_settings)


    )doc");

                    m.def(
                        "tabulated_from_files",
                        py::overload_cast<
                            const std::map<int, std::string>,
                            const std::map<int, std::string>, const double,
                            const double,
                            const std::vector<
                                ta::AerodynamicCoefficientsIndependentVariables>,
                            const ta::AerodynamicCoefficientFrames,
                            const ta::AerodynamicCoefficientFrames,
                            const Eigen::Vector3d &,
                            const std::shared_ptr<ti::InterpolatorSettings>>(
                            &tss::
                                readTabulatedAerodynamicCoefficientsFromFiles),
                        py::arg("force_coefficient_files"),
                        py::arg("moment_coefficient_files"),
                        py::arg("reference_length"), py::arg("reference_area"),
                        py::arg("independent_variable_names"),
                        py::arg("force_coefficients_frame") =
                            ta::negative_aerodynamic_frame_coefficients,
                        py::arg("moment_coefficients_frame") =
                            ta::body_fixed_frame_coefficients,
                        py::arg("moment_reference_point") =
                            Eigen::Vector3d::Constant(TUDAT_NAN),
                        py::arg("interpolator_settings") = nullptr,
                        R"doc(

Function for creating aerodynamic interface model settings from tabulated coefficients from files.

Function for settings object, defining aerodynamic interface model via user-defined, tabulated aerodynamic force and moment coefficients
(tabulated w.r.t. independent variable), obtained from data files.


Parameters
----------
force_coefficient_files : Dict[int, str]
    Path of the aerodynamic coefficient files corresponding to the force coefficient of the given dict key (0, 1 and 2 a are x-, y- and z-axis of force frame, respectively).
moment_coefficient_files : Dict[int, str]
    Path of the aerodynamic coefficient files corresponding to the moment coefficient of the given dict key (0, 1 and 2 a are x-, y- and z-axis of moment frame, respectively).
reference_length : float
    Reference length with which aerodynamic moments about x- and z- axes are non-dimensionalized.
reference_area : float
    Reference area with which aerodynamic forces and moments are non-dimensionalized.
independent_variable_names : list[environment.AerodynamicCoefficientsIndependentVariables]
    Vector with identifiers for the independent variable w.r.t. which the aerodynamic coefficients are defined.
force_coefficients_frame : AerodynamicCoefficientFrames, default = negative_aerodynamic_frame_coefficients
    Variable defining the frame in which the force coefficients are defined. By default, this is the negative aerodynamic
    frame, so that the coefficients are for drag, side force and lift (:math:`C_{D}, C_{S}, C_{L}`)

moment_coefficients_frame : AerodynamicCoefficientFrames, default = positive_body_frame_coefficients
    Variable defining the frame in which the moment coefficients are defined. By default, this is the positive body
    frame, so that the coefficients are roll, pitch yaw (:math:`C_{l}, C_{m}, C_{n}`)

moment_reference_point : numpy.ndarray[numpy.float64[3, 1]] = np.full([3, 1], np.nan)
    Point w.r.t. aerodynamic moment coefficients are defined. This variable is used to calculate the contribution of the aerodynamic
    force coefficients to the effective moment coefficients. See the ``add_force_contribution_to_moments`` attribute of the
    :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` for more details.
    If the present input is set to NaN (as is the default), the reference point is left undefined, and the aerodynamic moments are computed
    without computing any force coefficient contribution to the moment coefficients.

interpolator_settings : math.interpolators.InterpolatorSettings, default = None
    Interpolator settings object, where the conditions for interpolation of tabulated inputs are saved.

Returns
-------
TabulatedAerodynamicCoefficientSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.TabulatedAerodynamicCoefficientSettings` class





Examples
--------
This example is very similar to the one for `tabulated_force_only_from_files`, with the distinction that a pitching moment coefficient is added.

.. code-block:: python

  # Define the force coefficient files (leave C_S empty)
  force_coefficients_files = {0: "input/STS_CD.dat", 2:"input/STS_CL.dat"}
  # Define the moment coefficient files (leave C_S empty)
  moment_coefficients_files = {0: "input/STS_CM.dat"}
  # Setup the aerodynamic coefficients settings tabulated from the files
  coefficient_settings = environment_setup.aerodynamic_coefficients.tabulated_from_files(
      force_coefficient_files=force_coefficients_files,
      moment_coefficient_files=moment_coefficients_files,
      reference_length=11.9,
      reference_area=2690.0*0.3048*0.3048,
      independent_variable_names=[environment.angle_of_attack_dependent, environment.mach_number_dependent]
  )
  # Add the predefined aerodynamic coefficients database to the body
  environment_setup.add_aerodynamic_coefficient_interface(bodies, "STS", coefficient_settings)


    )doc");

                    m.def("scaled_by_constant",
                          py::overload_cast<
                              const std::shared_ptr<
                                  tss::AerodynamicCoefficientSettings>,
                              const double, const double, const bool>(
                              &tss::scaledAerodynamicCoefficientSettings),
                          py::arg("unscaled_coefficient_settings"),
                          py::arg("force_scaling_constant"),
                          py::arg("moment_scaling_constant"),
                          py::arg("is_scaling_absolute") = false,
                          R"doc(

Function for creating aerodynamic interface model settings by applying one constant scaling factor/value to all coefficients of an existing model settings object.

Function for settings object, defining aerodynamic interface based on scaling the coefficients of an existing model settings object by one constant factor or value.
Via the ``is_scaling_absolute``
boolean, the user can apply a constant scaling factor or an absolute value to the resulting force and moment coefficients (for instance for an uncertainty analysis).


Parameters
----------
unscaled_coefficient_settings : AerodynamicCoefficientSettings
    Existing aerodynamic interface model settings object that is used as the base for the scaled settings object.
force_scaling_constant : float
    Constant scaling factor to be applied to all aerodynamic force coefficients.
moment_scaling_constant : float
    Constant scaling factor to be applied to all aerodynamic moment coefficients.
is_scaling_absolute : bool, default = False
    Boolean indicating whether aerodynamic coefficient scaling is absolute.
    Setting this boolean to true will add the scaling value to the base value,
    instead of the default behaviour of multiplying the base value by the scaling factor.

Returns
-------
ScaledAerodynamicCoefficientInterfaceSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.ScaledAerodynamicCoefficientInterfaceSettings` class





Examples
--------
In this example, we first set constant aerodynamic coefficients, like in the earlier example.
Then, we use the `scaled_by_constant` function to scale the force coefficients by 1.1.
Since the `is_scaling_absolute` equals `False` by default, the force coefficients are then increased by 10%.

.. code-block:: python

  # Define the reference area and constant aerodynamic coefficients
  reference_area = 20.0
  drag_coefficient = 1.5
  lift_coefficient = 0.3
  # Create the aerodynamic interface settings
  aero_coefficient_settings = environment_setup.aerodynamic_coefficients.constant(
      reference_area,
      constant_force_coefficient=[drag_coefficient, 0, lift_coefficient]
  )
  # Define scaled aerodynamic coefficient to increase the force coefficients by 10%
  scaled_aero_coefficient_settings = environment_setup.aerodynamic_coefficients.scaled_by_constant(
      unscaled_coefficient_settings=aero_coefficient_settings,
      force_scaling_constant=1.1,
      moment_scaling_constant=1.0
  )
  # Assign aerodynamic interface to the vehicle
  environment_setup.add_aerodynamic_coefficient_interface(bodies, "Vehicle", scaled_aero_coefficient_settings)


    )doc");

                    m.def("scaled_by_vector",
                          py::overload_cast<
                              const std::shared_ptr<
                                  tss::AerodynamicCoefficientSettings>,
                              const Eigen::Vector3d, const Eigen::Vector3d,
                              const bool>(
                              &tss::scaledAerodynamicCoefficientSettings),
                          py::arg("unscaled_coefficient_settings"),
                          py::arg("force_scaling_vector"),
                          py::arg("moment_scaling_vector"),
                          py::arg("is_scaling_absolute") = false,
                          R"doc(

Function for creating aerodynamic interface model settings by applying constant scaling factors/values to the coefficients of an existing model settings object.

Function for settings object, defining aerodynamic interface based on scaling the coefficients of an existing model settings object by constant factors or values.
Via the ``is_scaling_absolute`` boolean, the user can apply one constant scaling factor or an absolute value to each resulting force and moment coefficient (for instance for an uncertainty analysis).


Parameters
----------
unscaled_coefficient_settings : AerodynamicCoefficientSettings
    Existing aerodynamic interface model settings object that is used as the base for the scaled settings object.
force_scaling_vector : numpy.ndarray[numpy.float64[3, 1]]
    Constant scaling factors to be applied to each aerodynamic force coefficient.
moment_scaling_vector : numpy.ndarray[numpy.float64[3, 1]]
    Constant scaling factors to be applied to each aerodynamic moment coefficient.
is_scaling_absolute : bool, default = False
    Boolean indicating whether aerodynamic coefficient scaling is absolute.
    Setting this boolean to true will add the scaling value to the base value,
    instead of the default behaviour of multiplying the base value by the scaling factor.

Returns
-------
ScaledAerodynamicCoefficientInterfaceSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.ScaledAerodynamicCoefficientInterfaceSettings` class





Examples
--------
In this example, we first set constant aerodynamic coefficients, like in the earlier example.
Then, we use the `scaled_by_vector` function to scale the drag coefficient by 2.

.. code-block:: python

  # Define the reference area and constant aerodynamic coefficients
  reference_area = 20.0
  drag_coefficient = 1.5
  lift_coefficient = 0.3
  # Create the aerodynamic interface settings
  aero_coefficient_settings = environment_setup.aerodynamic_coefficients.constant(
      reference_area,
      constant_force_coefficient=[drag_coefficient, 0, lift_coefficient]
  )
  # Define scaled aerodynamic coefficient to increase CD by a factor of 2
  scaled_aero_coefficient_settings = environment_setup.aerodynamic_coefficients.scaled_by_vector(
      unscaled_coefficient_settings=aero_coefficient_settings,
      force_scaling_vector=[2.0, 1.0, 1.0],
      moment_scaling_vector=[1.0, 1.0, 1.0]
  )
  # Assign aerodynamic interface to the vehicle
  environment_setup.add_aerodynamic_coefficient_interface(bodies, "Vehicle", scaled_aero_coefficient_settings)


    )doc");

                    m.def(
                        "scaled_by_vector_function",
                        py::overload_cast<
                            const std::shared_ptr<
                                tss::AerodynamicCoefficientSettings>,
                            const std::function<Eigen::Vector3d(const double)>,
                            const std::function<Eigen::Vector3d(const double)>,
                            const bool>(
                            &tss::scaledAerodynamicCoefficientSettings),
                        py::arg("unscaled_coefficient_settings"),
                        py::arg("force_scaling_vector_function"),
                        py::arg("moment_scaling_vector_function"),
                        py::arg("is_scaling_absolute") = false,
                        R"doc(

Function for creating aerodynamic interface model settings by applying custom scaling factors/values to the coefficients of an existing model settings object.

Function for settings object, defining aerodynamic interface based on scaling the coefficients of an existing model settings object by custom factors or values.
Via the ``is_scaling_absolute`` boolean, the user can apply the scaling factors or absolute values to each resulting force and moment coefficient (for instance for an uncertainty analysis).


Parameters
----------
unscaled_coefficient_settings : AerodynamicCoefficientSettings
    Existing aerodynamic interface model settings object that is used as the base for the scaled settings object.
force_scaling_vector_function : callable[[float], numpy.ndarray[numpy.float64[3, 1]]]
    Custom scaling factors to be applied to each aerodynamic force coefficient.
moment_scaling_vector_function : callable[[float], numpy.ndarray[numpy.float64[3, 1]]]
    Custom scaling factors to be applied to each aerodynamic moment coefficient.
is_scaling_absolute : bool, default = False
    Boolean indicating whether aerodynamic coefficient scaling is absolute.
    Setting this boolean to true will add the scaling value to the base value,
    instead of the default behaviour of multiplying the base value by the scaling factor.

Returns
-------
ScaledAerodynamicCoefficientInterfaceSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.ScaledAerodynamicCoefficientInterfaceSettings` class





Examples
--------
In this example, we first set constant aerodynamic coefficients, like in the earlier example.
Then, we use the `scaled_by_vector_function` function to scale the drag and lift coefficients according to a function that varies with time.
This scaling function essentially adds noise to the CD and CL following as a sin or cos function.

.. code-block:: python

  # Define the reference area and constant aerodynamic coefficients
  reference_area = 20.0
  drag_coefficient = 1.5
  lift_coefficient = 0.3
  # Create the aerodynamic interface settings
  aero_coefficient_settings = environment_setup.aerodynamic_coefficients.constant(
      reference_area,
      constant_force_coefficient=[drag_coefficient, 0, lift_coefficient]
  )
  # Define the aerodynamic coefficient scaling as a function of time
  def aero_coefficient_scaling(time):
      CD_scale = 1 + 0.25*np.sin(time/10)
      CL_scale = 1 + 0.25*np.cos(time/15)
      return [CD_scale, 1.0, CL_scale]
  # Define scaled aerodynamic coefficient to increase CD by a factor of 2
  scaled_aero_coefficient_settings = environment_setup.aerodynamic_coefficients.scaled_by_vector_function(
      unscaled_coefficient_settings=aero_coefficient_settings,
      force_scaling_vector_function=aero_coefficient_scaling,
      moment_scaling_vector_function=lambda x: [1.0, 1.0, 1.0]
  )
  # Assign aerodynamic interface to the vehicle
  environment_setup.add_aerodynamic_coefficient_interface(bodies, "Vehicle", scaled_aero_coefficient_settings)


    )doc");

                    m.def(
                        "custom_control_surface",
                        &tss::
                            customControlSurfaceIncrementAerodynamicCoefficientSettings,
                        py::arg("force_and_moment_coefficient_function"),
                        py::arg("independent_variable_names"),
                        R"doc(

Function for creating control surface aerodynamic model settings from custom coefficients.

Function for settings object, defining control surface aerodynamic interface model via a custom force and moment coefficient function
(function of independent variable). This function is essentically the control-surface equivalent of the
:func:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.custom_aerodynamic_force_and_moment_coefficients` function for body coefficient settings.


Parameters
----------
force_and_moment_coefficient_function : callable[[list[float]], numpy.ndarray[numpy.float64[6, 1]]]
    Function that is defining the aerodynamic force (first three entries) and moment (last three entries) coefficients as function of an independent variables (see  ``independent_variable_names``).
independent_variable_names : list[environment.AerodynamicCoefficientsIndependentVariables]
    Vector with identifiers for the independent variable w.r.t. which the control surface aerodynamic coefficients are defined. Typically, one entry from this list will be ``control_surface_deflection_dependent``
Returns
-------
ControlSurfaceIncrementAerodynamicCoefficientSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.ControlSurfaceIncrementAerodynamicCoefficientSettings` derived class





Examples
--------
In this example, we create :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.ControlSurfaceIncrementAerodynamicCoefficientSettings` for the artificial body "Vehicle", using a function based on the mach number:

.. code-block:: python

  def force_and_moment_coefficients(variables_list):
    # Extract the mach number
    mach_number = variables_list[0]
    # If the mach number is below 3, use fixed coefficients
    if mach_number <= 3:
        return [0.99, 0, 1.08, 0.94, 0.35, 0, 0]
    # Same if the mach number is above 10
    elif mach_number >= 10:
        return [0.82, 0, 0.88, 0.55, 0, 0]
    # Otherwise, vary linearly between the ones at M=3 and M=10
    CD = 1.0667-0.02457*mach_number
    CL = 1.1636-0.02786*mach_number
    Cl = 0.35 + 0.02857*mach_number
    return [CD, 0, CL, Cl, 0, 0]
  # Create the aerodynamic interface settings
  aero_coefficient_settings = environment_setup.aerodynamic_coefficients.custom_control_surface(
      force_and_moment_coefficients,
      independent_variable_names=[environment.AerodynamicCoefficientsIndependentVariables.mach_number_dependent]
  )
  # Assign the aerodynamic coefficient interface to the vehicle
  environment_setup.add_aerodynamic_coefficient_interface(bodies, "Vehicle", aero_coefficient_settings)



    )doc");

                    m.def(
                        "tabulated_from_files_control_surface",
                        py::overload_cast<
                            const std::map<int, std::string>,
                            const std::map<int, std::string>,
                            const std::vector<
                                ta::AerodynamicCoefficientsIndependentVariables>>(
                            &tss::
                                readTabulatedControlIncrementAerodynamicCoefficientsFromFiles),
                        py::arg("force_coefficient_files"),
                        py::arg("moment_coefficient_files"),
                        py::arg("independent_variable_names"),
                        R"doc(

Function for creating control surface aerodynamic model settings from tabulated coefficients from files.

Function for settings object, defining control surface aerodynamic interface model via user-defined, tabulated aerodynamic force and moment coefficients
(tabulated w.r.t. independent variable), obtained from data files.. This function is essentically the control-surface equivalent of the
:func:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.tabulated_from_files` function for body coefficient settings.

Returns
-------
ControlSurfaceIncrementAerodynamicCoefficientSettings
    Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.ControlSurfaceIncrementAerodynamicCoefficientSettings` derived class






    )doc");


                    /////////////////////////////////////////////////////////////////
                    //////////////// DEPRECATED
                    ////////////////////////////////////////
                    /////////////////////////////////////////////////////////////////

                    // m.def("custom",
                    //       &tss::customAerodynamicCoefficientSettingsDeprecatedPy,
                    //       py::arg("force_coefficient_function"),
                    //       py::arg("reference_area"),
                    //       py::arg("independent_variable_names"),
                    //       py::arg("are_coefficients_in_aerodynamic_frame") =
                    //       true,
                    //       py::arg("are_coefficients_in_negative_axis_direction")
                    //       = true );

                    // m.def("constant",
                    //       py::overload_cast<const double, const
                    //       Eigen::Vector3d &, const bool,
                    //               const bool>(
                    //               &tss::constantAerodynamicCoefficientSettingsDeprecated
                    //               ),
                    //       py::arg("reference_area"),
                    //       py::arg("constant_force_coefficient"),
                    //       py::arg("are_coefficients_in_aerodynamic_frame") =
                    //       true,
                    //       py::arg("are_coefficients_in_negative_axis_direction")
                    //       = true );

                    // m.def("custom_aerodynamic_force_coefficients",
                    //       py::overload_cast<
                    //               const std::function<Eigen::Vector3d(const
                    //               std::vector<double> &)>, const double,
                    //               const
                    //               std::vector<ta::AerodynamicCoefficientsIndependentVariables>,
                    //               const bool, const
                    //               bool>(&tss::customAerodynamicCoefficientSettingsDeprecated),
                    //       py::arg("force_coefficient_function"),
                    //       py::arg("reference_area"),
                    //       py::arg("independent_variable_names"),
                    //       py::arg("are_coefficients_in_aerodynamic_frame") =
                    //       true,
                    //       py::arg("are_coefficients_in_negative_axis_direction")
                    //       = true );

                    // m.def("custom_aerodynamic_force_and_moment_coefficients",
                    //       py::overload_cast<
                    //               const std::function< Eigen::Vector3d( const
                    //               std::vector< double >& ) >, const
                    //               std::function< Eigen::Vector3d( const
                    //               std::vector< double >& ) >, const double,
                    //               const double,
                    //               const Eigen::Vector3d&,
                    //               const std::vector<
                    //               ta::AerodynamicCoefficientsIndependentVariables
                    //               >, const bool, const bool
                    //               >(&tss::customAerodynamicCoefficientSettingsDeprecated),
                    //       py::arg("force_coefficient_function"),
                    //       py::arg("moment_coefficient_function"),
                    //       py::arg("reference_length"),
                    //       py::arg("reference_area"),
                    //       py::arg("moment_reference_point"),
                    //       py::arg("independent_variable_names"),
                    //       py::arg("are_coefficients_in_aerodynamic_frame") =
                    //       true,
                    //       py::arg("are_coefficients_in_negative_axis_direction")
                    //       = true );


                    // m.def("tabulated",
                    //       py::overload_cast<
                    //               const std::vector<double>,
                    //               const std::vector<Eigen::Vector3d>,
                    //               const std::vector<Eigen::Vector3d>,
                    //               const double,
                    //               const double,
                    //               const double,
                    //               const Eigen::Vector3d &,
                    //               const
                    //               ta::AerodynamicCoefficientsIndependentVariables,
                    //               const bool,
                    //               const bool,
                    //               const
                    //               std::shared_ptr<ti::InterpolatorSettings>>
                    //               (&tss::oneDimensionalTabulatedAerodynamicCoefficientSettingsDeprecated),
                    //       py::arg("independent_variables"),
                    //       py::arg("force_coefficients"),
                    //       py::arg("moment_coefficients"),
                    //       py::arg("reference_length"),
                    //       py::arg("reference_area"),
                    //       py::arg("lateral_reference_length"),
                    //       py::arg("moment_reference_point"),
                    //       py::arg("independent_variable_name"),
                    //       py::arg("are_coefficients_in_aerodynamic_frame") =
                    //       true,
                    //       py::arg("are_coefficients_in_negative_axis_direction")
                    //       = true, py::arg("interpolator_settings") =
                    //       nullptr);

                    // m.def("tabulated_force_only",
                    //       py::overload_cast<
                    //               const std::vector<double>,
                    //               const std::vector<Eigen::Vector3d>,
                    //               const double,
                    //               const
                    //               ta::AerodynamicCoefficientsIndependentVariables,
                    //               const bool,
                    //               const bool,
                    //               const
                    //               std::shared_ptr<ti::InterpolatorSettings>>
                    //               (&tss::oneDimensionalTabulatedAerodynamicCoefficientSettingsDeprecated),
                    //       py::arg("independent_variables"),
                    //       py::arg("force_coefficients"),
                    //       py::arg("reference_area"),
                    //       py::arg("independent_variable_name"),
                    //       py::arg("are_coefficients_in_aerodynamic_frame") =
                    //       true,
                    //       py::arg("are_coefficients_in_negative_axis_direction")
                    //       = true, py::arg("interpolator_settings") =
                    //       nullptr);

                    // m.def("tabulated_force_only_from_files",
                    //       py::overload_cast<
                    //               const std::map< int, std::string >,
                    //               const double,
                    //               const std::vector<
                    //               ta::AerodynamicCoefficientsIndependentVariables
                    //               >, const bool, const bool, const
                    //               std::shared_ptr< ti::InterpolatorSettings >
                    //               >
                    //               (&tss::readTabulatedAerodynamicCoefficientsFromFilesDeprecated),
                    //       py::arg("force_coefficient_files"),
                    //       py::arg("reference_area"),
                    //       py::arg("independent_variable_names"),
                    //       py::arg("are_coefficients_in_aerodynamic_frame") =
                    //       true,
                    //       py::arg("are_coefficients_in_negative_axis_direction")
                    //       = true, py::arg("interpolator_settings") = nullptr
                    //       );

                    // m.def("tabulated_from_files",
                    //       py::overload_cast<
                    //               const std::map< int, std::string >,
                    //               const std::map< int, std::string >,
                    //               const double,
                    //               const double,
                    //               const double,
                    //               const Eigen::Vector3d &,
                    //               const std::vector<
                    //               ta::AerodynamicCoefficientsIndependentVariables
                    //               >, const bool, const bool, const
                    //               std::shared_ptr< ti::InterpolatorSettings >
                    //               >
                    //               (&tss::readTabulatedAerodynamicCoefficientsFromFilesDeprecated),
                    //       py::arg("force_coefficient_files"),
                    //       py::arg("moment_coefficient_files"),
                    //       py::arg("reference_length"),
                    //       py::arg("reference_area"),
                    //       py::arg("lateral_reference_length"),
                    //       py::arg("moment_reference_point"),
                    //       py::arg("independent_variable_names"),
                    //       py::arg("are_coefficients_in_aerodynamic_frame") =
                    //       true,
                    //       py::arg("are_coefficients_in_negative_axis_direction")
                    //       = true, py::arg("interpolator_settings") = nullptr
                    //       );
                }


            }  // namespace aerodynamic_coefficients
        }      // namespace environment_setup
    }          // namespace numerical_simulation
}  // namespace tudatpy
