/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <tudat/astro/reference_frames/referenceFrameTransformations.h>
#include <tudat/basics/deprecationWarnings.h>
#include <tudat/simulation/environment_setup.h>


// #include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
// #include <pybind11/numpy.h>
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
            const std::function<Eigen::Vector3d(const std::vector<double>&)>
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

                PYBIND11_MODULE(expose_aerodynamic_coefficients, m) {
                    py::module_::import(
                        "tudatpy.numerical_simulation.environment");
                    /////////////////////////////////////////////////////////////////////////////
                    // createAerodynamicCoefficientInterface.h
                    /////////////////////////////////////////////////////////////////////////////
                    py::class_<
                        tss::AerodynamicCoefficientSettings,
                        std::shared_ptr<tss::AerodynamicCoefficientSettings>>(
                        m, "AerodynamicCoefficientSettings", "")
                        .def_property("add_force_contribution_to_moments",
                                      &tss::AerodynamicCoefficientSettings::
                                          getAddForceContributionToMoments,
                                      &tss::AerodynamicCoefficientSettings::
                                          setAddForceContributionToMoments,
                                      "")
                        .def_property("moment_reference_point",
                                      &tss::AerodynamicCoefficientSettings::
                                          getMomentReferencePoint,
                                      &tss::AerodynamicCoefficientSettings::
                                          setMomentReferencePoint,
                                      "")
                        .def("add_single_control_surface",
                             &tss::AerodynamicCoefficientSettings::
                                 addControlSurfaceSettings,
                             py::arg("control_surface_settings"),
                             py::arg("control_surface_name"), "");

                    py::class_<tss::ConstantAerodynamicCoefficientSettings,
                               std::shared_ptr<
                                   tss::ConstantAerodynamicCoefficientSettings>,
                               tss::AerodynamicCoefficientSettings>(
                        m, "ConstantAerodynamicCoefficientSettings",
                        R"doc(Class for defining model settings from constant aerodynamic coefficients.

	`AerodynamicCoefficientSettings` derived class for aerodynamic interface model settings using only constant aerodynamic coefficients.
)doc");

                    py::class_<
                        tss::
                            ControlSurfaceIncrementAerodynamicCoefficientSettings,
                        std::shared_ptr<
                            tss::
                                ControlSurfaceIncrementAerodynamicCoefficientSettings>>(
                        m,
                        "ControlSurfaceIncrementAerodynamicCoefficientSettings",
                        "");

                    m.def(
                        "constant",
                        py::overload_cast<
                            const double, const Eigen::Vector3d&,
                            const ta::AerodynamicCoefficientFrames>(
                            &tss::constantAerodynamicCoefficientSettings),
                        py::arg("reference_area"),
                        py::arg("constant_force_coefficient"),
                        py::arg("force_coefficients_frame") =
                            ta::negative_aerodynamic_frame_coefficients,
                        R"doc(Factory function for creating aerodynamic interface model settings entirely from constant coefficients.

	Factory function for settings object, defining aerodynamic interface model entirely from constant aerodynamic coefficients,
	i.e. coefficients are not a function of any independent variables.


	:param reference_area:
		Reference area with which aerodynamic forces and moments are non-dimensionalized.
	:param constant_force_coefficient:
		Constant force coefficients.
	:param are_coefficients_in_aerodynamic_frame:
		Boolean to define whether the aerodynamic coefficients are defined in the aerodynamic frame
		(drag, side, lift force) or in the body frame (typically denoted as Cx, Cy, Cz).

	:param are_coefficients_in_negative_axis_direction:
		Boolean to define whether the aerodynamic coefficients are positive along the positive axes of the body or
		aerodynamic frame (see arg are_coefficients_in_aerodynamic_frame).
		Note that for drag, side and lift force, the coefficients are typically defined in negative direction.

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.ConstantAerodynamicCoefficientSettings` class
)doc");


                    m.def(
                        "custom_aerodynamic_force_coefficients",
                        py::overload_cast<
                            const std::function<Eigen::Vector3d(
                                const std::vector<double>&)>,
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
                        R"doc(Factory function for creating aerodynamic interface model settings from custom coefficients.

	Factory function for settings object, defining aerodynamic interface model via a custom force coefficient function
	(function of independent variable).


	:param force_coefficient_function:
		Function that is defining the aerodynamic coefficients as function of an independent variable (see arg independent_variable_names).
	:param reference_area:
		Reference area with which aerodynamic forces and moments are non-dimensionalized.
	:param independent_variable_name:
		Vector with identifiers for the independent variable w.r.t. which the aerodynamic coefficients are defined.
	:param are_coefficients_in_aerodynamic_frame:
		Boolean to define whether the aerodynamic coefficients are defined in the aerodynamic frame
		(drag, side, lift force) or in the body frame (typically denoted as Cx, Cy, Cz).

	:param are_coefficients_in_negative_axis_direction:
		Boolean to define whether the aerodynamic coefficients are positive along the positive axes of the body or
		aerodynamic frame (see arg are_coefficients_in_aerodynamic_frame).
		Note that for drag, side and lift force, the coefficients are typically defined in negative direction.

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.CustomAerodynamicCoefficientSettings` class
)doc");

                    m.def(
                        "custom_aerodynamic_force_and_moment_coefficients",
                        py::overload_cast<
                            const std::function<Eigen::Vector3d(
                                const std::vector<double>&)>,
                            const std::function<Eigen::Vector3d(
                                const std::vector<double>&)>,
                            const double, const double,
                            const std::vector<
                                ta::AerodynamicCoefficientsIndependentVariables>,
                            const ta::AerodynamicCoefficientFrames,
                            const ta::AerodynamicCoefficientFrames,
                            const Eigen::Vector3d&>(
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
                        R"doc(Factory function for creating aerodynamic interface model settings from custom coefficients.

	Factory function for settings object, defining aerodynamic interface model via a custom force coefficient function
	(function of independent variable).


	:param force_coefficient_function:
		Function that is defining the aerodynamic force coefficients as function of an independent variable (see arg independent_variable_names).
	:param moment_coefficient_function:
		Function that is defining the aerodynamic moment coefficients as function of an independent variable (see arg independent_variable_names).
	:param reference_area:
		Reference area with which aerodynamic forces and moments are non-dimensionalized.
	:param reference_length:
		Reference length with which aerodynamic moments are non-dimensionalized.
	:param moment_reference_point:
		Reference point in the body-fixed frame w.r.t. which the moments are provided.
	:param independent_variable_name:
		Vector with identifiers for the independent variable w.r.t. which the aerodynamic coefficients are defined.
	:param are_coefficients_in_aerodynamic_frame:
		Boolean to define whether the aerodynamic coefficients are defined in the aerodynamic frame
		(drag, side, lift force) or in the body frame (typically denoted as Cx, Cy, Cz).

	:param are_coefficients_in_negative_axis_direction:
		Boolean to define whether the aerodynamic coefficients are positive along the positive axes of the body or
		aerodynamic frame (see arg are_coefficients_in_aerodynamic_frame).
		Note that for drag, side and lift force, the coefficients are typically defined in negative direction.

	:return:
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
                            const Eigen::Vector3d&,
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
                        R"doc(Factory function for creating aerodynamic interface model settings from user-defined, 1-d tabulated coefficients.

	Factory function for settings object, defining aerodynamic interface model via user-defined, 1-dimensional, tabulated aerodynamic force and moment coefficients
	(tabulated w.r.t. independent variable).


	:param independent_variables:
		Values of independent variables at which the coefficients in the input multi vector are defined (size 1).
	:param force_coefficients:
		Values of force coefficients at independent variables defined by independent_variables.
	:param moment_coefficients:
		Values of moment coefficients at independent variables defined by independent_variables.
	:param reference_length:
		Reference length with which aerodynamic moments about x- and z- axes are non-dimensionalized.
	:param reference_area:
		Reference area with which aerodynamic forces and moments are non-dimensionalized.
	:param lateral_reference_length:
		Reference length with which aerodynamic moment about y-axis is non-dimensionalized.
	:param moment_reference_point:
		Point w.r.t. aerodynamic moment is calculated.
	:param independent_variable_name:
		Vector with identifiers for the independent variable w.r.t. which the aerodynamic coefficients are defined.
	:param are_coefficients_in_aerodynamic_frame:
		Boolean to define whether the aerodynamic coefficients are defined in the aerodynamic frame
		(drag, side, lift force) or in the body frame (typically denoted as Cx, Cy, Cz).

	:param are_coefficients_in_negative_axis_direction:
		Boolean to define whether the aerodynamic coefficients are positive along the positive axes of the body or
		aerodynamic frame (see arg areCoefficientsInAerodynamicFrame).
		Note that for drag, side and lift force, the coefficients are typically defined in negative direction.

	:param interpolator_settings:
		Interpolator settings object, where the conditions for interpolation of tabulated inputs are saved.
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.TabulatedAerodynamicCoefficientSettings` class (via :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.TabulatedAerodynamicCoefficientSettingsBase` class)
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
                        R"doc(Factory function for creating aerodynamic interface model settings from user-defined, 1-d tabulated force coefficients.

	Factory function for settings object, defining aerodynamic interface model via user-defined, 1-dimensional, tabulated aerodynamic force coefficients
	(tabulated w.r.t. independent variable).


	:param independent_variables:
		Values of independent variables at which the coefficients in the input multi vector are defined (size 1)
	:param force_coefficients:
		Values of force coefficients at independent variables defined by independent_variables.
	:param reference_area:
		Reference area with which aerodynamic forces and moments are non-dimensionalized.
	:param independent_variable_name:
		Identifier of the independent variable w.r.t. which the aerodynamic coefficients are defined.
	:param are_coefficients_in_aerodynamic_frame:
		Boolean to define whether the aerodynamic coefficients are defined in the aerodynamic frame
		(drag, side, lift force) or in the body frame (typically denoted as Cx, Cy, Cz).

	:param are_coefficients_in_negative_axis_direction:
		Boolean to define whether the aerodynamic coefficients are positive along the positive axes of the body or
		aerodynamic frame (see arg areCoefficientsInAerodynamicFrame).
		Note that for drag, side and lift force, the coefficients are typically defined in negative direction.

	:param interpolator_settings:
		Interpolator settings object, where the conditions for interpolation of tabulated inputs are saved.
		Pointer to an interpolator settings object where the conditions for interpolation of tabulated inputs are saved.

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.TabulatedAerodynamicCoefficientSettings` class
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
                        R"doc(Factory function for creating aerodynamic interface model settings from tabulated force coefficients from files.

	Factory function for settings object, defining aerodynamic interface model via user-defined, tabulated aerodynamic force coefficients
	(tabulated w.r.t. independent variable), obtained from data files.


	:param force_coefficient_files:
		Path of the aerodynamic coefficient files corresponding to the force coefficient of the given dict key.
	:param reference_area:
		Reference area with which aerodynamic forces and moments are non-dimensionalized.
	:param independent_variable_names:
		Vector with identifiers for the independent variable w.r.t. which the aerodynamic coefficients are defined.
	:param are_coefficients_in_aerodynamic_frame:
		Boolean to define whether the aerodynamic coefficients are defined in the aerodynamic frame
		(drag, side, lift force) or in the body frame (typically denoted as Cx, Cy, Cz).

	:param are_coefficients_in_negative_axis_direction:
		Boolean to define whether the aerodynamic coefficients are positive along the positive axes of the body or
		aerodynamic frame (see arg areCoefficientsInAerodynamicFrame).
		Note that for drag, side and lift force, the coefficients are typically defined in negative direction.

	:param interpolator_settings:
		Interpolator settings object, where the conditions for interpolation of tabulated inputs are saved.

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.TabulatedAerodynamicCoefficientSettings` class
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
                            const Eigen::Vector3d&,
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
                        R"doc(Factory function for creating aerodynamic interface model settings from tabulated coefficients from files.

	Factory function for settings object, defining aerodynamic interface model via user-defined, tabulated aerodynamic force and moment coefficients
	(tabulated w.r.t. independent variable), obtained from data files.


	:param force_coefficient_files:
		Path of the aerodynamic coefficient files corresponding to the force coefficient of the given dict key.
	:param moment_coefficient_files:
		Path of the aerodynamic coefficient files corresponding to the moment coefficient of the given dict key.
	:param reference_length:
		Reference length with which aerodynamic moments about x- and z- axes are non-dimensionalized.
	:param reference_area:
		Reference area with which aerodynamic forces and moments are non-dimensionalized.
	:param lateral_reference_length:
		Reference length with which aerodynamic moment about y-axis is non-dimensionalized.
	:param moment_reference_point:
		Point w.r.t. aerodynamic moment is calculated.
	:param independent_variable_names:
		Vector with identifiers for the independent variable w.r.t. which the aerodynamic coefficients are defined.
	:param are_coefficients_in_aerodynamic_frame:
		Boolean to define whether the aerodynamic coefficients are defined in the aerodynamic frame
		(drag, side, lift force) or in the body frame (typically denoted as Cx, Cy, Cz).

	:param are_coefficients_in_negative_axis_direction:
		Boolean to define whether the aerodynamic coefficients are positive along the positive axes of the body or
		aerodynamic frame (see arg areCoefficientsInAerodynamicFrame).
		Note that for drag, side and lift force, the coefficients are typically defined in negative direction.

	:param interpolator_settings:
		Interpolator settings object, where the conditions for interpolation of tabulated inputs are saved.

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.TabulatedAerodynamicCoefficientSettings` class
)doc");

                    m.def(
                        "scaled_by_constant",
                        py::overload_cast<
                            const std::shared_ptr<
                                tss::AerodynamicCoefficientSettings>,
                            const double, const double, const bool>(
                            &tss::scaledAerodynamicCoefficientSettings),
                        py::arg("unscaled_coefficient_settings"),
                        py::arg("force_scaling_constant"),
                        py::arg("moment_scaling_constant"),
                        py::arg("is_scaling_absolute") = false,
                        R"doc(Factory function for creating aerodynamic interface model settings by applying one constant scaling factor/value to all coefficients of an existing model settings object.

	Factory function for settings object, defining aerodynamic interface based on scaling the coefficients of an existing model settings object by one constant factor or value.
	Via the ``is_scaling_absolute``
	boolean, the user can apply a constant scaling factor or an absolute value to the resulting force and moment coefficients (for instance for an uncertainty analysis).


	:param unscaled_coefficient_settings:
		Existing aerodynamic interface model settings object that is used as the base for the scaled settings object.
	:param force_scaling_constant:
		Constant scaling factor to be applied to all aerodynamic force coefficients.
	:param moment_scaling_constant:
		Constant scaling factor to be applied to all aerodynamic moment coefficients.
	:param is_scaling_absolute, default = False:
		Boolean indicating whether aerodynamic coefficient scaling is absolute.
		Setting this boolean to true will add the scaling value to the base value,
		instead of the default behaviour of multiplying the base value by the scaling factor.

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.ScaledAerodynamicCoefficientInterfaceSettings` class
)doc");

                    m.def(
                        "scaled_by_vector",
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
                        R"doc(Factory function for creating aerodynamic interface model settings by applying constant scaling factors/values to the coefficients of an existing model settings object.

	Factory function for settings object, defining aerodynamic interface based on scaling the coefficients of an existing model settings object by constant factors or values.
	Via the ``is_scaling_absolute`` boolean, the user can apply one constant scaling factor or an absolute value to each resulting force and moment coefficient (for instance for an uncertainty analysis).


	:param unscaled_coefficient_settings:
		Existing aerodynamic interface model settings object that is used as the base for the scaled settings object.
	:param force_scaling_vector:
		Constant scaling factors to be applied to each aerodynamic force coefficient.
	:param moment_scaling_vector:
		Constant scaling factors to be applied to each aerodynamic moment coefficient.
	:param is_scaling_absolute, default = False:
		Boolean indicating whether aerodynamic coefficient scaling is absolute.
		Setting this boolean to true will add the scaling value to the base value,
		instead of the default behaviour of multiplying the base value by the scaling factor.

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.ScaledAerodynamicCoefficientInterfaceSettings` class
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
                        R"doc(Factory function for creating aerodynamic interface model settings by applying custom scaling factors/values to the coefficients of an existing model settings object.

	Factory function for settings object, defining aerodynamic interface based on scaling the coefficients of an existing model settings object by custom factors or values.
	Via the ``is_scaling_absolute`` boolean, the user can apply the scaling factors or absolute values to each resulting force and moment coefficient (for instance for an uncertainty analysis).


	:param unscaled_coefficient_settings:
		Existing aerodynamic interface model settings object that is used as the base for the scaled settings object.
	:param force_scaling_vector_function:
		Custom scaling factors to be applied to each aerodynamic force coefficient.
	:param moment_scaling_vector_function:
		Custom scaling factors to be applied to each aerodynamic moment coefficient.
	:param is_scaling_absolute, default = False:
		Boolean indicating whether aerodynamic coefficient scaling is absolute.
		Setting this boolean to true will add the scaling value to the base value,
		instead of the default behaviour of multiplying the base value by the scaling factor.

	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.AerodynamicCoefficientSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.ScaledAerodynamicCoefficientInterfaceSettings` class
)doc");

                    m.def(
                        "custom_control_surface",
                        &tss::
                            customControlSurfaceIncrementAerodynamicCoefficientSettings,
                        py::arg("force_and_moment_coefficient_function"),
                        py::arg("independent_variable_names"), "");

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
                        py::arg("independent_variable_names"), "");


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
        }  // namespace environment_setup
    }  // namespace numerical_simulation
}  // namespace tudatpy
