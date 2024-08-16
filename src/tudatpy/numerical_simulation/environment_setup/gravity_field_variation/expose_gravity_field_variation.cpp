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
#include <tudat/simulation/environment_setup.h>

#include "tudatpy/docstrings.h"

// #include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
// #include <pybind11/numpy.h>
#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace tss = tudat::simulation_setup;
namespace tg = tudat::gravitation;

namespace tudat {

    namespace simulation_setup {

        inline std::shared_ptr<GravityFieldVariationSettings>
        degreeOrderVariableLoveNumberGravityFieldVariationSettingsPy(
            const std::vector<std::string>& deformingBodies,
            const std::map<int, std::vector<double>> loveNumber) {
            std::map<int, std::vector<std::complex<double>>> loveNumbers;
            for(auto loveNumberIt : loveNumber) {
                for(unsigned int i = 0; i < loveNumberIt.second.size(); i++) {
                    loveNumbers[loveNumberIt.first].push_back(
                        std::complex<double>(loveNumberIt.second.at(i), 0));
                }
            }
            return std::make_shared<
                BasicSolidBodyGravityFieldVariationSettings>(
                deformingBodies, loveNumbers, nullptr);
        }

        inline std::shared_ptr<GravityFieldVariationSettings>
        degreeOrderVariableLoveNumberGravityFieldVariationSettingsPy(
            const std::string deformingBody,
            const std::map<int, std::vector<double>> loveNumber) {
            std::map<int, std::vector<std::complex<double>>> loveNumbers;
            for(auto loveNumberIt : loveNumber) {
                for(unsigned int i = 0; i < loveNumberIt.second.size(); i++) {
                    loveNumbers[loveNumberIt.first].push_back(
                        std::complex<double>(loveNumberIt.second.at(i), 0));
                }
            }
            return std::make_shared<
                BasicSolidBodyGravityFieldVariationSettings>(
                std::vector<std::string>({deformingBody}), loveNumbers,
                nullptr);
        }

        inline std::shared_ptr<GravityFieldVariationSettings>
        degreeOrderVariableLoveNumberGravityFieldVariationSettingsPy(
            const std::string deformingBody,
            const std::map<int, std::vector<std::complex<double>>> loveNumber) {
            return std::make_shared<
                BasicSolidBodyGravityFieldVariationSettings>(
                std::vector<std::string>({deformingBody}), loveNumber, nullptr);
        }

    }  // namespace simulation_setup

}  // namespace tudat
namespace tudatpy {
    namespace numerical_simulation {
        namespace environment_setup {
            namespace gravity_field_variation {

                PYBIND11_MODULE(expose_gravity_field_variation, m) {
                    py::enum_<tg::BodyDeformationTypes>(
                        m, "BodyDeformationTypes", "<no_doc>")
                        .value("basic_solid_body", tg::basic_solid_body)
                        .value("tabulated_deformation", tg::tabulated_variation)
                        .export_values();

                    py::class_<
                        tss::GravityFieldVariationSettings,
                        std::shared_ptr<tss::GravityFieldVariationSettings>>(
                        m, "GravityFieldVariationSettings",
R"doc(Base class for providing settings for gravity field variations.

)doc");


                    py::class_<
                        tss::BasicSolidBodyGravityFieldVariationSettings,
                        std::shared_ptr<
                            tss::BasicSolidBodyGravityFieldVariationSettings>,
                        tss::GravityFieldVariationSettings>(
                        m, "BasicSolidBodyGravityFieldVariationSettings",
R"doc(Class for providing settings for solid body tidal gravity field variations, derived from GravityFieldVariationSettings.

)doc");

                    m.def(
                        "solid_body_tide",
                        py::overload_cast<const std::string, const double,
                                          const int>(
                            &tss::
                                fixedSingleDegreeLoveNumberGravityFieldVariationSettings),
                        py::arg("tide_raising_body"), py::arg("love_number"),
                        py::arg("degree"),
R"doc(Factory function for creating solid body tides.

	Factory function for creating solid body tides, using a single real Love number at a single degree (e.g. :math:`k_{2}`, :math:`k_{3}`, etc.). This function evaluates Eq. (6.6) from the IERS Conventions 2010, with real :math:`k_{l}=k_{lm}`, a single value of :math:`l` and a single tide-raising body :math:`j`.
	

	:param tide_raising_body:
		Name of body raising the tide.
	:param love_number:
		Constant real Love number to use for body undergoing deformation, at the spherical harmonic degree defined by 'degree' input.
	:param degree:
		Degree of the spherical harmonic gravity field, and associated Love number, that is to be considered.
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.BasicSolidBodyGravityFieldVariationSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.GravityFieldVariationSettings` class
)doc");

                    m.def(
                        "solid_body_tide_complex_k",
                        py::overload_cast<const std::string,
                                          const std::complex<double>,
                                          const int>(
                            &tss::
                                fixedSingleDegreeLoveNumberGravityFieldVariationSettings),
                        py::arg("tide_raising_body"), py::arg("love_number"),
                        py::arg("degree"),
R"doc(Factory function for creating solid body tides.

	As :func:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.solid_body_tide`, but with complex value for the Love number.
	

	:param tide_raising_body:
		Name of body raising the tide.
	:param love_number:
		Constant real Love number to use for body undergoing deformation, at the spherical harmonic degree defined by 'degree' input.
	:param degree:
		Degree of the spherical harmonic gravity field, and associated Love number, that is to be considered.
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.BasicSolidBodyGravityFieldVariationSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.GravityFieldVariationSettings` class
)doc");

                    m.def(
                        "solid_body_tide_degree_variable_k",
                        py::overload_cast<const std::string,
                                          std::map<int, double>>(
                            &tss::
                                fixedSingleDegreeLoveNumberGravityFieldVariationSettings),
                        py::arg("tide_raising_body"),
                        py::arg("love_number_per_degree"),
R"doc(Factory function for creating solid body tides.

	Factory function for creating solid body tides, using a set of real, separate, Love numbers at any number of degrees (e.g. :math:`k_{2}`, :math:`k_{3}`, etc.). This output of this function is effectively identical to a list of outputs to :func:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.solid_body_tide`, with differing degrees and associated Love numbers.  This function evaluates Eq. (6.6) from the IERS Conventions 2010, with a set of real values :math:`k_{l}=k_{lm}`, at a set of values of :math:`l` and a single tide-raising body :math:`j`.
	

	:param tide_raising_body:
		Name of body raising the tide.
	:param love_number_per_degree:
		Dictionary of Love numbers for each degree that is to be taken into account, with the key representing the degree :math:`l` of the Love number, and value containing the Love number :math:`k_{l}` itself
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.BasicSolidBodyGravityFieldVariationSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.GravityFieldVariationSettings` class
)doc");


                    m.def(
                        "solid_body_tide_degree_variable_complex_k",
                        py::overload_cast<const std::string,
                                          std::map<int, std::complex<double>>>(
                            &tss::
                                fixedSingleDegreeLoveNumberGravityFieldVariationSettings),
                        py::arg("tide_raising_body"),
                        py::arg("love_number_per_degree"),
R"doc(Factory function for creating solid body tides.

	As :func:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.solid_body_tide_degree_variable_k`, but with complex values for the Love numbers.
	

	:param tide_raising_body:
		Name of body raising the tide.
	:param love_number_per_degree:
		Dictionary of Love numbers for each degree that is to be taken into account, with the key representing the degree :math:`l` of the Love number, and value containing the Love number :math:`k_{l}` itself.
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.BasicSolidBodyGravityFieldVariationSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.GravityFieldVariationSettings` class
)doc");

                    m.def(
                        "solid_body_tide_degree_order_variable_k",
                        py::overload_cast<
                            const std::string,
                            const std::map<int, std::vector<double>>>(
                            &tss::
                                degreeOrderVariableLoveNumberGravityFieldVariationSettingsPy),
                        py::arg("tide_raising_body"),
                        py::arg("love_number_per_degree_and_order"),
R"doc(Factory function for creating solid body tides.

	Factory function for creating solid body tides, using a set of real, separate, Love numbers at any number of degrees and orders (e.g. :math:`k_{20}`, :math:`k_{21}`, :math:`k_{22}`, :math:`k_{30}`, etc.).  This function evaluates Eq. (6.6) from the IERS Conventions 2010, with a set of real values :math:`k_{lm}`, at a set of values of :math:`l` and a single tide-raising body :math:`j`.
	

	:param tide_raising_body:
		Name of body raising the tide.
	:param love_number_per_degree_and_order:
		Dictionary of Love numbers for each degree that is to be taken into account, with the key representing the degree :math:`l` of the Love number, and value containing the list of Love numbers :math:`k_{lm}` at this degree. Note that, for Love numbers at degree :math:`l`, the associated list should contain :math:`l+1` entries, representing the Love numbers (in order) :math:`k_{l0}`, :math:`k_{l1}`...:math:`k_{ll}`.
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.BasicSolidBodyGravityFieldVariationSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.GravityFieldVariationSettings` class
)doc");

                    m.def(
                        "solid_multi_body_tide_degree_order_variable_k",
                        py::overload_cast<
                            const std::vector<std::string>&,
                            const std::map<int, std::vector<double>>>(
                            &tss::
                                degreeOrderVariableLoveNumberGravityFieldVariationSettingsPy),
                        py::arg("tide_raising_bodies"),
                        py::arg("love_number_per_degree_and_order"),
get_docstring("solid_multi_body_tide_degree_order_variable_k").c_str());

                    m.def(
                        "solid_body_tide_degree_order_variable_complex_k",
                        py::overload_cast<
                            const std::string,
                            const std::map<int,
                                           std::vector<std::complex<double>>>>(
                            &tss::
                                degreeOrderVariableLoveNumberGravityFieldVariationSettingsPy),
                        py::arg("tide_raising_body"),
                        py::arg("love_number_per_degree_and_order"),
R"doc(Factory function for creating solid body tides.

	As :func:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.solid_body_tide_degree_order_variable_k`, but with complex values for the Love number.
	

	:param tide_raising_body:
		Name of body raising the tide.
	:param love_number_per_degree:
		Dictionary of Love numbers for each degree that is to be taken into account, with the key representing the degree :math:`l` of the Love number, and value containing the list of Love numbers :math:`k_{lm}` at this degree. Note that, for Love numbers at degree :math:`l`, the associated list should contain :math:`l+1` entries, representing the Love numbers (in order) :math:`k_{l0}`, :math:`k_{l1}`...:math:`k_{ll}`.
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.BasicSolidBodyGravityFieldVariationSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.GravityFieldVariationSettings` class
)doc");

                    m.def(
                        "solid_body_tide_degree_order_variable_complex_k",
                        py::overload_cast<
                            const std::string,
                            const std::map<int,
                                           std::vector<std::complex<double>>>>(
                            &tss::
                                degreeOrderVariableLoveNumberGravityFieldVariationSettingsPy),
                        py::arg("tide_raising_body"),
                        py::arg("love_number_per_degree_and_order"),
R"doc(Factory function for creating solid body tides.

	As :func:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.solid_body_tide_degree_order_variable_k`, but with complex values for the Love number.
	

	:param tide_raising_body:
		Name of body raising the tide.
	:param love_number_per_degree:
		Dictionary of Love numbers for each degree that is to be taken into account, with the key representing the degree :math:`l` of the Love number, and value containing the list of Love numbers :math:`k_{lm}` at this degree. Note that, for Love numbers at degree :math:`l`, the associated list should contain :math:`l+1` entries, representing the Love numbers (in order) :math:`k_{l0}`, :math:`k_{l1}`...:math:`k_{ll}`.
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.BasicSolidBodyGravityFieldVariationSettings` derived :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field_variation.GravityFieldVariationSettings` class
)doc");

                    m.def(
                        "single_period_periodic",
                        &tss::
                            periodicGravityFieldVariationsSettingsSingleFrequency,
                        py::arg("cosine_coefficient_amplitude_cosine_time"),
                        py::arg("cosine_coefficient_amplitude_sine_time"),
                        py::arg("sine_coefficient_amplitude_cosine_time"),
                        py::arg("sine_coefficient_amplitude_sine_time"),
                        py::arg("frequency"), py::arg("reference_epoch"),
                        py::arg("minimum_degree") = 2,
                        py::arg("minimum_order") = 0,
get_docstring("single_period_periodic").c_str());

                    m.def("periodic",
                          &tss::periodicGravityFieldVariationsSettings,
                          py::arg("cosine_coefficient_amplitudes_cosine_time"),
                          py::arg("cosine_coefficient_amplitudes_sine_time"),
                          py::arg("sine_coefficient_amplitudes_cosine_time"),
                          py::arg("sine_coefficient_amplitudes_sine_time"),
                          py::arg("frequencies"), py::arg("reference_epoch"),
                          py::arg("minimum_degree") = 2,
                          py::arg("minimum_order") = 0,
get_docstring("single_period_periodic").c_str());

                    m.def(
                        "single_power_polynomial",
                        &tss::
                            polynomialGravityFieldVariationsSettingsSinglePower,
                        py::arg("cosine_amplitudes"),
                        py::arg("sine_amplitudes"), py::arg("polynomial_power"),
                        py::arg("reference_epoch"),
                        py::arg("minimum_degree") = 2,
                        py::arg("minimum_order") = 0,
get_docstring("single_power_polynomial").c_str());

                    m.def("polynomial",
                          &tss::polynomialGravityFieldVariationsSettings,
                          py::arg("cosine_amplitudes_per_power"),
                          py::arg("sine_amplitudes_per_power"),
                          py::arg("reference_epoch"),
                          py::arg("minimum_degree") = 2,
                          py::arg("minimum_order") = 0,
get_docstring("polynomial").c_str());


                    m.def("tabulated",
                          &tss::tabulatedGravityFieldVariationSettings,
                          py::arg("cosine_variations_table"),
                          py::arg("sine_variations_table"),
                          py::arg("minimum_degree"), py::arg("minimum_order"),
                          py::arg("interpolation_settings"));
                }

            }  // namespace gravity_field_variation
        }  // namespace environment_setup
    }  // namespace numerical_simulation
}  // namespace tudatpy
