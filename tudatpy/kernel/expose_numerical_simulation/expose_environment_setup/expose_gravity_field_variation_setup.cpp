/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_gravity_field_variation_setup.h"

#include "tudatpy/docstrings.h"
#include <tudat/simulation/environment_setup.h>
#include <tudat/astro/reference_frames/referenceFrameTransformations.h>

//#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
//#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>

namespace py = pybind11;
namespace tss = tudat::simulation_setup;
namespace tg = tudat::gravitation;

namespace tudat
{

namespace simulation_setup
{

inline std::shared_ptr< GravityFieldVariationSettings > degreeOrderVariableLoveNumberGravityFieldVariationSettingsPy(
        const std::string deformingBody,
        const std::map< int, std::vector< double > > loveNumber )
{
    std::map< int, std::vector< std::complex< double > > > loveNumbers;
    for( auto loveNumberIt : loveNumber )
    {
        for( unsigned int i = 0; i < loveNumberIt.second.size( ); i++ )
        {
            loveNumbers[ loveNumberIt.first ].push_back( std::complex< double >( loveNumberIt.second.at( i ), 0 ) );
        }
    }
    return std::make_shared< BasicSolidBodyGravityFieldVariationSettings >(
                std::vector< std::string >( { deformingBody } ), loveNumbers, nullptr );
}

inline std::shared_ptr< GravityFieldVariationSettings > degreeOrderVariableLoveNumberGravityFieldVariationSettingsPy(
        const std::string deformingBody,
        const std::map< int, std::vector< std::complex< double > > > loveNumber )
{
    return std::make_shared< BasicSolidBodyGravityFieldVariationSettings >(
                std::vector< std::string >( { deformingBody } ), loveNumber, nullptr );
}

}

}
namespace tudatpy {
namespace numerical_simulation {
namespace environment_setup {
namespace gravity_field_variation {

void expose_gravity_field_variation_setup(py::module &m) {
    py::enum_<tg::BodyDeformationTypes>(
                m, "BodyDeformationTypes", "<no_doc>")
            .value("basic_solid_body", tg::basic_solid_body)
            .value("tabulated_deformation", tg::tabulated_variation)
            .export_values();

    py::class_<tss::GravityFieldVariationSettings,
            std::shared_ptr<tss::GravityFieldVariationSettings>>(
                m, "GravityFieldVariationSettings",
                get_docstring("GravityFieldVariationSettings").c_str( ) );


            py::class_<tss::BasicSolidBodyGravityFieldVariationSettings,
            std::shared_ptr<tss::BasicSolidBodyGravityFieldVariationSettings>,
            tss::GravityFieldVariationSettings>(
                m, "BasicSolidBodyGravityFieldVariationSettings",
                get_docstring("BasicSolidBodyGravityFieldVariationSettings").c_str( ) );

            m.def("solid_body_tide",
                  py::overload_cast<const std::string, const double, const int>(
                      &tss::fixedSingleDegreeLoveNumberGravityFieldVariationSettings),
                  py::arg("tide_raising_body"),
                  py::arg("love_number"),
                  py::arg("degree"),
                  get_docstring("solid_body_tide").c_str() );

    m.def("solid_body_tide_complex_k",
          py::overload_cast<const std::string, const std::complex<double>, const int>(
              &tss::fixedSingleDegreeLoveNumberGravityFieldVariationSettings),
          py::arg("tide_raising_body"),
          py::arg("love_number"),
          py::arg("degree"),
          get_docstring("solid_body_tide_complex_k").c_str() );

    m.def("solid_body_tide_degree_variable_k",
          py::overload_cast<const std::string, std::map<int, double> >(
              &tss::fixedSingleDegreeLoveNumberGravityFieldVariationSettings),
          py::arg("tide_raising_body"),
          py::arg("love_number_per_degree"),
          get_docstring("solid_body_tide_degree_variable_k").c_str() );


    m.def("solid_body_tide_degree_variable_complex_k",
          py::overload_cast<const std::string, std::map<int, std::complex<double> > >(
              &tss::fixedSingleDegreeLoveNumberGravityFieldVariationSettings),
          py::arg("tide_raising_body"),
          py::arg("love_number_per_degree"),
          get_docstring("solid_body_tide_degree_variable_complex_k").c_str() );

    m.def("solid_body_tide_degree_order_variable_k",
          py::overload_cast<const std::string, const std::map<int, std::vector<double> > >(
              &tss::degreeOrderVariableLoveNumberGravityFieldVariationSettingsPy),
          py::arg("tide_raising_body"),
          py::arg("love_number_per_degree_and_order"),
          get_docstring("solid_body_tide_degree_order_variable_k").c_str() );

    m.def("solid_body_tide_degree_order_variable_complex_k",
          py::overload_cast<const std::string, const std::map<int, std::vector<std::complex<double> > > >(
              &tss::degreeOrderVariableLoveNumberGravityFieldVariationSettingsPy),
          py::arg("tide_raising_body"),
          py::arg("love_number_per_degree_and_order"),
          get_docstring("solid_body_tide_degree_order_variable_complex_k").c_str() );

    m.def("solid_body_tide_degree_order_variable_complex_k",
          py::overload_cast<const std::string, const std::map<int, std::vector<std::complex<double> > > >(
              &tss::degreeOrderVariableLoveNumberGravityFieldVariationSettingsPy),
          py::arg("tide_raising_body"),
          py::arg("love_number_per_degree_and_order"),
          get_docstring("solid_body_tide_degree_order_variable_complex_k").c_str() );

    m.def("single_period_periodic",
          py::overload_cast< const Eigen::MatrixXd&, const Eigen::MatrixXd&, const Eigen::MatrixXd&,
                  const Eigen::MatrixXd&, const double, const double, const int, const int>(
              &tss::periodicGravityFieldVariationsSettingsSingleFrequency),
          py::arg("cosine_sh_amplitude_cosine_time"),
          py::arg("cosine_sh_amplitude_sine_time"),
          py::arg("sine_sh_amplitude_cosine_time"),
          py::arg("sine_sh_amplitude_sine_time"),
          py::arg("frequency"),
          py::arg("reference_epoch"),
          py::arg("minimum_degree") = 2,
          py::arg("minimum_order") = 0,
          get_docstring("single_period_periodic").c_str() );

    m.def("periodic",
          py::overload_cast< const std::vector< Eigen::MatrixXd >&, const std::vector< Eigen::MatrixXd >&,
                            const std::vector< Eigen::MatrixXd >&, const std::vector< Eigen::MatrixXd >&,
                             const std::vector< double >&, const double, const int, const int>(
              &tss::periodicGravityFieldVariationsSettings),
          py::arg("cosine_sh_amplitudes_cosine_time"),
          py::arg("cosine_sh_amplitudes_sine_time"),
          py::arg("sine_sh_amplitudes_cosine_time"),
          py::arg("sine_sh_amplitudes_sine_time"),
          py::arg("frequencies"),
          py::arg("reference_epoch"),
          py::arg("minimum_degree") = 2,
          py::arg("minimum_order") = 0,
          get_docstring("single_period_periodic").c_str() );

    m.def("single_period_polynomial",
          &tss::polynomialGravityFieldVariationsSettingsSinglePower,
          py::arg("cosine_amplitude_per_power"),
          py::arg("sine_amplitude_per_power"),
          py::arg("polynomial_power"),
          py::arg("reference_epoch"),
          py::arg("minimum_degree") = 2,
          py::arg("minimum_order") = 0,
          get_docstring("single_period_polynomial").c_str() );

    m.def("polynomial",
          &tss::polynomialGravityFieldVariationsSettings,
          py::arg("cosine_amplitude_per_power"),
          py::arg("sine_amplitude_per_power"),
          py::arg("reference_epoch"),
          py::arg("minimum_degree") = 2,
          py::arg("minimum_order") = 0,
          get_docstring("polynomial").c_str() );


    m.def("tabulated",
          &tss::tabulatedGravityFieldVariationSettings,
          py::arg("cosine_variations_table"),
          py::arg("sine_variations_table"),
          py::arg("minimum_degree"),
          py::arg("minimum_order"),
          py::arg("interpolation_settings"));
}

}// namespace gravity_field_variation
}// namespace environment_setup
}// namespace numerical_simulation
}// namespace tudatpy
