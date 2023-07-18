/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_io.h"

#include "tudatpy/docstrings.h"

#include "tudat/io/readHistoryFromFile.h"
#include "tudat/io/missileDatcomData.h"
#include "tudat/io/solarActivityData.h"
#include "tudat/io/readOdfFile.h"
#include "tudat/io/readTabulatedWeatherData.h"

#include <pybind11/stl.h>

namespace py = pybind11;
namespace tio = tudat::input_output;

namespace tudatpy {

namespace io {

void expose_io(py::module &m) {

      m.def("get_resource_path",
            &tudat::paths::get_resource_path,
            get_docstring("get_resource_path").c_str() );
      m.def("get_ephemeris_path",
            &tudat::paths::getEphemerisDataFilesPath,
            get_docstring("get_ephemeris_path").c_str() );
      m.def("get_earth_orientation_path",
            &tudat::paths::getEarthOrientationDataFilesPath,
            get_docstring("get_earth_orientation_path").c_str() );
      m.def("get_quadrature_path",
            &tudat::paths::getQuadratureDataPath,
            get_docstring("get_quadrature_path").c_str() );
      m.def("get_spice_kernel_path",
            &tudat::paths::getSpiceKernelPath,
            get_docstring("get_spice_kernel_path").c_str() );
      m.def("get_atmosphere_tables_path",
            &tudat::paths::getAtmosphereTablesPath,
            get_docstring("get_atmosphere_tables_path").c_str() );
      m.def("get_gravity_models_path",
            &tudat::paths::getGravityModelsPath,
            get_docstring("get_gravity_models_path").c_str() );
      m.def("get_space_weather_path",
            &tudat::paths::getSpaceWeatherDataPath,
            get_docstring("get_space_weather_path").c_str()  );

      m.def("read_vector_history_from_file",
            &tudat::input_output::readVectorHistoryFromFile< double, double >,
            py::arg("vector_size"),
            py::arg("file_name"),
            get_docstring("read_vector_history_from_file").c_str()
      );

      m.def("read_matrix_history_from_file",
            &tudat::input_output::readMatrixHistoryFromFile< double, double >,
            py::arg("matrix_rows"),
            py::arg("matrix_columns"),
            py::arg("file_name"),
            get_docstring("read_matrix_history_from_file").c_str()
      );

      py::class_<tudat::input_output::MissileDatcomData,
            std::shared_ptr<tudat::input_output::MissileDatcomData>>(m, "missile_DATCOM_data",
                                                                     get_docstring("missile_DATCOM_data").c_str())
            .def(py::init<
                  const std::string &>(),
                  py::arg("file_name_and_path"),
                  get_docstring("missile_DATCOM_data.ctor").c_str())
            .def("get_static_coefficient",
                  &tudat::input_output::MissileDatcomData::getStaticCoefficient,
                  py::arg("mach_index"),
                  py::arg("angle_of_attack_index"),
                  py::arg("coefficient_index"),
                  get_docstring("missile_DATCOM_data.get_static_coefficient").c_str())
            .def("get_dynamic_coefficient",
                  &tudat::input_output::MissileDatcomData::getDynamicCoefficient,
                  py::arg("mach_index"),
                  py::arg("angle_of_attack_index"),
                  py::arg("coefficient_index"),
                  get_docstring("missile_DATCOM_data.get_dynamic_coefficient").c_str())
            .def("get_angle_of_attacks",
                  &tudat::input_output::MissileDatcomData::getAngleOfAttacks,
                  get_docstring("missile_DATCOM_data.get_angle_of_attacks").c_str())
            .def("get_mach_numbers",
                  &tudat::input_output::MissileDatcomData::getMachNumbers,
                  get_docstring("missile_DATCOM_data.get_mach_numbers").c_str())
            .def("get_Reynolds_numbers",
                  &tudat::input_output::MissileDatcomData::getReynoldsNumbers,
                  get_docstring("missile_DATCOM_data.get_Reynolds_numbers").c_str())
            .def("write_all_coefficients_to_files",
                  &tudat::input_output::MissileDatcomData::writeAllCoefficientsToFiles,
                  py::arg("file_name_base"),
                  py::arg("base_precision") = 15,
                  py::arg("exponent_width") = 2,
                  get_docstring("missile_DATCOM_data.write_all_coefficients_to_files").c_str() )
            .def("write_force_and_moment_coefficients_to_files",
                  &tudat::input_output::MissileDatcomData::writeForceAndMomentCoefficientsToFiles,
                  py::arg("file_name_base"),
                  py::arg("base_precision") = 15,
                  py::arg("exponent_width") = 2,
                  get_docstring("missile_DATCOM_data.write_force_and_moment_coefficients_to_files").c_str() );

      py::enum_<tudat::input_output::MissileDatcomData::DynamicCoefficientNames>(m, "DynamicCoefficientNames",
                                                                                 get_docstring("DynamicCoefficientNames").c_str())
            .value("cnq", tudat::input_output::MissileDatcomData::DynamicCoefficientNames::cnq, get_docstring("DynamicCoefficientNames.cnq").c_str())
            .value("cmq", tudat::input_output::MissileDatcomData::DynamicCoefficientNames::cmq, get_docstring("DynamicCoefficientNames.cmq").c_str())
            .value("caq", tudat::input_output::MissileDatcomData::DynamicCoefficientNames::caq, get_docstring("DynamicCoefficientNames.caq").c_str())
            .value("cyq", tudat::input_output::MissileDatcomData::DynamicCoefficientNames::cyq, get_docstring("DynamicCoefficientNames.cyq").c_str())
            .value("clnq", tudat::input_output::MissileDatcomData::DynamicCoefficientNames::clnq, get_docstring("DynamicCoefficientNames.clnq").c_str())
            .value("cllq", tudat::input_output::MissileDatcomData::DynamicCoefficientNames::cllq, get_docstring("DynamicCoefficientNames.cllq").c_str())
            .value("cnr", tudat::input_output::MissileDatcomData::DynamicCoefficientNames::cnr, get_docstring("DynamicCoefficientNames.cnr").c_str())
            .value("cmr", tudat::input_output::MissileDatcomData::DynamicCoefficientNames::cmr, get_docstring("DynamicCoefficientNames.cmr").c_str())
            .value("car", tudat::input_output::MissileDatcomData::DynamicCoefficientNames::car, get_docstring("DynamicCoefficientNames.car").c_str())
            .value("cyr", tudat::input_output::MissileDatcomData::DynamicCoefficientNames::cyr, get_docstring("DynamicCoefficientNames.cyr").c_str())
            .value("clnr", tudat::input_output::MissileDatcomData::DynamicCoefficientNames::clnr, get_docstring("DynamicCoefficientNames.clnr").c_str())
            .value("cllr", tudat::input_output::MissileDatcomData::DynamicCoefficientNames::cllr, get_docstring("DynamicCoefficientNames.cllr").c_str())
            .value("cnp", tudat::input_output::MissileDatcomData::DynamicCoefficientNames::cnp, get_docstring("DynamicCoefficientNames.cnp").c_str())
            .value("cmp", tudat::input_output::MissileDatcomData::DynamicCoefficientNames::cmp, get_docstring("DynamicCoefficientNames.cmp").c_str())
            .value("cap", tudat::input_output::MissileDatcomData::DynamicCoefficientNames::cap, get_docstring("DynamicCoefficientNames.cap").c_str())
            .value("cyp", tudat::input_output::MissileDatcomData::DynamicCoefficientNames::cyp, get_docstring("DynamicCoefficientNames.cyp").c_str())
            .value("clnp", tudat::input_output::MissileDatcomData::DynamicCoefficientNames::clnp, get_docstring("DynamicCoefficientNames.clnp").c_str())
            .value("cllp", tudat::input_output::MissileDatcomData::DynamicCoefficientNames::cllp, get_docstring("DynamicCoefficientNames.cllp").c_str())
            .value("cnad", tudat::input_output::MissileDatcomData::DynamicCoefficientNames::cnad, get_docstring("DynamicCoefficientNames.cnad").c_str())
            .value("cmad", tudat::input_output::MissileDatcomData::DynamicCoefficientNames::cmad, get_docstring("DynamicCoefficientNames.cmad").c_str())
            .export_values();

      py::enum_<tudat::input_output::MissileDatcomData::StaticCoefficientNames>(m, "StaticCoefficientNames",
                                                                                get_docstring("StaticCoefficientNames").c_str())
            .value("cn", tudat::input_output::MissileDatcomData::StaticCoefficientNames::cn, get_docstring("StaticCoefficientNames.cn").c_str())
            .value("cm", tudat::input_output::MissileDatcomData::StaticCoefficientNames::cm, get_docstring("StaticCoefficientNames.cm").c_str())
            .value("ca", tudat::input_output::MissileDatcomData::StaticCoefficientNames::ca, get_docstring("StaticCoefficientNames.ca").c_str())
            .value("cy", tudat::input_output::MissileDatcomData::StaticCoefficientNames::cy, get_docstring("StaticCoefficientNames.cy").c_str())
            .value("cln", tudat::input_output::MissileDatcomData::StaticCoefficientNames::cln, get_docstring("StaticCoefficientNames.cln").c_str())
            .value("cll", tudat::input_output::MissileDatcomData::StaticCoefficientNames::cll, get_docstring("StaticCoefficientNames.cll").c_str())
            .value("cna", tudat::input_output::MissileDatcomData::StaticCoefficientNames::cna, get_docstring("StaticCoefficientNames.cna").c_str())
            .value("cma", tudat::input_output::MissileDatcomData::StaticCoefficientNames::cma, get_docstring("StaticCoefficientNames.cma").c_str())
            .value("cyb", tudat::input_output::MissileDatcomData::StaticCoefficientNames::cyb, get_docstring("StaticCoefficientNames.cyb").c_str())
            .value("cnb", tudat::input_output::MissileDatcomData::StaticCoefficientNames::cnb, get_docstring("StaticCoefficientNames.cnb").c_str())
            .value("clb", tudat::input_output::MissileDatcomData::StaticCoefficientNames::clb, get_docstring("StaticCoefficientNames.clb").c_str())
            .export_values();

      py::class_<
            tio::solar_activity::SolarActivityData,
            std::shared_ptr< tio::solar_activity::SolarActivityData > >(
                    m, "SolarActivityData", get_docstring("SolarActivityData").c_str() );

      py::class_<
            tio::OdfRawFileContents,
            std::shared_ptr< tio::OdfRawFileContents > >(
                    m, "OdfRawFileContents", get_docstring("OdfRawFileContents").c_str( ) )
            .def("write_to_text_file",
                 &tio::OdfRawFileContents::writeOdfToTextFile,
                 py::arg("output_file"),
                 get_docstring("OdfRawFileContents.write_to_text_file").c_str() );

      m.def("read_odf_file",
          &tio::readOdfFile,
          py::arg( "file_name" ),
          get_docstring("read_odf_file").c_str() );

      m.def("set_dsn_weather_data_in_ground_stations",
          py::overload_cast<
              tudat::simulation_setup::SystemOfBodies&,
              const std::vector< std::string >&,
              std::shared_ptr< tudat::interpolators::InterpolatorSettings >,
              const std::map< int, std::vector< std::string > >&,
              const std::string& >( &tio::setDsnWeatherDataInGroundStations ),
          py::arg( "bodies" ),
          py::arg( "weather_file_names" ),
          py::arg( "interpolator_settings" ) = tudat::interpolators::linearInterpolation( ),
          py::arg( "ground_stations_per_complex" ) = tudat::simulation_setup::getDefaultDsnStationNamesPerComplex( ),
          py::arg( "body_with_ground_stations_name" ) = "Earth",
          get_docstring("set_dsn_weather_data_in_ground_stations").c_str() );
};

}// namespace io
}// namespace tudatpy