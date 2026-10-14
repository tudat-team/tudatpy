#include "expose_ocean_loading.h"

#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include "tudat/io/readOceanLoadingData.h"
#include "tudat/io/readOceanPoleTideData.h"

namespace tudatpy
{
namespace data_input
{
namespace environment_data
{
namespace ocean_loading
{

void expose_ocean_loading( pybind11::module& module )
{
    namespace py = pybind11;
    namespace io = tudat::input_output;
    py::enum_< io::OceanLoadingReferenceOrigin >( module, "OceanLoadingReferenceOrigin" )
            .value( "unspecified", io::OceanLoadingReferenceOrigin::unspecified )
            .value( "solid_earth_centre", io::OceanLoadingReferenceOrigin::solidEarthCentre )
            .value( "whole_earth_centre_of_mass", io::OceanLoadingReferenceOrigin::wholeEarthCentreOfMass );
    py::class_< io::OceanLoadingStation >(
            module,
            "OceanLoadingStation",
            "BLQ source coefficients: rows up/radial, west, south; amplitudes in metres and Greenwich phase lags in degrees." )
            .def( py::init<>( ) )
            .def_readwrite( "station_identifier", &io::OceanLoadingStation::stationIdentifier )
            .def_readwrite( "amplitudes_metres", &io::OceanLoadingStation::amplitudesMetres )
            .def_readwrite( "phase_lags_degrees", &io::OceanLoadingStation::phaseLagsDegrees )
            .def_readwrite( "source", &io::OceanLoadingStation::source )
            .def_readwrite( "source_headers", &io::OceanLoadingStation::sourceHeaders )
            .def_readwrite( "ocean_tide_model", &io::OceanLoadingStation::oceanTideModel )
            .def_readwrite( "greens_function", &io::OceanLoadingStation::greensFunction )
            .def_readwrite( "reference_origin", &io::OceanLoadingStation::referenceOrigin )
            .def_readwrite( "provider_longitude_latitude_height", &io::OceanLoadingStation::providerLongitudeLatitudeHeight )
            .def_readwrite( "has_provider_coordinates", &io::OceanLoadingStation::hasProviderCoordinates );
    py::class_< io::OceanLoadingData >( module, "OceanLoadingData" )
            .def( py::init<>( ) )
            .def_readwrite( "source", &io::OceanLoadingData::source )
            .def_readwrite( "source_headers", &io::OceanLoadingData::sourceHeaders )
            .def_readwrite( "stations", &io::OceanLoadingData::stations )
            .def( "get_station",
                  &io::OceanLoadingData::getStation,
                  py::arg( "station_identifier" ),
                  py::return_value_policy::reference_internal )
            .def( "validate", []( const io::OceanLoadingData& data ) { io::validateOceanLoadingData( data ); } );
    module.def( "read_ocean_loading_data",
                py::overload_cast< const std::string& >( &io::readOceanLoadingData ),
                py::arg( "file_path" ),
                "Read validated BLQ data, preserving provider units, headers, station identifiers and origin declarations." );
    py::class_< io::OceanPoleTideGrid >( module,
                                         "OceanPoleTideGrid",
                                         "Dimensionless ocean pole-load grid: longitude-fast nodes, radial/north/east rows and "
                                         "real/imaginary columns; no interpolation." )
            .def( py::init<>( ) )
            .def_readwrite( "longitude_count", &io::OceanPoleTideGrid::longitudeCount )
            .def_readwrite( "latitude_count", &io::OceanPoleTideGrid::latitudeCount )
            .def_readwrite( "first_longitude_degrees", &io::OceanPoleTideGrid::firstLongitudeDegrees )
            .def_readwrite( "first_latitude_degrees", &io::OceanPoleTideGrid::firstLatitudeDegrees )
            .def_readwrite( "longitude_step_degrees", &io::OceanPoleTideGrid::longitudeStepDegrees )
            .def_readwrite( "latitude_step_degrees", &io::OceanPoleTideGrid::latitudeStepDegrees )
            .def_readwrite( "coefficients", &io::OceanPoleTideGrid::coefficients )
            .def_readwrite( "source", &io::OceanPoleTideGrid::source )
            .def_readwrite( "source_headers", &io::OceanPoleTideGrid::sourceHeaders )
            .def_readwrite( "reference_origin", &io::OceanPoleTideGrid::referenceOrigin )
            .def( "validate", &io::validateOceanPoleTideGrid );
    module.def(
            "read_ocean_pole_tide_grid",
            py::overload_cast< const std::string& >( &io::readOceanPoleTideGrid ),
            py::arg( "file_path" ),
            "Read an IERS plain-text ocean pole-load grid; decompress .gz files before reading. No interpolation or origin conversion." );
    module.def( "constituent_names", &io::getOceanLoadingConstituentNames );
}

}  // namespace ocean_loading
}  // namespace environment_data
}  // namespace data_input
}  // namespace tudatpy
