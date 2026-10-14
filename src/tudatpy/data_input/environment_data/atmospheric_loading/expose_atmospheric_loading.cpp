#include "expose_atmospheric_loading.h"
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include "tudat/io/readAtmosphericTideData.h"
namespace tudatpy
{
namespace data_input
{
namespace environment_data
{
namespace atmospheric_loading
{
void expose_atmospheric_loading( pybind11::module& module )
{
    namespace py = pybind11;
    namespace io = tudat::input_output;
    py::enum_< io::AtmosphericLoadingReferenceOrigin >( module, "AtmosphericLoadingReferenceOrigin" )
            .value( "unspecified", io::AtmosphericLoadingReferenceOrigin::unspecified )
            .value( "solid_earth_centre", io::AtmosphericLoadingReferenceOrigin::solidEarthCentre )
            .value( "whole_earth_centre_of_mass", io::AtmosphericLoadingReferenceOrigin::wholeEarthCentreOfMass );
    py::enum_< io::AtmosphericTideLatitudeConvention >( module, "AtmosphericTideLatitudeConvention" )
            .value( "unspecified", io::AtmosphericTideLatitudeConvention::unspecified )
            .value( "geocentric", io::AtmosphericTideLatitudeConvention::geocentric )
            .value( "geodetic", io::AtmosphericTideLatitudeConvention::geodetic );
    py::enum_< io::AtmosphericTideDisplacementBasis >( module, "AtmosphericTideDisplacementBasis" )
            .value( "unspecified", io::AtmosphericTideDisplacementBasis::unspecified )
            .value( "geocentric_radial_north_east", io::AtmosphericTideDisplacementBasis::geocentricRadialNorthEast )
            .value( "geodetic_up_north_east", io::AtmosphericTideDisplacementBasis::geodeticUpNorthEast );
    py::enum_< io::AtmosphericTideSourceProfile >( module, "AtmosphericTideSourceProfile" )
            .value( "unspecified", io::AtmosphericTideSourceProfile::unspecified )
            .value( "iers2010_rp03_no_ib_farrell", io::AtmosphericTideSourceProfile::iers2010Rp03NoIbFarrell );
    py::class_< io::AtmosphericTideMetadata >( module, "AtmosphericTideMetadata" )
            .def( py::init<>( ) )
            .def_readwrite( "source", &io::AtmosphericTideMetadata::source )
            .def_readwrite( "source_headers", &io::AtmosphericTideMetadata::sourceHeaders )
            .def_readwrite( "model_name", &io::AtmosphericTideMetadata::modelName )
            .def_readwrite( "greens_function", &io::AtmosphericTideMetadata::greensFunction )
            .def_readwrite( "reference_origin", &io::AtmosphericTideMetadata::referenceOrigin )
            .def_readwrite( "latitude_convention", &io::AtmosphericTideMetadata::latitudeConvention )
            .def_readwrite( "displacement_basis", &io::AtmosphericTideMetadata::displacementBasis )
            .def_readwrite( "source_profile", &io::AtmosphericTideMetadata::sourceProfile )
            .def_readwrite( "expected_longitude_count", &io::AtmosphericTideMetadata::expectedLongitudeCount )
            .def_readwrite( "expected_latitude_count", &io::AtmosphericTideMetadata::expectedLatitudeCount )
            .def_readwrite( "original_coefficient_units", &io::AtmosphericTideMetadata::originalCoefficientUnits )
            .def_readwrite( "original_coefficient_order", &io::AtmosphericTideMetadata::originalCoefficientOrder );
    py::class_< io::AtmosphericTideGrid >( module, "AtmosphericTideGrid" )
            .def( py::init<>( ) )
            .def_readwrite( "metadata", &io::AtmosphericTideGrid::metadata )
            .def_readwrite( "longitudes_degrees", &io::AtmosphericTideGrid::longitudesDegrees )
            .def_readwrite( "latitudes_degrees", &io::AtmosphericTideGrid::latitudesDegrees )
            .def_readwrite( "coefficients_metres", &io::AtmosphericTideGrid::coefficientsMetres )
            .def( "validate", &io::validateAtmosphericTideGrid );
    py::class_< io::AtmosphericTideStation >( module, "AtmosphericTideStation" )
            .def( py::init<>( ) )
            .def_readwrite( "station_identifier", &io::AtmosphericTideStation::stationIdentifier )
            .def_readwrite( "provider_station_identifier_line", &io::AtmosphericTideStation::providerStationIdentifierLine )
            .def_readwrite( "coefficients_metres", &io::AtmosphericTideStation::coefficientsMetres )
            .def_readwrite( "metadata", &io::AtmosphericTideStation::metadata )
            .def_readwrite( "has_provider_coordinates", &io::AtmosphericTideStation::hasProviderCoordinates )
            .def_readwrite( "provider_longitude_latitude_degrees", &io::AtmosphericTideStation::providerLongitudeLatitudeDegrees );
    py::class_< io::AtmosphericTideStationData >( module, "AtmosphericTideStationData" )
            .def( py::init<>( ) )
            .def_readwrite( "metadata", &io::AtmosphericTideStationData::metadata )
            .def_readwrite( "stations", &io::AtmosphericTideStationData::stations )
            .def( "get_station",
                  &io::AtmosphericTideStationData::getStation,
                  py::arg( "station_identifier" ),
                  py::return_value_policy::reference_internal )
            .def( "validate", &io::validateAtmosphericTideStationData );
    module.def( "read_atmospheric_tide_grid",
                py::overload_cast< const std::string&, const io::AtmosphericTideMetadata& >( &io::readAtmosphericTideGrid ),
                py::arg( "file_path" ),
                py::arg( "metadata" ) = io::AtmosphericTideMetadata( ) );
    module.def( "read_atmospheric_tide_station_data",
                py::overload_cast< const std::string&, const io::AtmosphericTideMetadata& >( &io::readAtmosphericTideStationData ),
                py::arg( "file_path" ),
                py::arg( "metadata" ) = io::AtmosphericTideMetadata( ) );
}
}  // namespace atmospheric_loading
}  // namespace environment_data
}  // namespace data_input
}  // namespace tudatpy
