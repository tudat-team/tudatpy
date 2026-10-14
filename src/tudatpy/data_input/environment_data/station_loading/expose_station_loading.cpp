#include "expose_station_loading.h"
#include <pybind11/stl.h>
#include "tudat/io/readStationLoadingData.h"

namespace tudatpy
{
namespace data_input
{
namespace environment_data
{
namespace station_loading
{
void expose_station_loading( pybind11::module& module )
{
    namespace py = pybind11;
    namespace io = tudat::input_output;
    py::enum_< io::StationLoadingTimeScale >( module, "StationLoadingTimeScale" )
            .value( "tai", io::StationLoadingTimeScale::tai )
            .value( "tt", io::StationLoadingTimeScale::tt )
            .value( "utc", io::StationLoadingTimeScale::utc )
            .value( "tdb", io::StationLoadingTimeScale::tdb );
    py::enum_< io::StationLoadingBasis >( module, "StationLoadingBasis" )
            .value( "unspecified", io::StationLoadingBasis::unspecified )
            .value( "geocentric_radial_east_north", io::StationLoadingBasis::geocentricRadialEastNorth )
            .value( "cartesian", io::StationLoadingBasis::cartesian );
    py::enum_< io::StationLoadingComponent >( module, "StationLoadingComponent" )
            .value( "atmosphere", io::StationLoadingComponent::atmosphere )
            .value( "hydrology", io::StationLoadingComponent::hydrology )
            .value( "ocean", io::StationLoadingComponent::ocean );
    py::enum_< io::StationLoadingSignalPresence >( module, "StationLoadingSignalPresence" )
            .value( "unknown", io::StationLoadingSignalPresence::unknown )
            .value( "absent", io::StationLoadingSignalPresence::absent )
            .value( "present", io::StationLoadingSignalPresence::present );
    py::enum_< io::StationLoadingOrigin >( module, "StationLoadingOrigin" )
            .value( "unspecified", io::StationLoadingOrigin::unspecified )
            .value( "whole_earth_centre_of_mass", io::StationLoadingOrigin::wholeEarthCentreOfMass )
            .value( "solid_earth_centre_of_mass", io::StationLoadingOrigin::solidEarthCentreOfMass )
            .value( "centre_of_figure", io::StationLoadingOrigin::centreOfFigure )
            .value( "custom", io::StationLoadingOrigin::custom );
    py::class_< io::StationLoadingEpoch >(
            module, "StationLoadingEpoch", "Integer MJD and fractional seconds in the declared source clock." )
            .def( py::init<>( ) )
            .def_readwrite( "modified_julian_day", &io::StationLoadingEpoch::modifiedJulianDay )
            .def_readwrite( "seconds_of_day", &io::StationLoadingEpoch::secondsOfDay );
    py::class_< io::StationLoadingMetadata >( module, "StationLoadingMetadata" )
            .def( py::init<>( ) )
            .def_readwrite( "source", &io::StationLoadingMetadata::source )
            .def_readwrite( "source_headers", &io::StationLoadingMetadata::sourceHeaders )
            .def_readwrite( "product", &io::StationLoadingMetadata::product )
            .def_readwrite( "release", &io::StationLoadingMetadata::release )
            .def_readwrite( "source_profile", &io::StationLoadingMetadata::sourceProfile )
            .def_readwrite( "reference_frame", &io::StationLoadingMetadata::referenceFrame )
            .def_readwrite( "reference_origin", &io::StationLoadingMetadata::referenceOrigin )
            .def_readwrite( "origin_definition", &io::StationLoadingMetadata::originDefinition )
            .def_readwrite( "ocean_response", &io::StationLoadingMetadata::oceanResponse )
            .def_readwrite( "components", &io::StationLoadingMetadata::components )
            .def_readwrite( "diurnal", &io::StationLoadingMetadata::diurnal )
            .def_readwrite( "semidiurnal", &io::StationLoadingMetadata::semidiurnal )
            .def_readwrite( "annual", &io::StationLoadingMetadata::annual )
            .def_readwrite( "semiannual", &io::StationLoadingMetadata::semiannual )
            .def_readwrite( "mean", &io::StationLoadingMetadata::mean )
            .def_readwrite( "trend", &io::StationLoadingMetadata::trend )
            .def_readwrite( "removed_harmonics", &io::StationLoadingMetadata::removedHarmonics )
            .def_readwrite( "original_displacement_units", &io::StationLoadingMetadata::originalDisplacementUnits );
    py::class_< io::StationLoadingSample >( module, "StationLoadingSample", "Source-basis displacement in metres." )
            .def( py::init<>( ) )
            .def_readwrite( "displacement", &io::StationLoadingSample::displacement )
            .def_readwrite( "source_line", &io::StationLoadingSample::sourceLine )
            .def_readwrite( "source_record", &io::StationLoadingSample::sourceRecord )
            .def_readwrite( "informational_epoch", &io::StationLoadingSample::informationalEpoch );
    py::class_< io::StationLoadingSite >(
            module, "StationLoadingSite", "Internal provider label and Cartesian metre location; no implicit name match." )
            .def( py::init<>( ) )
            .def_readwrite( "identifier", &io::StationLoadingSite::identifier )
            .def_readwrite( "cartesian_position", &io::StationLoadingSite::cartesianPosition )
            .def_readwrite( "validity_radius_metres", &io::StationLoadingSite::validityRadiusMetres )
            .def_readwrite( "source_line", &io::StationLoadingSite::sourceLine )
            .def_readwrite( "source_record", &io::StationLoadingSite::sourceRecord )
            .def_readwrite( "informational_coordinates", &io::StationLoadingSite::informationalCoordinates )
            .def_readwrite( "samples", &io::StationLoadingSite::samples );
    py::class_< io::StationLoadingData >(
            module, "StationLoadingData", "Preloaded source series; sample map indexes into epochs are zero based." )
            .def( py::init<>( ) )
            .def_readwrite( "metadata", &io::StationLoadingData::metadata )
            .def_readwrite( "time_scale", &io::StationLoadingData::timeScale )
            .def_readwrite( "basis", &io::StationLoadingData::basis )
            .def_readwrite( "epochs", &io::StationLoadingData::epochs )
            .def_readwrite( "declared_begin_epoch", &io::StationLoadingData::declaredBeginEpoch )
            .def_readwrite( "declared_end_epoch", &io::StationLoadingData::declaredEndEpoch )
            .def_readwrite( "nominal_sample_interval_seconds", &io::StationLoadingData::nominalSampleIntervalSeconds )
            .def_readwrite( "sites", &io::StationLoadingData::sites )
            .def_readwrite( "informational_begin_epoch", &io::StationLoadingData::informationalBeginEpoch )
            .def_readwrite( "informational_end_epoch", &io::StationLoadingData::informationalEndEpoch )
            .def( "validate", &io::validateStationLoadingData );
    module.def( "read_ephedisp_station_loading_data",
                py::overload_cast< const std::string&, const io::StationLoadingMetadata& >( &io::readEphedispStationLoadingData ),
                py::arg( "file_path" ),
                py::arg( "metadata" ) = io::StationLoadingMetadata( ) );
    module.def( "imls_merra2_atmospheric_loading_metadata",
                &io::imlsMerra2AtmosphericLoadingMetadata,
                py::arg( "release" ),
                "Explicit IMLS MERRA2 atmospheric source profile, based on the pinned provider description; release is caller supplied." );
}
}  // namespace station_loading
}  // namespace environment_data
}  // namespace data_input
}  // namespace tudatpy
