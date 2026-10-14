/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#include "expose_station_history.h"
#include <pybind11/stl.h>
#include <sstream>
#include "tudat/io/readSinexStationHistory.h"
namespace py = pybind11;
namespace tudatpy
{
namespace data_input
{
namespace environment_data
{
namespace station_history
{
namespace tio = tudat::input_output;
void expose_station_history( py::module& m )
{
    py::enum_< tio::SinexEpochMeaning >( m, "SinexEpochMeaning" )
            .value( "specified", tio::SinexEpochMeaning::specified )
            .value( "unspecified", tio::SinexEpochMeaning::unspecified )
            .value( "header_start_default", tio::SinexEpochMeaning::headerStartDefault )
            .value( "header_end_default", tio::SinexEpochMeaning::headerEndDefault )
            .value( "open_start", tio::SinexEpochMeaning::openStart )
            .value( "open_end", tio::SinexEpochMeaning::openEnd );
    py::enum_< tio::SinexNormalization >( m, "SinexNormalization" )
            .value( "si", tio::SinexNormalization::si )
            .value( "unresolved_extension", tio::SinexNormalization::unresolvedExtension )
            .value( "unsupported_unit", tio::SinexNormalization::unsupportedUnit );
    py::enum_< tio::SinexEccentricityBasis >( m, "SinexEccentricityBasis" )
            .value( "une", tio::SinexEccentricityBasis::une )
            .value( "xyz", tio::SinexEccentricityBasis::xyz );
    py::enum_< tio::SinexEquipmentKind >( m, "SinexEquipmentKind" )
            .value( "receiver", tio::SinexEquipmentKind::receiver )
            .value( "antenna", tio::SinexEquipmentKind::antenna );
    py::class_< tio::SinexSiteKey >( m, "SinexSiteKey" )
            .def( py::init( []( std::string site, std::string point ) { return tio::SinexSiteKey{ site, point }; } ),
                  py::arg( "site_code" ),
                  py::arg( "point_code" ) )
            .def_readonly( "site_code", &tio::SinexSiteKey::siteCode )
            .def_readonly( "point_code", &tio::SinexSiteKey::pointCode );
    py::class_< tio::SinexSolutionKey >( m, "SinexSolutionKey" )
            .def( py::init( []( std::string site, std::string point, std::string solution ) {
                      return tio::SinexSolutionKey{ { site, point }, solution };
                  } ),
                  py::arg( "site_code" ),
                  py::arg( "point_code" ),
                  py::arg( "solution_id" ) )
            .def_readonly( "site", &tio::SinexSolutionKey::site )
            .def_readonly( "solution_id", &tio::SinexSolutionKey::solutionId );
    py::class_< tio::SinexSourceRecord >( m, "SinexSourceRecord" )
            .def_readonly( "source", &tio::SinexSourceRecord::source )
            .def_readonly( "line", &tio::SinexSourceRecord::line )
            .def_readonly( "text", &tio::SinexSourceRecord::text );
    py::class_< tio::SinexCivilUtc >( m, "SinexCivilUtc" )
            .def_readonly( "year", &tio::SinexCivilUtc::year )
            .def_readonly( "day_of_year", &tio::SinexCivilUtc::dayOfYear )
            .def_readonly( "seconds_of_day", &tio::SinexCivilUtc::secondsOfDay );
    py::class_< tio::SinexEpoch >( m, "SinexEpoch" )
            .def_readonly( "original", &tio::SinexEpoch::original )
            .def_readonly( "meaning", &tio::SinexEpoch::meaning )
            .def_readonly( "civil_utc", &tio::SinexEpoch::civilUtc );
    py::class_< tio::SinexSiteRecord >( m, "SinexSiteRecord" )
            .def_readonly( "key", &tio::SinexSiteRecord::key )
            .def_readonly( "domes_id", &tio::SinexSiteRecord::domesId )
            .def_readonly( "observation_technique", &tio::SinexSiteRecord::observationTechnique )
            .def_readonly( "description", &tio::SinexSiteRecord::description )
            .def_readonly( "approximate_longitude_radians", &tio::SinexSiteRecord::approximateLongitudeRadians )
            .def_readonly( "approximate_latitude_radians", &tio::SinexSiteRecord::approximateLatitudeRadians )
            .def_readonly( "approximate_height_metres", &tio::SinexSiteRecord::approximateHeightMetres )
            .def_readonly( "provenance", &tio::SinexSiteRecord::provenance );
    py::class_< tio::SinexParameterRecord >( m, "SinexParameterRecord" )
            .def_readonly( "index", &tio::SinexParameterRecord::index )
            .def_readonly( "type", &tio::SinexParameterRecord::type )
            .def_readonly( "key", &tio::SinexParameterRecord::key )
            .def_readonly( "reference_epoch", &tio::SinexParameterRecord::referenceEpoch )
            .def_readonly( "original_unit", &tio::SinexParameterRecord::originalUnit )
            .def_readonly( "constraint_code", &tio::SinexParameterRecord::constraintCode )
            .def_readonly( "original_value", &tio::SinexParameterRecord::originalValue )
            .def_readonly( "original_standard_deviation", &tio::SinexParameterRecord::originalStandardDeviation )
            .def_readonly( "normalization", &tio::SinexParameterRecord::normalization )
            .def_readonly( "value_si", &tio::SinexParameterRecord::valueSi )
            .def_readonly( "standard_deviation_si", &tio::SinexParameterRecord::standardDeviationSi )
            .def_readonly( "si_unit", &tio::SinexParameterRecord::siUnit )
            .def_readonly( "provenance", &tio::SinexParameterRecord::provenance );
    py::class_< tio::SinexObservationSpan >( m, "SinexObservationSpan" )
            .def_readonly( "key", &tio::SinexObservationSpan::key )
            .def_readonly( "observation_technique", &tio::SinexObservationSpan::observationTechnique )
            .def_readonly( "start", &tio::SinexObservationSpan::start )
            .def_readonly( "end", &tio::SinexObservationSpan::end )
            .def_readonly( "mean", &tio::SinexObservationSpan::mean )
            .def_readonly( "provenance", &tio::SinexObservationSpan::provenance );
    py::class_< tio::SinexDiscontinuityRecord >( m, "SinexDiscontinuityRecord" )
            .def_readonly( "key", &tio::SinexDiscontinuityRecord::key )
            .def_readonly( "observation_technique", &tio::SinexDiscontinuityRecord::observationTechnique )
            .def_readonly( "kind", &tio::SinexDiscontinuityRecord::kind )
            .def_readonly( "start", &tio::SinexDiscontinuityRecord::start )
            .def_readonly( "end", &tio::SinexDiscontinuityRecord::end )
            .def_readonly( "description", &tio::SinexDiscontinuityRecord::description )
            .def_readonly( "provenance", &tio::SinexDiscontinuityRecord::provenance );
    py::class_< tio::SinexEccentricityRecord >( m, "SinexEccentricityRecord" )
            .def_readonly( "key", &tio::SinexEccentricityRecord::key )
            .def_readonly( "observation_technique", &tio::SinexEccentricityRecord::observationTechnique )
            .def_readonly( "start", &tio::SinexEccentricityRecord::start )
            .def_readonly( "end", &tio::SinexEccentricityRecord::end )
            .def_readonly( "basis", &tio::SinexEccentricityRecord::basis )
            .def_readonly( "components_metres", &tio::SinexEccentricityRecord::componentsMetres )
            .def_readonly( "provenance", &tio::SinexEccentricityRecord::provenance );
    py::class_< tio::SinexEquipmentRecord >( m, "SinexEquipmentRecord" )
            .def_readonly( "kind", &tio::SinexEquipmentRecord::kind )
            .def_readonly( "key", &tio::SinexEquipmentRecord::key )
            .def_readonly( "observation_technique", &tio::SinexEquipmentRecord::observationTechnique )
            .def_readonly( "start", &tio::SinexEquipmentRecord::start )
            .def_readonly( "end", &tio::SinexEquipmentRecord::end )
            .def_readonly( "equipment_type", &tio::SinexEquipmentRecord::equipmentType )
            .def_readonly( "serial_number", &tio::SinexEquipmentRecord::serialNumber )
            .def_readonly( "firmware_version", &tio::SinexEquipmentRecord::firmwareVersion )
            .def_readonly( "provenance", &tio::SinexEquipmentRecord::provenance );
    py::class_< tio::SinexOpaqueBlock >( m, "SinexOpaqueBlock" )
            .def_readonly( "name", &tio::SinexOpaqueBlock::name )
            .def_readonly( "records", &tio::SinexOpaqueBlock::records );
    py::class_< tio::SinexCartesianSolution >( m, "SinexCartesianSolution" )
            .def_readonly( "key", &tio::SinexCartesianSolution::key )
            .def_readonly( "site", &tio::SinexCartesianSolution::site )
            .def_readonly( "position_reference_epoch", &tio::SinexCartesianSolution::positionReferenceEpoch )
            .def_readonly( "position_metres", &tio::SinexCartesianSolution::positionMetres )
            .def_readonly( "position_parameter_indices", &tio::SinexCartesianSolution::positionParameterIndices )
            .def_readonly( "position_standard_deviations_metres", &tio::SinexCartesianSolution::positionStandardDeviationsMetres )
            .def_readonly( "velocity_metres_per_second", &tio::SinexCartesianSolution::velocityMetresPerSecond )
            .def_readonly( "velocity_parameter_indices", &tio::SinexCartesianSolution::velocityParameterIndices )
            .def_readonly( "velocity_standard_deviations_metres_per_second",
                           &tio::SinexCartesianSolution::velocityStandardDeviationsMetresPerSecond )
            .def_readonly( "velocity_reference_epochs", &tio::SinexCartesianSolution::velocityReferenceEpochs );
    py::class_< tio::SinexStationHistory >( m, "SinexStationHistory" )
            .def_readonly( "source", &tio::SinexStationHistory::source )
            .def_readonly( "original_header", &tio::SinexStationHistory::originalHeader )
            .def_readonly( "version", &tio::SinexStationHistory::version )
            .def_readonly( "agency", &tio::SinexStationHistory::agency )
            .def_readonly( "creation_epoch", &tio::SinexStationHistory::creationEpoch )
            .def_readonly( "header_start", &tio::SinexStationHistory::headerStart )
            .def_readonly( "header_end", &tio::SinexStationHistory::headerEnd )
            .def_readonly( "observation_technique", &tio::SinexStationHistory::observationTechnique )
            .def_readonly( "declared_parameter_count", &tio::SinexStationHistory::declaredParameterCount )
            .def_readwrite( "frame", &tio::SinexStationHistory::frame )
            .def_readwrite( "release", &tio::SinexStationHistory::release )
            .def_readwrite( "origin", &tio::SinexStationHistory::origin )
            .def_readwrite( "tide_system", &tio::SinexStationHistory::tideSystem )
            .def_readonly( "comments", &tio::SinexStationHistory::comments )
            .def_readonly( "parameters", &tio::SinexStationHistory::parameters )
            .def_readonly( "observation_spans", &tio::SinexStationHistory::observationSpans )
            .def_readonly( "discontinuities", &tio::SinexStationHistory::discontinuities )
            .def_readonly( "eccentricities", &tio::SinexStationHistory::eccentricities )
            .def_readonly( "equipment_histories", &tio::SinexStationHistory::equipmentHistories )
            .def_readonly( "opaque_blocks", &tio::SinexStationHistory::opaqueBlocks )
            .def_property_readonly( "sites", []( const tio::SinexStationHistory& history ) {
                // A Python dict would silently erase duplicate site/point records.
                std::vector< tio::SinexSiteRecord > records;
                for( const auto& row : history.sites ) records.push_back( row.second );
                return records;
            } );
    m.def( "read_sinex_station_history",
           py::overload_cast< const std::string& >( &tio::readSinexStationHistory ),
           py::arg( "file_name" ),
           "Read strict fixed-width SINEX station-product histories without constructing a model." );
    m.def(
            "read_sinex_station_history_text",
            []( const std::string& text, const std::string& source ) {
                std::istringstream stream( text );
                return tio::readSinexStationHistory( stream, source );
            },
            py::arg( "text" ),
            py::arg( "source" ) = "<text>",
            "Read SINEX text with the same validation and original-field preservation as the file reader." );
    m.def( "validate_sinex_station_history",
           &tio::validateSinexStationHistory,
           py::arg( "history" ),
           "Validate neutral histories without constructing or propagating an environment model." );
    m.def( "get_cartesian_solution_keys", &tio::getSinexCartesianSolutionKeys, py::arg( "history" ) );
    m.def( "get_cartesian_solution",
           &tio::getSinexCartesianSolution,
           py::arg( "history" ),
           py::arg( "key" ),
           py::arg( "require_complete_velocity" ) = false,
           "Select a complete Cartesian position explicitly; missing velocity components remain None." );
}
}  // namespace station_history
}  // namespace environment_data
}  // namespace data_input
}  // namespace tudatpy
