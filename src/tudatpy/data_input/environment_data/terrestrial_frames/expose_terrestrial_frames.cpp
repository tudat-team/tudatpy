#include "expose_terrestrial_frames.h"
#include <sstream>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include "tudat/io/readTerrestrialFrameTransformationData.h"
namespace tudatpy
{
namespace data_input
{
namespace environment_data
{
namespace terrestrial_frames
{
void expose_terrestrial_frames( pybind11::module& module )
{
    namespace py = pybind11;
    namespace io = tudat::input_output;
    py::enum_< io::HelmertRotationConvention >( module, "HelmertRotationConvention" )
            .value( "unspecified", io::HelmertRotationConvention::unspecified )
            .value( "position_vector", io::HelmertRotationConvention::positionVector )
            .value( "coordinate_frame", io::HelmertRotationConvention::coordinateFrame );
    py::class_< io::TerrestrialFrameIdentity >( module, "TerrestrialFrameIdentity" )
            .def( py::init<>( ) )
            .def_readwrite( "realization", &io::TerrestrialFrameIdentity::realization )
            .def_readwrite( "release", &io::TerrestrialFrameIdentity::release );
    py::class_< io::TerrestrialFrameTransformationSource >( module, "TerrestrialFrameTransformationSource" )
            .def( py::init<>( ) )
            .def_readwrite( "source", &io::TerrestrialFrameTransformationSource::source )
            .def_readwrite( "source_release", &io::TerrestrialFrameTransformationSource::sourceRelease )
            .def_readwrite( "source_sha256", &io::TerrestrialFrameTransformationSource::sourceSha256 );
    py::class_< io::TerrestrialFrameTransformationParameters >( module, "TerrestrialFrameTransformationParameters" )
            .def( py::init<>( ) )
            .def_readwrite( "source_frame", &io::TerrestrialFrameTransformationParameters::sourceFrame )
            .def_readwrite( "target_frame", &io::TerrestrialFrameTransformationParameters::targetFrame )
            .def_readwrite( "reference_decimal_year", &io::TerrestrialFrameTransformationParameters::referenceDecimalYear )
            .def_readwrite( "translation_metres", &io::TerrestrialFrameTransformationParameters::translationMetres )
            .def_readwrite( "rotation_radians", &io::TerrestrialFrameTransformationParameters::rotationRadians )
            .def_readwrite( "scale", &io::TerrestrialFrameTransformationParameters::scale )
            .def_readwrite( "translation_rate_metres_per_second",
                            &io::TerrestrialFrameTransformationParameters::translationRateMetresPerSecond )
            .def_readwrite( "rotation_rate_radians_per_second",
                            &io::TerrestrialFrameTransformationParameters::rotationRateRadiansPerSecond )
            .def_readwrite( "scale_rate_per_second", &io::TerrestrialFrameTransformationParameters::scaleRatePerSecond )
            .def_readwrite( "parameter_rate_year_seconds", &io::TerrestrialFrameTransformationParameters::parameterRateYearSeconds )
            .def_readwrite( "rotation_convention", &io::TerrestrialFrameTransformationParameters::rotationConvention )
            .def_readwrite( "original_translation_units", &io::TerrestrialFrameTransformationParameters::originalTranslationUnits )
            .def_readwrite( "original_rotation_units", &io::TerrestrialFrameTransformationParameters::originalRotationUnits )
            .def_readwrite( "original_scale_units", &io::TerrestrialFrameTransformationParameters::originalScaleUnits )
            .def_readwrite( "original_rate_time_units", &io::TerrestrialFrameTransformationParameters::originalRateTimeUnits )
            .def_readwrite( "provenance", &io::TerrestrialFrameTransformationParameters::provenance )
            .def_readwrite( "original_parameter_line", &io::TerrestrialFrameTransformationParameters::originalParameterLine )
            .def_readwrite( "original_rate_line", &io::TerrestrialFrameTransformationParameters::originalRateLine )
            .def_readwrite( "parameter_source_line", &io::TerrestrialFrameTransformationParameters::parameterSourceLine )
            .def_readwrite( "rate_source_line", &io::TerrestrialFrameTransformationParameters::rateSourceLine )
            .def( "validate", &io::validateTerrestrialFrameTransformationParameters );
    py::class_< io::TerrestrialFrameTransformationData >( module, "TerrestrialFrameTransformationData" )
            .def( py::init<>( ) )
            .def_readwrite( "reference_frame", &io::TerrestrialFrameTransformationData::referenceFrame )
            .def_readwrite( "provenance", &io::TerrestrialFrameTransformationData::provenance )
            .def_readwrite( "source_headers", &io::TerrestrialFrameTransformationData::sourceHeaders )
            .def_readwrite( "records", &io::TerrestrialFrameTransformationData::records )
            .def( "get_parameters",
                  &io::TerrestrialFrameTransformationData::getParameters,
                  py::arg( "target_realization" ),
                  py::return_value_policy::reference_internal )
            .def( "validate", &io::validateTerrestrialFrameTransformationData );
    module.def( "read_terrestrial_frame_transformation_data",
                py::overload_cast< const std::string&, const io::TerrestrialFrameTransformationSource& >(
                        &io::readTerrestrialFrameTransformationData ),
                py::arg( "file_path" ),
                py::arg( "source" ) = io::TerrestrialFrameTransformationSource( ) );
    module.def(
            "read_terrestrial_frame_transformation_data_text",
            []( const std::string& text, const std::string& source, io::TerrestrialFrameTransformationSource provenance ) {
                std::istringstream stream( text );
                provenance.source = source;
                return io::readTerrestrialFrameTransformationData( stream, provenance );
            },
            py::arg( "text" ),
            py::arg( "source" ) = "<text>",
            py::arg( "provenance" ) = io::TerrestrialFrameTransformationSource( ) );
    module.def( "official_itrf2014_parameters", &io::readOfficialItrf2014TransformationData );
    module.def( "official_itrf2020_parameters", &io::readOfficialItrf2020TransformationData );
}
}  // namespace terrestrial_frames
}  // namespace environment_data
}  // namespace data_input
}  // namespace tudatpy
