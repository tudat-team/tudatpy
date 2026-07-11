#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif

#include "expose_psf.h"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "scalarTypes.h"
#include "tudat/io/readPsfFile.h"

namespace py = pybind11;
namespace tio = tudat::input_output;

namespace tudatpy
{

namespace data_access
{

namespace tracking
{

namespace psf
{

void expose_psf( py::module& m )
{
    py::enum_< tio::psf::OpticalImageType >( m, "OpticalImageType" )
            .value( "star", tio::psf::OpticalImageType::star )
            .value( "planet", tio::psf::OpticalImageType::planet )
            .value( "satellite", tio::psf::OpticalImageType::satellite )
            .value( "rock", tio::psf::OpticalImageType::rock )
            .value( "end_marker", tio::psf::OpticalImageType::end_marker )
            .value( "unknown", tio::psf::OpticalImageType::unknown )
            .export_values( );

    py::class_< tio::psf::RawPsfMeasurement, std::shared_ptr< tio::psf::RawPsfMeasurement > >( m, "RawPsfMeasurement" )
            .def_readwrite( "optical_image_type", &tio::psf::RawPsfMeasurement::opticalImageType_ )
            .def_readwrite( "image_name", &tio::psf::RawPsfMeasurement::imageName_ )
            .def_readwrite( "image_id", &tio::psf::RawPsfMeasurement::imageId_ )
            .def_readwrite( "use_flag", &tio::psf::RawPsfMeasurement::useFlag_ )
            .def_readwrite( "observed_pixel_line", &tio::psf::RawPsfMeasurement::observedPixelLine_ )
            .def_readwrite( "local_correction", &tio::psf::RawPsfMeasurement::localCorrection_ )
            .def_readwrite( "sigma_pixel_line", &tio::psf::RawPsfMeasurement::sigmaPixelLine_ )
            .def( "get_effective_pixel_line", &tio::psf::RawPsfMeasurement::getEffectivePixelLine );

    py::class_< tio::psf::RawPsfStarMeasurement, std::shared_ptr< tio::psf::RawPsfStarMeasurement >, tio::psf::RawPsfMeasurement >(
            m, "RawPsfStarMeasurement" )
            .def_readwrite( "right_ascension_degrees", &tio::psf::RawPsfStarMeasurement::rightAscensionDegrees_ )
            .def_readwrite( "declination_degrees", &tio::psf::RawPsfStarMeasurement::declinationDegrees_ );

    py::class_< tio::psf::RawPsfFileImageContents >( m, "RawPsfFileImageContents" )
            .def_readwrite( "picture_name", &tio::psf::RawPsfFileImageContents::pictureName_ )
            .def_readwrite( "picture_number", &tio::psf::RawPsfFileImageContents::pictureNumber_ )
            .def_readwrite( "end_of_exposure_time_utc_string", &tio::psf::RawPsfFileImageContents::endOfExposureTimeUtcString_ )
            .def_readwrite( "camera_id", &tio::psf::RawPsfFileImageContents::cameraId_ )
            .def_readwrite( "exposure_time_seconds", &tio::psf::RawPsfFileImageContents::exposureTimeSeconds_ )
            .def_readwrite( "picture_deletion_flag", &tio::psf::RawPsfFileImageContents::pictureDeletionFlag_ )
            .def_readwrite( "right_ascension_degrees", &tio::psf::RawPsfFileImageContents::rightAscensionDegrees_ )
            .def_readwrite( "declination_degrees", &tio::psf::RawPsfFileImageContents::declinationDegrees_ )
            .def_readwrite( "twist_degrees", &tio::psf::RawPsfFileImageContents::twistDegrees_ )
            .def_readwrite( "measurements", &tio::psf::RawPsfFileImageContents::measurements_ );

    py::class_< tio::psf::RawPsfFileContents >( m, "RawPsfFileContents" )
            .def_readwrite( "spacecraft_id", &tio::psf::RawPsfFileContents::spacecraftId_ )
            .def_readwrite( "number_of_cameras", &tio::psf::RawPsfFileContents::numberOfCameras_ )
            .def_readwrite( "equinox", &tio::psf::RawPsfFileContents::equinox_ )
            .def_readwrite( "psf_id", &tio::psf::RawPsfFileContents::psfId_ )
            .def_readwrite( "psf_program", &tio::psf::RawPsfFileContents::psfProgram_ )
            .def_readwrite( "psf_generation_time_utc_string", &tio::psf::RawPsfFileContents::psfGenerationTimeUtcString_ )
            .def_readwrite( "psf_comments", &tio::psf::RawPsfFileContents::psfComments_ )
            .def_readwrite( "tracking_supplementary_data", &tio::psf::RawPsfFileContents::trackingSupplementaryData_ )
            .def_readwrite( "images", &tio::psf::RawPsfFileContents::images_ );

    m.def( "read_raw_psf_file_contents", &tio::psf::readRawPsfFile, py::arg( "psf_file" ) );
    m.def( "read_psf_files",
           &tio::psf::readPsfFiles< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "psf_file_names" ),
           py::arg( "receiver_body_name" ) = "",
           py::arg( "image_name_to_body_name" ) = std::map< std::string, std::string >( ),
           py::arg( "use_raw_image_name_as_body_name_if_unmapped" ) = true,
           py::arg( "use_corrected_pixel_line" ) = true,
           py::arg( "use_mid_exposure_time" ) = true,
           py::arg( "include_deleted_pictures" ) = false,
           py::arg( "include_end_marker_records" ) = false,
           py::arg( "filter_by_use_flag" ) = false,
           py::arg( "required_use_flag" ) = 0 );
    m.def( "read_psf_file",
           &tio::psf::readPsfFile< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "psf_file_name" ),
           py::arg( "receiver_body_name" ) = "",
           py::arg( "image_name_to_body_name" ) = std::map< std::string, std::string >( ),
           py::arg( "use_raw_image_name_as_body_name_if_unmapped" ) = true,
           py::arg( "use_corrected_pixel_line" ) = true,
           py::arg( "use_mid_exposure_time" ) = true,
           py::arg( "include_deleted_pictures" ) = false,
           py::arg( "include_end_marker_records" ) = false,
           py::arg( "filter_by_use_flag" ) = false,
           py::arg( "required_use_flag" ) = 0 );
}

}  // namespace psf

}  // namespace tracking

}  // namespace data_access

}  // namespace tudatpy
