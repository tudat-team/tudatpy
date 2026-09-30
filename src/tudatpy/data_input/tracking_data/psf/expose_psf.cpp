#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif

#include "expose_psf.h"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <optional>

#include "scalarTypes.h"
#include "tudat/io/readPsfFile.h"

namespace py = pybind11;
namespace tio = tudat::input_output;

namespace tudatpy
{

namespace data_input
{

namespace tracking_data
{

namespace psf
{

void expose_psf( py::module& m )
{
    py::enum_< tio::psf::OpticalImageType >( m, "OpticalImageType", R"doc(Image-object type encoded in a PSF file.)doc" )
            .value( "star", tio::psf::OpticalImageType::star )
            .value( "planet", tio::psf::OpticalImageType::planet )
            .value( "satellite", tio::psf::OpticalImageType::satellite )
            .value( "rock", tio::psf::OpticalImageType::rock )
            .value( "end_marker", tio::psf::OpticalImageType::end_marker )
            .value( "unknown", tio::psf::OpticalImageType::unknown )
            .export_values( );

    py::class_< tio::psf::RawPsfMeasurement, std::shared_ptr< tio::psf::RawPsfMeasurement > >(
            m, "RawPsfMeasurement", R"doc(Raw PSF image measurement record.)doc" )
            .def_readwrite( "optical_image_type", &tio::psf::RawPsfMeasurement::opticalImageType_ )
            .def_readwrite( "image_name", &tio::psf::RawPsfMeasurement::imageName_ )
            .def_readwrite( "image_id", &tio::psf::RawPsfMeasurement::imageId_ )
            .def_readwrite( "use_flag", &tio::psf::RawPsfMeasurement::useFlag_ )
            .def_readwrite( "observed_pixel_line", &tio::psf::RawPsfMeasurement::observedPixelLine_ )
            .def_readwrite( "local_correction", &tio::psf::RawPsfMeasurement::localCorrection_ )
            .def_readwrite( "sigma_pixel_line", &tio::psf::RawPsfMeasurement::sigmaPixelLine_ )
            .def( "get_effective_pixel_line",
                  &tio::psf::RawPsfMeasurement::getEffectivePixelLine,
                  R"doc(
         Return the corrected or effective pixel-line coordinates.

         Returns
         -------
         numpy.ndarray
             Effective pixel-line coordinates.
      )doc" );

    py::class_< tio::psf::RawPsfStarMeasurement, std::shared_ptr< tio::psf::RawPsfStarMeasurement >, tio::psf::RawPsfMeasurement >(
            m, "RawPsfStarMeasurement", R"doc(Raw PSF star measurement record.)doc" )
            .def_readwrite( "right_ascension_degrees", &tio::psf::RawPsfStarMeasurement::rightAscensionDegrees_ )
            .def_readwrite( "declination_degrees", &tio::psf::RawPsfStarMeasurement::declinationDegrees_ );

    py::class_< tio::psf::RawPsfFileImageContents >( m, "RawPsfFileImageContents", R"doc(Raw PSF contents for a single image.)doc" )
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

    py::class_< tio::psf::RawPsfFileContents >( m, "RawPsfFileContents", R"doc(Raw contents of a parsed PSF file.)doc" )
            .def_readwrite( "spacecraft_id", &tio::psf::RawPsfFileContents::spacecraftId_ )
            .def_readwrite( "number_of_cameras", &tio::psf::RawPsfFileContents::numberOfCameras_ )
            .def_readwrite( "equinox", &tio::psf::RawPsfFileContents::equinox_ )
            .def_readwrite( "psf_id", &tio::psf::RawPsfFileContents::psfId_ )
            .def_readwrite( "psf_program", &tio::psf::RawPsfFileContents::psfProgram_ )
            .def_readwrite( "psf_generation_time_utc_string", &tio::psf::RawPsfFileContents::psfGenerationTimeUtcString_ )
            .def_readwrite( "psf_comments", &tio::psf::RawPsfFileContents::psfComments_ )
            .def_readwrite( "tracking_supplementary_data", &tio::psf::RawPsfFileContents::trackingSupplementaryData_ )
            .def_readwrite( "images", &tio::psf::RawPsfFileContents::images_ );

    m.def( "read_raw_psf_file_contents",
           &tio::psf::readRawPsfFile,
           py::arg( "psf_file" ),
           R"doc(
         Read a PSF file without converting it to tracking data.

         This is a supporting reader for inspecting raw PSF contents. In the
         typical Tudat workflow, call :func:`read_psf_data` instead, so the file
         contents are converted to Tudat tracking-data objects.

         Parameters
         ----------
         psf_file : str
             Path to the PSF file.

         Returns
         -------
         RawPsfFileContents
             Raw parsed PSF file contents.
      )doc" );
    m.def(
            "read_psf_data",
            []( const std::vector< std::string >& psfFileNames,
                const std::string& receiverBodyName,
                const std::map< std::string, std::string >& imageNameToBodyName,
                const bool useRawImageNameAsBodyNameIfUnmapped,
                const std::optional< int >& requiredUseFlag ) {
                return tio::psf::readPsfFiles< STATE_SCALAR_TYPE, TIME_TYPE >( psfFileNames,
                                                                               receiverBodyName,
                                                                               imageNameToBodyName,
                                                                               useRawImageNameAsBodyNameIfUnmapped,
                                                                               true,
                                                                               true,
                                                                               false,
                                                                               false,
                                                                               requiredUseFlag.has_value( ),
                                                                               requiredUseFlag.value_or( 0 ) );
            },
            py::arg( "psf_file_names" ),
            py::arg( "receiver_body_name" ) = "",
            py::arg( "image_name_to_body_name" ) = std::map< std::string, std::string >( ),
            py::arg( "use_raw_image_name_as_body_name_if_unmapped" ) = true,
            py::arg( "required_use_flag" ) = py::none( ),
            R"doc(
         Read PSF files and convert them to tracking-data containers.

         This reader is intended for PSF files of the type distributed by the
         JPL Solar System Dynamics group as spacecraft optical observations of
         planetary satellites (https://ssd.jpl.nasa.gov/sats/obs_data.html).
         These files contain optical image measurements, camera metadata, and
         image orientation information. This reader converts the selected image
         records to pixel-line tracking data and corresponding supplementary
         data.
         The raw PSF contents can also be inspected with
         :func:`read_raw_psf_file_contents`; :func:`read_psf_data` applies the
         image-name mapping and optional use-flag filter before creating Tudat
         ``TrackingData`` and
         ``TrackingSupplementaryData`` objects.
         During conversion, deleted pictures and end-marker records are skipped,
         corrected pixel-line coordinates are used, and observation epochs are
         stored at mid-exposure.

         Parameters
         ----------
         psf_file_names : list[str]
             Paths to PSF files.
         receiver_body_name : str, default=""
             Optional receiving body name.
         image_name_to_body_name : dict[str, str], default={}
             Mapping from PSF image names to Tudat body names.
         use_raw_image_name_as_body_name_if_unmapped : bool, default=True
             Whether unmapped image names are used as body names.
         required_use_flag : int | None, default=None
             If not None, only measurements with this PSF ``USE`` flag are
             converted.

         Returns
         -------
         tuple[list[TrackingData], list[TrackingSupplementaryData]]
             Tracking data objects and supplementary data objects.
      )doc" );
}

}  // namespace psf

}  // namespace tracking_data

}  // namespace data_input

}  // namespace tudatpy
