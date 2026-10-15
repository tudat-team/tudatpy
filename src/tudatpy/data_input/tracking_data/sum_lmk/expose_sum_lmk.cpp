#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif

#include "expose_sum_lmk.h"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "scalarTypes.h"
#include "tudat/io/readSumLmkFiles.h"

namespace py = pybind11;
namespace tio = tudat::input_output;

namespace tudatpy
{

namespace data_input
{

namespace tracking_data
{

namespace sum_lmk
{

void expose_sum_lmk( py::module& m )
{
    py::class_< tio::sum_lmk::SumLandmarkObservation >( m, "SumLandmarkObservation", R"doc(One landmark measurement in a SUM file.)doc" )
            .def_readonly( "landmark_id", &tio::sum_lmk::SumLandmarkObservation::landmarkId_ )
            .def_readonly( "pixel_coordinates", &tio::sum_lmk::SumLandmarkObservation::pixelCoordinates_ );

    py::class_< tio::sum_lmk::SumImageData >( m, "SumImageData", R"doc(Raw contents of one SUM file, which describes one image.)doc" )
            .def_readonly( "image_id", &tio::sum_lmk::SumImageData::imageId_ )
            .def_readonly( "utc_epoch_string", &tio::sum_lmk::SumImageData::utcEpochString_ )
            .def_readonly( "image_size", &tio::sum_lmk::SumImageData::imageSize_ )
            .def_readonly( "focal_length_mm", &tio::sum_lmk::SumImageData::focalLengthMm_ )
            .def_readonly( "optical_center", &tio::sum_lmk::SumImageData::opticalCenter_ )
            .def_readonly( "spacecraft_object_vector", &tio::sum_lmk::SumImageData::spacecraftObjectVector_ )
            .def_readonly( "camera_axes", &tio::sum_lmk::SumImageData::cameraAxes_ )
            .def_readonly( "k_matrix", &tio::sum_lmk::SumImageData::kMatrix_ )
            .def_readonly( "pointing_sigma", &tio::sum_lmk::SumImageData::pointingSigma_ )
            .def_readonly( "landmark_observations", &tio::sum_lmk::SumImageData::landmarkObservations_ )
            .def_readonly( "source_file", &tio::sum_lmk::SumImageData::sourceFile_ );

    py::class_< tio::sum_lmk::LmkLandmarkData >(
            m, "LmkLandmarkData", R"doc(Raw contents of one LMK file, which describes one landmark.)doc" )
            .def_readonly( "landmark_id", &tio::sum_lmk::LmkLandmarkData::landmarkId_ )
            .def_readonly( "body_fixed_position", &tio::sum_lmk::LmkLandmarkData::bodyFixedPosition_ )
            .def_readonly( "landmark_position_sigma", &tio::sum_lmk::LmkLandmarkData::landmarkPositionSigma_ )
            .def_readonly( "local_x_axis", &tio::sum_lmk::LmkLandmarkData::localXAxis_ )
            .def_readonly( "local_y_axis", &tio::sum_lmk::LmkLandmarkData::localYAxis_ )
            .def_readonly( "local_z_axis", &tio::sum_lmk::LmkLandmarkData::localZAxis_ )
            .def_readonly( "source_file", &tio::sum_lmk::LmkLandmarkData::sourceFile_ );

    m.def( "read_raw_sum_files",
           &tio::sum_lmk::readSumFiles,
           py::arg( "sum_files" ),
           R"doc(
         Read SUM files without converting them to tracking data.

         This is a supporting reader for inspecting raw SUM contents, or for
         converting the same files several times with
         :func:`convert_raw_sum_lmk_data`. In the typical Tudat workflow, call
         :func:`read_sum_lmk_data` instead.

         Parameters
         ----------
         sum_files : list[str]
             Paths to the SUM files.

         Returns
         -------
         list[SumImageData]
             Raw contents of each image.
      )doc" );

    m.def( "read_raw_lmk_files",
           &tio::sum_lmk::readLmkFiles,
           py::arg( "lmk_files" ),
           R"doc(
         Read LMK files without converting them to tracking data.

         Parameters
         ----------
         lmk_files : list[str]
             Paths to the LMK files.

         Returns
         -------
         dict[str, LmkLandmarkData]
             Raw contents of each landmark, by landmark ID.
      )doc" );

    m.def( "convert_raw_sum_lmk_data",
           &tio::sum_lmk::convertSumLmkData< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "sum_images" ),
           py::arg( "landmarks" ),
           py::arg( "target_body_name" ),
           py::arg( "receiver_body_name" ),
           py::arg( "skip_observations_with_missing_landmarks" ) = false,
           py::arg( "validate_spacecraft_object_geometry" ) = true,
           R"doc(
         Convert raw SUM/LMK contents to tracking-data containers.

         Same as :func:`read_sum_lmk_data`, but starts from contents read with
         :func:`read_raw_sum_files` and :func:`read_raw_lmk_files`.

         Parameters
         ----------
         sum_images : list[SumImageData]
             Raw image contents.
         landmarks : dict[str, LmkLandmarkData]
             Raw landmark contents, by landmark ID.
         target_body_name : str
             Name of the body the landmarks are on.
         receiver_body_name : str
             Name of the body that carries the cameras.
         skip_observations_with_missing_landmarks : bool, default=False
             If True, measurements of landmarks without an LMK file are skipped.
             If False, they raise an error.
         validate_spacecraft_object_geometry : bool, default=True
             If True, check that every measured landmark lies in front of the
             camera, using the spacecraft position stored in the SUM file.

         Returns
         -------
         tuple[list[TrackingData], list[TrackingSupplementaryData]]
             Tracking data objects and supplementary data objects.
      )doc" );

    m.def( "read_sum_lmk_data",
           &tio::sum_lmk::readSumLmkData< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "sum_files" ),
           py::arg( "lmk_files" ),
           py::arg( "target_body_name" ),
           py::arg( "receiver_body_name" ),
           py::arg( "skip_observations_with_missing_landmarks" ) = false,
           py::arg( "validate_spacecraft_object_geometry" ) = true,
           R"doc(
         Read SUM and LMK files and convert them to tracking-data containers.

         SUM and LMK files are written by the stereophotoclinometry (SPC)
         software. Each SUM file describes one image: its epoch, camera
         calibration and pointing, and the pixel coordinates of the landmarks
         seen in it. Each LMK file describes one landmark: its position on the
         target body.

         Every (image, landmark) pair becomes one pixel-coordinate observation,
         sent by the landmark and received by the camera of the image. The
         supplementary data adds each landmark as a ground station on the
         target body, and one camera per image on the receiver body. Each
         camera is named ``Camera_<image ID>``, so its pointing can be
         estimated separately with
         :func:`~tudatpy.dynamics.parameters_setup.camera_pointing_correction`.
         The SUM ``SIGMA_PTG`` value is stored as the camera's
         ``pointing_uncertainty``.

         Apply the supplementary data with
         :func:`~tudatpy.estimation.observations.set_tracking_supplementary_data_in_bodies`
         before creating observations with
         :func:`~tudatpy.estimation.observations.create_observation_collection_from_tracking_data`.

         Parameters
         ----------
         sum_files : list[str]
             Paths to the SUM files.
         lmk_files : list[str]
             Paths to the LMK files.
         target_body_name : str
             Name of the body the landmarks are on.
         receiver_body_name : str
             Name of the body that carries the cameras.
         skip_observations_with_missing_landmarks : bool, default=False
             If True, measurements of landmarks without an LMK file are skipped.
             If False, they raise an error.
         validate_spacecraft_object_geometry : bool, default=True
             If True, check that every measured landmark lies in front of the
             camera, using the spacecraft position stored in the SUM file.

         Returns
         -------
         tuple[list[TrackingData], list[TrackingSupplementaryData]]
             Tracking data objects and supplementary data objects.
      )doc" );
}

}  // namespace sum_lmk

}  // namespace tracking_data

}  // namespace data_input

}  // namespace tudatpy
