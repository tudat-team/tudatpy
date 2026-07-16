#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif

#include "expose_tracking_data.h"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "fdets/expose_fdets.h"
#include "generic_text_file/expose_generic_text_file.h"
#include "ifms/expose_ifms.h"
#include "odf/expose_odf.h"
#include "psf/expose_psf.h"
#include "slr/expose_slr.h"
#include "scalarTypes.h"
#include "tudat/io/trackingData.h"
#include "tudat/io/trackingSupplementaryData.h"

namespace py = pybind11;
namespace tdat = tudat::data;

namespace tudatpy
{

namespace data_input
{

namespace tracking_data
{

void expose_tracking_data( py::module& m )
{
    py::class_< tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >, std::shared_ptr< tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE > > >(
            m,
            "TrackingData",
            R"doc(Container for source-loaded tracking observations from files and external sources. These objects are converted to observation collections by :func:`~tudatpy.estimation.observations.create_observation_collection_from_tracking_data`.)doc" )
            .def( py::init< const std::string,
                            const tdat::PlainLinkDefinition&,
                            const std::vector< Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >&,
                            const std::vector< TIME_TYPE >,
                            const std::string,
                            const std::string,
                            const std::string >( ),
                  py::arg( "observable_type" ),
                  py::arg( "link_ends" ),
                  py::arg( "observations" ),
                  py::arg( "epochs" ),
                  py::arg( "reference_link_end" ),
                  py::arg( "time_scale" ) = "TDB",
                  py::arg( "weighing_scheme" ) = "",
                  R"doc(
         Create a container for source-loaded tracking observations. Note that this constructor is typically not called directly by users,
         but that an instance of this class is created by one of the functions in this module's submodules.

         Parameters
         ----------
         observable_type : str
             Observable type represented by this data set.
         link_ends : list
             Link-end definition associated with the observations.
         observations : list
             Observation vectors.
         epochs : list
             Observation epochs.
         reference_link_end : str
             Reference link end for the observation epochs.
         time_scale : str, default "TDB"
             Time scale of the observation epochs.
         weighing_scheme : str, default ""
             Name of the weighing scheme assigned to the observations.
      )doc" )
            .def_property_readonly( "observable_type",
                                    &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getObservableType,
                                    R"doc(
         **read-only**

         Observable type represented by this data set.

         :type: str
      )doc" )
            .def_property_readonly( "link_ends",
                                    &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getLinkEnds,
                                    R"doc(
         **read-only**

         Link-end definition associated with the observations.

         :type: list
      )doc" )
            .def_property_readonly( "observations",
                                    &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getObservations,
                                    R"doc(
         **read-only**

         Observation vectors.

         :type: list
      )doc" )
            .def_property_readonly( "epochs",
                                    &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationEpochs,
                                    R"doc(
         **read-only**

         Observation epochs.

         :type: list
      )doc" )
            .def_property_readonly( "reference_link_end",
                                    &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getReferenceLinkEnd,
                                    R"doc(
         **read-only**

         Reference link end for the observation epochs.

         :type: str
      )doc" )
            .def_property_readonly( "time_scale",
                                    &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getTimeScale,
                                    R"doc(
         **read-only**

         Time scale of the observation epochs.

         :type: str
      )doc" )
            .def_property_readonly( "weighing_scheme",
                                    &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getWeighingScheme,
                                    R"doc(
         **read-only**

         Name of the weighing scheme assigned to the observations.

         :type: str
      )doc" )
            .def( "add_string_ancillary_setting",
                  py::overload_cast< const std::string, const std::string >(
                          &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::addAncillarySettings ),
                  py::arg( "ancillary_settings_type" ),
                  py::arg( "ancillary_settings_value" ),
                  R"doc(
         Add a string ancillary setting.

         Parameters
         ----------
         ancillary_settings_type : str
             Ancillary setting identifier.
         ancillary_settings_value : str
             Ancillary setting value.

         Returns
         -------
         None
      )doc" )
            .def( "add_string_vector_ancillary_setting",
                  py::overload_cast< const std::string, const std::vector< std::string > >(
                          &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::addAncillarySettings ),
                  py::arg( "ancillary_settings_type" ),
                  py::arg( "ancillary_settings_value" ),
                  R"doc(
         Add a string-vector ancillary setting.

         Parameters
         ----------
         ancillary_settings_type : str
             Ancillary setting identifier.
         ancillary_settings_value : list[str]
             Ancillary setting values.

         Returns
         -------
         None
      )doc" )
            .def( "add_double_ancillary_setting",
                  py::overload_cast< const std::string, const double >(
                          &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::addAncillarySettings ),
                  py::arg( "ancillary_settings_type" ),
                  py::arg( "ancillary_settings_value" ),
                  R"doc(
         Add a floating-point ancillary setting.

         Ancillary settings present in supported tracking-data files are loaded automatically by the readers.
         This function is intended only to override or extend those settings manually, which is not part of
         the typical data-loading workflow. See :ref:`ancillary_settings` for how ancillary settings are used
         by observation models.

         Parameters
         ----------
         ancillary_settings_type : str
             Ancillary setting identifier.
         ancillary_settings_value : float
             Ancillary setting value.

         Returns
         -------
         None
      )doc" )
            .def( "add_double_vector_ancillary_setting",
                  py::overload_cast< const std::string, const std::vector< double > >(
                          &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::addAncillarySettings ),
                  py::arg( "ancillary_settings_type" ),
                  py::arg( "ancillary_settings_value" ),
                  R"doc(
         Add a floating-point-vector ancillary setting.

         Ancillary settings present in supported tracking-data files are loaded automatically by the readers.
         This function is intended only to override or extend those settings manually, which is not part of
         the typical data-loading workflow. See :ref:`ancillary_settings` for how ancillary settings are used
         by observation models.

         Parameters
         ----------
         ancillary_settings_type : str
             Ancillary setting identifier.
         ancillary_settings_value : list[float]
             Ancillary setting values.

         Returns
         -------
         None
      )doc" )
            .def( "get_ancillary_settings_string_vector",
                  &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getAncillarySettingsStringVector,
                  R"doc(
         Return a string-vector ancillary setting.

         Parameters
         ----------
         ancillary_settings_type : str
             Ancillary setting identifier.

         Returns
         -------
         list[str]
             Ancillary setting values.
      )doc" )
            .def( "get_ancillary_settings_double",
                  &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getAncillarySettingsDouble,
                  R"doc(
         Return a floating-point ancillary setting.

         Parameters
         ----------
         ancillary_settings_type : str
             Ancillary setting identifier.

         Returns
         -------
         float
             Ancillary setting value.
      )doc" )
            .def( "get_ancillary_settings_double_vector",
                  &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getAncillarySettingsDoubleVector,
                  R"doc(
         Return a floating-point-vector ancillary setting.

         Parameters
         ----------
         ancillary_settings_type : str
             Ancillary setting identifier.

         Returns
         -------
         list[float]
             Ancillary setting values.
      )doc" )
            .def( "set_observation_weights",
                  &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::setObservationWeights,
                  py::arg( "observation_weights" ),
                  R"doc(
         Set per-observation weights.

         Parameters
         ----------
         observation_weights : list
             Observation weights.

         Returns
         -------
         None
      )doc" )
            .def( "reset_single_observation_weight",
                  &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::setSingleObservationWeight,
                  py::arg( "index" ),
                  py::arg( "observation_weight" ),
                  R"doc(
         Reset one observation weight.

         Parameters
         ----------
         index : int
             Observation index.
         observation_weight : float
             Replacement observation weight.

         Returns
         -------
         None
      )doc" )
            .def( "get_observation_weights",
                  &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationWeights,
                  R"doc(
         Return per-observation weights.

         Returns
         -------
         list
             Observation weights.
      )doc" )
            .def( "get_concatenated_observation_weights",
                  &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationWeightsVector,
                  R"doc(
         Return concatenated observation weights.

         Returns
         -------
         list[float]
             Concatenated observation weights.
      )doc" )
            .def( "set_observation_corrections",
                  &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::setObservationCorrections,
                  py::arg( "observation_corrections" ),
                  R"doc(
         Set per-observation corrections.

         Parameters
         ----------
         observation_corrections : list
             Observation correction vectors.

         Returns
         -------
         None
      )doc" )
            .def( "reset_single_observation_correction",
                  &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::setSingleObservationCorrection,
                  py::arg( "index" ),
                  py::arg( "observation_correction" ),
                  R"doc(
         Reset one observation correction.

         Parameters
         ----------
         index : int
             Observation index.
         observation_correction : numpy.ndarray
             Replacement observation correction vector.

         Returns
         -------
         None
      )doc" )
            .def( "get_observation_corrections",
                  &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationCorrections,
                  R"doc(
         Return per-observation corrections.

         Returns
         -------
         list
             Observation correction vectors.
      )doc" )
            .def( "remove_single_observation_entry",
                  &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::removeSingleObservationEntry,
                  py::arg( "index" ),
                  R"doc(
         Remove one observation entry.

         Parameters
         ----------
         index : int
             Observation index to remove.

         Returns
         -------
         None
      )doc" );

    py::class_< tdat::FrequencySupplementaryData, std::shared_ptr< tdat::FrequencySupplementaryData > >(
            m, "FrequencySupplementaryData", R"doc(
         Base class for frequency supplementary data.

         Derived classes store station transmission-frequency information. When
         applied with
         :func:`~tudatpy.estimation.observations.set_tracking_supplementary_data_in_bodies`,
         this data is used to set frequency histories on the associated
         :class:`~tudatpy.dynamics.environment.GroundStation`, or on the
         associated body's vehicle systems if no reference point is specified.
         Ground-station creation is described in :ref:`ground_station` and by
         :class:`~tudatpy.dynamics.environment_setup.ground_station.GroundStationSettings`;
         vehicle-system setup is described in :ref:`vehicle_systems`.
      )doc" );

    py::class_< tdat::RampedFrequencySupplementaryData::FrequencyRamp >( m, "FrequencyRamp", R"doc(
         Container for one ground-station transmission-frequency ramp.

         Frequency ramps are stored in
         :class:`RampedFrequencySupplementaryData` and are applied to the
         associated :class:`~tudatpy.dynamics.environment.GroundStation`, or to
         the associated body's vehicle systems if no reference point is
         specified. See :ref:`ground_station` and :ref:`vehicle_systems` for
         the corresponding environment models.
      )doc" )
            .def( py::init<>( ) )
            .def( py::init< const double, const double, const double, const double >( ),
                  py::arg( "start_time" ),
                  py::arg( "end_time" ),
                  py::arg( "start_frequency" ),
                  py::arg( "frequency_rate" ) )
            .def_readwrite( "start_time", &tdat::RampedFrequencySupplementaryData::FrequencyRamp::startTime_, R"doc(
         **read-only**

         Ramp start time.

         :type: float
      )doc" )
            .def_readwrite( "end_time", &tdat::RampedFrequencySupplementaryData::FrequencyRamp::endTime_, R"doc(
         **read-only**

         Ramp end time.

         :type: float
      )doc" )
            .def_readwrite( "start_frequency", &tdat::RampedFrequencySupplementaryData::FrequencyRamp::startFrequency_, R"doc(
         **read-only**

         Ramp start frequency.

         :type: float
      )doc" )
            .def_readwrite( "frequency_rate", &tdat::RampedFrequencySupplementaryData::FrequencyRamp::frequencyRate_, R"doc(
         **read-only**

         Ramp frequency rate.

         :type: float
      )doc" );

    py::class_< tdat::RampedFrequencySupplementaryData,
                std::shared_ptr< tdat::RampedFrequencySupplementaryData >,
                tdat::FrequencySupplementaryData >( m, "RampedFrequencySupplementaryData", R"doc(
         Container for ramped-frequency supplementary tracking data.

         This class stores piecewise-linear transmission-frequency ramps. When
         applied with
         :func:`~tudatpy.estimation.observations.set_tracking_supplementary_data_in_bodies`,
         the ramps are set on the associated
         :class:`~tudatpy.dynamics.environment.GroundStation`, or on the
         associated body's vehicle systems if no reference point is specified.
         See :ref:`ground_station`, :ref:`vehicle_systems`, and
         :class:`~tudatpy.dynamics.environment_setup.ground_station.GroundStationSettings`
         for the environment model to which this data is attached.
      )doc" )
            .def( py::init<>( ) )
            .def( py::init< const std::vector< tdat::RampedFrequencySupplementaryData::FrequencyRamp >& >( ), py::arg( "frequency_ramps" ) )
            .def_property_readonly( "frequency_ramps",
                                    &tdat::RampedFrequencySupplementaryData::getFrequencyRamps,
                                    R"doc(
         **read-only**

         Frequency-ramp records.

         :type: list[FrequencyRamp]
      )doc" )
            .def( "add_frequency_ramp",
                  &tdat::RampedFrequencySupplementaryData::addFrequencyRamp,
                  py::arg( "start_time" ),
                  py::arg( "end_time" ),
                  py::arg( "start_frequency" ),
                  py::arg( "frequency_rate" ),
                  R"doc(
         Add one frequency ramp.

         Parameters
         ----------
         start_time : float
             Ramp start time.
         end_time : float
             Ramp end time.
         start_frequency : float
             Ramp start frequency.
         frequency_rate : float
             Ramp frequency rate.

         Returns
         -------
         None
      )doc" );

    py::class_< tdat::PiecewiseConstantFrequencySupplementaryData,
                std::shared_ptr< tdat::PiecewiseConstantFrequencySupplementaryData >,
                tdat::FrequencySupplementaryData >( m, "PiecewiseConstantFrequencySupplementaryData", R"doc(
         Container for piecewise-constant frequency supplementary data.

         This class stores a transmission-frequency history that is constant
         between tabulated epochs. When applied with
         :func:`~tudatpy.estimation.observations.set_tracking_supplementary_data_in_bodies`,
         the history is set on the associated
         :class:`~tudatpy.dynamics.environment.GroundStation`, or on the
         associated body's vehicle systems if no reference point is specified.
         See :ref:`ground_station`, :ref:`vehicle_systems`, and
         :class:`~tudatpy.dynamics.environment_setup.ground_station.GroundStationSettings`
         for the environment model to which this data is attached.
      )doc" )
            .def( py::init<>( ) )
            .def( py::init< const std::map< double, double >& >( ), py::arg( "frequency_history" ) )
            .def_property_readonly( "frequency_history",
                                    &tdat::PiecewiseConstantFrequencySupplementaryData::getFrequencyHistory,
                                    R"doc(
         **read-only**

         Frequency history indexed by time.

         :type: dict[float, float]
      )doc" )
            .def( "set_frequency",
                  &tdat::PiecewiseConstantFrequencySupplementaryData::setFrequency,
                  py::arg( "time" ),
                  py::arg( "frequency" ),
                  R"doc(
         Set one piecewise-constant frequency entry.

         Parameters
         ----------
         time : float
             Time at which the frequency applies.
         frequency : float
             Frequency value.

         Returns
         -------
         None
      )doc" );

    py::class_< tdat::InstrumentSupplementaryData, std::shared_ptr< tdat::InstrumentSupplementaryData > >(
            m, "InstrumentSupplementaryData", R"doc(
         Base class for instrument supplementary data.

         Derived classes store instrument information that is converted to
         environment-system models when applied to the bodies. For camera data,
         this creates camera settings for the associated body; see
         :class:`~tudatpy.dynamics.environment_setup.vehicle_systems.CameraSettings`
         and :func:`~tudatpy.dynamics.environment_setup.vehicle_systems.pinhole_camera`.
      )doc" )
            .def_property_readonly( "instrument_supplementary_data_type",
                                    &tdat::InstrumentSupplementaryData::getInstrumentSupplementaryDataKind,
                                    R"doc(
         **read-only**

         Type of instrument supplementary data.

         :type: str
      )doc" );

    py::class_< tdat::CameraInstrumentSupplementaryData,
                std::shared_ptr< tdat::CameraInstrumentSupplementaryData >,
                tdat::InstrumentSupplementaryData >( m, "CameraInstrumentSupplementaryData", R"doc(
         Container for camera instrument supplementary data.

         This class stores camera calibration, mounting and optional orientation
         history. When applied with
         :func:`~tudatpy.estimation.observations.set_tracking_supplementary_data_in_bodies`,
         it is converted to camera data for the associated body. See
         :class:`~tudatpy.dynamics.environment_setup.vehicle_systems.CameraSettings`
         and :func:`~tudatpy.dynamics.environment_setup.vehicle_systems.pinhole_camera`
         for the corresponding environment setup interface.
      )doc" )
            .def( py::init<>( ) )
            .def( py::init< const std::string&,
                            const double,
                            const Eigen::Vector2d&,
                            const Eigen::Vector4d&,
                            const Eigen::Matrix< double, 2, 3 >&,
                            const Eigen::Matrix< double, 6, 1 >&,
                            const Eigen::Vector3d&,
                            const std::map< double, Eigen::Quaterniond >& >( ),
                  py::arg( "camera_id" ),
                  py::arg( "focal_length" ),
                  py::arg( "principal_point" ),
                  py::arg( "field_of_view_bounds" ),
                  py::arg( "k_matrix" ),
                  py::arg( "distortion_coefficients" ),
                  py::arg( "mounting_offsets" ),
                  py::arg( "rotation_from_inertial_to_camera_frame_history" ) = std::map< double, Eigen::Quaterniond >( ) )
            .def_property_readonly( "camera_id", &tdat::CameraInstrumentSupplementaryData::getCameraId, R"doc(
         **read-only**

         Camera identifier.

         :type: str
      )doc" )
            .def_property_readonly( "focal_length", &tdat::CameraInstrumentSupplementaryData::getFocalLength, R"doc(
         **read-only**

         Camera focal length.

         :type: float
      )doc" )
            .def_property_readonly( "principal_point", &tdat::CameraInstrumentSupplementaryData::getPrincipalPoint, R"doc(
         **read-only**

         Camera principal point.

         :type: numpy.ndarray
      )doc" )
            .def_property_readonly( "field_of_view_bounds",
                                    &tdat::CameraInstrumentSupplementaryData::getFieldOfViewBounds,
                                    R"doc(
         **read-only**

         Field-of-view bounds.

         :type: numpy.ndarray
      )doc" )
            .def_property_readonly( "k_matrix", &tdat::CameraInstrumentSupplementaryData::getKMatrix, R"doc(
         **read-only**

         Camera K matrix.

         :type: numpy.ndarray
      )doc" )
            .def_property_readonly( "distortion_coefficients",
                                    &tdat::CameraInstrumentSupplementaryData::getDistortionCoefficients,
                                    R"doc(
         **read-only**

         Camera distortion coefficients.

         :type: numpy.ndarray
      )doc" )
            .def_property_readonly( "mounting_offsets", &tdat::CameraInstrumentSupplementaryData::getMountingOffsets, R"doc(
         **read-only**

         Camera mounting offsets.

         :type: numpy.ndarray
      )doc" )
            .def_property_readonly( "rotation_from_inertial_to_camera_frame_history",
                                    &tdat::CameraInstrumentSupplementaryData::getRotationFromInertialToCameraFrameHistory,
                                    R"doc(
         **read-only**

         Camera-frame rotation history indexed by time.

         :type: dict[float, object]
      )doc" );

    py::class_< tdat::TranslationalStateSupplementaryData >( m, "TranslationalStateSupplementaryData", R"doc(
         Container for translational-state supplementary data.

         This class stores a tabulated Cartesian state history. When applied with
         :func:`~tudatpy.estimation.observations.set_tracking_supplementary_data_in_bodies`,
         the data is converted to a tabulated ephemeris for the associated
         body. Translational-state supplementary data is not applied to a
         reference point. See
         :class:`~tudatpy.dynamics.environment_setup.ephemeris.TabulatedEphemerisSettings`
         and :func:`~tudatpy.dynamics.environment_setup.ephemeris.tabulated`
         for the corresponding environment setup interface.
      )doc" )
            .def( py::init<>( ) )
            .def_property_readonly( "state_history",
                                    &tdat::TranslationalStateSupplementaryData::getStateHistory,
                                    R"doc(
         **read-only**

         Translational state history indexed by time.

         :type: dict[float, numpy.ndarray]
      )doc" )
            .def_property_readonly( "frame_origin",
                                    &tdat::TranslationalStateSupplementaryData::getFrameOrigin,
                                    R"doc(
         **read-only**

         Origin of the translational-state frame.

         :type: str
      )doc" )
            .def_property_readonly( "is_velocity_defined",
                                    &tdat::TranslationalStateSupplementaryData::isVelocityDefined,
                                    R"doc(
         **read-only**

         Boolean that defines whether velocity entries are defined, if not they are computed through finite differences when this object is processed to update the environment.

         :type: bool
      )doc" );

    py::class_< tdat::RotationalStateSupplementaryData >( m, "RotationalStateSupplementaryData", R"doc(
         Container for rotational-state supplementary data.

         This class stores a tabulated rotational state history. When applied
         with
         :func:`~tudatpy.estimation.observations.set_tracking_supplementary_data_in_bodies`,
         the data is converted to a rotational ephemeris for the associated
         body. Rotational-state supplementary data is not applied to a
         reference point. See
         :class:`~tudatpy.dynamics.environment_setup.rotation_model.RotationModelSettings`
         and :func:`~tudatpy.dynamics.environment_setup.rotation_model.custom_rotation_model`
         for the corresponding environment setup interface.
      )doc" )
            .def( py::init<>( ) )
            .def_property_readonly( "rotational_state_history",
                                    &tdat::RotationalStateSupplementaryData::getRotationalStateHistory,
                                    R"doc(
         **read-only**

         Rotational state history indexed by time (as quaternion and body-fixed angular velocity).

         :type: dict[float, numpy.ndarray]
      )doc" )
            .def_property_readonly( "base_frame_orientation",
                                    &tdat::RotationalStateSupplementaryData::getBaseFrameOrientation,
                                    R"doc(
         **read-only**

         Base frame of the rotational state.

         :type: str
      )doc" )
            .def_property_readonly( "is_angular_velocity_defined",
                                    &tdat::RotationalStateSupplementaryData::isAngularVelocityDefined,
                                    R"doc(
         **read-only**

         Whether angular-velocity entries are defined (if not, they are currently also not computed when processing this object to update the environment).

         :type: bool
      )doc" );

    py::class_< tdat::TrackingSupplementaryData, std::shared_ptr< tdat::TrackingSupplementaryData > >(
            m,
            "TrackingSupplementaryData",
            R"doc(Container for tracking supplementary data. These objects store auxiliary information loaded from tracking data files and are applied to a system of bodies by :func:`~tudatpy.estimation.observations.set_tracking_supplementary_data_in_bodies`.)doc" )
            .def( py::init<>( ) )
            .def( py::init< const std::string&, const std::string& >( ), py::arg( "body_name" ), py::arg( "reference_point_name" ) )
            .def_property( "body_name",
                           &tdat::TrackingSupplementaryData::getBodyName,
                           &tdat::TrackingSupplementaryData::setBodyName,
                           R"doc(
         **read-only**

         Body name.

         :type: str
      )doc" )
            .def_property( "reference_point_name",
                           &tdat::TrackingSupplementaryData::getReferencePointName,
                           &tdat::TrackingSupplementaryData::setReferencePointName,
                           R"doc(
         **read-only**

         Reference point name.

         :type: str
      )doc" )
            .def_property( "translational_state_supplementary_data",
                           &tdat::TrackingSupplementaryData::getTranslationalStateSupplementaryData,
                           &tdat::TrackingSupplementaryData::setTranslationalStateSupplementaryData,
                           R"doc(
         **read-only**

         Translational-state supplementary data.

         :type: TranslationalStateSupplementaryData
      )doc" )
            .def_property( "rotational_state_supplementary_data",
                           &tdat::TrackingSupplementaryData::getRotationalStateSupplementaryData,
                           &tdat::TrackingSupplementaryData::setRotationalStateSupplementaryData,
                           R"doc(
         **read-only**

         Rotational-state supplementary data.

         :type: RotationalStateSupplementaryData
      )doc" )
            .def_property( "frequency_supplementary_data",
                           &tdat::TrackingSupplementaryData::getFrequencySupplementaryData,
                           &tdat::TrackingSupplementaryData::setFrequencySupplementaryData,
                           R"doc(
         **read-only**

         Frequency supplementary data objects.

         :type: list[FrequencySupplementaryData]
      )doc" )
            .def_property( "instrument_supplementary_data",
                           &tdat::TrackingSupplementaryData::getInstrumentSupplementaryData,
                           &tdat::TrackingSupplementaryData::setInstrumentSupplementaryData,
                           R"doc(
         **read-only**

         Instrument supplementary data objects.

         :type: list[InstrumentSupplementaryData]
      )doc" )
            .def( "set_frequency_supplementary_data",
                  &tdat::TrackingSupplementaryData::setFrequencySupplementaryData,
                  py::arg( "frequency_supplementary_data" ),
                  R"doc(
         Set frequency supplementary data.

         Parameters
         ----------
         frequency_supplementary_data : list[FrequencySupplementaryData]
             Frequency supplementary data objects.

         Returns
         -------
         None
      )doc" );

    auto genericTextFile = m.def_submodule( "generic_text_file" );
    generic_text_file::expose_generic_text_file( genericTextFile );

    auto odfModule = m.def_submodule( "odf" );
    odf::expose_odf( odfModule );

    auto ifmsModule = m.def_submodule( "ifms" );
    ifms::expose_ifms( ifmsModule );

    auto fdetsModule = m.def_submodule( "fdets" );
    fdets::expose_fdets( fdetsModule );

    auto psfModule = m.def_submodule( "psf" );
    psf::expose_psf( psfModule );

    auto slrModule = m.def_submodule( "slr" );
    slr::expose_slr( slrModule );
}

}  // namespace tracking_data

}  // namespace data_input

}  // namespace tudatpy
