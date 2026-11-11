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
#include "landmarks/expose_landmarks.h"
#include "scalarTypes.h"
#include "tudat/io/trackingData.h"
#include "tudat/io/trackingDataTime.h"
#include "tudat/io/trackingSupplementaryData.h"

namespace py = pybind11;
namespace tdat = tudat::data;
namespace tom = tudat::observation_models;

namespace tudatpy
{

namespace data_input
{

namespace tracking_data
{

void expose_tracking_data( py::module& m )
{
    m.def( "get_tracking_data_epoch_bounds",
           &tdat::getTrackingDataEpochBounds< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "tracking_data" ),
           py::arg( "output_time_scale" ) = tudat::basic_astrodynamics::tdb_scale,
           R"doc(
         Return the earliest and latest observation epochs of a list of tracking data objects.

         The epochs of each tracking data object are converted from its own time scale to ``output_time_scale``. These time
         scale conversions are performed at the geocenter, so that no system of bodies is required. Tracking data objects
         without observations are skipped.

         Parameters
         ----------
         tracking_data : list[TrackingData]
             Tracking data objects, of which at least one must contain observations.
         output_time_scale : tudatpy.astro.time_representation.TimeScales, default = tdb_scale
             Time scale of the returned epochs (TAI, TT, TDB, UTC or UT1).

         Returns
         -------
         tuple[tudatpy.astro.time_representation.Time, tudatpy.astro.time_representation.Time]
             Earliest and latest observation epoch, in seconds since J2000 in ``output_time_scale``.
      )doc" );

    py::class_< tdat::AngularObservationCorrectionSettings >( m,
                                                              "AngularObservationCorrectionSettings",
                                                              R"doc(
         Class defining the corrections to angular position observations that are computed when converting tracking data to an
         observation dataset.

         Objects of this class are provided to :meth:`TrackingData.set_observation_correction_settings`. The corrections are
         computed and applied by :func:`~tudatpy.estimation.observations.create_observation_dataset_from_tracking_data` (if
         its ``apply_corrections`` input is ``True``), using the environment in the system of bodies, with the functions in
         :ref:`observation_corrections`. The observed body and the observer are taken from the link ends of the tracking data.
         The observations stored in the tracking data are not modified.
      )doc" )
            .def( py::init< const std::vector< std::string >&, const Eigen::VectorXd& >( ),
                  py::arg( "light_deflection_bodies" ) = std::vector< std::string >( ),
                  py::arg( "photocenter_body_dimensions" ) = Eigen::VectorXd( ),
                  R"doc(
         Class constructor.

         Parameters
         ----------
         light_deflection_bodies : list[str], default = []
             Names of the bodies for which the gravitational light deflection is corrected (no light deflection correction if
             empty).
         photocenter_body_dimensions : numpy.ndarray[numpy.float64[m, 1]], default = []
             Dimensions of the observed body, used to compute the photocenter correction: empty (no photocenter correction),
             a single radius, or the three semi-axes of an ellipsoid, in meters.
      )doc" )
            .def_readonly( "light_deflection_bodies", &tdat::AngularObservationCorrectionSettings::lightDeflectionBodies_ )
            .def_readonly( "photocenter_body_dimensions", &tdat::AngularObservationCorrectionSettings::photocenterBodyDimensions_ );

    py::class_< tom::ObservationWeightSettings >( m,
                                                  "ObservationWeightSettings",
                                                  R"doc(
         Class defining the weights of the observations in a :class:`TrackingData` object.

         Objects of this class are created with the static functions of this class, and provided to
         :meth:`TrackingData.set_observation_weight_settings` (or to
         :meth:`~tudatpy.estimation.observations.ObservationDataset.add_observation_set_with_weights`). The weights may be
         defined as a single value for all observations, as one value per observation or per observation component, as a matrix
         per observation (including correlations between its components), or as a single matrix for all observations (including
         correlations between different observations).
      )doc" )
            .def( py::init<>( ), R"doc(Class constructor, creating settings in which all weights are equal to 1.)doc" )
            .def_static( "default_weights",
                         &tom::ObservationWeightSettings::defaultWeights,
                         R"doc(
         Function to create settings in which all weights are equal to 1.

         Returns
         -------
         ObservationWeightSettings
             Weight settings.
      )doc" )
            .def_static( "constant_scalar",
                         &tom::ObservationWeightSettings::constantScalar,
                         py::arg( "weight" ),
                         R"doc(
         Function to create settings with the same weight for all components of all observations.

         Parameters
         ----------
         weight : float
             Weight of each observation component.

         Returns
         -------
         ObservationWeightSettings
             Weight settings.
      )doc" )
            .def_static( "scalar_per_observation",
                         &tom::ObservationWeightSettings::scalarPerObservation,
                         py::arg( "weights" ),
                         R"doc(
         Function to create settings with one weight per observation, used for all components of that observation.

         Parameters
         ----------
         weights : list[float]
             Weight of each observation.

         Returns
         -------
         ObservationWeightSettings
             Weight settings.
      )doc" )
            .def_static( "diagonal_per_observation",
                         &tom::ObservationWeightSettings::diagonalPerObservation,
                         py::arg( "weights" ),
                         R"doc(
         Function to create settings with a separate weight for each component of each observation.

         Parameters
         ----------
         weights : list[numpy.ndarray[numpy.float64[m, 1]]]
             Weight of each observation component, with one vector per observation.

         Returns
         -------
         ObservationWeightSettings
             Weight settings.
      )doc" )
            .def_static( "constant_block",
                         &tom::ObservationWeightSettings::constantBlock,
                         py::arg( "weight_block" ),
                         R"doc(
         Function to create settings with the same weight matrix (including correlations between the observation
         components) for all observations.

         Parameters
         ----------
         weight_block : numpy.ndarray[numpy.float64[m, m]]
             Symmetric weight matrix, with one row and column per observation component.

         Returns
         -------
         ObservationWeightSettings
             Weight settings.
      )doc" )
            .def_static( "block_per_observation",
                         &tom::ObservationWeightSettings::blockPerObservation,
                         py::arg( "weight_blocks" ),
                         R"doc(
         Function to create settings with a separate weight matrix (including correlations between the observation
         components) for each observation.

         Parameters
         ----------
         weight_blocks : list[numpy.ndarray[numpy.float64[m, m]]]
             Symmetric weight matrix of each observation, with one row and column per observation component.

         Returns
         -------
         ObservationWeightSettings
             Weight settings.
      )doc" )
            .def_static( "set_block",
                         &tom::ObservationWeightSettings::setBlock,
                         py::arg( "weight_block" ),
                         R"doc(
         Function to create settings with a single weight matrix for all observations, including correlations between
         different observations.

         Parameters
         ----------
         weight_block : numpy.ndarray[numpy.float64[n, n]]
             Symmetric weight matrix, with one row and column per component of each observation (ordered by observation,
             then by component).

         Returns
         -------
         ObservationWeightSettings
             Weight settings.
      )doc" );

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
             Name of the weighing scheme assigned to the observations (``"VFCC17"``, ``"ADES"`` or
             ``"ADES_NIGHT"``). The weights are computed when an observation collection is created
             from this object, unless weights have been set manually. For ``"ADES"`` and
             ``"ADES_NIGHT"``, the weights are computed from the rmsRA and rmsDec uncertainties
             stored in the metadata. For ``"ADES_NIGHT"``, the weights are additionally divided by
             the number of observations in the same local night (using the same
             definition of a local night as for ``"VFCC17"``).
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
         Return all string-vector ancillary settings, keyed by setting identifier.

         Returns
         -------
         dict[str, list[str]]
             Ancillary setting identifiers and their values.
      )doc" )
            .def( "get_ancillary_settings_double",
                  &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getAncillarySettingsDouble,
                  R"doc(
         Return all floating-point ancillary settings, keyed by setting identifier.

         Returns
         -------
         dict[str, float]
             Ancillary setting identifiers and their values.
      )doc" )
            .def( "get_ancillary_settings_double_vector",
                  &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getAncillarySettingsDoubleVector,
                  R"doc(
         Return all floating-point-vector ancillary settings, keyed by setting identifier.

         Returns
         -------
         dict[str, list[float]]
             Ancillary setting identifiers and their values.
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
            .def( "set_observation_weight_settings",
                  &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::setObservationWeightSettings,
                  py::arg( "observation_weight_settings" ),
                  R"doc(
         Set the weights of the observations.

         The weights are transferred to the observation dataset by
         :func:`~tudatpy.estimation.observations.create_observation_dataset_from_tracking_data`. Weights that were set before
         (also with :meth:`set_observation_weights`) are replaced.

         Parameters
         ----------
         observation_weight_settings : ObservationWeightSettings
             Settings defining the weights of all observations in this object.

         Returns
         -------
         None
      )doc" )
            .def( "get_observation_weight_settings",
                  &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationWeightSettings,
                  py::return_value_policy::copy,
                  R"doc(
         Return the settings defining the weights of the observations.

         Returns
         -------
         ObservationWeightSettings
             Weight settings of the observations. An error is raised if no weights are set.
      )doc" )
            .def_property_readonly( "has_observation_weight_settings",
                                    &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::hasObservationWeightSettings,
                                    R"doc(
         **read-only**

         Whether weights are set for the observations (with :meth:`set_observation_weight_settings` or
         :meth:`set_observation_weights`).

         :type: bool
      )doc" )
            .def( "clear_observation_weight_settings",
                  &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::clearObservationWeightSettings,
                  R"doc(
         Remove the weights of the observations, so that all weights are equal to 1 when converting to an observation dataset.

         Returns
         -------
         None
      )doc" )
            .def( "set_observation_correction_settings",
                  &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::setObservationCorrectionSettings,
                  py::arg( "settings" ),
                  R"doc(
         Set the corrections to the angular position observations that are computed when converting to an observation dataset.

         Only possible for two-component angular position observations with times defined at the receiver. See
         :class:`AngularObservationCorrectionSettings`.

         Parameters
         ----------
         settings : AngularObservationCorrectionSettings
             Settings defining the corrections.

         Returns
         -------
         None
      )doc" )
            .def_property_readonly( "observation_correction_settings",
                                    &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationCorrectionSettings,
                                    R"doc(
         **read-only**

         Settings defining the corrections to the angular position observations (``None`` if not set).

         :type: AngularObservationCorrectionSettings
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
            .def( "add_observation_metadata",
                  &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::addObservationMetadata,
                  py::arg( "key" ),
                  py::arg( "values" ),
                  "Attach one string per observation. Values remain aligned when observations are removed; "
                  "they are not interpreted as simulation ancillary settings." )
            .def( "add_numerical_observation_metadata",
                  &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::addNumericalObservationMetadata,
                  py::arg( "key" ),
                  py::arg( "values" ),
                  R"doc(
         Attach one numerical value per observation under a user-defined name.

         The values remain associated with the correct observations when observations are removed from this object, and are
         transferred to the observation dataset by
         :func:`~tudatpy.estimation.observations.create_observation_dataset_from_tracking_data` (see
         :meth:`~tudatpy.estimation.observations.ObservationDataset.get_numerical_observation_metadata`). They are not used by
         the observation models. Values that were stored under the same name are replaced.

         Parameters
         ----------
         key : str
             Name under which the values are stored.
         values : list[float]
             Value for each observation (one per observation, regardless of the number of observation components).

         Returns
         -------
         None
      )doc" )
            .def( "get_numerical_observation_metadata",
                  py::overload_cast< const std::string& >(
                          &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getNumericalObservationMetadata, py::const_ ),
                  py::arg( "key" ),
                  py::return_value_policy::copy,
                  R"doc(
         Return the numerical values stored under a given name (see :meth:`add_numerical_observation_metadata`).

         Parameters
         ----------
         key : str
             Name under which the values are stored.

         Returns
         -------
         list[float]
             Value for each observation. An ``IndexError`` is raised if no values are stored under ``key``.
      )doc" )
            .def_property_readonly(
                    "numerical_observation_metadata",
                    py::overload_cast<>( &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getNumericalObservationMetadata, py::const_ ),
                    py::return_value_policy::copy,
                    R"doc(
         **read-only**

         All numerical values stored with :meth:`add_numerical_observation_metadata`, with the name as key and the list
         of values (one per observation) as value.

         :type: dict[str, list[float]]
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
            .def( py::init< const tudat::Time, const tudat::Time, const double, const double >( ),
                  py::arg( "start_time" ),
                  py::arg( "end_time" ),
                  py::arg( "start_frequency" ),
                  py::arg( "frequency_rate" ) )
            .def_readwrite( "start_time", &tdat::RampedFrequencySupplementaryData::FrequencyRamp::startTime_, R"doc(
         **read-only**

         Ramp start time, preserving extended epoch precision. Float inputs are accepted.

         :type: tudatpy.astro.time_representation.Time
      )doc" )
            .def_readwrite( "end_time", &tdat::RampedFrequencySupplementaryData::FrequencyRamp::endTime_, R"doc(
         **read-only**

         Ramp end time, preserving extended epoch precision. Float inputs are accepted.

         :type: tudatpy.astro.time_representation.Time
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
                                    py::return_value_policy::copy,
                                    R"doc(
         **read-only**

         Independent copies of the frequency-ramp records. Appending ramps to
         the container does not invalidate previously retrieved records, and
         editing a retrieved record does not change the container.

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
         start_time : tudatpy.astro.time_representation.Time | float
             Ramp start time.
         end_time : tudatpy.astro.time_representation.Time | float
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
         Each entry takes effect at its epoch and is installed as a ramp with
         zero rate. The existing ramp calculator holds the first and last
         frequencies outside the tabulated epochs. Later entries replace earlier
         entries at the same epoch when histories or existing ramp tables are merged.
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

    py::class_< tdat::CameraCalibration >( m, "CameraCalibration", R"doc(
         Calibration of a camera, at a given time (see
         :meth:`CameraInstrumentSupplementaryData.get_calibration`).
      )doc" )
            .def_readonly( "focal_length", &tdat::CameraCalibration::focalLength_, R"doc(Camera focal length.)doc" )
            .def_readonly( "principal_point", &tdat::CameraCalibration::principalPoint_, R"doc(Camera principal point.)doc" )
            .def_readonly( "field_of_view_bounds", &tdat::CameraCalibration::fieldOfViewBounds_, R"doc(Field-of-view bounds.)doc" )
            .def_readonly( "k_matrix", &tdat::CameraCalibration::kMatrix_, R"doc(Camera K matrix.)doc" );

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
            .def( "get_calibration",
                  &tdat::CameraInstrumentSupplementaryData::getCalibration,
                  py::arg( "time" ) = TUDAT_NAN,
                  R"doc(
         Function to get the calibration of the camera at a given time.

         Returns the constant calibration of the camera, or, for a camera
         whose calibration changes between pictures (such as one created
         from SUM files), that of the picture closest in time.

         Parameters
         ----------
         time : float, default=nan
             Time, in UTC seconds since J2000 (the time scale of
             :attr:`rotation_from_inertial_to_camera_frame_history`). Only
             required when the calibration changes between pictures.

         Returns
         -------
         CameraCalibration
             Calibration of the camera at the given time.
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

         Camera-frame rotation history indexed by picture time, in UTC seconds
         since J2000 (as in the source files; the cameras created from this
         data use TDB instead). The rotations start from the inertial frame,
         or from the body-fixed frame of :attr:`pointing_reference_body_name`
         when that is set.

         :type: dict[float, object]
      )doc" )
            .def_property( "pointing_reference_body_name",
                           &tdat::CameraInstrumentSupplementaryData::getPointingReferenceBodyName,
                           &tdat::CameraInstrumentSupplementaryData::setPointingReferenceBodyName,
                           R"doc(
         Body whose body-fixed frame the camera rotation history starts from.
         An empty string means the inertial frame.

         :type: str
      )doc" )
            .def_property( "body_fixed_position",
                           &tdat::CameraInstrumentSupplementaryData::getBodyFixedPosition,
                           &tdat::CameraInstrumentSupplementaryData::setBodyFixedPosition,
                           R"doc(
         Position of the camera in the body-fixed frame of the body that
         carries it [m]. It is used when the camera is created from this
         supplementary data. Defaults to zero.

         :type: numpy.ndarray[numpy.float64[3, 1]]
      )doc" )
            .def_property( "pointing_uncertainty_history",
                           &tdat::CameraInstrumentSupplementaryData::getPointingUncertaintyHistory,
                           &tdat::CameraInstrumentSupplementaryData::setPointingUncertaintyHistory,
                           R"doc(
         A-priori uncertainty of the camera pointing about the three camera
         axes [rad], indexed by the picture times of the rotation history (in
         UTC seconds since J2000). Empty when the source gives no value.

         :type: dict[float, numpy.ndarray[numpy.float64[3, 1]]]
      )doc" );

    py::class_< tdat::ReferencePointSupplementaryData, std::shared_ptr< tdat::ReferencePointSupplementaryData > >(
            m, "ReferencePointSupplementaryData", R"doc(
         Container for the position of a reference point fixed to a body, such
         as a ground station or a surface landmark.

         When applied with
         :func:`~tudatpy.estimation.observations.set_tracking_supplementary_data_in_bodies`,
         a ground station with this position is created on the associated body,
         named after the reference point. If the ground station already exists,
         its position must match.
      )doc" )
            .def( py::init< const Eigen::Vector3d& >( ), py::arg( "body_fixed_position" ) )
            .def_property_readonly( "body_fixed_position",
                                    &tdat::ReferencePointSupplementaryData::getBodyFixedPosition,
                                    R"doc(
         **read-only**

         Cartesian position of the reference point in the body-fixed frame [m].

         :type: numpy.ndarray[numpy.float64[3, 1]]
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
            .def( py::init< const std::map< double, Eigen::Vector6d >&,
                            const std::string&,
                            const bool,
                            const std::string&,
                            const std::string& >( ),
                  py::arg( "state_history" ),
                  py::arg( "frame_origin" ),
                  py::arg( "is_velocity_defined" ),
                  py::arg( "time_scale" ) = "TDB",
                  py::arg( "frame_orientation" ) = "J2000",
                  R"doc(
         Create translational-state supplementary data from a tabulated state history.

         Parameters
         ----------
         state_history : dict[float, numpy.ndarray]
             Cartesian states (position and velocity) per epoch.
         frame_origin : str
             Origin of the frame in which the states are expressed.
         is_velocity_defined : bool
             Whether the velocity entries are valid. If False, velocities are
             computed from positions by finite differences when the data are applied.
         time_scale : str, default "TDB"
             Time scale of the epochs ("TDB" or "UTC").
         frame_orientation : str, default "J2000"
             Orientation of the frame in which the states are expressed.
      )doc" )
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
            .def_property_readonly( "frame_orientation",
                                    &tdat::TranslationalStateSupplementaryData::getFrameOrientation,
                                    R"doc(
         **read-only**

         Orientation of the translational-state frame.

         :type: str
      )doc" )
            .def_property_readonly( "is_velocity_defined",
                                    &tdat::TranslationalStateSupplementaryData::isVelocityDefined,
                                    R"doc(
         **read-only**

         Boolean that defines whether velocity entries are defined, if not they are computed through finite differences when this object is processed to update the environment.

         :type: bool
      )doc" )
            .def_property_readonly( "time_scale",
                                    &tdat::TranslationalStateSupplementaryData::getTimeScale,
                                    R"doc(
         **read-only**

         Time scale of the tabulated state-history epochs.

         :type: str
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

    py::class_< tdat::TrackingSupplementaryData, std::shared_ptr< tdat::TrackingSupplementaryData > >( m,
                                                                                                       "TrackingSupplementaryData",
                                                                                                       R"doc(
         Container for tracking supplementary data.

         These objects store auxiliary information loaded from tracking data
         files and are applied to a system of bodies by
         :func:`~tudatpy.estimation.observations.set_tracking_supplementary_data_in_bodies`.

         Parameters
         ----------
         body_name : str
             Name of the body to which the data apply.
         reference_point_name : str
             Name of the reference point on the body, or an empty string for
             body-level data.
      )doc" )
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
            .def_property( "is_passive_radar_reflector",
                           &tdat::TrackingSupplementaryData::isPassiveRadarReflector,
                           &tdat::TrackingSupplementaryData::setIsPassiveRadarReflector,
                           R"doc(
         Whether the named body is a passive radar reflector. Applying this
         supplementary data sets its turnaround ratio to one for every
         frequency-band pair.

         :type: bool
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
            .def_property( "reference_point_supplementary_data",
                           &tdat::TrackingSupplementaryData::getReferencePointSupplementaryData,
                           &tdat::TrackingSupplementaryData::setReferencePointSupplementaryData,
                           R"doc(
         Body-fixed position of the reference point, or None if not provided.

         :type: ReferencePointSupplementaryData | None
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

    auto landmarksModule = m.def_submodule( "landmarks" );
    landmarks::expose_landmarks( landmarksModule );
}

}  // namespace tracking_data

}  // namespace data_input

}  // namespace tudatpy
