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
            m, "TrackingData", R"doc(Container for source-loaded tracking observations.

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
                  py::arg( "weighing_scheme" ) = "" )
            .def_property_readonly( "observable_type",
                                    &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getObservableType,
                                    R"doc(str: Observable type represented by this data set.)doc" )
            .def_property_readonly( "link_ends",
                                    &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getLinkEnds,
                                    R"doc(list: Link-end definition associated with the observations.)doc" )
            .def_property_readonly( "observations",
                                    &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getObservations,
                                    R"doc(list: Observation vectors.)doc" )
            .def_property_readonly( "epochs",
                                    &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationEpochs,
                                    R"doc(list: Observation epochs.)doc" )
            .def_property_readonly( "reference_link_end",
                                    &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getReferenceLinkEnd,
                                    R"doc(str: Reference link end for the observation epochs.)doc" )
            .def_property_readonly( "time_scale",
                                    &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getTimeScale,
                                    R"doc(str: Time scale of the observation epochs.)doc" )
            .def_property_readonly( "weighing_scheme",
                                    &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getWeighingScheme,
                                    R"doc(str: Name of the weighing scheme assigned to the observations.)doc" )
            .def( "add_string_ancillary_setting",
                  py::overload_cast< const std::string, const std::string >(
                          &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::addAncillarySettings ),
                  py::arg( "ancillary_settings_type" ),
                  py::arg( "ancillary_settings_value" ),
                  R"doc(Add a string ancillary setting.

Parameters
----------
ancillary_settings_type : str
    Ancillary setting identifier.
ancillary_settings_value : str
    Ancillary setting value.

Returns
-------
None
    The setting is added in place.
)doc" )
            .def( "add_string_vector_ancillary_setting",
                  py::overload_cast< const std::string, const std::vector< std::string > >(
                          &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::addAncillarySettings ),
                  py::arg( "ancillary_settings_type" ),
                  py::arg( "ancillary_settings_value" ),
                  R"doc(Add a string-vector ancillary setting.

Parameters
----------
ancillary_settings_type : str
    Ancillary setting identifier.
ancillary_settings_value : list[str]
    Ancillary setting values.

Returns
-------
None
    The setting is added in place.
)doc" )
            .def( "add_double_ancillary_setting",
                  py::overload_cast< const std::string, const double >(
                          &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::addAncillarySettings ),
                  py::arg( "ancillary_settings_type" ),
                  py::arg( "ancillary_settings_value" ),
                  R"doc(Add a floating-point ancillary setting.

Parameters
----------
ancillary_settings_type : str
    Ancillary setting identifier.
ancillary_settings_value : float
    Ancillary setting value.

Returns
-------
None
    The setting is added in place.
)doc" )
            .def( "add_double_vector_ancillary_setting",
                  py::overload_cast< const std::string, const std::vector< double > >(
                          &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::addAncillarySettings ),
                  py::arg( "ancillary_settings_type" ),
                  py::arg( "ancillary_settings_value" ),
                  R"doc(Add a floating-point-vector ancillary setting.

Parameters
----------
ancillary_settings_type : str
    Ancillary setting identifier.
ancillary_settings_value : list[float]
    Ancillary setting values.

Returns
-------
None
    The setting is added in place.
)doc" )
            .def( "get_ancillary_settings_string_vector",
                  &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getAncillarySettingsStringVector,
                  R"doc(Return a string-vector ancillary setting.

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
                  R"doc(Return a floating-point ancillary setting.

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
                  R"doc(Return a floating-point-vector ancillary setting.

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
                  R"doc(Set per-observation weights.

Parameters
----------
observation_weights : list
    Observation weights.

Returns
-------
None
    The weights are set in place.
)doc" )
            .def( "reset_single_observation_weight",
                  &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::setSingleObservationWeight,
                  py::arg( "index" ),
                  py::arg( "observation_weight" ),
                  R"doc(Reset one observation weight.

Parameters
----------
index : int
    Observation index.
observation_weight : float
    Replacement observation weight.

Returns
-------
None
    The weight is set in place.
)doc" )
            .def( "get_observation_weights",
                  &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationWeights,
                  R"doc(Return per-observation weights.

Returns
-------
list
    Observation weights.
)doc" )
            .def( "get_concatenated_observation_weights",
                  &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationWeightsVector,
                  R"doc(Return concatenated observation weights.

Returns
-------
list[float]
    Concatenated observation weights.
)doc" )
            .def( "set_observation_corrections",
                  &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::setObservationCorrections,
                  py::arg( "observation_corrections" ),
                  R"doc(Set per-observation corrections.

Parameters
----------
observation_corrections : list
    Observation correction vectors.

Returns
-------
None
    The corrections are set in place.
)doc" )
            .def( "reset_single_observation_correction",
                  &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::setSingleObservationCorrection,
                  py::arg( "index" ),
                  py::arg( "observation_correction" ),
                  R"doc(Reset one observation correction.

Parameters
----------
index : int
    Observation index.
observation_correction : numpy.ndarray
    Replacement observation correction vector.

Returns
-------
None
    The correction is set in place.
)doc" )
            .def( "get_observation_corrections",
                  &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationCorrections,
                  R"doc(Return per-observation corrections.

Returns
-------
list
    Observation correction vectors.
)doc" )
            .def( "remove_single_observation_entry",
                  &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::removeSingleObservationEntry,
                  py::arg( "index" ),
                  R"doc(Remove one observation entry.

Parameters
----------
index : int
    Observation index to remove.

Returns
-------
None
    The observation entry is removed in place.
)doc" );

    py::class_< tdat::FrequencySupplementaryData, std::shared_ptr< tdat::FrequencySupplementaryData > >(
            m, "FrequencySupplementaryData", R"doc(Base class for frequency supplementary data.)doc" );

    py::class_< tdat::RampedFrequencySupplementaryData,
                std::shared_ptr< tdat::RampedFrequencySupplementaryData >,
                tdat::FrequencySupplementaryData >(
            m, "RampedFrequencySupplementaryData", R"doc(Container for ramped-frequency supplementary data.)doc" )
            .def( py::init<>( ) )
            .def( "add_frequency_ramp",
                  &tdat::RampedFrequencySupplementaryData::addFrequencyRamp,
                  py::arg( "start_time" ),
                  py::arg( "end_time" ),
                  py::arg( "start_frequency" ),
                  py::arg( "frequency_rate" ),
                  R"doc(Add one frequency ramp.

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
    The ramp is added in place.
)doc" );

    py::class_< tdat::TrackingSupplementaryData, std::shared_ptr< tdat::TrackingSupplementaryData > >(
            m, "TrackingSupplementaryData", R"doc(Container for tracking supplementary data.)doc" )
            .def( py::init<>( ) )
            .def( py::init< const std::string&, const std::string& >( ), py::arg( "body_name" ), py::arg( "reference_point_name" ) )
            .def_property_readonly( "body_name", &tdat::TrackingSupplementaryData::getBodyName, R"doc(str: Body name.)doc" )
            .def_property_readonly(
                    "reference_point_name", &tdat::TrackingSupplementaryData::getReferencePointName, R"doc(str: Reference point name.)doc" )
            .def( "set_frequency_supplementary_data",
                  &tdat::TrackingSupplementaryData::setFrequencySupplementaryData,
                  py::arg( "frequency_supplementary_data" ),
                  R"doc(Set frequency supplementary data.

Parameters
----------
frequency_supplementary_data : list[FrequencySupplementaryData]
    Frequency supplementary data objects.

Returns
-------
None
    The data are set in place.
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
