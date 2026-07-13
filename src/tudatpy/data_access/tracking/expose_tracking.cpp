#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif

#include "expose_tracking.h"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "expose_fdets.h"
#include "expose_generic_text_file.h"
#include "expose_ifms.h"
#include "expose_odf.h"
#include "expose_psf.h"
#include "scalarTypes.h"
#include "tudat/io/trackingData.h"
#include "tudat/io/trackingSupplementaryData.h"
#include "tudat/simulation/estimation_setup/createObservationCollection.h"

namespace py = pybind11;
namespace tdat = tudat::data;
namespace tom = tudat::observation_models;
namespace tss = tudat::simulation_setup;

namespace tudatpy
{

namespace data_access
{

namespace tracking
{

void expose_tracking( py::module& m )
{
    py::class_< tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >, std::shared_ptr< tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE > > >(
            m, "TrackingData", R"doc(Container for source-loaded tracking observations.)doc" )
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
            .def_property_readonly( "observable_type", &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getObservableType )
            .def_property_readonly( "link_ends", &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getLinkEnds )
            .def_property_readonly( "observations", &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getObservations )
            .def_property_readonly( "epochs", &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationEpochs )
            .def_property_readonly( "reference_link_end", &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getReferenceLinkEnd )
            .def_property_readonly( "time_scale", &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getTimeScale )
            .def_property_readonly( "weighing_scheme", &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getWeighingScheme )
            .def( "add_ancillary_settings",
                  py::overload_cast< const std::string, const std::vector< std::string > >(
                          &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::addAncillarySettings ),
                  py::arg( "ancillary_settings_type" ),
                  py::arg( "ancillary_settings_value" ) )
            .def( "add_ancillary_settings",
                  py::overload_cast< const std::string, const double >(
                          &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::addAncillarySettings ),
                  py::arg( "ancillary_settings_type" ),
                  py::arg( "ancillary_settings_value" ) )
            .def( "add_ancillary_settings",
                  py::overload_cast< const std::string, const std::vector< double > >(
                          &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::addAncillarySettings ),
                  py::arg( "ancillary_settings_type" ),
                  py::arg( "ancillary_settings_value" ) )
            .def( "get_ancillary_settings_string_vector",
                  &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getAncillarySettingsStringVector )
            .def( "get_ancillary_settings_double", &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getAncillarySettingsDouble )
            .def( "get_ancillary_settings_double_vector",
                  &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getAncillarySettingsDoubleVector )
            .def( "set_observation_weights",
                  &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::setObservationWeights,
                  py::arg( "observation_weights" ) )
            .def( "reset_single_observation_weight",
                  &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::setSingleObservationWeight,
                  py::arg( "index" ),
                  py::arg( "observation_weight" ) )
            .def( "get_observation_weights", &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationWeights )
            .def( "get_concatenated_observation_weights", &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationWeightsVector )
            .def( "set_observation_corrections",
                  &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::setObservationCorrections,
                  py::arg( "observation_corrections" ) )
            .def( "reset_single_observation_correction",
                  &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::setSingleObservationCorrection,
                  py::arg( "index" ),
                  py::arg( "observation_correction" ) )
            .def( "get_observation_corrections", &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationCorrections )
            .def( "remove_single_observation_entry",
                  &tdat::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::removeSingleObservationEntry,
                  py::arg( "index" ) );

    py::class_< tdat::FrequencySupplementaryData, std::shared_ptr< tdat::FrequencySupplementaryData > >( m, "FrequencySupplementaryData" );

    py::class_< tdat::RampedFrequencySupplementaryData,
                std::shared_ptr< tdat::RampedFrequencySupplementaryData >,
                tdat::FrequencySupplementaryData >( m, "RampedFrequencySupplementaryData" )
            .def( py::init<>( ) )
            .def( "add_frequency_ramp",
                  &tdat::RampedFrequencySupplementaryData::addFrequencyRamp,
                  py::arg( "start_time" ),
                  py::arg( "end_time" ),
                  py::arg( "start_frequency" ),
                  py::arg( "frequency_rate" ) );

    py::class_< tdat::TrackingSupplementaryData, std::shared_ptr< tdat::TrackingSupplementaryData > >( m, "TrackingSupplementaryData" )
            .def( py::init<>( ) )
            .def( py::init< const std::string&, const std::string& >( ), py::arg( "body_name" ), py::arg( "reference_point_name" ) )
            .def_property_readonly( "body_name", &tdat::TrackingSupplementaryData::getBodyName )
            .def_property_readonly( "reference_point_name", &tdat::TrackingSupplementaryData::getReferencePointName )
            .def( "set_frequency_supplementary_data",
                  &tdat::TrackingSupplementaryData::setFrequencySupplementaryData,
                  py::arg( "frequency_supplementary_data" ) );

    m.def( "create_observation_collection",
           &tom::createObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "tracking_data" ),
           py::arg( "bodies" ) );

    m.def( "set_tracking_supplementary_data_in_bodies",
           py::overload_cast< tss::SystemOfBodies&, const std::vector< std::shared_ptr< tdat::TrackingSupplementaryData > >& >(
                   &tom::setTrackingSupplementaryDataInBodies ),
           py::arg( "bodies" ),
           py::arg( "supplementary_data" ) );

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
}

}  // namespace tracking

}  // namespace data_access

}  // namespace tudatpy
