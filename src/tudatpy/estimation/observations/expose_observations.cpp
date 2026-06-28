/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif
#include "expose_observations.h"

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "scalarTypes.h"
#include "tudat/simulation/estimation_setup/observationDataset.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"
#include "observations_processing/expose_observations_processing.h"
#include "observations_geometry/expose_observations_geometry.h"

namespace py = pybind11;
namespace tss = tudat::simulation_setup;
namespace tom = tudat::observation_models;
namespace te = tudat::ephemerides;

namespace tudat
{

namespace simulation_setup
{

template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< tom::SingleObservationSet< ObservationScalarType, TimeType > > singleObservationSetWithoutDependentVariables(
        const tom::ObservableType observableType,
        const tom::LinkDefinition& linkEnds,
        const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
        const std::vector< TimeType > observationTimes,
        const tom::LinkEndType referenceLinkEnd,
        const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillarySettings = nullptr )
{
    std::cerr << "Function single_observation_set is deprecated. Use create_single_observation_set instead" << std::endl;
    return std::make_shared< tom::SingleObservationSet< ObservationScalarType, TimeType > >( observableType,
                                                                                             linkEnds,
                                                                                             observations,
                                                                                             observationTimes,
                                                                                             referenceLinkEnd,
                                                                                             std::vector< Eigen::VectorXd >( ),
                                                                                             nullptr,
                                                                                             ancillarySettings );
}

}  // namespace simulation_setup

}  // namespace tudat

namespace tudatpy
{
namespace estimation
{
namespace observations
{

void expose_observations( py::module& m )
{
    auto observations_processing = m.def_submodule( "observations_processing" );
    observations_processing::expose_observations_processing( observations_processing );

    auto observations_geometry = m.def_submodule( "observations_geometry" );
    observations_geometry::expose_observations_geometry( observations_geometry );

    // OBSERVATION DATASET

    py::class_< tom::ObservationSetMetadata< STATE_SCALAR_TYPE, TIME_TYPE > >( m,
                                                                               "ObservationSetMetadata",
                                                                               R"doc(
Metadata describing one logical observation set stored in an :class:`ObservationDataset`.

The metadata identifies the observable type, link definition, reference link end,
observable size and the registered ancillary/dependent-variable layouts used by
all observations in the set.
)doc" )
            .def_readonly( "observable_type",
                           &tom::ObservationSetMetadata< STATE_SCALAR_TYPE, TIME_TYPE >::observableType_,
                           R"doc(Observable type stored in this set.)doc" )
            .def_readonly( "link_definition_id",
                           &tom::ObservationSetMetadata< STATE_SCALAR_TYPE, TIME_TYPE >::linkDefinitionId_,
                           R"doc(Index of the link definition in the dataset link-definition registry.)doc" )
            .def_readonly( "reference_link_end",
                           &tom::ObservationSetMetadata< STATE_SCALAR_TYPE, TIME_TYPE >::referenceLinkEnd_,
                           R"doc(Reference link end used for all observations in this set.)doc" )
            .def_readonly( "observable_size",
                           &tom::ObservationSetMetadata< STATE_SCALAR_TYPE, TIME_TYPE >::observableSize_,
                           R"doc(Number of scalar components in one observation of this set.)doc" )
            .def_readonly( "ancillary_settings_id",
                           &tom::ObservationSetMetadata< STATE_SCALAR_TYPE, TIME_TYPE >::ancillarySettingsId_,
                           R"doc(Index of the ancillary settings in the dataset registry.)doc" )
            .def_readonly( "dependent_variable_layout_id",
                           &tom::ObservationSetMetadata< STATE_SCALAR_TYPE, TIME_TYPE >::dependentVariableLayoutId_,
                           R"doc(Index of the dependent-variable bookkeeping layout in the dataset registry.)doc" );

    py::class_< tom::ObservationDatasetRow< TIME_TYPE > >( m,
                                                           "ObservationDatasetRow",
                                                           R"doc(
Row-level storage metadata for one observation inside an :class:`ObservationDataset`.

Each row points to the first scalar component of the observation and records the
observation time, owning set and index within that set.
)doc" )
            .def_readonly( "time", &tom::ObservationDatasetRow< TIME_TYPE >::time_, R"doc(Observation time.)doc" )
            .def_readonly(
                    "set_id", &tom::ObservationDatasetRow< TIME_TYPE >::setId_, R"doc(Identifier of the owning observation set.)doc" )
            .def_readonly( "first_scalar_component",
                           &tom::ObservationDatasetRow< TIME_TYPE >::firstScalarComponent_,
                           R"doc(Index of the first scalar component belonging to this observation.)doc" )
            .def_readonly( "scalar_size",
                           &tom::ObservationDatasetRow< TIME_TYPE >::scalarSize_,
                           R"doc(Number of scalar components in this observation.)doc" )
            .def_readonly( "index_in_set",
                           &tom::ObservationDatasetRow< TIME_TYPE >::indexInSet_,
                           R"doc(Index of this observation within its observation set.)doc" )
            .def_readonly( "is_active",
                           &tom::ObservationDatasetRow< TIME_TYPE >::isActive_,
                           R"doc(Boolean flag indicating whether this row is active in projections.)doc" )
            .def_readonly( "rejection_reason",
                           &tom::ObservationDatasetRow< TIME_TYPE >::rejectionReason_,
                           R"doc(Optional text describing why this observation was rejected.)doc" );

    py::class_< tom::ObservationScalarComponentRow >( m,
                                                      "ObservationScalarComponentRow",
                                                      R"doc(
Storage metadata for one scalar component of an observation.

The row records the owning observation and component index inside that
observation.
)doc" )
            .def_readonly( "observation_id",
                           &tom::ObservationScalarComponentRow::observationId_,
                           R"doc(Identifier of the owning observation.)doc" )
            .def_readonly( "component_index",
                           &tom::ObservationScalarComponentRow::componentIndex_,
                           R"doc(Component index within the owning vector-valued observation.)doc" );

    py::class_< tom::ObservationWeightBlock >( m,
                                               "ObservationWeightBlock",
                                               R"doc(
Advanced dense weight block over selected scalar components.

This type is intended for rare off-diagonal correlations that are not naturally
represented as per-observation weights or as a full set-level weight block.
)doc" )
            .def( py::init<>( ), R"doc(Create an empty observation weight block.)doc" )
            .def_readwrite( "row_scalar_component_ids",
                            &tom::ObservationWeightBlock::rowScalarComponentIds_,
                            R"doc(Scalar component ids corresponding to the block rows.)doc" )
            .def_readwrite( "column_scalar_component_ids",
                            &tom::ObservationWeightBlock::columnScalarComponentIds_,
                            R"doc(Scalar component ids corresponding to the block columns.)doc" )
            .def_readwrite( "weight_block",
                            &tom::ObservationWeightBlock::weightBlock_,
                            R"doc(Dense weight block for the selected scalar components.)doc" );

    py::class_< tom::EstimationVectorProjection< STATE_SCALAR_TYPE, TIME_TYPE > >( m,
                                                                                   "EstimationVectorProjection",
                                                                                   R"doc(
Flat estimation-vector view of an :class:`ObservationDataset`.

This object contains the concatenated observation, residual and weight vectors,
together with the scalar-component provenance needed to map each entry back to
a dataset observation row and observation set. The diagonal weights are always
available through :attr:`weight_vector`. The full matrix is returned as a sparse
matrix and is only needed when off-diagonal terms are present.
)doc" )
            .def_property_readonly( "observation_vector",
                                    &tom::EstimationVectorProjection< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationVector,
                                    R"doc(Concatenated vector of observed values.)doc" )
            .def_property_readonly( "residual_vector",
                                    &tom::EstimationVectorProjection< STATE_SCALAR_TYPE, TIME_TYPE >::getResidualVector,
                                    R"doc(Concatenated vector of residual values.)doc" )
            .def_property_readonly( "weight_vector",
                                    &tom::EstimationVectorProjection< STATE_SCALAR_TYPE, TIME_TYPE >::getWeightVector,
                                    R"doc(Concatenated vector of scalar observation weights.)doc" )
            .def_property_readonly( "weight_matrix",
                                    &tom::EstimationVectorProjection< STATE_SCALAR_TYPE, TIME_TYPE >::getWeightMatrix,
                                    R"doc(
Sparse weight matrix in the same order as :attr:`observation_vector`.

For diagonal-only weights this matrix is generated from :attr:`weight_vector`.
For off-diagonal weights it contains the materialized sparse matrix assembled
from per-observation blocks, set-level blocks and advanced scalar-component
blocks.
)doc" )
            .def_property_readonly( "is_diagonal_weight_only",
                                    &tom::EstimationVectorProjection< STATE_SCALAR_TYPE, TIME_TYPE >::isDiagonalWeightOnly,
                                    R"doc(True when the projection weight matrix contains no off-diagonal entries.)doc" )
            .def_property_readonly( "has_off_diagonal_weights",
                                    &tom::EstimationVectorProjection< STATE_SCALAR_TYPE, TIME_TYPE >::hasOffDiagonalWeights,
                                    R"doc(True when the projection weight matrix contains off-diagonal entries.)doc" )
            .def_property_readonly( "times",
                                    &tom::EstimationVectorProjection< STATE_SCALAR_TYPE, TIME_TYPE >::getTimes,
                                    R"doc(Observation time associated with each scalar component.)doc" )
            .def_property_readonly( "observation_ids",
                                    &tom::EstimationVectorProjection< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationIds,
                                    R"doc(Observation row identifier associated with each scalar component.)doc" )
            .def_property_readonly( "set_ids",
                                    &tom::EstimationVectorProjection< STATE_SCALAR_TYPE, TIME_TYPE >::getSetIds,
                                    R"doc(Observation set identifier associated with each scalar component.)doc" )
            .def_property_readonly( "scalar_component_ids",
                                    &tom::EstimationVectorProjection< STATE_SCALAR_TYPE, TIME_TYPE >::getScalarComponentIds,
                                    R"doc(Scalar-component row identifier for each entry in the projection.)doc" );

    py::class_< tom::ObservationCondition< STATE_SCALAR_TYPE, TIME_TYPE > >( m,
                                                                             "ObservationCondition",
                                                                             R"doc(
Composable row-level condition used to select observations in an ObservationDataset.

Conditions operate on individual observation rows, not complete observation
sets. Combine conditions with ``&`` and ``|`` and negate them with ``~``.
)doc" )
            .def( py::init<>( ), R"doc(Create a condition that selects all observations.)doc" )
            .def_static( "all",
                         &tom::ObservationCondition< STATE_SCALAR_TYPE, TIME_TYPE >::all,
                         R"doc(Return a condition that selects all observations.)doc" )
            .def_static( "observable_type",
                         &tom::ObservationCondition< STATE_SCALAR_TYPE, TIME_TYPE >::observableType,
                         py::arg( "observable_type" ),
                         R"doc(Return a condition selecting observations of one observable type.)doc" )
            .def_static( "link_definition",
                         &tom::ObservationCondition< STATE_SCALAR_TYPE, TIME_TYPE >::linkDefinition,
                         py::arg( "link_definition" ),
                         R"doc(Return a condition selecting observations with a matching link definition.)doc" )
            .def_static( "link_end_type",
                         &tom::ObservationCondition< STATE_SCALAR_TYPE, TIME_TYPE >::linkEndType,
                         py::arg( "link_end_type" ),
                         R"doc(Return a condition selecting observations whose link definition contains a link-end type.)doc" )
            .def_static( "link_end",
                         &tom::ObservationCondition< STATE_SCALAR_TYPE, TIME_TYPE >::linkEnd,
                         py::arg( "link_end_type" ),
                         py::arg( "link_end_id" ),
                         R"doc(Return a condition selecting observations with a specific link end.)doc" )
            .def_static( "time_bounds",
                         &tom::ObservationCondition< STATE_SCALAR_TYPE, TIME_TYPE >::timeBounds,
                         py::arg( "start_time" ),
                         py::arg( "end_time" ),
                         R"doc(Return a condition selecting rows with start_time <= observation_time <= end_time.)doc" )
            .def_static( "active",
                         &tom::ObservationCondition< STATE_SCALAR_TYPE, TIME_TYPE >::active,
                         R"doc(Return a condition selecting active observations.)doc" )
            .def_static( "rejected",
                         &tom::ObservationCondition< STATE_SCALAR_TYPE, TIME_TYPE >::rejected,
                         R"doc(Return a condition selecting rejected observations.)doc" )
            .def_static( "residual_absolute_value_greater_than",
                         &tom::ObservationCondition< STATE_SCALAR_TYPE, TIME_TYPE >::residualAbsoluteValueGreaterThan,
                         py::arg( "limit" ),
                         R"doc(Return a condition selecting rows where any residual component exceeds the supplied absolute limit.)doc" )
            .def_static( "observation_absolute_value_greater_than",
                         &tom::ObservationCondition< STATE_SCALAR_TYPE, TIME_TYPE >::observationAbsoluteValueGreaterThan,
                         py::arg( "limit" ),
                         R"doc(Return a condition selecting rows where any observation component exceeds the supplied absolute limit.)doc" )
            .def_static( "dependent_variable_greater_than",
                         &tom::ObservationCondition< STATE_SCALAR_TYPE, TIME_TYPE >::dependentVariableGreaterThan,
                         py::arg( "dependent_variable_settings" ),
                         py::arg( "limit" ),
                         py::arg( "return_first_compatible_settings" ) = false,
                         R"doc(Return a condition selecting rows where a compatible dependent-variable component exceeds the limit.)doc" )
            .def(
                    "__and__",
                    []( const tom::ObservationCondition< STATE_SCALAR_TYPE, TIME_TYPE >& lhs,
                        const tom::ObservationCondition< STATE_SCALAR_TYPE, TIME_TYPE >& rhs ) { return lhs && rhs; },
                    py::is_operator( ),
                    R"doc(Return the logical AND of two conditions.)doc" )
            .def(
                    "__or__",
                    []( const tom::ObservationCondition< STATE_SCALAR_TYPE, TIME_TYPE >& lhs,
                        const tom::ObservationCondition< STATE_SCALAR_TYPE, TIME_TYPE >& rhs ) { return lhs || rhs; },
                    py::is_operator( ),
                    R"doc(Return the logical OR of two conditions.)doc" )
            .def(
                    "__invert__",
                    []( const tom::ObservationCondition< STATE_SCALAR_TYPE, TIME_TYPE >& condition ) { return !condition; },
                    py::is_operator( ),
                    R"doc(Return the logical negation of a condition.)doc" );

    py::class_< tom::ObservationDatasetViewer< STATE_SCALAR_TYPE, TIME_TYPE > >( m,
                                                                                 "ObservationDatasetViewer",
                                                                                 R"doc(
Read-only view on a selected subset of an ObservationDataset.

The viewer stores observation row identifiers selected from a parent dataset and
exposes only inspection and projection methods. It is invalidated if the parent
dataset is structurally modified.
)doc" )
            .def_property_readonly( "number_of_observations",
                                    &tom::ObservationDatasetViewer< STATE_SCALAR_TYPE, TIME_TYPE >::getNumberOfObservations,
                                    R"doc(Number of selected observation rows.)doc" )
            .def_property_readonly( "observation_ids",
                                    &tom::ObservationDatasetViewer< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationIds,
                                    R"doc(Selected observation row identifiers.)doc" )
            .def( "observation_row",
                  &tom::ObservationDatasetViewer< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationRow,
                  py::arg( "viewer_index" ),
                  R"doc(Return row metadata for one selected observation.)doc" )
            .def( "observation_value",
                  &tom::ObservationDatasetViewer< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationValue,
                  py::arg( "viewer_index" ),
                  R"doc(Return the vector-valued observation at the selected viewer index.)doc" )
            .def( "observation_time",
                  &tom::ObservationDatasetViewer< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationTime,
                  py::arg( "viewer_index" ),
                  R"doc(Return the observation time at the selected viewer index.)doc" )
            .def( "create_viewer",
                  &tom::ObservationDatasetViewer< STATE_SCALAR_TYPE, TIME_TYPE >::createViewer,
                  py::keep_alive< 0, 1 >( ),
                  py::arg( "condition" ),
                  R"doc(Return a narrower read-only viewer selected from this viewer.)doc" )
            .def( "estimation_vector_projection",
                  &tom::ObservationDatasetViewer< STATE_SCALAR_TYPE, TIME_TYPE >::createEstimationProjection,
                  py::arg( "include_rejected" ) = false,
                  R"doc(Return a projection for this viewer. Rejected observations are excluded by default.)doc" )
            .def( "legacy_vector_projection",
                  &tom::ObservationDatasetViewer< STATE_SCALAR_TYPE, TIME_TYPE >::createLegacyProjection,
                  py::arg( "include_inactive" ) = true,
                  R"doc(Return a legacy projection for this viewer.)doc" );

    py::class_< tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE > > >( m,
                                                                                              "ObservationDataset",
                                                                                              R"doc(
Backend storage object for observation data.

An ``ObservationDataset`` stores observation values, residuals, weights,
dependent variables, set metadata and link/ancillary registries in a single
dataset-centric representation. It is the new backend used by legacy
``SingleObservationSet`` and ``ObservationCollection`` compatibility wrappers.
)doc" )
            .def( py::init<>( ), R"doc(Create an empty observation dataset.)doc" )
            .def( "add_observation_set",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::addObservationSet,
                  py::arg( "observable_type" ),
                  py::arg( "link_definition" ),
                  py::arg( "observations" ),
                  py::arg( "times" ),
                  py::arg( "reference_link_end" ),
                  py::arg( "dependent_variables" ) = std::vector< Eigen::VectorXd >( ),
                  py::arg( "dependent_variable_bookkeeping" ) = nullptr,
                  py::arg( "ancillary_settings" ) = nullptr,
                  py::arg( "weights" ) = std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >( ),
                  py::arg( "residuals" ) = std::vector< Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >( ),
                  py::arg( "sort_observations" ) = false,
                  py::arg( "erase_duplicate_observations" ) = false,
                  R"doc(
Add a logical observation set and return its dataset set identifier.

Parameters are equivalent to the legacy ``SingleObservationSet`` constructor,
but the data are inserted directly into the dataset backend.
)doc" )
            .def( "add_observation_set_from_dataset",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::addObservationSetFromDataset,
                  py::arg( "source_dataset" ),
                  py::arg( "source_set_id" ),
                  R"doc(Copy one observation set from another dataset and return the new set identifier.)doc" )
            .def( "add_observation_set_with_scalar_weight",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::addObservationSetWithScalarWeight,
                  py::arg( "observable_type" ),
                  py::arg( "link_definition" ),
                  py::arg( "observations" ),
                  py::arg( "times" ),
                  py::arg( "reference_link_end" ),
                  py::arg( "weight" ),
                  py::arg( "dependent_variables" ) = std::vector< Eigen::VectorXd >( ),
                  py::arg( "dependent_variable_bookkeeping" ) = nullptr,
                  py::arg( "ancillary_settings" ) = nullptr,
                  py::arg( "residuals" ) = std::vector< Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >( ),
                  R"doc(Add a set using one compact scalar weight for every observation.)doc" )
            .def( "add_observation_set_with_scalar_weights",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::addObservationSetWithScalarWeights,
                  py::arg( "observable_type" ),
                  py::arg( "link_definition" ),
                  py::arg( "observations" ),
                  py::arg( "times" ),
                  py::arg( "reference_link_end" ),
                  py::arg( "weights" ),
                  py::arg( "dependent_variables" ) = std::vector< Eigen::VectorXd >( ),
                  py::arg( "dependent_variable_bookkeeping" ) = nullptr,
                  py::arg( "ancillary_settings" ) = nullptr,
                  py::arg( "residuals" ) = std::vector< Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >( ),
                  R"doc(Add a set using one compact scalar weight per observation.)doc" )
            .def( "add_observation_set_with_weight_block",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::addObservationSetWithWeightBlock,
                  py::arg( "observable_type" ),
                  py::arg( "link_definition" ),
                  py::arg( "observations" ),
                  py::arg( "times" ),
                  py::arg( "reference_link_end" ),
                  py::arg( "weight_block" ),
                  py::arg( "dependent_variables" ) = std::vector< Eigen::VectorXd >( ),
                  py::arg( "dependent_variable_bookkeeping" ) = nullptr,
                  py::arg( "ancillary_settings" ) = nullptr,
                  py::arg( "residuals" ) = std::vector< Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >( ),
                  R"doc(Add a set using one observable-size weight block for every observation.)doc" )
            .def( "add_observation_set_with_weight_blocks",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::addObservationSetWithWeightBlocks,
                  py::arg( "observable_type" ),
                  py::arg( "link_definition" ),
                  py::arg( "observations" ),
                  py::arg( "times" ),
                  py::arg( "reference_link_end" ),
                  py::arg( "weight_blocks" ),
                  py::arg( "dependent_variables" ) = std::vector< Eigen::VectorXd >( ),
                  py::arg( "dependent_variable_bookkeeping" ) = nullptr,
                  py::arg( "ancillary_settings" ) = nullptr,
                  py::arg( "residuals" ) = std::vector< Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >( ),
                  R"doc(Add a set using one observable-size weight block per observation.)doc" )
            .def( "add_observation_set_with_set_weight_block",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::addObservationSetWithSetWeightBlock,
                  py::arg( "observable_type" ),
                  py::arg( "link_definition" ),
                  py::arg( "observations" ),
                  py::arg( "times" ),
                  py::arg( "reference_link_end" ),
                  py::arg( "set_weight_block" ),
                  py::arg( "dependent_variables" ) = std::vector< Eigen::VectorXd >( ),
                  py::arg( "dependent_variable_bookkeeping" ) = nullptr,
                  py::arg( "ancillary_settings" ) = nullptr,
                  py::arg( "residuals" ) = std::vector< Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >( ),
                  R"doc(Add a set using one full M x M set-level weight block.)doc" )
            .def( "set_link_definition",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::setLinkDefinition,
                  py::arg( "set_id" ),
                  py::arg( "link_definition" ),
                  R"doc(Replace the link definition associated with an observation set.)doc" )
            .def( "set_observations_for_set",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::setObservationsForSet,
                  py::arg( "set_id" ),
                  py::arg( "observations" ),
                  R"doc(Replace all observation vectors in a set.)doc" )
            .def( "set_observation_vector_for_set",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::setObservationVectorForSet,
                  py::arg( "set_id" ),
                  py::arg( "observation_vector" ),
                  R"doc(Replace all observations in a set from one concatenated vector.)doc" )
            .def( "set_residuals_for_set",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::setResidualsForSet,
                  py::arg( "set_id" ),
                  py::arg( "residuals" ),
                  R"doc(Replace all residual vectors in a set.)doc" )
            .def( "set_residual_vector_for_set",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::setResidualVectorForSet,
                  py::arg( "set_id" ),
                  py::arg( "residual_vector" ),
                  R"doc(Replace all residuals in a set from one concatenated vector.)doc" )
            .def( "set_constant_weight_for_set",
                  py::overload_cast< const tom::ObservationSetId, const double >(
                          &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantWeightForSet ),
                  py::arg( "set_id" ),
                  py::arg( "weight" ),
                  R"doc(Set the same scalar weight for every scalar component in an observation set.)doc" )
            .def( "set_constant_weight_for_set",
                  py::overload_cast< const tom::ObservationSetId, const Eigen::VectorXd& >(
                          &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantWeightForSet ),
                  py::arg( "set_id" ),
                  py::arg( "weight" ),
                  R"doc(Set the same component-wise weight vector for every observation in a set.)doc" )
            .def( "set_weight_vector_for_set",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::setWeightVectorForSet,
                  py::arg( "set_id" ),
                  py::arg( "weight_vector" ),
                  R"doc(Replace all weights in a set from one concatenated vector.)doc" )
            .def( "set_weight_matrix_for_set",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::setWeightMatrixForSet,
                  py::arg( "set_id" ),
                  py::arg( "weight_matrix" ),
                  R"doc(
Store one full set-level weight matrix for an observation set.

The matrix must have size ``M x M``, where ``M`` is the total number of scalar
components in the set. For example, three angular-position observations require
a ``6 x 6`` matrix. A set-level matrix is used as the complete weight block for
that set in estimation projections and takes precedence over per-observation
blocks in that set.
)doc" )
            .def( "has_weight_matrix_for_set",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::hasWeightMatrixForSet,
                  py::arg( "set_id" ),
                  R"doc(Return whether a set has an explicitly stored set-level weight matrix.)doc" )
            .def( "set_weight_matrix_for_observation",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::setWeightMatrixForObservation,
                  py::arg( "observation_id" ),
                  py::arg( "weight_matrix" ),
                  R"doc(
Store one observable-size weight matrix for a single observation row.

The matrix must have size ``N x N``, where ``N`` is the scalar size of the
observable. This is the convenient interface for correlations between
components of one vector-valued observation, such as the two components of an
angular-position observation.
)doc" )
            .def( "has_weight_matrix_for_observation",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::hasWeightMatrixForObservation,
                  py::arg( "observation_id" ),
                  R"doc(Return whether an observation row has an explicitly stored observable-size weight matrix.)doc" )
            .def( "add_extra_weight_block",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::addExtraWeightBlock,
                  py::arg( "weight_block" ),
                  R"doc(
Add an advanced dense weight block over selected scalar components.

This method is intended for low-level use when a correlation cannot be expressed
as a per-observation block or as one full set-level block. The row and column
scalar-component ids define where the dense block is inserted in estimation
projections.
)doc" )
            .def( "set_weight_block",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::setWeightBlock,
                  py::arg( "row_observation_ids" ),
                  py::arg( "column_observation_ids" ),
                  py::arg( "weight_block" ),
                  py::arg( "row_components" ) = std::vector< unsigned int >( ),
                  py::arg( "column_components" ) = std::vector< unsigned int >( ),
                  py::arg( "make_symmetric" ) = false,
                  R"doc(
Store a dense weight block selected by observation ids.

Observation ids can be obtained with methods such as
:func:`observation_ids_for_set`, :func:`observation_ids_matching_condition`, or
from an :class:`ObservationDatasetViewer`. Empty component lists select all
scalar components of each selected observation. Non-empty component lists are
applied to every observation in the corresponding row or column selection.

The matrix size must match the expanded scalar selections:
``weight_block.shape == (number_of_selected_row_components,
number_of_selected_column_components)``.

If ``make_symmetric`` is ``True`` and the row and column selections differ, the
transposed block is stored automatically as well. If the selections are
identical, the supplied block must itself be symmetric.
)doc" )
            .def( "extra_weight_blocks",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getExtraWeightBlocks,
                  R"doc(Return the advanced scalar-component weight blocks stored on this dataset.)doc" )
            .def_property_readonly( "has_extra_weight_blocks",
                                    &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::hasExtraWeightBlocks,
                                    R"doc(True when the dataset stores advanced scalar-component weight blocks.)doc" )
            .def( "replace_observation_set_data",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::replaceObservationSetData,
                  py::arg( "set_id" ),
                  py::arg( "observations" ),
                  py::arg( "times" ),
                  py::arg( "dependent_variables" ) = std::vector< Eigen::VectorXd >( ),
                  py::arg( "weights" ) = std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >( ),
                  py::arg( "residuals" ) = std::vector< Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >( ),
                  R"doc(Replace all row-level data for an existing observation set.)doc" )
            .def( "add_observations_to_set",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::addObservationsToSet,
                  py::arg( "set_id" ),
                  py::arg( "observations" ),
                  py::arg( "times" ),
                  py::arg( "dependent_variables" ) = std::vector< Eigen::VectorXd >( ),
                  py::arg( "weights" ) = std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >( ),
                  py::arg( "residuals" ) = std::vector< Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >( ),
                  py::arg( "sort_observations" ) = true,
                  R"doc(Append observations to an existing set.)doc" )
            .def( "remove_observations_from_set",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::removeObservationsFromSet,
                  py::arg( "set_id" ),
                  py::arg( "indices_to_remove" ),
                  R"doc(Remove observations from a set by index within that set.)doc" )
            .def( "remove_observation_from_set",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::removeObservationFromSet,
                  py::arg( "set_id" ),
                  py::arg( "index_to_remove" ),
                  R"doc(Remove one observation from a set by index within that set.)doc" )
            .def( "filtered_observation_indices",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getFilteredObservationIndices,
                  py::arg( "set_id" ),
                  py::arg( "observation_filter" ),
                  R"doc(Return the indices selected by an observation filter for a given set.)doc" )
            .def(
                    "move_observations_to_set",
                    []( tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >& dataset,
                        const tom::ObservationSetId sourceSetId,
                        const std::shared_ptr< tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE > >& targetDataset,
                        const tom::ObservationSetId targetSetId,
                        const std::vector< unsigned int >& indices,
                        const bool removeFromSource ) {
                        if( targetDataset == nullptr )
                        {
                            throw std::runtime_error( "Error when moving observations, target dataset is None." );
                        }
                        dataset.moveObservationsToSet( sourceSetId, *targetDataset, targetSetId, indices, removeFromSource );
                    },
                    py::arg( "source_set_id" ),
                    py::arg( "target_dataset" ),
                    py::arg( "target_set_id" ),
                    py::arg( "indices" ),
                    py::arg( "remove_from_source" ) = true,
                    R"doc(Move or copy selected observations from one dataset set to another.)doc" )
            .def( "time_bounds_for_set",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getTimeBoundsForSet,
                  py::arg( "set_id" ),
                  R"doc(Return the minimum and maximum observation time in a set.)doc" )
            .def( "erase_duplicate_observations_from_set",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::eraseDuplicateObservationsFromSet,
                  py::arg( "set_id" ),
                  py::arg( "print_warning" ) = true,
                  R"doc(Remove duplicate consecutive observation epochs from a set.)doc" )
            .def( "computed_observations_for_set",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getComputedObservationsForSet,
                  py::arg( "set_id" ),
                  R"doc(Return computed observations, defined as observed values minus residuals.)doc" )
            .def( "computed_observation_vector_for_set",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getComputedObservationVectorForSet,
                  py::arg( "set_id" ),
                  R"doc(Return computed observations for a set as one concatenated vector.)doc" )
            .def( "rms_residuals_for_set",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getRmsResidualsForSet,
                  py::arg( "set_id" ),
                  R"doc(Return root-mean-square residuals per scalar component for a set.)doc" )
            .def( "mean_residuals_for_set",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getMeanResidualsForSet,
                  py::arg( "set_id" ),
                  R"doc(Return mean residuals per scalar component for a set.)doc" )
            .def_property_readonly( "number_of_observation_sets",
                                    &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getNumberOfObservationSets,
                                    R"doc(Number of logical observation sets stored in the dataset.)doc" )
            .def_property_readonly( "number_of_observations",
                                    &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getNumberOfObservations,
                                    R"doc(Number of observation rows stored in the dataset.)doc" )
            .def_property_readonly( "total_scalar_size",
                                    &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getTotalScalarSize,
                                    R"doc(Number of scalar observation components stored in the dataset.)doc" )
            .def_property_readonly(
                    "observation_set_metadata",
                    py::overload_cast<>( &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationSetMetadata, py::const_ ),
                    R"doc(List of metadata entries for all observation sets.)doc" )
            .def( "get_observation_set_metadata",
                  py::overload_cast< const tom::ObservationSetId >(
                          &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationSetMetadata, py::const_ ),
                  py::arg( "set_id" ),
                  R"doc(Return metadata for one observation set.)doc" )
            .def_property_readonly( "observation_rows",
                                    &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationRows,
                                    R"doc(Row metadata for all observations.)doc" )
            .def( "observation_row",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationRow,
                  py::arg( "observation_id" ),
                  R"doc(Return row metadata for one observation.)doc" )
            .def_property_readonly( "scalar_component_rows",
                                    &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getScalarComponentRows,
                                    R"doc(Row metadata for all scalar components.)doc" )
            .def( "scalar_component_row",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getScalarComponentRow,
                  py::arg( "scalar_component_id" ),
                  R"doc(Return metadata for one scalar component.)doc" )
            .def( "observation_ids_for_set",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationIdsForSet,
                  py::arg( "set_id" ),
                  R"doc(Return the observation row identifiers belonging to a set.)doc" )
            .def( "observations_for_set",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationsForSet,
                  py::arg( "set_id" ),
                  R"doc(Return all observation vectors in a set.)doc" )
            .def( "observation_vector_for_set",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationVectorForSet,
                  py::arg( "set_id" ),
                  R"doc(Return all observations in a set as one concatenated vector.)doc" )
            .def( "observation_value",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationValue,
                  py::arg( "observation_id" ),
                  R"doc(Return the vector-valued observation for one observation row.)doc" )
            .def( "observation_times_for_set",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationTimesForSet,
                  py::arg( "set_id" ),
                  R"doc(Return all observation times in a set.)doc" )
            .def( "observation_time",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationTime,
                  py::arg( "observation_id" ),
                  R"doc(Return the time of one observation row.)doc" )
            .def( "weights_for_set",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getWeightsForSet,
                  py::arg( "set_id" ),
                  R"doc(Return all weight vectors in a set.)doc" )
            .def( "weight_vector_for_set",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getWeightVectorForSet,
                  py::arg( "set_id" ),
                  R"doc(Return all weights in a set as one concatenated vector.)doc" )
            .def( "weight_matrix_for_set",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getWeightMatrixForSet,
                  py::arg( "set_id" ),
                  R"doc(
Return the full dense weight matrix for a set.

If a set-level matrix was stored, that matrix is returned. Otherwise the matrix
is assembled from the per-observation scalar weights or observable-size blocks
stored for the observations in the set.
)doc" )
            .def( "weight_value",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getWeightValue,
                  py::arg( "observation_id" ),
                  R"doc(Return the weight vector for one observation row.)doc" )
            .def( "weight_matrix_for_observation",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getWeightMatrixForObservation,
                  py::arg( "observation_id" ),
                  R"doc(
Return the observable-size dense weight matrix for one observation row.

Scalar or vector diagonal weights are materialized as a diagonal matrix.
)doc" )
            .def( "residuals_for_set",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getResidualsForSet,
                  py::arg( "set_id" ),
                  R"doc(Return all residual vectors in a set.)doc" )
            .def( "residual_vector_for_set",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getResidualVectorForSet,
                  py::arg( "set_id" ),
                  R"doc(Return all residuals in a set as one concatenated vector.)doc" )
            .def( "residual_value",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getResidualValue,
                  py::arg( "observation_id" ),
                  R"doc(Return the residual vector for one observation row.)doc" )
            .def( "dependent_variables_for_set",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getDependentVariablesForSet,
                  py::arg( "set_id" ),
                  R"doc(Return all dependent-variable vectors stored for a set.)doc" )
            .def( "dependent_variables",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getDependentVariables,
                  py::arg( "observation_id" ),
                  R"doc(Return the dependent-variable vector for one observation row.)doc" )
            .def( "single_dependent_variable_for_set",
                  py::overload_cast< const tom::ObservationSetId,
                                     const std::shared_ptr< tss::ObservationDependentVariableSettings >&,
                                     const bool >(
                          &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getSingleDependentVariableForSet, py::const_ ),
                  py::arg( "set_id" ),
                  py::arg( "dependent_variable_settings" ),
                  py::arg( "return_first_compatible_settings" ) = false,
                  R"doc(Return the values of one dependent variable stored for a set.)doc" )
            .def( "single_dependent_variable_for_set_by_index",
                  py::overload_cast< const tom::ObservationSetId, const std::pair< int, int >& >(
                          &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getSingleDependentVariableForSet, py::const_ ),
                  py::arg( "set_id" ),
                  py::arg( "dependent_variable_index_and_size" ),
                  R"doc(Return dependent-variable values by stored index and size.)doc" )
            .def( "compatible_dependent_variable_settings_for_set",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getCompatibleDependentVariableSettingsForSet,
                  py::arg( "set_id" ),
                  py::arg( "dependent_variable_settings" ),
                  R"doc(Return dependent-variable settings in a set compatible with the requested settings.)doc" )
            .def( "all_compatible_dependent_variables_for_set",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getAllCompatibleDependentVariablesForSet,
                  py::arg( "set_id" ),
                  py::arg( "dependent_variable_settings" ),
                  R"doc(Return all dependent-variable values in a set compatible with the requested settings.)doc" )
            .def( "set_dependent_variables_for_set",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::setDependentVariablesForSet,
                  py::arg( "set_id" ),
                  py::arg( "dependent_variables" ),
                  R"doc(Replace all dependent-variable vectors in a set.)doc" )
            .def( "clear_dependent_variables_for_set",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::clearDependentVariablesForSet,
                  py::arg( "set_id" ),
                  R"doc(Clear all dependent-variable vectors in a set.)doc" )
            .def( "number_of_observations_for_set",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getNumberOfObservationsForSet,
                  py::arg( "set_id" ),
                  R"doc(Return the number of observation rows in a set.)doc" )
            .def( "total_scalar_size_for_set",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getTotalScalarSizeForSet,
                  py::arg( "set_id" ),
                  R"doc(Return the number of scalar components in a set.)doc" )
            .def( "link_definition",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getLinkDefinition,
                  py::arg( "link_definition_id" ),
                  R"doc(Return a link definition from the dataset registry.)doc" )
            .def_property_readonly( "number_of_link_definitions",
                                    &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getNumberOfLinkDefinitions,
                                    R"doc(Number of unique link definitions registered in the dataset.)doc" )
            .def( "ancillary_settings",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getAncillarySettings,
                  py::arg( "ancillary_settings_id" ),
                  R"doc(Return ancillary settings from the dataset registry.)doc" )
            .def( "dependent_variable_bookkeeping",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getDependentVariableBookkeeping,
                  py::arg( "dependent_variable_layout_id" ),
                  R"doc(Return dependent-variable bookkeeping from the dataset registry.)doc" )
            .def( "observation_set_ids",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationSetIds,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(Return observation set identifiers selected by an observation parser.)doc" )
            .def( "observation_ids_matching_condition",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationIdsMatchingCondition,
                  py::arg( "condition" ),
                  R"doc(Return observation row identifiers selected by a row-level condition.)doc" )
            .def( "create_viewer",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::createViewer,
                  py::keep_alive< 0, 1 >( ),
                  py::arg( "condition" ),
                  R"doc(Return a read-only viewer over observations selected by a row-level condition.)doc" )
            .def( "create_new_and_keep",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::createNewAndKeep,
                  py::arg( "condition" ),
                  R"doc(Return a new dataset containing only observations selected by a row-level condition.)doc" )
            .def( "create_new_and_drop",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::createNewAndDrop,
                  py::arg( "condition" ),
                  R"doc(Return a new dataset excluding observations selected by a row-level condition.)doc" )
            .def( "reject_observations",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::rejectObservations,
                  py::arg( "condition" ),
                  py::arg( "reason" ) = "",
                  R"doc(Mark observations selected by a row-level condition as rejected without deleting them.)doc" )
            .def( "restore_observations",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::restoreObservations,
                  py::arg( "condition" ),
                  R"doc(Restore observations selected by a row-level condition to active status.)doc" )
            .def( "active_estimation_vector_projection",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::createEstimationProjection,
                  py::arg( "include_rejected" ) = false,
                  R"doc(Return an estimation projection. Rejected observations are excluded by default.)doc" )
            .def( "legacy_vector_projection",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::createLegacyProjection,
                  py::arg( "include_inactive" ) = true,
                  R"doc(Return the legacy flat projection. Inactive/rejected observations are included by default.)doc" )
            .def( "estimation_vector_projection",
                  &tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >::createLegacyProjection,
                  py::arg( "include_inactive" ) = true,
                  R"doc(Return a flat estimation-vector projection of the dataset.)doc" );

    m.def( "create_observation_dataset_from_single_observation_set",
           static_cast< std::shared_ptr< tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE > > ( * )(
                   const std::shared_ptr< tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE > >& ) >(
                   &tom::createObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE > ),
           py::arg( "observation_set" ),
           R"doc(Return the dataset backend used by a legacy ``SingleObservationSet`` wrapper.)doc" );

    m.def( "create_observation_dataset_from_collection",
           py::overload_cast< const std::shared_ptr< tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > >& >(
                   &tom::createObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE > ),
           py::arg( "observation_collection" ),
           R"doc(Return a dataset representation of a legacy ``ObservationCollection``.)doc" );

    m.def( "create_single_observation_set_from_dataset",
           static_cast< std::shared_ptr< tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE > > ( * )(
                   const std::shared_ptr< tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE > >&, const tom::ObservationSetId ) >(
                   &tom::createSingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE > ),
           py::arg( "observation_dataset" ),
           py::arg( "set_id" ) = 0,
           R"doc(Create a legacy ``SingleObservationSet`` wrapper around one dataset set.)doc" );

    m.def( "create_observation_collection_from_dataset",
           &tom::createObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "observation_dataset" ),
           R"doc(Create a legacy ``ObservationCollection`` wrapper around a dataset.)doc" );

    // SINGLE OBSERVATION SET

    py::class_< tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE > > >( m,
                                                                                                "SingleObservationSet",
                                                                                                R"doc(
        Class for storing and handling a single set of observations.

        This class stores a set of observations of a single observable type and a single link definition.
        The observations are stored as a vector of Eigen matrices, where each entry represents a single
        observation at a given time. The class also stores the observation times, reference link end,
        and other metadata.

        The pybind ``SingleObservationSet`` object is a list-like object, which can be iterated over, and from which specific
        observation data can be retrieved by index.

        When iterating, a tuple is returned with the following entries:
        - ``int``: index of the observation.
        - ``float``: time of the observation.
        - ``numpy.ndarray``: value of the observation.

        When using the ``[]`` operator, a tuple is returned with the following entries:
        - ``float``: time of the observation.
        - ``numpy.ndarray``: value of the observation.

        Parameters
        ----------
        observable_type : tudatpy.estimation.observable_models_setup.model_settings.ObservableType
            Type of observable.
        link_ends : tudatpy.estimation.observable_models_setup.links.LinkDefinition
            Definition of the link ends for the observation.
        observations : list[numpy.ndarray]
            List of observations. Each entry is a vector representing a single observation.
        observation_epochs : list[float]
            List of observation times.
        reference_link_end : tudatpy.estimation.observable_models_setup.links.LinkEndType
            Reference link end for the observation.
        observation_dependent_variables : list[numpy.ndarray], optional
            List of dependent variables for each observation.
        dependent_variable_calculator : tudatpy.estimation.observations.ObservationDependentVariableCalculator, optional
            Calculator for dependent variables.
        ancillary_settings : tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncillarySimulationSettings, optional
            Ancillary settings for the observation.
      )doc" )
            .def( py::init< const tom::ObservableType,
                            const tom::LinkDefinition,
                            const std::vector< Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >,
                            const std::vector< TIME_TYPE >,
                            const tom::LinkEndType,
                            const std::vector< Eigen::VectorXd >,
                            const std::shared_ptr< tss::ObservationDependentVariableBookkeeping >,
                            const std::shared_ptr< tom::ObservationAncillarySimulationSettings > >( ),
                  py::arg( "observable_type" ),
                  py::arg( "link_ends" ),
                  py::arg( "observations" ),
                  py::arg( "observation_epochs" ),
                  py::arg( "reference_link_end" ),
                  py::arg( "observation_dependent_variables" ) = std::vector< Eigen::VectorXd >( ),
                  py::arg( "dependent_variable_bookkeeping" ) = nullptr,
                  py::arg( "ancillary_settings" ) = nullptr )
            .def( "set_observations",
                  py::overload_cast< const std::vector< Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >& >(
                          &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::setObservations ),
                  py::arg( "observations" ),
                  R"doc(
Sets the observation values from a list of vectors.

Parameters
----------
observations : list[numpy.ndarray]
    The new list of observations.
)doc" )
            .def( "set_observations",
                  py::overload_cast< const Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 >& >(
                          &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::setObservations ),
                  py::arg( "observations" ),
                  R"doc(
Sets the observation values from a single concatenated vector.

Parameters
----------
observations : numpy.ndarray
    A single vector containing all observations concatenated.
)doc" )
            .def( "set_residuals",
                  py::overload_cast< const std::vector< Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >& >(
                          &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::setResiduals ),
                  py::arg( "residuals" ),
                  R"doc(
Sets the residuals for all observations from a list of vectors.

Parameters
----------
residuals : list[numpy.ndarray]
    The new list of residuals.
)doc" )
            .def( "set_residuals",
                  py::overload_cast< const Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 >& >(
                          &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::setResiduals ),
                  py::arg( "residuals" ),
                  R"doc(
Sets the residuals for all observations from a single concatenated vector.

Parameters
----------
residuals : numpy.ndarray
    A single vector containing all residuals concatenated.
)doc" )
            .def( "set_constant_weight",
                  py::overload_cast< const double >( &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantWeight ),
                  py::arg( "weight" ),
                  R"doc(
Sets a constant scalar weight for all observations.

Parameters
----------
weight : float
    The constant weight to apply.
)doc" )
            .def( "set_constant_weight",
                  py::overload_cast< const Eigen::Matrix< double, Eigen::Dynamic, 1 >& >(
                          &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantWeight ),
                  py::arg( "weight" ),
                  R"doc(
Sets a constant weight vector for all observations.

Parameters
----------
weight : numpy.ndarray
    The constant weight vector to apply to each observation.
)doc" )
            .def( "set_tabulated_weights",
                  py::overload_cast< const Eigen::VectorXd& >(
                          &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::setTabulatedWeights ),
                  py::arg( "weights" ),
                  R"doc(
Sets weights for all observations from a single concatenated vector.

Parameters
----------
weights : numpy.ndarray
    A single vector containing all weights concatenated.
)doc" )
            .def( "filter_observations",
                  &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::filterObservations,
                  py::arg( "filter" ),
                  py::arg( "save_filtered_obs" ) = true,
                  R"doc(
Filters observations based on a given filter criterion.

Parameters
----------
filter : tudatpy.estimation.observations.observations_processing.ObservationFilterBase
    The filter to apply.
save_filtered_obs : bool, optional
    If true, the filtered observations are stored in a separate set. Defaults to true.
)doc" )
            .def_property_readonly( "observable_type",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getObservableType,
                                    R"doc(

         **read-only**

         Type of observable for which the object stores observations

         :type: ObservableType
      )doc" )
            .def_property( "link_definition",
                           &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getLinkEnds,
                           &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::setLinkEnds,
                           R"doc(

         Definition of the link ends for which the object stores observations

         :type: LinkDefinition
      )doc" )
            .def_property_readonly( "reference_link_end",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getReferenceLinkEnd,
                                    R"doc(

         **read-only**

         Reference link end for all stored observations

         :type: LinkEndType
      )doc" )
            .def_property_readonly( "number_of_observables",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getNumberOfObservables,
                                    R"doc(
Returns the number of observations in the set.

Returns
-------
int
    The number of observations.
)doc" )
            .def_property_readonly( "single_observable_size",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getSingleObservableSize,
                                    R"doc(
Returns the size of a single observation value (e.g., 1 for range, 2 for angular position).

Returns
-------
int
    The size of a single observation.
)doc" )
            .def_property_readonly( "total_observation_set_size",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getTotalObservationSetSize,
                                    R"doc(
Returns the total size of all observation values in the set.

Returns
-------
int
    The total size of the observation set.
)doc" )
            .def_property_readonly( "time_bounds",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getTimeBounds,
                                    R"doc(
Returns the start and end time of the observation set.

Returns
-------
tuple[tudatpy.astro.time_representation.Time, tudatpy.astro.time_representation.Time]
    The start and end time of the observations.
)doc" )
            .def_property_readonly( "list_of_observations",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getObservations,
                                    R"doc(

         **read-only**

         List of separate stored observations. Each entry of this list is a vector containing a single observation. In cases where the observation is single-valued (range, Doppler), the vector is size 1, but for multi-valued observations such as angular position, each vector in the list will have size >1

         :type: list[ numpy.ndarray[numpy.float64[m, 1]] ]
      )doc" )
            .def_property_readonly( "observation_times",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationTimes,
                                    R"doc(

         **read-only**

         List of reference times for each of the observations in ``list_of_observations``

         :type: list[ tudatpy.astro.time_representation.Time ]
      )doc" )
            .def_property_readonly( "concatenated_observations",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationsVector,
                                    R"doc(

         **read-only**

         Concatenated vector of all stored observations

         :type: numpy.ndarray[numpy.float64[m, 1]]
      )doc" )
            .def_property_readonly( "computed_observations",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getComputedObservations,
                                    R"doc(
Returns the vector of computed observation values (observed - residual).

Returns
-------
list[numpy.ndarray]
    The list of computed observations.
)doc" )
            .def_property_readonly( "concatenated_computed_observations",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getComputedObservationsVector,
                                    R"doc(
Returns all computed observations concatenated into a single vector.

Returns
-------
numpy.ndarray
    A single vector containing all computed observations.
)doc" )
            .def_property_readonly( "residuals",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getResiduals,
                                    R"doc(
Returns the vector of residuals for all observations.

Returns
-------
list[numpy.ndarray]
    The list of residuals.
)doc" )
            .def_property_readonly( "concatenated_residuals",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getResidualsVector,
                                    R"doc(
Returns all residuals concatenated into a single vector.

Returns
-------
numpy.ndarray
    A single vector containing all residuals.
)doc" )
            .def_property_readonly( "rms_residuals",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getRmsResiduals,
                                    R"doc(
Returns the Root Mean Square (RMS) of the residuals.

Returns
-------
numpy.ndarray
    A vector containing the RMS of residuals for each component of the observable.
)doc" )
            .def_property_readonly( "mean_residuals",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getMeanResiduals,
                                    R"doc(
Returns the mean of the residuals.

Returns
-------
numpy.ndarray
    A vector containing the mean of residuals for each component of the observable.
)doc" )
            .def_property_readonly( "weights", &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getWeights, R"doc(
Returns the vector of weights for all observations.

Returns
-------
list[numpy.ndarray]
    The list of weights.
)doc" )
            .def_property_readonly( "concatenad_weights",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getWeightsVector,
                                    R"doc(
Returns all weights concatenated into a single vector.

Returns
-------
numpy.ndarray
    A single vector containing all weights.
)doc" )
            .def_property( "dependent_variables",
                           &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationsDependentVariables,
                           &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::setObservationsDependentVariables,
                           R"doc(List of dependent variables for all observations.)doc" )
            .def_property_readonly( "dependent_variables_history",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getDependentVariableHistory,
                                    R"doc(
Time history of all observation dependent variables.

Returns
-------
dict[tudatpy.astro.time_representation.Time, numpy.ndarray]
    A map from observation time to the vector of dependent variables.
)doc" )
            .def_property_readonly( "observations_history",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationsHistory,
                                    R"doc(

         **read-only**

         Dictionary of observations sorted by time. Created by making a dictionary with ``observation_times`` as keys and ``list_of_observations`` as values.

         :type: dict[ tudatpy.astro.time_representation.Time, numpy.ndarray[numpy.float64[m, 1]] ]
      )doc" )
            .def_property_readonly( "ancillary_settings",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getAncillarySettings,
                                    R"doc(

         **read-only**

         Ancillary settings for all stored observations

         :type: ObservationAncillarySimulationSettings
      )doc" )
            .def_property( "weights_vector",
                           &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getWeightsVector,
                           &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::setTabulatedWeights,
                           R"doc(Concatenated vector of weights for all observations.)doc" )
            .def_property_readonly( "filtered_observation_set",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getFilteredObservationSet,
                                    R"doc(
Returns the set of filtered observations.

Returns
-------
SingleObservationSet
    The observation set containing filtered observations.
)doc" )
            .def_property_readonly( "number_filtered_observations",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getNumberOfFilteredObservations,
                                    R"doc(
Returns the number of observations that have been filtered out.

Returns
-------
int
    The number of filtered observations.
)doc" )
            .def( "single_dependent_variable",
                  py::overload_cast< std::shared_ptr< tss::ObservationDependentVariableSettings >, const bool >(
                          &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getSingleDependentVariable ),
                  py::arg( "dependent_variable_settings" ),
                  py::arg( "return_first_compatible_settings" ) = false,
                  R"doc(
Returns the values of a single dependent variable (specified by dependent variable settings).

Parameters
----------
dependent_variable_settings : tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings
    Settings for the dependent variable to retrieve.
return_first_compatible_settings : bool, optional
    If true, returns the first compatible variable found. Defaults to false.

Returns
-------
numpy.ndarray
    A matrix of the dependent variable values over all observation times.
)doc" )
            .def( "compatible_dependent_variable_settings",
                  &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getCompatibleDependentVariablesSettingsList,
                  py::arg( "dependent_variable_settings" ),
                  R"doc(
Returns the list of all dependent variable settings compatible with the settings provided as inputs.

Parameters
----------
dependent_variable_settings : tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings
    The settings to check for compatibility.

Returns
-------
list[tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings]
    A list of compatible dependent variable settings.
)doc" )
            .def( "compatible_dependent_variables_list",
                  &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getAllCompatibleDependentVariables,
                  py::arg( "dependent_variable_settings" ),
                  R"doc(
Returns a vector containing the values of all dependent variables compatible with the settings provided as input.

Parameters
----------
dependent_variable_settings : tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings
    The settings to check for compatibility.

Returns
-------
list[numpy.ndarray]
    A list of matrices, each containing values of a compatible dependent variable.
)doc" )
            .def( "single_dependent_variable_history",
                  &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getSingleDependentVariableHistory,
                  py::arg( "dependent_variable_settings" ),
                  py::arg( "return_first_compatible_settings" ) = false,
                  R"doc(
Returns the time history of a single dependent variable (specified by settings).

Parameters
----------
dependent_variable_settings : tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings
    Settings for the dependent variable to retrieve.
return_first_compatible_settings : bool, optional
    If true, returns the first compatible variable found. Defaults to false.

Returns
-------
dict[float, numpy.ndarray]
    A map from observation time to the value of the specified dependent variable.
)doc" )
            .def_property_readonly( "dependent_variables_matrix",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationsDependentVariablesMatrix,
                                    R"doc(
Returns the dependent variables for all observations as a matrix.

Returns
-------
numpy.ndarray
    A matrix where each row corresponds to an observation and columns to dependent variables.
)doc" );

    m.def( "single_observation_set",
           &tss::singleObservationSetWithoutDependentVariables< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "observable_type" ),
           py::arg( "link_definition" ),
           py::arg( "observations" ),
           py::arg( "observation_times" ),
           py::arg( "reference_link_end" ),
           py::arg( "ancillary_settings" ) = nullptr,
           R"doc(

        Deprecated. Use :func:`~tudatpy.estimation.observations.create_single_observation_set` instead.

        )doc" );

    m.def( "create_single_observation_set",
           py::overload_cast< const tom::ObservableType,
                              const tom::LinkEnds&,
                              const std::vector< Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >&,
                              const std::vector< TIME_TYPE >,
                              const tom::LinkEndType,
                              const std::shared_ptr< tom::ObservationAncillarySimulationSettings > >(
                   &tom::createSingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE > ),
           py::arg( "observable_type" ),
           py::arg( "link_ends" ),
           py::arg( "observations" ),
           py::arg( "observation_times" ),
           py::arg( "reference_link_end" ),
           py::arg( "ancillary_settings" ) = nullptr,
           R"doc(
        Factory function to create a `SingleObservationSet` object.

        This function creates a `SingleObservationSet` object from a list of observations and their corresponding times.
        This is a simplified factory function that does not include dependent variables.

        Parameters
        ----------
        observable_type : :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservableType`
            Type of observable.
        link_ends : dict[:class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType`, :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndId`]
            Definition of the link ends for the observation.
        observations : list[numpy.ndarray]
            List of observations. Each entry is a vector representing a single observation.
        observation_times : list[float]
            List of observation times.
        reference_link_end : :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType`
            Reference link end for the observation.
        ancillary_settings : :class:`~tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncillarySimulationSettings`, optional
            Ancillary settings for the observation.

        Returns
        -------
        tudatpy.estimation.observations.SingleObservationSet
            A `SingleObservationSet` object.
        )doc" );

    // OBSERVATION COLLECTION

    py::class_< tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > > >( m,
                                                                                                 "ObservationCollection",
                                                                                                 R"doc(

         Class collecting all observations and associated data for use in an estimation.

         Class containing the full set of observations and associated data, typically for input into the estimation. When using simulated data,
         this class is instantiated via a call to the :func:`~tudatpy.estimation.observations_setup.observations_wrapper.simulate_observations` function. More information is provided
         on the `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`__





      )doc" )
            .def( py::init< std::vector< std::shared_ptr< tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE > > > >( ),
                  py::arg( "observation_sets" ),
                  R"doc(
Constructor for the ObservationCollection class.

Parameters
----------
observation_sets : list[tudatpy.estimation.observations.SingleObservationSet]
    List of single observation sets to be included in the collection.
)doc" )
            .def_property_readonly( "concatenated_times",
                                    &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedDoubleTimeVector,
                                    R"doc(

         **read-only**

         Vector containing concatenated observation times. See `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`__ for details on storage order

         :type: numpy.ndarray[numpy.float64[m, 1]]
      )doc" )
            .def_property_readonly( "concatenated_times_objects",
                                    &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedTimeVector,
                                    R"doc(

         **read-only**

         Vector containing concatenated observation times. See `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`__ for details on storage order

         :type: numpy.ndarray[numpy.float64[m, 1]]
      )doc" )
            .def_property_readonly( "concatenated_weights",
                                    &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getUnparsedConcatenatedWeights,
                                    R"doc(
         **read-only**

         Vector containing concatenated observation weights.

         :type: numpy.ndarray[numpy.float64[m, 1]]
)doc" )
            .def_property_readonly( "concatenated_observations",
                                    &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationVector,
                                    R"doc(

         **read-only**

         Vector containing concatenated observable values. See `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`__ for details on storage order

         :type: numpy.ndarray[numpy.float64[m, 1]]
      )doc" )
            .def_property_readonly(
                    "concatenated_link_definition_ids",
                    py::overload_cast<>( &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedLinkEndIds ),
                    R"doc(

         **read-only**

         Vector containing concatenated indices identifying the link ends. Each set of link ends is assigned a unique integer identifier (for a given instance of this class). The definition of a given integer identifier with the link ends is given by this class' :func:`link_definition_ids` function. See `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`__ for details on storage order of the present vector.

         :type: numpy.ndarray[ int ]
      )doc" )
            .def_property_readonly( "link_definition_ids",
                                    &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getInverseLinkEndIdentifierMap,
                                    R"doc(

         **read-only**

         Dictionary mapping a link end integer identifier to the specific link ends

         :type: dict[ int, dict[ LinkEndType, LinkEndId ] ]
      )doc" )
            .def_property_readonly( "observable_type_start_index_and_size",
                                    &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationTypeStartAndSize,
                                    R"doc(

         **read-only**

         Dictionary defining per observable type (dict key), the index in the full observation vector (:func:`concatenated_observations`) where the given observable type starts, and the number of subsequent entries in this vector containing a value of an observable of this type

         :type: dict[ ObservableType, [ int, int ] ]
      )doc" )
            .def_property_readonly(
                    "observation_set_start_index_and_size",
                    &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationSetStartAndSizePerLinkEndIndex,
                    R"doc(

         **read-only**

         The nested dictionary/list returned by this property mirrors the structure of the :func:`sorted_observation_sets` property of this class. The present function provides the start index and size of the observables in the full observation vector that come from the corresponding `SingleObservationSet` in the :func:`sorted_observation_sets` Consequently, the present property returns a nested dictionary defining per observable type, link end identifier, and `SingleObservationSet` index (for the given observable type and link end identifier), where the observables in the given `SingleObservationSet` starts, and the number of subsequent entries in this vector containing data from it.

         :type: dict[ ObservableType, dict[ int, list[ int, int ] ] ]
      )doc" )
            .def_property_readonly( "observation_vector_size",
                                    &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getTotalObservableSize,
                                    R"doc(

         **read-only**

         Length of the total vector of observations

         :type: int
      )doc" )
            .def_property_readonly( "sorted_observation_sets",
                                    &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getSortedObservationSets,
                                    R"doc(

         **read-only**

         The nested dictionary/list contains the list of `SingleObservationSet` objects, in the same method as they are stored internally in the present class. Specifics on the storage order are given in the `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`__

         :type: dict[ ObservableType, dict[ int, list[ SingleObservationSet ] ] ]
      )doc" )
            .def_property_readonly( "link_ends_per_observable_type",
                                    &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getLinkEndsPerObservableType,
                                    R"doc(
         **read-only**

         Dictionary mapping each observable type to a list of link ends for which observations are available.

         :type: dict[ ObservableType, list[LinkEnds] ]
)doc" )
            .def_property_readonly( "link_definitions_per_observable",
                                    &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getLinkDefinitionsPerObservable,
                                    R"doc(
         **read-only**

         Dictionary mapping each observable type to a list of link definitions.

         :type: dict[ ObservableType, list[LinkDefinition] ]
)doc" )
            .def_property_readonly( "time_bounds",
                                    &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getTimeBoundsDouble,
                                    R"doc(
         **read-only**

         Pair of minimum and maximum observation time in the collection.

         :type: tuple[float, float]
)doc" )
            .def_property_readonly( "time_bounds_time_object",
                                    &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getTimeBounds,
                                    R"doc(
         **read-only**

         Pair of minimum and maximum observation time in the collection, as Time objects.

         :type: tuple[Time, Time]
)doc" )
            .def_property_readonly( "sorted_per_set_time_bounds",
                                    &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getSortedObservationSetsTimeBounds,
                                    R"doc(
         **read-only**

         Nested dictionary providing the time bounds for each single observation set, sorted by observable type and link end ID.

         :type: dict[ ObservableType, dict[ int, list[tuple[float, float]] ] ]
)doc" )
            .def( "set_observations",
                  py::overload_cast< const Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 >& >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setObservations ),
                  py::arg( "observations" ),
                  R"doc(
Function to reset the full list of observable values. The order of the observations must be the same as for :attr:`~ObservationCollection.concatenated_observations`

Parameters
----------
observations : np.ndarray
    New list of observable values
     )doc" )
            .def( "set_observations",
                  py::overload_cast< const Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 >&,
                                     const std::shared_ptr< tom::ObservationCollectionParser > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setObservations ),
                  py::arg( "observations" ),
                  py::arg( "observation_parser" ),
                  R"doc(
Function to reset a subset of all observable values, with this subset defined by the ``observation_parser`` input.
The order and size of the new observation vector must be the same as when calling :attr:`~ObservationCollection.concatenated_observations` on
an ``ObservationCollection`` containing only the parsed observation.

Parameters
----------
observations : np.ndarray
    New list of observable values
observation_parser : ObservationCollectionParser
    Observation parser with which to select the subset of observations that is to be reset
     )doc" )
            .def( "set_observations",
                  py::overload_cast< const std::map< std::shared_ptr< tom::ObservationCollectionParser >,
                                                     Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >& >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setObservations ),
                  py::arg( "observations_per_parser" ),
                  R"doc(
Function to reset a subset of all observable values, with this subset defined by a list of observation parsers input.
Each observation parser is associated with a new set of observable values.
The order and size of the new observation vector for each parser must be the same as when calling :attr:`~ObservationCollection.concatenated_observations` on
an ``ObservationCollection`` containing only the parsed observation (from a single parser). NOTE: since the multiple parsers
are handled in order (iterating over the keys of ``observations_per_parser``) some observations may be reset several times,
in case.

Parameters
----------
observations : np.ndarray
    New list of observable values
observation_parser : ObservationCollectionParser
    Observation parser with which to select the subset of observations that is to be reset
     )doc" )
            .def( "set_residuals",
                  py::overload_cast< const Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 >& >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setResiduals ),
                  py::arg( "residuals" ),
                  R"doc(
Function to reset the full list of observation residuals. The order of the residuals must be the same as for :attr:`~ObservationCollection.concatenated_observations`

Parameters
----------
residuals : np.ndarray
    New list of observation residuals
     )doc" )
            .def( "set_residuals",
                  py::overload_cast< const Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 >&,
                                     const std::shared_ptr< tom::ObservationCollectionParser > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setResiduals ),
                  py::arg( "residuals" ),
                  py::arg( "observation_parser" ),
                  R"doc(
Function to reset a subset of all observation residuals, with this subset defined by the ``observation_parser`` input.

Parameters
----------
residuals : np.ndarray
    New list of observation residuals
observation_parser : ObservationCollectionParser
    Observation parser with which to select the subset of residuals that is to be reset
     )doc" )
            .def( "set_residuals",
                  py::overload_cast< const std::map< std::shared_ptr< tom::ObservationCollectionParser >,
                                                     Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >& >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setResiduals ),
                  py::arg( "residuals_per_parser" ),
                  R"doc(
Function to reset a subset of all observation residuals, with this subset defined by a list of observation parsers input.

Parameters
----------
residuals_per_parser : dict[ObservationCollectionParser, np.ndarray]
    Dictionary mapping observation parsers to new residual values.
     )doc" )
            .def( "get_link_definitions_for_observables",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getLinkDefinitionsForSingleObservable,
                  py::arg( "observable_type" ),
                  R"doc(
         Function to get all link definitions for a given observable type.

         Parameters
         ----------
         observable_type : :class:`tudatpy.estimation.observable_models_setup.model_settings.ObservableType`
             Observable type for which link definitions are to be retrieved.
         Returns
         -------
         list[ LinkDefinition ]
             List of link definitions for the given observable type.
     )doc" )
            .def( "get_single_link_and_type_observations",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getSingleLinkAndTypeObservationSets,
                  py::arg( "observable_type" ),
                  py::arg( "link_definition" ),
                  R"doc(
         Function to get all observation sets for a given observable type and link definition.


         Parameters
         ----------
         observable_type : :class:`tudatpy.estimation.observable_models_setup.model_settings.ObservableType`
             Observable type of which observations are to be simulated.
         link_ends : LinkDefinition
             Link ends for which observations are to be simulated.
         Returns
         -------
         list[ SingleObservationSet ]
             List of observation sets for given observable type and link definition.
     )doc" )
            .def( "get_observable_types",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getObservableTypes,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the observable types for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[ObservableType]
             A list of observable types present in the selected subset.
     )doc" )
            .def( "get_bodies_in_link_ends",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getBodiesInLinkEnds,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the names of bodies present in the link ends of a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[str]
             A list of body names.
     )doc" )
            .def( "get_reference_points_in_link_ends",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getReferencePointsInLinkEnds,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the names of reference points (e.g., ground stations) in the link ends of a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[str]
             A list of reference point names.
     )doc" )
            .def( "get_time_bounds_list",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getTimeBoundsListDouble,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the time bounds for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[tuple[float, float]]
             A list of time bound pairs (start_time, end_time).
     )doc" )
            .def( "get_time_bounds_list_time_object",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getTimeBoundsList,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the time bounds for a subset of observation sets as Time objects.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[tuple[Time, Time]]
             A list of time bound pairs (start_time, end_time) as Time objects.
     )doc" )
            .def( "get_time_bounds_per_set",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getTimeBoundsPerSetDouble,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the time bounds for each set in a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[tuple[float, float]]
             A list of time bound pairs for each observation set.
     )doc" )
            .def( "get_time_bounds_per_set_time_object",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getTimeBoundsPerSet,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the time bounds for each set in a subset of observation sets as Time objects.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[tuple[Time, Time]]
             A list of time bound pairs for each observation set as Time objects.
     )doc" )
            .def( "get_observations",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getObservations,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the observations for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[numpy.ndarray]
             A list of observation vectors, one for each selected observation set.
     )doc" )
            .def( "get_concatenated_observations",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedObservations,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the concatenated observations for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         numpy.ndarray
             The concatenated vector of observation values.
     )doc" )
            .def( "get_observation_times",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationTimesDouble,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the observation times for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[list[float]]
             A list of lists of observation times, one for each selected observation set.
     )doc" )
            .def( "get_observation_times_objects",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationTimes,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the observation times for a subset of observation sets as Time objects.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[list[Time]]
             A list of lists of observation times as Time objects.
     )doc" )
            .def( "get_concatenated_observation_times_objects",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedObservationTimes,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the concatenated observation times for a subset of observation sets as Time objects.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[Time]
             A concatenated list of observation times as Time objects.
     )doc" )
            .def( "get_concatenated_observation_times",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedDoubleObservationTimes,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the concatenated observation times for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[float]
             A concatenated list of observation times.
     )doc" )
            .def( "get_observations_and_times",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationsAndTimesDouble,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the observations and times for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         tuple[list[numpy.ndarray], list[list[float]]]
             A pair containing a list of observation vectors and a list of lists of observation times.
     )doc" )
            .def( "get_observations_and_times_objects",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationsAndTimes,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the observations and times for a subset of observation sets, with times as Time objects.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         tuple[list[numpy.ndarray], list[list[Time]]]
             A pair containing a list of observation vectors and a list of lists of observation times as Time objects.
     )doc" )
            .def( "get_concatenated_observations_and_times",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedObservationsAndTimesDouble,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the concatenated observations and times for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         tuple[numpy.ndarray, list[float]]
             A pair containing the concatenated observation vector and list of times.
     )doc" )
            .def( "get_concatenated_observations_and_times_objects",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedObservationsAndTimes,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the concatenated observations and times for a subset of observation sets, with times as Time objects.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         tuple[numpy.ndarray, list[Time]]
             A pair containing the concatenated observation vector and list of times as Time objects.
     )doc" )
            .def( "get_concatenated_link_definition_ids",
                  py::overload_cast< std::shared_ptr< tom::ObservationCollectionParser > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedLinkEndIds ),
                  py::arg( "observation_parser" ),
                  R"doc(
         Get the concatenated link definition IDs for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser
             Object that is used to select a subset of the observation sets.

         Returns
         -------
         list[int]
             A list of concatenated link end IDs.
     )doc" )
            .def( "get_weights",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getWeights,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the weights for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[numpy.ndarray]
             A list of weight vectors, one for each selected observation set.
     )doc" )
            .def( "get_concatenated_weights",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedWeights,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the concatenated weights for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         numpy.ndarray
             The concatenated vector of weights.
     )doc" )
            .def( "get_residuals",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getResiduals,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the residuals for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[numpy.ndarray]
             A list of residual vectors, one for each selected observation set.
     )doc" )
            .def( "get_concatenated_residuals",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedResiduals,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the concatenated residuals for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         numpy.ndarray
             The concatenated vector of residuals.
     )doc" )
            .def( "get_rms_residuals",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getRmsResiduals,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the RMS of residuals for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[numpy.ndarray]
             A list of RMS residual vectors, one for each selected observation set.
     )doc" )
            .def( "get_mean_residuals",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getMeanResiduals,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the mean of residuals for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[numpy.ndarray]
             A list of mean residual vectors, one for each selected observation set.
     )doc" )
            .def( "get_computed_observations",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getComputedObservations,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the computed observations for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[numpy.ndarray]
             A list of computed observation vectors, one for each selected observation set.
     )doc" )
            .def( "get_concatenated_computed_observations",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedComputedObservations,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the concatenated computed observations for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         numpy.ndarray
             The concatenated vector of computed observation values.
     )doc" )
            .def( "set_constant_weight",
                  py::overload_cast< const double, const std::shared_ptr< tom::ObservationCollectionParser > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantWeight ),
                  py::arg( "weight" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Set a constant weight for a subset of observation sets.

         Parameters
         ----------
         weight : float
             The constant weight to set.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, applying to all observation sets.
     )doc" )
            .def( "set_constant_weight",
                  py::overload_cast< const Eigen::VectorXd, const std::shared_ptr< tom::ObservationCollectionParser > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantWeight ),
                  py::arg( "weight" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Set a constant weight vector for a subset of observation sets.

         Parameters
         ----------
         weight : numpy.ndarray
             The constant weight vector to set.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, applying to all observation sets.
     )doc" )
            .def( "set_constant_weight_per_observation_parser",
                  py::overload_cast< std::map< std::shared_ptr< tom::ObservationCollectionParser >, double > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantWeightPerObservable ),
                  py::arg( "weights_per_observation_parser" ),
                  R"doc(
         Set a constant weight for multiple subsets of observation sets.

         Parameters
         ----------
         weights_per_observation_parser : dict[tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, float]
             A dictionary mapping observation parsers to constant weights.
     )doc" )
            .def( "set_constant_weight_per_observation_parser",
                  py::overload_cast< std::map< std::shared_ptr< tom::ObservationCollectionParser >, Eigen::VectorXd > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantWeightPerObservable ),
                  py::arg( "weights_per_observation_parser" ),
                  R"doc(
         Set a constant weight vector for multiple subsets of observation sets.

         Parameters
         ----------
         weights_per_observation_parser : dict[tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, numpy.ndarray]
             A dictionary mapping observation parsers to constant weight vectors.
     )doc" )
            .def( "set_tabulated_weights",
                  py::overload_cast< const Eigen::VectorXd, const std::shared_ptr< tom::ObservationCollectionParser > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setTabulatedWeights ),
                  py::arg( "tabulated_weights" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Set tabulated weights for a subset of observation sets.

         Parameters
         ----------
         tabulated_weights : numpy.ndarray
             The vector of tabulated weights to set.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, applying to all observation sets.
     )doc" )
            .def( "set_tabulated_weights",
                  py::overload_cast< std::map< std::shared_ptr< tom::ObservationCollectionParser >, Eigen::VectorXd > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setTabulatedWeights ),
                  py::arg( "tabulated_weights" ),
                  R"doc(
         Set tabulated weights for multiple subsets of observation sets.

         Parameters
         ----------
         tabulated_weights : dict[tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, numpy.ndarray]
             A dictionary mapping observation parsers to tabulated weight vectors.
     )doc" )
            .def( "append",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::appendObservationCollection,
                  py::arg( "observation_collection_to_append" ) )
            .def( "filter_observations",
                  py::overload_cast< const std::map< std::shared_ptr< tom::ObservationCollectionParser >,
                                                     std::shared_ptr< tom::ObservationFilterBase > >&,
                                     const bool >( &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::filterObservations ),
                  py::arg( "observation_filters" ),
                  py::arg( "save_filtered_observations" ) = true,
                  R"doc(
         Filter observations using a set of filters.

         This function filters the observations in the collection based on a map of observation filters, each associated with an observation parser.

         Parameters
         ----------
         observation_filters : dict[tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, tudatpy.estimation.observations.observations_processing.ObservationFilterBase]
             A dictionary mapping observation parsers to observation filters.
         save_filtered_observations : bool, optional
             If true, the filtered-out observations are saved within each observation set, by default True.
     )doc" )
            .def( "filter_observations",
                  py::overload_cast< std::shared_ptr< tom::ObservationFilterBase >,
                                     std::shared_ptr< tom::ObservationCollectionParser >,
                                     const bool >( &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::filterObservations ),
                  py::arg( "observation_filters" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  py::arg( "save_filtered_observations" ) = true,
                  R"doc(
         Filter observations using a single filter.

         This function filters a subset of observations (or all) using a single observation filter.

         Parameters
         ----------
         observation_filter : tudatpy.estimation.observations.observations_processing.ObservationFilterBase
             The observation filter to apply.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, applying to all observation sets.
         save_filtered_observations : bool, optional
             If true, the filtered-out observations are saved within each observation set, by default True.
     )doc" )
            .def( "split_observation_sets",
                  py::overload_cast< std::shared_ptr< tom::ObservationSetSplitterBase >,
                                     std::shared_ptr< tom::ObservationCollectionParser > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::splitObservationSets ),
                  py::arg( "observation_set_splitter" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Split observation sets based on a splitter.

         This function splits a subset of observation sets (or all) into smaller sets based on the criteria defined by the splitter.

         Parameters
         ----------
         observation_set_splitter : tudatpy.estimation.observations.observations_processing.ObservationSetSplitterBase
             The splitter to use for splitting the observation sets.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, applying to all observation sets.
     )doc" )
            .def( "get_single_observation_sets",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getSingleObservationSets,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get a subset of single observation sets.

         This function returns a list of the single observation sets that are selected by the provided observation parser.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[tudatpy.estimation.observations.SingleObservationSet]
             A list of the selected single observation sets.
     )doc" )
            .def( "print_observation_sets_start_and_size",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::printObservationSetsStartAndSize,
                  R"doc(Prints the structure of the observation collection, showing the start index and size of each individual observation set.)doc" )
            .def( "remove_single_observation_sets",
                  py::overload_cast< std::shared_ptr< tom::ObservationCollectionParser > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::removeSingleObservationSets ),
                  py::arg( "observation_parser" ),
                  R"doc(
         Remove a subset of single observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser
             Object that is used to select the subset of observation sets to remove.
     )doc" )
            .def( "set_reference_point",
                  py::overload_cast< tss::SystemOfBodies&,
                                     const Eigen::Vector3d&,
                                     const std::string&,
                                     const std::string&,
                                     const tom::LinkEndType,
                                     const std::shared_ptr< tom::ObservationCollectionParser > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setReferencePoint ),
                  py::arg( "bodies" ),
                  py::arg( "antenna_position" ),
                  py::arg( "antenna_name" ),
                  py::arg( "spacecraft_name" ),
                  py::arg( "link_end_type" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Set a fixed reference point for a subset of observations.

         Parameters
         ----------
         bodies : tudatpy.dynamics.environment.SystemOfBodies
             System of bodies.
         antenna_position : numpy.ndarray
             Position of the antenna in the spacecraft body-fixed frame.
         antenna_name : str
             Name of the antenna/reference point.
         spacecraft_name : str
             Name of the spacecraft body.
         link_end_type : LinkEndType
             Link end type to which the reference point should be applied.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Parser to select the observation sets to which the reference point should be applied.
     )doc" )
            .def( "set_reference_points",
                  py::overload_cast< tss::SystemOfBodies&,
                                     const std::map< double, Eigen::Vector3d >&,
                                     const std::string&,
                                     const tom::LinkEndType,
                                     const std::shared_ptr< tom::ObservationCollectionParser > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setReferencePoints ),
                  py::arg( "bodies" ),
                  py::arg( "antenna_switch_history" ),
                  py::arg( "spacecraft_name" ),
                  py::arg( "link_end_type" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Set multiple reference points based on an antenna switch history.

         Parameters
         ----------
         bodies : tudatpy.dynamics.environment.SystemOfBodies
             System of bodies.
         antenna_switch_history : dict[float, numpy.ndarray]
             Dictionary mapping time to antenna position.
         spacecraft_name : str
             Name of the spacecraft body.
         link_end_type : LinkEndType
             Link end type to which the reference points should be applied.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Parser to select the observation sets.
     )doc" )
            .def( "set_reference_point",
                  py::overload_cast< tss::SystemOfBodies&,
                                     const std::shared_ptr< te::Ephemeris >,
                                     const std::string&,
                                     const std::string&,
                                     const tom::LinkEndType,
                                     const std::shared_ptr< tom::ObservationCollectionParser > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setReferencePoint ),
                  py::arg( "bodies" ),
                  py::arg( "antenna_body_fixed_ephemeris" ),
                  py::arg( "antenna_name" ),
                  py::arg( "spacecraft_name" ),
                  py::arg( "link_end_type" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Set a time-varying reference point for a subset of observations.

         Parameters
         ----------
         bodies : tudatpy.dynamics.environment.SystemOfBodies
             System of bodies.
         antenna_body_fixed_ephemeris : tudatpy.dynamics.environment.Ephemeris
             Ephemeris of the antenna in the spacecraft body-fixed frame.
         antenna_name : str
             Name of the antenna/reference point.
         spacecraft_name : str
             Name of the spacecraft body.
         link_end_type : LinkEndType
             Link end type to which the reference point should be applied.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Parser to select the observation sets.
     )doc" )
            .def( "set_transponder_delay",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setTransponderDelay,
                  py::arg( "spacecraft_name" ),
                  py::arg( "transponder_delay" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Deprecated: set the transponder delay for a subset of observations by modifying the
         retransmission delay in their ancillary settings.

         For new simulations, set the default transponder delay on the spacecraft vehicle systems
         before creating the observation model:
         ``bodies.get_body(spacecraft_name).vehicle_systems.transponder_delay = transponder_delay``.

         Parameters
         ----------
         spacecraft_name : str
             Name of the spacecraft with the transponder.
         transponder_delay : float
             The transponder delay in seconds.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Parser to select the observation sets.
     )doc" )
            .def( "remove_empty_observation_sets",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::removeEmptySingleObservationSets,
                  R"doc(Remove all single observation sets that contain no observations.)doc" )
            .def( "add_dependent_variable",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::addDependentVariable,
                  py::arg( "dependent_variable_settings" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Add an observation dependent variable to a subset of the single observation sets.

         Parameters
         ----------
         dependent_variable_settings : tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings
             Settings for the dependent variable to add.
         bodies : tudatpy.dynamics.environment.SystemOfBodies
             System of bodies containing the environment.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Parser to select the observation sets to which the variable should be added.
         
         Returns
         -------
         tudatpy.estimation.observations.observations_processing.ObservationCollectionParser
             A parser that can be used to retrieve the added dependent variable.
     )doc" )
            .def( "dependent_variable",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getDependentVariables,
                  py::arg( "dependent_variable_settings" ),
                  py::arg( "first_compatible_settings" ) = false,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Retrieve the values of a given dependent variable, sorted per single observation set.

         Parameters
         ----------
         dependent_variable_settings : tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings
             Settings for the dependent variable to retrieve.
         first_compatible_settings : bool, optional
             If true, returns the first compatible variable found, by default False.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Parser to select a subset of observation sets.

         Returns
         -------
         tuple[list[numpy.ndarray], tudatpy.estimation.observations.observations_processing.ObservationCollectionParser]
             A pair containing a list of matrices with the dependent variable values and the parser used.
     )doc" )
            .def( "concatenated_dependent_variable",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedDependentVariables,
                  py::arg( "dependent_variable_settings" ),
                  py::arg( "first_compatible_settings" ) = false,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Retrieve the concatenated values of a given dependent variable.

         Parameters
         ----------
         dependent_variable_settings : tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings
             Settings for the dependent variable to retrieve.
         first_compatible_settings : bool, optional
             If true, returns the first compatible variable found, by default False.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Parser to select a subset of observation sets.

         Returns
         -------
         tuple[numpy.ndarray, tudatpy.estimation.observations.observations_processing.ObservationCollectionParser]
             A pair containing a matrix with the concatenated dependent variable values and the parser used.
     )doc" )
            .def( "compatible_dependent_variable_settings",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getCompatibleDependentVariablesSettingsList,
                  py::arg( "dependent_variable_settings" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get a list of all dependent variable settings compatible with the input settings.

         Parameters
         ----------
         dependent_variable_settings : tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings
             Settings for the dependent variable.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Parser to select a subset of observation sets.

         Returns
         -------
         tuple[list[list[tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings]], tudatpy.estimation.observations.observations_processing.ObservationCollectionParser]
             A pair containing a list of lists of compatible settings and the parser used.
     )doc" )
            .def( "compatible_dependent_variables_list",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getAllCompatibleDependentVariables,
                  py::arg( "dependent_variable_settings" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get all dependent variables compatible with the input settings.

         Parameters
         ----------
         dependent_variable_settings : tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings
             Settings for the dependent variable.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Parser to select a subset of observation sets.

         Returns
         -------
         tuple[list[list[numpy.ndarray]], tudatpy.estimation.observations.observations_processing.ObservationCollectionParser]
             A pair containing a list of lists of dependent variable values and the parser used.
     )doc" )
            .def( "dependent_variable_history_per_set",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getDependentVariableHistoryPerObservationSetDouble,
                  py::arg( "dependent_variable_settings" ),
                  py::arg( "first_compatible_settings" ) = false,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Retrieve the time history of a given dependent variable, sorted per observation set.

         Parameters
         ----------
         dependent_variable_settings : tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings
             Settings for the dependent variable to retrieve.
         first_compatible_settings : bool, optional
             If true, returns the first compatible variable found, by default False.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Parser to select a subset of observation sets.

         Returns
         -------
         list[dict[float, numpy.ndarray]]
             A list of maps from time to dependent variable value, one for each set.
     )doc" )
            .def( "dependent_variable_history_per_set_time_object",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getDependentVariableHistoryPerObservationSet,
                  py::arg( "dependent_variable_settings" ),
                  py::arg( "first_compatible_settings" ) = false,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Retrieve the time history of a given dependent variable, sorted per observation set, with times as Time objects.

         Parameters
         ----------
         dependent_variable_settings : tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings
             Settings for the dependent variable to retrieve.
         first_compatible_settings : bool, optional
             If true, returns the first compatible variable found, by default False.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Parser to select a subset of observation sets.

         Returns
         -------
         list[dict[Time, numpy.ndarray]]
             A list of maps from time to dependent variable value, one for each set, with times as Time objects.
     )doc" )
            .def( "dependent_variable_history",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getDependentVariableHistoryDouble,
                  py::arg( "dependent_variable_settings" ),
                  py::arg( "first_compatible_settings" ) = false,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Retrieve the concatenated time history of a given dependent variable.

         Parameters
         ----------
         dependent_variable_settings : tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings
             Settings for the dependent variable to retrieve.
         first_compatible_settings : bool, optional
             If true, returns the first compatible variable found, by default False.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Parser to select a subset of observation sets.

         Returns
         -------
         dict[float, numpy.ndarray]
             A map from time to dependent variable value.
     )doc" )
            .def( "dependent_variable_history_time_object",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getDependentVariableHistory,
                  py::arg( "dependent_variable_settings" ),
                  py::arg( "first_compatible_settings" ) = false,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Retrieve the concatenated time history of a given dependent variable, with times as Time objects.

         Parameters
         ----------
         dependent_variable_settings : tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings
             Settings for the dependent variable to retrieve.
         first_compatible_settings : bool, optional
             If true, returns the first compatible variable found, by default False.
         observation_parser : tudatpy.estimation.observations.observations_processing.ObservationCollectionParser, optional
             Parser to select a subset of observation sets.

         Returns
         -------
         dict[Time, numpy.ndarray]
             A map from time to dependent variable value, with times as Time objects.
     )doc" );

    m.def( "compute_residuals_and_dependent_variables",
           static_cast< void ( * )( std::shared_ptr< tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > >,
                                    const std::vector< std::shared_ptr< tom::ObservationSimulatorBase< STATE_SCALAR_TYPE, TIME_TYPE > > >&,
                                    const tss::SystemOfBodies& ) >(
                   &tss::computeResidualsAndDependentVariables< STATE_SCALAR_TYPE, TIME_TYPE > ),
           py::arg( "observation_collection" ),
           py::arg( "observation_simulators" ),
           py::arg( "bodies" ),
           R"doc(
        Computes residuals and dependent variables for a given observation collection.

        This function simulates observations based on the settings of the input `observation_collection`
        (which typically contains real data). It then computes the residuals by subtracting the
        simulated observations from the original observations. The computed residuals and any
        observation-dependent variables are then stored in the input `observation_collection`.

        Parameters
        ----------
        observation_collection : tudatpy.estimation.observations.ObservationCollection
            The collection of observations for which to compute residuals and dependent variables.
            This object is modified in-place.
        observation_simulators : list[tudatpy.estimation.observable_models.observables_simulation.ObservationSimulator]
            List of observation simulators to be used for computing the observations.
        bodies : tudatpy.dynamics.environment.SystemOfBodies
            The system of bodies required for the observation simulation.
        )doc" );

    m.def( "observation_simulation_settings_from_dataset",
           &tss::getObservationSimulationSettingsFromObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "observation_dataset" ),
           py::arg( "bodies" ),
           R"doc(
        Create observation simulation settings from a dataset.

        The returned settings reproduce the observable type, link definition,
        observation epochs, reference link end, ancillary settings and dependent
        variable bookkeeping stored in each set of the dataset. They can be
        passed to :func:`tudatpy.estimation.observations_setup.observations_wrapper.simulate_observation_dataset`.

        Parameters
        ----------
        observation_dataset : tudatpy.estimation.observations.ObservationDataset
            Dataset from which simulation settings should be reconstructed.
        bodies : tudatpy.dynamics.environment.SystemOfBodies
            System of bodies used for consistency with observation-dependent
            variable setup.

        Returns
        -------
        list[tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings]
            Simulation settings matching the dataset observation sets.
        )doc" );

    m.def(
            "compute_residuals_and_dependent_variables_for_dataset",
            []( const std::shared_ptr< tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE > >& observationDataset,
                const std::vector< std::shared_ptr< tom::ObservationSimulatorBase< STATE_SCALAR_TYPE, TIME_TYPE > > >&
                        observationSimulators,
                const tss::SystemOfBodies& bodies ) {
                if( observationDataset == nullptr )
                {
                    throw std::runtime_error(
                            "Error when computing residuals and dependent variables for dataset, input dataset is None." );
                }
                tss::computeResidualsAndDependentVariables< STATE_SCALAR_TYPE, TIME_TYPE >(
                        observationDataset, observationSimulators, bodies );
            },
            py::arg( "observation_dataset" ),
            py::arg( "observation_simulators" ),
            py::arg( "bodies" ),
            R"doc(
        Compute residuals and dependent variables for a dataset in-place.

        This is the dataset-backed counterpart of
        :func:`compute_residuals_and_dependent_variables`. Residuals and
        dependent variables are written back into the supplied
        :class:`ObservationDataset`.

        Parameters
        ----------
        observation_dataset : tudatpy.estimation.observations.ObservationDataset
            Dataset containing observed values and metadata. This object is
            modified in-place.
        observation_simulators : list[tudatpy.estimation.observable_models.observables_simulation.ObservationSimulator]
            Observation simulators used to compute modeled values.
        bodies : tudatpy.dynamics.environment.SystemOfBodies
            System of bodies required for observation simulation.
        )doc" );

    m.def( "filter_observations",
           py::overload_cast< const std::shared_ptr< tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE > >,
                              const std::shared_ptr< tom::ObservationFilterBase >,
                              const bool >( &tom::filterObservations< STATE_SCALAR_TYPE, TIME_TYPE > ),
           py::arg( "original_observation_set" ),
           py::arg( "observation_filter" ),
           py::arg( "save_filtered_observations" ) = false,
           R"doc(

Deprecated. Use :func:`~tudatpy.estimation.observations.create_filtered_observation_set` instead.


        )doc" );

    m.def( "create_filtered_observation_set",
           py::overload_cast< const std::shared_ptr< tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE > >,
                              const std::shared_ptr< tom::ObservationFilterBase >,
                              const bool >( &tom::filterObservations< STATE_SCALAR_TYPE, TIME_TYPE > ),
           py::arg( "original_observation_set" ),
           py::arg( "observation_filter" ),
           py::arg( "save_filtered_observations" ) = false,
           R"doc(

        Filters a single observation set and returns a new set containing the filtered observations.

        This function creates a copy of the input observation set and applies the given filter.
        The original observation set is not modified.

        Parameters
        ----------
        original_observation_set : tudatpy.estimation.observations.SingleObservationSet
            The observation set to filter.
        observation_filter : tudatpy.estimation.observations.observations_processing.ObservationFilterBase
            The filter to apply.
        save_filtered_observations : bool, optional
            If true, the observations that are filtered out are stored within the new observation set. Defaults to false.

        Returns
        -------
        tudatpy.estimation.observations.SingleObservationSet
            A new observation set containing only the observations that passed the filter.
        )doc" );

    m.def( "split_observation_set",
           py::overload_cast< const std::shared_ptr< tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE > >,
                              const std::shared_ptr< tom::ObservationSetSplitterBase >,
                              const bool >( &tom::splitObservationSet< STATE_SCALAR_TYPE, TIME_TYPE > ),
           py::arg( "original_observation_set" ),
           py::arg( "observation_splitter" ),
           py::arg( "print_warning" ) = true,
           R"doc(
        Splits a single observation set into multiple sets based on a splitter.

        This function takes an observation set and divides it into a list of smaller observation sets
        according to the criteria defined in the observation splitter.

        Parameters
        ----------
        original_observation_set : tudatpy.estimation.observations.SingleObservationSet
            The observation set to split.
        observation_splitter : tudatpy.estimation.observations.observations_processing.ObservationSetSplitterBase
            The splitter defining the splitting criteria.
        print_warning : bool, optional
            If true, a warning is printed if the original set contains filtered observations that will be lost. Defaults to true.

        Returns
        -------
        list[tudatpy.estimation.observations.SingleObservationSet]
            A list of new observation sets resulting from the split.
        )doc" );

    m.def( "merge_observation_collections",
           &tss::mergeObservationCollections< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "observation_collection_list" ) );

    // The following functions create a new ObservationCollection object from an existing one

    m.def( "create_filtered_observation_collection",
           py::overload_cast<
                   const std::shared_ptr< tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > >,
                   const std::map< std::shared_ptr< tom::ObservationCollectionParser >, std::shared_ptr< tom::ObservationFilterBase > >& >(
                   &tom::filterObservations< STATE_SCALAR_TYPE, TIME_TYPE > ),
           py::arg( "original_observation_collection" ),
           py::arg( "observation_filters_map" ),
           R"doc(

        Creates a new, filtered observation collection from an existing one using multiple filters.

        This function applies a set of filters, each associated with a parser, to an observation collection
        and returns a new collection with the filtered observations. The original collection is not modified.

        Parameters
        ----------
        original_observation_collection : tudatpy.estimation.observations.ObservationCollection
            The observation collection to filter.
        observation_filters_map : dict[tudatpy.estimation.observations.ObservationCollectionParser, tudatpy.estimation.observations.observations_processing.ObservationFilterBase]
            A dictionary mapping parsers to filters.

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollection
            A new observation collection containing the filtered observations.

        )doc" );

    m.def( "create_filtered_observation_collection",
           py::overload_cast< const std::shared_ptr< tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > >,
                              const std::shared_ptr< tom::ObservationFilterBase >,
                              const std::shared_ptr< tom::ObservationCollectionParser > >(
                   &tom::filterObservations< STATE_SCALAR_TYPE, TIME_TYPE > ),
           py::arg( "original_observation_collection" ),
           py::arg( "observation_filter" ),
           py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
           R"doc(

        Creates a new, filtered observation collection from an existing one using a single filter.

        This function applies a single filter to a subset of an observation collection (selected by a parser)
        and returns a new collection with the filtered observations. The original collection is not modified.

        Parameters
        ----------
        original_observation_collection : tudatpy.estimation.observations.ObservationCollection
            The observation collection to filter.
        observation_filter : tudatpy.estimation.observations.observations_processing.ObservationFilterBase
            The filter to apply.
        observation_parser : tudatpy.estimation.observations.ObservationCollectionParser, optional
            Parser to select the subset of observations to filter. Defaults to an empty parser (all observations).

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollection
            A new observation collection containing the filtered observations.

        )doc" );

    m.def( "split_observation_collection",
           &tom::splitObservationSets< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "original_observation_collection" ),
           py::arg( "observation_set_splitter" ),
           py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
           R"doc(
        Creates a new observation collection by splitting sets from an existing collection.

        This function splits observation sets from the original collection based on a splitter and
        returns a new collection containing the resulting sets. The original collection is not modified.

        Parameters
        ----------
        original_observation_collection : tudatpy.estimation.observations.ObservationCollection
            The observation collection from which to split sets.
        observation_set_splitter : tudatpy.estimation.observations.observations_processing.ObservationSetSplitterBase
            The splitter defining how to split the sets.
        observation_parser : tudatpy.estimation.observations.ObservationCollectionParser, optional
            Parser to select which observation sets to split. Defaults to an empty parser (all sets).

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollection
            A new observation collection with the split observation sets.
        )doc" );

    m.def( "create_new_observation_collection",
           &tom::createNewObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "original_observation_collection" ),
           py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
           R"doc(

        Creates a new observation collection containing a subset of an existing collection.

        This function selects a subset of observation sets from an original collection using a parser
        and creates a new collection containing only those sets.

        Parameters
        ----------
        original_observation_collection : tudatpy.estimation.observations.ObservationCollection
            The collection from which to extract a subset.
        observation_parser : tudatpy.estimation.observations.ObservationCollectionParser, optional
            Parser to select the observation sets to include in the new collection. Defaults to an empty parser (all sets).

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollection
            A new observation collection containing the selected subset of observation sets.
        )doc" );
}

}  // namespace observations
}  // namespace estimation
}  // namespace tudatpy
