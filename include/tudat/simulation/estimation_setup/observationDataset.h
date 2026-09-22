/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OBSERVATION_DATASET_H
#define TUDAT_OBSERVATION_DATASET_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <initializer_list>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/basics/basicTypedefs.h"
#include "tudat/simulation/estimation_setup/observationVectorData.h"
#include "tudat/simulation/estimation_setup/observationDatasetRows.h"
#include "tudat/simulation/estimation_setup/observationOutput.h"
#include "tudat/simulation/estimation_setup/observationCondition.h"
#include "tudat/simulation/estimation_setup/observationWeights.h"
#include "tudat/simulation/estimation_setup/observationsProcessing.h"
#include "tudat/basics/tudatTypeTraits.h"

namespace tudat
{

namespace observation_models
{
//! Return whether one dataset set is selected by a legacy collection parser.
template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
bool isObservationSetSelectedByLegacyParser( const ObservationDataset< ObservationScalarType, TimeType >& dataset,
                                             const unsigned int setId,
                                             const std::shared_ptr< ObservationCollectionParser >& observationParser );

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type >
class SingleObservationSet;

//! Queryable backend for Tudat observation data.
/*!
 * ObservationDataset stores observations as event rows plus scalar-component
 * arrays. A vector observable has one observation row and N scalar components:
 * observedValues_ is flat by scalar component, while observationRows_ records
 * the event boundary. Scalar reverse mappings are derived from those boundaries.
 * Set-level metadata is stored once in registries and referenced by id. Flat
 * estimator vectors are derived by explicit observation-vector-data builders, not used as the
 * primary data model. Inspection getters copy only requested fields into independent
 * values in internal or estimation order.
 */
template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
class ObservationDataset : public std::enable_shared_from_this< ObservationDataset< ObservationScalarType, TimeType > >
{
public:
    //! Create an empty observation dataset.
    ObservationDataset( ) = default;
    //! Copy an observation dataset while creating an independent lifetime identity.
    ObservationDataset( const ObservationDataset& other );
    //! Move an observation dataset into a new owner.
    ObservationDataset( ObservationDataset&& ) = default;
    //! Move-assign an observation dataset.
    ObservationDataset& operator=( ObservationDataset&& ) = default;
    //! Copy-assign an observation dataset while refreshing its lifetime identity.
    ObservationDataset& operator=( const ObservationDataset& other )
    {
        if( this != &other )
        {
            ObservationDataset copy( other );
            *this = std::move( copy );
        }
        return *this;
    }

    //! Return whether two datasets contain equivalent observation data and metadata.
    bool operator==( const ObservationDataset& rhs ) const
    {
        const auto pointedObjectsEqual = []( const auto& lhs, const auto& rhs ) {
            return static_cast< bool >( lhs ) == static_cast< bool >( rhs ) && ( !lhs || *lhs == *rhs );
        };
        return observationRows_ == rhs.observationRows_ && setMetadata_ == rhs.setMetadata_ &&
                observationIdsBySet_ == rhs.observationIdsBySet_ && linkDefinitionRegistry_ == rhs.linkDefinitionRegistry_ &&
                observedValues_ == rhs.observedValues_ && residualValues_ == rhs.residualValues_ &&
                observationWeights_ == rhs.observationWeights_ &&
                std::equal( ancillarySettingsRegistry_.begin( ),
                            ancillarySettingsRegistry_.end( ),
                            rhs.ancillarySettingsRegistry_.begin( ),
                            rhs.ancillarySettingsRegistry_.end( ),
                            pointedObjectsEqual ) &&
                std::equal( dependentVariableLayoutRegistry_.begin( ),
                            dependentVariableLayoutRegistry_.end( ),
                            rhs.dependentVariableLayoutRegistry_.begin( ),
                            rhs.dependentVariableLayoutRegistry_.end( ),
                            pointedObjectsEqual );
    }

    //////////////////////////////////////////////////////////
    /////////////////       SET CREATION            //////////
    //////////////////////////////////////////////////////////

    //! Add one logical observation set and register its metadata.
    /*!
     * Inputs are supplied as per-observation vectors. The method validates
     * observable dimensions, optionally sorts by time and removes duplicates,
     * registers shared metadata, then appends event rows and scalar components.
     */
    int addObservationSet(
            const ObservableType observableType,
            const LinkDefinition& linkDefinition,
            const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
            const std::vector< TimeType >& times,
            const LinkEndType referenceLinkEnd,
            const std::vector< Eigen::VectorXd >& dependentVariables = std::vector< Eigen::VectorXd >( ),
            const std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping >& dependentVariableBookkeeping = nullptr,
            const std::shared_ptr< ObservationAncillarySimulationSettings >& ancillarySettings = nullptr,
            const std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >& weights =
                    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >( ),
            const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& residuals =
                    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >( ),
            const bool sortObservations = false,
            const bool eraseDuplicateObservations = false );

    //! Copy one set from another dataset into this dataset.
    int addObservationSetFromDataset( const ObservationDataset< ObservationScalarType, TimeType >& sourceDataset,
                                      const unsigned int sourceSetId );

    //! Add a set and apply the requested compact/block weight policy.
    int addObservationSetWithWeights(
            const ObservableType observableType,
            const LinkDefinition& linkDefinition,
            const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
            const std::vector< TimeType >& times,
            const LinkEndType referenceLinkEnd,
            const ObservationWeightSettings& weightSettings,
            const std::vector< Eigen::VectorXd >& dependentVariables = std::vector< Eigen::VectorXd >( ),
            const std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping >& dependentVariableBookkeeping = nullptr,
            const std::shared_ptr< ObservationAncillarySimulationSettings >& ancillarySettings = nullptr,
            const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& residuals =
                    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >( ),
            const bool sortObservations = false );

    //////////////////////////////////////////////////////////
    /////////////////       VALUE MUTATION          //////////
    //////////////////////////////////////////////////////////

    //! Replace the vector-valued measurements for all observation events in one set.
    void setObservationsForSet( const unsigned int setId,
                                const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations );

    //! Replace the vector-valued residuals for all observation events in one set.
    void setResidualsForSet( const unsigned int setId,
                             const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& residuals );

    //////////////////////////////////////////////////////////
    /////////////////       WEIGHTS FUNCTIONALITY   //////////
    //////////////////////////////////////////////////////////

    //! Apply one compact scalar weight to every observation event in one set.
    void setConstantSingleObservationScalarWeightForSet( const unsigned int setId, const double weight );

    //! Apply one component-wise weight vector to every observation event in one set.
    void setConstantSingleObservationDiagonalWeightForSet( const unsigned int setId, const Eigen::VectorXd& weight );

    //! Apply one observable-size N x N weight block to every observation event in one set.
    void setConstantSingleObservationMatrixWeightForSet( const unsigned int setId, const Eigen::MatrixXd& weight );

    //! Set a scalar constant weight on all observations matching a row-level condition.
    void setConstantSingleObservationScalarWeight( const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition,
                                                   const double weight );

    //! Set a vector-valued constant weight on all observations matching a row-level condition.
    void setConstantSingleObservationDiagonalWeight( const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition,
                                                     const Eigen::VectorXd& weight );

    //! Set an observable-size N x N weight block on all observations matching a row-level condition.
    void setConstantSingleObservationMatrixWeight( const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition,
                                                   const Eigen::MatrixXd& weight );

    //! Store one full M x M weight block for an observation set.
    void setWeightMatrixForSet( const unsigned int setId, const Eigen::MatrixXd& weightMatrix );

    //! Return whether the effective set principal block contains off-diagonal weights.
    bool hasWeightMatrixForSet( const unsigned int setId ) const;

    //! Store one observable-size N x N weight block for an observation event.
    void setWeightMatrixForObservation( const unsigned int observationId, const Eigen::MatrixXd& weightMatrix );

    //! Return whether the effective observation block contains off-diagonal weights.
    bool hasWeightMatrixForObservation( const unsigned int observationId ) const;

    //! Store a dense weight block selected by observation ids.
    /*!
     * This is the public interface for correlations between unrelated
     * observations. Empty component lists select all scalar components of each
     * observation. Non-empty component lists are applied to every observation in
     * the corresponding row or column selection. Observation ids are converted
     * to scalar-component ids immediately, so the block remains valid under
     * later observation vector data as long as the dataset structure is unchanged.
     * The public interface always preserves a symmetric total weight matrix:
     * a transposed block is added for different row and column selections, and
     * identical selections must already define a symmetric block.
     */
    void setWeightBlock( const std::vector< unsigned int >& rowObservationIds,
                         const std::vector< unsigned int >& columnObservationIds,
                         const Eigen::MatrixXd& weightBlock,
                         const std::vector< unsigned int >& rowComponents = std::vector< unsigned int >( ),
                         const std::vector< unsigned int >& columnComponents = std::vector< unsigned int >( ) );

    //! Return whether advanced scalar-component weight blocks are stored on the dataset.
    bool hasExtraWeightBlocks( ) const;

    //////////////////////////////////////////////////////////
    /////////////////       SET MUTATION            //////////
    //////////////////////////////////////////////////////////

    //! Append per-observation data to an existing set.
    /*!
     * Missing weights and residuals are filled with unit weights and zero
     * residuals. If requested, the complete set is sorted by time after append.
     */
    void addObservationsToSet( const unsigned int setId,
                               const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
                               const std::vector< TimeType >& times,
                               const std::vector< Eigen::VectorXd >& dependentVariables = std::vector< Eigen::VectorXd >( ),
                               const std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >& weights =
                                       std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >( ),
                               const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& residuals =
                                       std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >( ),
                               const bool sortObservations = true );

    //! Remove observation events from a set by per-set observation index.
    /*!
     * This rebuilds the affected dataset storage once per call. Prefer passing
     * all rows to remove in one call instead of looping over single indices.
     */
    void removeObservationsFromSet( const unsigned int setId, std::vector< unsigned int > indicesToRemove );

    //! Remove all observation events matching a row-level condition.
    void removeObservations( const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition );

    //! Physically remove all currently rejected observation events.
    void removeRejectedObservations( );

    //! Compatibility alias for physically removing all rejected observations.
    void deleteRejectedObservations( )
    {
        removeRejectedObservations( );
    }

    //! Return the time bounds of one observation set.
    std::pair< TimeType, TimeType > getTimeBoundsForSet( const unsigned int setId ) const;

    //! Return time bounds for all observations in the dataset.
    std::pair< TimeType, TimeType > getTimeBounds( ) const;

    //////////////////////////////////////////////////////////
    /////////////////       DATA ACCESSORS          //////////
    //////////////////////////////////////////////////////////

    //! Metadata by set ID: (set metadata, link definition, ancillary settings, dependent-variable layout).
    //! The layout maps (start, size) to cloned settings. Null ancillary settings and empty layouts are preserved.
    using InspectionMetadata = std::map<
            unsigned int,
            std::tuple< ObservationSetMetadata< ObservationScalarType, TimeType >,
                        LinkDefinition,
                        std::shared_ptr< ObservationAncillarySimulationSettings >,
                        std::map< std::pair< int, int >, std::shared_ptr< simulation_setup::ObservationDependentVariableSettings > > > >;

    //! Inspection results are independent, mutable value snapshots in both directions.
    //! All inspection getters include rejected rows by default and default to internal order.
    //! Ordering never changes membership. Selection and ordering do not mutate the dataset.
    //! Event-aligned time values, preserving TimeType.
    std::vector< TimeType > getTimes( const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition =
                                              ObservationSelectionCondition< ObservationScalarType, TimeType >::all( ),
                                      const ObservationOrdering ordering = ObservationOrdering::internal ) const;

    //! One owning observation vector per event; mixed dimensions are supported.
    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > getObservations(
            const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition =
                    ObservationSelectionCondition< ObservationScalarType, TimeType >::all( ),
            const ObservationOrdering ordering = ObservationOrdering::internal ) const;

    //! One owning residual vector per event; omitted residuals retain the stored zeros.
    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > getResiduals(
            const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition =
                    ObservationSelectionCondition< ObservationScalarType, TimeType >::all( ),
            const ObservationOrdering ordering = ObservationOrdering::internal ) const;

    //! Stable observation IDs, one per event.
    std::vector< unsigned int > getObservationIds( const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition =
                                                           ObservationSelectionCondition< ObservationScalarType, TimeType >::all( ),
                                                   const ObservationOrdering ordering = ObservationOrdering::internal ) const;

    //! Set IDs, one per event (including repeated IDs).
    std::vector< unsigned int > getSetIds( const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition =
                                                   ObservationSelectionCondition< ObservationScalarType, TimeType >::all( ),
                                           const ObservationOrdering ordering = ObservationOrdering::internal ) const;

    //! Event rows copied recursively, including dependent values and rejection status.
    std::vector< ObservationDatasetRow< TimeType > > getRows(
            const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition =
                    ObservationSelectionCondition< ObservationScalarType, TimeType >::all( ),
            const ObservationOrdering ordering = ObservationOrdering::internal ) const;

    //! One dependent-variable vector per event; missing values are empty vectors.
    std::vector< Eigen::VectorXd > getDependentVariableValues(
            const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition =
                    ObservationSelectionCondition< ObservationScalarType, TimeType >::all( ),
            const ObservationOrdering ordering = ObservationOrdering::internal ) const;

    //! Scalar-aligned (observation ID, component index) pairs, in component order within each event.
    std::vector< std::pair< unsigned int, unsigned int > > getScalarComponents(
            const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition =
                    ObservationSelectionCondition< ObservationScalarType, TimeType >::all( ),
            const ObservationOrdering ordering = ObservationOrdering::internal ) const;

    //! Effective weight diagonal, aligned with getScalarComponents; no correlations are gathered.
    Eigen::VectorXd getWeightDiagonal( const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition =
                                               ObservationSelectionCondition< ObservationScalarType, TimeType >::all( ),
                                       const ObservationOrdering ordering = ObservationOrdering::internal ) const;

    //! Effective selected principal weight matrix, with both axes aligned with getScalarComponents.
    Eigen::SparseMatrix< double > getWeightMatrix( const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition =
                                                           ObservationSelectionCondition< ObservationScalarType, TimeType >::all( ),
                                                   const ObservationOrdering ordering = ObservationOrdering::internal ) const;

    //! Detached metadata for selected sets, keyed by stable set ID; see InspectionMetadata.
    InspectionMetadata getMetadata( const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition =
                                            ObservationSelectionCondition< ObservationScalarType, TimeType >::all( ),
                                    const ObservationOrdering ordering = ObservationOrdering::internal ) const;

    //! Return computed observations as observed-minus-residual values for one set.
    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > getComputedObservationsForSet(
            const unsigned int setId ) const;

    //! Compute component-wise RMS residuals over all observations in one set.
    Eigen::VectorXd getRmsResidualsForSet( const unsigned int setId ) const;

    //! Compute component-wise mean residuals over all observations in one set.
    Eigen::VectorXd getMeanResidualsForSet( const unsigned int setId ) const;

    //! Return the number of logical observation sets.
    std::size_t getNumberOfObservationSets( ) const;

    //! Return the number of observation events.
    std::size_t getNumberOfObservations( ) const;

    //! Return the total number of scalar observation components.
    std::size_t getTotalScalarSize( ) const;

    //! Return metadata for every observation set.
    const std::vector< ObservationSetMetadata< ObservationScalarType, TimeType > >& getObservationSetMetadata( ) const;

    //! Return metadata for one observation set.
    const ObservationSetMetadata< ObservationScalarType, TimeType >& getObservationSetMetadata( const unsigned int setId ) const;

    //! Return all event rows in internal storage order.
    const std::vector< ObservationDatasetRow< TimeType > >& getObservationRows( ) const;

    //! Return the event row with a stable observation identity.
    const ObservationDatasetRow< TimeType >& getObservationRow( const unsigned int observationId ) const;

    //! Return derived rows for every scalar observation component.
    std::vector< ObservationScalarComponentRow > getScalarComponentRows( ) const;

    //! Return the derived row for one scalar observation component.
    ObservationScalarComponentRow getScalarComponentRow( const unsigned int scalarComponentId ) const;

    //! Return stable observation identities belonging to one set.
    const std::vector< unsigned int >& getObservationIdsForSet( const unsigned int setId ) const;

    //! Return one vector-valued measurement per observation event in a set.
    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > getObservationsForSet( const unsigned int setId ) const;

    //! Reconstruct the vector-valued measurement for one observation event.
    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getObservationValue( const unsigned int observationId ) const;

    //! Return one reference-link-end time per observation event in a set.
    std::vector< TimeType > getObservationTimesForSet( const unsigned int setId ) const;

    //! Return the reference-link-end time for one observation event.
    TimeType getObservationTime( const unsigned int observationId ) const;

    //! Return one vector of scalar-component weights per observation event.
    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > getWeightsForSet( const unsigned int setId ) const;

    //! Reconstruct the scalar-component weight vector for one observation event.
    Eigen::Matrix< double, Eigen::Dynamic, 1 > getWeightValue( const unsigned int observationId ) const;

    //! Reconstruct the weight matrix for one observation event.
    Eigen::MatrixXd getWeightMatrixForObservation( const unsigned int observationId ) const;

    //! Materialize the effective principal submatrix for one set, including correlations.
    Eigen::MatrixXd getWeightMatrixForSet( const unsigned int setId ) const;

    //! Return one vector-valued residual per observation event in a set.
    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > getResidualsForSet( const unsigned int setId ) const;

    //! Reconstruct the vector-valued residual for one observation event.
    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getResidualValue( const unsigned int observationId ) const;

    //////////////////////////////////////////////////////////
    /////////////////       DEPENDENT VARIABLES    ///////////
    //////////////////////////////////////////////////////////

    //! Return per-observation dependent-variable vectors for one set.
    std::vector< Eigen::VectorXd > getDependentVariablesForSet( const unsigned int setId ) const;

    //! Return the dependent-variable vector for one observation event.
    Eigen::VectorXd getDependentVariables( const unsigned int observationId ) const;

    //! Extract one dependent-variable block by column start and size.
    Eigen::MatrixXd getSingleDependentVariableForSet( const unsigned int setId,
                                                      const std::pair< int, int >& dependentVariableIndexAndSize ) const;

    //! Extract one dependent-variable block compatible with the requested settings.
    Eigen::MatrixXd getSingleDependentVariableForSet(
            const unsigned int setId,
            const std::shared_ptr< simulation_setup::ObservationDependentVariableSettings >& dependentVariableSettings,
            const bool returnFirstCompatibleSettings = false ) const;

    //! List dependent-variable settings in a set that match the requested settings.
    std::vector< std::shared_ptr< simulation_setup::ObservationDependentVariableSettings > > getCompatibleDependentVariableSettingsForSet(
            const unsigned int setId,
            const std::shared_ptr< simulation_setup::ObservationDependentVariableSettings >& dependentVariableSettings ) const;

    //! Extract every dependent-variable block compatible with the requested settings.
    std::vector< Eigen::MatrixXd > getAllCompatibleDependentVariablesForSet(
            const unsigned int setId,
            const std::shared_ptr< simulation_setup::ObservationDependentVariableSettings >& dependentVariableSettings ) const;

    //! Replace per-observation dependent-variable vectors for one set.
    void setDependentVariablesForSet( const unsigned int setId, const std::vector< Eigen::VectorXd >& dependentVariables );

    //! Clear dependent-variable values for every observation event in one set.
    void clearDependentVariablesForSet( const unsigned int setId );

    //! Return the number of observation events in one set.
    std::size_t getNumberOfObservationsForSet( const unsigned int setId ) const;

    //! Return the scalar observation size of one set.
    std::size_t getTotalScalarSizeForSet( const unsigned int setId ) const;

    //! Return a registered link definition by registry identity.
    const LinkDefinition& getLinkDefinition( const unsigned int linkDefinitionId ) const;

    //! Return the number of registered link definitions.
    std::size_t getNumberOfLinkDefinitions( ) const;

    //! Replace a named reference point for matching link ends in selected observation sets.
    void setLinkEndReferencePoint( const std::string& bodyName,
                                   const std::string& referencePointName,
                                   const LinkEndType linkEndType,
                                   const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition =
                                           ObservationSelectionCondition< ObservationScalarType, TimeType >::all( ) );

    //! Return registered ancillary settings by registry identity.
    const std::shared_ptr< ObservationAncillarySimulationSettings >& getAncillarySettings( const unsigned int ancillarySettingsId ) const;

    //! Return ancillary settings associated with one observation set.
    const std::shared_ptr< ObservationAncillarySimulationSettings >& getAncillarySettingsForSet( const unsigned int setId ) const;

    //! Return registered dependent-variable bookkeeping by registry identity.
    const std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping >& getDependentVariableBookkeeping(
            const unsigned int dependentVariableLayoutId ) const;

    //////////////////////////////////////////////////////////
    /////////////////       SELECTION AND STATUS    //////////
    //////////////////////////////////////////////////////////

    //! Return observation ids for rows that satisfy a new row-level condition.
    std::vector< unsigned int > getObservationIdsMatchingCondition(
            const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition ) const;

    //! Create a new dataset containing only observations satisfying a condition.
    std::shared_ptr< ObservationDataset< ObservationScalarType, TimeType > > createNewAndKeep(
            const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition ) const;

    //! Create a new dataset containing all observations except those satisfying a condition.
    std::shared_ptr< ObservationDataset< ObservationScalarType, TimeType > > createNewAndDrop(
            const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition ) const;

    //! Mark matching observation rows as rejected without physically deleting data.
    void rejectObservations( const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition,
                             const std::string& reason = "" );

    //! Restore matching observation rows to active status.
    void restoreObservations( const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition );

    //////////////////////////////////////////////////////////
    /////////////////    OBSERVATION VECTOR DATA    //////////
    //////////////////////////////////////////////////////////

    //! Materialize scalar-aligned observation data in canonical vector order.
    ObservationVectorData< ObservationScalarType, TimeType > createObservationVectorData( const bool includeRejected = false ) const;

    //! Materialize observation vector data used by residual/dependent-variable computation.
    ObservationVectorData< ObservationScalarType, TimeType > createComputationObservationVectorData(
            const bool includeRejected = true ) const;

    //! Materialize observation vector data in observable-type/link-end order.
    /*!
     * Row order is observable type, link ends, set index within that
     * group, observation row within the set, and scalar component within the
     * observation. Times and ids are repeated per scalar component, matching
     * concatenated-vector conventions.
     */
    ObservationVectorData< ObservationScalarType, TimeType > createOrderedObservationVectorData( const bool includeInactive = true ) const;

    //! Return set ids in the ordered observable-type/link-ends/index ordering.
    std::vector< unsigned int > getSetIdsInObservationVectorOrder( ) const;

    //! Return the revision used to invalidate structural observation mappings.
    std::size_t getStructuralVersion( ) const;

    //////////////////////////////////////////////////////////
    /////////////////       LEGACY INTERFACES       //////////
    //////////////////////////////////////////////////////////

    //! Replace the legacy concatenated observation vector for one set.
    void setObservationVectorForSet( const unsigned int setId,
                                     const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& observationVector );

    //! Replace the legacy concatenated residual vector for one set.
    void setResidualVectorForSet( const unsigned int setId,
                                  const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& residualVector );

    //! Replace the legacy concatenated diagonal weight vector for one set.
    void setWeightVectorForSet( const unsigned int setId, const Eigen::VectorXd& weightVector );

    //! Set a constant scalar weight on sets selected by a legacy parser.
    void setConstantWeight(
            const double weight = 1.0,
            const std::shared_ptr< ObservationCollectionParser > observationParser = std::make_shared< ObservationCollectionParser >( ) );

    //! Set constant component weights on sets selected by a legacy parser.
    void setConstantWeight(
            const Eigen::VectorXd weight,
            const std::shared_ptr< ObservationCollectionParser > observationParser = std::make_shared< ObservationCollectionParser >( ) );

    //! Set scalar weights for multiple legacy parser selections.
    void setConstantWeightPerObservable(
            const std::map< std::shared_ptr< ObservationCollectionParser >, double > weightsPerObservationParser );

    //! Set component weights for multiple legacy parser selections.
    void setConstantWeightPerObservable(
            const std::map< std::shared_ptr< ObservationCollectionParser >, Eigen::VectorXd > weightsPerObservationParser );

    //! Set a tabulated diagonal weight vector on a legacy parser selection.
    void setTabulatedWeights(
            const Eigen::VectorXd tabulatedWeights,
            const std::shared_ptr< ObservationCollectionParser > observationParser = std::make_shared< ObservationCollectionParser >( ) );

    //! Set tabulated diagonal weights for multiple legacy parser selections.
    void setTabulatedWeights(
            const std::map< std::shared_ptr< ObservationCollectionParser >, Eigen::VectorXd > weightsPerObservationParser );

    //! Return per-set indices accepted by a legacy observation filter.
    std::vector< unsigned int > getFilteredObservationIndices( const unsigned int setId,
                                                               const std::shared_ptr< ObservationFilterBase >& observationFilter ) const;

    //! Move or copy selected observations between sets through the backwards-compatible dataset interface.
    void moveObservationsToSet( const unsigned int sourceSetId,
                                ObservationDataset< ObservationScalarType, TimeType >& targetDataset,
                                const unsigned int targetSetId,
                                const std::vector< unsigned int >& indices,
                                const bool removeFromSource = true );

    //! Remove duplicate time tags from one set through the legacy interface.
    void eraseDuplicateObservationsFromSet( const unsigned int setId, const bool printWarning = true );

    //! Return the legacy concatenated computed-observation vector for one set.
    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getComputedObservationVectorForSet( const unsigned int setId ) const;

    //! Return the legacy concatenated observed-value vector for one set.
    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getObservationVectorForSet( const unsigned int setId ) const;

    //! Return the legacy concatenated diagonal weight vector for one set.
    Eigen::VectorXd getWeightVectorForSet( const unsigned int setId ) const;

    //! Return the legacy concatenated residual vector for one set.
    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getResidualVectorForSet( const unsigned int setId ) const;

    //! Return set identities selected by a legacy parser.
    std::vector< unsigned int > getObservationSetIds( const std::shared_ptr< ObservationCollectionParser >& observationParser =
                                                              std::make_shared< ObservationCollectionParser >( ) ) const;

    //! Return legacy start-and-size entries for all sets in vector order.
    std::vector< std::pair< int, int > > getObservationSetStartAndSize( ) const;

    //! Return legacy set start-and-size entries grouped by observable and link ends.
    std::map< ObservableType, std::map< LinkEnds, std::vector< std::pair< int, int > > > > getObservationSetStartAndSizeByLink( ) const;

    //! Return aggregate legacy start-and-size entries grouped by observable and link ends.
    std::map< ObservableType, std::map< LinkEnds, std::pair< int, int > > > getObservationTypeAndLinkEndStartAndSize( ) const;

    //! Return aggregate legacy start-and-size entries grouped by observable type.
    std::map< ObservableType, std::pair< int, int > > getObservableTypeStartAndSize( ) const;

    //! Return legacy link-end groups for each observable type.
    std::map< ObservableType, std::vector< LinkEnds > > getLinkEndsPerObservableType( ) const;

    //! Return concatenated observed values for one observable type.
    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getObservationVectorForObservableType(
            const ObservableType observableType ) const;

    //! Replace concatenated observed values for one observable type.
    void setObservationVectorForObservableType( const ObservableType observableType,
                                                const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& observations );

    //! Replace concatenated diagonal weights for one observable type.
    void setWeightVectorForObservableType( const ObservableType observableType, const Eigen::VectorXd& weights );

    //! Replace the complete residual vector in observation-vector order.
    void setResidualVector( const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& residualVector );

private:
    friend class cereal::access;

    //! Serialize dataset rows, registries, scalar values, and weights.
    template< class Archive >
    void save( Archive& ar ) const
    {
        ar( observationRows_,
            nextObservationId_,
            setMetadata_,
            observationIdsBySet_,
            linkDefinitionRegistry_,
            ancillarySettingsRegistry_,
            dependentVariableLayoutRegistry_,
            observedValues_,
            residualValues_,
            observationWeights_ );
    }

    //! Deserialize dataset storage and rebuild derived row indexing.
    template< class Archive >
    void load( Archive& ar )
    {
        ar( observationRows_,
            nextObservationId_,
            setMetadata_,
            observationIdsBySet_,
            linkDefinitionRegistry_,
            ancillarySettingsRegistry_,
            dependentVariableLayoutRegistry_,
            observedValues_,
            residualValues_,
            observationWeights_ );
        rebuildRowIndex( );
        ++structuralVersion_;
    }

    //////////////////////////////////////////////////////////
    /////////////////       FACADE ACCESS           //////////
    //////////////////////////////////////////////////////////

    friend struct detail::ObservationSelectionIndices< ObservationScalarType, TimeType >;

    template< typename SetObservationScalarType,
              typename SetTimeType,
              typename std::enable_if< is_state_scalar_and_time_type< SetObservationScalarType, SetTimeType >::value, int >::type >
    friend class SingleObservationSet;

    template< typename CollectionObservationScalarType,
              typename CollectionTimeType,
              typename std::enable_if< is_state_scalar_and_time_type< CollectionObservationScalarType, CollectionTimeType >::value,
                                       int >::type >
    friend class ObservationCollection;

    //! Non-copyable identity value used to detect mappings from another dataset instance.
    struct LifetimeToken {
        //! Create a fresh lifetime identity.
        LifetimeToken( ): value_( std::make_shared< const int >( 0 ) ) {}
        //! Copy construction deliberately creates a fresh identity.
        LifetimeToken( const LifetimeToken& ): LifetimeToken( ) {}
        //! Move construction invalidates both the old and new identities.
        LifetimeToken( LifetimeToken&& other ): LifetimeToken( )
        {
            other.value_ = std::make_shared< const int >( 0 );
        }
        //! Copy assignment replaces the current identity.
        LifetimeToken& operator=( const LifetimeToken& )
        {
            value_ = std::make_shared< const int >( 0 );
            return *this;
        }
        //! Move assignment refreshes both source and destination identities.
        LifetimeToken& operator=( LifetimeToken&& other )
        {
            value_ = std::make_shared< const int >( 0 );
            other.value_ = std::make_shared< const int >( 0 );
            return *this;
        }

        //! Shared marker observed weakly by derived vector-data mappings.
        std::shared_ptr< const int > value_;
    };

    //! Return a weak view of this dataset's current lifetime identity.
    std::weak_ptr< const int > getLifetimeToken( ) const
    {
        return lifetimeToken_.value_;
    }

    //! Replace the registered link definition associated with one set.
    void resetLinkDefinitionForSet( const unsigned int setId, const LinkDefinition& linkDefinition );

    //! Replace dependent-variable bookkeeping associated with one set.
    void resetDependentVariableBookkeepingForSet(
            const unsigned int setId,
            const std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping >& dependentVariableBookkeeping );

    //! Return all observation ids in dataset row order.
    std::vector< unsigned int > getAllObservationIds( ) const;

    //! Return observation ids in observation-vector set order and row order within each set.
    std::vector< unsigned int > getObservationIdsInObservationVectorOrder( ) const;

    //! Materialize observation vector data for selected observation ids.
    ObservationVectorData< ObservationScalarType, TimeType > createObservationVectorDataFromObservationIds(
            const std::vector< unsigned int >& selectedObservationIds,
            const bool includeInactive = true ) const;

public:
    //////////////////////////////////////////////////////////
    /////////////////       ORDERED VECTOR HELPERS  //////////
    //////////////////////////////////////////////////////////

    //! Return ordered set ids for a single observable type.
    std::vector< unsigned int > getObservationSetIdsForObservableType( const ObservableType observableType ) const;

    //! Return total scalar size for all sets of one observable type.
    std::size_t getTotalScalarSizeForObservableType( const ObservableType observableType ) const;

    //! Return concatenated observations for a single observable/link-definition pair.
    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getSingleLinkObservations( const ObservableType observableType,
                                                                                         const LinkDefinition& linkDefinition ) const;

    //! Return scalar-aligned times for a single observable/link-definition pair.
    std::vector< TimeType > getSingleLinkTimes( const ObservableType observableType, const LinkDefinition& linkDefinition ) const;

    //! Return concatenated observations and scalar-aligned times for a single observable/link-definition pair.
    std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, std::vector< TimeType > > getSingleLinkObservationsAndTimes(
            const ObservableType observableType,
            const LinkDefinition& linkDefinition ) const;

    //! Return set ids whose metadata can support the requested dependent-variable settings.
    std::vector< unsigned int > getObservationSetIdsForDependentVariableSettings(
            const std::shared_ptr< simulation_setup::ObservationDependentVariableSettings >& dependentVariableSettings ) const;

    //! Return set ids that already store values for the requested dependent-variable settings.
    std::vector< unsigned int > getObservationSetIdsWithDependentVariableValues(
            const std::shared_ptr< simulation_setup::ObservationDependentVariableSettings >& dependentVariableSettings ) const;

    //! Return compatible dependent-variable settings grouped by observation set.
    std::vector< std::vector< std::shared_ptr< simulation_setup::ObservationDependentVariableSettings > > >
    getCompatibleDependentVariableSettingsPerSet(
            const std::shared_ptr< simulation_setup::ObservationDependentVariableSettings >& dependentVariableSettings ) const;

    //! Return all compatible dependent-variable values grouped by observation set.
    std::vector< std::vector< Eigen::MatrixXd > > getAllCompatibleDependentVariablesPerSet(
            const std::shared_ptr< simulation_setup::ObservationDependentVariableSettings >& dependentVariableSettings ) const;

    //! Return concatenated time history for one dependent-variable setting.
    std::map< TimeType, Eigen::VectorXd > getDependentVariableHistory(
            const std::shared_ptr< simulation_setup::ObservationDependentVariableSettings >& dependentVariableSettings,
            const bool returnFirstCompatibleSettings = false ) const;

    //! Add dependent-variable settings to all compatible sets matching an optional condition.
    void addDependentVariableToSets(
            const std::shared_ptr< simulation_setup::ObservationDependentVariableSettings >& dependentVariableSettings,
            const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition =
                    ObservationSelectionCondition< ObservationScalarType, TimeType >::all( ) );

    //! Fail before using vector data from another dataset or an invalidated mapping.
    void validateObservationVectorData( const ObservationVectorData< ObservationScalarType, TimeType >& observationVectorData ) const
    {
        if( observationVectorData.source_.lock( ) != lifetimeToken_.value_ ||
            observationVectorData.structuralVersion_ != structuralVersion_ ||
            observationVectorData.vectorDataVersion_ != vectorDataVersion_ )
        {
            throw std::runtime_error( "Observation vector data belong to another dataset or have been invalidated." );
        }
        // Legacy APIs expose mutable ancillary pointers. Check their values because
        // changes through such pointers cannot increment the dataset revision.
        for( const auto& entry : observationVectorData.ancillaryBySet_ )
        {
            const auto& current = getAncillarySettingsForSet( entry.first );
            if( static_cast< bool >( current ) != static_cast< bool >( entry.second ) || ( current && !( *current == *entry.second ) ) )
            {
                throw std::runtime_error( "Observation vector data ancillary settings have changed." );
            }
        }
    }

    //! Set residuals for scalar rows described by observation vector data.
    void setResidualVector( const ObservationVectorData< ObservationScalarType, TimeType >& observationVectorData,
                            const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& residualVector );

private:
    //////////////////////////////////////////////////////////
    /////////////////       PRIVATE HELPERS         //////////
    //////////////////////////////////////////////////////////

    //! Record the correlations introduced by a nonzero entry between firstScalar and secondScalar in observationWeights_.
    //! Update both affected sets without reducing either set's existing within-set block size.
    void promoteWeightStructureForEntry( const unsigned int firstScalar, const unsigned int secondScalar );

    //! Recheck correlations after changing set membership, removing observations, or replacing cross-set weights.
    //! Preserve known within-set block sizes; clear cross-set flags only where no connecting entries remain.
    void refreshWeightStructures( );

    //! Single selection/order route shared by inspection and numerical preparation.
    std::vector< unsigned int > resolveObservationIds( const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition,
                                                       const ObservationOrdering ordering,
                                                       const bool includeRejected = true ) const;

    //! Compact all selected rows/scalars together, preserving event and metadata identities.
    void retainObservationRows( const std::vector< unsigned int >& retainedIds );
    //! Rebuild the derived observation-id to row-position lookup.
    void rebuildRowIndex( );
    //! Sort one set's observation identities by observation time.
    void sortObservationIdsForSet( const unsigned int setId );
    //! Return mutable access to one event row by stable identity.
    ObservationDatasetRow< TimeType >& mutableObservationRow( const unsigned int id )
    {
        return observationRows_.at( rowPositionById_.at( id ) );
    }

    //! Validate per-observation vectors before replacing/appending set data.
    void validateObservationSetData( const unsigned int setId,
                                     const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
                                     const std::vector< TimeType >& times,
                                     const std::vector< Eigen::VectorXd >& dependentVariables,
                                     const std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >& weights,
                                     const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& residuals ) const;

    //! Return a permutation that sorts observation events by their time tags.
    static std::vector< std::size_t > getTimeSortingPermutation( const std::vector< TimeType >& observationTimes );

    //! Apply a precomputed observation-event permutation to one parallel vector.
    template< typename VectorType >
    static void reorderVector( std::vector< VectorType >& data, const std::vector< std::size_t >& permutation );

    //! Remove observation-event entries from a parallel vector by index.
    template< typename VectorType >
    static void removeEntries( std::vector< VectorType >& data, const std::vector< unsigned int >& indicesToRemove );

    //! Gather scalar values for one set into a legacy flat vector.
    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > createSetVector(
            const unsigned int setId,
            const std::vector< ObservationScalarType >& scalarValues ) const;

    //! Overwrite all scalar components belonging to one observation event.
    void setObservationValue( const unsigned int observationId,
                              const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& observation );

    //! Overwrite all residual scalar components belonging to one observation event.
    void setResidualValue( const unsigned int observationId, const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& residual );

    //! Overwrite all weight scalar components belonging to one observation event.
    void setWeightValue( const unsigned int observationId, const Eigen::VectorXd& weight );

    //! Expand observation ids and optional component indices to scalar-component ids.
    std::vector< unsigned int > getScalarComponentIdsForObservationSelection( const std::vector< unsigned int >& observationIds,
                                                                              const std::vector< unsigned int >& components ) const;

    //! Return the registry id for a link definition, inserting it if needed.
    int registerLinkDefinition( const LinkDefinition& linkDefinition );

    //! Return the registry id for ancillary settings, inserting them if needed.
    int registerAncillarySettings( const std::shared_ptr< ObservationAncillarySimulationSettings >& ancillarySettings );

    //! Return the registry id for dependent-variable bookkeeping, inserting it if needed.
    int registerDependentVariableLayout( const std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping >& bookkeeping );

    //////////////////////////////////////////////////////////
    /////////////////       MEMBER VARIABLES        //////////
    //////////////////////////////////////////////////////////

    //! One row per observation event; vector observables occupy one row, not N rows.
    std::vector< ObservationDatasetRow< TimeType > > observationRows_;

    //! For each ObservationDatasetRow::observationId_, store its current index in observationRows_.
    //! The observation ID remains unchanged, while the index may change when observations are reordered or removed.
    std::unordered_map< unsigned int, std::size_t > rowPositionById_;

    //! Next stable identity assigned to a newly added observation event.
    unsigned int nextObservationId_ = 0;

    //!  For each set id, the metadata record belonging to that set.
    std::vector< ObservationSetMetadata< ObservationScalarType, TimeType > > setMetadata_;

    //! ObservationDatasetRow::observationId_ values grouped and ordered within each observation set.
    //! The outer vector index is the setId; each inner vector is ordered by ObservationDatasetRow::indexInSet_.
    std::vector< std::vector< unsigned int > > observationIdsBySet_;

    //! Full link definitions referenced by the observation sets.
    //! For a given setId, setMetadata_.at(setId).linkDefinitionId_ gives the index in this vector.
    std::vector< LinkDefinition > linkDefinitionRegistry_;

    //! Ancillary settings referenced by the observation sets.
    //! For a given setId, setMetadata_.at(setId).ancillarySettingsId_ gives the index in this vector.
    std::vector< std::shared_ptr< ObservationAncillarySimulationSettings > > ancillarySettingsRegistry_;

    //! Dependent-variable bookkeeping referenced by the observation sets.
    //! For a given setId, setMetadata_.at(setId).dependentVariableLayoutId_ gives the index in this vector.
    std::vector< std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping > > dependentVariableLayoutRegistry_;

    //! Scalar observed values; vector observables contribute one entry per component.
    std::vector< ObservationScalarType > observedValues_;

    //! Scalar residual values aligned one-to-one with observedValues_.
    std::vector< ObservationScalarType > residualValues_;

    //! Compact observation weight storage; materialized into vectors/matrices only on request.
    ObservationWeights observationWeights_;

    //! Version number for the dataset structure recorded in ObservationVectorData.
    //! It increases when observation rows or set metadata change; writeback is rejected if the stored number differs.
    std::size_t structuralVersion_ = 0;

    //! Version number for the observed values and each row's active or inactive state.
    //! It increases when either changes; ObservationVectorData created before that change cannot be written back.
    std::size_t vectorDataVersion_ = 0;

    //! Identifier for this specific ObservationDataset object.
    //! ObservationVectorData and compatibility wrappers retain it to detect data from another or already-destroyed dataset.
    LifetimeToken lifetimeToken_;
};
}  // namespace observation_models

}  // namespace tudat

#include "tudat/simulation/estimation_setup/observationLegacyParserAdapter.h"
#include "tudat/simulation/estimation_setup/observationDatasetInspectionImplementation.h"

#include "tudat/simulation/estimation_setup/observationDatasetMutationImplementation.h"
#include "tudat/simulation/estimation_setup/observationDatasetWeightsImplementation.h"
#include "tudat/simulation/estimation_setup/observationDatasetAccessorsImplementation.h"
#include "tudat/simulation/estimation_setup/observationDatasetSelectionImplementation.h"
#include "tudat/simulation/estimation_setup/observationDatasetVectorDataImplementation.h"
#include "tudat/simulation/estimation_setup/observationDatasetPrivateImplementation.h"
#include "tudat/simulation/estimation_setup/observationDatasetLegacyImplementation.h"

#endif  // TUDAT_OBSERVATION_DATASET_H
