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
#include "tudat/basics/tudatTypeTraits.h"

namespace tudat
{

namespace observation_models
{

//! Class holding a set of observations and associated data, for use in an estimation.
/*!
 * Each observation is stored as one row, while its observed values, residuals and weights are stored with one entry
 * per observation component. Metadata shared by all observations in a set is stored once, and referenced by id. The
 * concatenated vectors used in the estimation are created with createObservationVectorData. The functions that
 * retrieve the data of selected observations return copies, in internal or estimation order.
 */
template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
class ObservationDataset
{
public:
    //! Create an empty observation dataset.
    ObservationDataset( ) = default;
    //! Copy constructor; the copy is independent of the original dataset.
    ObservationDataset( const ObservationDataset& other );
    //! Move constructor.
    ObservationDataset( ObservationDataset&& ) = default;
    //! Move-assign an observation dataset.
    ObservationDataset& operator=( ObservationDataset&& ) = default;
    //! Copy assignment operator; the copy is independent of the original dataset.
    ObservationDataset& operator=( const ObservationDataset& other )
    {
        if( this != &other )
        {
            ObservationDataset copy( other );
            *this = std::move( copy );
        }
        return *this;
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

    //! Add an observation set, with weights defined by weight settings.
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

    //! Replace the concatenated diagonal weight vector for one set.
    void setWeightVectorForSet( const unsigned int setId, const Eigen::VectorXd& weightVector );

    //! Return whether the effective set principal block contains off-diagonal weights.
    bool hasWeightMatrixForSet( const unsigned int setId ) const
    {
        return observationWeights_.restricted( getScalarComponentIdsForObservationSelection( observationIdsBySet_.at( setId ), {} ) )
                .hasOffDiagonalWeights( );
    }

    //! Store one observable-size N x N weight block for an observation event.
    void setWeightMatrixForObservation( const unsigned int observationId, const Eigen::MatrixXd& weightMatrix );

    //! Return whether the effective observation block contains off-diagonal weights.
    bool hasWeightMatrixForObservation( const unsigned int observationId ) const
    {
        return observationWeights_.restricted( getScalarComponentIdsForObservationSelection( { observationId }, {} ) )
                .hasOffDiagonalWeights( );
    }

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

    //! Return whether the weights of the dataset contain off-diagonal entries.
    bool hasExtraWeightBlocks( ) const
    {
        return observationWeights_.hasOffDiagonalWeights( );
    }

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

    //! Remove all observation events matching a row-level condition.
    void removeObservations( const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition )
    {
        retainObservationRows( getObservationIdsMatchingCondition( !condition ) );
    }

    //! Physically remove all currently rejected observation events.
    void deleteRejectedObservations( )
    {
        removeObservations( ObservationSelectionCondition< ObservationScalarType, TimeType >::rejected( ) );
    }

    //! Return time bounds for all observations in the dataset.
    std::pair< TimeType, TimeType > getTimeBounds( ) const;

    //////////////////////////////////////////////////////////
    /////////////////       DATA ACCESSORS          //////////
    //////////////////////////////////////////////////////////

    //! Metadata per set id: set metadata, link definition, ancillary settings and dependent-variable layout.
    //! The layout maps the (start index, size) of each dependent variable to a copy of its settings.
    using InspectionMetadata = std::map<
            unsigned int,
            std::tuple< ObservationSetMetadata< ObservationScalarType, TimeType >,
                        LinkDefinition,
                        std::shared_ptr< ObservationAncillarySimulationSettings >,
                        std::map< std::pair< int, int >, std::shared_ptr< simulation_setup::ObservationDependentVariableSettings > > > >;

    //! The following functions return copies of the data of the observations selected by a condition (by default all
    //! observations, including rejected ones), in internal or estimation order. The order does not affect the selection.
    //! Function to retrieve the time of each selected observation.
    std::vector< TimeType > getTimes( const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition =
                                              ObservationSelectionCondition< ObservationScalarType, TimeType >::all( ),
                                      const ObservationOrdering ordering = ObservationOrdering::internal ) const;

    //! Function to retrieve the value of each selected observation.
    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > getObservations(
            const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition =
                    ObservationSelectionCondition< ObservationScalarType, TimeType >::all( ),
            const ObservationOrdering ordering = ObservationOrdering::internal ) const;

    //! Function to retrieve the residual of each selected observation (zero if not computed).
    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > getResiduals(
            const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition =
                    ObservationSelectionCondition< ObservationScalarType, TimeType >::all( ),
            const ObservationOrdering ordering = ObservationOrdering::internal ) const;

    //! Function to retrieve the id of each selected observation.
    std::vector< unsigned int > getObservationIds( const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition =
                                                           ObservationSelectionCondition< ObservationScalarType, TimeType >::all( ),
                                                   const ObservationOrdering ordering = ObservationOrdering::internal ) const
    {
        return resolveObservationIds( condition, ordering );
    }

    //! Function to retrieve the set id of each selected observation.
    std::vector< unsigned int > getSetIds( const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition =
                                                   ObservationSelectionCondition< ObservationScalarType, TimeType >::all( ),
                                           const ObservationOrdering ordering = ObservationOrdering::internal ) const;

    //! Function to retrieve the row (including dependent variables and rejection status) of each selected observation.
    std::vector< ObservationDatasetRow< TimeType > > getRows(
            const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition =
                    ObservationSelectionCondition< ObservationScalarType, TimeType >::all( ),
            const ObservationOrdering ordering = ObservationOrdering::internal ) const;

    //! Function to retrieve the dependent variables of each selected observation (empty if none are stored).
    std::vector< Eigen::VectorXd > getDependentVariableValues(
            const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition =
                    ObservationSelectionCondition< ObservationScalarType, TimeType >::all( ),
            const ObservationOrdering ordering = ObservationOrdering::internal ) const;

    //! Function to retrieve the (observation id, component index) of each component of the selected observations.
    std::vector< std::pair< unsigned int, unsigned int > > getScalarComponents(
            const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition =
                    ObservationSelectionCondition< ObservationScalarType, TimeType >::all( ),
            const ObservationOrdering ordering = ObservationOrdering::internal ) const;

    //! Function to retrieve the diagonal of the weight matrix of the selected observations, ordered as getScalarComponents.
    Eigen::VectorXd getWeightDiagonal( const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition =
                                               ObservationSelectionCondition< ObservationScalarType, TimeType >::all( ),
                                       const ObservationOrdering ordering = ObservationOrdering::internal ) const
    {
        return observationWeights_.getDiagonal(
                getScalarComponentIdsForObservationSelection( resolveObservationIds( condition, ordering ), {} ) );
    }

    //! Function to retrieve the weight matrix of the selected observations, with rows and columns ordered as getScalarComponents.
    Eigen::SparseMatrix< double > getWeightMatrix( const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition =
                                                           ObservationSelectionCondition< ObservationScalarType, TimeType >::all( ),
                                                   const ObservationOrdering ordering = ObservationOrdering::internal ) const
    {
        return observationWeights_
                .restricted( getScalarComponentIdsForObservationSelection( resolveObservationIds( condition, ordering ), {} ) )
                .sparseMatrix( );
    }

    //! Function to retrieve the metadata of the sets containing the selected observations, with the set id as key.
    InspectionMetadata getMetadata( const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition =
                                            ObservationSelectionCondition< ObservationScalarType, TimeType >::all( ),
                                    const ObservationOrdering ordering = ObservationOrdering::internal ) const;

    //! Return the number of logical observation sets.
    std::size_t getNumberOfObservationSets( ) const
    {
        return setMetadata_.size( );
    }

    //! Return the number of observation events.
    std::size_t getNumberOfObservations( ) const
    {
        return observationRows_.size( );
    }

    //! Return the total number of scalar observation components.
    std::size_t getTotalScalarSize( ) const
    {
        return observedValues_.size( );
    }

    //! Return metadata for every observation set.
    const std::vector< ObservationSetMetadata< ObservationScalarType, TimeType > >& getObservationSetMetadata( ) const
    {
        return setMetadata_;
    }

    //! Return metadata for one observation set.
    const ObservationSetMetadata< ObservationScalarType, TimeType >& getObservationSetMetadata( const unsigned int setId ) const
    {
        return setMetadata_.at( setId );
    }

    //! Return all event rows in internal storage order.
    const std::vector< ObservationDatasetRow< TimeType > >& getObservationRows( ) const
    {
        return observationRows_;
    }

    //! Return the row of one observation.
    const ObservationDatasetRow< TimeType >& getObservationRow( const unsigned int observationId ) const
    {
        return observationRows_.at( rowPositionById_.at( observationId ) );
    }

    //! Return the observation id and component index of one entry in the component storage.
    ObservationScalarComponentRow getScalarComponentRow( const unsigned int scalarComponentId ) const;

    //! Return the ids of the observations in one set.
    const std::vector< unsigned int >& getObservationIdsForSet( const unsigned int setId ) const
    {
        return observationIdsBySet_.at( setId );
    }

    //! Store user-defined numerical values for a list of observations, under a given key.
    void addNumericalObservationMetadata( const std::string& key,
                                          const std::vector< unsigned int >& observationIds,
                                          const std::vector< double >& values );

    //! Return metadata values in the supplied observation-ID order; missing values throw.
    std::vector< double > getNumericalObservationMetadata( const std::string& key,
                                                           const std::vector< unsigned int >& observationIds ) const;

    //! Return all named numerical metadata, including values belonging to rejected observations.
    const std::map< std::string, std::map< unsigned int, double > >& getNumericalObservationMetadata( ) const
    {
        return numericalObservationMetadata_;
    }

    //! Return one vector-valued measurement per observation event in a set.
    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > getObservationsForSet( const unsigned int setId ) const;

    //! Return the concatenated observed-value vector for one set.
    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getObservationVectorForSet( const unsigned int setId ) const
    {
        return createSetVector( setId, observedValues_ );
    }

    //! Reconstruct the vector-valued measurement for one observation event.
    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getObservationValue( const unsigned int observationId ) const;

    //! Return one reference-link-end time per observation event in a set.
    std::vector< TimeType > getObservationTimesForSet( const unsigned int setId ) const;

    //! Return the reference-link-end time for one observation event.
    TimeType getObservationTime( const unsigned int observationId ) const
    {
        return getObservationRow( observationId ).time_;
    }

    //! Return one vector of scalar-component weights per observation event.
    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > getWeightsForSet( const unsigned int setId ) const;

    //! Reconstruct the scalar-component weight vector for one observation event.
    Eigen::Matrix< double, Eigen::Dynamic, 1 > getWeightValue( const unsigned int observationId ) const
    {
        return observationWeights_.getDiagonal( getScalarComponentIdsForObservationSelection( { observationId }, {} ) );
    }

    //! Reconstruct the weight matrix for one observation event.
    Eigen::MatrixXd getWeightMatrixForObservation( const unsigned int observationId ) const
    {
        return observationWeights_.restricted( getScalarComponentIdsForObservationSelection( { observationId }, {} ) )
                .sparseMatrix( )
                .toDense( );
    }

    //! Return the weight matrix of one set, including correlations.
    Eigen::MatrixXd getWeightMatrixForSet( const unsigned int setId ) const
    {
        return observationWeights_.restricted( getScalarComponentIdsForObservationSelection( observationIdsBySet_.at( setId ), {} ) )
                .sparseMatrix( )
                .toDense( );
    }

    //! Return one vector-valued residual per observation event in a set.
    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > getResidualsForSet( const unsigned int setId ) const;

    //! Return the concatenated residual vector for one set.
    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getResidualVectorForSet( const unsigned int setId ) const
    {
        return createSetVector( setId, residualValues_ );
    }

    //! Reconstruct the vector-valued residual for one observation event.
    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getResidualValue( const unsigned int observationId ) const;

    //////////////////////////////////////////////////////////
    /////////////////       DEPENDENT VARIABLES    ///////////
    //////////////////////////////////////////////////////////

    //! Return per-observation dependent-variable vectors for one set.
    std::vector< Eigen::VectorXd > getDependentVariablesForSet( const unsigned int setId ) const;

    //! Return the dependent-variable vector for one observation event.
    Eigen::VectorXd getDependentVariables( const unsigned int observationId ) const
    {
        return getObservationRow( observationId ).dependentVariableValues_;
    }

    //! Replace per-observation dependent-variable vectors for one set.
    void setDependentVariablesForSet( const unsigned int setId, const std::vector< Eigen::VectorXd >& dependentVariables );

    //! Clear dependent-variable values for every observation event in one set.
    void clearDependentVariablesForSet( const unsigned int setId );

    //! Add dependent-variable settings to all compatible sets matching an optional condition.
    void addDependentVariableToSets(
            const std::shared_ptr< simulation_setup::ObservationDependentVariableSettings >& dependentVariableSettings,
            const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition =
                    ObservationSelectionCondition< ObservationScalarType, TimeType >::all( ) );

    //! Return the number of observation events in one set.
    std::size_t getNumberOfObservationsForSet( const unsigned int setId ) const
    {
        return observationIdsBySet_.at( setId ).size( );
    }

    //! Return the scalar observation size of one set.
    std::size_t getTotalScalarSizeForSet( const unsigned int setId ) const
    {
        return observationIdsBySet_.at( setId ).size( ) * getObservationSetMetadata( setId ).observableSize_;
    }

    //! Return a registered link definition by registry identity.
    const LinkDefinition& getLinkDefinition( const unsigned int linkDefinitionId ) const
    {
        return linkDefinitionRegistry_.at( linkDefinitionId );
    }

    //! Replace a named reference point for matching link ends in selected observation sets.
    void setLinkEndReferencePoint( const std::string& bodyName,
                                   const std::string& referencePointName,
                                   const LinkEndType linkEndType,
                                   const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition =
                                           ObservationSelectionCondition< ObservationScalarType, TimeType >::all( ) );

    //! Return registered ancillary settings by registry identity.
    const std::shared_ptr< ObservationAncillarySimulationSettings >& getAncillarySettings( const unsigned int ancillarySettingsId ) const
    {
        return ancillarySettingsRegistry_.at( ancillarySettingsId );
    }

    //! Return ancillary settings associated with one observation set.
    const std::shared_ptr< ObservationAncillarySimulationSettings >& getAncillarySettingsForSet( const unsigned int setId ) const
    {
        return getAncillarySettings( getObservationSetMetadata( setId ).ancillarySettingsId_ );
    }

    //! Return registered dependent-variable bookkeeping by registry identity.
    const std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping >& getDependentVariableBookkeeping(
            const unsigned int dependentVariableLayoutId ) const
    {
        return dependentVariableLayoutRegistry_.at( dependentVariableLayoutId );
    }

    //////////////////////////////////////////////////////////
    /////////////////       SELECTION AND STATUS    //////////
    //////////////////////////////////////////////////////////

    //! Return the ids of the observations that satisfy a condition.
    std::vector< unsigned int > getObservationIdsMatchingCondition(
            const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition ) const;

    //! Create a new dataset containing only observations satisfying a condition.
    std::shared_ptr< ObservationDataset< ObservationScalarType, TimeType > > createNewAndKeep(
            const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition ) const;

    //! Create a new dataset containing all observations except those satisfying a condition.
    std::shared_ptr< ObservationDataset< ObservationScalarType, TimeType > > createNewAndDrop(
            const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition ) const
    {
        return createNewAndKeep( !condition );
    }

    //! Mark matching observation rows as rejected without physically deleting data.
    void rejectObservations( const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition,
                             const std::string& reason = "" );

    //! Reject active selected rows if any residual exceeds sigmaLimit times its measurement uncertainty.
    //! Covariance includes the complete correlated set; residuals and weights remain unchanged.
    void rejectObservations( const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition,
                             const double sigmaLimit,
                             const std::string& reason = "" );

    //! Restore matching observation rows to active status and clear their rejection reason.
    void restoreObservations( const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition );

    //////////////////////////////////////////////////////////
    /////////////////    OBSERVATION VECTOR DATA    //////////
    //////////////////////////////////////////////////////////

    //! Create the concatenated observation data, in the order used in the estimation.
    ObservationVectorData< ObservationScalarType, TimeType > createObservationVectorData( const bool includeRejected = false ) const;

    //! Return set ids in the ordered observable-type/link-ends/index ordering.
    std::vector< unsigned int > getSetIdsInObservationVectorOrder( ) const;

    //! Check that observation vector data was created from this dataset, and that the dataset has not changed since.
    void validateObservationVectorData( const ObservationVectorData< ObservationScalarType, TimeType >& observationVectorData ) const
    {
        if( observationVectorData.source_.lock( ) != lifetimeToken_ || observationVectorData.structuralVersion_ != structuralVersion_ ||
            observationVectorData.vectorDataVersion_ != vectorDataVersion_ )
        {
            throw std::runtime_error( "Observation vector data belong to another dataset or have been invalidated." );
        }
        // Ancillary settings are accessible through mutable pointers. Check their values because
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

    //! Replace the registered link definition associated with one set.
    void resetLinkDefinitionForSet( const unsigned int setId, const LinkDefinition& linkDefinition )
    {
        setMetadata_.at( setId ).linkDefinitionId_ = registerEntry( linkDefinitionRegistry_, linkDefinition );
        ++structuralVersion_;
    }

    //! Replace dependent-variable bookkeeping associated with one set.
    void resetDependentVariableBookkeepingForSet(
            const unsigned int setId,
            const std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping >& dependentVariableBookkeeping );

    //! Return all observation ids in dataset row order.
    std::vector< unsigned int > getAllObservationIds( ) const;

    //! Return observation ids in observation-vector set order and row order within each set.
    std::vector< unsigned int > getObservationIdsInObservationVectorOrder( ) const;

    //! Record the correlations introduced by a nonzero entry between firstScalar and secondScalar in observationWeights_.
    //! Update both affected sets without reducing either set's existing within-set block size.
    void promoteWeightStructureForEntry( const unsigned int firstScalar, const unsigned int secondScalar );

    //! Recheck correlations after changing set membership, removing observations, or replacing cross-set weights.
    //! Preserve known within-set block sizes; clear cross-set flags only where no connecting entries remain.
    void refreshWeightStructures( );

    //! Return the ids of the observations selected by a condition, in the requested order.
    std::vector< unsigned int > resolveObservationIds( const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition,
                                                       const ObservationOrdering ordering,
                                                       const bool includeRejected = true ) const;

    //! Remove all observations that are not in retainedIds, keeping the ids of the remaining observations and sets.
    void retainObservationRows( const std::vector< unsigned int >& retainedIds );
    //! Sort one set's observation identities by observation time.
    void sortObservationIdsForSet( const unsigned int setId );
    //! Return a modifiable reference to the row of one observation.
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

    //! Gather scalar values for one set into a flat vector.
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

    //! Return the index of the first registry entry equal to entry, appending entry if none is equal.
    template< typename EntryType >
    static unsigned int registerEntry( std::vector< EntryType >& registry, const EntryType& entry );

    //////////////////////////////////////////////////////////
    /////////////////       MEMBER VARIABLES        //////////
    //////////////////////////////////////////////////////////

    //! One row per observation event; vector observables occupy one row, not N rows.
    std::vector< ObservationDatasetRow< TimeType > > observationRows_;

    //! For each ObservationDatasetRow::observationId_, store its current index in observationRows_.
    //! The observation ID remains unchanged, while the index may change when observations are reordered or removed.
    std::unordered_map< unsigned int, std::size_t > rowPositionById_;

    //! Id that is assigned to the next observation that is added.
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

    //! Weights of all observation components, including off-diagonal entries.
    ObservationWeights observationWeights_;

    //! User-defined numerical metadata, with the name and observation id as keys.
    std::map< std::string, std::map< unsigned int, double > > numericalObservationMetadata_;

    //! Copy of the complete weights used by ObservationVectorData; recreated when the weights or observations change.
    mutable std::shared_ptr< const typename ObservationVectorData< ObservationScalarType, TimeType >::CompleteWeightData >
            completeWeightData_;

    //! Version number for the dataset structure recorded in ObservationVectorData.
    //! It increases when observation rows or set metadata change; residuals cannot be written back from ObservationVectorData with a different number.
    std::size_t structuralVersion_ = 0;

    //! Version number for the observed values and each row's active or inactive state.
    //! It increases when either changes; ObservationVectorData created before that change cannot be written back.
    std::size_t vectorDataVersion_ = 0;

    //! Identifier for this specific ObservationDataset object; copies receive a new one and moves transfer it.
    //! ObservationVectorData retain it weakly to detect data from another or already-destroyed dataset.
    std::shared_ptr< const int > lifetimeToken_ = std::make_shared< const int >( 0 );
};
}  // namespace observation_models

}  // namespace tudat

#include "tudat/simulation/estimation_setup/observationDatasetInspectionImplementation.h"

#include "tudat/simulation/estimation_setup/observationDatasetMutationImplementation.h"
#include "tudat/simulation/estimation_setup/observationDatasetWeightsImplementation.h"
#include "tudat/simulation/estimation_setup/observationDatasetAccessorsImplementation.h"
#include "tudat/simulation/estimation_setup/observationDatasetSelectionImplementation.h"
#include "tudat/simulation/estimation_setup/observationDatasetVectorDataImplementation.h"
#include "tudat/simulation/estimation_setup/observationDatasetPrivateImplementation.h"

#endif  // TUDAT_OBSERVATION_DATASET_H
