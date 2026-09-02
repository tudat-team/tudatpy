#ifndef TUDAT_SINGLE_OBSERVATION_SET_FILTERING_AND_APPENDING_H
#define TUDAT_SINGLE_OBSERVATION_SET_FILTERING_AND_APPENDING_H

#include "tudat/simulation/estimation_setup/singleObservationSet.h"

namespace tudat
{

namespace observation_models
{

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type
                  IsStateScalarAndTimeType >
//! Applies an observation filter by moving matching entries between active and filtered storage.
void SingleObservationSet< ObservationScalarType, TimeType, IsStateScalarAndTimeType >::filterObservations(
        const std::shared_ptr< ObservationFilterBase > observationFilter,
        const bool saveFilteredObservations )
{
    // Lazily create filtered storage when filtering-out is requested for the first time.
    if( observationFilter->filterOut( ) && filteredObservationSet_ == nullptr )
    {
        // Initialise empty filtered observation set
        filteredObservationSet_ = std::make_shared< SingleObservationSet< ObservationScalarType, TimeType > >(
                observableType_,
                linkEnds_,
                std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >( ),
                std::vector< TimeType >( ),
                referenceLinkEnd_,
                std::vector< Eigen::VectorXd >( ),
                dependentVariableBookkeeping_,
                ancillarySettings_ );
    }
    if( !observationFilter->filterOut( ) && filteredObservationSet_ == nullptr )
    {
        throw std::runtime_error(
                "Error when attempting to un-filter observations, filtered observation set is "
                "empty." );
    }

    unsigned int nbObservationsToTest =
            ( observationFilter->filterOut( ) ? numberOfObservations_ : filteredObservationSet_->getNumberOfObservables( ) );
    bool useOppositeCondition = observationFilter->useOppositeCondition( );

    const bool filterOutObservations = observationFilter->filterOut( );
    const auto getObservationTimeAt = [ & ]( const unsigned int index ) -> TimeType {
        return filterOutObservations ? observationTimes_.at( index ) : filteredObservationSet_->getObservationTime( index );
    };

    // Build a predicate that marks observations to be moved for the configured filter type.
    std::function< bool( unsigned int ) > shouldRemoveObservation;
    switch( observationFilter->getFilterType( ) )
    {
        case residual_filtering: {
            Eigen::VectorXd residualCutOff = Eigen::VectorXd::Zero( singleObservationSize_ );
            const std::shared_ptr< ObservationFilter< double > > scalarFilter =
                    std::dynamic_pointer_cast< ObservationFilter< double > >( observationFilter );
            const std::shared_ptr< ObservationFilter< Eigen::VectorXd > > vectorFilter =
                    std::dynamic_pointer_cast< ObservationFilter< Eigen::VectorXd > >( observationFilter );
            if( scalarFilter != nullptr )
            {
                residualCutOff = scalarFilter->getFilterValue( ) * Eigen::VectorXd::Ones( singleObservationSize_ );
            }
            else if( vectorFilter != nullptr )
            {
                if( vectorFilter->getFilterValue( ).size( ) != singleObservationSize_ )
                {
                    throw std::runtime_error(
                            "Error when performing residual filtering, size of the residual "
                            "cut off vector inconsistent with observable size." );
                }
                residualCutOff = vectorFilter->getFilterValue( );
            }

            shouldRemoveObservation = [ &, residualCutOff ]( const unsigned int index ) {
                const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > residualVector =
                        filterOutObservations ? residuals_.at( index ) : filteredObservationSet_->getResidual( index );
                for( unsigned int k = 0; k < singleObservationSize_; ++k )
                {
                    const bool passesFilter = std::fabs( residualVector[ k ] ) > residualCutOff[ k ];
                    if( ( !useOppositeCondition && passesFilter ) || ( useOppositeCondition && !passesFilter ) )
                    {
                        return true;
                    }
                }
                return false;
            };
            break;
        }
        case absolute_value_filtering: {
            Eigen::VectorXd absoluteValueCutOff = Eigen::VectorXd::Zero( singleObservationSize_ );
            const std::shared_ptr< ObservationFilter< double > > scalarFilter =
                    std::dynamic_pointer_cast< ObservationFilter< double > >( observationFilter );
            const std::shared_ptr< ObservationFilter< Eigen::VectorXd > > vectorFilter =
                    std::dynamic_pointer_cast< ObservationFilter< Eigen::VectorXd > >( observationFilter );
            if( scalarFilter != nullptr )
            {
                absoluteValueCutOff = scalarFilter->getFilterValue( ) * Eigen::VectorXd::Ones( singleObservationSize_ );
            }
            else if( vectorFilter != nullptr )
            {
                if( vectorFilter->getFilterValue( ).size( ) != singleObservationSize_ )
                {
                    throw std::runtime_error(
                            "Error when performing observation value filtering, size of the "
                            "filter value inconsistent with observable size." );
                }
                absoluteValueCutOff = vectorFilter->getFilterValue( );
            }

            shouldRemoveObservation = [ &, absoluteValueCutOff ]( const unsigned int index ) {
                const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > observationVector =
                        filterOutObservations ? observations_.at( index ) : filteredObservationSet_->getObservation( index );
                for( unsigned int k = 0; k < singleObservationSize_; ++k )
                {
                    const bool passesFilter = observationVector[ k ] > absoluteValueCutOff[ k ];
                    if( ( !useOppositeCondition && passesFilter ) || ( useOppositeCondition && !passesFilter ) )
                    {
                        return true;
                    }
                }
                return false;
            };
            break;
        }
        case epochs_filtering: {
            const std::vector< double > filterEpochs =
                    std::dynamic_pointer_cast< ObservationFilter< std::vector< double > > >( observationFilter )->getFilterValue( );
            shouldRemoveObservation = [ &, filterEpochs ]( const unsigned int index ) {
                const TimeType observationTime = getObservationTimeAt( index );
                const bool passesFilter = std::count( filterEpochs.begin( ), filterEpochs.end( ), observationTime ) > 0;
                return ( !useOppositeCondition && passesFilter ) || ( useOppositeCondition && !passesFilter );
            };
            break;
        }
        case time_bounds_filtering: {
            const std::pair< double, double > timeBounds =
                    std::dynamic_pointer_cast< ObservationFilter< std::pair< double, double > > >( observationFilter )->getFilterValue( );
            shouldRemoveObservation = [ &, timeBounds ]( const unsigned int index ) {
                const TimeType observationTime = getObservationTimeAt( index );
                const bool passesFilter = ( observationTime >= timeBounds.first ) && ( observationTime <= timeBounds.second );
                return ( !useOppositeCondition && passesFilter ) || ( useOppositeCondition && !passesFilter );
            };
            break;
        }
        case dependent_variable_filtering: {
            const std::shared_ptr< ObservationDependentVariableFilter > dependentVariableFilter =
                    std::dynamic_pointer_cast< ObservationDependentVariableFilter >( observationFilter );
            if( dependentVariableFilter == nullptr )
            {
                throw std::runtime_error(
                        "Error when performing dependent variable filtering, inconsistent "
                        "filter input (should be ObservationDependentVariableFilter object)." );
            }

            const std::shared_ptr< ObservationDependentVariableSettings > settings =
                    dependentVariableFilter->getDependentVariableSettings( );
            const unsigned int dependentVariableSize = getObservationDependentVariableSize( settings, linkEnds_.linkEnds_ );

            const Eigen::VectorXd dependentVariableCutOff = dependentVariableFilter->getFilterValue( );
            if( dependentVariableCutOff.size( ) != dependentVariableSize )
            {
                throw std::runtime_error(
                        "Error when performing dependent variable filtering, size of the "
                        "dependent variable cut off vector inconsistent with dependent "
                        "variable size." );
            }

            const Eigen::MatrixXd singleDependentVariableValues = filterOutObservations
                    ? getSingleDependentVariable( settings )
                    : filteredObservationSet_->getSingleDependentVariable( settings );
            if( ( singleDependentVariableValues.rows( ) != nbObservationsToTest ) ||
                ( singleDependentVariableValues.cols( ) != dependentVariableSize ) )
            {
                throw std::runtime_error(
                        "Error when performing dependent variable filtering, size of "
                        "observation dependent variables is inconsistent with the "
                        "number of observations and presupposed dependent variable size." );
            }

            shouldRemoveObservation =
                    [ &, dependentVariableCutOff, singleDependentVariableValues, dependentVariableSize ]( const unsigned int index ) {
                        for( unsigned int k = 0; k < dependentVariableSize; ++k )
                        {
                            const bool passesFilter = singleDependentVariableValues( index, k ) > dependentVariableCutOff[ k ];
                            if( ( !useOppositeCondition && passesFilter ) || ( useOppositeCondition && !passesFilter ) )
                            {
                                return true;
                            }
                        }
                        return false;
                    };
            break;
        }
        default:
            throw std::runtime_error( "Observation filter type not recognised." );
    }

    std::vector< unsigned int > indicesToRemove;
    indicesToRemove.reserve( nbObservationsToTest );
    // Evaluate the predicate over the currently active/filter-targeted observations.
    for( unsigned int j = 0; j < nbObservationsToTest; ++j )
    {
        if( shouldRemoveObservation( j ) )
        {
            indicesToRemove.push_back( j );
        }
    }

    // Move selected entries to filtered storage or back to active storage.
    if( observationFilter->filterOut( ) )
    {
        moveObservationsInOutFilteredSet( indicesToRemove, true, saveFilteredObservations );
    }
    else
    {
        moveObservationsInOutFilteredSet( indicesToRemove, false, true );
    }
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type
                  IsStateScalarAndTimeType >
//! Appends observations with diagonal weights, enforcing consistency with the current weight model.
void SingleObservationSet< ObservationScalarType, TimeType, IsStateScalarAndTimeType >::addObservations(
        const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
        const std::vector< TimeType >& times,
        const std::vector< Eigen::VectorXd >& dependentVariables,
        const std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >& weights,
        const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& residuals,
        const bool sortObservations )
{
    validateOptionalBatchSize( weights, observations.size( ), "weights", "adding observations to SingleObservationSet" );

    // Full-matrix correlations cannot be safely expanded by this append API.
    if( weightState_.fullWeights.rows( ) > 0 || fullWeightCrossCorrelationWithFilteredSet_.rows( ) > 0 ||
        ( filteredObservationSet_ != nullptr &&
          ( filteredObservationSet_->weightState_.fullWeights.rows( ) > 0 ||
            filteredObservationSet_->fullWeightCrossCorrelationWithFilteredSet_.rows( ) > 0 ) ) )
    {
        throw std::runtime_error(
                "Error when adding observations to SingleObservationSet: this operation is not supported when full "
                "weights are defined." );
    }

    if( weightState_.matrixType != diagonal_weights_matrix )
    {
        std::cerr << "Warning when adding observations to SingleObservationSet: resetting off-diagonal weights to diagonal defaults."
                  << std::endl;
        resetWeightsToUnitDiagonal( );
    }

    // Convert input diagonal weights into the internal append payload format.
    WeightState weightSubsetData;
    weightSubsetData.matrixType = diagonal_weights_matrix;
    weightSubsetData.diagonalWeights.reserve( observations.size( ) );

    for( unsigned int k = 0; k < observations.size( ); ++k )
    {
        if( !weights.empty( ) && weights.at( k ).size( ) != singleObservationSize_ )
        {
            throw std::runtime_error(
                    "Error when adding observations to SingleObservationSet, new weight "
                    "size is inconsistent." );
        }
        if( weights.empty( ) )
        {
            weightSubsetData.diagonalWeights.push_back( Eigen::Matrix< double, Eigen::Dynamic, 1 >::Ones( singleObservationSize_, 1 ) );
        }
        else
        {
            weightSubsetData.diagonalWeights.push_back( weights.at( k ) );
        }
    }

    addObservationsWithWeightData( observations, times, dependentVariables, residuals, weightSubsetData, sortObservations );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type
                  IsStateScalarAndTimeType >
Eigen::MatrixXd
//! Extracts the full-weight submatrix induced by a subset of observation indices.
SingleObservationSet< ObservationScalarType, TimeType, IsStateScalarAndTimeType >::extractFullWeightMatrixSubset(
        const std::vector< unsigned int >& indices ) const
{
    if( weightState_.fullWeights.rows( ) == 0 )
    {
        return Eigen::MatrixXd( 0, 0 );
    }

    const int blockSize = static_cast< int >( singleObservationSize_ );
    const int subsetSize = static_cast< int >( indices.size( ) ) * blockSize;
    Eigen::MatrixXd subsetMatrix = Eigen::MatrixXd::Zero( subsetSize, subsetSize );
    for( unsigned int i = 0; i < indices.size( ); i++ )
    {
        const int sourceRow = static_cast< int >( indices.at( i ) ) * blockSize;
        for( unsigned int j = 0; j < indices.size( ); j++ )
        {
            const int sourceCol = static_cast< int >( indices.at( j ) ) * blockSize;
            subsetMatrix.block( static_cast< int >( i ) * blockSize, static_cast< int >( j ) * blockSize, blockSize, blockSize ) =
                    weightState_.fullWeights.block( sourceRow, sourceCol, blockSize, blockSize );
        }
    }
    return subsetMatrix;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type
                  IsStateScalarAndTimeType >
typename SingleObservationSet< ObservationScalarType, TimeType, IsStateScalarAndTimeType >::WeightState
//! Extracts base/full weight payload for the selected observation indices.
SingleObservationSet< ObservationScalarType, TimeType, IsStateScalarAndTimeType >::extractWeightSubsetData(
        const std::vector< unsigned int >& indices ) const
{
    WeightState weightSubset;
    weightSubset.matrixType = weightState_.matrixType;

    if( weightState_.matrixType == diagonal_weights_matrix )
    {
        weightSubset.diagonalWeights.reserve( indices.size( ) );
        for( const unsigned int index : indices )
        {
            weightSubset.diagonalWeights.push_back( weightState_.diagonalWeights.at( index ) );
        }
    }
    else if( weightState_.matrixType == block_diagonal_weights_matrix )
    {
        weightSubset.blockWeights.reserve( indices.size( ) );
        for( const unsigned int index : indices )
        {
            weightSubset.blockWeights.push_back( weightState_.blockWeights.at( index ) );
        }
    }

    weightSubset.fullWeights = extractFullWeightMatrixSubset( indices );
    return weightSubset;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type
                  IsStateScalarAndTimeType >
//! Core append routine that inserts observations together with explicit weight payload data.
void SingleObservationSet< ObservationScalarType, TimeType, IsStateScalarAndTimeType >::addObservationsWithWeightData(
        const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
        const std::vector< TimeType >& times,
        const std::vector< Eigen::VectorXd >& dependentVariables,
        const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& residuals,
        const WeightState& weightSubsetData,
        const bool sortObservations )
{
    if( observations.empty( ) )
    {
        return;
    }

    // Validate observation/residual/dependent-variable dimensions before mutating state.
    validateObservationBatchSizes(
            observations, times, residuals, dependentVariables, "adding observations with explicit weight data to SingleObservationSet" );
    validateObservationDimensionsAgainstSingleSize( observations, "adding observations to SingleObservationSet" );
    validatePerObservationVectorSize( residuals, "residual", "adding observations to SingleObservationSet" );
    validateDependentVariableBatch( dependentVariables, observations.size( ), false, "adding observations to SingleObservationSet" );

    if( numberOfObservations_ > 0 && weightSubsetData.matrixType != weightState_.matrixType )
    {
        throw std::runtime_error(
                "Error when adding observations to SingleObservationSet, base weight matrix type is inconsistent with existing data." );
    }

    // Validate weight payload shape against the selected matrix representation.
    if( weightSubsetData.matrixType == diagonal_weights_matrix )
    {
        if( weightSubsetData.diagonalWeights.size( ) != observations.size( ) )
        {
            throw std::runtime_error(
                    "Error when adding observations to SingleObservationSet, number of diagonal weight vectors is inconsistent." );
        }
    }
    else if( weightSubsetData.matrixType == block_diagonal_weights_matrix )
    {
        if( weightSubsetData.blockWeights.size( ) != observations.size( ) )
        {
            throw std::runtime_error(
                    "Error when adding observations to SingleObservationSet, number of block weight matrices is inconsistent." );
        }
    }
    else
    {
        throw std::runtime_error( "Error when adding observations to SingleObservationSet, unsupported weight matrix type." );
    }

    const int oldObservationSetSize = static_cast< int >( getTotalObservationSetSize( ) );
    const int newObservationSetSize = static_cast< int >( observations.size( ) * singleObservationSize_ );
    if( weightSubsetData.fullWeights.rows( ) > 0 &&
        ( weightSubsetData.fullWeights.rows( ) != newObservationSetSize || weightSubsetData.fullWeights.cols( ) != newObservationSetSize ) )
    {
        throw std::runtime_error(
                "Error when adding observations to SingleObservationSet, full weight matrix contribution size is incompatible." );
    }

    if( numberOfObservations_ == 0 )
    {
        weightState_.matrixType = weightSubsetData.matrixType;
        weightState_.diagonalWeights.clear( );
        weightState_.blockWeights.clear( );
    }

    // Append base metadata and per-observation base weights.
    for( unsigned int k = 0; k < observations.size( ); k++ )
    {
        observations_.push_back( observations.at( k ) );
        observationTimes_.push_back( times.at( k ) );

        if( residuals.size( ) > 0 )
        {
            residuals_.push_back( residuals.at( k ) );
        }
        else
        {
            residuals_.push_back( Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( singleObservationSize_, 1 ) );
        }

        if( weightSubsetData.matrixType == diagonal_weights_matrix )
        {
            if( weightSubsetData.diagonalWeights.at( k ).size( ) != singleObservationSize_ )
            {
                throw std::runtime_error(
                        "Error when adding observations to SingleObservationSet, new diagonal weights size is inconsistent." );
            }
            weightState_.diagonalWeights.push_back( weightSubsetData.diagonalWeights.at( k ) );
        }
        else
        {
            if( weightSubsetData.blockWeights.at( k ).rows( ) != static_cast< int >( singleObservationSize_ ) ||
                weightSubsetData.blockWeights.at( k ).cols( ) != static_cast< int >( singleObservationSize_ ) )
            {
                throw std::runtime_error(
                        "Error when adding observations to SingleObservationSet, new block weights size is inconsistent." );
            }
            weightState_.blockWeights.push_back( weightSubsetData.blockWeights.at( k ) );
        }

        if( ( observationsDependentVariables_.size( ) > 0 || numberOfObservations_ == 0 ) && dependentVariables.size( ) > 0 )
        {
            observationsDependentVariables_.push_back( dependentVariables.at( k ) );
        }

        numberOfObservations_ += 1;
    }

    if( weightState_.fullWeights.rows( ) > 0 || weightSubsetData.fullWeights.rows( ) > 0 )
    {
        // Expand/merge full-weight contribution blocks for old and appended observations.
        Eigen::MatrixXd combinedFullWeights =
                Eigen::MatrixXd::Zero( oldObservationSetSize + newObservationSetSize, oldObservationSetSize + newObservationSetSize );
        if( weightState_.fullWeights.rows( ) > 0 )
        {
            combinedFullWeights.block( 0, 0, oldObservationSetSize, oldObservationSetSize ) = weightState_.fullWeights;
        }
        if( weightSubsetData.fullWeights.rows( ) > 0 )
        {
            combinedFullWeights.block( oldObservationSetSize, oldObservationSetSize, newObservationSetSize, newObservationSetSize ) =
                    weightSubsetData.fullWeights;
        }
        weightState_.fullWeights = combinedFullWeights;
    }

    if( sortObservations )
    {
        // Keep active observations in canonical time-sorted order if requested.
        orderObservationsAndMetadata( );
    }

    updateTimeBounds( );
    validateCombinedWeightsMatrix( );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type
                  IsStateScalarAndTimeType >
//! Moves observations between active and filtered sets while preserving full-weight block consistency.
void SingleObservationSet< ObservationScalarType, TimeType, IsStateScalarAndTimeType >::moveObservationsInOutFilteredSet(
        const std::vector< unsigned int >& indices,
        const bool moveInFilteredSet,
        const bool saveFilteredObservations )
{
    if( filteredObservationSet_ == nullptr )
    {
        throw std::runtime_error( "Error when moving observations in/out filtered set, filtered observation set pointer is empty." );
    }

    const unsigned int currentObservationCountBeforeMove = numberOfObservations_;
    const unsigned int filteredObservationCountBeforeMove = filteredObservationSet_->getNumberOfObservables( );

    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > observations;
    std::vector< TimeType > times;
    std::vector< Eigen::VectorXd > dependentVariables;
    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > residuals;

    // Move observations from active storage to filtered storage.
    if( moveInFilteredSet )
    {
        const Eigen::MatrixXd combinedFullWeightsBeforeMove =
                getCombinedFullWeightContributionMatrix( currentObservationCountBeforeMove, filteredObservationCountBeforeMove );
        const std::vector< unsigned int > retainedCurrentIndices =
                getComplementObservationIndices( currentObservationCountBeforeMove, indices );

        WeightState weightSubsetData;
        if( saveFilteredObservations )
        {
            weightSubsetData = extractWeightSubsetData( indices );
        }

        std::vector< TimeType > filteredTimesBeforeSorting = filteredObservationSet_->observationTimes_;
        for( auto index : indices )
        {
            if( index >= currentObservationCountBeforeMove )
            {
                throw std::runtime_error(
                        "Error when moving observation to filtered observation set, index "
                        "incompatible with number of observations." );
            }

            observations.push_back( observations_.at( index ) );
            times.push_back( observationTimes_.at( index ) );
            residuals.push_back( residuals_.at( index ) );
            filteredTimesBeforeSorting.push_back( observationTimes_.at( index ) );

            // If dependent variables not empty
            if( observationsDependentVariables_.size( ) > 0 )
            {
                dependentVariables.push_back( observationsDependentVariables_.at( index ) );
            }
        }

        if( saveFilteredObservations )
        {
            filteredObservationSet_->addObservationsWithWeightData(
                    observations, times, dependentVariables, residuals, weightSubsetData, true );
        }

        removeObservations( indices );

        const std::vector< unsigned int > currentIndicesAfterMove = retainedCurrentIndices;
        std::vector< unsigned int > filteredIndicesAfterMove;
        if( saveFilteredObservations )
        {
            const std::vector< unsigned int > oldFilteredCombinedIndices =
                    createSequentialObservationIndicesWithOffset( filteredObservationCountBeforeMove, currentObservationCountBeforeMove );
            const std::vector< unsigned int > movedCombinedIndices = indices;
            const std::vector< unsigned int > filteredCombinedConcatenatedIndices =
                    concatenateObservationIndices( oldFilteredCombinedIndices, movedCombinedIndices );
            const std::vector< std::size_t > filteredSortPermutation = getTimeSortingPermutation( filteredTimesBeforeSorting );
            filteredIndicesAfterMove = applyPermutationToObservationIndices( filteredCombinedConcatenatedIndices, filteredSortPermutation );
        }
        else
        {
            filteredIndicesAfterMove =
                    createSequentialObservationIndicesWithOffset( filteredObservationCountBeforeMove, currentObservationCountBeforeMove );
        }

        const std::vector< unsigned int > combinedIndicesAfterMove =
                concatenateObservationIndices( currentIndicesAfterMove, filteredIndicesAfterMove );
        const Eigen::MatrixXd combinedFullWeightsAfterMove =
                extractObservationBlockMatrix( combinedFullWeightsBeforeMove, combinedIndicesAfterMove, combinedIndicesAfterMove );
        setCombinedFullWeightContributionMatrix(
                combinedFullWeightsAfterMove,
                currentObservationCountBeforeMove - static_cast< unsigned int >( indices.size( ) ),
                filteredObservationCountBeforeMove + ( saveFilteredObservations ? static_cast< unsigned int >( indices.size( ) ) : 0 ) );
    }
    // Move observations from filtered storage back to active storage.
    else
    {
        const Eigen::MatrixXd combinedFullWeightsBeforeMove =
                getCombinedFullWeightContributionMatrix( currentObservationCountBeforeMove, filteredObservationCountBeforeMove );
        const std::vector< unsigned int > retainedFilteredIndices =
                getComplementObservationIndices( filteredObservationCountBeforeMove, indices );

        WeightState weightSubsetData = filteredObservationSet_->extractWeightSubsetData( indices );
        std::vector< TimeType > currentTimesBeforeSorting = observationTimes_;
        for( auto index : indices )
        {
            if( getNumberOfFilteredObservations( ) == 0 )
            {
                throw std::runtime_error(
                        "Error when moving observation back from filtered observation set, "
                        "filtered observation set is empty." );
            }
            if( index >= getNumberOfFilteredObservations( ) )
            {
                throw std::runtime_error(
                        "Error when moving observation back from filtered observation set, "
                        "index incompatible with number of observations." );
            }

            observations.push_back( filteredObservationSet_->getObservations( ).at( index ) );
            times.push_back( filteredObservationSet_->getObservationTimes( ).at( index ) );
            residuals.push_back( filteredObservationSet_->getResiduals( ).at( index ) );
            currentTimesBeforeSorting.push_back( filteredObservationSet_->getObservationTimes( ).at( index ) );

            // If dependent variables are set
            if( filteredObservationSet_->getObservationsDependentVariables( ).size( ) > 0 )
            {
                dependentVariables.push_back( filteredObservationSet_->getObservationsDependentVariables( ).at( index ) );
            }
        }

        addObservationsWithWeightData( observations, times, dependentVariables, residuals, weightSubsetData, true );
        filteredObservationSet_->removeObservations( indices );

        const std::vector< unsigned int > currentCombinedIndicesBeforeMove =
                createSequentialObservationIndices( currentObservationCountBeforeMove );
        const std::vector< unsigned int > movedCombinedIndices =
                applyOffsetToObservationIndices( indices, currentObservationCountBeforeMove );
        std::vector< unsigned int > currentCombinedConcatenatedIndices =
                concatenateObservationIndices( currentCombinedIndicesBeforeMove, movedCombinedIndices );
        const std::vector< std::size_t > currentSortPermutationBeforeMove = getTimeSortingPermutation( currentTimesBeforeSorting );
        const std::vector< unsigned int > currentIndicesAfterMove =
                applyPermutationToObservationIndices( currentCombinedConcatenatedIndices, currentSortPermutationBeforeMove );
        const std::vector< unsigned int > filteredIndicesAfterMove =
                applyOffsetToObservationIndices( retainedFilteredIndices, currentObservationCountBeforeMove );
        const std::vector< unsigned int > combinedIndicesAfterMove =
                concatenateObservationIndices( currentIndicesAfterMove, filteredIndicesAfterMove );
        const Eigen::MatrixXd combinedFullWeightsAfterMove =
                extractObservationBlockMatrix( combinedFullWeightsBeforeMove, combinedIndicesAfterMove, combinedIndicesAfterMove );
        setCombinedFullWeightContributionMatrix( combinedFullWeightsAfterMove,
                                                 currentObservationCountBeforeMove + static_cast< unsigned int >( indices.size( ) ),
                                                 filteredObservationCountBeforeMove - static_cast< unsigned int >( indices.size( ) ) );
    }
}

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_SINGLE_OBSERVATION_SET_FILTERING_AND_APPENDING_H
