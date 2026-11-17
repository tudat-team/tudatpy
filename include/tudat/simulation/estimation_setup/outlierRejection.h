/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OUTLIERREJECTION_H
#define TUDAT_OUTLIERREJECTION_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <Eigen/Core>

#include "tudat/simulation/estimation_setup/observationVectorData.h"
#include "tudat/simulation/estimation_setup/observationCondition.h"
#include "tudat/simulation/estimation_setup/observationDataset.h"
#include "tudat/simulation/estimation_setup/outlierRejectionSettings.h"
#include <tudat/astro/observation_models/observableTypes.h>

namespace tudat
{

namespace simulation_setup
{

using ObservableType = tudat::observation_models::ObservableType;

//! Data from a single estimation iteration that is provided to an outlier rejection algorithm.
/*!
 * Bundles information from the estimation that may be used in the outlier rejection algorithms, such as the observation data,
 * residuals, parameter covariance, design matrix, iteration number etc.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
struct OutlierRejectionInput {
    OutlierRejectionInput( const int iterationNumber,
                           const observation_models::ObservationVectorData< ObservationScalarType, TimeType >& flattenedObservationData,
                           const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& residuals,
                           const Eigen::MatrixXd& designMatrix,
                           const Eigen::MatrixXd& parameterCovariance,
                           const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& parameterCorrection ):
        iterationNumber_( iterationNumber ), flattenedObservationData_( flattenedObservationData ), residuals_( residuals ),
        designMatrix_( designMatrix ), parameterCovariance_( parameterCovariance ), parameterCorrection_( parameterCorrection )
    {}

    //! Return the rows of the design matrix H that belong to a single observation.
    /*!
     * The returned block has as many rows as the observable has components, and one column per estimated parameter.
     * This is what a criterion that is evaluated per observation needs.
     */
    Eigen::MatrixXd getDesignMatrixBlock( const unsigned int observationId ) const
    {
        const std::pair< int, unsigned int > rowRange = getRowRange( observationId );
        return designMatrix_.block( rowRange.first, 0, rowRange.second, designMatrix_.cols( ) );
    }

    //! Return the residuals of a single observation (one entry per component of the observable).
    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getResidualBlock( const unsigned int observationId ) const
    {
        const std::pair< int, unsigned int > rowRange = getRowRange( observationId );
        return residuals_.segment( rowRange.first, rowRange.second );
    }

    //! Return the covariance of a single observation (of size observable size by observable size).
    /*!
     * This is the diagonal block of the complete observation covariance, which for a right ascension/declination
     * pair, for instance, is the 2 x 2 covariance of that pair.
     */
    Eigen::MatrixXd getObservationCovariance( const unsigned int observationId ) const
    {
        // Preserve the input contract: only observations represented by this iteration may be queried, even though
        // ObservationVectorData retains the complete weight snapshot needed to recover the covariance block.
        static_cast< void >( getRowRange( observationId ) );
        return flattenedObservationData_.getInverseWeightMatrixForObservation( observationId );
    }

    //! Iteration of the estimation for which this data was computed, counted from zero.
    const int iterationNumber_;

    //! Flattened data over all observations, defining the row order of the residuals and the design matrix.
    const observation_models::ObservationVectorData< ObservationScalarType, TimeType >& flattenedObservationData_;

    //! Residual vector of all observations.
    const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& residuals_;

    //! Design matrix (partials of the observations w.r.t. the estimated parameters) of all observations.
    const Eigen::MatrixXd& designMatrix_;

    //! Covariance of the estimated parameters, as computed in the current iteration.
    const Eigen::MatrixXd& parameterCovariance_;

    //! Correction to the parameter vector computed in the current iteration.
    const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& parameterCorrection_;

private:
    //! Return the first row and the number of rows that one observation occupies in the residuals/design matrix.
    std::pair< int, unsigned int > getRowRange( const unsigned int observationId ) const
    {
        const int firstRow = flattenedObservationData_.getFirstVectorRowForObservation( observationId );
        const unsigned int numberOfRows = flattenedObservationData_.getScalarSizeForObservation( observationId );
        if( firstRow < 0 || numberOfRows == 0 )
        {
            throw std::runtime_error( "Error when retrieving outlier rejection input for observation " + std::to_string( observationId ) +
                                      ", observation is not present in the current iteration data." );
        }
        return std::make_pair( firstRow, numberOfRows );
    }
};

//! Return the outlier rejection settings, throwing if they are null.
template< typename SettingsType >
const std::shared_ptr< SettingsType >& getNonNullOutlierRejectionSettings( const std::shared_ptr< SettingsType >& outlierRejectionSettings )
{
    if( outlierRejectionSettings == nullptr )
    {
        throw std::runtime_error( "Error when creating outlier rejection object, settings are null." );
    }
    return outlierRejectionSettings;
}

//! Base class for algorithms that reject and recover outlying observations during an estimation.
/*!
 * One object of this type is created at the start of an estimation, from the settings that the user provided, and is
 * called once per iteration of the least-squares process. The object keeps the rejection status of every observation
 * and writes that status through to the ObservationDataset. The base class is not instantiated, the derived classes are.
 * Derived classes implement one function only: computeRejectionStatus, which fills the isRejected_ vector for the
 * current iteration. The surrounding steps (initializing the status from the dataset, and writing the new status back
 * to it) are identical for every algorithm and are handled by this base class.
 *
 * Each observation is assigned to (at most) one threshold group, defined by the observation selection conditions of the settings, which
 * are evaluated once, when this object is created. Observations that match no condition are not affected by the outlier rejection, and
 * an error is thrown if an observation matches more than one condition.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
class OutlierRejection
{
public:
    //! Constructor.
    /*!
     * \param outlierRejectionSettings Settings of the algorithm implemented by the derived class. They define the algorithm type (used
     * for logging and error messages), the first iteration with rejection, and the thresholds per observation selection condition; the
     * index of a condition in the settings is the index of its threshold group.
     * \param observationDataset Dataset holding the observations of the estimation. Observations that are already
     * rejected in this dataset (for instance by the user, before the estimation) start out as rejected.
     */
    OutlierRejection(
            const std::shared_ptr< OutlierRejectionSettings< ObservationScalarType, TimeType > >& outlierRejectionSettings,
            const std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > >& observationDataset ):
        outlierRejectionType_( getNonNullOutlierRejectionSettings( outlierRejectionSettings )->getOutlierRejectionType( ) ),
        firstIterationWithRejection_( outlierRejectionSettings->getFirstIterationWithRejection( ) ),
        numberOfThresholdGroups_( static_cast< unsigned int >( outlierRejectionSettings->getThresholdsPerCondition( ).size( ) ) ),
        observationDataset_( observationDataset ), outlierRejectionSettings_( outlierRejectionSettings )
    {
        if( observationDataset_ == nullptr )
        {
            throw std::runtime_error( "Error when creating outlier rejection object, observation dataset is null." );
        }

        for( const auto& row : observationDataset_->getObservationRows( ) )
        {
            observationIds_.push_back( row.observationId_ );
        }
        std::sort( observationIds_.begin( ), observationIds_.end( ) );
        isRejected_.resize( observationIds_.size( ) );
        for( unsigned int observationIndex = 0; observationIndex < observationIds_.size( ); observationIndex++ )
        {
            const unsigned int observationId = observationIds_.at( observationIndex );
            indexByObservationId_.emplace( observationId, observationIndex );
            isRejected_.at( observationIndex ) = !observationDataset_->getObservationRow( observationId ).isActive_;
        }
        isDisabled_ = isRejected_;
        numberOfDisabledObservations_ = getNumberOfRejectedObservations( );

        assignThresholdGroups( outlierRejectionSettings_->getConditions( ) );
    }

    //! Destructor
    virtual ~OutlierRejection( ) = default;

    //! Update the rejection status of all observations, using the data of the current estimation iteration.
    /*!
     * This function is called once per iteration by the estimation. It first lets the algorithm decide which
     * observations are outliers, and then applies that decision to the observation dataset. Before the first iteration with
     * rejection, the status is left unchanged.
     */
    void updateRejectionStatus( const OutlierRejectionInput< ObservationScalarType, TimeType >& outlierRejectionInput )
    {
        if( outlierRejectionInput.iterationNumber_ >= firstIterationWithRejection_ )
        {
            computeRejectionStatus( outlierRejectionInput );
        }

        if( isRejected_.size( ) != observationDataset_->getNumberOfObservations( ) )
        {
            throw std::runtime_error( "Error after computing outlier rejection status, size of rejection status vector (" +
                                      std::to_string( isRejected_.size( ) ) + ") is not equal to the number of observations (" +
                                      std::to_string( observationDataset_->getNumberOfObservations( ) ) + ")." );
        }

        applyRejectionStatusToObservationDataset( );
    }

    //! Return rejection status, in order of increasing observation id.
    const std::vector< bool >& getRejectionStatus( ) const
    {
        return isRejected_;
    }

    //! Reset the rejection status (in order of increasing observation id), and apply it to the observation dataset.
    void resetRejectionStatus( const std::vector< bool >& isRejected )
    {
        if( isRejected.size( ) != isRejected_.size( ) )
        {
            throw std::runtime_error( "Error when resetting outlier rejection status, size of rejection status vector (" +
                                      std::to_string( isRejected.size( ) ) + ") is not equal to the number of observations (" +
                                      std::to_string( isRejected_.size( ) ) + ")." );
        }
        isRejected_ = isRejected;
        applyRejectionStatusToObservationDataset( );
    }

    //! Return the number of observations that are currently rejected.
    unsigned int getNumberOfRejectedObservations( ) const
    {
        unsigned int numberOfRejectedObservations = 0;
        for( const bool observationIsRejected : isRejected_ )
        {
            if( observationIsRejected )
            {
                numberOfRejectedObservations++;
            }
        }
        return numberOfRejectedObservations;
    }

    int getFirstIterationWithRejection( ) const
    {
        return firstIterationWithRejection_;
    }

    //! Return the threshold group (index of the matching condition) of each observation, in order of increasing observation id; -1 if
    //! the observation matches no condition.
    const std::vector< int >& getThresholdGroupIndices( ) const
    {
        return thresholdGroupIndex_;
    }

    //! Return the number of threshold groups (equal to the number of conditions).
    unsigned int getNumberOfThresholdGroups( ) const
    {
        return numberOfThresholdGroups_;
    }

protected:
    //! A method that must update the isRejected_ member of OutlierRejection
    virtual void computeRejectionStatus( const OutlierRejectionInput< ObservationScalarType, TimeType >& outlierRejectionInput ) = 0;

    //! Write the current rejection status to the observation dataset.
    void applyRejectionStatusToObservationDataset( )
    {
        // Condition that selects exactly those observations that this object has marked as rejected. The condition is
        // evaluated by the dataset for each of its observations in turn.
        const observation_models::ObservationSelectionCondition< ObservationScalarType, TimeType > isRejectedCondition(
                [ this ]( const observation_models::ObservationDataset< ObservationScalarType, TimeType >&, const int observationId ) {
                    return isRejected_.at( indexByObservationId_.at( observationId ) );
                } );

        // Only reject active observations, so that already rejected ones (e.g. by the user) keep their rejection reason.
        observationDataset_->rejectObservations(
                isRejectedCondition && observation_models::ObservationSelectionCondition< ObservationScalarType, TimeType >::active( ),
                getOutlierRejectionTypeString( outlierRejectionType_ ) );
        observationDataset_->restoreObservations( !isRejectedCondition );
    }

    //! Return whether the algorithm may change the rejection status of an observation: it was not rejected at the start of the
    //! estimation, and it matches one of the conditions.
    bool isManaged( const unsigned int observationIndex ) const
    {
        return !isDisabled_.at( observationIndex ) && thresholdGroupIndex_.at( observationIndex ) >= 0;
    }

    //! Evaluate the conditions for each observation, and store the index of the matching condition in thresholdGroupIndex_.
    void assignThresholdGroups(
            const std::vector< observation_models::ObservationSelectionCondition< ObservationScalarType, TimeType > >& conditions )
    {
        thresholdGroupIndex_.assign( observationIds_.size( ), -1 );
        for( unsigned int observationIndex = 0; observationIndex < observationIds_.size( ); observationIndex++ )
        {
            const unsigned int observationId = observationIds_.at( observationIndex );
            for( unsigned int conditionIndex = 0; conditionIndex < conditions.size( ); conditionIndex++ )
            {
                if( !conditions.at( conditionIndex )( *observationDataset_, observationId ) )
                {
                    continue;
                }
                if( thresholdGroupIndex_.at( observationIndex ) >= 0 )
                {
                    throw std::runtime_error( "Error when creating " + getOutlierRejectionTypeString( outlierRejectionType_ ) +
                                              " object, observation " + std::to_string( observationId ) + " matches both condition " +
                                              std::to_string( thresholdGroupIndex_.at( observationIndex ) ) + " and condition " +
                                              std::to_string( conditionIndex ) +
                                              " of the outlier rejection settings; each observation may match at most one condition." );
                }
                thresholdGroupIndex_.at( observationIndex ) = static_cast< int >( conditionIndex );
            }
        }
    }

    //! Algorithm implemented by the derived class.
    const OutlierRejectionType outlierRejectionType_;

    //! Index of the first estimation iteration in which observations may be rejected.
    const int firstIterationWithRejection_;

    //! Observation ids in the order used by the rejection status vectors.
    std::vector< unsigned int > observationIds_;

    //! Map a persistent observation id to its position in the compact status vectors.
    std::unordered_map< unsigned int, unsigned int > indexByObservationId_;

    //! Rejection status of each observation in observationIds_ order; true means rejected.
    std::vector< bool > isRejected_;

    //! Vector that tracks if the observation was inactive at the start of estimation
    std::vector< bool > isDisabled_;

    unsigned int numberOfDisabledObservations_;

    //! Number of threshold groups (equal to the number of conditions).
    const unsigned int numberOfThresholdGroups_;

    //! Threshold group of each observation in observationIds_ order; -1 if the observation matches no condition.
    std::vector< int > thresholdGroupIndex_;

    //! Dataset holding the observations of the estimation.
    std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > > observationDataset_;

    //! Settings defining the thresholds of the algorithm.
    const std::shared_ptr< OutlierRejectionSettings< ObservationScalarType, TimeType > > outlierRejectionSettings_;
};

//! Outlier rejection algorithm of Carpino et al. (2003).
/*!
 * The algorithm computes, for each observation, a chi-squared value from the residual of that observation and the
 * covariance of that residual, and compares this value against the rejection and recovery thresholds in the
 * settings. Computing the residual covariance requires the observation covariance, the design matrix and the
 * covariance of the estimated parameters, all of which are available from the OutlierRejectionInput object. The adaptive rejection
 * threshold and the limit on the number of new rejections per iteration are computed separately for each threshold group.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
class CarpinoOutlierRejection : public OutlierRejection< ObservationScalarType, TimeType >
{
public:
    CarpinoOutlierRejection(
            const std::shared_ptr< CarpinoOutlierRejectionSettings< ObservationScalarType, TimeType > >& outlierRejectionSettings,
            const std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > >& observationDataset ):
        OutlierRejection< ObservationScalarType, TimeType >( outlierRejectionSettings, observationDataset )
    {
        // If the weights of the observation dataset are all equal to 1 (the default), they are most likely not set, in which case the
        // Carpino outlier rejection algorithm gives meaningless results. Unit weights are correct for an uncertainty of 1 (in
        // the unit of the observations), so this is a warning.
        ObservationVectorData< ObservationScalarType, TimeType > flattenedData = observationDataset->createObservationVectorData( true );
        const Eigen::VectorXd& observationWeightVector = flattenedData.getWeightVector( );
        const bool hasOffDiagonalWeights = flattenedData.hasOffDiagonalWeights( );
        const bool hasDefaultWeights = ( observationWeightVector.array( ) == 1.0 ).all( );
        if( !hasOffDiagonalWeights && hasDefaultWeights )
        {
            std::cerr << "Warning in outlier rejection: the Carpino algorithm requires observation weights equal to the inverse of "
                         "the squared observation uncertainty, but all weights of the ObservationDataset are equal to 1 (the "
                         "default). Set the weights, unless an uncertainty of 1 is intended."
                      << std::endl;
        }
    }

protected:
    //! Determine which observations are outliers according to the algorithm of Carpino et al. (2003).
    //! Sets the rejection flags inside isRejected_ based on the chi-squared values of the residuals in the current iteration

    void computeRejectionStatus( const OutlierRejectionInput< ObservationScalarType, TimeType >& outlierRejectionInput ) override
    {
        // Compute chi-squared for each observation, using its persistent dataset id to find the matching rows. Observations that are
        // not managed by the algorithm (rejected at the start of the estimation, or matching no condition) get no chi-squared value.
        std::vector< double > chiSquaredPerObservation( this->isRejected_.size( ), 0.0 );
        for( unsigned int observationIndex = 0; observationIndex < this->isRejected_.size( ); observationIndex++ )
        {
            if( !this->isManaged( observationIndex ) )
            {
                continue;
            }
            const unsigned int observationId = this->observationIds_.at( observationIndex );
            chiSquaredPerObservation.at( observationIndex ) =
                    computeChiSquared( outlierRejectionInput.getDesignMatrixBlock( observationId ),
                                       outlierRejectionInput.getResidualBlock( observationId ).template cast< double >( ),
                                       outlierRejectionInput.parameterCorrection_.template cast< double >( ),
                                       outlierRejectionInput.parameterCovariance_,
                                       outlierRejectionInput.getObservationCovariance( observationId ),
                                       this->getRejectionStatus( ).at( observationIndex ) );
        }

        const std::vector< double > rejectionThresholdPerGroup = getRejectionThresholdPerGroup( chiSquaredPerObservation );

        // Rejection status for the next iteration. The isRejected_ member keeps the status of the current
        // iteration until the end of this function, so that every step below can compare against it.
        std::vector< bool > newRejectionStatus = this->isRejected_;
        for( unsigned int observationIndex = 0; observationIndex < newRejectionStatus.size( ); observationIndex++ )
        {
            if( !this->isManaged( observationIndex ) )
            {
                continue;
            }
            const int groupIndex = this->thresholdGroupIndex_.at( observationIndex );
            const bool isCurrentlyRejected = this->getRejectionStatus( ).at( observationIndex );
            newRejectionStatus.at( observationIndex ) =
                    decideRejectionStatus( isCurrentlyRejected,
                                           chiSquaredPerObservation.at( observationIndex ),
                                           this->outlierRejectionSettings_->getThresholds( groupIndex ).first,
                                           rejectionThresholdPerGroup.at( groupIndex ) );
        }

        // Limit the number of observations that are newly rejected in this iteration
        applyMaximumRejectionsPerIteration( newRejectionStatus, chiSquaredPerObservation );

        this->isRejected_ = newRejectionStatus;
    }

    //! No more than 20% of the number of available observations (within a group) may be rejected per iteration, although at least
    //! one observation per group may always be rejected
    static constexpr double maximumRejectedFractionPerIteration = 0.2;

    //! Limit the number of observations that are newly rejected in this iteration, per threshold group.
    /*!
     * At most maximumRejectedFractionPerIteration of the observations of a threshold group that are in the fit at the start of this
     * iteration may be newly rejected, rounded down, but with a minimum of one, so that an outlier can also be rejected from a group
     * with fewer than 1 / maximumRejectedFractionPerIteration observations. If more observations of that group meet the rejection
     * criterion, only those with the largest
     * chi-squared values are rejected, and the others stay in the fit. Observations that were already rejected, observations that are
     * recovered, and observations that are not managed by the algorithm, are not touched and do not count towards the limit.
     *
     * \param newRejectionStatus Rejection status for the next iteration, modified in place by this function.
     * \param chiSquaredPerObservation Chi-squared value of each observation, in order of increasing observation id.
     */
    void applyMaximumRejectionsPerIteration( std::vector< bool >& newRejectionStatus,
                                             const std::vector< double >& chiSquaredPerObservation ) const
    {
        // Count the observations of each group that are in the fit at the start of this iteration, and collect those that are newly
        // rejected
        std::vector< std::size_t > numberOfAcceptedObservations( this->numberOfThresholdGroups_, 0 );
        std::vector< std::vector< unsigned int > > newlyRejectedObservationIndices( this->numberOfThresholdGroups_ );
        for( unsigned int observationIndex = 0; observationIndex < newRejectionStatus.size( ); observationIndex++ )
        {
            if( !this->isManaged( observationIndex ) || this->isRejected_.at( observationIndex ) )
            {
                continue;
            }
            const int groupIndex = this->thresholdGroupIndex_.at( observationIndex );
            numberOfAcceptedObservations.at( groupIndex )++;
            if( newRejectionStatus.at( observationIndex ) )
            {
                newlyRejectedObservationIndices.at( groupIndex ).push_back( observationIndex );
            }
        }

        std::size_t numberOfKeptObservations = 0;
        for( unsigned int groupIndex = 0; groupIndex < this->numberOfThresholdGroups_; groupIndex++ )
        {
            std::vector< unsigned int >& newlyRejectedInGroup = newlyRejectedObservationIndices.at( groupIndex );
            const std::size_t allowedNumberOfNewRejections = std::max(
                    static_cast< std::size_t >( 1 ),
                    static_cast< std::size_t >( std::floor( maximumRejectedFractionPerIteration *
                                                            static_cast< double >( numberOfAcceptedObservations.at( groupIndex ) ) ) ) );

            if( newlyRejectedInGroup.size( ) <= allowedNumberOfNewRejections )
            {
                continue;
            }

            // Order the new rejections by decreasing chi-squared, so that the worst observations are the ones that are rejected
            std::stable_sort( newlyRejectedInGroup.begin( ),
                              newlyRejectedInGroup.end( ),
                              [ &chiSquaredPerObservation ]( const unsigned int firstIndex, const unsigned int secondIndex ) {
                                  return chiSquaredPerObservation.at( firstIndex ) > chiSquaredPerObservation.at( secondIndex );
                              } );

            for( std::size_t i = allowedNumberOfNewRejections; i < newlyRejectedInGroup.size( ); i++ )
            {
                newRejectionStatus.at( newlyRejectedInGroup.at( i ) ) = false;
            }
            numberOfKeptObservations += newlyRejectedInGroup.size( ) - allowedNumberOfNewRejections;
        }

        if( numberOfKeptObservations > 0 )
        {
            std::cerr << "Warning during outlier rejection, " << numberOfKeptObservations
                      << " observations that met the rejection criterion were kept in the fit, since at most "
                      << 100.0 * maximumRejectedFractionPerIteration
                      << "% of the observations in the fit may be rejected per iteration, per outlier rejection condition." << std::endl;
        }
    }

    //! Compute the rejection threshold of each threshold group for the current iteration.
    /*!
     * Carpino et al. (2003), Section 2.4: the threshold of a group is the maximum of a quarter of the largest chi-squared value of the
     * accepted observations in that group, and the rejection threshold from the settings plus a fudge term that increases when few
     * observations are left in the fit. The fudge term reflects how well the solution is determined, so it uses the number of all
     * observations in the fit, including those of other groups and those that match no condition.
     */
    std::vector< double > getRejectionThresholdPerGroup( const std::vector< double >& chiSquaredPerObservation ) const
    {
        // Maximum chi-squared of the accepted observations per group, and number of observations in the fit (in all groups, and
        // including observations that match no condition)
        std::vector< double > chiSquaredMax( this->numberOfThresholdGroups_, 0.0 );
        int numberOfObservationsInFit = 0;
        for( unsigned int observationIndex = 0; observationIndex < this->isRejected_.size( ); observationIndex++ )
        {
            if( this->isRejected_.at( observationIndex ) )
            {
                continue;
            }
            numberOfObservationsInFit++;
            if( this->isManaged( observationIndex ) )
            {
                const int groupIndex = this->thresholdGroupIndex_.at( observationIndex );
                chiSquaredMax.at( groupIndex ) =
                        std::max( chiSquaredMax.at( groupIndex ), chiSquaredPerObservation.at( observationIndex ) );
            }
        }

        // Increases rejection threshold when small number of observations are left in the fit
        const double fudgeTerm = 400 * std::pow( 1.2, -numberOfObservationsInFit );

        std::vector< double > rejectionThresholdPerGroup( this->numberOfThresholdGroups_ );
        for( unsigned int groupIndex = 0; groupIndex < this->numberOfThresholdGroups_; groupIndex++ )
        {
            rejectionThresholdPerGroup.at( groupIndex ) =
                    std::max( 0.25 * chiSquaredMax.at( groupIndex ),
                              this->outlierRejectionSettings_->getThresholds( groupIndex ).second + fudgeTerm );
        }
        return rejectionThresholdPerGroup;
    }

    bool decideRejectionStatus( const bool isCurrentlyRejected,
                                const double chiSquared,
                                const double chiSquaredRecoveryThreshold,
                                const double chiSquaredRejectionThreshold ) const
    {
        if( chiSquared >= 0 && std::isfinite( chiSquared ) )
        {
            if( isCurrentlyRejected )
            {
                // Observation was rejected but is now below recovery threshold
                if( chiSquared < chiSquaredRecoveryThreshold )
                {
                    return false;
                }
            }
            else
            {
                // Observation was accepted but is now above rejection threshold
                if( chiSquared > chiSquaredRejectionThreshold )
                {
                    return true;
                }
            }
        }
        else
        {
            std::cerr << "Warning: in Carpino outlier rejection, one or more chi-squared values are negative/NaN. Leaving rejection flag "
                      << "unchanged" << std::endl;
        }
        return isCurrentlyRejected;  // Return unchanged
    }

    //! Compute Chi2 for one observation. See Carpino (2003) for details
    double computeChiSquared( const Eigen::MatrixXd& partialsForObservation,
                              const Eigen::VectorXd& residualVector,
                              const Eigen::VectorXd& parameterCorrection,
                              const Eigen::MatrixXd& parameterCovariance,
                              const Eigen::MatrixXd& singleObservationCovariance,
                              const bool isRejected ) const
    {
        // Get linear approximation of post-fit residuals
        const Eigen::VectorXd residualsVectorPostFit = residualVector - partialsForObservation * parameterCorrection;

        // Residual covariance. Note this is different from pure observation covariance
        Eigen::MatrixXd residualCovariance;
        if( isRejected )
        {
            residualCovariance =
                    singleObservationCovariance + partialsForObservation * parameterCovariance * partialsForObservation.transpose( );
        }
        else
        {
            residualCovariance =
                    singleObservationCovariance - partialsForObservation * parameterCovariance * partialsForObservation.transpose( );
        }

        // Chi-squared value
        return ( residualsVectorPostFit.transpose( ) * residualCovariance.inverse( ) * residualsVectorPostFit ).value( );
    }
};

//! Simple outlier rejection algorithm
/*!
 * Rejects an accepted observation if its largest absolute residual exceeds the rejection threshold of its threshold group, and recovers a
 * rejected observation if it drops below the recovery threshold of its threshold group.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
class SimpleOutlierRejection : public OutlierRejection< ObservationScalarType, TimeType >
{
public:
    SimpleOutlierRejection(
            const std::shared_ptr< SimpleOutlierRejectionSettings< ObservationScalarType, TimeType > >& outlierRejectionSettings,
            const std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > >& observationDataset ):
        OutlierRejection< ObservationScalarType, TimeType >( outlierRejectionSettings, observationDataset )
    {}

protected:
    void computeRejectionStatus( const OutlierRejectionInput< ObservationScalarType, TimeType >& outlierRejectionInput ) override
    {
        for( unsigned int observationIndex = 0; observationIndex < this->isRejected_.size( ); observationIndex++ )
        {
            // Skip observations that were excluded at the start of the estimation, or that match no condition
            if( !this->isManaged( observationIndex ) )
            {
                continue;
            }
            const unsigned int observationId = this->observationIds_.at( observationIndex );

            // (recovery, rejection) thresholds for this observation
            const std::pair< double, double >& residualThresholds =
                    this->outlierRejectionSettings_->getThresholds( this->thresholdGroupIndex_.at( observationIndex ) );

            const double maximumAbsoluteResidual =
                    outlierRejectionInput.getResidualBlock( observationId ).template cast< double >( ).cwiseAbs( ).maxCoeff( );
            if( this->isRejected_.at( observationIndex ) )
            {
                if( maximumAbsoluteResidual < residualThresholds.first )
                {
                    this->isRejected_.at( observationIndex ) = false;
                }
            }
            else if( maximumAbsoluteResidual > residualThresholds.second )
            {
                this->isRejected_.at( observationIndex ) = true;
            }
        }
    }
};

//! Function to create the outlier rejection object defined by a set of outlier rejection settings.
/*!
 * \param outlierRejectionSettings Settings provided by the user; a null pointer means that no outlier rejection is to
 * be performed, in which case a null pointer is returned.
 * \param observationDataset Dataset holding the observations of the estimation.
 * \return Outlier rejection object to be used during the estimation, or a null pointer if no settings were provided.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< OutlierRejection< ObservationScalarType, TimeType > > createOutlierRejection(
        const std::shared_ptr< OutlierRejectionSettings< ObservationScalarType, TimeType > >& outlierRejectionSettings,
        const std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > >& observationDataset )
{
    if( outlierRejectionSettings == nullptr )
    {
        return nullptr;
    }

    switch( outlierRejectionSettings->getOutlierRejectionType( ) )
    {
        case OutlierRejectionType::carpino_outlier_rejection: {
            // The settings are held as a pointer to the base class, so they must be converted to a pointer to the
            // derived class before the algorithm-specific parameters can be read. dynamic_pointer_cast performs this
            // conversion, and returns a null pointer if the object is not of the requested type.
            const std::shared_ptr< CarpinoOutlierRejectionSettings< ObservationScalarType, TimeType > > carpinoSettings =
                    std::dynamic_pointer_cast< CarpinoOutlierRejectionSettings< ObservationScalarType, TimeType > >(
                            outlierRejectionSettings );
            if( carpinoSettings == nullptr )
            {
                throw std::runtime_error(
                        "Error when creating outlier rejection object, settings are of type carpino_outlier_rejection, but are not "
                        "CarpinoOutlierRejectionSettings." );
            }
            return std::make_shared< CarpinoOutlierRejection< ObservationScalarType, TimeType > >( carpinoSettings, observationDataset );
        }

        case OutlierRejectionType::simple_outlier_rejection: {
            const std::shared_ptr< SimpleOutlierRejectionSettings< ObservationScalarType, TimeType > > simpleSettings =
                    std::dynamic_pointer_cast< SimpleOutlierRejectionSettings< ObservationScalarType, TimeType > >(
                            outlierRejectionSettings );
            if( simpleSettings == nullptr )
            {
                throw std::runtime_error(
                        "Error when creating outlier rejection object, settings are of type simple_outlier_rejection, but are not "
                        "SimpleOutlierRejectionSettings." );
            }
            return std::make_shared< SimpleOutlierRejection< ObservationScalarType, TimeType > >( simpleSettings, observationDataset );
        }
        default:
            throw std::runtime_error( "Error when creating outlier rejection object, type " +
                                      getOutlierRejectionTypeString( outlierRejectionSettings->getOutlierRejectionType( ) ) +
                                      " not recognized." );
    }
}

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_OUTLIERREJECTION_H
