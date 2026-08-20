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
#include <vector>

#include <Eigen/Core>

#include "tudat/simulation/estimation_setup/flattenedObservationData.h"
#include "tudat/simulation/estimation_setup/observationCondition.h"
#include "tudat/simulation/estimation_setup/observationDataset.h"
#include "tudat/simulation/estimation_setup/outlierRejectionSettings.h"

namespace tudat
{

namespace simulation_setup
{

//! Data from a single estimation iteration that is provided to an outlier rejection algorithm.
/*!
 * This object bundles the quantities that the least-squares iteration has computed anyway. Which of these an
 * algorithm uses is up to the algorithm: an algorithm that only inspects residual values simply ignores the
 * remaining entries. The object stores references only, so creating it is cheap.
 *
 * All quantities are unnormalized (that is, they are the physical quantities, not the internally scaled ones that
 * the least-squares solution uses), and they cover *all* observations in the dataset, including those that are
 * currently rejected and that therefore did not contribute to the estimation. The latter is required for algorithms
 * to be able to recover observations that were rejected in an earlier iteration.
 *
 * The observation covariance is the inverse of the complete weight matrix. The complete matrix is inverted, rather
 * than the weight block of each observation, because the weight matrix is not assumed to be block diagonal:
 * correlations between arbitrary observations may be defined, and in that case the covariance of a single
 * observation is a block of the inverse of the complete matrix, which is not the same as the inverse of that
 * observation's own weight block.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
struct OutlierRejectionInput {
    OutlierRejectionInput(
            const int iterationNumber,
            const observation_models::FlattenedObservationData< ObservationScalarType, TimeType >& flattenedObservationData,
            const Eigen::MatrixXd& observationCovariance,
            const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& residuals,
            const Eigen::MatrixXd& designMatrix,
            const Eigen::MatrixXd& parameterCovariance,
            const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& parameterCorrection ):
        iterationNumber_( iterationNumber ), flattenedObservationData_( flattenedObservationData ),
        observationCovariance_( observationCovariance ), residuals_( residuals ), designMatrix_( designMatrix ),
        parameterCovariance_( parameterCovariance ), parameterCorrection_( parameterCorrection )
    { }

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
        const std::pair< int, unsigned int > rowRange = getRowRange( observationId );
        return observationCovariance_.block( rowRange.first, rowRange.first, rowRange.second, rowRange.second );
    }

    //! Iteration of the estimation for which this data was computed, counted from zero.
    const int iterationNumber_;

    //! Flattened data over all observations, defining the row order of the residuals and the design matrix.
    const observation_models::FlattenedObservationData< ObservationScalarType, TimeType >& flattenedObservationData_;

    //! Covariance of all observations: the inverse of the complete weight matrix, in the row order above.
    const Eigen::MatrixXd& observationCovariance_;

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
        const int firstRow = flattenedObservationData_.getFirstFlattenedRowForObservation( observationId );
        const unsigned int numberOfRows = flattenedObservationData_.getScalarSizeForObservation( observationId );
        if( firstRow < 0 || numberOfRows == 0 )
        {
            throw std::runtime_error( "Error when retrieving outlier rejection input for observation " +
                                      std::to_string( observationId ) + ", observation is not present in the current iteration data." );
        }
        return std::make_pair( firstRow, numberOfRows );
    }
};

//! Base class for algorithms that reject and recover outlying observations during an estimation.
/*!
 * One object of this type is created at the start of an estimation, from the settings that the user provided, and is
 * called once per iteration of the least-squares process. The object keeps the rejection status of every observation
 * and writes that status through to the ObservationDataset, which is what the estimation reads to decide which
 * observations to use. The dataset therefore remains the single place where the rejection status is stored.
 *
 * Derived classes implement one function only: computeRejectionStatus, which fills the isRejected_ vector for the
 * current iteration. The surrounding steps (initializing the status from the dataset, and writing the new status back
 * to it) are identical for every algorithm and are handled by this base class.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
class OutlierRejection
{
public:
    //! Constructor.
    /*!
     * \param outlierRejectionType Algorithm implemented by the derived class, used for logging and error messages.
     * \param observationDataset Dataset holding the observations of the estimation. Observations that are already
     * rejected in this dataset (for instance by the user, before the estimation) start out as rejected.
     */
    OutlierRejection( const OutlierRejectionType outlierRejectionType,
                      const std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > >& observationDataset ):
        outlierRejectionType_( outlierRejectionType ), observationDataset_( observationDataset )
    {
        if( observationDataset_ == nullptr )
        {
            throw std::runtime_error( "Error when creating outlier rejection object, observation dataset is null." );
        }

        isRejected_.resize( observationDataset_->getNumberOfObservations( ) );
        for( unsigned int observationId = 0; observationId < isRejected_.size( ); observationId++ )
        {
            isRejected_.at( observationId ) = !observationDataset_->getObservationRow( observationId ).isActive_;
        }
    }

    //! Destructor (virtual, since objects of derived types are deleted through a pointer to this base class).
    virtual ~OutlierRejection( ) = default;

    //! Update the rejection status of all observations, using the data of the current estimation iteration.
    /*!
     * This function is called once per iteration by the estimation. It first lets the algorithm decide which
     * observations are outliers, and then applies that decision to the observation dataset.
     */
    void updateRejectionStatus( const OutlierRejectionInput< ObservationScalarType, TimeType >& outlierRejectionInput )
    {
        computeRejectionStatus( outlierRejectionInput );

        if( isRejected_.size( ) != observationDataset_->getNumberOfObservations( ) )
        {
            throw std::runtime_error( "Error after computing outlier rejection status, size of rejection status vector (" +
                                      std::to_string( isRejected_.size( ) ) + ") is not equal to the number of observations (" +
                                      std::to_string( observationDataset_->getNumberOfObservations( ) ) + ")." );
        }

        applyRejectionStatusToObservationDataset( );
    }

    //! Return the rejection status of each observation, indexed by observation id.
    const std::vector< bool >& getRejectionStatus( ) const
    {
        return isRejected_;
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

    OutlierRejectionType getOutlierRejectionType( ) const
    {
        return outlierRejectionType_;
    }

protected:
    //! Determine which observations are outliers, by filling the isRejected_ vector.
    /*!
     * The '= 0' makes this function pure virtual: this class provides no implementation, every derived class must
     * provide one, and an object of this base class can no longer be created directly.
     *
     * Implementations must set an entry for every observation in the dataset, both for observations that are
     * currently used in the estimation (which may become rejected) and for observations that are currently rejected
     * (which may be recovered).
     */
    virtual void computeRejectionStatus( const OutlierRejectionInput< ObservationScalarType, TimeType >& outlierRejectionInput ) = 0;

    //! Write the current rejection status to the observation dataset.
    void applyRejectionStatusToObservationDataset( )
    {
        // Condition that selects exactly those observations that this object has marked as rejected. The condition is
        // evaluated by the dataset for each of its observations in turn.
        const observation_models::ObservationSelectionCondition< ObservationScalarType, TimeType > isRejectedCondition(
                [ this ]( const observation_models::ObservationDataset< ObservationScalarType, TimeType >&, const int observationId ) {
                    return isRejected_.at( observationId );
                } );

        observationDataset_->rejectObservations( isRejectedCondition, getOutlierRejectionTypeString( outlierRejectionType_ ) );
        observationDataset_->restoreObservations( !isRejectedCondition );
    }

    //! Algorithm implemented by the derived class.
    const OutlierRejectionType outlierRejectionType_;

    //! Rejection status of each observation, indexed by observation id; true means the observation is rejected.
    std::vector< bool > isRejected_;

    //! Dataset holding the observations of the estimation.
    std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > > observationDataset_;
};

//! Outlier rejection algorithm of Carpino et al. (2003).
/*!
 * The algorithm computes, for each observation, a chi-squared value from the residual of that observation and the
 * covariance of that residual, and compares this value against the rejection and recovery thresholds in the
 * settings. Computing the residual covariance requires the observation covariance, the design matrix and the
 * covariance of the estimated parameters, all of which are available from the OutlierRejectionInput object.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
class CarpinoOutlierRejection : public OutlierRejection< ObservationScalarType, TimeType >
{
public:
    CarpinoOutlierRejection(
            const std::shared_ptr< CarpinoOutlierRejectionSettings >& outlierRejectionSettings,
            const std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > >& observationDataset ):
        OutlierRejection< ObservationScalarType, TimeType >( OutlierRejectionType::carpino_outlier_rejection, observationDataset ),
        outlierRejectionSettings_( outlierRejectionSettings )
    {
        if( outlierRejectionSettings_ == nullptr )
        {
            throw std::runtime_error( "Error when creating Carpino outlier rejection object, settings are null." );
        }

        //! If the weights of the observation dataset are all equal to 1 (i.e. default), they are not meaningful and the Carpino outlier
        //! rejection algorithm produces garbage results
        FlattenedObservationData<ObservationScalarType, TimeType> flattenedData =
            observationDataset->createComputationFlattenedObservationData(  );
        const Eigen::VectorXd& observationWeightVector = flattenedData.getWeightVector(  );
        const bool hasOffDiagonalWeights = flattenedData.hasOffDiagonalWeights(  );
        const bool hasDefaultWeights = (observationWeightVector.array() == 1.0).all();
        if(!hasOffDiagonalWeights && hasDefaultWeights)
        {
            throw std::runtime_error("Error in outlier rejection: Carpino algorithm expects observation weights to be compatible with the"
                                     "per-observation uncertainty, but default weights were found on the ObservationDataset");
        }

    }

protected:
    //! Determine which observations are outliers according to the algorithm of Carpino et al. (2003).
    //! Sets the rejection flags inside isRejected_ based on the chi-squared values of the residuals in the current iteration

    void computeRejectionStatus( const OutlierRejectionInput< ObservationScalarType, TimeType >& outlierRejectionInput ) override
    {
        // Only update rejections if later than first iteration with rejection
        if(outlierRejectionInput.iterationNumber_ >= outlierRejectionSettings_->getFirstIterationWithRejection(  ))
        {
            std::vector<double> chiSquaredPerObservation;

            // Compute chi-squared for each observation ID
            for(unsigned int observationId = 0; observationId < this->isRejected_.size(); observationId++)
            {
                const double chiSquared = computeChiSquared(
                    outlierRejectionInput.getDesignMatrixBlock( observationId ),
                    outlierRejectionInput.getResidualBlock( observationId ).template cast<double>(),
                    outlierRejectionInput.parameterCorrection_.template cast<double>() ,
                    outlierRejectionInput.parameterCovariance_,
                    outlierRejectionInput.getObservationCovariance( observationId ),
                    this->getRejectionStatus(  ).at( observationId ));

                chiSquaredPerObservation.push_back( chiSquared );
            }

            const double rejectionThreshold = getRejectionThreshold( chiSquaredPerObservation );

            // Rejection status for the next iteration. The isRejected_ member keeps the status of the current
            // iteration until the end of this function, so that every step below can compare against it.
            std::vector< bool > newRejectionStatus = this->isRejected_;
            for(unsigned int observationId = 0; observationId < newRejectionStatus.size(); observationId++)
            {
                const bool isCurrentlyRejected = this->getRejectionStatus(  ).at( observationId );
                newRejectionStatus.at(observationId) = decideRejectionStatus( isCurrentlyRejected,
                    chiSquaredPerObservation.at( observationId ), rejectionThreshold);
            }

            // Prevent too many observations being rejected at once during one iteration
            applyMaximumRejectedFraction( newRejectionStatus, chiSquaredPerObservation );

            this->isRejected_ = newRejectionStatus;
        }
    }

    //! Limit the number of observations rejected in this iteration to the fraction that the settings allow.
    /*!
     * Observations that were already rejected before this iteration take up part of the allowed number of rejected
     * observations, but are not reconsidered here: their chi-squared is computed with the observation outside the
     * fit, and is therefore not comparable with the chi-squared of the observations that are rejected in this
     * iteration. Of the newly rejected observations, only those with the highest chi-squared remain rejected, and
     * the remainder is accepted again. Observations with an equal chi-squared are ordered by observation id, so that
     * the result does not depend on the order in which the observations happen to be stored.
     *
     * \param newRejectionStatus Rejection status for the next iteration, modified in place by this function.
     * \param chiSquaredPerObservation Chi-squared value of each observation, indexed by observation id.
     */
    void applyMaximumRejectedFraction( std::vector< bool >& newRejectionStatus,
                                       const std::vector< double >& chiSquaredPerObservation ) const
    {
        const std::size_t maximumNumberRejected = static_cast< std::size_t >( std::floor(
                outlierRejectionSettings_->getMaximumRejectedFraction( ) * static_cast< double >( newRejectionStatus.size( ) ) ) );

        // Split the observations that are to be rejected into those that were already rejected before this
        // iteration, and those that are rejected by the criteria of this iteration
        std::size_t numberOfRetainedRejections = 0;
        std::vector< unsigned int > newlyRejectedObservationIds;
        for( unsigned int observationId = 0; observationId < newRejectionStatus.size( ); observationId++ )
        {
            if( newRejectionStatus.at( observationId ) )
            {
                if( this->isRejected_.at( observationId ) )
                {
                    numberOfRetainedRejections++; // Number of observations that were already rejected
                }
                else
                {
                    newlyRejectedObservationIds.push_back( observationId ); // Observations that are now, but were not rejected
                }
            }
        }

        if( numberOfRetainedRejections + newlyRejectedObservationIds.size( ) <= maximumNumberRejected )
        {
            return;
        }

        // Number of new rejections that fits in what is left of the allowed number of rejected observations
        const std::size_t allowedNumberOfNewRejections =
                ( maximumNumberRejected > numberOfRetainedRejections ) ? maximumNumberRejected - numberOfRetainedRejections : 0;

        // Order the new rejections by decreasing chi-squared, so that the worst observations are the ones that are kept
        std::stable_sort( newlyRejectedObservationIds.begin( ),
                          newlyRejectedObservationIds.end( ),
                          [ &chiSquaredPerObservation ]( const unsigned int firstObservationId, const unsigned int secondObservationId ) {
                              return chiSquaredPerObservation.at( firstObservationId ) >
                                      chiSquaredPerObservation.at( secondObservationId );
                          } );

        for( std::size_t i = allowedNumberOfNewRejections; i < newlyRejectedObservationIds.size( ); i++ )
        {
            newRejectionStatus.at( newlyRejectedObservationIds.at( i ) ) = false;
        }

        std::cerr << "Warning during outlier rejection, " << newlyRejectedObservationIds.size( ) - allowedNumberOfNewRejections
                  << " observations that met the rejection criterion were kept, since the maximum rejected fraction of "
                  << outlierRejectionSettings_->getMaximumRejectedFraction( ) << " was reached." << std::endl;
    }

    double getRejectionThreshold(const std::vector<double>& chiSquaredPerObservation) const
    {
        // Get Maximum chi-squared of all accepted observations
        double chiSquaredMax = 0.0;
        for(unsigned int observationId = 0; observationId < this->isRejected_.size(); observationId++)
        {
            if(!this->isRejected_.at( observationId ))
            {
                chiSquaredMax = std::max(chiSquaredMax, chiSquaredPerObservation.at(observationId));
            }
        }

        const int numberOfAcceptedObservations = this->observationDataset_->getNumberOfObservations(  ) -
            this->getNumberOfRejectedObservations(  );

        // Increases rejection threshold when small number of observations are left
        const double fudgeTerm = 400 * std::pow(1.2, -numberOfAcceptedObservations);

        double rejectionThreshold =
            std::max(0.25 * chiSquaredMax, this->outlierRejectionSettings_->getChi2RejectionThreshold(  )) + fudgeTerm;

        return rejectionThreshold;

    }

    bool decideRejectionStatus( const bool isCurrentlyRejected, const double chiSquared, const double chiSquaredRejectionThreshold ) const
    {
        if(chiSquared > 0 && std::isfinite(chiSquared))
        {
            if( isCurrentlyRejected )
            {
                // Observation was rejected but is now below recovery threshold
                if(chiSquared < outlierRejectionSettings_->getChi2RecoveryThreshold(  ))
                {
                    return false;
                }
            }
            else
            {
                // Observation was accepted but is now above rejection threshold
                if(chiSquared > chiSquaredRejectionThreshold)
                {
                    return true;
                }
            }
        }
        else
        {
            std::cerr << "Warning: in Carpino outlier rejection, one or more chi-squared values are negative/NaN. Leaving rejection flag"
            << "unchanged" << std::endl;
        }
        return isCurrentlyRejected; // Return unchanged
    }

    // Compute Chi2 for one observation
    double computeChiSquared(const Eigen::MatrixXd& partialsForObservation,
                             const Eigen::VectorXd& residualVector,
                             const Eigen::VectorXd& parameterCorrection,
                             const Eigen::MatrixXd& parameterCovariance,
                             const Eigen::MatrixXd& singleObservationCovariance,
                             const bool isRejected) const
    {
        // Get linear approximation of post-fit residuals
        Eigen::VectorXd residualsVectorPostFit;
        residualsVectorPostFit = residualVector -  partialsForObservation * parameterCorrection;

        // Residual covariance
        Eigen::MatrixXd residualCovariance;
        if( isRejected )
        {
            residualCovariance = singleObservationCovariance + partialsForObservation * parameterCovariance * partialsForObservation.transpose();
        }
        else
        {
            residualCovariance = singleObservationCovariance - partialsForObservation * parameterCovariance * partialsForObservation.transpose();
        }

        // Chi-squared value
        double chiSquared;
        chiSquared = (residualsVectorPostFit.transpose() * residualCovariance.inverse( ) * residualsVectorPostFit).value();

        return chiSquared;

    }

    //! Settings defining the thresholds of the algorithm.
    std::shared_ptr< CarpinoOutlierRejectionSettings > outlierRejectionSettings_;
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
        const std::shared_ptr< OutlierRejectionSettings >& outlierRejectionSettings,
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
            const std::shared_ptr< CarpinoOutlierRejectionSettings > carpinoSettings =
                    std::dynamic_pointer_cast< CarpinoOutlierRejectionSettings >( outlierRejectionSettings );
            if( carpinoSettings == nullptr )
            {
                throw std::runtime_error(
                        "Error when creating outlier rejection object, settings are of type carpino_outlier_rejection, but are not "
                        "CarpinoOutlierRejectionSettings." );
            }
            return std::make_shared< CarpinoOutlierRejection< ObservationScalarType, TimeType > >( carpinoSettings, observationDataset );
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
