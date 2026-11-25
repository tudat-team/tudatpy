/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OUTLIERREJECTIONSETTINGS_H
#define TUDAT_OUTLIERREJECTIONSETTINGS_H

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <tudat/astro/observation_models/observableTypes.h>
#include "tudat/simulation/estimation_setup/observationCondition.h"

namespace tudat
{

namespace simulation_setup
{

using ObservableType = tudat::observation_models::ObservableType;

//! Outlier rejection algorithms that can be used during an estimation.
enum class OutlierRejectionType {
    carpino_outlier_rejection,
    simple_outlier_rejection,
};

inline std::string getOutlierRejectionTypeString( const OutlierRejectionType outlierRejectionType )
{
    switch( outlierRejectionType )
    {
        case OutlierRejectionType::carpino_outlier_rejection:
            return "carpino_outlier_rejection";
        case OutlierRejectionType::simple_outlier_rejection:
            return "simple_outlier_rejection";
        default:
            throw std::runtime_error( "Error when retrieving outlier rejection type string, type not recognized." );
    }
}

//! Check that a (recovery, rejection) pair of thresholds is valid: 0 <= recovery <= rejection.
inline void checkOutlierRejectionThresholds( const std::pair< double, double >& thresholds,
                                             const OutlierRejectionType outlierRejectionType )
{
    if( !( thresholds.first >= 0.0 ) || !( thresholds.second >= 0.0 ) )
    {
        throw std::runtime_error( "Error when creating " + getOutlierRejectionTypeString( outlierRejectionType ) +
                                  " settings, the recovery and rejection thresholds must not be negative." );
    }
    if( thresholds.first > thresholds.second )
    {
        throw std::runtime_error( "Error when creating " + getOutlierRejectionTypeString( outlierRejectionType ) +
                                  " settings, the recovery threshold must not exceed the rejection threshold." );
    }
}

//! Base class for the settings of an outlier rejection algorithm.
/*!
 * The base class is never instantiated itself, the derived classes are. They hold the configuration settings for the outlier rejection
 * mechanisms. The settings object is used to create the actual OutlierRejection object when the estimation is performed.
 *
 * Every algorithm uses a (recovery, rejection) pair of thresholds per observation selection condition. Each condition defines a group of
 * observations that share one pair of thresholds. When the estimation starts, every condition is evaluated once for every observation.
 * Observations that match no condition are not affected by the outlier rejection, and an error is thrown if an observation matches more
 * than one condition. Since the conditions are evaluated once, conditions that depend on the residuals or the rejection status use the
 * state of the observations before the estimation.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
class OutlierRejectionSettings
{
public:
    using Condition = observation_models::ObservationSelectionCondition< ObservationScalarType, TimeType >;
    using ThresholdsPerCondition = std::vector< std::pair< Condition, std::pair< double, double > > >;

    //! Constructor.
    /*!
     * \param outlierRejectionType Algorithm that these settings define.
     * \param thresholdsPerCondition (recovery, rejection) thresholds per observation selection condition. Must not be empty, and
     * for each pair 0 <= recovery <= rejection must hold.
     * \param firstIterationWithRejection Index of the first estimation iteration (counted from zero) in which observations may be
     * rejected.
     */
    OutlierRejectionSettings( const OutlierRejectionType outlierRejectionType,
                              const ThresholdsPerCondition& thresholdsPerCondition,
                              const int firstIterationWithRejection ):
        outlierRejectionType_( outlierRejectionType ), thresholdsPerCondition_( thresholdsPerCondition ),
        firstIterationWithRejection_( firstIterationWithRejection )
    {
        if( firstIterationWithRejection_ < 0 )
        {
            throw std::runtime_error(
                    "Error when creating outlier rejection settings, the first iteration with rejection must not be "
                    "negative." );
        }
        if( thresholdsPerCondition_.empty( ) )
        {
            throw std::runtime_error( "Error when creating " + getOutlierRejectionTypeString( outlierRejectionType ) +
                                      " settings, the list of thresholds per condition is empty." );
        }
        for( const auto& thresholdsForCondition : thresholdsPerCondition_ )
        {
            checkOutlierRejectionThresholds( thresholdsForCondition.second, outlierRejectionType );
        }
    }

    //! Destructor.
    virtual ~OutlierRejectionSettings( ) = default;

    //! Return the type of the outlier rejection algorithm
    OutlierRejectionType getOutlierRejectionType( ) const
    {
        return outlierRejectionType_;
    }

    int getFirstIterationWithRejection( ) const
    {
        return firstIterationWithRejection_;
    }

    //! Return the (recovery, rejection) thresholds per observation selection condition.
    const ThresholdsPerCondition& getThresholdsPerCondition( ) const
    {
        return thresholdsPerCondition_;
    }

    //! Return the observation selection conditions, in the order in which they were provided.
    std::vector< Condition > getConditions( ) const
    {
        std::vector< Condition > conditions;
        for( const auto& thresholdsForCondition : thresholdsPerCondition_ )
        {
            conditions.push_back( thresholdsForCondition.first );
        }
        return conditions;
    }

    //! Return the (recovery, rejection) thresholds of the condition with the given index.
    const std::pair< double, double >& getThresholds( const unsigned int conditionIndex ) const
    {
        return thresholdsPerCondition_.at( conditionIndex ).second;
    }

protected:
    //! Algorithm that these settings define, used to select the outlier rejection object that is to be created.
    const OutlierRejectionType outlierRejectionType_;

    //! (recovery, rejection) thresholds per observation selection condition.
    const ThresholdsPerCondition thresholdsPerCondition_;

    //! Index of the first estimation iteration in which observations may be rejected.
    const int firstIterationWithRejection_;
};

//! Settings for the outlier rejection algorithm of Carpino et al. (2003).
/*!
 * This algorithm rejects and recovers observations based on a per-observation chi-squared value. For scalar observations, this is the
 * square of the ratio of the residual to the residual uncertainty. For vector observations, this is the squared Mahalanobis distance of the
 * residual vector. Note that the residual covariance, that is used in this algorithm, is different from the typical observation
 * covariance. The former also takes into account whether an observation was used in the least-squares inversion or not. An observation
 * that is currently used in the estimation is rejected when its chi-squared value exceeds the rejection threshold. An observation that
 * was rejected in an earlier iteration is recovered when its chi-squared value drops below the recovery threshold. Thresholds are
 * given as a (recovery, rejection) pair. Using a recovery threshold that is lower than the rejection threshold prevents observations
 * from oscillating between the rejected and the accepted state in successive iterations; the two may also be equal.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
class CarpinoOutlierRejectionSettings : public OutlierRejectionSettings< ObservationScalarType, TimeType >
{
public:
    using Base = OutlierRejectionSettings< ObservationScalarType, TimeType >;
    using typename Base::Condition;
    using typename Base::ThresholdsPerCondition;

    //! Constructor.
    /*!
     * \param chi2ThresholdsPerCondition (recovery, rejection) chi-squared thresholds per observation selection condition. An accepted
     * observation is rejected when its chi-squared value exceeds the rejection threshold, and a rejected observation is recovered when
     * its chi-squared value drops below the recovery threshold. Both must not be negative, and the recovery threshold must not exceed
     * the rejection threshold.
     * \param firstIterationWithRejection Index of the first estimation iteration in which observations may be
     * rejected. Iterations are counted from zero, so the default value of 1 leaves the first iteration (which uses
     * the a priori parameter values, and can therefore have large residuals for all observations) untouched.
     */
    CarpinoOutlierRejectionSettings( const ThresholdsPerCondition& chi2ThresholdsPerCondition, const int firstIterationWithRejection = 1 ):
        Base( OutlierRejectionType::carpino_outlier_rejection, chi2ThresholdsPerCondition, firstIterationWithRejection )
    {}

    //! Constructor, with one (recovery, rejection) pair of chi-squared thresholds for all observations.
    CarpinoOutlierRejectionSettings( const std::pair< double, double >& chi2Thresholds, const int firstIterationWithRejection = 1 ):
        CarpinoOutlierRejectionSettings( ThresholdsPerCondition{ { Condition::all( ), chi2Thresholds } }, firstIterationWithRejection )
    {}
};

//! Settings for the simple outlier rejection algorithm.
/*!
 * This algorithm rejects and recovers observations based on the largest absolute O-C residual of an observation. Thresholds are given as
 * a (recovery, rejection) pair: an accepted observation is rejected when its residual exceeds the rejection threshold, and a rejected
 * observation is recovered when its residual drops below the recovery threshold. One pair is used for all observations, or one pair
 * per observation selection condition is provided (see OutlierRejectionSettings). Units are assumed SI.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
class SimpleOutlierRejectionSettings : public OutlierRejectionSettings< ObservationScalarType, TimeType >
{
public:
    using Base = OutlierRejectionSettings< ObservationScalarType, TimeType >;
    using typename Base::Condition;
    using typename Base::ThresholdsPerCondition;

    //! Constructor.
    /*!
     * \param residualThresholdsPerCondition (recovery, rejection) residual thresholds per observation selection condition. Both must not
     * be negative, and the recovery threshold must not exceed the rejection threshold. A recovery threshold of 0 means rejected
     * observations are never recovered.
     * \param firstIterationWithRejection Index of the first estimation iteration in which observations may be rejected.
     */
    SimpleOutlierRejectionSettings( const ThresholdsPerCondition& residualThresholdsPerCondition,
                                    const int firstIterationWithRejection = 1 ):
        Base( OutlierRejectionType::simple_outlier_rejection, residualThresholdsPerCondition, firstIterationWithRejection )
    {}

    //! Constructor, with one (recovery, rejection) pair of residual thresholds for all observations.
    SimpleOutlierRejectionSettings( const std::pair< double, double >& residualThresholds, const int firstIterationWithRejection = 1 ):
        SimpleOutlierRejectionSettings( ThresholdsPerCondition{ { Condition::all( ), residualThresholds } }, firstIterationWithRejection )
    {}
};

//! Function to create settings for the outlier rejection algorithm of Carpino et al. (2003), with one (recovery, rejection) pair for all
//! observations.
template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< OutlierRejectionSettings< ObservationScalarType, TimeType > > carpinoOutlierRejectionSettings(
        const std::pair< double, double >& chi2Thresholds,
        const int firstIterationWithRejection = 1 )
{
    return std::make_shared< CarpinoOutlierRejectionSettings< ObservationScalarType, TimeType > >( chi2Thresholds,
                                                                                                   firstIterationWithRejection );
}

//! Function to create settings for the outlier rejection algorithm of Carpino et al. (2003), with a (recovery, rejection) pair per
//! observation selection condition.
template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< OutlierRejectionSettings< ObservationScalarType, TimeType > > carpinoOutlierRejectionSettings(
        const typename CarpinoOutlierRejectionSettings< ObservationScalarType, TimeType >::ThresholdsPerCondition&
                chi2ThresholdsPerCondition,
        const int firstIterationWithRejection = 1 )
{
    return std::make_shared< CarpinoOutlierRejectionSettings< ObservationScalarType, TimeType > >( chi2ThresholdsPerCondition,
                                                                                                   firstIterationWithRejection );
}

//! Function to create settings for the simple outlier rejection algorithm, with one (recovery, rejection) pair for all observations.
template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< OutlierRejectionSettings< ObservationScalarType, TimeType > > simpleOutlierRejectionSettings(
        const std::pair< double, double >& residualThresholds,
        const int firstIterationWithRejection = 1 )
{
    return std::make_shared< SimpleOutlierRejectionSettings< ObservationScalarType, TimeType > >( residualThresholds,
                                                                                                  firstIterationWithRejection );
}

//! Function to create settings for the simple outlier rejection algorithm, with a (recovery, rejection) pair per observation selection
//! condition.
template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< OutlierRejectionSettings< ObservationScalarType, TimeType > > simpleOutlierRejectionSettings(
        const typename SimpleOutlierRejectionSettings< ObservationScalarType, TimeType >::ThresholdsPerCondition&
                residualThresholdsPerCondition,
        const int firstIterationWithRejection = 1 )
{
    return std::make_shared< SimpleOutlierRejectionSettings< ObservationScalarType, TimeType > >( residualThresholdsPerCondition,
                                                                                                  firstIterationWithRejection );
}

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_OUTLIERREJECTIONSETTINGS_H
