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
#include <map>
#include <utility>
#include <tudat/astro/observation_models/observableTypes.h>

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

//! Base class for the settings of an outlier rejection algorithm.
/*!
 * The base class is never instantiated itself, the derived classes are. They hold the configuration settings for the outlier rejection
 * mechanisms. The settings object is used to create the actual OutlierRejection object when the estimation is performed.
 */
class OutlierRejectionSettings
{
public:
    //! Constructor.
    /*!
     * \param outlierRejectionType Algorithm that these settings define.
     * \param firstIterationWithRejection Index of the first estimation iteration (counted from zero) in which observations may be
     * rejected.
     */
    OutlierRejectionSettings( const OutlierRejectionType outlierRejectionType, const int firstIterationWithRejection ):
        outlierRejectionType_( outlierRejectionType ), firstIterationWithRejection_( firstIterationWithRejection )
    {
        if( firstIterationWithRejection_ < 0 )
        {
            throw std::runtime_error(
                    "Error when creating outlier rejection settings, the first iteration with rejection must not be "
                    "negative." );
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

protected:
    //! Algorithm that these settings define, used to select the outlier rejection object that is to be created.
    const OutlierRejectionType outlierRejectionType_;

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
 *
 * In each iteration, at most 20% of the observations that are in the fit at the start of that iteration may be newly rejected. If
 * more observations meet the rejection criterion, only those with the largest chi-squared values are rejected.
 */
class CarpinoOutlierRejectionSettings : public OutlierRejectionSettings
{
public:
    //! Constructor.
    /*!
     * \param chi2Thresholds (recovery, rejection) chi-squared thresholds. An accepted observation is rejected when its chi-squared value
     * exceeds the rejection threshold, and a rejected observation is recovered when its chi-squared value drops below the recovery
     * threshold. Both must be positive, and the recovery threshold must not exceed the rejection threshold.
     * \param firstIterationWithRejection Index of the first estimation iteration in which observations may be
     * rejected. Iterations are counted from zero, so the default value of 1 leaves the first iteration (which uses
     * the a priori parameter values, and can therefore have large residuals for all observations) untouched.
     */
    CarpinoOutlierRejectionSettings( const std::pair< double, double >& chi2Thresholds, const int firstIterationWithRejection = 1 ):
        OutlierRejectionSettings( OutlierRejectionType::carpino_outlier_rejection, firstIterationWithRejection ),
        chi2Thresholds_( chi2Thresholds )
    {
        if( !( chi2Thresholds_.first > 0.0 ) || !( chi2Thresholds_.second > 0.0 ) )
        {
            throw std::runtime_error( "Error when creating Carpino outlier rejection settings, chi-squared thresholds must be positive." );
        }

        if( chi2Thresholds_.first > chi2Thresholds_.second )
        {
            throw std::runtime_error(
                    "Error when creating Carpino outlier rejection settings, the chi-squared recovery threshold must not exceed the "
                    "chi-squared rejection threshold." );
        }
    }

    //! Return the (recovery, rejection) chi-squared thresholds.
    const std::pair< double, double >& getChi2Thresholds( ) const
    {
        return chi2Thresholds_;
    }

    double getChi2RecoveryThreshold( ) const
    {
        return chi2Thresholds_.first;
    }

    double getChi2RejectionThreshold( ) const
    {
        return chi2Thresholds_.second;
    }

protected:
    //! (recovery, rejection) chi-squared thresholds.
    const std::pair< double, double > chi2Thresholds_;
};

//! Settings for the simple outlier rejection algorithm.
/*!
 * This algorithm rejects and recovers observations based on the largest absolute O-C residual of an observation. Thresholds are given as
 * a (recovery, rejection) pair: an accepted observation is rejected when its residual exceeds the rejection threshold, and a rejected
 * observation is recovered when its residual drops below the recovery threshold. One pair is used for all observable types, or one pair
 * per observable type is provided. Units are SI.
 */
class SimpleOutlierRejectionSettings : public OutlierRejectionSettings
{
public:
    //! Constructor.
    /*!
     * \param residualThresholds (recovery, rejection) residual thresholds, used for all observable types. A recovery threshold of 0
     * means rejected observations are never recovered.
     * \param firstIterationWithRejection Index of the first estimation iteration in which observations may be rejected.
     */
    SimpleOutlierRejectionSettings( const std::pair< double, double >& residualThresholds, const int firstIterationWithRejection = 1 ):
        OutlierRejectionSettings( OutlierRejectionType::simple_outlier_rejection, firstIterationWithRejection ),
        residualThresholds_( residualThresholds )
    {
        checkResidualThresholds( residualThresholds_ );
    }

    //! Constructor.
    /*!
     * \param residualThresholdsMap (recovery, rejection) residual thresholds per observable type. An error is thrown during the
     * estimation if an observable type is missing.
     * \param firstIterationWithRejection Index of the first estimation iteration in which observations may be rejected.
     */
    SimpleOutlierRejectionSettings( const std::map< ObservableType, std::pair< double, double > >& residualThresholdsMap,
                                    const int firstIterationWithRejection = 1 ):
        OutlierRejectionSettings( OutlierRejectionType::simple_outlier_rejection, firstIterationWithRejection ),
        residualThresholdsMap_( residualThresholdsMap ), residualThresholds_( -1.0, -1.0 )
    {
        if( residualThresholdsMap_.empty( ) )
        {
            throw std::runtime_error(
                    "Error when creating simple outlier rejection settings, the map of residual thresholds per observable type is "
                    "empty." );
        }
        for( const auto& thresholdsPerObservable : residualThresholdsMap_ )
        {
            checkResidualThresholds( thresholdsPerObservable.second );
        }
    }

    const std::map< ObservableType, std::pair< double, double > >& getResidualThresholdsMap( ) const
    {
        return residualThresholdsMap_;
    }

    const std::pair< double, double >& getResidualThresholds( ) const
    {
        return residualThresholds_;
    }

protected:
    //! Check that 0 <= recovery <= rejection, and rejection > 0.
    static void checkResidualThresholds( const std::pair< double, double >& residualThresholds )
    {
        if( !( residualThresholds.second > 0.0 ) )
        {
            throw std::runtime_error( "Error when creating simple outlier rejection settings, the rejection threshold must be positive." );
        }
        if( !( residualThresholds.first >= 0.0 ) || residualThresholds.first > residualThresholds.second )
        {
            throw std::runtime_error(
                    "Error when creating simple outlier rejection settings, the recovery threshold must be between 0 and the rejection "
                    "threshold." );
        }
    }

    //! (recovery, rejection) residual thresholds per observable type; empty if one pair is used for all types.
    const std::map< ObservableType, std::pair< double, double > > residualThresholdsMap_;

    //! (recovery, rejection) residual thresholds for all observable types; (-1, -1) if the map is used.
    const std::pair< double, double > residualThresholds_;
};

//! Function to create settings for the outlier rejection algorithm of Carpino et al. (2003).
inline std::shared_ptr< OutlierRejectionSettings > carpinoOutlierRejectionSettings( const std::pair< double, double >& chi2Thresholds,
                                                                                    const int firstIterationWithRejection = 1 )
{
    return std::make_shared< CarpinoOutlierRejectionSettings >( chi2Thresholds, firstIterationWithRejection );
}

//! Function to create settings for the simple outlier rejection algorithm, with one (recovery, rejection) pair for all observables.
inline std::shared_ptr< OutlierRejectionSettings > simpleOutlierRejectionSettings( const std::pair< double, double >& residualThresholds,
                                                                                   const int firstIterationWithRejection = 1 )
{
    return std::make_shared< SimpleOutlierRejectionSettings >( residualThresholds, firstIterationWithRejection );
}

//! Function to create settings for the simple outlier rejection algorithm, with a (recovery, rejection) pair per observable type.
inline std::shared_ptr< OutlierRejectionSettings > simpleOutlierRejectionSettings(
        const std::map< ObservableType, std::pair< double, double > >& residualThresholdsMap,
        const int firstIterationWithRejection = 1 )
{
    return std::make_shared< SimpleOutlierRejectionSettings >( residualThresholdsMap, firstIterationWithRejection );
}

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_OUTLIERREJECTIONSETTINGS_H
