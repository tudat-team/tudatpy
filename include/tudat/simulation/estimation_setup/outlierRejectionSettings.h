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
#include <tudat/astro/observation_models/observableTypes.h>

namespace tudat
{

namespace simulation_setup
{

using ObservableType = tudat::observation_models::ObservableType;

//! Outlier rejection algorithms that can be used during an estimation.
enum class OutlierRejectionType
{
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
    OutlierRejectionSettings( const OutlierRejectionType outlierRejectionType ): outlierRejectionType_( outlierRejectionType ) { }

    //! Destructor.
    virtual ~OutlierRejectionSettings( ) = default;

    //! Return the type of the outlier rejection algorithm
    OutlierRejectionType getOutlierRejectionType( ) const
    {
        return outlierRejectionType_;
    }

protected:
    //! Algorithm that these settings define, used to select the outlier rejection object that is to be created.
    const OutlierRejectionType outlierRejectionType_;
};

//! Settings for the outlier rejection algorithm of Carpino et al. (2003).
/*!
 * This algorithm rejects and recovers observations based on a per-observation chi-squared value. For scalar observations, this is the ratio
 * of the residual to the residual uncertainty. For vector observations, this is the Mahalanobis distance of the residual vector. Note that the
 * residual covariance, that is used in this algorithm, is different from the typical observation covariance. The former also takes into account
 * whether an observation was used in the least-squares inversion or not. An observation that is currently used in the estimation
 * is rejected when its chi-squared value exceeds the rejection threshold. An observation that was rejected in an
 * earlier iteration is recovered when its chi-squared value drops below the recovery threshold. Using a recovery
 * threshold that is lower than the rejection threshold prevents observations from oscillating between the rejected
 * and the accepted state in successive iterations.
 */
class CarpinoOutlierRejectionSettings : public OutlierRejectionSettings
{
public:
    //! Constructor.
    /*!
     * \param chi2RejectionThreshold Chi-squared value above which an accepted observation is rejected.
     * \param chi2RecoveryThreshold Chi-squared value below which a rejected observation is recovered. Must be smaller
     * than the rejection threshold.
     * \param maximumRejectedFraction Maximum fraction (0 to 1) of all observations that may be in the rejected state.
     * \param firstIterationWithRejection Index of the first estimation iteration in which observations may be
     * rejected. Iterations are counted from zero, so the default value of 1 leaves the first iteration (which uses
     * the a priori parameter values, and can therefore have large residuals for all observations) untouched.
     */
    CarpinoOutlierRejectionSettings( const double chi2RejectionThreshold = 9.0,
                                     const double chi2RecoveryThreshold = 8.0,
                                     const double maximumRejectedFraction = 0.25,
                                     const int firstIterationWithRejection = 1 ):
        OutlierRejectionSettings( OutlierRejectionType::carpino_outlier_rejection ),
        chi2RejectionThreshold_( chi2RejectionThreshold ), chi2RecoveryThreshold_( chi2RecoveryThreshold ),
        maximumRejectedFraction_( maximumRejectedFraction ), firstIterationWithRejection_( firstIterationWithRejection )
    {
        if( chi2RejectionThreshold_ <= 0.0 || chi2RecoveryThreshold_ <= 0.0 )
        {
            throw std::runtime_error( "Error when creating Carpino outlier rejection settings, chi-squared thresholds must be positive." );
        }

        if( chi2RecoveryThreshold_ >= chi2RejectionThreshold_ )
        {
            throw std::runtime_error(
                    "Error when creating Carpino outlier rejection settings, the chi-squared recovery threshold must be smaller than the "
                    "chi-squared rejection threshold." );
        }

        if( maximumRejectedFraction_ < 0.0 || maximumRejectedFraction_ > 1.0 )
        {
            throw std::runtime_error(
                    "Error when creating Carpino outlier rejection settings, the maximum rejected fraction must be between 0 and 1." );
        }

        if( firstIterationWithRejection_ < 0 )
        {
            throw std::runtime_error(
                    "Error when creating Carpino outlier rejection settings, the first iteration with rejection must not be negative." );
        }
    }

    double getChi2RejectionThreshold( ) const
    {
        return chi2RejectionThreshold_;
    }

    double getChi2RecoveryThreshold( ) const
    {
        return chi2RecoveryThreshold_;
    }

    double getMaximumRejectedFraction( ) const
    {
        return maximumRejectedFraction_;
    }

    int getFirstIterationWithRejection( ) const
    {
        return firstIterationWithRejection_;
    }

protected:
    //! Chi-squared value above which an accepted observation is rejected.
    const double chi2RejectionThreshold_;

    //! Chi-squared value below which a rejected observation is recovered.
    const double chi2RecoveryThreshold_;

    //! Maximum fraction of all observations that may be in the rejected state.
    const double maximumRejectedFraction_;

    //! Index of the first estimation iteration in which observations may be rejected.
    const int firstIterationWithRejection_;
};

//! Create settings for the Simple outlier rejection algorithm
/*!
 * This algorithm rejects and recovers observations based on the absolute value of the O-C residuals. If a residual exceeds the threshold,
 * the observation is rejected. There are two options: a scalar residual threshold is provided, which holds for all types of observation
 * used in the estimation, or a map (dictionary) of thresholds is provided per observable type. Units are assumed SI.
*/
class SimpleOutlierRejectionSettings : public OutlierRejectionSettings
{
public:
    //! Constructor
    /*!
     *\param maximumAllowedResidualValue: the absolute threshold for residuals
     *\param firstIterationWithRejection: iteration at which outlier rejection kicks in
     *\param allowRestore: bool if observations are allowed to be restored in later iterations or not
     */
    SimpleOutlierRejectionSettings(
        const double maximumAllowedResidualValue,
        const int firstIterationWithRejection = 1,
        const bool allowRestore = true) : OutlierRejectionSettings( OutlierRejectionType::simple_outlier_rejection ),
            firstIterationWithRejection_( firstIterationWithRejection ), allowRestore_( allowRestore ),
             maximumAllowedResidualValue_( maximumAllowedResidualValue )
    {
        if(firstIterationWithRejection < 0)
        {
            throw std::runtime_error("Error when creating simple outlier rejection settings, the first iteration with rejection must not be negative.");
        }
        if(maximumAllowedResidualValue <= 0.0)
        {
            throw std::runtime_error("Error when creating simple outlier rejection settings. The maximum allowed residual value must be a positive value");
        }
    }

    //! Constructor
    /*!
     *\param maximumAllowedResidualValueMap: absolute residual threshold for each observable type that may be used in the estimation. Will
     *throw an error downstream if an entry is missing. Units are assumed SI
     *\param firstIterationWithRejection: iteration at which outlier rejection kicks in
     *\param allowRestore: bool if observations are allowed to be restored in later iterations or not
     */
    SimpleOutlierRejectionSettings(
        const std::map<tudat::observation_models::ObservableType, double>& maximumAllowedResidualValueMap,
        const int firstIterationWithRejection = 1,
        const bool allowRestore = true) : OutlierRejectionSettings( OutlierRejectionType::simple_outlier_rejection ),
            firstIterationWithRejection_( firstIterationWithRejection ), allowRestore_( allowRestore ),
            maximumAllowedResidualValueMap_(maximumAllowedResidualValueMap), maximumAllowedResidualValue_( -1.0 )
    {
        if(firstIterationWithRejection < 0)
        {
            throw std::runtime_error("Error when creating simple outlier rejection settings, the first iteration with rejection must not be negative");
        }
        if(maximumAllowedResidualValueMap.empty())
        {
            throw std::runtime_error("Error when creating simple outlier rejection settings, the map for the maximum residual values per "
                                     "observable type was left empty");
        }
        for(const auto& pair : maximumAllowedResidualValueMap)
        {
            if(pair.second <= 0.0)
            {
                throw std::runtime_error("Error when creating simple outlier rejection settings. The maximum allowed residual value must be a positive value");
            }
        }
    }

    const std::map<ObservableType, double>& getMaximumAllowedResidualValueMap( ) const
    {
        return maximumAllowedResidualValueMap_;
    }

    double getMaximumAllowedResidualValue( ) const
    {
        return maximumAllowedResidualValue_;
    }

    int getFirstIterationWithRejection( ) const
    {
        return firstIterationWithRejection_;
    }

    bool getAllowRestore( ) const
    {
        return allowRestore_;
    }


protected:
    // First iteration at which to start outlier rejection process
    const int firstIterationWithRejection_;

    // If true, allow rejected observations to be restored in later iterations
    const bool allowRestore_;

    // A map of the maximum allowed residual value per observable type
    const std::map<ObservableType, double> maximumAllowedResidualValueMap_;

    // One maximum allowed value per observable type
    const double maximumAllowedResidualValue_;
};

//! Function to create settings for the outlier rejection algorithm of Carpino et al. (2003).
inline std::shared_ptr< OutlierRejectionSettings > carpinoOutlierRejectionSettings( const double chi2RejectionThreshold = 9.0,
                                                                                    const double chi2RecoveryThreshold = 8.0,
                                                                                    const double maximumRejectedFraction = 0.25,
                                                                                    const int firstIterationWithRejection = 1 )
{
    return std::make_shared< CarpinoOutlierRejectionSettings >(
            chi2RejectionThreshold, chi2RecoveryThreshold, maximumRejectedFraction, firstIterationWithRejection );
}

//! Function to create settings the simple residual based outlier rejection settings
inline std::shared_ptr< OutlierRejectionSettings > simpleOutlierRejectionSettings( const double maximumAllowedResidualValue,
                                                                                    const int firstIterationWithRejection = 1,
                                                                                    const bool allowRestore = true)
{
    return std::make_shared< SimpleOutlierRejectionSettings>(maximumAllowedResidualValue, firstIterationWithRejection, allowRestore );
}

inline std::shared_ptr< OutlierRejectionSettings > simpleOutlierRejectionSettings( const std::map<ObservableType, double>& maximumAllowedResidualValueMap_,
                                                                                    const int firstIterationWithRejection = 1,
                                                                                    const bool allowRestore = true)
{
    return std::make_shared< SimpleOutlierRejectionSettings>(maximumAllowedResidualValueMap_, firstIterationWithRejection, allowRestore );
}



}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_OUTLIERREJECTIONSETTINGS_H
