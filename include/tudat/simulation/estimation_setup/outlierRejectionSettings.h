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

namespace tudat
{

namespace simulation_setup
{

//! Outlier rejection algorithms that can be used during an estimation.
enum class OutlierRejectionType { carpino_outlier_rejection };

inline std::string getOutlierRejectionTypeString( const OutlierRejectionType outlierRejectionType )
{
    switch( outlierRejectionType )
    {
        case OutlierRejectionType::carpino_outlier_rejection:
            return "carpino_outlier_rejection";
        default:
            throw std::runtime_error( "Error when retrieving outlier rejection type string, type not recognized." );
    }
}

//! Base class for the settings of an outlier rejection algorithm.
/*!
 * Objects of this type only carry the user's choices; they perform no computation. They are provided to an
 * EstimationInput object, and are converted into the corresponding OutlierRejection object (which does the actual
 * work) when the estimation starts. Each algorithm defines its own derived settings class, holding the parameters
 * that are meaningful for that algorithm only.
 */
class OutlierRejectionSettings
{
public:
    OutlierRejectionSettings( const OutlierRejectionType outlierRejectionType ): outlierRejectionType_( outlierRejectionType ) { }

    //! Destructor.
    /*!
     * This destructor is virtual because objects of derived types are stored and deleted through a pointer to this
     * base class. Without the virtual keyword, deleting such a pointer would only destroy the base class part of the
     * object, and the members added by the derived class would be leaked.
     */
    virtual ~OutlierRejectionSettings( ) = default;

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
 * The algorithm computes a chi-squared value for each observation from its residual and the covariance of that
 * residual, and compares this value against two thresholds. An observation that is currently used in the estimation
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

//! Function to create settings for the outlier rejection algorithm of Carpino et al. (2003).
inline std::shared_ptr< OutlierRejectionSettings > carpinoOutlierRejectionSettings( const double chi2RejectionThreshold = 9.0,
                                                                                    const double chi2RecoveryThreshold = 8.0,
                                                                                    const double maximumRejectedFraction = 0.25,
                                                                                    const int firstIterationWithRejection = 1 )
{
    return std::make_shared< CarpinoOutlierRejectionSettings >(
            chi2RejectionThreshold, chi2RecoveryThreshold, maximumRejectedFraction, firstIterationWithRejection );
}

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_OUTLIERREJECTIONSETTINGS_H
