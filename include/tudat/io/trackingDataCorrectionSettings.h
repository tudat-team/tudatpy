/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_TRACKING_DATA_CORRECTION_SETTINGS_H
#define TUDAT_TRACKING_DATA_CORRECTION_SETTINGS_H

#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "tudat/io/angularObservationCorrectionSettings.h"

namespace tudat
{

namespace data
{

//! Corrections that are known before creating an observation dataset, one per observation
template< typename ObservationScalarType = double >
class ExplicitCorrectionSettings
{
public:
    ExplicitCorrectionSettings( const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& corrections ):
        corrections_( corrections )
    {}

    void removeObservation( const unsigned int observationIndex )
    {
        corrections_.erase( corrections_.begin( ) + observationIndex );
    }

    //! Corrections that are added to the observations, one per observation
    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > corrections_;
};

//! Settings that define the corrections to the observations in a TrackingData object (applied when creating an observation dataset
//! with corrections applied)
template< typename ObservationScalarType = double >
class TrackingDataCorrectionSettings
{
public:
    //! Constructor, for the observations of a TrackingData object
    /*!
     * \param numberOfObservations Number of observations
     * \param singleObservationSize Size of a single observation
     * \param observableType Observable type of the observations (as in TrackingData)
     * \param referenceLinkEnd Reference link end of the observations (as in TrackingData)
     */
    TrackingDataCorrectionSettings( const unsigned int numberOfObservations,
                                    const unsigned int singleObservationSize,
                                    const std::string& observableType,
                                    const std::string& referenceLinkEnd ):
        numberOfObservations_( numberOfObservations ), singleObservationSize_( singleObservationSize ), observableType_( observableType ),
        referenceLinkEnd_( referenceLinkEnd )
    {}

    //! Set corrections to the observations (optional)
    void setObservationCorrections( const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observationCorrections )
    {
        // Check size consistency (for the total number of observations)
        if( observationCorrections.size( ) != numberOfObservations_ )
        {
            throw std::runtime_error( "Error when adding observation corrections to tracking data object, size of corrections (" +
                                      std::to_string( observationCorrections.size( ) ) + ") does not match number of observations (" +
                                      std::to_string( numberOfObservations_ ) + ")." );
        }

        // Check size consistency (for a single observable)
        for( auto correction : observationCorrections )
        {
            if( correction.size( ) != singleObservationSize_ )
            {
                throw std::runtime_error( "Error when adding observation corrections to tracking data object, size of single correction (" +
                                          std::to_string( correction.size( ) ) + ") does not match single observable size (should be " +
                                          std::to_string( singleObservationSize_ ) + ")." );
            }
        }

        // Warn only after the complete replacement has passed validation. Rejected
        // replacements must leave the existing corrections untouched.
        if( explicitCorrections_ != nullptr )
        {
            std::cerr << "Warning when adding observation corrections to tracking data object, corrections already existed and are "
                         "overwritten ."
                      << std::endl;
        }

        // If all sizes are consistent, store observation corrections
        explicitCorrections_ = std::make_shared< ExplicitCorrectionSettings< ObservationScalarType > >( observationCorrections );
    }

    //! Function that (re-)sets a single observation correction
    void setSingleObservationCorrection( const unsigned int index,
                                         const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& observationCorrection )
    {
        // Check that observation corrections are already defined
        if( explicitCorrections_ == nullptr )
        {
            throw std::runtime_error(
                    "Error when resetting single observation correction in TrackingData object, corrections not yet defined." );
        }

        // Check that the index of the observation correction that needs overwriting does not exceed the size of the observation vector
        if( index >= numberOfObservations_ )
        {
            throw std::runtime_error(
                    "Error when resetting single observation correction in TrackingData object, index exceeds number of observations "
                    "contained in TrackingData object." );
        }

        // Check size consistency
        if( observationCorrection.size( ) != singleObservationSize_ )
        {
            throw std::runtime_error( "Error when resetting single observation correction in TrackingData object, new correction size (" +
                                      std::to_string( observationCorrection.size( ) ) + ") inconsistent with observation size (" +
                                      std::to_string( singleObservationSize_ ) + ")." );
        }

        // Overwrites specific correction entry
        explicitCorrections_->corrections_.at( index ) = observationCorrection;
    }

    //! Function that return a vector of observation corrections (empty of observation corrections are not provided)
    const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& getObservationCorrections( ) const
    {
        if( explicitCorrections_ != nullptr )
        {
            return explicitCorrections_->corrections_;
        }
        static const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > emptyCorrections;
        return emptyCorrections;
    }

    //! Store corrections to be evaluated when constructing an observation dataset.
    void setObservationCorrectionSettings( const AngularObservationCorrectionSettings& settings )
    {
        if( observableType_ != "AngularPosition" || singleObservationSize_ != 2 || referenceLinkEnd_ != "receiver" )
        {
            throw std::runtime_error( "Angular correction settings require two-component reception-time angular observations." );
        }
        angularCorrectionSettings_ = std::make_shared< AngularObservationCorrectionSettings >( settings );
    }

    //! Return the angular correction settings (nullptr if not provided)
    std::shared_ptr< AngularObservationCorrectionSettings > getObservationCorrectionSettings( ) const
    {
        return angularCorrectionSettings_;
    }

    bool hasObservationCorrectionSettings( ) const
    {
        return angularCorrectionSettings_ != nullptr && angularCorrectionSettings_->hasCorrections( );
    }

    //! Remove a single observation from all settings
    void removeObservation( const unsigned int observationIndex )
    {
        if( explicitCorrections_ != nullptr )
        {
            explicitCorrections_->removeObservation( observationIndex );
        }
        numberOfObservations_--;
    }

private:
    //! Number of observations
    unsigned int numberOfObservations_;

    //! Size of a single observation
    unsigned int singleObservationSize_;

    //! Observable type of the observations (as in TrackingData)
    std::string observableType_;

    //! Reference link end of the observations (as in TrackingData)
    std::string referenceLinkEnd_;

    //! Corrections that are known before creating the observation dataset (nullptr if none)
    std::shared_ptr< ExplicitCorrectionSettings< ObservationScalarType > > explicitCorrections_;

    //! Angular corrections computed from the environment when creating the observation dataset (nullptr if none)
    std::shared_ptr< AngularObservationCorrectionSettings > angularCorrectionSettings_;
};

}  // namespace data

}  // namespace tudat

#endif  // TUDAT_TRACKING_DATA_CORRECTION_SETTINGS_H
