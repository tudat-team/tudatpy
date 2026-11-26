#ifndef TUDAT_ANGULAR_OBSERVATION_CORRECTION_SETTINGS_H
#define TUDAT_ANGULAR_OBSERVATION_CORRECTION_SETTINGS_H

#include <algorithm>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace tudat
{
namespace data
{

//! Plain settings for angular corrections that require an environment during dataset conversion.
struct AngularObservationCorrectionSettings {
    AngularObservationCorrectionSettings( const std::vector< std::string >& lightDeflectionBodies = {},
                                          const bool photocenterCorrection = false,
                                          const std::vector< unsigned int >& observationIndices = {} ):
        lightDeflectionBodies_( lightDeflectionBodies ), photocenterCorrection_( photocenterCorrection ),
        observationIndices_( observationIndices )
    {
        if( std::set< unsigned int >( observationIndices_.begin( ), observationIndices_.end( ) ).size( ) != observationIndices_.size( ) )
        {
            throw std::runtime_error( "Observation indices of angular corrections must be unique." );
        }
        for( const auto& name : lightDeflectionBodies_ )
        {
            if( name.empty( ) )
            {
                throw std::runtime_error( "Light-deflection body names must not be empty." );
            }
        }
    }

    bool hasCorrections( ) const
    {
        return !lightDeflectionBodies_.empty( ) || photocenterCorrection_;
    }

    //! Remove a single observation, return whether corrections for any observation remain
    bool removeObservation( const unsigned int observationIndex )
    {
        if( observationIndices_.empty( ) )
        {
            return true;
        }
        observationIndices_.erase( std::remove( observationIndices_.begin( ), observationIndices_.end( ), observationIndex ),
                                   observationIndices_.end( ) );
        for( unsigned int& index : observationIndices_ )
        {
            if( index > observationIndex )
            {
                index--;
            }
        }
        return !observationIndices_.empty( );
    }

    std::vector< std::string > lightDeflectionBodies_;
    bool photocenterCorrection_;

    //! Indices of the observations that are corrected (empty for all observations)
    std::vector< unsigned int > observationIndices_;
};

}  // namespace data
}  // namespace tudat

#endif
