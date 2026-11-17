#ifndef TUDAT_ANGULAR_OBSERVATION_CORRECTION_SETTINGS_H
#define TUDAT_ANGULAR_OBSERVATION_CORRECTION_SETTINGS_H

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
                                          const bool photocenterCorrection = false ):
        lightDeflectionBodies_( lightDeflectionBodies ), photocenterCorrection_( photocenterCorrection )
    {
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

    std::vector< std::string > lightDeflectionBodies_;
    bool photocenterCorrection_;
};

}  // namespace data
}  // namespace tudat

#endif
