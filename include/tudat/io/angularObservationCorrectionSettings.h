#ifndef TUDAT_ANGULAR_OBSERVATION_CORRECTION_SETTINGS_H
#define TUDAT_ANGULAR_OBSERVATION_CORRECTION_SETTINGS_H

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Core>

namespace tudat
{
namespace data
{

//! Plain settings for angular corrections that require an environment during dataset conversion.
struct AngularObservationCorrectionSettings {
    AngularObservationCorrectionSettings( const std::vector< std::string >& lightDeflectionBodies = {},
                                          const Eigen::VectorXd& photocenterBodyDimensions = Eigen::VectorXd( ) ):
        lightDeflectionBodies_( lightDeflectionBodies ), photocenterBodyDimensions_( photocenterBodyDimensions )
    {
        const auto size = photocenterBodyDimensions_.size( );
        if( size != 0 && size != 1 && size != 3 )
        {
            throw std::runtime_error( "Photocenter dimensions must contain a radius or three ellipsoid semi-axes [m]." );
        }
        for( Eigen::Index i = 0; i < size; ++i )
        {
            if( !std::isfinite( photocenterBodyDimensions_( i ) ) || photocenterBodyDimensions_( i ) <= 0.0 )
            {
                throw std::runtime_error( "Photocenter dimensions must be finite and positive [m]." );
            }
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
        return !lightDeflectionBodies_.empty( ) || photocenterBodyDimensions_.size( ) != 0;
    }

    std::vector< std::string > lightDeflectionBodies_;
    Eigen::VectorXd photocenterBodyDimensions_;
};

}  // namespace data
}  // namespace tudat

#endif
