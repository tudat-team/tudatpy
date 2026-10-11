#ifndef TUDAT_BARYCENTER_ORIGIN_ACCELERATION_H
#define TUDAT_BARYCENTER_ORIGIN_ACCELERATION_H

#include "tudat/astro/basic_astro/accelerationModel.h"

namespace tudat
{
namespace basic_astrodynamics
{

//! Acceleration of coordinates relative to a prescribed member of a propagated system.
class BarycenterOriginAcceleration : public AccelerationModel3d
{
public:
    BarycenterOriginAcceleration( const std::string& barycenter,
                                  const SingleBodyAccelerationMap& externalAccelerations,
                                  const std::function< Eigen::Vector3d( double ) >& memberOffsetAcceleration ):
        barycenter_( barycenter ), externalAccelerations_( externalAccelerations ), memberOffsetAcceleration_( memberOffsetAcceleration )
    {}

    void updateMembers( const double time = TUDAT_NAN ) override
    {
        if( time != currentTime_ )
        {
            currentAcceleration_ = -memberOffsetAcceleration_( time );
            for( const auto& entry : externalAccelerations_ )
            {
                for( const auto& model : entry.second )
                {
                    model->updateMembers( time );
                    currentAcceleration_ -= model->getAcceleration( );
                }
            }
            currentTime_ = time;
        }
    }
    void resetCurrentTime( ) override
    {
        AccelerationModel3d::resetCurrentTime( );
        for( const auto& entry : externalAccelerations_ )
        {
            for( const auto& model : entry.second )
            {
                model->resetCurrentTime( );
            }
        }
    }
    std::string getBarycenterName( ) const
    {
        return barycenter_;
    }
    const SingleBodyAccelerationMap& getExternalAccelerations( ) const
    {
        return externalAccelerations_;
    }

private:
    std::string barycenter_;
    SingleBodyAccelerationMap externalAccelerations_;
    std::function< Eigen::Vector3d( double ) > memberOffsetAcceleration_;
};

}  // namespace basic_astrodynamics
}  // namespace tudat
#endif
