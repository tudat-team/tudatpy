#ifndef THRUSTACCELERATIONMODEL_H
#define THRUSTACCELERATIONMODEL_H

#include <limits>

#include <boost/function.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"


namespace tudat
{

namespace basic_astrodynamics
{

class ThrustAcceleration : public basic_astrodynamics::AccelerationModel< Eigen::Vector3d >
{
public:
    ThrustAcceleration(
            const boost::function< double( ) > thrustMagnitudeFunction,
            const boost::function< Eigen::Vector3d( ) > thrustDirectionFunction,
            const boost::function< double( ) > bodyMassFunction,
            const boost::function< void( const double ) > thrustUpdateFunction ):
        AccelerationModel< Eigen::Vector3d >( ),
    thrustMagnitudeFunction_( thrustMagnitudeFunction ),
    thrustDirectionFunction_( thrustDirectionFunction ),
    bodyMassFunction_( bodyMassFunction ),
    thrustUpdateFunction_( thrustUpdateFunction ){ }

    ~ThrustAcceleration( ){ }

    Eigen::Vector3d getAcceleration( )
    {
        return currentAcceleration_;
    }

    void updateMembers( const double currentTime = TUDAT_NAN )
    {
        thrustUpdateFunction_( currentTime );

        currentAccelerationDirection_ = thrustDirectionFunction_( );

        if( ( std::fabs( currentAccelerationDirection_.norm( ) ) - 1.0 ) > 5.0 * std::numeric_limits< double >::epsilon( ) )
        {
            throw std::runtime_error( "Error in thrust acceleration, direction is not a unit vector"  );
        }

        currentAcceleration_ = currentAccelerationDirection_ * thrustMagnitudeFunction_( ) / bodyMassFunction_( );
    }

private:

    boost::function< double( ) > thrustMagnitudeFunction_;

    boost::function< Eigen::Vector3d( ) > thrustDirectionFunction_;

    boost::function< double( ) > bodyMassFunction_;

    boost::function< void( const double ) > thrustUpdateFunction_;

    Eigen::Vector3d currentAcceleration_;

    Eigen::Vector3d currentAccelerationDirection_;
};

}

}
#endif // THRUSTACCELERATIONMODEL_H
