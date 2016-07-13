#ifndef THRUSTACCELERATIONMODEL_H
#define THRUSTACCELERATIONMODEL_H

#include <limits>

#include <boost/function.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/massRateModel.h"

namespace tudat
{

namespace propulsion
{

class ThrustAcceleration : public basic_astrodynamics::AccelerationModel< Eigen::Vector3d >
{
public:
    ThrustAcceleration(
            const boost::function< double( ) > thrustMagnitudeFunction,
            const boost::function< Eigen::Vector3d( ) > thrustDirectionFunction,
            const boost::function< double( ) > bodyMassFunction,
            const boost::function< double( ) > massRateFunction,
            const std::string associatedThroustSource,
            const boost::function< void( const double ) > thrustUpdateFunction ):
        AccelerationModel< Eigen::Vector3d >( ),
        thrustMagnitudeFunction_( thrustMagnitudeFunction ),
        thrustDirectionFunction_( thrustDirectionFunction ),
        bodyMassFunction_( bodyMassFunction ),
        massRateFunction_( massRateFunction ),
        associatedThroustSource_( associatedThroustSource ),
        thrustUpdateFunction_( thrustUpdateFunction ){ }

    ~ThrustAcceleration( ){ }

    Eigen::Vector3d getAcceleration( )
    {
        return currentAcceleration_;
    }

    void updateMembers( const double currentTime = TUDAT_NAN )
    {
        if( !( currentTime_ == currentTime ) )
        {
            thrustUpdateFunction_( currentTime );

            currentAccelerationDirection_ = thrustDirectionFunction_( );

            if( ( std::fabs( currentAccelerationDirection_.norm( ) ) - 1.0 ) > 5.0 * std::numeric_limits< double >::epsilon( ) )
            {
                throw std::runtime_error( "Error in thrust acceleration, direction is not a unit vector"  );
            }

            currentThrustMagnitude_ = thrustMagnitudeFunction_( );
            currentMassRate_ = -massRateFunction_( );

            currentAcceleration_ = currentAccelerationDirection_ * currentThrustMagnitude_ / bodyMassFunction_( );
            currentTime_ = currentTime;
        }

    }

    double getCurrentMassRate( )
    {
        return currentMassRate_;
    }

    std::string getAssociatedThroustSource( )
    {
        return associatedThroustSource_;
    }

private:

    boost::function< double( ) > thrustMagnitudeFunction_;

    boost::function< Eigen::Vector3d( ) > thrustDirectionFunction_;

    boost::function< double( ) > bodyMassFunction_;

    boost::function< double( ) > massRateFunction_;

    std::string associatedThroustSource_;

    boost::function< void( const double ) > thrustUpdateFunction_;

    Eigen::Vector3d currentAcceleration_;

    Eigen::Vector3d currentAccelerationDirection_;

    double currentThrustMagnitude_;

    double currentMassRate_;
};

}

}
#endif // THRUSTACCELERATIONMODEL_H
