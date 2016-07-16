/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef THRUSTACCELERATIONMODEL_H
#define THRUSTACCELERATIONMODEL_H

#include <limits>

#include <boost/function.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/massRateModel.h"
#include "Tudat/Astrodynamics/Propagators/environmentUpdateTypes.h"

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
            const boost::function< void( const double ) > thrustUpdateFunction,
            const std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >& requiredModelUpdates =
            std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >( ) ):
        AccelerationModel< Eigen::Vector3d >( ),
        thrustMagnitudeFunction_( thrustMagnitudeFunction ),
        thrustDirectionFunction_( thrustDirectionFunction ),
        bodyMassFunction_( bodyMassFunction ),
        massRateFunction_( massRateFunction ),
        associatedThroustSource_( associatedThroustSource ),
        thrustUpdateFunction_( thrustUpdateFunction ),
        requiredModelUpdates_( requiredModelUpdates ){ }

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

            if( ( std::fabs( currentAccelerationDirection_.norm( ) ) - 1.0 ) > 10.0 * std::numeric_limits< double >::epsilon( ) )
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

    std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > getRequiredModelUpdates( )
    {
        return requiredModelUpdates_;
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

    std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > requiredModelUpdates_;
};

} // namespace propulsion

} // namespace tudat

#endif // THRUSTACCELERATIONMODEL_H
