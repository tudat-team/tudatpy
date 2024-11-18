/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "tudat/astro/gravitation/ringGravityModel.h"

namespace tudat
{
namespace gravitation
{

//void RingGravitationalAccelerationModel::updateMembers( const double currentTime )
//{
//    if( !( this->currentTime_ == currentTime ) )
//    {
//
//        rotationToIntegrationFrame_ = rotationFromBodyFixedToIntegrationFrameFunction_( );
//
//        subjectPositionFunction_( positionOfBodySubjectToAcceleration_ );
//        sourcePositionFunction_( positionOfBodyExertingAcceleration_ );
//        currentInertialRelativePosition_ =
//                positionOfBodySubjectToAcceleration_ - positionOfBodyExertingAcceleration_ ;
//
//        currentRelativePosition_ = rotationToIntegrationFrame_.inverse( ) * currentInertialRelativePosition_;
//
//        ringCache_->setRingRadius( ringRadiusFunction_( ) );
//        ringCache_->update( currentRelativePosition_ );
//
//        // Compute the current acceleration
//        currentAccelerationInBodyFixedFrame_ = computeRingGravitationalAcceleration(
//                currentRelativePosition_,
//                ringRadiusFunction_( ),
//                gravitationalParameterFunction_(),
//                ringCache_->getEllipticIntegralB( ),
//                ringCache_->getEllipticIntegralE( ),
//                ringCache_->getEllipticIntegralS( ) );
//
//        currentAcceleration_ = rotationToIntegrationFrame_ * currentAccelerationInBodyFixedFrame_;
//
//        // Compute the current gravitational potential
//        if ( updatePotential_ )
//        {
//            currentPotential_ = computeRingGravitationalPotential(
//                currentRelativePosition_,
//                ringRadiusFunction_( ),
//                gravitationalParameterFunction_(),
//                ringCache_->getEllipticIntegralK( ) );
//        }
//
//    }
//}


} // namespace gravitation

} // namespace tudat