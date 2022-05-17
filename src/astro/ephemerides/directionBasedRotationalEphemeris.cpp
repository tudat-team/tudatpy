/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "tudat/astro/ephemerides/directionBasedRotationalEphemeris.h"

namespace tudat
{


namespace ephemerides
{

void DirectionBasedRotationalEphemeris::update( const double currentTime )
{
    if( !( currentTime == currentTime_ ) )
    {
        if( currentTime == currentTime )
        {
            currentTime_ = currentTime;
            currentInertialDirection_ = inertialBodyAxisDirectionFunction_( currentTime_ ).normalized( );
        }
        else
        {
            resetCurrentTime( );
        }
    }
}


void DirectionBasedRotationalEphemeris::resetCurrentTime(  )
{
    currentTime_ = TUDAT_NAN;
    currentInertialDirection_.setConstant( TUDAT_NAN );

    currentEulerAnglesTime_ = TUDAT_NAN;
    eulerAngles_.setConstant( TUDAT_NAN );
}

} // namespace tudat
} // namespace ephemerides
