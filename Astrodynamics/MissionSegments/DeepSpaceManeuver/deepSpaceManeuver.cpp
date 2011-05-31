/*! \file deepSpaceManeuver.cpp
 *    This source file contains a base class for deep space maneuver elements.
 *
 *    Path              : /Astrodynamics/MissionSegments/DeepSpaceManeuver/
 *    Version           : 1
 *    Check status      : Checked
 *
 *    Author            : E. Iorfida
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : elisabetta_iorfida@yahoo.it
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Date created      : 24 February, 2011
 *    Last modified     : 6 April, 2011
 *
 *    References
 *
 *    Notes
 *
 *    Copyright (c) 2010 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110224    E. Iorfida        First creation of code.
 *      110406    K. Kumar          Minor modifications.
 */

// Include statements.
#include "deepSpaceManeuver.h"

//! Default constructor.
DeepSpaceManeuver::DeepSpaceManeuver( ):
        deltaV_( -1.0 ),
        timeOfDeepSpaceManeuver_( 0.0 )
{
}

//! Default destructor.
DeepSpaceManeuver::~DeepSpaceManeuver( )
{
}

//! Set time of deep space maneuver event.
void DeepSpaceManeuver::setTime( const double &timeOfDeepSpaceManeuver )
{
    timeOfDeepSpaceManeuver_ = timeOfDeepSpaceManeuver;
}

//! Set state at deep space maneuver event.
void DeepSpaceManeuver::setState( State* pointerToState )
{
    pointerToState_ = pointerToState;
}

//! Set delta-V of deep space maneuver event.
void DeepSpaceManeuver::setDeltaV( const double &deltaV )
{
    deltaV_ = deltaV;
}

//! Get time of deep space maneuver event.
double& DeepSpaceManeuver::getTime( )
{
    return timeOfDeepSpaceManeuver_;
}

//! Get state at deep space maneuver event.
State* DeepSpaceManeuver::getState( )
{
    return pointerToState_;
}

//! Get delta-V of deep space maneuver event.
double& DeepSpaceManeuver::getDeltaV( )
{
    return deltaV_;
}

// End of file.
