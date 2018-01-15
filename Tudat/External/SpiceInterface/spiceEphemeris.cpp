/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include <stdexcept>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

#include "Tudat/External/SpiceInterface/spiceEphemeris.h"

namespace tudat
{

namespace ephemerides
{

//! Constructor
SpiceEphemeris::SpiceEphemeris( const std::string& targetBodyName,
                                const std::string& observerBodyName,
                                const bool correctForStellarAberration,
                                const bool correctForLightTimeAberration,
                                const bool convergeLighTimeAberration,
                                const std::string& referenceFrameName,
                                const double referenceJulianDay )
    : Ephemeris( observerBodyName, referenceFrameName ),
      targetBodyName_( targetBodyName )
{
    referenceDayOffSet_ = ( referenceJulianDay - basic_astrodynamics::JULIAN_DAY_ON_J2000 ) * physical_constants::JULIAN_DAY;

    // Check consistency of input.
    if ( correctForLightTimeAberration == 0 && convergeLighTimeAberration == 1 )
    {
        throw std::runtime_error(
                    "Error, requested multiple iterations for light time correction, but not light time correction itself." );
    }

    if ( correctForLightTimeAberration == 0 && correctForStellarAberration ==  1 )
    {

        throw std::runtime_error(
                    "Error, requested stellar aberration, but not light-time correction, this is not supprted by spice." );
    }

    // Set aberration corrections variable.
    aberrationCorrections_ = "";
    if ( correctForLightTimeAberration && !convergeLighTimeAberration )
    {
        aberrationCorrections_.append( "LT" );
    }

    else if ( correctForLightTimeAberration && convergeLighTimeAberration )
    {
        aberrationCorrections_.append( "CN" );
    }

    else if ( !correctForLightTimeAberration && !correctForStellarAberration )
    {
        aberrationCorrections_.append( "NONE" );
    }

    if ( correctForLightTimeAberration && correctForStellarAberration )
    {
        aberrationCorrections_.append( " +S" );
    }
}

//! Get Cartesian state from ephemeris.
Eigen::Vector6d SpiceEphemeris::getCartesianState(
        const double secondsSinceEpoch )
{
    using namespace basic_astrodynamics;

    // Retrieve body state at given ephemeris time, using settings passed to constructor of this
    // object.

    // Calculate ephemeris time at which cartesian state is to be determind.
    const double ephemerisTime = secondsSinceEpoch;

    // Retrieve Cartesian state from spice.
    const Eigen::Vector6d cartesianStateAtEpoch =
            spice_interface::getBodyCartesianStateAtEpoch(
                targetBodyName_, referenceFrameOrigin_, referenceFrameOrientation_,
                aberrationCorrections_, ephemerisTime + referenceDayOffSet_ );

    return cartesianStateAtEpoch;
}

} // namespace ephemerides
} // namespace tudat
