/*    Copyright (c) 2010-2017, Delft University of Technology
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

#include <boost/lexical_cast.hpp>
#include <boost/exception/all.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

#include "Tudat/External/SpiceInterface/spiceEphemeris.h"

namespace tudat
{

namespace ephemerides
{

//! Constructor
SpiceEphemeris::SpiceEphemeris( const std::string& targetBodyName,
                                const std::string& observerBodyName,
                                const bool correctForStellarAbberation,
                                const bool correctForLightTimeAbberation,
                                const bool convergeLighTimeAbberation,
                                const std::string& referenceFrameName,
                                const double referenceJulianDay )
    : Ephemeris( observerBodyName, referenceFrameName ),
      targetBodyName_( targetBodyName )
{
    referenceDayOffSet_ = ( referenceJulianDay - basic_astrodynamics::JULIAN_DAY_ON_J2000 ) * physical_constants::JULIAN_DAY;

    // Check consistency of input.
    if ( correctForLightTimeAbberation == 0 && convergeLighTimeAbberation == 1 )
    {
        throw std::runtime_error(
                    "Error, requested multiple iterations for light time correction, but not light time correction itself." );
    }

    if ( correctForLightTimeAbberation == 0 && correctForStellarAbberation ==  1 )
    {

        throw std::runtime_error(
                    "Error, requested stellar aberration, but not light-time correction, this is not supprted by spice." );
    }

    // Set aberration corrections variable.
    abberationCorrections_ = "";
    if ( correctForLightTimeAbberation && !convergeLighTimeAbberation )
    {
        abberationCorrections_.append( "LT" );
    }

    else if ( correctForLightTimeAbberation && convergeLighTimeAbberation )
    {
        abberationCorrections_.append( "CN" );
    }

    else if ( !correctForLightTimeAbberation && !correctForStellarAbberation )
    {
        abberationCorrections_.append( "NONE" );
    }

    if ( correctForLightTimeAbberation && correctForStellarAbberation )
    {
        abberationCorrections_.append( "+S" );
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
                abberationCorrections_, ephemerisTime + referenceDayOffSet_ );

    return cartesianStateAtEpoch;
}

} // namespace ephemerides
} // namespace tudat
