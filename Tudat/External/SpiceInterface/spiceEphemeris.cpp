/*    Copyright (c) 2010-2015, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120717    D. Dirkx          Creation of file.
 *      130120    D. Dirkx          Updated with new Julian day + seconds since Julian day input.
 *      130226    K. Kumar          Updated return-type for getCartesianState().
 *
 *    References
 *
 *    Notes
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
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error(
                            "Warning, requested multiple iterations for light time correction, "
                            "but not light time correction itself." ) ) );
    }

    if ( correctForLightTimeAbberation == 0 && correctForStellarAbberation ==  1 )
    {

        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error
                        ( "Warning, requested stellar aberration, but not light-time correction, "
                          "this is not supprted by spice." ) ) );
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
basic_mathematics::Vector6d SpiceEphemeris::getCartesianState(
        const double secondsSinceEpoch )
{
    using namespace basic_astrodynamics;

    // Retrieve body state at given ephemeris time, using settings passed to constructor of this
    // object.

    // Calculate ephemeris time at which cartesian state is to be determind.
    const double ephemerisTime = secondsSinceEpoch;

    // Retrieve Cartesian state from spice.
    const basic_mathematics::Vector6d cartesianStateAtEpoch =
            spice_interface::getBodyCartesianStateAtEpoch(
                targetBodyName_, referenceFrameOrigin_, referenceFrameOrientation_,
                abberationCorrections_, ephemerisTime + referenceDayOffSet_ );

    return cartesianStateAtEpoch;
}

} // namespace ephemerides
} // namespace tudat
