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
 *      150501    D. Dirkx          Ported from personal code
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_DEFAULTBODIES_H
#define TUDAT_DEFAULTBODIES_H

#include "Tudat/SimulationSetup/createBodies.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to create default settings for a body's atmosphere model.
/*!
 *  Function to create default settings for a body's atmosphere model. Currently set to no
 *  atmosphere, except for Earth, for which a tabulated version of the 1976 Standard Atmosphere is
 *  set.
 *  \param bodyName Name of body for which default atmosphere settings are to be retrieved.
 *  \param initialTime Start time at which environment models in body are to be created
 *  (not currently used by this function, but included for consistency).
 *  \param finalTime End time up to which environment models in body are to be created
 *  (not currently used by this function, but included for consistency).
 */
boost::shared_ptr< AtmosphereSettings > getDefaultAtmosphereModelSettings(
        const std::string& bodyName,
        const double initialTime,
        const double finalTime );

//! Function to create default settings for a body's ephemeris.
/*!
 *  Function to create default settings for a body's ephemeris. Currently set to a
 *  creating a 6th order Lagrange interpolator from Spice, with a 300 s time step.
 *  \param bodyName Name of body for which default ephemeris settings are to be retrieved.
 *  \param initialTime Start time at which ephemeris is to be created
 *  \param finalTime End time up to which ephemeris is to be created
 */
boost::shared_ptr< EphemerisSettings > getDefaultEphemerisSettings(
        const std::string& bodyName,
        const double initialTime,
        const double finalTime );

//! Function to create default settings for a body's gravity field model.
/*!
 *  Function to create default settings for a body's gravity field model. Currently set to
 *  a point mass gravty field, with the gravitational parameter obtained from Spice.
 *  \param bodyName Name of body for which default gravity field settings are to be retrieved.
 *  \param initialTime Start time at which environment models in body are to be created
 *  (not currently used by this function, but included for consistency).
 *  \param finalTime End time up to which environment models in body are to be created
 *  (not currently used by this function, but included for consistency).
 */
boost::shared_ptr< GravityFieldSettings > getDefaultGravityFieldSettings(
        const std::string& bodyName,
        const double initialTime,
        const double finalTime );

//! Function to create default settings for a body's rotation model.
/*!
 *  Function to create default settings for a body's rotation model. Currently set to
 *  a rotation model taken directly from Spice
 *  \param bodyName Name of body for which default rotation model settings are to be retrieved.
 *  \param initialTime Start time at which environment models in body are to be created
 *  (not currently used by this function, but included for consistency).
 *  \param finalTime End time up to which environment models in body are to be created
 *  (not currently used by this function, but included for consistency).
 */
boost::shared_ptr< RotationModelSettings > getDefaultRotationModelSettings(
        const std::string& bodyName,
        const double initialTime,
        const double finalTime );

//! Function to create default settings for a body's shape model.
/*!
 *  Function to create default settings for a body's shape model. Currently set to
 *  a spherical model, with the radius taken from Spice
 *  \param body Name of body for which default shape model settings are to be retrieved.
 *  \param initialTime Start time at which environment models in body are to be created
 *  (not currently used by this function, but included for consistency).
 *  \param finalTime End time up to which environment models in body are to be created
 *  (not currently used by this function, but included for consistency).
 */
boost::shared_ptr< BodyShapeSettings > getDefaultBodyShapeSettings(
        const std::string& body,
        const double initialTime, const double finalTime );


//! Function to create default settings from which to create a single body object.
/*!
 *  Function to create default settings from which to create a single body object using
 *  the code in createBodies.h/.cpp. This function is included to streamline and simplify the
 *  creation of typical celestial bodies. The default settings for the various
 *  environment models of the body are defined in the various functions defined in this file.
 *  \param bodyName Name of body for which default settings are to be retrieved.
 *  \param initialTime Start time at which environment models in body are to be created
 *  (included as some environment models require e.g., interpolators to be created over
 *  a certain time period).
 *  \param finalTime End time up to which environment models in body are to be created
 *  (included as some environment models require e.g., interpolators to be created over
 *  a certain time period).
 */
boost::shared_ptr< RotationModelSettings > getDefaultRotationModelSettings(
        const std::string& bodyName,
        const double initialTime,
        const double finalTime );

//! Function to create default settings for a single for body.
/*!
 *  Function to create default settings for a single body from which to create a object using
 *  the code in createBodies.h/.cpp. This function is included to streamline and simplify the
 *  creation of typical celestial bodies. The default settings for the various
 *  environment models of the body are defined in the various functions defined in this file.
 *  \param body Name of body for which default settings are to be retrieved.
 *  \param initialTime Start time at which environment models in body are to be created
 *  (included as some environment models require e.g., interpolators to be created over
 *  a certain time period).
 *  \param finalTime End time up to which environment models in body are to be created
 *  (included as some environment models require e.g., interpolators to be created over
 *  a certain time period).
 */
boost::shared_ptr< BodySettings > getDefaultSingleBodySettings(
        const std::string& body,
        const double initialTime,
        const double finalTime );

//! Function to create default settings from which to create a set of body objects.
/*!
 *  Function to create default settings from which to create a set of body objects using
 *  the code in createBodies.h/.cpp. This function is included to streamline and simplify the
 *  creation of typical celestial bodies. The default settings for the various
 *  environment models of the body are defined in the various functions defined in this file.
 *  \param bodies List of bodies for which default settings are to be retrieved.
 *  \param initialTime Start time at which environment models in body are to be created
 *  (included as some environment models require e.g., interpolators to be created over
 *  a certain time period).
 *  \param finalTime End time up to which environment models in body are to be created
 *  (included as some environment models require e.g., interpolators to be created over
 *  a certain time period).
 */
std::map< std::string, boost::shared_ptr< BodySettings > > getDefaultBodySettings(
        const std::vector< std::string >& bodies,
        const double initialTime,
        const double finalTime );

}

}

#endif // TUDAT_DEFAULTBODIES_H
