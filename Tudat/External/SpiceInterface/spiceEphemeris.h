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

#ifndef TUDAT_SPICE_EPHEMERIS_H
#define TUDAT_SPICE_EPHEMERIS_H

#include <string>

#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/Astrodynamics/Ephemerides/ephemeris.h"

#include "Tudat/External/SpiceInterface/spiceInterface.h"

#include "Tudat/Basics/basicTypedefs.h"

namespace tudat
{
namespace ephemerides
{

//! Ephemeris derived class which retrieves the state of a body directly from the SPICE library.
/*!
 *  Ephemeris derived class which retrieves the state of a body directly from the SPICE library.
 *  The body of which the ephemeris is to be retrieved, as well as the origin and orientation
 *  of the reference frame in which the states are returned, and any corrections that are
 *  applied, are defined once during object construction.
 */
class SpiceEphemeris : public Ephemeris
{
public:

    using Ephemeris::getCartesianState;

    //! Constructor.
    /*!
     * Constructor, sets the input variables for the calls to the spice function to retrieve state.
     * \param targetBodyName Name of body of which the ephemeris is to be calculated.
     * \param observerBodyName Name of body relative to which the ephemeris is to be calculated.
     * \param correctForStellarAberration Boolean whether to correct for stellar Aberration in
     *          retrieved values of (observed state).
     * \param correctForLightTimeAberration Boolean whether to correct for light time in
     *          retrieved values of (observed state).
     * \param convergeLighTimeAberration Boolean whether to use single iteration or max. 3
     *          iterations for calculating light time.
     * \param referenceFrameName Name of the reference frame in which the epehemeris is to be
     *          calculated.
     * \param referenceJulianDay Reference julian day w.r.t. which ephemeris is evaluated.
     */
    SpiceEphemeris( const std::string& targetBodyName, const std::string& observerBodyName,
                    const bool correctForStellarAberration = true,
                    const bool correctForLightTimeAberration = true,
                    const bool convergeLighTimeAberration = false,
                    const std::string& referenceFrameName = "ECLIPJ2000",
                    const double referenceJulianDay = basic_astrodynamics::JULIAN_DAY_ON_J2000 );

    //! Get Cartesian state from ephemeris.
    /*!
     * Returns Cartesian state from ephemeris at given Julian day.
     * \param secondsSinceEpoch Seconds since epoch at which ephemeris is to be evaluated..
     * \return State from ephemeris.
     */
    Eigen::Vector6d getCartesianState(
            const double secondsSinceEpoch );

private:

    //! Name of body of which ephemeris is to be determined
    /*!
     * Name of body of which ephemeris is to be determined. Name can be either normal name
     * (example Jupiter, not case-sensitive), or NAIF ID number used internally in Spice.
     */
    std::string targetBodyName_;

    //! Name of body w.r.t. which ephemeris is to be determined.
    /*!
     * Name of body w.r.t. which ephemeris is to be determined (i.e. where reference frame
     * is to be centered). Name can be either normal name (example Jupiter, not case-sensitive),
     * or NAIF ID number used internally in Spice.
     */
    std::string observerBodyName_;

    //! Name of the reference frame in which ephemeris is to be expressed.
    /*!
     * Name of the reference frame in which ephemeris is to be expressed. This identifier gives
     * only the orientation of the reference frame, the origin is defined by the observerBodyName_
     * variable.
     */
    std::string referenceFrameName_;

    //! Aberration corrections that are to be applied when retrieving ephemeris.
    /*!
     * Aberration corrections that are to be applied when retrieving ephemeris. Stellar and light-
     * time aberration can be corrected for. Note that only the finite speed of light is included
     * in the light-time aberration, not general relativistic delay or path bending. The variable
     * defines the corrections the same as the required input to the spkezr spice function
     * ( see corresponding spice documentation ).
     */
    std::string aberrationCorrections_;

    //! Offset of reference julian day (from J2000) w.r.t. which ephemeris is evaluated.
    double referenceDayOffSet_;
};

} // namespace ephemerides
} // namespace tudat

#endif // TUDAT_SPICE_EPHEMERIS_H
