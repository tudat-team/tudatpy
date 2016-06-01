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
 *      110128    K. Kumar          File created.
 *      110221    K. Kumar          Updated code to work as base class for derived ephemeris
 *                                  classes.
 *      120322    D. Dirkx          Modified to new Ephemeris interfaces.
 *      130120    K. Kumar          Updated VectorXd to Vector6d; added shared-ptr typedef.
 *      130120    D. Dirkx          Updated with new Julian day + seconds since Julian day input.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_EPHEMERIS_H
#define TUDAT_EPHEMERIS_H

#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>

#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"

namespace tudat
{
namespace ephemerides
{


//! Ephemeris base class.
/*!
 * Ephemeris base class.
 */
class Ephemeris
{
public:

    //! Constructor
    /*!
     *  Constructor, sets reference frame associated with ephemeris (default empty).
     *  \param referenceFrameOrigin Origin of reference frame (string identifier).
     *  \param referenceFrameOrientation Orientation of reference frame (string identifier).
     */
    Ephemeris( const std::string& referenceFrameOrigin = "",
               const std::string& referenceFrameOrientation = "" ):
        referenceFrameOrigin_( referenceFrameOrigin ),
        referenceFrameOrientation_( referenceFrameOrientation )
    { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~Ephemeris( ) { }

    //! Get state from ephemeris.
    /*!
     * Returns state from ephemeris at given Julian date.
     * \param secondsSinceEpoch Seconds since epoch.
     * \param julianDayAtEpoch Reference epoch in Julian day.
     * \return State from ephemeris.
     */
    virtual basic_mathematics::Vector6d getCartesianStateFromEphemeris(
            const double secondsSinceEpoch,
            const double julianDayAtEpoch = basic_astrodynamics::JULIAN_DAY_ON_J2000 ) = 0;

    //! Get state from ephemeris (with long double as state scalar).
    /*!
     * Returns state from ephemeris with long double as state scalar at given time. By default, this
     * function casts the double getCartesianStateFromEphemeris to long double. It may be overridden
     * by derived classes to make use of full long double computations.
     * \param secondsSinceEpoch Seconds since epoch at which ephemeris is to be evaluated.
     * \param julianDayAtEpoch Reference epoch in Julian day.
     * \return State from ephemeris with long double as state scalar
     */
    virtual Eigen::Matrix< long double, 6, 1 > getCartesianLongStateFromEphemeris(
            const double secondsSinceEpoch,
            const double julianDayAtEpoch = basic_astrodynamics::JULIAN_DAY_ON_J2000 )
    {
        return getCartesianStateFromEphemeris( secondsSinceEpoch, julianDayAtEpoch ).cast< long double >( );
    }

    //! Get state from ephemeris, with state scalar as template type.
    /*!
     * Returns state from ephemeris (state scalar as template type) at given time.
     * \param time Time at which ephemeris is to be evaluated (JULIAN_DAY_ON_J2000 used as reference
              julian day when needed).
     * \return State from ephemeris with requested state scalar type.
     */
    template< typename StateScalarType, typename TimeType >
    Eigen::Matrix< StateScalarType, 6, 1 > getTemplatedStateFromEphemeris( const TimeType& time );

    //! Get reference frame origin.
    /*!
     * Returns reference frame origin as a string.
     * \return Reference frame origin.
     */
    std::string getReferenceFrameOrigin( ) { return referenceFrameOrigin_; }

    //! Get reference frame orientation.
    /*!
     * Returns the reference frame orientation as a string.
     * \return Reference frame orientation
     */
    std::string getReferenceFrameOrientation( ) { return referenceFrameOrientation_; }

protected:

    //! Reference frame origin.
    /*!
     * Reference frame origin. This identifier gives only the origin of the reference frame,
     * its orientation is defined by the referenceFrameOrientation_ variable.
     */
    std::string referenceFrameOrigin_;

    //! Reference frame orientation
    /*!
     * Reference frame orientation. This identifier gives only the orientation of the
     * reference frame, the origin is defined by the referenceFrameOrigin_ variable.
     */
    std::string referenceFrameOrientation_;
};

//! Typedef for shared-pointer to Ephemeris object.
typedef boost::shared_ptr< Ephemeris > EphemerisPointer;

} // namespace ephemerides
} // namespace tudat

#endif // TUDAT_EPHEMERIS_H
