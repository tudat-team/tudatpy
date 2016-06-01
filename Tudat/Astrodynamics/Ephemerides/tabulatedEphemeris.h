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
 *      YYMMDD    Author              Comment
 *      141105    D. Dirkx            File created.
 *
 *    References
 *
 *    Notes
 */

#ifndef TUDAT_TABULATEDEPHEMERIS_H
#define TUDAT_TABULATEDEPHEMERIS_H

#include <map>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/Astrodynamics/Ephemerides/ephemeris.h"
#include "Tudat/Mathematics/Interpolators/oneDimensionalInterpolator.h"


namespace tudat
{

namespace ephemerides
{

//! Class that determines an ephemeris from tabulated data.
/*!
 *  Class that determines an ephemeris from tabulated data, by using numerical interpolation of
 *  this data. Required input to this class is a OneDimensionalInterpolator, which may be reset.
 *  This class may for instance be used for setting the numerically integrated state of a body
 *  as its 'new' ephemeris
 */
template< typename StateScalarType = double, typename TimeType = double >
class TabulatedCartesianEphemeris : public Ephemeris
{
public:

    typedef Eigen::Matrix< StateScalarType, 6, 1 > StateType;

    //! Typedef for state interpolator
    typedef boost::shared_ptr< interpolators::OneDimensionalInterpolator
    < TimeType, StateType  > > StateInterpolatorPointer;

    //! Constructor, sets data interpolator and frame data.
    /*!
     *  Constructor, sets data interpolator and frame data.
     *  \param interpolator Interpolator that returns the interpolated state as a function of time.
     *  \param referenceFrameOrigin Origin of reference frame in which state is defined.
     *  \param referenceFrameOrientation Orientation of reference frame in which state is defined.
     *  \param julianDayAtEpoch Julian day for times to be used as input to interpolator.
     */
    TabulatedCartesianEphemeris(
            const StateInterpolatorPointer interpolator,
            const std::string referenceFrameOrigin = "SSB",
            const std::string referenceFrameOrientation = "ECLIPJ2000",
            const double julianDayAtEpoch = basic_astrodynamics::JULIAN_DAY_ON_J2000 ):
        Ephemeris( referenceFrameOrigin, referenceFrameOrientation ), interpolator_( interpolator ),
            julianDayAtEpoch_( julianDayAtEpoch )
    {  }

    //! Destructor
    /*!
     *  Destructor
     */
    ~TabulatedCartesianEphemeris( ){ }

    //! Function to reset the state interpolator.
    /*!
     *  Function to reset the state interpolator, for instance following an update of the states of
     *  the body after a new numerical integration.
     *  \param interpolator New interpolator that returns the interpolated state as a function of
     *  time.
     *  \param julianDayAtEpoch New reference Julian day for times to be used as input to
     *  interpolator.
     */
    void resetInterpolator( const StateInterpolatorPointer interpolator,
                            const double julianDayAtEpoch =
                                basic_astrodynamics::JULIAN_DAY_ON_J2000 )
    {
        interpolator_ = interpolator;
        julianDayAtEpoch_ = julianDayAtEpoch;
    }

    //! Get cartesian state from ephemeris.
    /*!
     * Returns cartesian state from ephemeris, as calculated from interpolator_.
     * \param secondsSinceEpoch Seconds since epoch.
     * \param julianDayAtEpoch Reference epoch in Julian day.
     * \return State in Cartesian elements from ephemeris.
     */
    basic_mathematics::Vector6d getCartesianStateFromEphemeris(
            const double secondsSinceEpoch,
            const double julianDayAtEpoch = basic_astrodynamics::JULIAN_DAY_ON_J2000 );

    //! Get cartesian state from ephemeris (in long double precision).
    /*!
     * Returns cartesian state from ephemeris  (in long double precision), as calculated from interpolator_. For
     * double StateScalarType class template argument, this function returns the double precision interpolated values,
     * cast to long double. Only for long double StateScalarType argument is this function used to its fullest.
     * \param secondsSinceEpoch Seconds since epoch.
     * \param julianDayAtEpoch Reference epoch in Julian day.
     * \return State in Cartesian elements from ephemeris.
     */
    Eigen::Matrix< long double, 6, 1 > getCartesianLongStateFromEphemeris(
            const double secondsSinceEpoch,
            const double julianDayAtEpoch = basic_astrodynamics::JULIAN_DAY_ON_J2000 );


    //! Function to return the interpolator
    /*!
     *  Function to return the interpolator that is to be used to calculate the state.
     *  \return Interpolator that is to be used to calculate the state.
     */
    StateInterpolatorPointer getInterpolator( )
    {
        return interpolator_;
    }


private:

    //! Interpolator that returns body state as a function of time.
    /*!
     *  Interpolator that returns body state as a function of time by calling the interpolate
     *  function (i.e. time as independent variable and states as dependent variables ).
     */
    StateInterpolatorPointer interpolator_;


    //! Reference Julian day for times to be used as input to interpolator_.
    /*!
     *  Reference Julian day for times to be used as input to interpolator_, i.e. the input
     *  argument to the interpolate() function of interpolator_ should be in seconds since
     *  the Julian day given by this variable.
     */
    double julianDayAtEpoch_;

};

} // namespace ephemerides

} // namespace tudat
#endif // TUDAT_TABULATEDEPHEMERIS_H
