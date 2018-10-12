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

#ifndef TUDAT_EPHEMERIS_H
#define TUDAT_EPHEMERIS_H

#include <memory>
#include <functional>

#include "Tudat/Astrodynamics/Ephemerides/ephemeris.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"
#include "Tudat/Basics/basicTypedefs.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/Basics/timeType.h"

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
     * \param secondsSinceEpoch Seconds since epoch at which ephemeris is to be evaluated.
     * \return State from ephemeris.
     */
    virtual Eigen::Vector6d getCartesianState(
            const double secondsSinceEpoch ) = 0;

    //! Get state from ephemeris (with long double as state scalar).
    /*!
     * Returns state from ephemeris with long double as state scalar at given time. By default, this
     * function casts the double getCartesianState to long double. It may be overridden
     * by derived classes to make use of full long double computations.
     * \param secondsSinceEpoch Seconds since epoch at which ephemeris is to be evaluated.
     * \return State from ephemeris with long double as state scalar
     */
    virtual Eigen::Matrix< long double, 6, 1 > getCartesianLongState(
            const double secondsSinceEpoch )
    {
        return getCartesianState( secondsSinceEpoch ).cast< long double >( );
    }


    //! Get state from ephemeris (with double as state scalar and Time as time type).
    /*!
     * Returns state from ephemeris with double as state scalar at given time (as custom Time type). By default, this
     * function casts the double getCartesianState to double. It may be overridden
     * by derived classes.
     * \param currentTime Time at which state is to be evaluated
     * \return State from ephemeris with double as state scalar
     */
    virtual Eigen::Vector6d getCartesianStateFromExtendedTime(
            const Time& currentTime )
    {
        return getCartesianState( currentTime.getSeconds< double >( ) );
    }

    //! Get state from ephemeris (with long double as state scalar and Time as time type).
    /*!
     * Returns state from ephemeris with long double as state scalar at given time (as custom Time type). By default, this
     * function casts the double getCartesianState to long double. It may be overridden
     * by derived classes to make use of full long double computations.
     * \param currentTime Time at which state is to be evaluated
     * \return State from ephemeris with long double as state scalar
     */
    virtual Eigen::Matrix< long double, 6, 1 > getCartesianLongStateFromExtendedTime(
            const Time& currentTime )
    {
        return getCartesianLongState( currentTime.getSeconds< double >( ) );
    }

    //! Get state from ephemeris, with state scalar as template type.
    /*!
     * Returns state from ephemeris (state scalar as template type) at given time.
     * \param time Time at which ephemeris is to be evaluated
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
typedef std::shared_ptr< Ephemeris > EphemerisPointer;

//! Function to compute the relative state from two state functions.
/*!
 *  Function to compute the relative state from two state functions.
 *  \param relativeState Relative state returned by stateFunctionOfBody w.r.t. state returned by stateFunctionOfCentralBody
 *  (returned by reference).
 *  \param stateFunctionOfBody Function returning state of body for which relative state is to be computed.
 *  \param stateFunctionOfCentralBody Function returning state of central body w.r.t. which the relative state is to be
 *  computed.
 */
void getRelativeState(
        Eigen::Vector6d& relativeState,
        const std::function< Eigen::Vector6d( ) > stateFunctionOfBody,
        const std::function< Eigen::Vector6d( ) > stateFunctionOfCentralBody );

} // namespace ephemerides

} // namespace tudat

#endif // TUDAT_EPHEMERIS_H
