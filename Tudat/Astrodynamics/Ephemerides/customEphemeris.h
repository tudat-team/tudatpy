/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#ifndef TUDAT_CUSTOMEPHEMERIS_H
#define TUDAT_CUSTOMEPHEMERIS_H

#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>

#include "Tudat/Astrodynamics/Ephemerides/ephemeris.h"


namespace tudat
{

namespace ephemerides
{

//! Ephemeris class that gives a custom (i.e. arbitrarily defined as a function of time) state.
class CustomEphemeris : public Ephemeris
{
public:

    using Ephemeris::getCartesianState;

    //! Constructor.
    /*!
     *  Constructor of an Ephemeris object from a custom function
     *  \param stateFunction Function returning the state as a function of time
     *  \param referenceFrameOrigin Origin of reference frame in which state is defined.
     *  \param referenceFrameOrientation Orientation of reference frame in which state is defined.
     */
    CustomEphemeris( const boost::function< Eigen::Vector6d( const double ) > stateFunction,
                     const std::string& referenceFrameOrigin = "SSB",
                     const std::string& referenceFrameOrientation = "ECLIPJ2000" ):
        Ephemeris( referenceFrameOrigin, referenceFrameOrientation ),
        stateFunction_( stateFunction ) { }

    //! Get state from ephemeris according to custom function
    /*!
     * Returns state from ephemeris at given time.
     * \param seconsSinceEpoch Seconds since epoch at which ephemeris is to be evaluated
              (not used in this derived class)
     * \return State given by stateFunction_
     */
    Eigen::Vector6d getCartesianState(
            const double seconsSinceEpoch = 0.0 )
    {
        return stateFunction_( seconsSinceEpoch );
    }

private:

    //! Time-independent state function.
    /*!
     *  Function that returns a constant cartesian state.
     */
    boost::function< Eigen::Vector6d( const double ) > stateFunction_;

};

} // namespace ephemerides

} // namespace tudat

#endif // TUDAT_CUSTOMEPHEMERIS_H
